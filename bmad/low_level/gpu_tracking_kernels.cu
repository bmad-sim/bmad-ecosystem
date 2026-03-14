/*
 * gpu_tracking_kernels.cu
 *
 * CUDA kernels for GPU-accelerated particle tracking through supported
 * elements (drift, quadrupole, sbend, lcavity).  Called from Fortran
 * via iso_c_binding wrappers defined in gpu_tracking_mod.f90.
 *
 * Build requirements:
 *   - CUDA Toolkit (cuda_runtime.h)
 *   - Compile with nvcc: -DUSE_GPU_TRACKING
 *   - Link with: -lcudart
 *
 * The wrapper caches device memory allocations to avoid re-allocation
 * on repeated calls.
 */

#ifdef USE_GPU_TRACKING

#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/* --------------------------------------------------------------------------
 * CUDA error checking macro.  Returns -1 from the calling function on error.
 * Host wrappers (void) check the return value and issue a bare "return;".
 * -------------------------------------------------------------------------- */
#define CUDA_CHECK(call) do { \
    cudaError_t err_ = (call); \
    if (err_ != cudaSuccess) { \
        fprintf(stderr, "[gpu_tracking] CUDA error: %s at %s:%d\n", \
                cudaGetErrorString(err_), __FILE__, __LINE__); \
        return -1; \
    } \
} while(0)

/* Same as CUDA_CHECK but for void functions (host wrappers) */
#define CUDA_CHECK_VOID(call) do { \
    cudaError_t err_ = (call); \
    if (err_ != cudaSuccess) { \
        fprintf(stderr, "[gpu_tracking] CUDA error: %s at %s:%d\n", \
                cudaGetErrorString(err_), __FILE__, __LINE__); \
        return; \
    } \
} while(0)

/* Physical constants (must match Bmad values exactly) */
#define C_LIGHT   2.99792458e8
#define ALIVE_ST  1   /* alive$ */
#define LOST_PZ   8   /* lost_pz$ */

/* --------------------------------------------------------------------------
 * Cached device buffers
 * -------------------------------------------------------------------------- */
static double *d_vec[6]  = {NULL,NULL,NULL,NULL,NULL,NULL};
static int    *d_state   = NULL;
static double *d_beta    = NULL;
static double *d_p0c     = NULL;
static double *d_s       = NULL;
static double *d_t       = NULL;
static int     d_cap     = 0;          /* allocated capacity (particles) */

/* Cached device buffers for multipole data */
static double *d_a2  = NULL;
static double *d_b2  = NULL;
static double *d_ea2 = NULL;
static double *d_eb2 = NULL;
static double *d_cm  = NULL;

/* --------------------------------------------------------------------------
 * ensure_buffers — (re-)allocate device arrays when size changes
 * -------------------------------------------------------------------------- */
static int ensure_buffers(int n)
{
    if (n <= 0) return (n == 0) ? 0 : -1;
    if (n <= d_cap) return 0;

    /* Free old */
    for (int k = 0; k < 6; k++) { if (d_vec[k]) cudaFree(d_vec[k]); d_vec[k] = NULL; }
    if (d_state) cudaFree(d_state); d_state = NULL;
    if (d_beta)  cudaFree(d_beta);  d_beta  = NULL;
    if (d_p0c)   cudaFree(d_p0c);   d_p0c   = NULL;
    if (d_s)     cudaFree(d_s);     d_s     = NULL;
    if (d_t)     cudaFree(d_t);     d_t     = NULL;

    size_t db = (size_t)n * sizeof(double);
    size_t ib = (size_t)n * sizeof(int);

    for (int k = 0; k < 6; k++) {
        if (cudaMalloc((void**)&d_vec[k], db) != cudaSuccess) goto fail;
    }
    if (cudaMalloc((void**)&d_state, ib) != cudaSuccess) goto fail;
    if (cudaMalloc((void**)&d_beta,  db) != cudaSuccess) goto fail;
    if (cudaMalloc((void**)&d_p0c,   db) != cudaSuccess) goto fail;
    if (cudaMalloc((void**)&d_s,     db) != cudaSuccess) goto fail;
    if (cudaMalloc((void**)&d_t,     db) != cudaSuccess) goto fail;

    d_cap = n;
    return 0;

fail:
    fprintf(stderr, "[gpu_tracking] cudaMalloc failed for %d particles\n", n);
    /* Clean up any partially allocated buffers */
    for (int k = 0; k < 6; k++) { if (d_vec[k]) cudaFree(d_vec[k]); d_vec[k] = NULL; }
    if (d_state) cudaFree(d_state); d_state = NULL;
    if (d_beta)  cudaFree(d_beta);  d_beta  = NULL;
    if (d_p0c)   cudaFree(d_p0c);   d_p0c   = NULL;
    if (d_s)     cudaFree(d_s);     d_s     = NULL;
    if (d_t)     cudaFree(d_t);     d_t     = NULL;
    d_cap = 0;
    return -1;
}

/* --------------------------------------------------------------------------
 * gpu_tracking_available — query whether a CUDA GPU is present
 * -------------------------------------------------------------------------- */
extern "C" int gpu_tracking_available_(void)
{
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    return (err == cudaSuccess && count > 0) ? 1 : 0;
}

/* --------------------------------------------------------------------------
 * Numerically stable sqrt(1+x) - 1
 * -------------------------------------------------------------------------- */
__device__ __forceinline__ double sqrt_one_dev(double x)
{
    double sq = sqrt(1.0 + x);
    return x / (sq + 1.0);
}

/* --------------------------------------------------------------------------
 * Low energy z correction (matches low_energy_z_correction.f90).
 * Returns the z increment for one integration step.
 * -------------------------------------------------------------------------- */
__device__ __forceinline__ double low_energy_z_correction_dev(
    double pz_val, double step_len, double beta_val, double beta_ref,
    double mc2, double e_tot_ele)
{
    if (mc2 * (beta_ref * pz_val) * (beta_ref * pz_val) < 3e-7 * e_tot_ele) {
        /* Taylor expansion for small pz — avoids precision loss */
        double mr = mc2 / e_tot_ele;
        double b02 = beta_ref * beta_ref;
        double f_tay = b02 * (2.0 * b02 - mr * mr * 0.5);
        return step_len * pz_val * (1.0 - 1.5 * pz_val * b02 + pz_val * pz_val * f_tay) * mr * mr;
    } else {
        return step_len * (beta_val - beta_ref) / beta_ref;
    }
}

/* --------------------------------------------------------------------------
 * drift_body_dev — core drift physics for a single particle
 *
 * Updates position (x, y), longitudinal (z), and time (t) for a drift of
 * length ds.  Returns 0 on success, 1 if particle is lost (pxy2 >= 1).
 * Does NOT update s_pos — callers handle that themselves.
 * -------------------------------------------------------------------------- */
__device__ int drift_body_dev(
    double *x, double *px, double *y, double *py, double *z, double *pz,
    double *beta, double *t, double mc2, double p0c, double ds)
{
    double delta  = *pz;
    double rel_pc = 1.0 + delta;
    double px_rel = *px / rel_pc;
    double py_rel = *py / rel_pc;
    double pxy2   = px_rel * px_rel + py_rel * py_rel;

    if (pxy2 >= 1.0) return 1;
    if (*beta <= 0.0) return 1;

    double ps_rel = sqrt(1.0 - pxy2);

    *x += ds * px_rel / ps_rel;
    *y += ds * py_rel / ps_rel;

    double p_tot = p0c * rel_pc;
    double A = (mc2 * mc2 * (2.0 * delta + delta * delta)) /
               (p_tot * p_tot + mc2 * mc2);
    *z += ds * (sqrt_one_dev(A) + sqrt_one_dev(-pxy2) / ps_rel);

    *t += ds / (*beta * ps_rel * C_LIGHT);

    return 0;
}

/* =========================================================================
 * DRIFT KERNEL
 * Replicates the physics of track_a_drift (forward, include_ref_motion=true)
 * ========================================================================= */
__global__ void drift_kernel(
    double *vx, double *vpx, double *vy, double *vpy, double *vz, double *vpz,
    int *state, double *beta, double *p0c, double *s_pos, double *t_time,
    double mc2, double length, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    if (state[i] != ALIVE_ST) return;

    if (drift_body_dev(&vx[i], &vpx[i], &vy[i], &vpy[i], &vz[i], &vpz[i],
                        &beta[i], &t_time[i], mc2, p0c[i], length)) {
        state[i] = LOST_PZ;
        return;
    }

    /* s update (direction = +1) */
    s_pos[i] += length;
}

/* --------------------------------------------------------------------------
 * Device helper: compute x^p for small non-negative integer p
 * -------------------------------------------------------------------------- */
__device__ __forceinline__ double ipow(double x, int p)
{
    if (p == 0) return 1.0;
    double r = 1.0;
    for (int k = 0; k < p; k++) r *= x;
    return r;
}

/* --------------------------------------------------------------------------
 * Device helper: apply magnetic multipole kicks (ab_multipole_kick equivalent)
 *
 * Precomputed c_multi coefficients are passed in cm[N_MULTI*N_MULTI].
 * Scaled a2[n], b2[n] arrays include charge/orientation/scale factors.
 * -------------------------------------------------------------------------- */
/* Must match Bmad's n_pole_maxx + 1 (defined in bmad_struct.f90) */
#define N_MULTI 22

__device__ void multipole_kick_dev(
    const double *a2, const double *b2, int ix_max,
    const double *cm,
    double x, double y, double *kx_out, double *ky_out)
{
    double kx = 0.0, ky = 0.0;
    for (int nn = 0; nn <= ix_max; nn++) {
        if (a2[nn] == 0.0 && b2[nn] == 0.0) continue;
        /* even m — cm is Fortran column-major: cm(nn, m) at offset m*N_MULTI+nn */
        for (int m = 0; m <= nn; m += 2) {
            double f = cm[m * N_MULTI + nn] * ipow(x, nn - m) * ipow(y, m);
            kx += b2[nn] * f;
            ky -= a2[nn] * f;
        }
        /* odd m */
        for (int m = 1; m <= nn; m += 2) {
            double f = cm[m * N_MULTI + nn] * ipow(x, nn - m) * ipow(y, m);
            kx += a2[nn] * f;
            ky += b2[nn] * f;
        }
    }
    *kx_out = kx;
    *ky_out = ky;
}

/* --------------------------------------------------------------------------
 * apply_electric_kick_dev — apply a scaled electric multipole kick and
 * update pz, beta, and z.  Caller pre-computes the scaled kick (kx_s, ky_s)
 * including any element-specific factors (1/beta, (1+g*x)/ps, etc.).
 * Returns 1 if particle is lost (alpha < -1), 0 on success.
 * -------------------------------------------------------------------------- */
__device__ int apply_electric_kick_dev(
    double kx_s, double ky_s,
    double *px, double *py, double *pz, double *z,
    double *beta_val, double *beta_arr_i,
    double mc2, double p0c_val, int *state_i)
{
    double px_old = *px, py_old = *py, pz_old = *pz;
    *px += kx_s;
    *py += ky_s;
    double alpha = (kx_s * (2.0*px_old + kx_s) + ky_s * (2.0*py_old + ky_s))
                   / ((1.0 + pz_old) * (1.0 + pz_old));
    if (alpha < -1.0) { *state_i = LOST_PZ; return 1; }
    *pz = pz_old + (1.0 + pz_old) * sqrt_one_dev(alpha);
    double new_beta = (1.0 + *pz) / sqrt((1.0 + *pz) * (1.0 + *pz)
                      + (mc2 / p0c_val) * (mc2 / p0c_val));
    *z = *z * new_beta / *beta_val;
    *beta_val = new_beta;
    *beta_arr_i = new_beta;
    return 0;
}

/* --------------------------------------------------------------------------
 * quad_mat2_calc_dev — compute 2x2 transfer matrix and z-correction terms
 *
 * Given focusing strength k_val and step length, computes:
 *   c, s: cosine-like and sine-like matrix elements
 *   zc1, zc2, zc3: z-correction coefficients for the plane
 * -------------------------------------------------------------------------- */
__device__ void quad_mat2_calc_dev(
    double k_val, double step_len, double rel_p,
    double *c_out, double *s_out,
    double *zc1, double *zc2, double *zc3)
{
    double abs_k = fabs(k_val);
    double sqrt_k = sqrt(abs_k);
    double sk_l = sqrt_k * step_len;

    if (fabs(sk_l) < 1e-10) {
        double kl2 = k_val * step_len * step_len;
        *c_out = 1.0 + kl2 * 0.5;
        *s_out = (1.0 + kl2 / 6.0) * step_len;
    } else if (k_val < 0.0) {
        *c_out = cos(sk_l);
        *s_out = sin(sk_l) / sqrt_k;
    } else {
        *c_out = cosh(sk_l);
        *s_out = sinh(sk_l) / sqrt_k;
    }

    double c = *c_out, s = *s_out;
    *zc1 = k_val * (-c * s + step_len) / 4.0;
    *zc2 = -k_val * s * s / (2.0 * rel_p);
    *zc3 = -(c * s + step_len) / (4.0 * rel_p * rel_p);
}

/* =========================================================================
 * QUADRUPOLE KERNEL
 * Replicates the full body of track_a_quadrupole including:
 *   - Split-step integration with n_step steps
 *   - Magnetic and electric multipole kicks (interleaved)
 *   - low_energy_z_correction
 *   - Forward tracking only (direction=1, time_dir=1)
 *
 * When ix_mag_max < 0, ix_elec_max < 0, and n_step == 1, reduces to simple case.
 * ========================================================================= */
__global__ void quad_kernel(
    double *vx, double *vpx, double *vy, double *vpy, double *vz, double *vpz,
    int *state, double *beta_arr, double *p0c_arr, double *t_arr,
    double mc2, double b1, double ele_length,
    double delta_ref_time, double e_tot_ele, double charge_dir,
    int n_particles,
    /* Multipole parameters (may be NULL/unused if ix < 0) */
    const double *d_a2, const double *d_b2, const double *d_cm,
    int ix_mag_max, int n_step,
    /* Electric multipole parameters */
    const double *d_ea2, const double *d_eb2, int ix_elec_max)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_particles) return;
    if (state[i] != ALIVE_ST) return;

    int has_mag = (ix_mag_max >= 0);
    int has_elec = (ix_elec_max >= 0);
    double step_len = ele_length / (double)n_step;
    double z_start = vz[i];
    double t_start = t_arr[i];
    double beta_val = beta_arr[i];
    double p0c_val = p0c_arr[i];
    double beta_ref = p0c_val / e_tot_ele;

    /* Entrance half magnetic multipole kick (scale = r_step/2, built into d_a2/d_b2) */
    if (has_mag) {
        double kx, ky;
        multipole_kick_dev(d_a2, d_b2, ix_mag_max, d_cm,
                           vx[i], vy[i], &kx, &ky);
        vpx[i] += 0.5 * kx;
        vpy[i] += 0.5 * ky;
    }

    /* Entrance half electric multipole kick */
    if (has_elec) {
        double kx, ky;
        multipole_kick_dev(d_ea2, d_eb2, ix_elec_max, d_cm,
                           vx[i], vy[i], &kx, &ky);
        if (apply_electric_kick_dev(0.5 * kx / beta_val, 0.5 * ky / beta_val,
                &vpx[i], &vpy[i], &vpz[i], &vz[i],
                &beta_val, &beta_arr[i], mc2, p0c_val, &state[i])) return;
    }

    /* Body: n_step integration steps */
    for (int istep = 1; istep <= n_step; istep++) {

        double rel_p = 1.0 + vpz[i];
        double k1 = charge_dir * b1 / (ele_length * rel_p);

        /* quad_mat2_calc for x plane (k_x = -k1) and y plane (k_y = +k1) */
        double cx, sx, zc_x1, zc_x2, zc_x3;
        double cy, sy, zc_y1, zc_y2, zc_y3;
        quad_mat2_calc_dev(-k1, step_len, rel_p, &cx, &sx, &zc_x1, &zc_x2, &zc_x3);
        quad_mat2_calc_dev( k1, step_len, rel_p, &cy, &sy, &zc_y1, &zc_y2, &zc_y3);

        /* Save pre-matrix coords for z update */
        double x0 = vx[i], px0 = vpx[i];
        double y0 = vy[i], py0 = vpy[i];

        /* z update from quad focusing */
        vz[i] += zc_x1*x0*x0 + zc_x2*x0*px0 + zc_x3*px0*px0 +
                 zc_y1*y0*y0 + zc_y2*y0*py0 + zc_y3*py0*py0;

        /* Apply 2x2 matrices */
        double k1_x = -k1, k1_y = k1;
        vx[i]  = cx * x0 + (sx / rel_p) * px0;
        vpx[i] = (k1_x * sx * rel_p) * x0 + cx * px0;
        vy[i]  = cy * y0 + (sy / rel_p) * py0;
        vpy[i] = (k1_y * sy * rel_p) * y0 + cy * py0;

        /* Low energy z correction */
        vz[i] += low_energy_z_correction_dev(vpz[i], step_len, beta_val, beta_ref, mc2, e_tot_ele);

        /* Magnetic multipole kick (half at last step, full otherwise) */
        if (has_mag) {
            double kx, ky;
            multipole_kick_dev(d_a2, d_b2, ix_mag_max, d_cm,
                               vx[i], vy[i], &kx, &ky);
            double scl = (istep == n_step) ? 0.5 : 1.0;
            vpx[i] += scl * kx;
            vpy[i] += scl * ky;
        }

        /* Electric multipole kick (half at last step, full otherwise) */
        if (has_elec) {
            double kx, ky;
            multipole_kick_dev(d_ea2, d_eb2, ix_elec_max, d_cm,
                               vx[i], vy[i], &kx, &ky);
            double scl = (istep == n_step) ? 0.5 : 1.0;
            if (apply_electric_kick_dev(scl * kx / beta_val, scl * ky / beta_val,
                    &vpx[i], &vpy[i], &vpz[i], &vz[i],
                    &beta_val, &beta_arr[i], mc2, p0c_val, &state[i])) return;
        }
    }

    /* Time update */
    t_arr[i] = t_start + delta_ref_time + (z_start - vz[i]) / (beta_val * C_LIGHT);
}

/* --------------------------------------------------------------------------
 * upload_particle_data — H→D transfer of core particle arrays
 * -------------------------------------------------------------------------- */
static int upload_particle_data(int n,
    double *h_vx, double *h_vpx, double *h_vy, double *h_vpy,
    double *h_vz, double *h_vpz,
    int *h_state, double *h_beta, double *h_p0c, double *h_t)
{
    size_t db = (size_t)n * sizeof(double);
    size_t ib = (size_t)n * sizeof(int);
    CUDA_CHECK(cudaMemcpy(d_vec[0], h_vx,    db, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vec[1], h_vpx,   db, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vec[2], h_vy,    db, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vec[3], h_vpy,   db, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vec[4], h_vz,    db, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_vec[5], h_vpz,   db, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_state,  h_state,  ib, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_beta,   h_beta,   db, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_p0c,    h_p0c,    db, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_t,      h_t,      db, cudaMemcpyHostToDevice));
    return 0;
}

/* --------------------------------------------------------------------------
 * download_particle_data — D→H transfer of core particle arrays
 *
 * copy_beta/copy_p0c: set to 1 to also download beta/p0c arrays.
 * Drift: both 0.  Quad/bend: copy_beta only when electric multipoles.
 * Lcavity: both 1.
 * -------------------------------------------------------------------------- */
static int download_particle_data(int n,
    double *h_vx, double *h_vpx, double *h_vy, double *h_vpy,
    double *h_vz, double *h_vpz,
    int *h_state, double *h_beta, double *h_p0c, double *h_t,
    int copy_beta, int copy_p0c)
{
    size_t db = (size_t)n * sizeof(double);
    size_t ib = (size_t)n * sizeof(int);
    CUDA_CHECK(cudaMemcpy(h_vx,    d_vec[0], db, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vpx,   d_vec[1], db, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vy,    d_vec[2], db, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vpy,   d_vec[3], db, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vz,    d_vec[4], db, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_vpz,   d_vec[5], db, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_state,  d_state,  ib, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_t,      d_t,      db, cudaMemcpyDeviceToHost));
    if (copy_beta) { CUDA_CHECK(cudaMemcpy(h_beta, d_beta, db, cudaMemcpyDeviceToHost)); }
    if (copy_p0c)  { CUDA_CHECK(cudaMemcpy(h_p0c, d_p0c,  db, cudaMemcpyDeviceToHost)); }
    return 0;
}

/* --------------------------------------------------------------------------
 * upload_multipole_data — H→D transfer of multipole coefficient arrays
 * -------------------------------------------------------------------------- */
static int upload_multipole_data(
    double *h_a2, double *h_b2, double *h_cm,
    double *h_ea2, double *h_eb2,
    int ix_mag_max, int ix_elec_max)
{
    size_t multi_sz = N_MULTI * sizeof(double);
    size_t cm_sz    = N_MULTI * N_MULTI * sizeof(double);

    if (ix_mag_max >= 0 || ix_elec_max >= 0) {
        if (!d_cm) { CUDA_CHECK(cudaMalloc((void**)&d_cm, cm_sz)); }
        CUDA_CHECK(cudaMemcpy(d_cm, h_cm, cm_sz, cudaMemcpyHostToDevice));
    }
    if (ix_mag_max >= 0) {
        if (!d_a2) { CUDA_CHECK(cudaMalloc((void**)&d_a2, multi_sz)); }
        if (!d_b2) { CUDA_CHECK(cudaMalloc((void**)&d_b2, multi_sz)); }
        CUDA_CHECK(cudaMemcpy(d_a2, h_a2, multi_sz, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_b2, h_b2, multi_sz, cudaMemcpyHostToDevice));
    }
    if (ix_elec_max >= 0) {
        if (!d_ea2) { CUDA_CHECK(cudaMalloc((void**)&d_ea2, multi_sz)); }
        if (!d_eb2) { CUDA_CHECK(cudaMalloc((void**)&d_eb2, multi_sz)); }
        CUDA_CHECK(cudaMemcpy(d_ea2, h_ea2, multi_sz, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_eb2, h_eb2, multi_sz, cudaMemcpyHostToDevice));
    }
    return 0;
}

/* =========================================================================
 * HOST WRAPPER: gpu_track_drift
 *
 * Fortran signature (iso_c_binding):
 *   subroutine gpu_track_drift(vec_x, vec_px, vec_y, vec_py, vec_z, vec_pz,
 *              state, beta, p0c, s_pos, t_time, mc2, length, n) bind(C)
 * ========================================================================= */
extern "C" void gpu_track_drift_(
    double *h_vx, double *h_vpx, double *h_vy, double *h_vpy,
    double *h_vz, double *h_vpz,
    int *h_state, double *h_beta, double *h_p0c,
    double *h_s, double *h_t,
    double mc2, double length, int n)
{
    if (ensure_buffers(n) != 0) return;

    size_t db = (size_t)n * sizeof(double);

    /* Host -> Device */
    if (upload_particle_data(n, h_vx, h_vpx, h_vy, h_vpy, h_vz, h_vpz,
                             h_state, h_beta, h_p0c, h_t) != 0) return;
    CUDA_CHECK_VOID(cudaMemcpy(d_s, h_s, db, cudaMemcpyHostToDevice));

    /* Launch kernel */
    int threads = 256;
    int blocks  = (n + threads - 1) / threads;
    drift_kernel<<<blocks, threads>>>(
        d_vec[0], d_vec[1], d_vec[2], d_vec[3], d_vec[4], d_vec[5],
        d_state, d_beta, d_p0c, d_s, d_t, mc2, length, n);
    CUDA_CHECK_VOID(cudaGetLastError());
    CUDA_CHECK_VOID(cudaDeviceSynchronize());

    /* Device -> Host */
    if (download_particle_data(n, h_vx, h_vpx, h_vy, h_vpy, h_vz, h_vpz,
                               h_state, h_beta, h_p0c, h_t, 0, 0) != 0) return;
    CUDA_CHECK_VOID(cudaMemcpy(h_s, d_s, db, cudaMemcpyDeviceToHost));
}

/* =========================================================================
 * HOST WRAPPER: gpu_track_quad
 * ========================================================================= */
extern "C" void gpu_track_quad_(
    double *h_vx, double *h_vpx, double *h_vy, double *h_vpy,
    double *h_vz, double *h_vpz,
    int *h_state, double *h_beta, double *h_p0c, double *h_t,
    double mc2, double b1, double ele_length,
    double delta_ref_time, double e_tot_ele, double charge_dir,
    int n_particles,
    double *h_a2, double *h_b2, double *h_cm,
    int ix_mag_max, int n_step,
    double *h_ea2, double *h_eb2, int ix_elec_max)
{
    if (ensure_buffers(n_particles) != 0) return;

    /* Host -> Device */
    if (upload_particle_data(n_particles, h_vx, h_vpx, h_vy, h_vpy, h_vz, h_vpz,
                             h_state, h_beta, h_p0c, h_t) != 0) return;
    if (upload_multipole_data(h_a2, h_b2, h_cm, h_ea2, h_eb2,
                              ix_mag_max, ix_elec_max) != 0) return;

    /* Launch kernel */
    int threads = 256;
    int blocks  = (n_particles + threads - 1) / threads;
    quad_kernel<<<blocks, threads>>>(
        d_vec[0], d_vec[1], d_vec[2], d_vec[3], d_vec[4], d_vec[5],
        d_state, d_beta, d_p0c, d_t,
        mc2, b1, ele_length, delta_ref_time, e_tot_ele, charge_dir,
        n_particles, d_a2, d_b2, d_cm, ix_mag_max, n_step,
        d_ea2, d_eb2, ix_elec_max);
    CUDA_CHECK_VOID(cudaGetLastError());
    CUDA_CHECK_VOID(cudaDeviceSynchronize());

    /* Device -> Host (electric kicks can modify beta) */
    if (download_particle_data(n_particles, h_vx, h_vpx, h_vy, h_vpy, h_vz, h_vpz,
                               h_state, h_beta, h_p0c, h_t,
                               (ix_elec_max >= 0), 0) != 0) return;
}

/* =========================================================================
 * BEND KERNEL
 * Replicates the body of track_a_bend including:
 *   - General nonlinear bend map (g != 0)
 *   - sbend_body_with_k1_map (b1 != 0, quad component in bend)
 *   - Pure drift fallback (g=0, dg=0)
 *   - Linear approximation for near-axis particles
 *   - Split-step magnetic/electric multipole kicks
 *   - Low energy z correction (for k1 map)
 *   - Forward tracking only (direction=1, time_dir=1)
 * ========================================================================= */

/* Device helper: sinc(x) = sin(x)/x, numerically stable for small x */
__device__ __forceinline__ double sinc_dev(double x)
{
    if (fabs(x) < 1e-8) return 1.0;
    return sin(x) / x;
}

/* Device helper: cosc(x) = (1-cos(x))/x^2, numerically stable */
__device__ __forceinline__ double cosc_dev(double x)
{
    if (fabs(x) < 1e-8) return 0.5;
    double h = x * 0.5;
    double s = sinc_dev(h);
    return 0.5 * s * s;
}

/* Device helper: sincc(x) = (x - sin(x))/x^3, numerically stable */
__device__ __forceinline__ double sincc_dev(double x)
{
    double x2 = x * x;
    if (fabs(x) < 0.1) {
        return (1.0/6.0) + x2 * ((-1.0/120.0) + x2 * ((1.0/5040.0) + x2 * (-1.0/362880.0)));
    }
    return (x - sin(x)) / (x * x2);
}

__global__ void bend_kernel(
    double *vx, double *vpx, double *vy, double *vpy, double *vz, double *vpz,
    int *state, double *beta_arr, double *p0c_arr, double *t_arr,
    double mc2, double g, double g_tot, double dg, double b1,
    double ele_length, double delta_ref_time, double e_tot_ele,
    double rel_charge_dir,
    double p0c_ele,
    int n_particles,
    /* Multipole parameters */
    const double *d_a2, const double *d_b2, const double *d_cm,
    int ix_mag_max, int n_step,
    /* Electric multipole parameters */
    const double *d_ea2, const double *d_eb2, int ix_elec_max)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_particles) return;
    if (state[i] != ALIVE_ST) return;

    int has_mag = (ix_mag_max >= 0);
    int has_elec = (ix_elec_max >= 0);
    double step_len = ele_length / (double)n_step;
    double angle = g * step_len;
    double z_start = vz[i];
    double t_start = t_arr[i];
    double beta_val = beta_arr[i];
    double p0c_val = p0c_arr[i];
    double beta_ref = p0c_ele / e_tot_ele;

    /* Entrance half magnetic multipole kick */
    if (has_mag) {
        double kx, ky;
        multipole_kick_dev(d_a2, d_b2, ix_mag_max, d_cm,
                           vx[i], vy[i], &kx, &ky);
        /* For bends, magnetic kick includes (1+g*x) factor and uses field formulation.
         * The precomputed a2/b2 arrays include all element-level scaling.
         * The (1+g*x) factor is applied here per-particle. */
        double f_gx = 1.0 + g * vx[i];
        vpx[i] += 0.5 * kx * f_gx;
        vpy[i] += 0.5 * ky * f_gx;
    }

    /* Entrance half electric multipole kick */
    if (has_elec) {
        double kx, ky;
        multipole_kick_dev(d_ea2, d_eb2, ix_elec_max, d_cm,
                           vx[i], vy[i], &kx, &ky);
        double f_gx = 1.0 + g * vx[i];
        double rel_p0 = 1.0 + vpz[i];
        double ps = sqrt(rel_p0*rel_p0 - vpx[i]*vpx[i] - vpy[i]*vpy[i]) / rel_p0;
        double f_scale = 0.5 * f_gx / (ps * beta_val);
        if (apply_electric_kick_dev(kx * f_scale, ky * f_scale,
                &vpx[i], &vpy[i], &vpz[i], &vz[i],
                &beta_val, &beta_arr[i], mc2, p0c_val, &state[i])) return;
    }

    /* Body: n_step integration steps */
    for (int istep = 1; istep <= n_step; istep++) {

        double pz = vpz[i];
        double rel_p = 1.0 + pz;

        /* ---- Branch 1: b1 != 0 → sbend_body_with_k1_map ---- */
        if (b1 != 0.0) {
            double k1 = b1 / ele_length;
            double k_x = k1 + g * g_tot;
            double x_c = (g * rel_p - g_tot) / k_x;
            double om_x = sqrt(fabs(k_x) / rel_p);
            double om_y = sqrt(fabs(k1) / rel_p);
            double tau_x = (k_x > 0) ? -1.0 : 1.0;
            double tau_y = (k1 > 0) ? 1.0 : -1.0;
            /* Note: k_x > 0 means focusing in x → sin/cos; k_x < 0 → sinh/cosh */
            /* But tau_x = -sign(1, k_x), so tau_x < 0 when k_x > 0 */

            double arg_x = om_x * step_len;
            double s_x, c_x, z2;
            if (arg_x < 1e-6) {
                s_x = (1.0 + tau_x * arg_x * arg_x / 6.0) * step_len;
                c_x = 1.0 + tau_x * arg_x * arg_x / 2.0;
                z2 = g * step_len * step_len / (2.0 * rel_p);
            } else if (k_x > 0) {
                s_x = sin(arg_x) / om_x;
                c_x = cos(arg_x);
                z2 = tau_x * g * (1.0 - c_x) / (rel_p * om_x * om_x);
            } else {
                s_x = sinh(arg_x) / om_x;
                c_x = cosh(arg_x);
                z2 = tau_x * g * (1.0 - c_x) / (rel_p * om_x * om_x);
            }

            double arg_y = om_y * step_len;
            double s_y, c_y;
            if (arg_y < 1e-6) {
                s_y = (1.0 + tau_y * arg_y * arg_y / 6.0) * step_len;
                c_y = 1.0 + tau_y * arg_y * arg_y / 2.0;
            } else if (k1 < 0) {
                s_y = sin(om_y * step_len) / om_y;
                c_y = cos(om_y * step_len);
            } else {
                s_y = sinh(om_y * step_len) / om_y;
                c_y = cosh(om_y * step_len);
            }

            double r1 = vx[i] - x_c;
            double r2 = vpx[i];
            double r3 = vy[i];
            double r4 = vpy[i];

            double z0  = -g * x_c * step_len;
            double z1  = -g * s_x;
            double z11 = tau_x * om_x*om_x * (step_len - c_x*s_x) / 4.0;
            double z12 = -tau_x * om_x*om_x * s_x*s_x / (2.0 * rel_p);
            double z22 = -(step_len + c_x*s_x) / (4.0 * rel_p * rel_p);
            double z33 = tau_y * om_y*om_y * (step_len - c_y*s_y) / 4.0;
            double z34 = -tau_y * om_y*om_y * s_y*s_y / (2.0 * rel_p);
            double z44 = -(step_len + c_y*s_y) / (4.0 * rel_p * rel_p);

            vx[i]  = c_x * r1 + s_x * r2 / rel_p + x_c;
            vpx[i] = tau_x * om_x*om_x * rel_p * s_x * r1 + c_x * r2;
            vy[i]  = c_y * r3 + s_y * r4 / rel_p;
            vpy[i] = tau_y * om_y*om_y * rel_p * s_y * r3 + c_y * r4;
            /* orientation*direction = 1 for forward tracking */
            vz[i] += z0 + z1*r1 + z2*r2 +
                     z11*r1*r1 + z12*r1*r2 + z22*r2*r2 +
                     z33*r3*r3 + z34*r3*r4 + z44*r4*r4;

            /* Low energy z correction for k1 map */
            vz[i] += low_energy_z_correction_dev(vpz[i], step_len, beta_val, beta_ref, mc2, e_tot_ele);

        /* ---- Branch 2: g=0 and dg=0 → pure drift ---- */
        } else if ((g == 0.0 && dg == 0.0) || step_len == 0.0) {
            double t_dummy = t_arr[i];
            if (drift_body_dev(&vx[i], &vpx[i], &vy[i], &vpy[i], &vz[i], &vpz[i],
                                &beta_val, &t_dummy, mc2, p0c_val, step_len)) {
                state[i] = LOST_PZ; return;
            }
            /* time handled at end (not from drift_body_dev) */

        /* ---- Branch 3: General bend (g != 0, b1 = 0) ---- */
        } else {
            double x  = vx[i];
            double px = vpx[i];
            double y  = vy[i];
            double py = vpy[i];
            double z  = vz[i];
            double rel_p2 = rel_p * rel_p;

            /* Linear approximation for near-axis particles */
            if (dg == 0.0 && fabs(x*g) < 1e-9 && fabs(px) < 1e-9 &&
                fabs(py) < 1e-9 && fabs(pz) < 1e-9) {
                double ll = step_len;
                double cos_a = cos(angle);
                double sin_a = sin(angle);
                double sinc_a = sinc_dev(angle);
                double cosc_a = cosc_dev(angle);
                double gam2 = mc2*mc2 / (rel_p2 * p0c_val*p0c_val + mc2*mc2);
                double m56 = ll * (gam2 - (g*ll)*(g*ll) * sincc_dev(angle));

                vx[i]  = cos_a * x      + ll*sinc_a * px + g*ll*ll*cosc_a * pz;
                vpx[i] = -g*sin_a * x   + cos_a * px     + g*ll*sinc_a * pz;
                vy[i]  = y + ll * py;
                vz[i]  = -g*ll*sinc_a*x - g*ll*ll*cosc_a*px + z + m56*pz;

            /* General nonlinear case */
            } else {
                double sinc_a = sinc_dev(angle);
                double pt = sqrt(rel_p2 - py*py);
                if (fabs(px) > pt) { state[i] = LOST_PZ; return; }
                double g_p = g_tot / pt;
                double phi_1 = asin(px / pt);
                double cos_a = cos(angle);
                double sin_a = sin(angle);
                double cosc_a = cosc_dev(angle);
                double cos_plus = cos(angle + phi_1);
                double sin_plus = sin(angle + phi_1);
                double alpha_b = 2.0*(1.0+g*x)*sin_plus*step_len*sinc_a -
                                 g_p*((1.0+g*x)*step_len*sinc_a)*((1.0+g*x)*step_len*sinc_a);
                double r_val = cos_plus*cos_plus + g_p*alpha_b;

                if (r_val < 0.0 || (fabs(g_p) < 1e-5 && fabs(cos_plus) < 1e-5)) {
                    state[i] = LOST_PZ;
                    return;
                }

                double rad = sqrt(r_val);
                double xi;
                if (cos_plus > 0.0) {
                    double denom = rad + cos_plus;
                    xi = alpha_b / denom;
                } else {
                    if (fabs(g_p) < 1e-30) { state[i] = LOST_PZ; return; }
                    xi = (rad - cos_plus) / g_p;
                }
                vx[i] = x*cos_a - step_len*step_len*g*cosc_a + xi;

                /* Check aperture limit */
                if (fabs(vx[i]) > 1.0) {
                    state[i] = LOST_PZ;
                    return;
                }

                double L_u = xi;
                double L_v = -(step_len*sinc_a + x*sin_a);  /* time_dir=1 */
                double L_c = sqrt(L_v*L_v + L_u*L_u);
                double angle_p = 2.0*(angle + phi_1 - atan2(L_u, -L_v));  /* time_dir=1 */
                double L_p = L_c / sinc_dev(angle_p*0.5);  /* time_dir=1 */
                vpx[i] = pt * sin(phi_1 + angle - angle_p);
                vy[i]  = y + py * L_p / pt;
                vz[i]  = z + beta_val * step_len / beta_ref - rel_p * L_p / pt;
            }
        }

        /* Multipole kick after each step */
        if (has_mag) {
            double kx, ky;
            multipole_kick_dev(d_a2, d_b2, ix_mag_max, d_cm,
                               vx[i], vy[i], &kx, &ky);
            double scl = (istep == n_step) ? 0.5 : 1.0;
            double f_gx = 1.0 + g * vx[i];
            vpx[i] += scl * kx * f_gx;
            vpy[i] += scl * ky * f_gx;
        }

        if (has_elec) {
            double kx, ky;
            multipole_kick_dev(d_ea2, d_eb2, ix_elec_max, d_cm,
                               vx[i], vy[i], &kx, &ky);
            double f_gx = 1.0 + g * vx[i];
            double rel_p0 = 1.0 + vpz[i];
            double ps = sqrt(rel_p0*rel_p0 - vpx[i]*vpx[i] - vpy[i]*vpy[i]) / rel_p0;
            double scl = (istep == n_step) ? 0.5 : 1.0;
            double f_scale = scl * f_gx / (ps * beta_val);
            if (apply_electric_kick_dev(kx * f_scale, ky * f_scale,
                    &vpx[i], &vpy[i], &vpz[i], &vz[i],
                    &beta_val, &beta_arr[i], mc2, p0c_val, &state[i])) return;
        }
    }

    /* Time update */
    t_arr[i] = t_start + delta_ref_time + (z_start - vz[i]) / (beta_val * C_LIGHT);
}

/* =========================================================================
 * HOST WRAPPER: gpu_track_bend
 * ========================================================================= */
extern "C" void gpu_track_bend_(
    double *h_vx, double *h_vpx, double *h_vy, double *h_vpy,
    double *h_vz, double *h_vpz,
    int *h_state, double *h_beta, double *h_p0c, double *h_t,
    double mc2, double g, double g_tot, double dg, double b1,
    double ele_length, double delta_ref_time, double e_tot_ele,
    double rel_charge_dir,
    double p0c_ele,
    int n_particles,
    double *h_a2, double *h_b2, double *h_cm,
    int ix_mag_max, int n_step,
    double *h_ea2, double *h_eb2, int ix_elec_max)
{
    if (ensure_buffers(n_particles) != 0) return;

    /* Host -> Device */
    if (upload_particle_data(n_particles, h_vx, h_vpx, h_vy, h_vpy, h_vz, h_vpz,
                             h_state, h_beta, h_p0c, h_t) != 0) return;
    if (upload_multipole_data(h_a2, h_b2, h_cm, h_ea2, h_eb2,
                              ix_mag_max, ix_elec_max) != 0) return;

    /* Launch kernel */
    int threads = 256;
    int blocks  = (n_particles + threads - 1) / threads;
    bend_kernel<<<blocks, threads>>>(
        d_vec[0], d_vec[1], d_vec[2], d_vec[3], d_vec[4], d_vec[5],
        d_state, d_beta, d_p0c, d_t,
        mc2, g, g_tot, dg, b1, ele_length, delta_ref_time, e_tot_ele,
        rel_charge_dir, p0c_ele,
        n_particles, d_a2, d_b2, d_cm, ix_mag_max, n_step,
        d_ea2, d_eb2, ix_elec_max);
    CUDA_CHECK_VOID(cudaGetLastError());
    CUDA_CHECK_VOID(cudaDeviceSynchronize());

    /* Device -> Host (electric kicks can modify beta) */
    if (download_particle_data(n_particles, h_vx, h_vpx, h_vy, h_vpy, h_vz, h_vpz,
                               h_state, h_beta, h_p0c, h_t,
                               (ix_elec_max >= 0), 0) != 0) return;
}

/* =========================================================================
 * LCAVITY KERNEL
 * Replicates track_a_lcavity stair-step RF approximation.
 * Each particle independently loops through drift + energy kick steps.
 * Handles ponderomotive transverse kicks for standing-wave cavities.
 * Forward tracking only, relative time tracking.
 * ========================================================================= */

#define TWOPI 6.283185307179586476925286766559
#define STANDING_WAVE 1
#define TRAVELING_WAVE 2

/* Device helper: dpc_given_dE — momentum change from energy change.
 * Uses rationalization to avoid catastrophic cancellation for small dE. */
__device__ __forceinline__ double dpc_given_dE_dev(double pc_old, double mc2, double dE)
{
    double del2 = dE * dE + 2.0 * sqrt(pc_old * pc_old + mc2 * mc2) * dE;
    double pc_new = sqrt(pc_old * pc_old + del2);
    /* Rationalize: pc_new - pc_old = del2 / (pc_new + pc_old) to avoid cancellation */
    return del2 / (pc_new + pc_old);
}

/* Device helper: lcavity fringe kick (entrance or exit)
 * Matches track_a_lcavity.f90 fringe_kick subroutine (lines 232-301)
 * for GPU case: body_dir=1, time_dir=1.
 * edge: +1 for entrance, -1 for exit.
 */
__device__ void lcavity_fringe_kick_dev(
    double &x, double &px, double &y, double &py, double &z, double &pz,
    double &beta_val, double p0c, double mc2,
    int edge, double gradient_tot, double charge_ratio,
    double rf_frequency, double phi0_total)
{
    double particle_time = -z / (beta_val * C_LIGHT);
    double phase = TWOPI * (phi0_total + particle_time * rf_frequency);

    double ez_field = gradient_tot * cos(phase);
    double rf_omega = TWOPI * rf_frequency / C_LIGHT;
    double dez_dz_field = gradient_tot * sin(phase) * rf_omega;

    double ff = edge * charge_ratio;
    double f = ff / p0c;

    double dE = -ff * 0.5 * dez_dz_field * (x * x + y * y);

    double pc = p0c * (1.0 + pz);
    double pz_end = pz + dpc_given_dE_dev(pc, mc2, dE) / p0c;

    /* to_energy_coords */
    z = z / beta_val;
    pz = (1.0 + pz) / beta_val;

    /* kicks in energy coords */
    px = px - f * ez_field * x;
    py = py - f * ez_field * y;
    pz = pz + dE / p0c;

    /* to_momentum_coords */
    double pc_new = (1.0 + pz_end) * p0c;
    double beta_new = pc_new / (p0c * pz);
    z = z * beta_new;
    pz = pz_end;
    beta_val = beta_new;
}

/* --------------------------------------------------------------------------
 * ponderomotive_kick_dev — standing-wave ponderomotive transverse kick
 *
 * Applied symmetrically before and after each RF energy kick step.
 * -------------------------------------------------------------------------- */
__device__ void ponderomotive_kick_dev(
    double *px, double *py, double *z, double *t,
    double x, double y, double pz, double beta_val,
    double grad, double l_active, double p0c, int n_rf_steps)
{
    double rel_p = 1.0 + pz;
    double coef = grad * grad * l_active / (16.0 * p0c * p0c * rel_p * n_rf_steps);
    *px -= coef * x;
    *py -= coef * y;
    double dzp = -0.5 * coef * (x * x + y * y) / rel_p;
    *z += dzp;
    *t -= dzp / (C_LIGHT * beta_val);
}

__global__ void lcavity_kernel(
    double *vx, double *vpx, double *vy, double *vpy, double *vz, double *vpz,
    int *state, double *beta_arr, double *p0c_arr, double *t_arr,
    double mc2,
    /* Step data arrays: n_steps_total entries (indices 0..n_rf_steps+1) */
    const double *step_s0, const double *step_s,
    const double *step_p0c, const double *step_p1c,
    const double *step_scale, const double *step_time,
    int n_rf_steps,
    /* Element parameters */
    double voltage, double voltage_err, double field_autoscale,
    double rf_frequency,
    double phi0_total,      /* phi0 + phi0_err + phi0_multipass (precomputed) */
    double voltage_tot, double l_active,
    int cavity_type,        /* 1=standing_wave, 2=traveling_wave */
    int fringe_at,          /* 0=none, 1=entrance, 2=exit, 3=both */
    double charge_ratio,    /* charge_of(species) / (2 * charge_of(ref_species)) */
    int n_particles)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_particles) return;
    if (state[i] != ALIVE_ST) return;

    double x = vx[i], px = vpx[i], y = vy[i], py = vpy[i], z = vz[i], pz = vpz[i];
    double beta_val = beta_arr[i], p0c = p0c_arr[i], t = t_arr[i];

    double s_now = step_s0[0];
    int ix_step_end = n_rf_steps + 1;  /* phantom step index */

    /* Ponderomotive kicks only for standing wave + forward tracking */
    int do_ponderomotive = (cavity_type == STANDING_WAVE) && (l_active > 0.0);

    for (int ix = 0; ix <= ix_step_end; ix++) {
        /* ---- Drift to step boundary ---- */
        double ds = step_s[ix] - s_now;
        s_now = step_s[ix];

        if (ds != 0.0) {
            if (drift_body_dev(&x, &px, &y, &py, &z, &pz,
                                &beta_val, &t, mc2, p0c, ds)) {
                state[i] = LOST_PZ; return;
            }
        }

        /* ---- Entrance fringe: step 0, after drift, before energy kick ---- */
        if (ix == 0 && (fringe_at & 1) && l_active > 0.0) {
            double grad_tot = voltage_tot * field_autoscale / l_active;
            lcavity_fringe_kick_dev(x, px, y, py, z, pz, beta_val, p0c, mc2,
                +1, grad_tot, charge_ratio, rf_frequency, phi0_total);
        }

        /* ---- Stair-step kick (not at phantom step) ---- */
        if (ix != ix_step_end) {
            /* Upstream ponderomotive kick (skip at step 0) */
            if (do_ponderomotive && ix > 0) {
                double grad = field_autoscale * voltage_tot / l_active;
                ponderomotive_kick_dev(&px, &py, &z, &t, x, y, pz, beta_val,
                                       grad, l_active, p0c, n_rf_steps);
            }

            /* ---- Energy kick ---- */
            /* Compute RF phase (relative time tracking) */
            double particle_time = -z / (beta_val * C_LIGHT);
            double phase = TWOPI * (phi0_total + particle_time * rf_frequency);

            /* Energy change */
            double dE_amp = (voltage + voltage_err) * step_scale[ix] * field_autoscale;
            double dE = dE_amp * cos(phase);

            /* Compute pz_end before coordinate conversion (avoid round-off) */
            double rel_p = 1.0 + pz;
            double pc = rel_p * p0c;
            double pz_end = pz + dpc_given_dE_dev(pc, mc2, dE) / p0c;

            /* to_energy_coords: (z, pz) -> (c*(t0-t), E/p0c) */
            z = z / beta_val;
            pz = (1.0 + pz) / beta_val;  /* now pz = E/p0c in energy coords */

            /* Apply energy kick in energy coords */
            pz = pz + dE / p0c;

            /* to_momentum_coords: convert back */
            double pc_new = (1.0 + pz_end) * p0c;
            double beta_new = pc_new / (p0c * pz);  /* pz here is E/p0c */
            z = z * beta_new;
            pz = pz_end;
            beta_val = beta_new;

            /* orbit_reference_energy_correction */
            double p1c = step_p1c[ix];
            double p_rel = p0c / p1c;
            px = px * p_rel;
            py = py * p_rel;
            pz = (pz * p0c - (p1c - p0c)) / p1c;
            p0c = p1c;

            /* Update beta after reference energy change */
            pc_new = (1.0 + pz) * p0c;
            beta_val = pc_new / sqrt(pc_new * pc_new + mc2 * mc2);

            /* Downstream ponderomotive kick (skip at step n_rf_steps) */
            if (do_ponderomotive && ix < n_rf_steps) {
                double grad = field_autoscale * voltage_tot / l_active;
                ponderomotive_kick_dev(&px, &py, &z, &t, x, y, pz, beta_val,
                                       grad, l_active, p0c, n_rf_steps);
            }
        }

        /* ---- Exit fringe: step n_rf_steps, after energy kick ---- */
        if (ix == n_rf_steps && (fringe_at & 2) && l_active > 0.0) {
            double grad_tot = voltage_tot * field_autoscale / l_active;
            lcavity_fringe_kick_dev(x, px, y, py, z, pz, beta_val, p0c, mc2,
                -1, grad_tot, charge_ratio, rf_frequency, phi0_total);
        }
    }

    /* Write back */
    vx[i] = x;   vpx[i] = px;  vy[i] = y;   vpy[i] = py;
    vz[i] = z;   vpz[i] = pz;
    beta_arr[i] = beta_val;
    p0c_arr[i]  = p0c;
    t_arr[i]    = t;
}

/* =========================================================================
 * HOST WRAPPER: gpu_track_lcavity
 * ========================================================================= */

/* Step data device buffers */
static double *d_step_s0   = NULL;
static double *d_step_s    = NULL;
static double *d_step_p0c  = NULL;
static double *d_step_p1c  = NULL;
static double *d_step_scl  = NULL;
static double *d_step_time = NULL;
static int     d_step_cap  = 0;

static int ensure_step_buffers(int n_steps_total)
{
    if (n_steps_total <= d_step_cap) return 0;
    if (d_step_s0)   cudaFree(d_step_s0);
    if (d_step_s)    cudaFree(d_step_s);
    if (d_step_p0c)  cudaFree(d_step_p0c);
    if (d_step_p1c)  cudaFree(d_step_p1c);
    if (d_step_scl)  cudaFree(d_step_scl);
    if (d_step_time) cudaFree(d_step_time);
    d_step_s0 = d_step_s = d_step_p0c = d_step_p1c = d_step_scl = d_step_time = NULL;

    size_t sz = (size_t)n_steps_total * sizeof(double);
    if (cudaMalloc((void**)&d_step_s0,   sz) != cudaSuccess) goto sfail;
    if (cudaMalloc((void**)&d_step_s,    sz) != cudaSuccess) goto sfail;
    if (cudaMalloc((void**)&d_step_p0c,  sz) != cudaSuccess) goto sfail;
    if (cudaMalloc((void**)&d_step_p1c,  sz) != cudaSuccess) goto sfail;
    if (cudaMalloc((void**)&d_step_scl,  sz) != cudaSuccess) goto sfail;
    if (cudaMalloc((void**)&d_step_time, sz) != cudaSuccess) goto sfail;
    d_step_cap = n_steps_total;
    return 0;
sfail:
    fprintf(stderr, "[gpu_tracking] cudaMalloc failed for %d step buffers\n", n_steps_total);
    /* Clean up any partially allocated step buffers */
    if (d_step_s0)   cudaFree(d_step_s0);   d_step_s0   = NULL;
    if (d_step_s)    cudaFree(d_step_s);    d_step_s    = NULL;
    if (d_step_p0c)  cudaFree(d_step_p0c);  d_step_p0c  = NULL;
    if (d_step_p1c)  cudaFree(d_step_p1c);  d_step_p1c  = NULL;
    if (d_step_scl)  cudaFree(d_step_scl);  d_step_scl  = NULL;
    if (d_step_time) cudaFree(d_step_time); d_step_time = NULL;
    d_step_cap = 0;
    return -1;
}

extern "C" void gpu_track_lcavity_(
    double *h_vx, double *h_vpx, double *h_vy, double *h_vpy,
    double *h_vz, double *h_vpz,
    int *h_state, double *h_beta, double *h_p0c, double *h_t,
    double mc2,
    double *h_step_s0, double *h_step_s,
    double *h_step_p0c, double *h_step_p1c,
    double *h_step_scale, double *h_step_time,
    int n_rf_steps,
    double voltage, double voltage_err, double field_autoscale,
    double rf_frequency, double phi0_total,
    double voltage_tot, double l_active,
    int cavity_type,
    int fringe_at, double charge_ratio,
    int n_particles)
{
    if (ensure_buffers(n_particles) != 0) return;

    int n_steps_total = n_rf_steps + 2;  /* indices 0..n_rf_steps+1 */
    if (ensure_step_buffers(n_steps_total) != 0) return;

    size_t sb = (size_t)n_steps_total * sizeof(double);

    /* Host -> Device: particle data */
    if (upload_particle_data(n_particles, h_vx, h_vpx, h_vy, h_vpy, h_vz, h_vpz,
                             h_state, h_beta, h_p0c, h_t) != 0) return;

    /* Host -> Device: step data */
    CUDA_CHECK_VOID(cudaMemcpy(d_step_s0,   h_step_s0,    sb, cudaMemcpyHostToDevice));
    CUDA_CHECK_VOID(cudaMemcpy(d_step_s,    h_step_s,     sb, cudaMemcpyHostToDevice));
    CUDA_CHECK_VOID(cudaMemcpy(d_step_p0c,  h_step_p0c,   sb, cudaMemcpyHostToDevice));
    CUDA_CHECK_VOID(cudaMemcpy(d_step_p1c,  h_step_p1c,   sb, cudaMemcpyHostToDevice));
    CUDA_CHECK_VOID(cudaMemcpy(d_step_scl,  h_step_scale, sb, cudaMemcpyHostToDevice));
    CUDA_CHECK_VOID(cudaMemcpy(d_step_time, h_step_time,  sb, cudaMemcpyHostToDevice));

    /* Launch kernel */
    int threads = 256;
    int blocks  = (n_particles + threads - 1) / threads;
    lcavity_kernel<<<blocks, threads>>>(
        d_vec[0], d_vec[1], d_vec[2], d_vec[3], d_vec[4], d_vec[5],
        d_state, d_beta, d_p0c, d_t, mc2,
        d_step_s0, d_step_s, d_step_p0c, d_step_p1c, d_step_scl, d_step_time,
        n_rf_steps,
        voltage, voltage_err, field_autoscale,
        rf_frequency, phi0_total, voltage_tot, l_active,
        cavity_type,
        fringe_at, charge_ratio,
        n_particles);
    CUDA_CHECK_VOID(cudaGetLastError());
    CUDA_CHECK_VOID(cudaDeviceSynchronize());

    /* Device -> Host (lcavity changes beta and p0c) */
    if (download_particle_data(n_particles, h_vx, h_vpx, h_vy, h_vpy, h_vz, h_vpz,
                               h_state, h_beta, h_p0c, h_t, 1, 1) != 0) return;
}

/* --------------------------------------------------------------------------
 * gpu_tracking_cleanup — release cached device buffers
 * -------------------------------------------------------------------------- */
extern "C" void gpu_tracking_cleanup_(void)
{
    for (int k = 0; k < 6; k++) { if (d_vec[k]) cudaFree(d_vec[k]); d_vec[k] = NULL; }
    if (d_state) cudaFree(d_state); d_state = NULL;
    if (d_beta)  cudaFree(d_beta);  d_beta  = NULL;
    if (d_p0c)   cudaFree(d_p0c);   d_p0c   = NULL;
    if (d_s)     cudaFree(d_s);     d_s     = NULL;
    if (d_t)     cudaFree(d_t);     d_t     = NULL;
    if (d_a2)    cudaFree(d_a2);    d_a2    = NULL;
    if (d_b2)    cudaFree(d_b2);    d_b2    = NULL;
    if (d_ea2)   cudaFree(d_ea2);   d_ea2   = NULL;
    if (d_eb2)   cudaFree(d_eb2);   d_eb2   = NULL;
    if (d_cm)    cudaFree(d_cm);    d_cm    = NULL;
    if (d_step_s0)   cudaFree(d_step_s0);   d_step_s0   = NULL;
    if (d_step_s)    cudaFree(d_step_s);    d_step_s    = NULL;
    if (d_step_p0c)  cudaFree(d_step_p0c);  d_step_p0c  = NULL;
    if (d_step_p1c)  cudaFree(d_step_p1c);  d_step_p1c  = NULL;
    if (d_step_scl)  cudaFree(d_step_scl);  d_step_scl  = NULL;
    if (d_step_time) cudaFree(d_step_time); d_step_time = NULL;
    d_step_cap = 0;
    d_cap = 0;
}

#endif /* USE_GPU_TRACKING */
