/*
 * gpu_tracking_kernels.cu
 *
 * CUDA kernels for GPU-accelerated particle tracking through drift and
 * quadrupole elements.  Called from Fortran via iso_c_binding wrappers
 * defined in gpu_tracking_mod.f90.
 *
 * Build requirements:
 *   - CUDA Toolkit (cuda_runtime.h)
 *   - Compile with nvcc: -DUSE_GPU_TRACKING
 *   - Link with: -lcudart
 *
 * The wrapper caches device memory allocations to avoid re-allocation
 * on repeated calls (similar pattern to cufft_wrapper.c).
 */

#ifdef USE_GPU_TRACKING

#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

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

/* --------------------------------------------------------------------------
 * ensure_buffers — (re-)allocate device arrays when size changes
 * -------------------------------------------------------------------------- */
static int ensure_buffers(int n)
{
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

    double delta  = vpz[i];
    double rel_pc = 1.0 + delta;
    double px_rel = vpx[i] / rel_pc;
    double py_rel = vpy[i] / rel_pc;
    double pxy2   = px_rel * px_rel + py_rel * py_rel;

    if (pxy2 >= 1.0) {
        state[i] = LOST_PZ;
        return;
    }

    double ps_rel = sqrt(1.0 - pxy2);

    /* Position update */
    vx[i] += length * px_rel / ps_rel;
    vy[i] += length * py_rel / ps_rel;

    /* z update: dz = length * (sqrt_one(A) + sqrt_one(-pxy2)/ps_rel) */
    double p_tot = p0c[i] * rel_pc;
    double A = (mc2 * mc2 * (2.0 * delta + delta * delta)) /
               (p_tot * p_tot + mc2 * mc2);
    double dz = length * (sqrt_one_dev(A) + sqrt_one_dev(-pxy2) / ps_rel);
    vz[i] += dz;

    /* Time update */
    double dt = length / (beta[i] * ps_rel * C_LIGHT);
    t_time[i] += dt;

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
#define N_MULTI 22   /* n_pole_maxx + 1 = 22 */

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

    /* Entrance half electric multipole kick (scale = step_len/2, built into d_ea2/d_eb2) */
    /* d_ea2/d_eb2 are precomputed WITHOUT 1/beta; we apply 1/beta per-particle here */
    if (has_elec) {
        double kx, ky;
        multipole_kick_dev(d_ea2, d_eb2, ix_elec_max, d_cm,
                           vx[i], vy[i], &kx, &ky);
        kx /= beta_val;
        ky /= beta_val;
        double px_old = vpx[i], py_old = vpy[i], pz_old = vpz[i];
        double kx_h = 0.5 * kx, ky_h = 0.5 * ky;
        vpx[i] += kx_h;
        vpy[i] += ky_h;
        /* pz update for electric kick */
        double alpha = (kx_h * (2.0*px_old + kx_h) + ky_h * (2.0*py_old + ky_h)) / ((1.0 + pz_old) * (1.0 + pz_old));
        if (alpha < -1.0) { state[i] = LOST_PZ; return; }
        double dpz = (1.0 + pz_old) * sqrt_one_dev(alpha);
        vpz[i] = pz_old + dpz;
        double new_beta = (1.0 + vpz[i]) / sqrt((1.0 + vpz[i])*(1.0 + vpz[i]) + (mc2/p0c_val)*(mc2/p0c_val));
        vz[i] = vz[i] * new_beta / beta_val;
        beta_val = new_beta;
        beta_arr[i] = beta_val;
    }

    /* Body: n_step integration steps */
    for (int istep = 1; istep <= n_step; istep++) {

        double rel_p = 1.0 + vpz[i];
        double k1 = charge_dir * b1 / (ele_length * rel_p);

        /* quad_mat2_calc for x plane (k_x = -k1) */
        double k1_x = -k1;
        double abs_k_x = fabs(k1_x);
        double sqrt_k = sqrt(abs_k_x);
        double sk_l = sqrt_k * step_len;
        double cx, sx;

        if (fabs(sk_l) < 1e-10) {
            double kl2 = k1_x * step_len * step_len;
            cx = 1.0 + kl2 * 0.5;
            sx = (1.0 + kl2 / 6.0) * step_len;
        } else if (k1_x < 0.0) {
            cx = cos(sk_l);
            sx = sin(sk_l) / sqrt_k;
        } else {
            cx = cosh(sk_l);
            sx = sinh(sk_l) / sqrt_k;
        }

        double zc_x1 = k1_x * (-cx * sx + step_len) / 4.0;
        double zc_x2 = -k1_x * sx * sx / (2.0 * rel_p);
        double zc_x3 = -(cx * sx + step_len) / (4.0 * rel_p * rel_p);

        /* quad_mat2_calc for y plane (k_y = +k1) */
        double k1_y = k1;
        double abs_k_y = fabs(k1_y);
        double sqrt_ky = sqrt(abs_k_y);
        double sk_ly = sqrt_ky * step_len;
        double cy, sy;

        if (fabs(sk_ly) < 1e-10) {
            double kl2y = k1_y * step_len * step_len;
            cy = 1.0 + kl2y * 0.5;
            sy = (1.0 + kl2y / 6.0) * step_len;
        } else if (k1_y < 0.0) {
            cy = cos(sk_ly);
            sy = sin(sk_ly) / sqrt_ky;
        } else {
            cy = cosh(sk_ly);
            sy = sinh(sk_ly) / sqrt_ky;
        }

        double zc_y1 = k1_y * (-cy * sy + step_len) / 4.0;
        double zc_y2 = -k1_y * sy * sy / (2.0 * rel_p);
        double zc_y3 = -(cy * sy + step_len) / (4.0 * rel_p * rel_p);

        /* Save pre-matrix coords for z update */
        double x0 = vx[i], px0 = vpx[i];
        double y0 = vy[i], py0 = vpy[i];

        /* z update from quad focusing */
        vz[i] += zc_x1*x0*x0 + zc_x2*x0*px0 + zc_x3*px0*px0 +
                 zc_y1*y0*y0 + zc_y2*y0*py0 + zc_y3*py0*py0;

        /* Apply 2x2 matrices */
        vx[i]  = cx * x0 + (sx / rel_p) * px0;
        vpx[i] = (k1_x * sx * rel_p) * x0 + cx * px0;
        vy[i]  = cy * y0 + (sy / rel_p) * py0;
        vpy[i] = (k1_y * sy * rel_p) * y0 + cy * py0;

        /* Low energy z correction (matches low_energy_z_correction.f90) */
        {
            double pz_val = vpz[i];
            if (mc2 * (beta_ref * pz_val) * (beta_ref * pz_val) < 3e-7 * e_tot_ele) {
                /* Taylor expansion for small pz — avoids precision loss */
                double mr = mc2 / e_tot_ele;  /* mass / e_tot */
                double b02 = beta_ref * beta_ref;
                double f_tay = b02 * (2.0 * b02 - mr * mr * 0.5);
                vz[i] += step_len * pz_val * (1.0 - 1.5 * pz_val * b02 + pz_val * pz_val * f_tay) * mr * mr;
            } else {
                vz[i] += step_len * (beta_val - beta_ref) / beta_ref;
            }
        }

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
        /* d_ea2/d_eb2 are precomputed WITHOUT 1/beta; apply per-particle */
        if (has_elec) {
            double kx, ky;
            multipole_kick_dev(d_ea2, d_eb2, ix_elec_max, d_cm,
                               vx[i], vy[i], &kx, &ky);
            kx /= beta_val;
            ky /= beta_val;
            double scl = (istep == n_step) ? 0.5 : 1.0;
            double px_old = vpx[i], py_old = vpy[i], pz_old = vpz[i];
            double kx_s = scl * kx, ky_s = scl * ky;
            vpx[i] += kx_s;
            vpy[i] += ky_s;
            /* pz update for electric kick */
            double alpha = (kx_s * (2.0*px_old + kx_s) + ky_s * (2.0*py_old + ky_s)) / ((1.0 + pz_old) * (1.0 + pz_old));
            if (alpha < -1.0) { state[i] = LOST_PZ; return; }
            double dpz = (1.0 + pz_old) * sqrt_one_dev(alpha);
            vpz[i] = pz_old + dpz;
            double new_beta = (1.0 + vpz[i]) / sqrt((1.0 + vpz[i])*(1.0 + vpz[i]) + (mc2/p0c_val)*(mc2/p0c_val));
            vz[i] = vz[i] * new_beta / beta_val;
            beta_val = new_beta;
            beta_arr[i] = beta_val;
        }
    }

    /* Time update */
    t_arr[i] = t_start + delta_ref_time + (z_start - vz[i]) / (beta_val * C_LIGHT);
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
    size_t ib = (size_t)n * sizeof(int);

    /* Host -> Device */
    cudaMemcpy(d_vec[0], h_vx,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[1], h_vpx, db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[2], h_vy,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[3], h_vpy, db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[4], h_vz,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[5], h_vpz, db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_state,  h_state, ib, cudaMemcpyHostToDevice);
    cudaMemcpy(d_beta,   h_beta,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_p0c,    h_p0c,   db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_s,      h_s,     db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_t,      h_t,     db, cudaMemcpyHostToDevice);

    /* Launch kernel */
    int threads = 256;
    int blocks  = (n + threads - 1) / threads;
    drift_kernel<<<blocks, threads>>>(
        d_vec[0], d_vec[1], d_vec[2], d_vec[3], d_vec[4], d_vec[5],
        d_state, d_beta, d_p0c, d_s, d_t, mc2, length, n);
    cudaDeviceSynchronize();

    /* Device -> Host */
    cudaMemcpy(h_vx,    d_vec[0], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vpx,   d_vec[1], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vy,    d_vec[2], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vpy,   d_vec[3], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vz,    d_vec[4], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vpz,   d_vec[5], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_state,  d_state,  ib, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_s,      d_s,      db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_t,      d_t,      db, cudaMemcpyDeviceToHost);
}

/* Cached device buffers for multipole data */
static double *d_a2  = NULL;
static double *d_b2  = NULL;
static double *d_ea2 = NULL;
static double *d_eb2 = NULL;
static double *d_cm  = NULL;

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

    size_t db = (size_t)n_particles * sizeof(double);
    size_t ib = (size_t)n_particles * sizeof(int);
    size_t multi_sz = N_MULTI * sizeof(double);
    size_t cm_sz    = N_MULTI * N_MULTI * sizeof(double);

    /* Host -> Device: particle data */
    cudaMemcpy(d_vec[0], h_vx,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[1], h_vpx, db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[2], h_vy,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[3], h_vpy, db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[4], h_vz,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[5], h_vpz, db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_state,  h_state, ib, cudaMemcpyHostToDevice);
    cudaMemcpy(d_beta,   h_beta,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_p0c,    h_p0c,   db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_t,      h_t,     db, cudaMemcpyHostToDevice);

    /* Ensure c_multi table is always uploaded if any multipoles present */
    if (ix_mag_max >= 0 || ix_elec_max >= 0) {
        if (!d_cm) cudaMalloc((void**)&d_cm, cm_sz);
        cudaMemcpy(d_cm, h_cm, cm_sz, cudaMemcpyHostToDevice);
    }

    /* Host -> Device: magnetic multipole data */
    if (ix_mag_max >= 0) {
        if (!d_a2) cudaMalloc((void**)&d_a2, multi_sz);
        if (!d_b2) cudaMalloc((void**)&d_b2, multi_sz);
        cudaMemcpy(d_a2, h_a2, multi_sz, cudaMemcpyHostToDevice);
        cudaMemcpy(d_b2, h_b2, multi_sz, cudaMemcpyHostToDevice);
    }

    /* Host -> Device: electric multipole data */
    if (ix_elec_max >= 0) {
        if (!d_ea2) cudaMalloc((void**)&d_ea2, multi_sz);
        if (!d_eb2) cudaMalloc((void**)&d_eb2, multi_sz);
        cudaMemcpy(d_ea2, h_ea2, multi_sz, cudaMemcpyHostToDevice);
        cudaMemcpy(d_eb2, h_eb2, multi_sz, cudaMemcpyHostToDevice);
    }

    /* Launch kernel */
    int threads = 256;
    int blocks  = (n_particles + threads - 1) / threads;
    quad_kernel<<<blocks, threads>>>(
        d_vec[0], d_vec[1], d_vec[2], d_vec[3], d_vec[4], d_vec[5],
        d_state, d_beta, d_p0c, d_t,
        mc2, b1, ele_length, delta_ref_time, e_tot_ele, charge_dir,
        n_particles, d_a2, d_b2, d_cm, ix_mag_max, n_step,
        d_ea2, d_eb2, ix_elec_max);
    cudaDeviceSynchronize();

    /* Device -> Host */
    cudaMemcpy(h_vx,    d_vec[0], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vpx,   d_vec[1], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vy,    d_vec[2], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vpy,   d_vec[3], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vz,    d_vec[4], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vpz,   d_vec[5], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_state,  d_state,  ib, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_t,      d_t,      db, cudaMemcpyDeviceToHost);
    /* Electric kicks can modify beta and pz — copy back */
    if (ix_elec_max >= 0) {
        cudaMemcpy(h_beta, d_beta, db, cudaMemcpyDeviceToHost);
    }
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
    double rel_charge_dir, double charge_dir_for_multipole,
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
        double ps2 = rel_p0*rel_p0 - vpx[i]*vpx[i] - vpy[i]*vpy[i];
        double ps = sqrt(ps2) / rel_p0;
        kx *= f_gx / (ps * beta_val);
        ky *= f_gx / (ps * beta_val);
        double px_old = vpx[i], py_old = vpy[i], pz_old = vpz[i];
        double kx_h = 0.5 * kx, ky_h = 0.5 * ky;
        vpx[i] += kx_h;
        vpy[i] += ky_h;
        double alpha_e = (kx_h*(2.0*px_old+kx_h) + ky_h*(2.0*py_old+ky_h)) / (rel_p0*rel_p0);
        if (alpha_e < -1.0) { state[i] = LOST_PZ; return; }
        vpz[i] = pz_old + rel_p0 * sqrt_one_dev(alpha_e);
        double new_beta = (1.0+vpz[i]) / sqrt((1.0+vpz[i])*(1.0+vpz[i]) + (mc2/p0c_val)*(mc2/p0c_val));
        vz[i] = vz[i] * new_beta / beta_val;
        beta_val = new_beta;
        beta_arr[i] = beta_val;
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
            {
                double pz_v = vpz[i];
                if (mc2 * (beta_ref * pz_v) * (beta_ref * pz_v) < 3e-7 * e_tot_ele) {
                    double mr = mc2 / e_tot_ele;
                    double b02 = beta_ref * beta_ref;
                    double f_tay = b02 * (2.0*b02 - mr*mr*0.5);
                    vz[i] += step_len * pz_v * (1.0 - 1.5*pz_v*b02 + pz_v*pz_v*f_tay) * mr*mr;
                } else {
                    vz[i] += step_len * (beta_val - beta_ref) / beta_ref;
                }
            }

        /* ---- Branch 2: g=0 and dg=0 → pure drift ---- */
        } else if ((g == 0.0 && dg == 0.0) || step_len == 0.0) {
            double delta = vpz[i];
            double rel_pc = 1.0 + delta;
            double px_rel = vpx[i] / rel_pc;
            double py_rel = vpy[i] / rel_pc;
            double pxy2 = px_rel*px_rel + py_rel*py_rel;
            if (pxy2 >= 1.0) { state[i] = LOST_PZ; return; }
            double ps_rel = sqrt(1.0 - pxy2);
            vx[i] += step_len * px_rel / ps_rel;
            vy[i] += step_len * py_rel / ps_rel;
            double p_tot = p0c_val * rel_pc;
            double A = (mc2*mc2*(2.0*delta + delta*delta)) / (p_tot*p_tot + mc2*mc2);
            vz[i] += step_len * (sqrt_one_dev(A) + sqrt_one_dev(-pxy2)/ps_rel);
            /* time handled at end */

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
                    state[i] = 8; /* lost$ mapped to lost_pz for simplicity */
                    return;
                }

                double rad = sqrt(r_val);
                double xi;
                if (cos_plus > 0.0) {
                    double denom = rad + cos_plus;
                    xi = alpha_b / denom;
                } else {
                    xi = (rad - cos_plus) / g_p;
                }
                vx[i] = x*cos_a - step_len*step_len*g*cosc_a + xi;

                /* Check aperture limit */
                if (fabs(vx[i]) > 1.0) {
                    state[i] = 8; /* lost$ */
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
            double ps2 = rel_p0*rel_p0 - vpx[i]*vpx[i] - vpy[i]*vpy[i];
            double ps = sqrt(ps2) / rel_p0;
            kx *= f_gx / (ps * beta_val);
            ky *= f_gx / (ps * beta_val);
            double scl = (istep == n_step) ? 0.5 : 1.0;
            double px_old = vpx[i], py_old = vpy[i], pz_old = vpz[i];
            double kx_s = scl * kx, ky_s = scl * ky;
            vpx[i] += kx_s;
            vpy[i] += ky_s;
            double alpha_e = (kx_s*(2.0*px_old+kx_s) + ky_s*(2.0*py_old+ky_s)) / ((1.0+pz_old)*(1.0+pz_old));
            if (alpha_e < -1.0) { state[i] = LOST_PZ; return; }
            vpz[i] = pz_old + (1.0+pz_old) * sqrt_one_dev(alpha_e);
            double new_beta = (1.0+vpz[i]) / sqrt((1.0+vpz[i])*(1.0+vpz[i]) + (mc2/p0c_val)*(mc2/p0c_val));
            vz[i] = vz[i] * new_beta / beta_val;
            beta_val = new_beta;
            beta_arr[i] = beta_val;
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
    double rel_charge_dir, double charge_dir_for_multipole,
    double p0c_ele,
    int n_particles,
    double *h_a2, double *h_b2, double *h_cm,
    int ix_mag_max, int n_step,
    double *h_ea2, double *h_eb2, int ix_elec_max)
{
    if (ensure_buffers(n_particles) != 0) return;

    size_t db = (size_t)n_particles * sizeof(double);
    size_t ib = (size_t)n_particles * sizeof(int);
    size_t multi_sz = N_MULTI * sizeof(double);
    size_t cm_sz    = N_MULTI * N_MULTI * sizeof(double);

    /* Host -> Device: particle data */
    cudaMemcpy(d_vec[0], h_vx,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[1], h_vpx, db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[2], h_vy,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[3], h_vpy, db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[4], h_vz,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_vec[5], h_vpz, db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_state,  h_state, ib, cudaMemcpyHostToDevice);
    cudaMemcpy(d_beta,   h_beta,  db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_p0c,    h_p0c,   db, cudaMemcpyHostToDevice);
    cudaMemcpy(d_t,      h_t,     db, cudaMemcpyHostToDevice);

    /* Upload multipole data if present */
    if (ix_mag_max >= 0 || ix_elec_max >= 0) {
        if (!d_cm) cudaMalloc((void**)&d_cm, cm_sz);
        cudaMemcpy(d_cm, h_cm, cm_sz, cudaMemcpyHostToDevice);
    }
    if (ix_mag_max >= 0) {
        if (!d_a2) cudaMalloc((void**)&d_a2, multi_sz);
        if (!d_b2) cudaMalloc((void**)&d_b2, multi_sz);
        cudaMemcpy(d_a2, h_a2, multi_sz, cudaMemcpyHostToDevice);
        cudaMemcpy(d_b2, h_b2, multi_sz, cudaMemcpyHostToDevice);
    }
    if (ix_elec_max >= 0) {
        if (!d_ea2) cudaMalloc((void**)&d_ea2, multi_sz);
        if (!d_eb2) cudaMalloc((void**)&d_eb2, multi_sz);
        cudaMemcpy(d_ea2, h_ea2, multi_sz, cudaMemcpyHostToDevice);
        cudaMemcpy(d_eb2, h_eb2, multi_sz, cudaMemcpyHostToDevice);
    }

    /* Launch kernel */
    int threads = 256;
    int blocks  = (n_particles + threads - 1) / threads;
    bend_kernel<<<blocks, threads>>>(
        d_vec[0], d_vec[1], d_vec[2], d_vec[3], d_vec[4], d_vec[5],
        d_state, d_beta, d_p0c, d_t,
        mc2, g, g_tot, dg, b1, ele_length, delta_ref_time, e_tot_ele,
        rel_charge_dir, charge_dir_for_multipole, p0c_ele,
        n_particles, d_a2, d_b2, d_cm, ix_mag_max, n_step,
        d_ea2, d_eb2, ix_elec_max);
    cudaDeviceSynchronize();

    /* Device -> Host */
    cudaMemcpy(h_vx,    d_vec[0], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vpx,   d_vec[1], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vy,    d_vec[2], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vpy,   d_vec[3], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vz,    d_vec[4], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_vpz,   d_vec[5], db, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_state,  d_state,  ib, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_t,      d_t,      db, cudaMemcpyDeviceToHost);
    if (ix_elec_max >= 0) {
        cudaMemcpy(h_beta, d_beta, db, cudaMemcpyDeviceToHost);
    }
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
    d_cap = 0;
}

#endif /* USE_GPU_TRACKING */
