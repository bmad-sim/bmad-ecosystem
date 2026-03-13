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

        /* Low energy z correction */
        vz[i] += step_len * (beta_val - beta_ref) / beta_ref;

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
