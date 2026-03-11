/*
 * cufft_wrapper.c
 *
 * C wrapper around NVIDIA cuFFT for 3D complex-to-complex (Z2Z) double-precision FFTs.
 * Called from Fortran via iso_c_binding, replacing the FFTW backend in fft_interface_mod.f90.
 *
 * Build requirements:
 *   - CUDA Toolkit (cuda_runtime.h, cufft.h)
 *   - Link with: -lcufft -lcudart
 *   - Compile with: -DUSE_CUFFT
 *
 * The wrapper caches a single cuFFT plan and device memory allocation for the most recent
 * grid size, avoiding costly re-planning on every call (the main solver calls ccfft3d 7+
 * times per step with the same dimensions).
 */

#ifdef USE_CUFFT

#include <cuda_runtime.h>
#include <cufft.h>
#include <stdio.h>
#include <string.h>

/* --------------------------------------------------------------------------
 * Cached state: one plan + one pair of device buffers
 * -------------------------------------------------------------------------- */
static cufftHandle cached_plan    = 0;
static int         cached_n1      = 0;
static int         cached_n2      = 0;
static int         cached_n3      = 0;
static int         plan_valid     = 0;

static cufftDoubleComplex *d_in   = NULL;
static cufftDoubleComplex *d_out  = NULL;
static size_t              d_nelems = 0;   /* number of complex elements */

/* --------------------------------------------------------------------------
 * cufft_gpu_available  –  query whether a CUDA GPU is present
 * Returns 1 if at least one GPU is found, 0 otherwise.
 * -------------------------------------------------------------------------- */
int cufft_gpu_available_(void)          /* trailing underscore for Fortran/gfortran */
{
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    return (err == cudaSuccess && count > 0) ? 1 : 0;
}

/* --------------------------------------------------------------------------
 * ensure_resources  –  (re-)create plan and device buffers when grid size changes
 * Returns 0 on success, -1 on failure.
 * -------------------------------------------------------------------------- */
static int ensure_resources(int n1, int n2, int n3)
{
    size_t nelems = (size_t)n1 * (size_t)n2 * (size_t)n3;

    if (plan_valid && cached_n1 == n1 && cached_n2 == n2 && cached_n3 == n3)
        return 0;   /* nothing to do */

    /* ------- tear down old resources ------- */
    if (plan_valid) {
        cufftDestroy(cached_plan);
        plan_valid = 0;
    }

    if (nelems != d_nelems) {
        if (d_in)  { cudaFree(d_in);  d_in  = NULL; }
        if (d_out) { cudaFree(d_out); d_out = NULL; }

        size_t bytes = nelems * sizeof(cufftDoubleComplex);
        cudaError_t e1 = cudaMalloc((void **)&d_in,  bytes);
        cudaError_t e2 = cudaMalloc((void **)&d_out, bytes);
        if (e1 != cudaSuccess || e2 != cudaSuccess) {
            fprintf(stderr, "[cufft_wrapper] cudaMalloc failed (%zu bytes): %s / %s\n",
                    bytes, cudaGetErrorString(e1), cudaGetErrorString(e2));
            if (d_in)  { cudaFree(d_in);  d_in  = NULL; }
            if (d_out) { cudaFree(d_out); d_out = NULL; }
            d_nelems = 0;
            return -1;
        }
        d_nelems = nelems;
    }

    /* ------- create plan ------- */
    /*
     * Fortran array(n1,n2,n3) is column-major, which is identical in memory
     * layout to C array[n3][n2][n1].  cufftPlan3d(plan, nx, ny, nz) expects
     * the outermost, middle, innermost C-order dimensions, so we pass (n3,n2,n1).
     */
    cufftResult res = cufftPlan3d(&cached_plan, n3, n2, n1, CUFFT_Z2Z);
    if (res != CUFFT_SUCCESS) {
        fprintf(stderr, "[cufft_wrapper] cufftPlan3d(%d,%d,%d) failed: error %d\n",
                n3, n2, n1, (int)res);
        return -1;
    }

    plan_valid = 1;
    cached_n1  = n1;
    cached_n2  = n2;
    cached_n3  = n3;
    return 0;
}

/* --------------------------------------------------------------------------
 * cufft_ccfft3d  –  3D complex-to-complex double-precision FFT
 *
 * Fortran signature (iso_c_binding):
 *   subroutine cufft_ccfft3d(a, b, idir1, n1, n2, n3) bind(C)
 *     complex(C_DOUBLE_COMPLEX), intent(in)  :: a(n1,n2,n3)
 *     complex(C_DOUBLE_COMPLEX), intent(out) :: b(n1,n2,n3)
 *     integer(C_INT), value, intent(in)      :: idir1, n1, n2, n3
 *   end subroutine
 *
 * idir1 semantics (matching the existing FFTW convention in fft_interface_mod):
 *   idir1 ==  1  →  "backward" / inverse  →  CUFFT_INVERSE  (+1 exponent sign)
 *   idir1 == -1  →  "forward"             →  CUFFT_FORWARD  (-1 exponent sign)
 * -------------------------------------------------------------------------- */
void cufft_ccfft3d_(const cufftDoubleComplex *a,
                          cufftDoubleComplex *b,
                    int idir1, int n1, int n2, int n3)
{
    size_t nelems = (size_t)n1 * (size_t)n2 * (size_t)n3;
    size_t bytes  = nelems * sizeof(cufftDoubleComplex);

    /* Map direction: idir(1)==1 → backward/inverse, else forward */
    int direction = (idir1 == 1) ? CUFFT_INVERSE : CUFFT_FORWARD;

    /* Set up plan + device memory */
    if (ensure_resources(n1, n2, n3) != 0) {
        fprintf(stderr, "[cufft_wrapper] Resource setup failed – copying input to output as fallback\n");
        if ((const void *)a != (const void *)b)
            memcpy(b, a, bytes);
        return;
    }

    /* Host → Device */
    cudaMemcpy(d_in, a, bytes, cudaMemcpyHostToDevice);

    /* Execute */
    cufftResult res = cufftExecZ2Z(cached_plan, d_in, d_out, direction);
    if (res != CUFFT_SUCCESS) {
        fprintf(stderr, "[cufft_wrapper] cufftExecZ2Z failed: error %d\n", (int)res);
    }

    /* Device → Host */
    cudaMemcpy(b, d_out, bytes, cudaMemcpyDeviceToHost);
}

/* --------------------------------------------------------------------------
 * cufft_cleanup  –  release all cached GPU resources
 * Call from Fortran at program shutdown (optional but tidy).
 * -------------------------------------------------------------------------- */
void cufft_cleanup_(void)
{
    if (plan_valid) {
        cufftDestroy(cached_plan);
        plan_valid = 0;
    }
    if (d_in)  { cudaFree(d_in);  d_in  = NULL; }
    if (d_out) { cudaFree(d_out); d_out = NULL; }
    d_nelems = 0;
    cached_n1 = cached_n2 = cached_n3 = 0;
}

#endif /* USE_CUFFT */
