#include "tms_lib.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#define GL2_CUDA_CHECK(call)                                                    \
    do {                                                                        \
        cudaError_t err__ = (call);                                             \
        if (err__ != cudaSuccess) {                                             \
            std::fprintf(stderr,                                                \
                "CUDA error %s:%d: %s\n",                                      \
                __FILE__, __LINE__, cudaGetErrorString(err__));                 \
            std::abort();                                                       \
        }                                                                       \
    } while (0)

static inline double gl2_clamp_discrepancy_square(double x) {
    return x < 0.0 && x > -1e-10 ? 0.0 : x;
}

__global__
static void gl2_build_Y_soa_kernel(const double* __restrict__ points,
                                   double* __restrict__ Y,
                                   int npts,
                                   int dim,
                                   int stride_dim) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int total = npts * dim;

    if (tid >= total) return;

    const int i = tid / dim;
    const int d = tid - i * dim;

    // SoA layout on GPU:
    // Y[d * npts + i] = 2 - points[i * stride_dim + d].
    Y[(size_t)d * (size_t)npts + (size_t)i] =
        2.0 - points[(size_t)i * (size_t)stride_dim + (size_t)d];
}

__global__
static void gl2_term1_diag_kernel(const double* __restrict__ Y,
                                  double* __restrict__ partial_term1,
                                  double* __restrict__ partial_diag,
                                  int npts,
                                  int dim) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= npts) return;

    double prod1 = 1.0;
    double prod_diag = 1.0;

    for (int d = 0; d < dim; ++d) {
        const double y = Y[(size_t)d * (size_t)npts + (size_t)i];
        const double x = 2.0 - y;

        prod1 *= 0.5 * (3.0 - x * x);
        prod_diag *= y;
    }

    partial_term1[i] = prod1;
    partial_diag[i] = prod_diag;
}

__global__
static void gl2_term2_upper_per_i_kernel(const double* __restrict__ Y,
                                         double* __restrict__ partial_upper,
                                         int npts,
                                         int dim,
                                         int i_begin,
                                         int i_count) {
    const int local_i = blockIdx.x;
    if (local_i >= i_count) return;

    const int i = i_begin + local_i;

    extern __shared__ double shared[];
    double* yi = shared;
    double* reduce = shared + dim;

    if (i >= npts - 1) {
        if (threadIdx.x == 0) partial_upper[i] = 0.0;
        return;
    }

    // Cache point i in shared memory. This avoids re-reading Y[d*npts+i]
    // for every k handled by the block.
    for (int d = threadIdx.x; d < dim; d += blockDim.x) {
        yi[d] = Y[(size_t)d * (size_t)npts + (size_t)i];
    }

    __syncthreads();

    double local_sum = 0.0;

    // Threads within the block cover k = i+1 .. npts-1.
    // For fixed d, threads in a warp read consecutive Y[d*npts+k], so reads
    // are coalesced.
    for (int k = i + 1 + threadIdx.x; k < npts; k += blockDim.x) {
        double prod = 1.0;

        for (int d = 0; d < dim; ++d) {
            const double a = yi[d];
            const double b = Y[(size_t)d * (size_t)npts + (size_t)k];
            prod *= a < b ? a : b;
        }

        local_sum += prod;
    }

    reduce[threadIdx.x] = local_sum;
    __syncthreads();

    // blockDim.x must be a power of two for this simple reduction.
    for (int offset = blockDim.x >> 1; offset > 0; offset >>= 1) {
        if (threadIdx.x < offset) {
            reduce[threadIdx.x] += reduce[threadIdx.x + offset];
        }
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        partial_upper[i] = reduce[0];
    }
}

static bool gl2_is_power_of_two(int x) {
    return x > 0 && (x & (x - 1)) == 0;
}

extern "C"
double generalized_l2_discrepancy_squared_cuda(
    const double* h_points,
    int npts,
    int dim,
    int stride_dim,
    int threads_per_block,
    int i_chunk_size) {

    assert(h_points != nullptr);
    assert(npts > 0);
    assert(dim > 0);

    if (stride_dim == 0) stride_dim = dim;
    if (threads_per_block <= 0) threads_per_block = 256;
    if (i_chunk_size <= 0) i_chunk_size = 1024;

    assert(stride_dim >= dim);
    assert(gl2_is_power_of_two(threads_per_block));

    const double invN = 1.0 / (double)npts;

    const size_t points_size =
        (size_t)npts * (size_t)stride_dim * sizeof(double);
    const size_t Y_size =
        (size_t)npts * (size_t)dim * sizeof(double);
    const size_t partial_size =
        (size_t)npts * sizeof(double);

    double* d_points = nullptr;
    double* d_Y = nullptr;
    double* d_partial_term1 = nullptr;
    double* d_partial_diag = nullptr;
    double* d_partial_upper = nullptr;

    GL2_CUDA_CHECK(cudaMalloc(&d_points, points_size));
    GL2_CUDA_CHECK(cudaMalloc(&d_Y, Y_size));
    GL2_CUDA_CHECK(cudaMalloc(&d_partial_term1, partial_size));
    GL2_CUDA_CHECK(cudaMalloc(&d_partial_diag, partial_size));
    GL2_CUDA_CHECK(cudaMalloc(&d_partial_upper, partial_size));

    GL2_CUDA_CHECK(cudaMemcpy(d_points,
                              h_points,
                              points_size,
                              cudaMemcpyHostToDevice));

    // Build Y = 2 - x in SoA form.
    {
        const int total = npts * dim;
        const int block = 256;
        const int grid = (total + block - 1) / block;

        gl2_build_Y_soa_kernel<<<grid, block>>>(
            d_points,
            d_Y,
            npts,
            dim,
            stride_dim);

        GL2_CUDA_CHECK(cudaGetLastError());
    }

    // First linear term and diagonal of the quadratic term.
    {
        const int block = 256;
        const int grid = (npts + block - 1) / block;

        gl2_term1_diag_kernel<<<grid, block>>>(
            d_Y,
            d_partial_term1,
            d_partial_diag,
            npts,
            dim);

        GL2_CUDA_CHECK(cudaGetLastError());
    }

    // Upper triangular quadratic term.
    // A kernel launch handles a chunk of i-values. This avoids extremely long
    // single launches on Windows display GPUs.
    {
        const int shared_bytes =
            (int)(((size_t)dim + (size_t)threads_per_block) * sizeof(double));

        for (int i_begin = 0; i_begin < npts; i_begin += i_chunk_size) {
            const int i_count = std::min(i_chunk_size, npts - i_begin);

            gl2_term2_upper_per_i_kernel<<<i_count,
                                           threads_per_block,
                                           shared_bytes>>>(
                d_Y,
                d_partial_upper,
                npts,
                dim,
                i_begin,
                i_count);

            GL2_CUDA_CHECK(cudaGetLastError());
        }
    }

    GL2_CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<double> h_term1(npts);
    std::vector<double> h_diag(npts);
    std::vector<double> h_upper(npts);

    GL2_CUDA_CHECK(cudaMemcpy(h_term1.data(),
                              d_partial_term1,
                              partial_size,
                              cudaMemcpyDeviceToHost));
    GL2_CUDA_CHECK(cudaMemcpy(h_diag.data(),
                              d_partial_diag,
                              partial_size,
                              cudaMemcpyDeviceToHost));
    GL2_CUDA_CHECK(cudaMemcpy(h_upper.data(),
                              d_partial_upper,
                              partial_size,
                              cudaMemcpyDeviceToHost));

    GL2_CUDA_CHECK(cudaFree(d_points));
    GL2_CUDA_CHECK(cudaFree(d_Y));
    GL2_CUDA_CHECK(cudaFree(d_partial_term1));
    GL2_CUDA_CHECK(cudaFree(d_partial_diag));
    GL2_CUDA_CHECK(cudaFree(d_partial_upper));

    long double term0 = 1.0L;
    for (int d = 0; d < dim; ++d) {
        term0 *= 4.0L / 3.0L;
    }

    long double term1_sum = 0.0L;
    long double term2_diag = 0.0L;
    long double term2_upper = 0.0L;

    for (int i = 0; i < npts; ++i) {
        term1_sum += (long double)h_term1[i];
        term2_diag += (long double)h_diag[i];
        term2_upper += (long double)h_upper[i];
    }

    const long double term2_sum = term2_diag + 2.0L * term2_upper;

    const long double d2 =
        term0
        - 2.0L * (long double)invN * term1_sum
        + (long double)invN * (long double)invN * term2_sum;

    return gl2_clamp_discrepancy_square((double)d2);
}

extern "C"
double generalized_l2_discrepancy_cuda(
    const double* h_points,
    int npts,
    int dim,
    int stride_dim,
    int threads_per_block,
    int i_chunk_size) {

    return std::sqrt(generalized_l2_discrepancy_squared_cuda(
        h_points,
        npts,
        dim,
        stride_dim,
        threads_per_block,
        i_chunk_size));
}

extern "C"
double generalized_l2_discrepancy_squared_cuda_default(
    const double* h_points,
    int npts,
    int dim) {

    return generalized_l2_discrepancy_squared_cuda(
        h_points,
        npts,
        dim,
        dim,
        256,
        1024);
}

extern "C"
double generalized_l2_discrepancy_cuda_default(
    const double* h_points,
    int npts,
    int dim) {

    return std::sqrt(generalized_l2_discrepancy_squared_cuda_default(
        h_points,
        npts,
        dim));
}
