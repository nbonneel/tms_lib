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





#define CUDA_CHECK(call)                                                   \
    do {                                                                   \
        cudaError_t err__ = (call);                                        \
        if (err__ != cudaSuccess) {                                        \
            std::fprintf(stderr,                                           \
                "CUDA error %s:%d: %s\n",                                  \
                __FILE__, __LINE__, cudaGetErrorString(err__));            \
            std::abort();                                                  \
        }                                                                  \
    } while (0)

static inline double sanitize_unit_coord_host(double x) {
    if (!std::isfinite(x)) return 0.0;
    if (x < 0.0) return 0.0;
    if (x >= 1.0) return std::nextafter(1.0, 0.0);
    return x;
}

static unsigned long long checked_mul_u64(unsigned long long a,
                                          unsigned long long b) {
    if (a != 0 && b > ~0ULL / a) return ~0ULL;
    return a * b;
}

static void build_candidate_coordinates_host(
    const double* points,
    int npts,
    int dim,
    int stride_dim,
    std::vector<double>& coords_flat,
    std::vector<int>& offsets,
    std::vector<int>& sizes,
    unsigned long long& n_candidates) {
    if (stride_dim == 0) stride_dim = dim;

    offsets.resize(dim + 1);
    sizes.resize(dim);
    coords_flat.clear();

    n_candidates = 1ULL;

    for (int d = 0; d < dim; ++d) {
        std::vector<double> c;
        c.reserve(size_t(npts) + 1);

        for (int i = 0; i < npts; ++i) {
            c.push_back(sanitize_unit_coord_host(points[size_t(i) * stride_dim + d]));
        }
        c.push_back(1.0);

        std::sort(c.begin(), c.end());
        c.erase(std::unique(c.begin(), c.end()), c.end());

        offsets[d] = static_cast<int>(coords_flat.size());
        sizes[d] = static_cast<int>(c.size());
        coords_flat.insert(coords_flat.end(), c.begin(), c.end());

        n_candidates = checked_mul_u64(n_candidates, static_cast<unsigned long long>(sizes[d]));
    }

    offsets[dim] = static_cast<int>(coords_flat.size());
}

int star_discrepancy_cuda_exhaustive_feasible(
    const double* points,
    int npts,
    int dim,
    int stride_dim,
    unsigned long long max_candidate_boxes) {
    if (!points || npts <= 0 || dim <= 0) return 0;
    if (stride_dim == 0) stride_dim = dim;
    if (stride_dim < dim) return 0;

    std::vector<double> coords_flat;
    std::vector<int> offsets;
    std::vector<int> sizes;
    unsigned long long n_candidates = 0;

    build_candidate_coordinates_host(points, npts, dim, stride_dim,
                                     coords_flat, offsets, sizes,
                                     n_candidates);

    return n_candidates <= max_candidate_boxes ? 1 : 0;
}

__device__ inline double atomicMaxDoubleCuda(double* address, double val) {
    unsigned long long int* address_as_ull =
        reinterpret_cast<unsigned long long int*>(address);
    unsigned long long int old = *address_as_ull;
    unsigned long long int assumed;

    do {
        assumed = old;
        double old_val = __longlong_as_double(static_cast<long long>(assumed));
        if (old_val >= val) break;
        old = atomicCAS(address_as_ull,
                        assumed,
                        static_cast<unsigned long long int>(__double_as_longlong(val)));
    } while (assumed != old);

    return __longlong_as_double(static_cast<long long>(old));
}

__global__ void star_exhaustive_kernel(
    const double* __restrict__ points,
    const double* __restrict__ coords_flat,
    const int* __restrict__ offsets,
    const int* __restrict__ sizes,
    int npts,
    int dim,
    int stride_dim,
    unsigned long long n_candidates,
    double invN,
    double* best_out) {

    unsigned long long tid =
        static_cast<unsigned long long>(blockIdx.x) * blockDim.x + threadIdx.x;
    unsigned long long step =
        static_cast<unsigned long long>(gridDim.x) * blockDim.x;

    // Supports runtime dimension. For each candidate rank, decode mixed radix.
    for (unsigned long long rank = tid; rank < n_candidates; rank += step) {
        unsigned long long r = rank;

        // Small fixed upper bound for stack-local candidates. Increase if needed.
        double u[64];
        if (dim > 64) return;

        double volume = 1.0;

        for (int d = dim - 1; d >= 0; --d) {
            int sd = sizes[d];
            int idx = static_cast<int>(r % static_cast<unsigned long long>(sd));
            r /= static_cast<unsigned long long>(sd);

            double x = coords_flat[offsets[d] + idx];
            u[d] = x;
            volume *= x;
        }

        int count_less = 0;
        int count_leq = 0;

        for (int i = 0; i < npts; ++i) {
            bool less = true;
            bool leq = true;

            for (int d = 0; d < dim; ++d) {
                double p = points[static_cast<size_t>(i) * stride_dim + d];

                // points were sanitized on the host for candidate generation,
                // but the original input is used here. Apply equivalent clipping.
                if (p < 0.0) p = 0.0;
                if (p >= 1.0) p = nextafter(1.0, 0.0);

                if (!(p < u[d])) less = false;
                if (!(p <= u[d])) leq = false;

                if (!less && !leq) break;
            }

            count_less += less ? 1 : 0;
            count_leq += leq ? 1 : 0;
        }

        double dplus = static_cast<double>(count_leq) * invN - volume;
        double dminus = volume - static_cast<double>(count_less) * invN;
        double b = dplus > dminus ? dplus : dminus;

        if (b > 0.0) {
            atomicMaxDoubleCuda(best_out, b);
        }
    }
}

double star_discrepancy_exact_cuda_exhaustive(
    const double* points,
    int npts,
    int dim,
    int stride_dim,
    unsigned long long max_candidate_boxes) {

    assert(points != nullptr);
    assert(npts > 0);
    assert(dim > 0);
    if (stride_dim == 0) stride_dim = dim;
    assert(stride_dim >= dim);

    std::vector<double> coords_flat;
    std::vector<int> offsets;
    std::vector<int> sizes;
    unsigned long long n_candidates = 0;

    build_candidate_coordinates_host(points, npts, dim, stride_dim,
                                     coords_flat, offsets, sizes,
                                     n_candidates);

    if (n_candidates > max_candidate_boxes) {
        std::fprintf(stderr,
                     "star_discrepancy_exact_cuda_exhaustive: too many candidate boxes: %llu > %llu.\n"
                     "Use the CPU branch-and-bound path for this case.\n",
                     n_candidates, max_candidate_boxes);
        return -1.0;
    }

    double* d_points = nullptr;
    double* d_coords = nullptr;
    int* d_offsets = nullptr;
    int* d_sizes = nullptr;
    double* d_best = nullptr;

    const size_t points_bytes = size_t(npts) * size_t(stride_dim) * sizeof(double);
    const size_t coords_bytes = coords_flat.size() * sizeof(double);
    const size_t offsets_bytes = offsets.size() * sizeof(int);
    const size_t sizes_bytes = sizes.size() * sizeof(int);

    CUDA_CHECK(cudaMalloc(&d_points, points_bytes));
    CUDA_CHECK(cudaMalloc(&d_coords, coords_bytes));
    CUDA_CHECK(cudaMalloc(&d_offsets, offsets_bytes));
    CUDA_CHECK(cudaMalloc(&d_sizes, sizes_bytes));
    CUDA_CHECK(cudaMalloc(&d_best, sizeof(double)));

    double zero = 0.0;
    CUDA_CHECK(cudaMemcpy(d_points, points, points_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_coords, coords_flat.data(), coords_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offsets, offsets.data(), offsets_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_sizes, sizes.data(), sizes_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_best, &zero, sizeof(double), cudaMemcpyHostToDevice));

    int block = 128;
    int grid = 0;
    if (n_candidates < 65535ULL * block) {
        grid = static_cast<int>((n_candidates + block - 1) / block);
        if (grid < 1) grid = 1;
    } else {
        grid = 65535;
    }

    star_exhaustive_kernel<<<grid, block>>>(
        d_points,
        d_coords,
        d_offsets,
        d_sizes,
        npts,
        dim,
        stride_dim,
        n_candidates,
        1.0 / static_cast<double>(npts),
        d_best);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    double best = 0.0;
    CUDA_CHECK(cudaMemcpy(&best, d_best, sizeof(double), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_coords));
    CUDA_CHECK(cudaFree(d_offsets));
    CUDA_CHECK(cudaFree(d_sizes));
    CUDA_CHECK(cudaFree(d_best));

    return best;
}
