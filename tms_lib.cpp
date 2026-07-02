#include "tms_lib.h"

#include <limits>
#include <cstdio>
#include <random>
#include <map>
#include <unordered_map>

#include <cstdint>
#include <cstdio>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif



inline int popcount_u64(uint64_t x) {
#if defined(_MSC_VER)
    return static_cast<int>(__popcnt64(x));
#else
    return __builtin_popcountll(x);
#endif
}


char invGalois[MAX_GF + 1][MAX_GF] = { // we set inv(0)=0
    {}, {}, // bases 0 and 1
    {0, 1}, // base 2
    {0, 1, 2}, // base 3
    {0, 1, 3, 2}, // base 4 : GF(2)[x]/(x^2+x+1)
    {0, 1, 3, 2, 4}, // base 5
    {}, // base 6 : not a Galois Field, not p^r
    {0, 1, 4, 5, 2, 3, 6}, // base 7
    {0, 1, 5, 6, 7, 2, 3, 4}, // base 8 : GF(2)[x]/(x^3+x+1)
    {0, 1, 2, 6, 5, 4, 3, 8, 7}, // base 9 : GF(3)[x]/(x^2+1)
    {}, // base 10 : not a Galois Field, not p^r
    {0, 1, 6, 4, 3, 9, 2, 8, 7, 5, 10}, // base 11
    {}, // base 12 : not a Galois Field, not p^r
    {0, 1, 7, 9, 10, 8, 11, 2, 5, 3, 4, 6, 12}, // base 13
    {}, // base 14 : not a Galois Field, not p^r
    {}, // base 15 : not a Galois Field, not p^r
    {0, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8} // base 16 : GF(2)[x]/(x^4+x+1)
};



char add_non_prime[MAX_GF + 1][MAX_GF][MAX_GF] = {
    {{}}, {{}},     // bases 0 and 1
    {{}},           // base 2 : prime basis, use modulo
    {{}},           // base 3 : prime basis, use modulo

    {   // base 4 : GF(2)[x]/(x^2+x+1)
        {0,1,2,3},
        {1,0,3,2},
        {2,3,0,1},
        {3,2,1,0}
    },

    {{}},           // base 5 : prime basis, use modulo
    {{}},           // base 6 : not a Galois Field, not p^r
    {{}},           // base 7 : prime basis, use modulo

    {   // base 8 : GF(2)[x]/(x^3+x+1)
        {0,1,2,3,4,5,6,7},
        {1,0,3,2,5,4,7,6},
        {2,3,0,1,6,7,4,5},
        {3,2,1,0,7,6,5,4},
        {4,5,6,7,0,1,2,3},
        {5,4,7,6,1,0,3,2},
        {6,7,4,5,2,3,0,1},
        {7,6,5,4,3,2,1,0}
    },

    {   // base 9 : GF(3)[x]/(x^2+1)
        {0,1,2,3,4,5,6,7,8},
        {1,2,0,4,5,3,7,8,6},
        {2,0,1,5,3,4,8,6,7},
        {3,4,5,6,7,8,0,1,2},
        {4,5,3,7,8,6,1,2,0},
        {5,3,4,8,6,7,2,0,1},
        {6,7,8,0,1,2,3,4,5},
        {7,8,6,1,2,0,4,5,3},
        {8,6,7,2,0,1,5,3,4}
    },

    {{}},           // base 10 : not a Galois Field, not p^r
    {{}},           // base 11 : prime basis, use modulo
    {{}},           // base 12 : not a Galois Field, not p^r
    {{}},           // base 13 : prime basis, use modulo
    {{}},           // base 14 : not a Galois Field, not p^r
    {{}},           // base 15 : not a Galois Field, not p^r

    {   // base 16 : GF(2)[x]/(x^4+x+1)
        {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15},
        {1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14},
        {2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13},
        {3,2,1,0,7,6,5,4,11,10,9,8,15,14,13,12},
        {4,5,6,7,0,1,2,3,12,13,14,15,8,9,10,11},
        {5,4,7,6,1,0,3,2,13,12,15,14,9,8,11,10},
        {6,7,4,5,2,3,0,1,14,15,12,13,10,11,8,9},
        {7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8},
        {8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7},
        {9,8,11,10,13,12,15,14,1,0,3,2,5,4,7,6},
        {10,11,8,9,14,15,12,13,2,3,0,1,6,7,4,5},
        {11,10,9,8,15,14,13,12,3,2,1,0,7,6,5,4},
        {12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3},
        {13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2},
        {14,15,12,13,10,11,8,9,6,7,4,5,2,3,0,1},
        {15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0}
    }
};

char neg_non_prime[MAX_GF + 1][MAX_GF] = {
    {}, {}, // bases 0 and 1
    {0, 1}, // base 2
    {0, 2, 1}, // base 3
    {0, 1, 2, 3}, // base 4 : GF(2)[x]/(x^2+x+1)
    {0, 4, 3, 2, 1}, // base 5
    {}, // base 6 : not a Galois Field, not p^r
    {0, 6, 5, 4, 3, 2, 1}, // base 7
    {0, 1, 2, 3, 4, 5, 6, 7}, // base 8 : GF(2)[x]/(x^3+x+1)
    {0, 2, 1, 6, 8, 7, 3, 5, 4}, // base 9 : GF(3)[x]/(x^2+1)
    {}, // base 10 : not a Galois Field, not p^r
    {0, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1}, // base 11
    {}, // base 12 : not a Galois Field, not p^r
    {0, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1}, // base 13
    {}, // base 14 : not a Galois Field, not p^r
    {}, // base 15 : not a Galois Field, not p^r
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15} // base 16 : GF(2)[x]/(x^4+x+1)
};

char sub_non_prime[MAX_GF + 1][MAX_GF][MAX_GF] = {
    {{}}, {{}},     // bases 0 and 1
    {{}},           // base 2 : prime basis, use modulo
    {{}},           // base 3 : prime basis, use modulo

    {   // base 4
        {0,1,2,3},
        {1,0,3,2},
        {2,3,0,1},
        {3,2,1,0}
    },

    {{}},           // base 5 : prime basis, use modulo
    {{}},           // base 6 : not a Galois Field
    {{}},           // base 7 : prime basis, use modulo

    {   // base 8
        {0,1,2,3,4,5,6,7},
        {1,0,3,2,5,4,7,6},
        {2,3,0,1,6,7,4,5},
        {3,2,1,0,7,6,5,4},
        {4,5,6,7,0,1,2,3},
        {5,4,7,6,1,0,3,2},
        {6,7,4,5,2,3,0,1},
        {7,6,5,4,3,2,1,0}
    },

    {   // base 9
        {0,2,1,6,8,7,3,5,4},
        {1,0,2,7,6,8,4,3,5},
        {2,1,0,8,7,6,5,4,3},
        {3,5,4,0,2,1,6,8,7},
        {4,3,5,1,0,2,7,6,8},
        {5,4,3,2,1,0,8,7,6},
        {6,8,7,3,5,4,0,2,1},
        {7,6,8,4,3,5,1,0,2},
        {8,7,6,5,4,3,2,1,0}
    },

    {{}},           // base 10 : not a Galois Field
    {{}},           // base 11 : prime basis, use modulo
    {{}},           // base 12 : not a Galois Field
    {{}},           // base 13 : prime basis, use modulo
    {{}},           // base 14 : not a Galois Field
    {{}},           // base 15 : not a Galois Field

    {   // base 16
        {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15},
        {1,0,3,2,5,4,7,6,9,8,11,10,13,12,15,14},
        {2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13},
        {3,2,1,0,7,6,5,4,11,10,9,8,15,14,13,12},
        {4,5,6,7,0,1,2,3,12,13,14,15,8,9,10,11},
        {5,4,7,6,1,0,3,2,13,12,15,14,9,8,11,10},
        {6,7,4,5,2,3,0,1,14,15,12,13,10,11,8,9},
        {7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8},
        {8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7},
        {9,8,11,10,13,12,15,14,1,0,3,2,5,4,7,6},
        {10,11,8,9,14,15,12,13,2,3,0,1,6,7,4,5},
        {11,10,9,8,15,14,13,12,3,2,1,0,7,6,5,4},
        {12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3},
        {13,12,15,14,9,8,11,10,5,4,7,6,1,0,3,2},
        {14,15,12,13,10,11,8,9,6,7,4,5,2,3,0,1},
        {15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0}
    }
};

char mul_non_prime[MAX_GF + 1][MAX_GF][MAX_GF] = {
    {{}}, {{}},     // bases 0 and 1
    {{}},           // base 2 : prime basis, use modulo
    {{}},           // base 3 : prime basis, use modulo

    {   // base 4 : GF(2)[x]/(x^2+x+1)
        {0,0,0,0},
        {0,1,2,3},
        {0,2,3,1},
        {0,3,1,2}
    },

    {{}},           // base 5 : prime basis, use modulo
    {{}},           // base 6 : not a Galois Field, not p^r
    {{}},           // base 7 : prime basis, use modulo

    {   // base 8 : GF(2)[x]/(x^3+x+1)
        {0,0,0,0,0,0,0,0},
        {0,1,2,3,4,5,6,7},
        {0,2,4,6,3,1,7,5},
        {0,3,6,5,7,4,1,2},
        {0,4,3,7,6,2,5,1},
        {0,5,1,4,2,7,3,6},
        {0,6,7,1,5,3,2,4},
        {0,7,5,2,1,6,4,3}
    },

    {   // base 9 : GF(3)[x]/(x^2+1)
        {0,0,0,0,0,0,0,0,0},
        {0,1,2,3,4,5,6,7,8},
        {0,2,1,6,8,7,3,5,4},
        {0,3,6,2,5,8,1,4,7},
        {0,4,8,5,6,1,7,2,3},
        {0,5,7,8,1,3,4,6,2},
        {0,6,3,1,7,4,2,8,5},
        {0,7,5,4,2,6,8,3,1},
        {0,8,4,7,3,2,5,1,6}
    },

    {{}},           // base 10 : not a Galois Field, not p^r
    {{}},           // base 11 : prime basis, use modulo
    {{}},           // base 12 : not a Galois Field, not p^r
    {{}},           // base 13 : prime basis, use modulo
    {{}},           // base 14 : not a Galois Field, not p^r
    {{}},           // base 15 : not a Galois Field, not p^r

    {   // base 16 : GF(2)[x]/(x^4+x+1)
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15},
        {0,2,4,6,8,10,12,14,3,1,7,5,11,9,15,13},
        {0,3,6,5,12,15,10,9,11,8,13,14,7,4,1,2},
        {0,4,8,12,3,7,11,15,6,2,14,10,5,1,13,9},
        {0,5,10,15,7,2,13,8,14,11,4,1,9,12,3,6},
        {0,6,12,10,11,13,7,1,5,3,9,15,14,8,2,4},
        {0,7,14,9,15,8,1,6,13,10,3,4,2,5,12,11},
        {0,8,3,11,6,14,5,13,12,4,15,7,10,2,9,1},
        {0,9,1,8,2,11,3,10,4,13,5,12,6,15,7,14},
        {0,10,7,13,14,4,9,3,15,5,8,2,1,11,6,12},
        {0,11,5,14,10,1,15,4,7,12,2,9,13,6,8,3},
        {0,12,11,7,5,9,14,2,10,6,1,13,15,3,4,8},
        {0,13,9,4,1,12,8,5,2,15,11,6,3,14,10,7},
        {0,14,15,1,13,3,2,12,9,7,6,8,4,10,11,5},
        {0,15,13,2,9,6,4,11,1,14,12,3,8,7,5,10}
    }
};

char div_non_prime[MAX_GF + 1][MAX_GF][MAX_GF] = {
    {{}}, {{}},     // bases 0 and 1
    {{}},           // base 2 : prime basis, use modulo
    {{}},           // base 3 : prime basis, use modulo

    {   // base 4
        {0,0,0,0},
        {0,1,3,2},
        {0,2,1,3},
        {0,3,2,1}
    },

    {{}},           // base 5 : prime basis, use modulo
    {{}},           // base 6 : not a Galois Field
    {{}},           // base 7 : prime basis, use modulo

    {   // base 8
        {0,0,0,0,0,0,0,0},
        {0,1,5,6,7,2,3,4},
        {0,2,1,7,5,4,6,3},
        {0,3,4,1,2,6,5,7},
        {0,4,2,5,1,3,7,6},
        {0,5,7,3,6,1,4,2},
        {0,6,3,2,4,7,1,5},
        {0,7,6,4,3,5,2,1}
    },

    {   // base 9
        {0,0,0,0,0,0,0,0,0},
        {0,1,2,6,5,4,3,8,7},
        {0,2,1,3,7,8,6,4,5},
        {0,3,6,1,8,5,2,7,4},
        {0,4,8,7,1,6,5,3,2},
        {0,5,7,4,3,1,8,2,6},
        {0,6,3,2,4,7,1,5,8},
        {0,7,5,8,6,2,4,1,3},
        {0,8,4,5,2,3,7,6,1}
    },

    {{}},           // base 10 : not a Galois Field
    {{}},           // base 11 : prime basis, use modulo
    {{}},           // base 12 : not a Galois Field
    {{}},           // base 13 : prime basis, use modulo
    {{}},           // base 14 : not a Galois Field
    {{}},           // base 15 : not a Galois Field

    {   // base 16
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,1,9,14,13,11,7,6,15,2,12,5,10,4,3,8},
        {0,2,1,15,9,5,14,12,13,4,11,10,7,8,6,3},
        {0,3,8,1,4,14,9,10,2,6,7,15,13,12,5,11},
        {0,4,2,13,1,10,15,11,9,8,5,7,14,3,12,6},
        {0,5,11,3,12,1,8,13,6,10,9,2,4,7,15,14},
        {0,6,3,2,8,15,1,7,4,12,14,13,9,11,10,5},
        {0,7,10,12,5,4,6,1,11,14,2,8,3,15,9,13},
        {0,8,4,9,2,7,13,5,1,3,10,14,15,6,11,12},
        {0,9,13,7,15,12,10,3,14,1,6,11,5,2,8,4},
        {0,10,5,6,11,2,3,9,12,7,1,4,8,14,13,15},
        {0,11,12,8,6,9,4,15,3,5,13,1,2,10,14,7},
        {0,12,6,4,3,13,2,14,8,11,15,9,1,5,7,10},
        {0,13,15,10,14,6,5,8,7,9,3,12,11,1,4,2},
        {0,14,7,11,10,8,12,2,5,15,4,3,6,13,1,9},
        {0,15,14,5,7,3,11,4,10,13,8,6,12,9,2,1}
    }
};




inline double clamp_discrepancy_square(double x) {
    return x < 0.0 && x > -1e-12 ? 0.0 : x;
}




inline double sanitize_unit_coord(double x) {
    assert(std::isfinite(x));

    if (x < 0.0) {
        return 0.0;
    }

    if (x >= 1.0) {
        return std::nextafter(1.0, 0.0);
    }

    return x;
}




double generalized_l2_discrepancy_squared_exact_runtime_aos(
    const double* points,
    int npts,
    size_t dim,
    int block_size = 64
) {
    assert(points != 0);
    assert(npts > 0);
    assert(dim > 0);
    assert(block_size > 0);

    const double invN = 1.0 / double(npts);

    // AoS layout:
    // Y[i * dim + d] = 2 - points[i * dim + d]
    std::vector<double> Y(std::size_t(npts) * std::size_t(dim));

#pragma omp parallel for
    for (int i = 0; i < npts; ++i) {
        const double* p = points + std::size_t(i) * std::size_t(dim);
        double* y = Y.data() + std::size_t(i) * std::size_t(dim);

        for (int d = 0; d < dim; ++d) {
            y[d] = 2.0 - p[d];
        }
    }

    double term0 = 1.0;
    for (int d = 0; d < dim; ++d) {
        term0 *= 4.0 / 3.0;
    }

    double term1_sum = 0.0;

#pragma omp parallel for reduction(+:term1_sum)
    for (int i = 0; i < npts; ++i) {
        const double* y = Y.data() + std::size_t(i) * std::size_t(dim);

        double prod = 1.0;

        for (int d = 0; d < dim; ++d) {
            const double x = 2.0 - y[d];
            prod *= 0.5 * (3.0 - x * x);
        }

        term1_sum += prod;
    }

    double term2_diag = 0.0;

#pragma omp parallel for reduction(+:term2_diag)
    for (int i = 0; i < npts; ++i) {
        const double* y = Y.data() + std::size_t(i) * dim;

        double prod = 1.0;

        for (int d = 0; d < dim; ++d) {
            prod *= y[d];
        }

        term2_diag += prod;
    }

    const int nb = (npts + block_size - 1) / block_size;

    double term2_upper = 0.0;

#pragma omp parallel for schedule(dynamic, 1) reduction(+:term2_upper)
    for (int bi = 0; bi < nb; ++bi) {
        const int i0 = bi * block_size;
        const int i1 = std::min(i0 + block_size, npts);

        double local = 0.0;


		for (size_t i = i0; i < i1; ++i) {
			const double* yi = Y.data() + i * dim;

			for (size_t k = i + 1; k < i1; ++k) {
				const double* yk = Y.data() + k * dim;

				double prod = 1.0;
				for (size_t d = 0; d < dim; ++d) {
					const double a = yi[d];
					const double b = yk[d];
					prod *= a < b ? a : b;
				}
				local += prod;
			}
		}
        

        for (int bj = bi+1; bj < nb; ++bj) {
            const int k0 = bj * block_size;
            const int k1 = std::min(k0 + block_size, npts);

                for (size_t i = i0; i < i1; ++i) {
                    const double* yi = Y.data() + i *dim;

                    for (size_t k = k0; k < k1; ++k) {
                        const double* yk = Y.data() + k * dim;

                        double prod = 1.0;
                        for (size_t d = 0; d < dim; ++d) { 
                            const double a = yi[d];
                            const double b = yk[d];
                            prod *= a < b ? a : b;
                        }

                        local += prod;
                    }
                }            
        }

        term2_upper += local;
    }

    const double term2_sum = term2_diag + 2.0 * term2_upper;

    const double d2 =
        term0
        - 2.0 * invN * term1_sum
        + invN * invN * term2_sum;

    return clamp_discrepancy_square(d2);
}


double generalized_l2_discrepancy_squared_auto(
    const double* points,
    int npts,
    int dim,
    int block_size, 
    bool use_gpu
) {
#ifdef TMS_USE_CUDA
    if (use_gpu)
        return generalized_l2_discrepancy_squared_cuda(points, npts, dim, dim);
    else
        return generalized_l2_discrepancy_squared_exact_runtime_aos(points, npts, dim, block_size);
#else
    return generalized_l2_discrepancy_squared_exact_runtime_aos(points, npts, dim, block_size);
#endif
}


// O(n^2 d)
double generalized_l2_discrepancy(const double* points, int npts, int dim, int block_size, bool use_gpu) {
    return std::sqrt(  generalized_l2_discrepancy_squared_auto(points, npts, dim, block_size, use_gpu)  );
}


inline double gl2_pair_kernel_y_aos(const double* yi,
    const double* yk,
    int dim) {
    // yi[d] = 2 - x_i[d]
    // yk[d] = 2 - x_k[d]
    //
    // prod_d (2 - max(x_i[d], x_k[d]))
    // = prod_d min(2 - x_i[d], 2 - x_k[d])
    double prod = 1.0;

    for (int d = 0; d < dim; ++d) {
        const double a = yi[d];
        const double b = yk[d];
        prod *= a < b ? a : b;
    }

    return prod;
}

void generalized_l2_discrepancy_curve_squared_exact_runtime_aos(
    const double* points,
    long long max_npts,
    int dim,
    const long long* sample_counts,
    long long n_samples,
    double* out_d2,
    int stride_dim = 0,
    int block_size = 128,
    int internal_num_threads = 1
) {
    assert(points != 0);
    assert(max_npts > 0);
    assert(dim > 0);
    assert(sample_counts != 0);
    assert(n_samples > 0);
    assert(out_d2 != 0);
    assert(block_size > 0);

    if (stride_dim == 0) {
        stride_dim = dim;
    }

    assert(stride_dim >= dim);

    for (int s = 0; s < n_samples; ++s) {
        assert(sample_counts[s] > 0);
        assert(sample_counts[s] <= max_npts);

        if (s > 0) {
            assert(sample_counts[s] >= sample_counts[s - 1]);
        }
    }

    bool do_omp = (internal_num_threads>1);

#ifndef _OPENMP
    do_omp = false;
#endif

    // Compact AoS layout:
    // Y[i * dim + d] = 2 - points[i * stride_dim + d].
    std::vector<double> Y(max_npts * dim);

    std::vector<double> incA(max_npts, 0.0);
    std::vector<double> incB(max_npts, 0.0);

#pragma omp parallel for schedule(static) if(do_omp)
    for (long long i = 0; i < max_npts; ++i) {
        const double* p =
            points + i * stride_dim;

        double* y =
            Y.data() + i * dim;

        double a_i = 1.0;
        double b_ii = 1.0;

        for (int d = 0; d < dim; ++d) {
            const double x = p[d];
            const double yy = 2.0 - x;

            y[d] = yy;

            a_i *= 0.5 * (3.0 - x * x);
            b_ii *= yy;
        }

        incA[i] = a_i;

        // Start incB with the diagonal contribution b_ii.
        incB[i] = b_ii;
    }

    const int nb =
        (max_npts + block_size - 1) / block_size;

    // Compute:
    //
    // incB[i] = b_ii + 2 * sum_{k<i} b_ik
    //
    // but with block locality:
    //
    // for each output i-block bi:
    //     accumulate all previous k-blocks bj <= bi
    //
    // This is race-free because each bi writes a disjoint range of incB.
#pragma omp parallel for schedule(dynamic, 1) if(do_omp) num_threads(internal_num_threads)
    for (int bi = 0; bi < nb; ++bi) {
        const long long i0 = bi * block_size;
        const long long i1 = std::min(i0 + block_size, max_npts);

        std::vector<double> local_pair_sum(i1 - i0, 0.0);

        for (int bj = 0; bj <= bi; ++bj) {
            const long long k0 = bj * block_size;
            const long long k1 = std::min(k0 + block_size, max_npts);

            if (bj < bi) {
                // Full rectangular block: k < i automatically.
                for (long long i = i0; i < i1; ++i) {
                    const double* yi =
                        Y.data() + i * dim;

                    double sum_i = 0.0;

                    for (long long k = k0; k < k1; ++k) {
                        const double* yk =
                            Y.data() + k *dim;

                        sum_i += gl2_pair_kernel_y_aos(yi, yk, dim);
                    }

                    local_pair_sum[i - i0] += sum_i;
                }
            }
            else {
                // Diagonal block: keep only k < i.
                for (int i = i0; i < i1; ++i) {
                    const double* yi =
                        Y.data() + std::size_t(i) * std::size_t(dim);

                    double sum_i = 0.0;

                    for (int k = i0; k < i; ++k) {
                        const double* yk =
                            Y.data() + std::size_t(k) * std::size_t(dim);

                        sum_i += gl2_pair_kernel_y_aos(yi, yk, dim);
                    }

                    local_pair_sum[std::size_t(i - i0)] += sum_i;
                }
            }
        }

        for (int i = i0; i < i1; ++i) {
            incB[i] += 2.0 * local_pair_sum[std::size_t(i - i0)];
        }
    }

    double term0 = 1.0;

    for (int d = 0; d < dim; ++d) {
        term0 *= 4.0 / 3.0;
    }

    // Prefix accumulation:
    //
    // A_N = sum_{i<N} incA[i]
    // B_N = sum_{i<N} incB[i]
    double A = 0.0;
    double B = 0.0;

    int sample_id = 0;

    for (int i = 0; i < max_npts && sample_id < n_samples; ++i) {
        A += incA[i];
        B += incB[i];

        const int n = i + 1;

        while (sample_id < n_samples &&
            sample_counts[sample_id] == n) {
            const double invN = 1.0 / double(n);

            const double d2 =
                term0
                - 2.0 * invN * A
                + invN * invN * B;

            out_d2[sample_id] = clamp_discrepancy_square(d2);

            ++sample_id;
        }
    }

    assert(sample_id == n_samples);
}


void generalized_l2_discrepancy_curve(const double* points, long long max_npts, int dim, const long long* sample_counts,
    long long n_samples, double* out_D, int stride_dim, int block_size, int internal_num_threads) {
    generalized_l2_discrepancy_curve_squared_exact_runtime_aos(points, max_npts, dim, sample_counts, n_samples, out_D,
        stride_dim, block_size, internal_num_threads);

    for (long long s = 0; s < n_samples; ++s) {
        out_D[s] = std::sqrt(out_D[s]);
    }
}


struct StarDiscrepancyBB {
    const double* points;
    int npts;
    int dim;

    int words;
    double invN;

    std::vector<double> P;
    std::vector<std::vector<double> > coords;

    std::vector<std::vector<uint64_t> > less_masks;
    std::vector<std::vector<uint64_t> > leq_masks;
    std::vector<uint64_t> all_ones;
    std::vector<uint64_t> tmp_mask;

    double best;

    StarDiscrepancyBB(const double* points_,
        int npts_,
        int dim_)
        : points(points_),
        npts(npts_),
        dim(dim_),
        words((npts_ + 63) >> 6),
        invN(1.0 / double(npts_)),
        best(0.0) {
        assert(points != 0);
        assert(npts > 0);
        assert(dim > 0);

        build();
    }

    void build() {
        P.resize(npts * dim);

        for (int i = 0; i < npts; ++i) {
            for (int j = 0; j < dim; ++j) {
                P[i * dim + j] =
                    sanitize_unit_coord(points[i * dim + j]);
            }
        }

        coords.resize(dim);

        for (int j = 0; j < dim; ++j) {
            coords[j].reserve(npts + 1);

            for (int i = 0; i < npts; ++i) {
                coords[j].push_back(P[i * dim + j]);
            }

            coords[j].push_back(1.0);

            std::sort(coords[j].begin(), coords[j].end());

            coords[j].erase(
                std::unique(coords[j].begin(), coords[j].end()),
                coords[j].end()
            );
        }

        all_ones.assign(words, ~uint64_t(0));

        const int extra_bits = words * 64 - npts;
        if (extra_bits > 0) {
            all_ones[words - 1] >>= extra_bits;
        }

        tmp_mask.resize(words);

        less_masks.resize(dim);
        leq_masks.resize(dim);

        for (int j = 0; j < dim; ++j) {
            const int nj = static_cast<int>(coords[j].size());

            less_masks[j].assign(nj * words, uint64_t(0));
            leq_masks[j].assign(nj * words, uint64_t(0));

            for (int c = 0; c < nj; ++c) {
                const double x = coords[j][c];

                uint64_t* less = &less_masks[j][c * words];
                uint64_t* leq = &leq_masks[j][c * words];

                for (int i = 0; i < npts; ++i) {
                    const double p = P[i * dim + j];

                    const int w = i >> 6;
                    const uint64_t bit = uint64_t(1) << (i & 63);

                    if (p < x) {
                        less[w] |= bit;
                    }

                    if (p <= x) {
                        leq[w] |= bit;
                    }
                }
            }
        }
    }

    int count_intersection(const std::vector<int>& idx,
        bool use_leq) {
        for (int w = 0; w < words; ++w) {
            tmp_mask[w] = all_ones[w];
        }

        for (int j = 0; j < dim; ++j) {
            const std::vector<uint64_t>& masks =
                use_leq ? leq_masks[j] : less_masks[j];

            const uint64_t* mj = &masks[idx[j] * words];

            for (int w = 0; w < words; ++w) {
                tmp_mask[w] &= mj[w];
            }
        }

        int count = 0;

        for (int w = 0; w < words; ++w) {
            count += popcount_u64(tmp_mask[w]);
        }

        return count;
    }

    double volume_from_indices(const std::vector<int>& idx) const {
        double vol = 1.0;

        for (int j = 0; j < dim; ++j) {
            vol *= coords[j][idx[j]];
        }

        return vol;
    }

    void evaluate_single_box(const std::vector<int>& idx) {
        const double volume = volume_from_indices(idx);

        const int count_less = count_intersection(idx, false);
        const int count_leq = count_intersection(idx, true);

        const double dplus =
            double(count_leq) * invN - volume;

        const double dminus =
            volume - double(count_less) * invN;

        if (dplus > best) {
            best = dplus;
        }

        if (dminus > best) {
            best = dminus;
        }
    }

    void recurse(std::vector<int>& lo,
        std::vector<int>& hi) {
        // Compute a rigorous upper bound for all boxes in this block.
        const double vol_min = volume_from_indices(lo);
        const double vol_max = volume_from_indices(hi);

        const int count_leq_max =
            count_intersection(hi, true);

        const int count_less_min =
            count_intersection(lo, false);

        const double ub_plus =
            double(count_leq_max) * invN - vol_min;

        const double ub_minus =
            vol_max - double(count_less_min) * invN;

        const double ub =
            ub_plus > ub_minus ? ub_plus : ub_minus;

        if (ub <= best + 1e-15) {
            return;
        }

        // Check whether the block is a single grid point.
        int split_dim = -1;
        int largest_range = 0;

        for (int j = 0; j < dim; ++j) {
            const int r = hi[j] - lo[j];

            if (r > largest_range) {
                largest_range = r;
                split_dim = j;
            }
        }

        if (split_dim < 0) {
            evaluate_single_box(lo);
            return;
        }

        const int d = split_dim;
        const int mid = (lo[d] + hi[d]) >> 1;

        // Visit the side with larger coordinates first.
        // It often improves best quickly via Dminus.
        {
            const int old_lo = lo[d];
            lo[d] = mid + 1;
            recurse(lo, hi);
            lo[d] = old_lo;
        }

        {
            const int old_hi = hi[d];
            hi[d] = mid;
            recurse(lo, hi);
            hi[d] = old_hi;
        }
    }

    double compute() {
        std::vector<int> lo(dim);
        std::vector<int> hi(dim);

        for (int j = 0; j < dim; ++j) {
            lo[j] = 0;
            hi[j] = static_cast<int>(coords[j].size()) - 1;
        }

        // Evaluate the full upper corner first.
        // This gives a nonzero lower bound immediately.
        evaluate_single_box(hi);

        recurse(lo, hi);

        return best;
    }
};

double star_discrepancy_exact_bb(const double* points,
    int npts,
    int dim) {
    StarDiscrepancyBB solver(points, npts, dim);
    return solver.compute();
}



struct StarDiscrepancy4D_BB {
    int npts;
    int words;
    double invN;

    std::vector<double> P;
    std::vector<double> coords[4];

    std::vector<uint64_t> less_masks[4];
    std::vector<uint64_t> leq_masks[4];
    std::vector<uint64_t> all_ones;
    std::vector<uint64_t> tmp;

    double best;

    StarDiscrepancy4D_BB(const double* points, int npts_)
        : npts(npts_),
        words((npts_ + 63) >> 6),
        invN(1.0 / double(npts_)),
        best(0.0) {
        P.resize(npts * 4);

        for (int i = 0; i < npts; ++i) {
            for (int j = 0; j < 4; ++j) {
                P[4 * i + j] = sanitize_unit_coord(points[4 * i + j]);
            }
        }

        build_coords();
        build_masks();

        tmp.resize(words);

        count_cache_leq.resize(COUNT_CACHE_SIZE);
        count_cache_less.resize(COUNT_CACHE_SIZE);

        for (int i = 0; i < COUNT_CACHE_SIZE; ++i) {
            count_cache_leq[i].key = ~0ULL;
            count_cache_leq[i].value = 0;

            count_cache_less[i].key = ~0ULL;
            count_cache_less[i].value = 0;
        }
    }

    struct CountCacheEntry {
        unsigned long long key;
        int value;
    };

    static const int COUNT_CACHE_BITS = 16;
    static const int COUNT_CACHE_SIZE = 1 << COUNT_CACHE_BITS;

    std::vector<CountCacheEntry> count_cache_leq;
    std::vector<CountCacheEntry> count_cache_less;

    void build_coords() {
        for (int j = 0; j < 4; ++j) {
            coords[j].reserve(npts + 1);

            for (int i = 0; i < npts; ++i) {
                coords[j].push_back(P[4 * i + j]);
            }

            coords[j].push_back(1.0);

            std::sort(coords[j].begin(), coords[j].end());
            coords[j].erase(
                std::unique(coords[j].begin(), coords[j].end()),
                coords[j].end()
            );
        }
    }

    void build_masks() {
        all_ones.assign(words, ~uint64_t(0));

        const int extra_bits = words * 64 - npts;
        if (extra_bits > 0) {
            all_ones[words - 1] >>= extra_bits;
        }

        for (int d = 0; d < 4; ++d) {
            const int nd = static_cast<int>(coords[d].size());

            less_masks[d].assign(nd * words, uint64_t(0));
            leq_masks[d].assign(nd * words, uint64_t(0));

            for (int c = 0; c < nd; ++c) {
                const double x = coords[d][c];

                uint64_t* less = &less_masks[d][c * words];
                uint64_t* leq = &leq_masks[d][c * words];

                for (int i = 0; i < npts; ++i) {
                    const double p = P[4 * i + d];
                    const int w = i >> 6;
                    const uint64_t bit = uint64_t(1) << (i & 63);

                    if (p < x)  less[w] |= bit;
                    if (p <= x) leq[w] |= bit;
                }
            }
        }
    }

    unsigned long long pack_key4(int a, int b, int c, int d) const {
        // Works if each coordinate index is < 65536.
        // For 1400 points, this is fine.
        return
            (unsigned long long)(unsigned int)a |
            ((unsigned long long)(unsigned int)b << 16) |
            ((unsigned long long)(unsigned int)c << 32) |
            ((unsigned long long)(unsigned int)d << 48);
    }

    int count4_cached(std::vector<uint64_t>* masks,
        std::vector<CountCacheEntry>& cache,
        int a, int b, int c, int d) {
        const unsigned long long key = pack_key4(a, b, c, d);

        unsigned int h = static_cast<unsigned int>(key ^ (key >> 32));
        h *= 2654435761u;
        h &= COUNT_CACHE_SIZE - 1;

        CountCacheEntry& e = cache[h];

        if (e.key == key) {
            return e.value;
        }

        const int v = count4(masks, a, b, c, d);

        e.key = key;
        e.value = v;

        return v;
    }

    int count4(const std::vector<uint64_t>* masks,
        int a, int b, int c, int d) {
        const uint64_t* m0 = &masks[0][a * words];
        const uint64_t* m1 = &masks[1][b * words];
        const uint64_t* m2 = &masks[2][c * words];
        const uint64_t* m3 = &masks[3][d * words];

        int count = 0;

        for (int w = 0; w < words; ++w) {
            const uint64_t x = m0[w] & m1[w] & m2[w] & m3[w];
            count += popcount_u64(x);
        }

        return count;
    }

    double volume4(int a, int b, int c, int d) const {
        return coords[0][a] * coords[1][b] * coords[2][c] * coords[3][d];
    }

    void eval_box(int a, int b, int c, int d) {
        const double vol = volume4(a, b, c, d);

        const int cless = count4(less_masks, a, b, c, d);
        const int cleq = count4(leq_masks, a, b, c, d);

        const double dplus = double(cleq) * invN - vol;
        const double dminus = vol - double(cless) * invN;

        if (dplus > best)  best = dplus;
        if (dminus > best) best = dminus;
    }

    void recurse(int a0, int a1,
        int b0, int b1,
        int c0, int c1,
        int d0, int d1) {
        const double vol_min = volume4(a0, b0, c0, d0);
        const double vol_max = volume4(a1, b1, c1, d1);

        const double eps = 1e-15;

        if (1.0 - vol_min <= best + eps && vol_max <= best + eps) {
            return;
        }

        bool must_recurse = false;

        if (1.0 - vol_min > best + eps) {
            const int cleq_max =
                count4_cached(leq_masks, count_cache_leq, a1, b1, c1, d1);

            const double ub_plus =
                double(cleq_max) * invN - vol_min;

            if (ub_plus > best + eps) {
                must_recurse = true;
            }
        }

        if (!must_recurse && vol_max > best + eps) {
            const int cless_min =
                count4_cached(less_masks, count_cache_less, a0, b0, c0, d0);

            const double ub_minus =
                vol_max - double(cless_min) * invN;

            if (ub_minus > best + eps) {
                must_recurse = true;
            }
        }

        if (!must_recurse) {
            return;
        }

        const int ra = a1 - a0;
        const int rb = b1 - b0;
        const int rc = c1 - c0;
        const int rd = d1 - d0;

        if ((ra | rb | rc | rd) == 0) {
            eval_box(a0, b0, c0, d0);
            return;
        }

        if (ra >= rb && ra >= rc && ra >= rd) {
            const int mid = (a0 + a1) >> 1;
            recurse(mid + 1, a1, b0, b1, c0, c1, d0, d1);
            recurse(a0, mid, b0, b1, c0, c1, d0, d1);
        }
        else if (rb >= rc && rb >= rd) {
            const int mid = (b0 + b1) >> 1;
            recurse(a0, a1, mid + 1, b1, c0, c1, d0, d1);
            recurse(a0, a1, b0, mid, c0, c1, d0, d1);
        }
        else if (rc >= rd) {
            const int mid = (c0 + c1) >> 1;
            recurse(a0, a1, b0, b1, mid + 1, c1, d0, d1);
            recurse(a0, a1, b0, b1, c0, mid, d0, d1);
        }
        else {
            const int mid = (d0 + d1) >> 1;
            recurse(a0, a1, b0, b1, c0, c1, mid + 1, d1);
            recurse(a0, a1, b0, b1, c0, c1, d0, mid);
        }
    }

    struct StarNode4D {
        int a0, a1;
        int b0, b1;
        int c0, c1;
        int d0, d1;
    };

    void recurse_iterative(int a0_init, int a1_init,
        int b0_init, int b1_init,
        int c0_init, int c1_init,
        int d0_init, int d1_init) {
        const double eps = 1e-15;

        std::vector<StarNode4D> stack;
        stack.reserve(4096);

        stack.push_back({
            a0_init, a1_init,
            b0_init, b1_init,
            c0_init, c1_init,
            d0_init, d1_init
            });

        while (!stack.empty()) {
            const StarNode4D node = stack.back();
            stack.pop_back();

            const int a0 = node.a0;
            const int a1 = node.a1;
            const int b0 = node.b0;
            const int b1 = node.b1;
            const int c0 = node.c0;
            const int c1 = node.c1;
            const int d0 = node.d0;
            const int d1 = node.d1;

            const double vol_min = volume4(a0, b0, c0, d0);
            const double vol_max = volume4(a1, b1, c1, d1);

            // Cheap count-free pruning.
            if (1.0 - vol_min <= best + eps &&
                vol_max <= best + eps) {
                continue;
            }

            bool must_recurse = false;

            // Lazy upper-bound evaluation.
            if (1.0 - vol_min > best + eps) {
                const int cleq_max =
                    count4_cached(leq_masks, count_cache_leq, a1, b1, c1, d1);

                const double ub_plus =
                    double(cleq_max) * invN - vol_min;

                if (ub_plus > best + eps) {
                    must_recurse = true;
                }
            }

            if (!must_recurse && vol_max > best + eps) {
                const int cless_min =
                    count4_cached(less_masks, count_cache_less, a0, b0, c0, d0);

                const double ub_minus =
                    vol_max - double(cless_min) * invN;

                if (ub_minus > best + eps) {
                    must_recurse = true;
                }
            }

            if (!must_recurse) {
                continue;
            }

            const int ra = a1 - a0;
            const int rb = b1 - b0;
            const int rc = c1 - c0;
            const int rd = d1 - d0;

            if ((ra | rb | rc | rd) == 0) {
                eval_box(a0, b0, c0, d0);
                continue;
            }

            // Split along the largest remaining interval.
            //
            // We push the low side first and the high side second,
            // because the stack is LIFO. This means the high side is visited first.
            // This matches the recursive ordering that often improves best earlier.
            if (ra >= rb && ra >= rc && ra >= rd) {
                const int mid = (a0 + a1) >> 1;

                stack.push_back({
                    a0, mid,
                    b0, b1,
                    c0, c1,
                    d0, d1
                    });

                stack.push_back({
                    mid + 1, a1,
                    b0, b1,
                    c0, c1,
                    d0, d1
                    });
            }
            else if (rb >= rc && rb >= rd) {
                const int mid = (b0 + b1) >> 1;

                stack.push_back({
                    a0, a1,
                    b0, mid,
                    c0, c1,
                    d0, d1
                    });

                stack.push_back({
                    a0, a1,
                    mid + 1, b1,
                    c0, c1,
                    d0, d1
                    });
            }
            else if (rc >= rd) {
                const int mid = (c0 + c1) >> 1;

                stack.push_back({
                    a0, a1,
                    b0, b1,
                    c0, mid,
                    d0, d1
                    });

                stack.push_back({
                    a0, a1,
                    b0, b1,
                    mid + 1, c1,
                    d0, d1
                    });
            }
            else {
                const int mid = (d0 + d1) >> 1;

                stack.push_back({
                    a0, a1,
                    b0, b1,
                    c0, c1,
                    d0, mid
                    });

                stack.push_back({
                    a0, a1,
                    b0, b1,
                    c0, c1,
                    mid + 1, d1
                    });
            }
        }
    }

    void initialize_best() {
        const int a = static_cast<int>(coords[0].size()) - 1;
        const int b = static_cast<int>(coords[1].size()) - 1;
        const int c = static_cast<int>(coords[2].size()) - 1;
        const int d = static_cast<int>(coords[3].size()) - 1;

        eval_box(a, b, c, d);

        // Cheap warm start: evaluate boxes anchored at actual points.
        for (int i = 0; i < npts; ++i) {
            int idx[4];

            for (int j = 0; j < 4; ++j) {
                const double v = P[4 * i + j];

                std::vector<double>::const_iterator it =
                    std::lower_bound(coords[j].begin(), coords[j].end(), v);

                idx[j] = static_cast<int>(it - coords[j].begin());
            }

            eval_box(idx[0], idx[1], idx[2], idx[3]);
        }
    }

    double compute() {
        initialize_best();

        /*recurse(
            0, static_cast<int>(coords[0].size()) - 1,
            0, static_cast<int>(coords[1].size()) - 1,
            0, static_cast<int>(coords[2].size()) - 1,
            0, static_cast<int>(coords[3].size()) - 1
        );*/
        recurse_iterative(
            0, static_cast<int>(coords[0].size()) - 1,
            0, static_cast<int>(coords[1].size()) - 1,
            0, static_cast<int>(coords[2].size()) - 1,
            0, static_cast<int>(coords[3].size()) - 1
        );

        return best;
    }
};

double star_discrepancy_4d_exact_bb_fast(const double* points, int npts) {
    StarDiscrepancy4D_BB solver(points, npts);
    return solver.compute();
}


struct StarPoint2D {
    double x;
    double y;
    int y_index;
};


// O(n^2)
double star_discrepancy_2d_exact_sweep(const double* points,
    int npts) {
    assert(points != 0);
    assert(npts > 0);

    const double invN = 1.0 / double(npts);

    std::vector<StarPoint2D> P(npts);
    std::vector<double> xs;
    std::vector<double> ys;

    xs.reserve(npts + 1);
    ys.reserve(npts + 1);

    // Copy and sanitize points.
    for (int i = 0; i < npts; ++i) {
        const double x = sanitize_unit_coord(points[2 * i + 0]);
        const double y = sanitize_unit_coord(points[2 * i + 1]);

        P[i].x = x;
        P[i].y = y;
        P[i].y_index = -1;

        xs.push_back(x);
        ys.push_back(y);
    }

    // The candidate value 1 is required.
    xs.push_back(1.0);
    ys.push_back(1.0);

    std::sort(xs.begin(), xs.end());
    std::sort(ys.begin(), ys.end());

    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());

    const int nx = static_cast<int>(xs.size());
    const int ny = static_cast<int>(ys.size());

    // Assign a y-rank to every point.
    for (int i = 0; i < npts; ++i) {
        std::vector<double>::iterator it =
            std::lower_bound(ys.begin(), ys.end(), P[i].y);

        assert(it != ys.end());
        assert(*it == P[i].y);

        P[i].y_index = static_cast<int>(it - ys.begin());
    }

    // Sort points by x-coordinate.
    std::sort(P.begin(), P.end(),
        [](const StarPoint2D& a, const StarPoint2D& b) {
            if (a.x < b.x) return true;
            if (b.x < a.x) return false;
            return a.y < b.y;
        });

    std::vector<int> active_y_count(ny, 0);

    double best = 0.0;
    int next_point = 0;

    for (int ix = 0; ix < nx; ++ix) {
        const double x = xs[ix];

        // ------------------------------------------------------------
        // Dminus:
        // active set contains points with p.x < x.
        //
        // For a candidate y, the relevant count is:
        //     #{ p : p.x < x and p.y < y }
        //
        // During the scan over y, prefix_count stores points with
        // y-index strictly smaller than the current candidate.
        // ------------------------------------------------------------
        {
            int prefix_count = 0;

            for (int iy = 0; iy < ny; ++iy) {
                const double y = ys[iy];
                const double volume = x * y;

                const double dminus =
                    volume - double(prefix_count) * invN;

                if (dminus > best) {
                    best = dminus;
                }

                // Make points with p.y == current y active for larger y.
                prefix_count += active_y_count[iy];
            }
        }

        // Add all points with p.x == x.
        // After this, active set contains points with p.x <= x.
        while (next_point < npts && P[next_point].x == x) {
            ++active_y_count[P[next_point].y_index];
            ++next_point;
        }

        // ------------------------------------------------------------
        // Dplus:
        // active set contains points with p.x <= x.
        //
        // For a candidate y, the relevant count is:
        //     #{ p : p.x <= x and p.y <= y }
        //
        // During the scan over y, prefix_count is updated before
        // evaluating the discrepancy.
        // ------------------------------------------------------------
        {
            int prefix_count = 0;

            for (int iy = 0; iy < ny; ++iy) {
                prefix_count += active_y_count[iy];

                const double y = ys[iy];
                const double volume = x * y;

                const double dplus =
                    double(prefix_count) * invN - volume;

                if (dplus > best) {
                    best = dplus;
                }
            }
        }
    }

    return best;
}

double star_discrepancy(const double* points,
    int npts,
    int dim, bool gpu) {

#ifdef TMS_USE_CUDA
    const unsigned long long max_candidate_boxes = 20000000ULL;

    if (dim <= 3 &&
        star_discrepancy_cuda_exhaustive_feasible(
            points,
            npts,
            dim,
            dim,
            max_candidate_boxes)) {
        const double r =
            star_discrepancy_exact_cuda_exhaustive(
                points,
                npts,
                dim,
                dim,
                max_candidate_boxes);

        if (r >= 0.0) {
            return r;
        }
    }
#endif

    if (dim == 2) {
        return star_discrepancy_2d_exact_sweep(points, npts);
    }

    if (dim == 4) {
        return star_discrepancy_4d_exact_bb_fast(points, npts);
    }

    return star_discrepancy_exact_bb(points, npts, dim);
}

void print_point_range(const double* points, int npts, int dim) {
    double mn = points[0];
    double mx = points[0];

    for (int i = 0; i < npts * dim; ++i) {
        mn = std::min(mn, points[i]);
        mx = std::max(mx, points[i]);
    }

    std::printf("point range = [%g, %g]\n", mn, mx);
}






inline const char* discrepancy_metric_name(DiscrepancyPlotMetric metric) {
    switch (metric) {
    case DISCREPANCY_GENERALIZED_L2:
        return "Generalized L2 discrepancy";
    case DISCREPANCY_STAR:
        return "Star discrepancy";
    default:
        return "Discrepancy";
    }
}

inline double safe_log10(double x) {
    const double eps = 1e-300;
    return std::log10(x > eps ? x : eps);
}

inline long long integer_power_clamped(int base, int exponent) {
    assert(base >= 2);
    assert(exponent >= 0);

    long long r = 1;

    for (int i = 0; i < exponent; ++i) {
        if (r > std::numeric_limits<long long>::max() / base) {
            return std::numeric_limits<long long>::max();
        }

        r *= base;
    }

    return r;
}

inline void format_integer_label(long long v, char* out, int out_size) {
    if (v < 1000000LL) {
        std::snprintf(out, out_size, "%lld", v);
    }
    else {
        std::snprintf(out, out_size, "%.3g", double(v));
    }
}

inline void svg_text(FILE* f,
    double x,
    double y,
    const char* text,
    int font_size,
    const char* anchor,
    const char* fill) {
    std::fprintf(f,
        "<text x=\"%.3f\" y=\"%.3f\" "
        "font-family=\"Arial, Helvetica, sans-serif\" "
        "font-size=\"%d\" text-anchor=\"%s\" fill=\"%s\">%s</text>\n",
        x, y, font_size, anchor, fill, text);
}

inline void svg_line(FILE* f,
    double x1,
    double y1,
    double x2,
    double y2,
    const char* stroke,
    double width,
    const char* extra) {
    std::fprintf(f,
        "<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" "
        "stroke=\"%s\" stroke-width=\"%.3f\" %s/>\n",
        x1, y1, x2, y2, stroke, width, extra);
}

inline void svg_circle(FILE* f,
    double cx,
    double cy,
    double r,
    const char* fill,
    const char* stroke, double opacity) {
    std::fprintf(f,
        "<circle cx=\"%.3f\" cy=\"%.3f\" r=\"%.3f\" "
        "fill=\"%s\" stroke=\"%s\" stroke-width=\"1\" fill-opacity=\"%.3f\"/>\n",
        cx, cy, r, fill, stroke, opacity);
}



uint64_t ipow_u64_checked(int base, int exp) {
    assert(base >= 2);
    assert(exp >= 0);

    uint64_t r = 1;

    for (int i = 0; i < exp; ++i) {
        assert(r <= UINT64_MAX / uint64_t(base));
        r *= uint64_t(base);
    }

    return r;
}

inline int integer_log_exact_u64(uint64_t n, int base) {
    assert(n > 0);
    assert(base >= 2);

    uint64_t v = 1;
    int m = 0;

    while (v < n) {
        if (v > UINT64_MAX / uint64_t(base)) {
            return -1;
        }

        v *= uint64_t(base);
        ++m;
    }

    return v == n ? m : -1;
}


inline uint64_t box_cell_from_double(double x,
    uint64_t scale) {
    assert(scale > 0);
    assert(std::isfinite(x));

    if (x <= 0.0) {
        return 0;
    }

    if (x >= 1.0) {
        return scale - 1;
    }

    const double y = x * double(scale);

    // For digital points, y is often theoretically an integer.
    // Due to binary floating-point, it may be represented as n - eps.
    // In that case, floor(y) would send it to the previous box.
    const double nearest = std::floor(y + 0.5);

    const double tol =
        64.0 * std::numeric_limits<double>::epsilon()
        * std::max(1.0, std::fabs(y));

    uint64_t cell;

    if (std::fabs(y - nearest) <= tol) {
        cell = uint64_t(nearest);
    }
    else {
        cell = uint64_t(std::floor(y));
    }

    if (cell >= scale) {
        cell = scale - 1;
    }

    return cell;
}

bool test_t_value_pointset_real_sorted(const double* points,
    int npts,
    int dim,
    int stride_dim,
    int base,
    int t_factor,
    bool dbg = false) {
    assert(points != 0);
    assert(npts > 0);
    assert(dim > 0);
    assert(stride_dim >= dim);
    assert(base >= 2);

    const int m = integer_log_exact_u64(uint64_t(npts), base);

    if (m < 0) {
        return false;
    }

    assert(t_factor >= 0);
    assert(t_factor <= m);

    const int q = m - t_factor;

    const uint64_t expected_count = ipow_u64_checked(base, t_factor);
    const uint64_t expected_nboxes = ipow_u64_checked(base, q);

    std::vector<uint64_t> codes;
    codes.resize(npts);

    bool all_ok = true;

    enumerate_compositions(dim, q, [&](const std::vector<int>& k) -> bool {
        for (int ip = 0; ip < npts; ++ip) {
            uint64_t code = 0;

            for (int d = 0; d < dim; ++d) {
                const int kd = k[d];

                uint64_t cell = 0;

                if (kd > 0) {
                    const double x =
                        sanitize_unit_coord(points[ip * stride_dim + d]);

                    const uint64_t scale = ipow_u64_checked(base, kd);

                    cell = box_cell_from_double(x, scale);

                    if (cell >= scale) {
                        cell = scale - 1;
                    }
                }

                const uint64_t radix = ipow_u64_checked(base, kd);
                code = code * radix + cell;
            }

            codes[ip] = code;
        }

        std::sort(codes.begin(), codes.end());

        uint64_t n_distinct_boxes = 0;
        int pos = 0;

        while (pos < npts) {
            const uint64_t c = codes[pos];

            int next = pos + 1;

            while (next < npts && codes[next] == c) {
                ++next;
            }

            const uint64_t count = uint64_t(next - pos);

            if (count != expected_count) {
                if (dbg) {
                    std::printf("Fail t=%d: box %llu count=%llu expected=%llu\n",
                        t_factor,
                        (unsigned long long)c,
                        (unsigned long long)count,
                        (unsigned long long)expected_count);

                    std::printf("composition k = ");
                    for (int d = 0; d < dim; ++d) {
                        std::printf("%d ", k[d]);
                    }
                    std::printf("\n");
                }

                all_ok = false;
                return false;
            }

            ++n_distinct_boxes;
            pos = next;
        }

        if (n_distinct_boxes != expected_nboxes) {
            if (dbg) {
                std::printf("Fail t=%d: distinct boxes=%llu expected boxes=%llu\n",
                    t_factor,
                    (unsigned long long)n_distinct_boxes,
                    (unsigned long long)expected_nboxes);

                std::printf("composition k = ");
                for (int d = 0; d < dim; ++d) {
                    std::printf("%d ", k[d]);
                }
                std::printf("\n");
            }

            all_ok = false;
            return false;
        }

        return true;
        });

    return all_ok;
}


int t_value_pointset(const double* points,
    int npts,
    int dim,
    int stride_dim,
    int base,
    bool dbg = false) {
    assert(points != 0);
    assert(npts > 0);
    assert(dim > 0);
    assert(stride_dim >= dim);
    assert(base >= 2);

    const int m = integer_log_exact_u64(uint64_t(npts), base);

    if (m < 0) {
        if (dbg) {
            std::printf("npts=%d is not a power of base=%d\n", npts, base);
        }

        return -1;
    }

    for (int t = 0; t <= m; ++t) {
        if (dbg) {
            std::printf("testing t=%d\n", t);
        }

        if (test_t_value_pointset_real_sorted(points,
            npts,
            dim,
            stride_dim,
            base,
            t,
            dbg)) {
            return t;
        }
    }

    return m;
}

int t_value_pointset(const double* points,
    int npts,
    int dim,
    int base,
    bool dbg) {
    return t_value_pointset(points, npts, dim, dim, base, dbg);
}



inline void random_permutation(int* p,
    int n,
    FastRNG& rng) {
    for (int i = 0; i < n; ++i) {
        p[i] = i;
    }

    for (int i = n - 1; i > 0; --i) {
        const int j = static_cast<int>(rng.bounded(uint64_t(i + 1)));
        std::swap(p[i], p[j]);
    }
}

inline uint64_t double_coord_to_grid_integer(double x,
    int base,
    int m) {
    const uint64_t scale = ipow_u64_checked(base, m);

    if (x <= 0.0) return 0;
    if (x >= 1.0) return scale - 1;

    const double y = x * double(scale);

    // Important: use nearest integer, not floor.
    // The input point is supposed to be exactly on the base^m grid.
    uint64_t I = uint64_t(std::floor(y + 0.5));

    if (I >= scale) {
        I = scale - 1;
    }

    return I;
}

bool check_owen_tree_1d(const OwenTree1D& tree) {
    if (tree.base < 2) return false;
    if (tree.depth < 0) return false;

    uint64_t n_nodes = 1;

    for (int level = 0; level < tree.depth; ++level) {
        const std::vector<int>& perms = tree.level_perm[level];

        if (perms.size() != std::size_t(n_nodes) * std::size_t(tree.base)) {
            return false;
        }

        for (uint64_t node = 0; node < n_nodes; ++node) {
            std::vector<int> seen(tree.base, 0);

            for (int a = 0; a < tree.base; ++a) {
                const int b = perms[std::size_t(node) * tree.base + a];

                if (b < 0 || b >= tree.base) {
                    return false;
                }

                ++seen[b];
            }

            for (int b = 0; b < tree.base; ++b) {
                if (seen[b] != 1) {
                    return false;
                }
            }
        }

        n_nodes *= uint64_t(tree.base);
    }

    return true;
}

bool check_owen_tree_nd(const OwenTreeND& tree) {
    if (tree.dim <= 0) return false;
    if (tree.base < 2) return false;

    if (tree.trees.size() != std::size_t(tree.dim)) {
        return false;
    }

    for (int d = 0; d < tree.dim; ++d) {
        if (tree.trees[d].base != tree.base) return false;

        if (!check_owen_tree_1d(tree.trees[d])) {
            return false;
        }
    }

    return true;
}

OwenTree1D make_random_owen_tree_1d(int base,
    int depth,
    FastRNG& rng) {
    assert(base >= 2);
    assert(depth >= 0);

    OwenTree1D tree(base, depth);

    long long n_nodes = 1;

    for (int level = 0; level < depth; ++level) {
        tree.level_perm[level].resize(
            static_cast<std::size_t>(n_nodes) * static_cast<std::size_t>(base)
        );

        for (long long node = 0; node < n_nodes; ++node) {
            int* perm =
                &tree.level_perm[level][static_cast<std::size_t>(node) * base];

            random_permutation(perm, base, rng);
        }

        n_nodes *= base;
    }

    return tree;
}

OwenTreeND* make_random_owen_tree_nd(int dim,
    int base,
    int depth,
    uint64_t seed) {
    assert(dim > 0);
    assert(base >= 2);
    assert(depth >= 0);

    FastRNG rng(seed);

    OwenTreeND* tree = new OwenTreeND(dim, base, depth);

    for (int d = 0; d < dim; ++d) {
        tree->trees[d] = make_random_owen_tree_1d(base, depth, rng);
    }

    return tree;
}


double apply_owen_1d_real(double x,
    const OwenTree1D& tree,
    int m) {
    assert(tree.base >= 2);
    assert(m >= 0);
    //assert(check_owen_tree_1d(tree));

    const int base = tree.base;
    const uint64_t scale = ipow_u64_checked(base, m);

    uint64_t I = double_coord_to_grid_integer(x, base, m);

    std::vector<int> digits(m, 0);

    // Extract m base-b digits, most significant first.
    for (int j = m - 1; j >= 0; --j) {
        digits[j] = int(I % uint64_t(base));
        I /= uint64_t(base);
    }

    uint64_t prefix = 0;
    uint64_t J = 0;

    const int active_depth =
        tree.depth < m ? tree.depth : m;

    for (int level = 0; level < m; ++level) {
        const int digit = digits[level];
        int new_digit = digit;

        if (level < active_depth) {
            const std::vector<int>& perms = tree.level_perm[level];

            new_digit =
                perms[std::size_t(prefix) * std::size_t(base)
                + std::size_t(digit)];
        }

        J = J * uint64_t(base) + uint64_t(new_digit);

        // Important: prefix is the original prefix, as in your Mathematica code.
        prefix = prefix * uint64_t(base) + uint64_t(digit);
    }

    return double(J) / double(scale);
}



void apply_owen_permutation_real(const double* points_in,
    double* points_out,
    int npts,
    int dim,
    int stride_in,
    int stride_out,
    int m,
    const OwenTreeND& tree) {
    assert(points_in != 0);
    assert(points_out != 0);
    assert(npts >= 0);
    assert(dim > 0);
    assert(stride_in >= dim);
    assert(stride_out >= dim);
    assert(tree.dim == dim);
    //assert(check_owen_tree_nd(tree));

    for (int i = 0; i < npts; ++i) {
        for (int d = 0; d < dim; ++d) {
            const double x = points_in[i * stride_in + d];

            points_out[i * stride_out + d] =
                apply_owen_1d_real(x, tree.trees[d], m);
        }
    }
}



void apply_owen_permutation_real(const double* points_in,
    double* points_out,
    int npts,
    int dim,
    int m,
    const OwenTreeND& tree) {
    apply_owen_permutation_real(points_in,
        points_out,
        npts,
        dim,
        dim,
        dim,
        m,
        tree);
}

void padd_least_significant_digits(double* pts, long long n_pts, int dim, int base, int m, long long seed) {

    std::mt19937_64 engine;
    engine.seed(seed);
    std::uniform_real_distribution<double> unif(0, 1);

    double scale = std::pow((double)base, -m);

    for (int i = 0; i < n_pts * dim; i++) {
        pts[i] += unif(engine) * scale;
    }

}


///////////////////// plotting routines

inline const ProjectionHighlight* find_projection_highlight(
    int d0,
    int d1,
    const ProjectionHighlight* highlights,
    int n_highlights
) {
    for (int i = 0; i < n_highlights; ++i) {
        const int a = highlights[i].dim0;
        const int b = highlights[i].dim1;

        if ((a == d0 && b == d1) || (a == d1 && b == d0)) {
            return &highlights[i];
        }
    }
    return 0;
}

inline void svg_begin(FILE* f, int width, int height) {
    std::fprintf(f,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        "version=\"1.1\" width=\"%d\" height=\"%d\" "
        "viewBox=\"0 0 %d %d\">\n",
        width, height, width, height);

    std::fprintf(f,
        "<rect x=\"0\" y=\"0\" width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n",
        width, height);
}

inline void svg_end(FILE* f) {
    std::fprintf(f, "</svg>\n");
}

inline void svg_rect(FILE* f,
    double x,
    double y,
    double w,
    double h,
    const char* fill,
    const char* stroke,
    double stroke_width) {
    std::fprintf(f,
        "<rect x=\"%.3f\" y=\"%.3f\" width=\"%.3f\" height=\"%.3f\" "
        "fill=\"%s\" stroke=\"%s\" stroke-width=\"%.3f\"/>\n",
        x, y, w, h, fill, stroke, stroke_width);
}



bool draw_2d_projections_svg(
    double* points, int dim, long long n_points, const char* filename, const ProjectionHighlight* highlights,
    int n_highlights, int cell_size, int cell_inner_margin, int outer_margin, int label_band, double point_radius,
    const char* point_color, double point_opacity, const char* default_border_color, double default_border_width,
    int label_font_size
) {
    assert(dim >= 2);
    assert(filename != 0);

    if (point_radius < 0) point_radius = sqrt(std::abs(point_radius) * cell_size * cell_size / (3.1416 * n_points));
    const int ncols = dim - 1;
    const int nrows = dim - 1;

    const int width = 2 * outer_margin + label_band + ncols * cell_size;
    const int height = 2 * outer_margin + label_band + nrows * cell_size;

    FILE* f = std::fopen(filename, "w");
    if (!f) {
        return false;
    }

    svg_begin(f, width, height);

    for (int d = 0; d < dim - 1; ++d) {
        char txt[32];
        std::snprintf(txt, sizeof(txt), "%d", d);
        svg_text(f,
            outer_margin + label_band + d * cell_size + cell_size * 0.5,
            height - outer_margin + 4.0,
            txt,
            label_font_size,
            "middle",
            "#444444");
    }

    for (int d = 1; d < dim; ++d) {
        char txt[32];
        std::snprintf(txt, sizeof(txt), "%d", d);
        svg_text(f,
            outer_margin + label_band - 6.0,
            outer_margin + (d - 1) * cell_size + cell_size * 0.5 + 4.0,
            txt,
            label_font_size,
            "end",
            "#444444");
    }

    for (int d1 = 1; d1 < dim; ++d1) {
        for (int d0 = 0; d0 < d1; ++d0) {

            const double x0 = outer_margin + label_band + d0 * cell_size;
            const double y0 = outer_margin + (d1 - 1) * cell_size;

            const ProjectionHighlight* hl =
                find_projection_highlight(d0, d1, highlights, n_highlights);

            if (!hl) {
                svg_rect(f,
                    x0, y0,
                    (double)cell_size, (double)cell_size,
                    "#ffffff",
                    default_border_color,
                    default_border_width);
            }
        }
    }

    // first normal squares and *then* highlights
    for (int d1 = 1; d1 < dim; ++d1) {
        for (int d0 = 0; d0 < d1; ++d0) {

            const double x0 = outer_margin + label_band + d0 * cell_size;
            const double y0 = outer_margin + (d1 - 1) * cell_size;

            const ProjectionHighlight* hl =
                find_projection_highlight(d0, d1, highlights, n_highlights);

            if (hl) {
                svg_rect(f,
                    x0, y0,
                    (double)cell_size, (double)cell_size,
                    "#ffffff",
                    hl->color,
                    hl->stroke_width);
            }
        }
    }

    for (int d1 = 1; d1 < dim; ++d1) {
        for (int d0 = 0; d0 < d1; ++d0) {

            const double x0 = outer_margin + label_band + d0 * cell_size;
            const double y0 = outer_margin + (d1 - 1) * cell_size;

            const ProjectionHighlight* hl =
                find_projection_highlight(d0, d1, highlights, n_highlights);


            const double px0 = x0 + cell_inner_margin;
            const double py0 = y0 + cell_inner_margin;
            const double pw = cell_size - 2.0 * cell_inner_margin;
            const double ph = cell_size - 2.0 * cell_inner_margin;

            for (long long p = 0; p < n_points; ++p) {
                double u = points[(size_t)p * (size_t)dim + (size_t)d0];
                double v = points[(size_t)p * (size_t)dim + (size_t)d1];

                if (u < 0.0) u = 0.0;
                if (u > 1.0) u = 1.0;
                if (v < 0.0) v = 0.0;
                if (v > 1.0) v = 1.0;

                const double cx = px0 + u * pw;
                const double cy = py0 + (1.0 - v) * ph;

                svg_circle(f, cx, cy, point_radius, point_color, point_color, point_opacity);
            }
        }
    }

    svg_end(f);
    std::fclose(f);
    return true;
}



//// plotting curves




struct PlotRect {
    double x0;
    double y0;
    double x1;
    double y1;
};

struct RGB {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

struct ResolvedCurve {
    const DiscrepancyCurve* curve;
    std::string color;
    double stroke_width;
    double opacity;
    double point_radius;
    double point_opacity;
    std::string group_key;
    std::string legend_label;
};

struct ResolvedReference {
    const ReferenceCurveSpec* ref;
    std::string color;
    double stroke_width;
    double scale;
};

struct LegendEntry {
    std::string label;
    std::string color;
    double width;
    double opacity;
    bool dashed, show_point;
    std::string dash_array;
};



double clamp_positive(double x) {
    return x > 1e-300 ? x : 1e-300;
}

bool nearly_integer(double x, double eps = 1e-9) {
    return std::fabs(x - std::round(x)) <= eps;
}


std::string format_number(double value, PlotTickNumberFormat fmt) {
    char buf[128];

    switch (fmt) {
    case PLOT_TICK_INTEGER:
        std::snprintf(buf, sizeof(buf), "%.0f", value);
        break;

    case PLOT_TICK_SCIENTIFIC:
        std::snprintf(buf, sizeof(buf), "%.2e", value);
        break;

    case PLOT_TICK_COMPACT:
        if (std::fabs(value) < 1000000.0 && std::fabs(value) >= 0.001) {
            if (std::fabs(value - std::round(value)) < 1e-8) {
                std::snprintf(buf, sizeof(buf), "%.0f", value);
            }
            else {
                std::snprintf(buf, sizeof(buf), "%.4g", value);
            }
        }
        else {
            std::snprintf(buf, sizeof(buf), "%.3g", value);
        }
        break;

    case PLOT_TICK_LOG_VALUE:
    default:
        std::snprintf(buf, sizeof(buf), "%.3g", value);
        break;
    }

    return std::string(buf);
}

std::string group_key_for_curve(const DiscrepancyCurve& c, int idx) {
    if (!c.legend_group.empty()) {
        return c.legend_group;
    }

    if (!c.label.empty()) {
        return c.label;
    }

    char buf[64];
    std::snprintf(buf, sizeof(buf), "__curve_%d", idx);
    return std::string(buf);
}

bool has_symmetric_or_asymmetric_errors(const DiscrepancyCurve& c) {
    return !c.error_minus.empty() || !c.error_plus.empty();
}

double curve_error_minus(const DiscrepancyCurve& c, int i) {
    if (!c.error_minus.empty() && i < (int)c.error_minus.size()) {
        return c.error_minus[i];
    }

    if (!c.error_plus.empty() && i < (int)c.error_plus.size()) {
        return c.error_plus[i];
    }

    return 0.0;
}

double curve_error_plus(const DiscrepancyCurve& c, int i) {
    if (!c.error_plus.empty() && i < (int)c.error_plus.size()) {
        return c.error_plus[i];
    }

    if (!c.error_minus.empty() && i < (int)c.error_minus.size()) {
        return c.error_minus[i];
    }

    return 0.0;
}

bool curve_show_point(const DiscrepancyCurve& c, int i) {
    if (!c.show_point_mask.empty()) {
        return i < (int)c.show_point_mask.size() && c.show_point_mask[i] != 0;
    }

    return c.show_points;
}

bool curve_show_error(const DiscrepancyCurve& c, int i) {
    if (c.error_style == PLOT_ERROR_NONE || !has_symmetric_or_asymmetric_errors(c)) {
        return false;
    }

    if (!c.error_mask.empty()) {
        return i < (int)c.error_mask.size() && c.error_mask[i] != 0;
    }

    return true;
}

RGB parse_color(const std::string& s, RGB fallback) {
    if (s.size() != 7 || s[0] != '#') {
        return fallback;
    }

    char* end = 0;
    const long v = std::strtol(s.c_str() + 1, &end, 16);

    if (end == 0 || *end != '\0') {
        return fallback;
    }

    RGB c;
    c.r = (unsigned char)((v >> 16) & 255);
    c.g = (unsigned char)((v >> 8) & 255);
    c.b = (unsigned char)(v & 255);
    return c;
}

const char* default_ref_color(int idx) {
    static const char* colors[4] = {
        "#D36C6C",
        "#6C8CD3",
        "#6CA87A",
        "#888888"
    };

    return colors[idx % 4];
}

void svg_text(
    FILE* f,
    double x,
    double y,
    const std::string& text,
    int font_size,
    const char* anchor,
    const char* fill,
    double rotation_degrees = 0.0)
{
    if (std::fabs(rotation_degrees) > 1e-12) {
        std::fprintf(
            f,
            "  <text x=\"%.3f\" y=\"%.3f\" font-family=\"Arial, Helvetica, sans-serif\" "
            "font-size=\"%d\" text-anchor=\"%s\" fill=\"%s\" "
            "transform=\"rotate(%.3f %.3f %.3f)\">%s</text>\n",
            x,
            y,
            font_size,
            anchor,
            fill,
            rotation_degrees,
            x,
            y,
            text.c_str());
    }
    else {
        std::fprintf(
            f,
            "  <text x=\"%.3f\" y=\"%.3f\" font-family=\"Arial, Helvetica, sans-serif\" "
            "font-size=\"%d\" text-anchor=\"%s\" fill=\"%s\">%s</text>\n",
            x,
            y,
            font_size,
            anchor,
            fill,
            text.c_str());
    }
}

void svg_line(
    FILE* f,
    double x1,
    double y1,
    double x2,
    double y2,
    const std::string& stroke,
    double width,
    double opacity,
    const std::string& extra = std::string())
{
    std::fprintf(
        f,
        "  <line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" "
        "stroke=\"%s\" stroke-width=\"%.3f\" stroke-opacity=\"%.3f\" %s/>\n",
        x1,
        y1,
        x2,
        y2,
        stroke.c_str(),
        width,
        opacity,
        extra.c_str());
}

void svg_circle(
    FILE* f,
    double cx,
    double cy,
    double r,
    const std::string& fill,
    const std::string& stroke,
    double fill_opacity,
    double stroke_opacity)
{
    std::fprintf(
        f,
        "  <circle cx=\"%.3f\" cy=\"%.3f\" r=\"%.3f\" fill=\"%s\" fill-opacity=\"%.3f\" "
        "stroke=\"%s\" stroke-opacity=\"%.3f\" stroke-width=\"1\"/>\n",
        cx,
        cy,
        r,
        fill.c_str(),
        fill_opacity,
        stroke.c_str(),
        stroke_opacity);
}

double estimate_svg_text_width(const std::string& s, int font_size) {
    // Rough but robust enough for Arial/Helvetica-like fonts.
    // This intentionally overestimates a bit to avoid legend text overflow.
    double w = 0.0;

    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(s[i]);

        if (c == ' ') {
            w += 0.32 * font_size;
        }
        else if (c == 'i' || c == 'l' || c == 'I' || c == '!' || c == '.' || c == ',') {
            w += 0.30 * font_size;
        }
        else if (c == 'm' || c == 'w' || c == 'M' || c == 'W') {
            w += 0.90 * font_size;
        }
        else if (c >= 128) {
            // UTF-8 bytes: overestimate slightly rather than trying to decode.
            w += 0.35 * font_size;
        }
        else {
            w += 0.58 * font_size;
        }
    }

    return w;
}

class Bitmap {
public:
    Bitmap(int w, int h)
        : width(w),
        height(h),
        pixels((size_t)w* (size_t)h * 3, 255) {
    }

    void blend_pixel(int x, int y, RGB c, double a) {
        if (x < 0 || x >= width || y < 0 || y >= height) {
            return;
        }

        if (a <= 0.0) {
            return;
        }

        if (a > 1.0) {
            a = 1.0;
        }

        const size_t idx = ((size_t)y * (size_t)width + (size_t)x) * 3;

        pixels[idx + 0] = (unsigned char)(pixels[idx + 0] * (1.0 - a) + c.r * a);
        pixels[idx + 1] = (unsigned char)(pixels[idx + 1] * (1.0 - a) + c.g * a);
        pixels[idx + 2] = (unsigned char)(pixels[idx + 2] * (1.0 - a) + c.b * a);
    }

    void draw_line(double x0, double y0, double x1, double y1, RGB c, double a, double width_px) {
        const double dx = x1 - x0;
        const double dy = y1 - y0;
        const double len = std::sqrt(dx * dx + dy * dy);

        if (len <= 0.0) {
            draw_disc((int)std::round(x0), (int)std::round(y0), 0.5 * width_px, c, a);
            return;
        }

        const int steps = std::max(1, (int)std::ceil(len));
        const double radius = std::max(0.5, 0.5 * width_px);

        for (int i = 0; i <= steps; ++i) {
            const double t = double(i) / double(steps);
            const double x = x0 + t * dx;
            const double y = y0 + t * dy;
            draw_disc((int)std::round(x), (int)std::round(y), radius, c, a);
        }
    }

    void draw_dashed_line(
        double x0,
        double y0,
        double x1,
        double y1,
        RGB c,
        double a,
        double width_px,
        double dash_len,
        double gap_len)
    {
        const double dx = x1 - x0;
        const double dy = y1 - y0;
        const double len = std::sqrt(dx * dx + dy * dy);

        if (len <= 0.0) {
            return;
        }

        double pos = 0.0;

        while (pos < len) {
            const double a0 = pos / len;
            const double a1 = std::min(len, pos + dash_len) / len;

            draw_line(
                x0 + a0 * dx,
                y0 + a0 * dy,
                x0 + a1 * dx,
                y0 + a1 * dy,
                c,
                a,
                width_px);

            pos += dash_len + gap_len;
        }
    }

    void draw_polyline(const std::vector<double>& xs,
        const std::vector<double>& ys,
        RGB c,
        double a,
        double width_px,
        bool dashed)
    {
        if (xs.size() < 2 || xs.size() != ys.size()) {
            return;
        }

        for (int i = 1; i < (int)xs.size(); ++i) {
            if (dashed) {
                draw_dashed_line(xs[i - 1], ys[i - 1], xs[i], ys[i], c, a, width_px, 7.0, 5.0);
            }
            else {
                draw_line(xs[i - 1], ys[i - 1], xs[i], ys[i], c, a, width_px);
            }
        }
    }

    void draw_disc(int cx, int cy, double r, RGB c, double a) {
        const int rr = (int)std::ceil(r);

        for (int y = cy - rr; y <= cy + rr; ++y) {
            for (int x = cx - rr; x <= cx + rr; ++x) {
                const double dx = double(x) - double(cx);
                const double dy = double(y) - double(cy);

                if (dx * dx + dy * dy <= r * r) {
                    blend_pixel(x, y, c, a);
                }
            }
        }
    }

    void fill_polygon(const std::vector<double>& xs, const std::vector<double>& ys, RGB c, double a) {
        if (xs.size() < 3 || xs.size() != ys.size()) {
            return;
        }

        double ymin = ys[0];
        double ymax = ys[0];

        for (int i = 1; i < (int)ys.size(); ++i) {
            ymin = std::min(ymin, ys[i]);
            ymax = std::max(ymax, ys[i]);
        }

        int y0 = std::max(0, (int)std::floor(ymin));
        int y1 = std::min(height - 1, (int)std::ceil(ymax));

        std::vector<double> intersections;

        for (int y = y0; y <= y1; ++y) {
            const double yy = double(y) + 0.5;
            intersections.clear();

            for (int i = 0; i < (int)xs.size(); ++i) {
                const int j = (i + 1) % (int)xs.size();
                const double yA = ys[i];
                const double yB = ys[j];

                if ((yA <= yy && yy < yB) || (yB <= yy && yy < yA)) {
                    const double t = (yy - yA) / (yB - yA);
                    intersections.push_back(xs[i] + t * (xs[j] - xs[i]));
                }
            }

            std::sort(intersections.begin(), intersections.end());

            for (int k = 0; k + 1 < (int)intersections.size(); k += 2) {
                const int x0 = std::max(0, (int)std::floor(intersections[k]));
                const int x1 = std::min(width - 1, (int)std::ceil(intersections[k + 1]));

                for (int x = x0; x <= x1; ++x) {
                    blend_pixel(x, y, c, a);
                }
            }
        }
    }

    bool save_bmp(const char* filename) const {
        FILE* f = std::fopen(filename, "wb");

        if (!f) {
            return false;
        }

        const int row_stride = (width * 3 + 3) & ~3;
        const int image_size = row_stride * height;
        const int file_size = 54 + image_size;

        unsigned char header[54];
        std::memset(header, 0, sizeof(header));

        header[0] = 'B';
        header[1] = 'M';
        write_le32(header + 2, file_size);
        write_le32(header + 10, 54);
        write_le32(header + 14, 40);
        write_le32(header + 18, width);
        write_le32(header + 22, height);
        write_le16(header + 26, 1);
        write_le16(header + 28, 24);
        write_le32(header + 34, image_size);

        std::fwrite(header, 1, 54, f);

        std::vector<unsigned char> row(row_stride, 255);

        for (int y = height - 1; y >= 0; --y) {
            std::memset(row.data(), 255, row.size());

            for (int x = 0; x < width; ++x) {
                const size_t src = ((size_t)y * (size_t)width + (size_t)x) * 3;
                const size_t dst = (size_t)x * 3;

                row[dst + 0] = pixels[src + 2];
                row[dst + 1] = pixels[src + 1];
                row[dst + 2] = pixels[src + 0];
            }

            std::fwrite(row.data(), 1, row.size(), f);
        }

        std::fclose(f);
        return true;
    }

private:
    int width;
    int height;
    std::vector<unsigned char> pixels;

    static void write_le16(unsigned char* p, int v) {
        p[0] = (unsigned char)(v & 255);
        p[1] = (unsigned char)((v >> 8) & 255);
    }

    static void write_le32(unsigned char* p, int v) {
        p[0] = (unsigned char)(v & 255);
        p[1] = (unsigned char)((v >> 8) & 255);
        p[2] = (unsigned char)((v >> 16) & 255);
        p[3] = (unsigned char)((v >> 24) & 255);
    }
};

class PlotGeometry {
public:
    PlotGeometry(const std::vector<DiscrepancyCurve>& curves_,
        const std::vector<ReferenceCurveSpec>& refs_,
        const std::vector<PlotPowerGuideSpec>& power_guides_,
        const DiscrepancyPlotOptions& opt_)
        : curves(curves_),
        refs(refs_),
        power_guides(power_guides_),
        opt(opt_),
        xmin(0.0),
        xmax(1.0),
        ymin(1e-6),
        ymax(1.0) {
    }

    bool prepare() {
        if (curves.empty() || opt.base < 2) {
            return false;
        }

        resolve_curves();

        bool found_data = false;

        xmin = std::numeric_limits<double>::infinity();
        xmax = -std::numeric_limits<double>::infinity();
        ymin = std::numeric_limits<double>::infinity();
        ymax = -std::numeric_limits<double>::infinity();

        for (int i = 0; i < (int)curves.size(); ++i) {
            const DiscrepancyCurve& c = curves[i];

            if (c.n_points.size() != c.values.size()) {
                return false;
            }

            for (int k = 0; k < (int)c.n_points.size(); ++k) {
                const long long N = c.n_points[k];
                const double y = c.values[k];

                if (N <= 0 || !std::isfinite(y)) {
                    continue;
                }

                const double xlog = log_base_N(N);

                xmin = std::min(xmin, xlog);
                xmax = std::max(xmax, xlog);

                update_y_range(y);
                found_data = true;

                if (curve_show_error(c, k)) {
                    const double ym = y - curve_error_minus(c, k);
                    const double yp = y + curve_error_plus(c, k);
                    update_y_range(ym);
                    update_y_range(yp);
                }
            }
        }

        if (!found_data || !std::isfinite(xmin) || !std::isfinite(xmax)) {
            return false;
        }

        if (xmax <= xmin) {
            xmax = xmin + 1.0;
        }

        resolve_references();
        include_reference_ranges();
        pad_y_range();

        plot_rect.x0 = opt.left_margin;
        plot_rect.y0 = opt.top_margin;
        plot_rect.x1 = double(opt.width) - opt.right_margin;
        plot_rect.y1 = double(opt.height) - opt.bottom_margin;

        return true;
    }

    bool write_svg(const char* filename) const {
        FILE* f = std::fopen(filename, "w");

        if (!f) {
            return false;
        }

        std::fprintf(
            f,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" "
            "width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n",
            opt.width,
            opt.height,
            opt.width,
            opt.height);

        std::fprintf(
            f,
            "  <rect x=\"0\" y=\"0\" width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n",
            opt.width,
            opt.height);

        std::fprintf(
            f,
            "  <rect x=\"%.3f\" y=\"%.3f\" width=\"%.3f\" height=\"%.3f\" "
            "fill=\"#fbfbfb\" stroke=\"#cccccc\" stroke-width=\"1\"/>\n",
            plot_rect.x0,
            plot_rect.y0,
            plot_rect.x1 - plot_rect.x0,
            plot_rect.y1 - plot_rect.y0);

        draw_svg_title(f);
        draw_svg_grid(f);
        draw_svg_power_guides(f);
        draw_svg_axes_and_ticks(f);
        draw_svg_error_fills(f);
        draw_svg_references(f);
        draw_svg_curves(f);
        draw_svg_error_bars(f);
        draw_svg_points(f);
        draw_svg_legend(f);

        std::fprintf(f, "</svg>\n");
        std::fclose(f);

        return true;
    }

    bool write_bmp(const char* filename) const {
        Bitmap bmp(opt.width, opt.height);

        draw_bmp_background(bmp);
        draw_bmp_grid(bmp);
        draw_bmp_power_guides(bmp);
        draw_bmp_error_fills(bmp);
        draw_bmp_references(bmp);
        draw_bmp_curves(bmp);
        draw_bmp_error_bars(bmp);
        draw_bmp_points(bmp);

        return bmp.save_bmp(filename);
    }

private:
    const std::vector<DiscrepancyCurve>& curves;
    const std::vector<ReferenceCurveSpec>& refs;
    const std::vector<PlotPowerGuideSpec>& power_guides;
    const DiscrepancyPlotOptions& opt;

    std::vector<ResolvedCurve> resolved_curves;
    std::vector<ResolvedReference> resolved_refs;

    double xmin;
    double xmax;
    double ymin;
    double ymax;
    PlotRect plot_rect;

    double log_base_N(long long N) const {
        return std::log(double(N)) / std::log(double(opt.base));
    }

    double x_to_pixel(double x) const {
        return plot_rect.x0 + (x - xmin) / (xmax - xmin) * (plot_rect.x1 - plot_rect.x0);
    }

    double y_to_pixel(double y) const {
        if (opt.y_scale == PLOT_Y_LOG10) {
            const double ly = safe_log10(clamp_positive(y));
            const double lymin = safe_log10(ymin);
            const double lymax = safe_log10(ymax);
            return plot_rect.y1 - (ly - lymin) / (lymax - lymin) * (plot_rect.y1 - plot_rect.y0);
        }

        return plot_rect.y1 - (y - ymin) / (ymax - ymin) * (plot_rect.y1 - plot_rect.y0);
    }

    void update_y_range(double y) {
        if (!std::isfinite(y)) {
            return;
        }

        if (opt.y_scale == PLOT_Y_LOG10 && y <= 0.0) {
            return;
        }

        ymin = std::min(ymin, y);
        ymax = std::max(ymax, y);
    }

    void resolve_curves() {
        std::map<std::string, std::string> group_color;

        resolved_curves.clear();
        resolved_curves.reserve(curves.size());

        int next_color = 0;

        for (int i = 0; i < (int)curves.size(); ++i) {
            const DiscrepancyCurve& c = curves[i];

            const std::string key = group_key_for_curve(c, i);

            if (group_color.find(key) == group_color.end()) {
                if (!c.color.empty()) {
                    group_color[key] = c.color;
                }
                else {
                    group_color[key] = opt.default_palette[next_color % opt.default_palette.size()];
                }

                ++next_color;
            }

            ResolvedCurve r;
            r.curve = &c;
            r.group_key = key;
            r.legend_label = !c.legend_group.empty() ? c.legend_group : c.label;
            r.color = !c.color.empty() ? c.color : group_color[key];
            r.stroke_width = c.stroke_width > 0.0 ? c.stroke_width : opt.default_curve_width;
            r.opacity = c.opacity;
            r.point_radius = c.point_radius;
            r.point_opacity = c.point_opacity;

            resolved_curves.push_back(r);
        }
    }

    void resolve_references() {
        resolved_refs.clear();
        resolved_refs.reserve(refs.size());

        double anchor_N = -1.0;
        double anchor_y = -1.0;

        for (int i = 0; i < (int)curves.size() && anchor_N <= 0.0; ++i) {
            const DiscrepancyCurve& c = curves[i];

            for (int k = 0; k < (int)c.n_points.size(); ++k) {
                if (c.n_points[k] > 0 && c.values[k] > 0.0) {
                    anchor_N = double(c.n_points[k]);
                    anchor_y = c.values[k];
                    break;
                }
            }
        }

        if (anchor_N <= 0.0 || anchor_y <= 0.0) {
            anchor_N = std::pow(double(opt.base), xmin);
            anchor_y = 1.0;
        }

        for (int i = 0; i < (int)refs.size(); ++i) {
            const ReferenceCurveSpec& ref = refs[i];

            ResolvedReference rr;
            rr.ref = &ref;
            rr.color = !ref.color.empty() ? ref.color : default_ref_color(i);
            rr.stroke_width = ref.stroke_width > 0.0 ? ref.stroke_width : opt.default_ref_width;

            if (ref.scale > 0.0) {
                rr.scale = ref.scale;
            }
            else {
                rr.scale = anchor_y * std::pow(anchor_N, ref.exponent);
            }

            resolved_refs.push_back(rr);
        }
    }

    void include_reference_ranges() {
        const int nref_samples = 128;

        for (int i = 0; i < (int)resolved_refs.size(); ++i) {
            const ResolvedReference& r = resolved_refs[i];

            for (int t = 0; t < nref_samples; ++t) {
                const double alpha = double(t) / double(nref_samples - 1);
                const double xlog = xmin + alpha * (xmax - xmin);
                const double N = std::pow(double(opt.base), xlog);
                const double y = r.scale * std::pow(N, -r.ref->exponent);

                update_y_range(y);
            }
        }
    }

    void pad_y_range() {
        if (opt.y_scale == PLOT_Y_LOG10) {
            ymin = clamp_positive(ymin);
            ymax = clamp_positive(ymax);

            double lymin = safe_log10(ymin);
            double lymax = safe_log10(ymax);

            if (lymax <= lymin) {
                lymax = lymin + 1.0;
            }

            const double pad = 0.08 * (lymax - lymin);
            lymin -= pad;
            lymax += pad;

            ymin = std::pow(10.0, lymin);
            ymax = std::pow(10.0, lymax);
        }
        else {
            if (ymax <= ymin) {
                ymax = ymin + 1.0;
            }

            const double pad = 0.08 * (ymax - ymin);
            ymin -= pad;
            ymax += pad;
        }
    }

    void draw_svg_title(FILE* f) const {
        if (!opt.title.empty()) {
            svg_text(
                f,
                opt.width * 0.5,
                28.0,
                opt.title,
                opt.title_font_size,
                "middle",
                "#222222");
        }
    }

    void draw_svg_grid(FILE* f) const {
        if (!opt.draw_grid) {
            return;
        }

        if (opt.draw_grid_vertical) {
            for_each_x_tick([&](double xlog, bool) {
                const double x = x_to_pixel(xlog);
                svg_line(f, x, plot_rect.y0, x, plot_rect.y1, "#e0e0e0", 1.0, 1.0);
                });
        }

        if (opt.draw_grid_horizontal) {
            for_each_y_tick([&](double yvalue) {
                const double y = y_to_pixel(yvalue);
                svg_line(f, plot_rect.x0, y, plot_rect.x1, y, "#e0e0e0", 1.0, 1.0);
                });
        }
    }

    void draw_svg_power_guides(FILE* f) const {
        for (int i = 0; i < (int)power_guides.size(); ++i) {
            const PlotPowerGuideSpec& g = power_guides[i];

            if (g.base <= 1) {
                continue;
            }

            const int e0 = guide_min_exp(g);
            const int e1 = guide_max_exp(g);

            const std::string extra = g.dashed
                ? std::string("stroke-dasharray=\"") + g.dash_array + "\""
                : std::string();

            for (int e = e0; e <= e1; ++e) {
                const double N = std::pow(double(g.base), double(e));
                const double xlog = std::log(N) / std::log(double(opt.base));

                if (xlog < xmin - 1e-12 || xlog > xmax + 1e-12) {
                    continue;
                }

                const double x = x_to_pixel(xlog);
                svg_line(f, x, plot_rect.y0, x, plot_rect.y1, g.color, g.stroke_width, g.opacity, extra);
            }
        }
    }

    void draw_svg_axes_and_ticks(FILE* f) const {
        svg_line(f, plot_rect.x0, plot_rect.y1, plot_rect.x1, plot_rect.y1, "#333333", 1.4, 1.0);
        svg_line(f, plot_rect.x0, plot_rect.y0, plot_rect.x0, plot_rect.y1, "#333333", 1.4, 1.0);

        for_each_x_tick([&](double xlog, bool draw_label) {
            const double x = x_to_pixel(xlog);

            svg_line(f, x, plot_rect.y1, x, plot_rect.y1 + 5.0, "#333333", 1.0, 1.0);

            if (!draw_label) {
                return;
            }

            std::string label;

            if (opt.x_tick_label_as_counts) {
                const double N = std::pow(double(opt.base), xlog);
                label = format_number(N, opt.x_tick_format);
            }
            else {
                label = format_number(xlog, PLOT_TICK_LOG_VALUE);
            }

            svg_text(
                f,
                x,
                plot_rect.y1 + 22.0,
                label,
                opt.x_tick_font_size,
                std::fabs(opt.x_tick_rotation_degrees) > 1e-12 ? "end" : "middle",
                "#333333",
                opt.x_tick_rotation_degrees);
            });

        for_each_y_tick([&](double yvalue) {
            const double y = y_to_pixel(yvalue);

            svg_line(f, plot_rect.x0 - 5.0, y, plot_rect.x0, y, "#333333", 1.0, 1.0);

            const std::string label = format_number(yvalue, opt.y_tick_format);
            svg_text(f, plot_rect.x0 - 10.0, y + 4.0, label, opt.y_tick_font_size, "end", "#333333");
            });

        if (!opt.x_label.empty()) {
            svg_text(
                f,
                0.5 * (plot_rect.x0 + plot_rect.x1),
                opt.height - 28.0,
                opt.x_label,
                opt.axis_font_size,
                "middle",
                "#222222");
        }

        if (!opt.y_label.empty()) {
            std::fprintf(
                f,
                "  <text x=\"24\" y=\"%.3f\" font-family=\"Arial, Helvetica, sans-serif\" "
                "font-size=\"%d\" text-anchor=\"middle\" fill=\"#222222\" "
                "transform=\"rotate(-90 24 %.3f)\">%s</text>\n",
                0.5 * (plot_rect.y0 + plot_rect.y1),
                opt.axis_font_size,
                0.5 * (plot_rect.y0 + plot_rect.y1),
                opt.y_label.c_str());
        }
    }

    void draw_svg_error_fills(FILE* f) const {
        for (int i = 0; i < (int)resolved_curves.size(); ++i) {
            const ResolvedCurve& rc = resolved_curves[i];
            const DiscrepancyCurve& c = *rc.curve;

            if (!(c.error_style == PLOT_ERROR_FILL || c.error_style == PLOT_ERROR_BARS_AND_FILL)) {
                continue;
            }

            if (!has_symmetric_or_asymmetric_errors(c)) {
                continue;
            }

            std::vector<std::pair<double, double> > upper;
            std::vector<std::pair<double, double> > lower;

            for (int k = 0; k < (int)c.n_points.size(); ++k) {
                if (c.n_points[k] <= 0 || !curve_show_error(c, k)) {
                    continue;
                }

                const double x = x_to_pixel(log_base_N(c.n_points[k]));
                const double ym = c.values[k] - curve_error_minus(c, k);
                const double yp = c.values[k] + curve_error_plus(c, k);

                if (opt.y_scale == PLOT_Y_LOG10 && ym <= 0.0) {
                    continue;
                }

                upper.push_back(std::make_pair(x, y_to_pixel(yp)));
                lower.push_back(std::make_pair(x, y_to_pixel(ym)));
            }

            if (upper.size() < 2 || upper.size() != lower.size()) {
                continue;
            }

            std::fprintf(
                f,
                "  <polygon fill=\"%s\" fill-opacity=\"%.3f\" stroke=\"none\" points=\"",
                rc.color.c_str(),
                c.error_fill_opacity);

            for (int k = 0; k < (int)upper.size(); ++k) {
                std::fprintf(f, "%.3f,%.3f ", upper[k].first, upper[k].second);
            }

            for (int k = (int)lower.size() - 1; k >= 0; --k) {
                std::fprintf(f, "%.3f,%.3f ", lower[k].first, lower[k].second);
            }

            std::fprintf(f, "\"/>\n");
        }
    }

    void draw_svg_references(FILE* f) const {
        for (int i = 0; i < (int)resolved_refs.size(); ++i) {
            const ResolvedReference& r = resolved_refs[i];
            const int nref_samples = 128;

            std::fprintf(
                f,
                "  <polyline fill=\"none\" stroke=\"%s\" stroke-width=\"%.3f\" "
                "stroke-opacity=\"%.3f\" ",
                r.color.c_str(),
                r.stroke_width,
                r.ref->opacity);

            if (r.ref->dashed) {
                std::fprintf(f, "stroke-dasharray=\"%s\" ", r.ref->dash_array.c_str());
            }

            std::fprintf(f, "stroke-linejoin=\"round\" stroke-linecap=\"round\" points=\"");

            for (int t = 0; t < nref_samples; ++t) {
                const double alpha = double(t) / double(nref_samples - 1);
                const double xlog = xmin + alpha * (xmax - xmin);
                const double N = std::pow(double(opt.base), xlog);
                const double y = r.scale * std::pow(N, -r.ref->exponent);

                std::fprintf(f, "%.3f,%.3f ", x_to_pixel(xlog), y_to_pixel(y));
            }

            std::fprintf(f, "\"/>\n");
        }
    }

    void draw_svg_curves(FILE* f) const {
        for (int i = 0; i < (int)resolved_curves.size(); ++i) {
            const ResolvedCurve& rc = resolved_curves[i];
            const DiscrepancyCurve& c = *rc.curve;

            if (c.n_points.empty()) {
                continue;
            }

            std::fprintf(
                f,
                "  <polyline fill=\"none\" stroke=\"%s\" stroke-width=\"%.3f\" "
                "stroke-opacity=\"%.3f\" ",
                rc.color.c_str(),
                rc.stroke_width,
                rc.opacity);

            if (c.dashed) {
                std::fprintf(f, "stroke-dasharray=\"%s\" ", c.dash_array.c_str());
            }

            std::fprintf(f, "stroke-linejoin=\"round\" stroke-linecap=\"round\" points=\"");

            for (int k = 0; k < (int)c.n_points.size(); ++k) {
                if (c.n_points[k] <= 0) {
                    continue;
                }

                const double xlog = log_base_N(c.n_points[k]);
                const double y = c.values[k];

                if (opt.y_scale == PLOT_Y_LOG10 && y <= 0.0) {
                    continue;
                }

                std::fprintf(f, "%.3f,%.3f ", x_to_pixel(xlog), y_to_pixel(y));
            }

            std::fprintf(f, "\"/>\n");
        }
    }

    void draw_svg_error_bars(FILE* f) const {
        for (int i = 0; i < (int)resolved_curves.size(); ++i) {
            const ResolvedCurve& rc = resolved_curves[i];
            const DiscrepancyCurve& c = *rc.curve;

            if (!(c.error_style == PLOT_ERROR_BARS || c.error_style == PLOT_ERROR_BARS_AND_FILL)) {
                continue;
            }

            for (int k = 0; k < (int)c.n_points.size(); ++k) {
                if (c.n_points[k] <= 0 || !curve_show_error(c, k)) {
                    continue;
                }

                const double x = x_to_pixel(log_base_N(c.n_points[k]));
                const double ym = c.values[k] - curve_error_minus(c, k);
                const double yp = c.values[k] + curve_error_plus(c, k);

                if (opt.y_scale == PLOT_Y_LOG10 && ym <= 0.0) {
                    continue;
                }

                const double y0 = y_to_pixel(ym);
                const double y1 = y_to_pixel(yp);
                const double cap = 0.5 * c.error_bar_cap_width_pixels;

                svg_line(f, x, y0, x, y1, rc.color, c.error_bar_width_pixels, rc.opacity);
                svg_line(f, x - cap, y0, x + cap, y0, rc.color, c.error_bar_width_pixels, rc.opacity);
                svg_line(f, x - cap, y1, x + cap, y1, rc.color, c.error_bar_width_pixels, rc.opacity);
            }
        }
    }

    void draw_svg_points(FILE* f) const {
        for (int i = 0; i < (int)resolved_curves.size(); ++i) {
            const ResolvedCurve& rc = resolved_curves[i];
            const DiscrepancyCurve& c = *rc.curve;

            for (int k = 0; k < (int)c.n_points.size(); ++k) {
                if (c.n_points[k] <= 0 || !curve_show_point(c, k)) {
                    continue;
                }

                const double xlog = log_base_N(c.n_points[k]);
                const double y = c.values[k];

                if (opt.y_scale == PLOT_Y_LOG10 && y <= 0.0) {
                    continue;
                }

                std::fprintf(
                    f,
                    "  <circle cx=\"%.3f\" cy=\"%.3f\" r=\"%.3f\" "
                    "fill=\"%s\" fill-opacity=\"%.3f\" stroke=\"#2f2f2f\" "
                    "stroke-opacity=\"%.3f\" stroke-width=\"1\">\n",
                    x_to_pixel(xlog),
                    y_to_pixel(y),
                    rc.point_radius,
                    rc.color.c_str(),
                    rc.point_opacity,
                    std::min(1.0, rc.opacity + 0.1));

                std::fprintf(
                    f,
                    "    <title>N=%lld, value=%.12g</title>\n",
                    c.n_points[k],
                    c.values[k]);

                std::fprintf(f, "  </circle>\n");
            }
        }
    }

    void draw_svg_legend(FILE* f) const {
        if (!opt.draw_legend) {
            return;
        }

        const std::vector<LegendEntry> entries = build_legend_entries();

        if (entries.empty()) {
            return;
        }

        const double row_h = 19.0;
        //const double legend_w = 300.0;

        double max_label_w = 0.0;

        for (int i = 0; i < (int)entries.size(); ++i) {
            max_label_w = std::max(
                max_label_w,
                estimate_svg_text_width(entries[i].label, opt.legend_font_size)
            );
        }

        // 58 px = left padding + line sample + gap before text.
        // 18 px = right padding.
        // Keep a minimum width for short legends.
        const double legend_w = std::max(180.0, 58.0 + max_label_w + 18.0);

        const double legend_h = 16.0 + row_h * entries.size();

        double lx = plot_rect.x1 - legend_w - 12.0;
        double ly = plot_rect.y0 + 12.0;

        if (opt.legend_position == PLOT_LEGEND_TOP_LEFT) {
            lx = plot_rect.x0 + 12.0;
            ly = plot_rect.y0 + 12.0;
        }
        else if (opt.legend_position == PLOT_LEGEND_BOTTOM_RIGHT) {
            lx = plot_rect.x1 - legend_w - 12.0;
            ly = plot_rect.y1 - legend_h - 12.0;
        }
        else if (opt.legend_position == PLOT_LEGEND_BOTTOM_LEFT) {
            lx = plot_rect.x0 + 12.0;
            ly = plot_rect.y1 - legend_h - 12.0;
        }

        std::fprintf(
            f,
            "  <rect x=\"%.3f\" y=\"%.3f\" width=\"%.3f\" height=\"%.3f\" "
            "rx=\"6\" ry=\"6\" fill=\"#ffffff\" fill-opacity=\"0.92\" stroke=\"#d0d0d0\"/>\n",
            lx,
            ly,
            legend_w,
            legend_h);

        for (int i = 0; i < (int)entries.size(); ++i) {
            const LegendEntry& e = entries[i];
            const double yy = ly + 16.0 + i * row_h;

            const std::string extra = e.dashed
                ? std::string("stroke-dasharray=\"") + e.dash_array + "\""
                : std::string();

            svg_line(f, lx + 12.0, yy, lx + 48.0, yy, e.color, e.width, e.opacity, extra);
            if (e.show_point)
                svg_circle(f, lx + 30.0, yy, 3.0, e.color, "#2f2f2f", std::min(1.0, e.opacity + 0.1), 0.75);

            svg_text(
                f,
                lx + 58.0,
                yy + 4.0,
                e.label,
                opt.legend_font_size,
                "start",
                "#333333");
        }
    }

    std::vector<LegendEntry> build_legend_entries() const {
        std::vector<LegendEntry> entries;
        std::vector<std::string> keys;

        for (int i = 0; i < (int)resolved_curves.size(); ++i) {
            const ResolvedCurve& rc = resolved_curves[i];
            const DiscrepancyCurve& c = *rc.curve;

            if (!c.show_in_legend) {
                continue;
            }

            if (std::find(keys.begin(), keys.end(), rc.group_key) != keys.end()) {
                continue;
            }

            keys.push_back(rc.group_key);

            LegendEntry e;
            e.label = !rc.legend_label.empty() ? rc.legend_label : rc.group_key;
            e.color = rc.color;
            e.show_point = c.show_points;
            e.width = rc.stroke_width;
            e.opacity = rc.opacity;
            e.dashed = c.dashed;
            e.dash_array = c.dash_array;

            entries.push_back(e);
        }

        for (int i = 0; i < (int)resolved_refs.size(); ++i) {
            const ResolvedReference& rr = resolved_refs[i];

            if (!rr.ref->show_in_legend) {
                continue;
            }

            LegendEntry e;
            e.label = rr.ref->label.empty() ? "Reference" : rr.ref->label;
            e.color = rr.color;
            e.width = rr.stroke_width;
            e.opacity = rr.ref->opacity;
            e.dashed = rr.ref->dashed;
            e.dash_array = rr.ref->dash_array;

            entries.push_back(e);
        }

        for (int i = 0; i < (int)power_guides.size(); ++i) {
            const PlotPowerGuideSpec& g = power_guides[i];

            if (!g.show_in_legend) {
                continue;
            }

            LegendEntry e;
            e.label = !g.label.empty() ? g.label : std::string("powers of ") + format_number(g.base, PLOT_TICK_INTEGER);
            e.color = g.color;
            e.width = g.stroke_width;
            e.opacity = g.opacity;
            e.dashed = g.dashed;
            e.dash_array = g.dash_array;

            entries.push_back(e);
        }

        return entries;
    }

    template<class Fn>
    void for_each_x_tick(Fn fn) const {
        const double step = opt.x_tick_octave_step > 0.0 ? opt.x_tick_octave_step : 1.0;
        const double first = std::ceil(xmin / step) * step;

        for (double xt = first; xt <= xmax + 1e-12; xt += step) {
            bool draw_label = true;

            if (opt.x_tick_label_only_integer_octaves && !nearly_integer(xt)) {
                draw_label = false;
            }

            fn(xt, draw_label);
        }
    }

    template<class Fn>
    void for_each_y_tick(Fn fn) const {
        if (opt.y_scale == PLOT_Y_LOG10) {
            const double lymin = safe_log10(ymin);
            const double lymax = safe_log10(ymax);

            const int step = std::max(1, opt.y_log10_tick_step);
            const int emin = (int)std::ceil(lymin);
            const int emax = (int)std::floor(lymax);

            for (int e = emin; e <= emax; e += step) {
                fn(std::pow(10.0, double(e)));
            }
        }
        else {
            if (opt.y_linear_tick_step <= 0.0) {
                return;
            }

            const double first = std::ceil(ymin / opt.y_linear_tick_step) * opt.y_linear_tick_step;

            for (double yt = first; yt <= ymax + 1e-12; yt += opt.y_linear_tick_step) {
                fn(yt);
            }
        }
    }

    int guide_min_exp(const PlotPowerGuideSpec& g) const {
        if (g.min_exponent != std::numeric_limits<int>::min()) {
            return g.min_exponent;
        }

        const double Nmin = std::pow(double(opt.base), xmin);
        return (int)std::ceil(std::log(Nmin) / std::log(double(g.base)) - 1e-12);
    }

    int guide_max_exp(const PlotPowerGuideSpec& g) const {
        if (g.max_exponent >= g.min_exponent) {
            return g.max_exponent;
        }

        const double Nmax = std::pow(double(opt.base), xmax);
        return (int)std::floor(std::log(Nmax) / std::log(double(g.base)) + 1e-12);
    }

    void get_curve_pixels(const DiscrepancyCurve& c,
        std::vector<double>& xs,
        std::vector<double>& ys) const
    {
        xs.clear();
        ys.clear();

        for (int k = 0; k < (int)c.n_points.size(); ++k) {
            if (c.n_points[k] <= 0) {
                continue;
            }

            const double y = c.values[k];

            if (opt.y_scale == PLOT_Y_LOG10 && y <= 0.0) {
                continue;
            }

            xs.push_back(x_to_pixel(log_base_N(c.n_points[k])));
            ys.push_back(y_to_pixel(y));
        }
    }

    void draw_bmp_background(Bitmap& bmp) const {
        const RGB border = parse_color("#cccccc", RGB{ 204, 204, 204 });
        bmp.draw_line(plot_rect.x0, plot_rect.y0, plot_rect.x1, plot_rect.y0, border, 1.0, 1.0);
        bmp.draw_line(plot_rect.x1, plot_rect.y0, plot_rect.x1, plot_rect.y1, border, 1.0, 1.0);
        bmp.draw_line(plot_rect.x1, plot_rect.y1, plot_rect.x0, plot_rect.y1, border, 1.0, 1.0);
        bmp.draw_line(plot_rect.x0, plot_rect.y1, plot_rect.x0, plot_rect.y0, border, 1.0, 1.0);
    }

    void draw_bmp_grid(Bitmap& bmp) const {
        if (!opt.draw_grid) {
            return;
        }

        const RGB grid = parse_color("#e0e0e0", RGB{ 224, 224, 224 });

        if (opt.draw_grid_vertical) {
            for_each_x_tick([&](double xlog, bool) {
                const double x = x_to_pixel(xlog);
                bmp.draw_line(x, plot_rect.y0, x, plot_rect.y1, grid, 1.0, 1.0);
                });
        }

        if (opt.draw_grid_horizontal) {
            for_each_y_tick([&](double yvalue) {
                const double y = y_to_pixel(yvalue);
                bmp.draw_line(plot_rect.x0, y, plot_rect.x1, y, grid, 1.0, 1.0);
                });
        }

        const RGB axis = parse_color("#333333", RGB{ 51, 51, 51 });
        bmp.draw_line(plot_rect.x0, plot_rect.y1, plot_rect.x1, plot_rect.y1, axis, 1.0, 1.4);
        bmp.draw_line(plot_rect.x0, plot_rect.y0, plot_rect.x0, plot_rect.y1, axis, 1.0, 1.4);
    }

    void draw_bmp_power_guides(Bitmap& bmp) const {
        for (int i = 0; i < (int)power_guides.size(); ++i) {
            const PlotPowerGuideSpec& g = power_guides[i];

            if (g.base <= 1) {
                continue;
            }

            const RGB color = parse_color(g.color, RGB{ 180, 180, 180 });
            const int e0 = guide_min_exp(g);
            const int e1 = guide_max_exp(g);

            for (int e = e0; e <= e1; ++e) {
                const double N = std::pow(double(g.base), double(e));
                const double xlog = std::log(N) / std::log(double(opt.base));

                if (xlog < xmin - 1e-12 || xlog > xmax + 1e-12) {
                    continue;
                }

                const double x = x_to_pixel(xlog);

                if (g.dashed) {
                    bmp.draw_dashed_line(x, plot_rect.y0, x, plot_rect.y1, color, g.opacity, g.stroke_width, 4.0, 5.0);
                }
                else {
                    bmp.draw_line(x, plot_rect.y0, x, plot_rect.y1, color, g.opacity, g.stroke_width);
                }
            }
        }
    }

    void draw_bmp_error_fills(Bitmap& bmp) const {
        for (int i = 0; i < (int)resolved_curves.size(); ++i) {
            const ResolvedCurve& rc = resolved_curves[i];
            const DiscrepancyCurve& c = *rc.curve;

            if (!(c.error_style == PLOT_ERROR_FILL || c.error_style == PLOT_ERROR_BARS_AND_FILL)) {
                continue;
            }

            RGB color = parse_color(rc.color, RGB{ 100, 100, 100 });

            std::vector<double> xs;
            std::vector<double> ys;
            std::vector<double> lx;
            std::vector<double> ly;

            for (int k = 0; k < (int)c.n_points.size(); ++k) {
                if (c.n_points[k] <= 0 || !curve_show_error(c, k)) {
                    continue;
                }

                const double ym = c.values[k] - curve_error_minus(c, k);
                const double yp = c.values[k] + curve_error_plus(c, k);

                if (opt.y_scale == PLOT_Y_LOG10 && ym <= 0.0) {
                    continue;
                }

                xs.push_back(x_to_pixel(log_base_N(c.n_points[k])));
                ys.push_back(y_to_pixel(yp));

                lx.push_back(x_to_pixel(log_base_N(c.n_points[k])));
                ly.push_back(y_to_pixel(ym));
            }

            if (xs.size() < 2 || xs.size() != lx.size()) {
                continue;
            }

            for (int k = (int)lx.size() - 1; k >= 0; --k) {
                xs.push_back(lx[k]);
                ys.push_back(ly[k]);
            }

            bmp.fill_polygon(xs, ys, color, c.error_fill_opacity);
        }
    }

    void draw_bmp_references(Bitmap& bmp) const {
        const int nref_samples = 256;

        for (int i = 0; i < (int)resolved_refs.size(); ++i) {
            const ResolvedReference& r = resolved_refs[i];
            const RGB color = parse_color(r.color, RGB{ 100, 100, 100 });

            std::vector<double> xs;
            std::vector<double> ys;

            for (int t = 0; t < nref_samples; ++t) {
                const double alpha = double(t) / double(nref_samples - 1);
                const double xlog = xmin + alpha * (xmax - xmin);
                const double N = std::pow(double(opt.base), xlog);
                const double y = r.scale * std::pow(N, -r.ref->exponent);

                xs.push_back(x_to_pixel(xlog));
                ys.push_back(y_to_pixel(y));
            }

            bmp.draw_polyline(xs, ys, color, r.ref->opacity, r.stroke_width, r.ref->dashed);
        }
    }

    void draw_bmp_curves(Bitmap& bmp) const {
        for (int i = 0; i < (int)resolved_curves.size(); ++i) {
            const ResolvedCurve& rc = resolved_curves[i];
            const RGB color = parse_color(rc.color, RGB{ 100, 100, 100 });

            std::vector<double> xs;
            std::vector<double> ys;
            get_curve_pixels(*rc.curve, xs, ys);

            bmp.draw_polyline(xs, ys, color, rc.opacity, rc.stroke_width, rc.curve->dashed);
        }
    }

    void draw_bmp_error_bars(Bitmap& bmp) const {
        for (int i = 0; i < (int)resolved_curves.size(); ++i) {
            const ResolvedCurve& rc = resolved_curves[i];
            const DiscrepancyCurve& c = *rc.curve;
            const RGB color = parse_color(rc.color, RGB{ 100, 100, 100 });

            if (!(c.error_style == PLOT_ERROR_BARS || c.error_style == PLOT_ERROR_BARS_AND_FILL)) {
                continue;
            }

            for (int k = 0; k < (int)c.n_points.size(); ++k) {
                if (c.n_points[k] <= 0 || !curve_show_error(c, k)) {
                    continue;
                }

                const double x = x_to_pixel(log_base_N(c.n_points[k]));
                const double ym = c.values[k] - curve_error_minus(c, k);
                const double yp = c.values[k] + curve_error_plus(c, k);

                if (opt.y_scale == PLOT_Y_LOG10 && ym <= 0.0) {
                    continue;
                }

                const double y0 = y_to_pixel(ym);
                const double y1 = y_to_pixel(yp);
                const double cap = 0.5 * c.error_bar_cap_width_pixels;

                bmp.draw_line(x, y0, x, y1, color, rc.opacity, c.error_bar_width_pixels);
                bmp.draw_line(x - cap, y0, x + cap, y0, color, rc.opacity, c.error_bar_width_pixels);
                bmp.draw_line(x - cap, y1, x + cap, y1, color, rc.opacity, c.error_bar_width_pixels);
            }
        }
    }

    void draw_bmp_points(Bitmap& bmp) const {
        for (int i = 0; i < (int)resolved_curves.size(); ++i) {
            const ResolvedCurve& rc = resolved_curves[i];
            const DiscrepancyCurve& c = *rc.curve;
            const RGB color = parse_color(rc.color, RGB{ 100, 100, 100 });

            for (int k = 0; k < (int)c.n_points.size(); ++k) {
                if (c.n_points[k] <= 0 || !curve_show_point(c, k)) {
                    continue;
                }

                const double y = c.values[k];

                if (opt.y_scale == PLOT_Y_LOG10 && y <= 0.0) {
                    continue;
                }

                bmp.draw_disc(
                    (int)std::round(x_to_pixel(log_base_N(c.n_points[k]))),
                    (int)std::round(y_to_pixel(y)),
                    rc.point_radius,
                    color,
                    rc.point_opacity);
            }
        }
    }
};




ReferenceCurveSpec make_reference_n_minus_half() {
    ReferenceCurveSpec r;
    r.label = "N^{-1/2}";
    r.exponent = 0.5;
    r.dashed = true;
    r.dash_array = "7 5";
    return r;
}

ReferenceCurveSpec make_reference_n_minus_one() {
    ReferenceCurveSpec r;
    r.label = "N^{-1}";
    r.exponent = 1.0;
    r.dashed = true;
    r.dash_array = "7 5";
    return r;
}

bool plot_discrepancy_curves(
    const std::vector<DiscrepancyCurve>& curves,
    const std::vector<ReferenceCurveSpec>& refs,
    const std::vector<PlotPowerGuideSpec>& power_guides,
    const DiscrepancyPlotOptions& opt,
    const char* filename)
{
    PlotGeometry plot(curves, refs, power_guides, opt);

    if (!plot.prepare()) {
        return false;
    }

    if (opt.output_format == PLOT_OUTPUT_BMP) {
        return plot.write_bmp(filename);
    }

    return plot.write_svg(filename);
}

bool plot_discrepancy_curves_svg(
    const std::vector<DiscrepancyCurve>& curves,
    const std::vector<ReferenceCurveSpec>& refs,
    const DiscrepancyPlotOptions& opt,
    const char* filename)
{
    const std::vector<PlotPowerGuideSpec> no_power_guides;
    DiscrepancyPlotOptions svg_opt = opt;
    svg_opt.output_format = PLOT_OUTPUT_SVG;

    return plot_discrepancy_curves(curves, refs, no_power_guides, svg_opt, filename);
}
