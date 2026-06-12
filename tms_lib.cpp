#include "tms_lib.h"

#include <limits>
#include <cstdio>
#include <random>
#include <unordered_map>

#include <cstdint>
#include <cstdio>

#if defined(_MSC_VER)
#include <intrin.h>
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
    int dim) {
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
    const char* stroke) {
    std::fprintf(f,
        "<circle cx=\"%.3f\" cy=\"%.3f\" r=\"%.3f\" "
        "fill=\"%s\" stroke=\"%s\" stroke-width=\"1\"/>\n",
        cx, cy, r, fill, stroke);
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

OwenTreeND make_random_owen_tree_nd(int dim,
    int base,
    int depth,
    uint64_t seed) {
    assert(dim > 0);
    assert(base >= 2);
    assert(depth >= 0);

    FastRNG rng(seed);

    OwenTreeND tree(dim, base, depth);

    for (int d = 0; d < dim; ++d) {
        tree.trees[d] = make_random_owen_tree_1d(base, depth, rng);
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



inline double plot_clamp_positive(double x) {
    return x > 1e-300 ? x : 1e-300;
}

inline void plot_format_integer_label(long long v, char* out, int out_size) {
    if (v < 1000000LL) {
        std::snprintf(out, out_size, "%lld", v);
    }
    else {
        std::snprintf(out, out_size, "%.3g", double(v));
    }
}

inline bool plot_is_nearly_integer(double x, double eps = 1e-9) {
    return std::fabs(x - std::round(x)) <= eps;
}

inline void plot_svg_text(FILE* f,
    double x,
    double y,
    const char* text,
    int font_size = 12,
    const char* anchor = "middle",
    const char* fill = "#333333") {
    std::fprintf(f,
        "<text x=\"%.3f\" y=\"%.3f\" "
        "font-family=\"Arial, Helvetica, sans-serif\" "
        "font-size=\"%d\" text-anchor=\"%s\" fill=\"%s\">%s</text>\n",
        x, y, font_size, anchor, fill, text);
}

inline void plot_svg_line(FILE* f,
    double x1,
    double y1,
    double x2,
    double y2,
    const char* stroke,
    double width = 1.0,
    double opacity = 1.0,
    const char* extra = "") {
    std::fprintf(f,
        "<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" "
        "stroke=\"%s\" stroke-width=\"%.3f\" stroke-opacity=\"%.3f\" %s/>\n",
        x1, y1, x2, y2, stroke, width, opacity, extra);
}

inline void plot_svg_circle(FILE* f,
    double cx,
    double cy,
    double r,
    const char* fill,
    const char* stroke = "#333333",
    double fill_opacity = 1.0,
    double stroke_opacity = 1.0) {
    std::fprintf(f,
        "<circle cx=\"%.3f\" cy=\"%.3f\" r=\"%.3f\" "
        "fill=\"%s\" fill-opacity=\"%.3f\" "
        "stroke=\"%s\" stroke-opacity=\"%.3f\" stroke-width=\"1\"/>\n",
        cx, cy, r, fill, fill_opacity, stroke, stroke_opacity);
}

inline const char* plot_default_ref_color(int idx) {
    static const char* colors[4] = {
        "#D36C6C",
        "#6C8CD3",
        "#6CA87A",
        "#888888"
    };
    return colors[idx % 4];
}

inline std::string plot_resolve_group_key(const DiscrepancyCurve& c, int idx) {
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

inline int plot_find_string(const std::vector<std::string>& arr, const std::string& s) {
    for (int i = 0; i < (int)arr.size(); ++i) {
        if (arr[i] == s) {
            return i;
        }
    }
    return -1;
}

bool plot_discrepancy_curves_svg(
    const std::vector<DiscrepancyCurve>& curves,
    const std::vector<ReferenceCurveSpec>& refs,
    const DiscrepancyPlotOptions& opt,
    const char* filename)
{
    assert(filename != 0);
    assert(opt.base >= 2);

    if (curves.empty()) {
        return false;
    }

    struct ResolvedCurve {
        const DiscrepancyCurve* src;
        std::string color;
        double stroke_width;
        double opacity;
        double point_radius;
        double point_opacity;
        std::string group_key;
        std::string legend_label;
    };

    std::vector<ResolvedCurve> rc;
    rc.reserve(curves.size());

    std::vector<std::string> group_keys_seen;
    std::vector<std::string> group_colors;

    for (int i = 0; i < (int)curves.size(); ++i) {
        const DiscrepancyCurve& c = curves[i];

        if (c.n_points.size() != c.values.size()) {
            return false;
        }

        std::string group_key = plot_resolve_group_key(c, i);
        int group_idx = plot_find_string(group_keys_seen, group_key);

        if (group_idx < 0) {
            group_idx = (int)group_keys_seen.size();
            group_keys_seen.push_back(group_key);

            if (!c.color.empty()) {
                group_colors.push_back(c.color);
            }
            else {
                group_colors.push_back(opt.default_palette[group_idx % opt.default_palette.size()]);
            }
        }

        ResolvedCurve r;
        r.src = &c;
        r.group_key = group_key;
        r.legend_label = !c.legend_group.empty() ? c.legend_group : c.label;
        r.color = !c.color.empty() ? c.color : group_colors[group_idx];
        r.stroke_width = (c.stroke_width > 0.0) ? c.stroke_width : opt.default_curve_width;
        r.opacity = c.opacity;
        r.point_radius = c.point_radius;
        r.point_opacity = c.point_opacity;

        rc.push_back(r);
    }

    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();

    bool found_data = false;

    for (int i = 0; i < (int)curves.size(); ++i) {
        const DiscrepancyCurve& c = curves[i];

        for (int k = 0; k < (int)c.n_points.size(); ++k) {
            const long long N = c.n_points[k];
            const double y = c.values[k];

            if (N <= 0) {
                continue;
            }

            const double x = std::log(double(N)) / std::log(double(opt.base));

            xmin = std::min(xmin, x);
            xmax = std::max(xmax, x);

            if (opt.y_scale == PLOT_Y_LOG10) {
                if (y > 0.0) {
                    ymin = std::min(ymin, y);
                    ymax = std::max(ymax, y);
                    found_data = true;
                }
            }
            else {
                ymin = std::min(ymin, y);
                ymax = std::max(ymax, y);
                found_data = true;
            }
        }
    }

    if (!found_data || !std::isfinite(xmin) || !std::isfinite(xmax)) {
        return false;
    }

    if (xmax <= xmin) {
        xmax = xmin + 1.0;
    }

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

    struct ResolvedReference {
        const ReferenceCurveSpec* src;
        std::string color;
        double stroke_width;
        double scale;
    };

    std::vector<ResolvedReference> rr;
    rr.reserve(refs.size());

    for (int i = 0; i < (int)refs.size(); ++i) {
        const ReferenceCurveSpec& r = refs[i];

        ResolvedReference z;
        z.src = &r;
        z.color = !r.color.empty() ? r.color : plot_default_ref_color(i);
        z.stroke_width = (r.stroke_width > 0.0) ? r.stroke_width : opt.default_ref_width;

        if (r.scale > 0.0) {
            z.scale = r.scale;
        }
        else {
            z.scale = anchor_y * std::pow(anchor_N, r.exponent);
        }

        rr.push_back(z);

        const int nref_samples = 128;
        for (int t = 0; t < nref_samples; ++t) {
            const double alpha = (nref_samples == 1) ? 0.0 : double(t) / double(nref_samples - 1);
            const double x = xmin + alpha * (xmax - xmin);
            const double N = std::pow(double(opt.base), x);
            const double y = z.scale * std::pow(N, -r.exponent);

            if (opt.y_scale == PLOT_Y_LOG10) {
                if (y > 0.0) {
                    ymin = std::min(ymin, y);
                    ymax = std::max(ymax, y);
                }
            }
            else {
                ymin = std::min(ymin, y);
                ymax = std::max(ymax, y);
            }
        }
    }

    if (opt.y_scale == PLOT_Y_LOG10) {
        ymin = plot_clamp_positive(ymin);
        ymax = plot_clamp_positive(ymax);

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

    const double plot_x0 = opt.left_margin;
    const double plot_y0 = opt.top_margin;
    const double plot_w = double(opt.width) - opt.left_margin - opt.right_margin;
    const double plot_h = double(opt.height) - opt.top_margin - opt.bottom_margin;
    const double plot_x1 = plot_x0 + plot_w;
    const double plot_y1 = plot_y0 + plot_h;

    auto map_x = [&](double x) -> double {
        return plot_x0 + (x - xmin) / (xmax - xmin) * plot_w;
        };

    auto map_y = [&](double y) -> double {
        if (opt.y_scale == PLOT_Y_LOG10) {
            const double ly = safe_log10(plot_clamp_positive(y));
            const double lymin = safe_log10(ymin);
            const double lymax = safe_log10(ymax);
            return plot_y1 - (ly - lymin) / (lymax - lymin) * plot_h;
        }
        else {
            return plot_y1 - (y - ymin) / (ymax - ymin) * plot_h;
        }
        };

    FILE* f = std::fopen(filename, "w");
    if (!f) {
        return false;
    }

    std::fprintf(f,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        "version=\"1.1\" width=\"%d\" height=\"%d\" "
        "viewBox=\"0 0 %d %d\">\n",
        opt.width, opt.height, opt.width, opt.height);

    std::fprintf(f,
        "<rect x=\"0\" y=\"0\" width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n",
        opt.width, opt.height);

    std::fprintf(f,
        "<rect x=\"%.3f\" y=\"%.3f\" width=\"%.3f\" height=\"%.3f\" "
        "fill=\"#fbfbfb\" stroke=\"#cccccc\" stroke-width=\"1\"/>\n",
        plot_x0, plot_y0, plot_w, plot_h);

    if (!opt.title.empty()) {
        plot_svg_text(f,
            opt.width * 0.5,
            28.0,
            opt.title.c_str(),
            opt.title_font_size,
            "middle",
            "#222222");
    }

    if (opt.draw_grid) {
        for (double xt = std::ceil(xmin / opt.x_tick_octave_step) * opt.x_tick_octave_step;
            xt <= xmax + 1e-12;
            xt += opt.x_tick_octave_step)
        {
            const double x = map_x(xt);
            plot_svg_line(f, x, plot_y0, x, plot_y1, "#e0e0e0", 1.0, 1.0);
        }

        if (opt.y_scale == PLOT_Y_LOG10) {
            const double lymin = safe_log10(ymin);
            const double lymax = safe_log10(ymax);

            const int emin = (int)std::ceil(lymin);
            const int emax = (int)std::floor(lymax);

            for (int e = emin; e <= emax; e += std::max(1, opt.y_log10_tick_step)) {
                const double yv = std::pow(10.0, double(e));
                const double y = map_y(yv);
                plot_svg_line(f, plot_x0, y, plot_x1, y, "#e0e0e0", 1.0, 1.0);
            }
        }
        else {
            if (opt.y_linear_tick_step > 0.0) {
                const double first_tick = std::ceil(ymin / opt.y_linear_tick_step) * opt.y_linear_tick_step;
                for (double yt = first_tick; yt <= ymax + 1e-12; yt += opt.y_linear_tick_step) {
                    const double y = map_y(yt);
                    plot_svg_line(f, plot_x0, y, plot_x1, y, "#e0e0e0", 1.0, 1.0);
                }
            }
        }
    }

    plot_svg_line(f, plot_x0, plot_y1, plot_x1, plot_y1, "#333333", 1.4, 1.0);
    plot_svg_line(f, plot_x0, plot_y0, plot_x0, plot_y1, "#333333", 1.4, 1.0);

    for (double xt = std::ceil(xmin / opt.x_tick_octave_step) * opt.x_tick_octave_step;
        xt <= xmax + 1e-12;
        xt += opt.x_tick_octave_step)
    {
        const double x = map_x(xt);

        plot_svg_line(f, x, plot_y1, x, plot_y1 + 5.0, "#333333", 1.0, 1.0);

        bool draw_label = true;
        if (opt.x_tick_label_only_integer_octaves && !plot_is_nearly_integer(xt)) {
            draw_label = false;
        }

        if (draw_label) {
            char label[64];

            if (opt.x_tick_label_as_counts) {
                const long long N = plot_is_nearly_integer(xt)
                    ? integer_power_clamped(opt.base, (int)std::llround(xt))
                    : (long long)std::floor(std::pow(double(opt.base), xt) + 0.5);

                plot_format_integer_label(N, label, sizeof(label));
            }
            else {
                std::snprintf(label, sizeof(label), "%.3g", xt);
            }

            plot_svg_text(f, x, plot_y1 + 22.0, label,
                opt.tick_font_size, "middle", "#333333");
        }
    }

    if (opt.y_scale == PLOT_Y_LOG10) {
        const double lymin = safe_log10(ymin);
        const double lymax = safe_log10(ymax);

        const int emin = (int)std::ceil(lymin);
        const int emax = (int)std::floor(lymax);

        for (int e = emin; e <= emax; e += std::max(1, opt.y_log10_tick_step)) {
            const double yv = std::pow(10.0, double(e));
            const double y = map_y(yv);

            plot_svg_line(f, plot_x0 - 5.0, y, plot_x0, y, "#333333", 1.0, 1.0);

            char label[64];
            std::snprintf(label, sizeof(label), "1e%d", e);

            plot_svg_text(f, plot_x0 - 10.0, y + 4.0, label,
                opt.tick_font_size, "end", "#333333");
        }
    }
    else if (opt.y_linear_tick_step > 0.0) {
        const double first_tick = std::ceil(ymin / opt.y_linear_tick_step) * opt.y_linear_tick_step;

        for (double yt = first_tick; yt <= ymax + 1e-12; yt += opt.y_linear_tick_step) {
            const double y = map_y(yt);

            plot_svg_line(f, plot_x0 - 5.0, y, plot_x0, y, "#333333", 1.0, 1.0);

            char label[64];
            std::snprintf(label, sizeof(label), "%.3g", yt);

            plot_svg_text(f, plot_x0 - 10.0, y + 4.0, label,
                opt.tick_font_size, "end", "#333333");
        }
    }

    if (!opt.x_label.empty()) {
        plot_svg_text(f,
            plot_x0 + plot_w * 0.5,
            opt.height - 28.0,
            opt.x_label.c_str(),
            opt.axis_font_size,
            "middle",
            "#222222");
    }

    if (!opt.y_label.empty()) {
        std::fprintf(f,
            "<text x=\"24\" y=\"%.3f\" "
            "font-family=\"Arial, Helvetica, sans-serif\" "
            "font-size=\"%d\" text-anchor=\"middle\" fill=\"#222222\" "
            "transform=\"rotate(-90 24 %.3f)\">%s</text>\n",
            plot_y0 + plot_h * 0.5,
            opt.axis_font_size,
            plot_y0 + plot_h * 0.5,
            opt.y_label.c_str());
    }

    for (int i = 0; i < (int)rr.size(); ++i) {
        const ResolvedReference& r = rr[i];
        const int nref_samples = 128;

        std::fprintf(f,
            "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"%.3f\" "
            "stroke-opacity=\"%.3f\" ",
            r.color.c_str(),
            r.stroke_width,
            r.src->opacity);

        if (r.src->dashed) {
            std::fprintf(f, "stroke-dasharray=\"%s\" ",
                r.src->dash_array.c_str());
        }

        std::fprintf(f, "points=\"");

        for (int t = 0; t < nref_samples; ++t) {
            const double alpha = (nref_samples == 1) ? 0.0 : double(t) / double(nref_samples - 1);
            const double xlog = xmin + alpha * (xmax - xmin);
            const double N = std::pow(double(opt.base), xlog);
            const double y = r.scale * std::pow(N, -r.src->exponent);

            std::fprintf(f, "%.3f,%.3f ", map_x(xlog), map_y(y));
        }

        std::fprintf(f, "\"/>\n");
    }

    for (int i = 0; i < (int)rc.size(); ++i) {
        const ResolvedCurve& r = rc[i];
        const DiscrepancyCurve& c = *r.src;

        if (c.n_points.empty()) {
            continue;
        }

        std::fprintf(f,
            "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"%.3f\" "
            "stroke-opacity=\"%.3f\" ",
            r.color.c_str(),
            r.stroke_width,
            r.opacity);

        if (c.dashed) {
            std::fprintf(f, "stroke-dasharray=\"%s\" ", c.dash_array.c_str());
        }

        std::fprintf(f, "stroke-linejoin=\"round\" stroke-linecap=\"round\" points=\"");

        for (int k = 0; k < (int)c.n_points.size(); ++k) {
            if (c.n_points[k] <= 0) {
                continue;
            }

            const double xlog = std::log(double(c.n_points[k])) / std::log(double(opt.base));
            const double y = c.values[k];

            std::fprintf(f, "%.3f,%.3f ", map_x(xlog), map_y(y));
        }

        std::fprintf(f, "\"/>\n");

        if (c.show_points) {
            for (int k = 0; k < (int)c.n_points.size(); ++k) {
                if (c.n_points[k] <= 0) {
                    continue;
                }

                const double xlog = std::log(double(c.n_points[k])) / std::log(double(opt.base));
                const double y = c.values[k];

                std::fprintf(f,
                    "<circle cx=\"%.3f\" cy=\"%.3f\" r=\"%.3f\" "
                    "fill=\"%s\" fill-opacity=\"%.3f\" "
                    "stroke=\"#2f2f2f\" stroke-opacity=\"%.3f\" stroke-width=\"1\">",
                    map_x(xlog), map_y(y),
                    r.point_radius,
                    r.color.c_str(),
                    r.point_opacity,
                    std::min(1.0, r.opacity + 0.1));

                std::fprintf(f,
                    "<title>N=%lld, value=%.12g</title>",
                    c.n_points[k], c.values[k]);

                std::fprintf(f, "</circle>\n");
            }
        }
    }

    if (opt.draw_legend) {
        struct LegendEntry {
            std::string label;
            std::string color;
            double width;
            double opacity;
            bool dashed;
            std::string dash_array;
        };

        std::vector<LegendEntry> legend_entries;
        std::vector<std::string> legend_keys;

        for (int i = 0; i < (int)rc.size(); ++i) {
            const ResolvedCurve& r = rc[i];
            const DiscrepancyCurve& c = *r.src;

            if (!c.show_in_legend) {
                continue;
            }

            const std::string key = r.group_key;
            if (plot_find_string(legend_keys, key) >= 0) {
                continue;
            }

            legend_keys.push_back(key);

            LegendEntry e;
            e.label = !r.legend_label.empty() ? r.legend_label : key;
            e.color = r.color;
            e.width = r.stroke_width;
            e.opacity = r.opacity;
            e.dashed = c.dashed;
            e.dash_array = c.dash_array;

            legend_entries.push_back(e);
        }

        for (int i = 0; i < (int)rr.size(); ++i) {
            const ResolvedReference& r = rr[i];

            if (!r.src->show_in_legend) {
                continue;
            }

            LegendEntry e;
            e.label = r.src->label.empty() ? "Reference" : r.src->label;
            e.color = r.color;
            e.width = r.stroke_width;
            e.opacity = r.src->opacity;
            e.dashed = r.src->dashed;
            e.dash_array = r.src->dash_array;

            legend_entries.push_back(e);
        }

        if (!legend_entries.empty()) {
            const double row_h = 19.0;
            const double legend_w = 290.0;
            const double legend_h = 16.0 + row_h * legend_entries.size();

            const double lx = plot_x1 - legend_w - 12.0;
            const double ly = plot_y0 + 12.0;

            std::fprintf(f,
                "<rect x=\"%.3f\" y=\"%.3f\" width=\"%.3f\" height=\"%.3f\" "
                "rx=\"6\" ry=\"6\" fill=\"#ffffff\" stroke=\"#d0d0d0\"/>\n",
                lx, ly, legend_w, legend_h);

            for (int i = 0; i < (int)legend_entries.size(); ++i) {
                const LegendEntry& e = legend_entries[i];
                const double yy = ly + 16.0 + i * row_h;

                plot_svg_line(f, lx + 12.0, yy, lx + 48.0, yy,
                    e.color.c_str(), e.width, e.opacity,
                    e.dashed ? (std::string("stroke-dasharray=\"") + e.dash_array + "\"").c_str() : "");

                plot_svg_circle(f, lx + 30.0, yy, 3.0,
                    e.color.c_str(), "#2f2f2f",
                    std::min(1.0, e.opacity + 0.1),
                    std::min(1.0, e.opacity + 0.1));

                plot_svg_text(f, lx + 58.0, yy + 4.0,
                    e.label.c_str(),
                    opt.legend_font_size,
                    "start",
                    "#333333");
            }
        }
    }

    std::fprintf(f, "</svg>\n");
    std::fclose(f);

    return true;
}