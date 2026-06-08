#include "tms_lib.h"

#include <cmath>
#include <limits>
#include <cstdio>
#include <cstring>

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

double generalized_l2_discrepancy_squared(const double* points,
    int npts,
    int dim) {
    assert(points != 0);
    assert(npts > 0);
    assert(dim > 0);

    const double invN = 1.0 / double(npts);

    double term0 = 1.0;
    for (int j = 0; j < dim; ++j) {
        term0 *= 4.0 / 3.0;
    }

    double term1_sum = 0.0;

    for (int i = 0; i < npts; ++i) {
        double prod = 1.0;

        for (int j = 0; j < dim; ++j) {
            const double z = points[i * dim + j];
            prod *= (3.0 - z * z) * 0.5;
        }

        term1_sum += prod;
    }

    double term2_sum = 0.0;

    for (int i = 0; i < npts; ++i) {
        const double* pi = points + i * dim;

        for (int k = 0; k < npts; ++k) {
            const double* pk = points + k * dim;

            double prod = 1.0;

            for (int j = 0; j < dim; ++j) {
                const double m = pi[j] > pk[j] ? pi[j] : pk[j];
                prod *= 2.0 - m;
            }

            term2_sum += prod;
        }
    }

    const double d2 =
        term0
        - 2.0 * invN * term1_sum
        + invN * invN * term2_sum;

    return clamp_discrepancy_square(d2);
}

double generalized_l2_discrepancy(const double* points,
    int npts,
    int dim) {
    const double d2 =
        generalized_l2_discrepancy_squared(points, npts, dim);

    return std::sqrt(d2);
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
    }

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

        const int cleq_max =
            count4(leq_masks, a1, b1, c1, d1);

        const int cless_min =
            count4(less_masks, a0, b0, c0, d0);

        const double ub_plus =
            double(cleq_max) * invN - vol_min;

        const double ub_minus =
            vol_max - double(cless_min) * invN;

        const double ub = ub_plus > ub_minus ? ub_plus : ub_minus;

        if (ub <= best + 1e-15) {
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

        // Split along the largest remaining interval.
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

        recurse(
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