#pragma once
#include <type_traits>
#include <iostream>
#include <ostream>
#include <vector>
#include <cassert>
#include <algorithm>

#define MAX_GF 16

extern char add_non_prime[MAX_GF + 1][MAX_GF][MAX_GF];
extern char sub_non_prime[MAX_GF + 1][MAX_GF][MAX_GF];
extern char mul_non_prime[MAX_GF + 1][MAX_GF][MAX_GF];
extern char div_non_prime[MAX_GF + 1][MAX_GF][MAX_GF];
extern char neg_non_prime[MAX_GF + 1][MAX_GF];
extern char invGalois[MAX_GF + 1][MAX_GF];

double generalized_l2_discrepancy(const double* points, int npts, int dim, int block_size = 64);
double star_discrepancy(const double* points, int npts, int dim);
void print_point_range(const double* points, int npts, int dim);

template<int p, int r>
struct GFCardinality {
    enum { value = p * GFCardinality<p, r - 1>::value };
};

template<int p>
struct GFCardinality<p, 0> {
    enum { value = 1 };
};

template<int p, int r>
struct GF_ext {
    int v;
};

template<int p>
using GF_prime = int;

template<int p, int r>
using GF = std::conditional_t<r == 1, int, GF_ext<p, r>>;


template<int p, int r>
inline GF_ext<p, r> operator+(GF_ext<p, r> a, GF_ext<p, r> b) {
    enum { Q = GFCardinality<p, r>::value };
    return GF_ext<p, r>{ static_cast<unsigned char>(add_non_prime[Q][a.v][b.v]) };
}
template<int p, int r>
inline GF_ext<p, r>& operator+=(GF_ext<p, r>& a, GF_ext<p, r> b) {
    enum { Q = GFCardinality<p, r>::value };

    a.v = static_cast<unsigned char>(add_non_prime[Q][a.v][b.v]);
    return a;
}
template<int p, int r>
inline GF_ext<p, r> operator-(GF_ext<p, r> a, GF_ext<p, r> b) {
    enum { Q = GFCardinality<p, r>::value };
    return GF_ext<p, r>{ static_cast<unsigned char>(sub_non_prime[Q][a.v][b.v]) };
}
template<int p, int r>
inline GF_ext<p, r>& operator-=(GF_ext<p, r>& a, GF_ext<p, r> b) {
    enum { Q = GFCardinality<p, r>::value };

    a.v = static_cast<unsigned char>(sub_non_prime[Q][a.v][b.v]);
    return a;
}
template<int p, int r>
inline GF_ext<p, r> operator*(GF_ext<p, r> a, GF_ext<p, r> b) {
    enum { Q = GFCardinality<p, r>::value };
    return GF_ext<p, r>{ static_cast<unsigned char>(mul_non_prime[Q][a.v][b.v]) };
}
template<int p, int r>
inline GF_ext<p, r>& operator*=(GF_ext<p, r>& a, GF_ext<p, r> b) {
    enum { Q = GFCardinality<p, r>::value };

    a.v = static_cast<unsigned char>(mul_non_prime[Q][a.v][b.v]);
    return a;
}


template<int p, int r>
inline GF_ext<p, r> operator/(GF_ext<p, r> a, GF_ext<p, r> b) {
    enum { Q = GFCardinality<p, r>::value };
    return GF_ext<p, r>{ static_cast<unsigned char>(div_non_prime[Q][a.v][b.v]) };
}
template<int p, int r>
inline GF_ext<p, r>& operator/=(GF_ext<p, r>& a, GF_ext<p, r> b) {
    enum { Q = GFCardinality<p, r>::value };

    a.v = static_cast<unsigned char>(div_non_prime[Q][a.v][b.v]);
    return a;
}

template<int p, int r>
inline GF_ext<p, r> operator-(GF_ext<p, r> a) {
    enum { Q = GFCardinality<p, r>::value };
    return GF_ext<p, r>{ static_cast<unsigned char>(neg_non_prime[Q][a.v]) };
}

template<int p, int r>
inline bool operator==(GF_ext<p, r> a, GF_ext<p, r> b) {
    return a.v == b.v;
}

template<int p, int r>
inline bool operator!=(GF_ext<p, r> a, GF_ext<p, r> b) {
    return a.v != b.v;
}

template<int p, int r>
std::ostream& operator<<(std::ostream& out, GF_ext<p, r> a) {
    return out << a.v;
}


template<int p, int r>
inline typename std::enable_if<r == 1, int>::type
gf_div(int a, int b) {
    return (a * static_cast<unsigned char>(invGalois[p][b])) % p;
}

template<int p, int r>
inline typename std::enable_if<r != 1, GF_ext<p, r>>::type
gf_div(GF_ext<p, r> a, GF_ext<p, r> b) {
    enum { Q = GFCardinality<p, r>::value };
    return GF_ext<p, r>{ static_cast<unsigned char>(div_non_prime[Q][a.v][b.v]) };
}

template<int p, int r, typename T>
inline typename std::enable_if<r == 1, T>::type
gf_reduce(T x) {
    T m = x % p;
    return m>=0 ? m:m+p;
}

template<int p, int r, typename T>
inline typename std::enable_if<r != 1, T>::type
gf_reduce(T x) {
    return x;
}

template<int p_, int r_>
struct Field {
    enum { p = p_, r = r_ };
    typedef GF<p, r> T;

    static inline T zero() {
        return T{ 0 };
    }

    static inline T one() {
        return T{ 1 };
    }

    static inline T add(T a, T b) { return a+b; }
    static inline T sub(T a, T b) { return a - b; }
    static inline T mul(T a, T b) { return a * b; }
    static inline T div(T a, T b) { return gf_div<p, r>(a, b); }
    static inline T neg(T a) { return -a; }
};
template<int p, int r>
inline int positive_mod_p(int x) {
    int y = x % p;
    return y < 0 ? y + p : y;
}
template<int p, int r>
inline typename std::enable_if<r == 1, int>::type
gf_from_int(int x) {
    return positive_mod_p<p, r>(x);
}
template<int p, int r>
inline typename std::enable_if<r != 1, GF_ext<p, r>>::type
gf_from_int(int x) {
    return GF_ext<p, r>{ positive_mod_p<p, r>(x) };
}

template<class F>
typename std::enable_if<F::r == 1, int>::type
gf_to_raw_index(typename F::T x) {
    return x;
}

template<class F>
typename std::enable_if<F::r != 1, int>::type
gf_to_raw_index(typename F::T x) {
    return x.v;
}

// assumes digits has a sufficiently large allocated size ! typically std::floor(std::log(max_val) / std::log(Q))+1;
template<class F>
void decompose_integer_into_base(long long val, typename F::T* digits, int max_digits) {
    typedef typename F::T T;
    enum { Q = GFCardinality<F::p, F::r>::value };

    assert(Q >= 2);
    assert(val >= 0);

    for (int i = 0; i < max_digits; ++i) {
        digits[i] = T{ static_cast<int>(val % Q) };
        val /= Q;
    }
}



template<class F>
struct MatrixView;

template<class F>
struct Matrix;

template<class Derived, class F>
struct MatrixBase {
    typedef typename F::T T;

    Derived& derived() {
        return static_cast<Derived&>(*this);
    }

    const Derived& derived() const {
        return static_cast<const Derived&>(*this);
    }

    T& operator()(int i, int j) {
        return derived().values[i * derived().n + j];
    }

    const T& operator()(int i, int j) const {
        return derived().values[i * derived().n + j];
    }

    T* data() {
        return derived().values;
    }

    const T* data() const {
        return derived().values;
    }

    int rows() const {
        return derived().m;
    }

    int cols() const {
        return derived().n;
    }

    void set_zero() {
        memset(data(), 0, rows() * cols() * sizeof(T));
    }
    void set_id() {
        int n = cols();
        memset(data(), 0, rows() * cols() * sizeof(T));
        for (int i = 0; i < n; i++) {
            data()[i * n + i] = T{ 1 };
        }
    }

    T operator[](int i) const {
        return derived().values[i];
    }
    T& operator[](int i) {
        return derived().values[i];
    }

    // test dets of square submatrices anchored *at the top left*, and return the row when det==0 mod GF. Returns -1 if invertible
    // tmpMatT and tmpInv should be allocated of size m x m, tmpiAc of size m x 1, and tmpriA of size 1 x m
    // first 5x5 dets in closed form ; larger sizes with Sherman-Morrison formula/det lemma
    // assumes no pivot is needed (typical from Sobol' matrices and their combination) -- eg won't work for an antidiagonal matrix
    int checkDet(MatrixView<F> tmpMatT, MatrixView<F> tmpInv, MatrixView<F> tmpiAc, MatrixView<F> tmpriA) const;

    // same but slower : need to allocate memory instead of using preallocated tmp storage
    int checkDet() const;

    // assumes no pivot is needed (typical from Sobol' matrices and their combination) -- eg won't work for an antidiagonal matrix
    Matrix<F> get_inverse(MatrixView<F> tmpMatT, MatrixView<F> tmpiAc, MatrixView<F> tmpriA) const;

    // same but slower : need to allocate memory instead of using preallocated tmp storage
    Matrix<F> get_inverse() const;

    // in case the matrix is missing modulos
    void reduce() {
        for (int i = 0; i < rows() * cols(); i++) {
            data()[i] = gf_reduce<F::p, F::r>(data()[i]);
        }
    }

    // rank of the matrix : includes both column and row pivoting, so most robust
    int rank(MatrixView<F> tmpPivoted) const;

    // same but slower : need to allocate memory instead of using preallocated tmp storage
    int rank() const;


    void get_point_coordinates(long long nb_points, double* coords, int stride = 1) const {

        typedef typename F::T T;
        enum { Q = GFCardinality<F::p, F::r>::value };

        int n_digits = cols();
        if (std::floor(log((double)nb_points) / log((double)Q)) + 1 > n_digits) {
            std::cout << "request too many points : matrix is too small" << std::endl;
        }

        Matrix<F> digits(n_digits, 1);
        Matrix<F> result_mul(n_digits, 1);

        for (long long i = 0; i < nb_points; i++) {
            decompose_integer_into_base<F>(i, digits.data(), n_digits);
            result_mul = (*this)* digits;
            double val = 0;
            for (int j = 0; j < n_digits; j++) {
                val += gf_to_raw_index<F>(result_mul[j]) / (double)std::pow((double)Q, (double)j+1);
            }
            coords[i*stride] = val;
        }
    }

    std::vector<double> get_point_coordinates(long long nb_points) const {

        std::vector<double> result(nb_points);
        get_point_coordinates(nb_points, &result[0]);
        return result;
    }

    bool save_svg(const char* filename,
        int cell_size = 24,
        int margin = 6,
        bool draw_grid = true) const;
};

template<class F>
void draw_2D_points(MatrixView<F> M1, MatrixView<F> M2, int nb_points, const std::string &svg_filename) {

    std::vector<double> x = M1.get_point_coordinates(nb_points);
    std::vector<double> y = M2.get_point_coordinates(nb_points);


    FILE* f = fopen(svg_filename.c_str(), "w+");
	fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" viewBox=\"0 0 1 1\" width = \"1024\" height = \"1024\">\n");

	fprintf(f, "<g>\n");
	for (int i = 0; i < nb_points; i++) {
		fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r = \"%3.3f\" />\n", x[i], 1.-y[i], sqrt(0.3/(nb_points*3.14)));
	}
	fprintf(f, "</g>\n");

	fprintf(f, "</svg>\n");
    fclose(f);

}



template<class F>
struct MatrixView : MatrixBase<MatrixView<F>, F> {
    typedef typename F::T T;

    T* values;
    int m, n;

    MatrixView() : values(0), m(0), n(0) {}

    MatrixView(T* values_, int m_, int n_)
        : values(values_), m(m_), n(n_) {
    }
    T operator[](int i) const {
        return values[i];
    }
    T& operator[](int i) {
        return values[i];
    }

};

template<class F>
struct Matrix : MatrixBase<Matrix<F>, F> {
    typedef typename F::T T;

    T* values;
    int m, n;

    Matrix() : values(0), m(0), n(0) {}

    Matrix(int m_, int n_)
        : values(new T[m_ * n_]), m(m_), n(n_) {
    }

    Matrix(const MatrixView<F>& view)
        : values(new T[view.m * view.n]), m(view.m), n(view.n) {
        for (int i = 0; i < m * n; ++i) {
            values[i] = view.values[i];
        }
    }

    Matrix(int m_, int n_, std::initializer_list<int> init)
        : values(new T[m_ * n_]), m(m_), n(n_) {
        assert(static_cast<int>(init.size()) <= m * n);

        int k = 0;
        for (int x : init) {
            values[k++] = gf_from_int<F::p, F::r>(x);
        }

        for (; k < m * n; ++k) {
            values[k] = T{ 0 };
        }
    }

    ~Matrix() {
        delete[] values;
    }

    Matrix(const Matrix&) = delete;
    Matrix& operator=(const Matrix&) = delete;

    Matrix(Matrix&& other) noexcept
        : values(other.values), m(other.m), n(other.n) {
        other.values = 0;
        other.m = 0;
        other.n = 0;
    }

    Matrix& operator=(Matrix&& other) noexcept {
        if (this != &other) {
            delete[] values;

            values = other.values;
            m = other.m;
            n = other.n;

            other.values = 0;
            other.m = 0;
            other.n = 0;
        }
        return *this;
    }

    MatrixView<F> view() {
        return MatrixView<F>{ values, m, n };
    }

    const MatrixView<F> view() const {
        return MatrixView<F>{ values, m, n };
    }

    T operator[](int i) const {
        return values[i];
    }
    T& operator[](int i) {
        return values[i];
    }
};


template<typename T>
void printMat(const Matrix<T>& mat, int n) {
    for (int i = 0; i < mat.m; i++) {
        for (int j = 0; j < mat.n; j++) {
            std::cout << (int)mat.values[i * mat.n + j] << "\t";
        }
        std::cout << std::endl;
    }
}


template<class F>
inline typename std::enable_if<F::r == 1, int>::type
gf_svg_index(typename F::T x) {
    enum { Q = GFCardinality<F::p, F::r>::value };
    int v = x % Q;
    if (v < 0) v += Q;
    return v;
}
template<class F>
inline typename std::enable_if<F::r != 1, int>::type
gf_svg_index(typename F::T x) {
    return x.v;
}

inline const char* gf_svg_color(int idx) {
    static const char* colors[16] = {
        "#FFFFFF", // 0
        "#A8DADC", // 1  pastel cyan
        "#FFCAD4", // 2  pastel pink
        "#CDEAC0", // 3  pastel green
        "#FFF1A8", // 4  pastel yellow
        "#CDB4DB", // 5  pastel lavender
        "#BDE0FE", // 6  pastel sky
        "#FFD6A5", // 7  pastel peach
        "#CAFFBF", // 8  pastel mint
        "#F1C0E8", // 9  pastel mauve
        "#A0C4FF", // 10 pastel blue
        "#FDFFB6", // 11 pastel lemon
        "#FFC6FF", // 12 pastel magenta
        "#D0F4DE", // 13 pastel aqua
        "#E4C1F9", // 14 pastel violet
        "#B8F2E6"  // 15 pastel turquoise
    };

    if (idx < 0) idx = 0;
    if (idx > 15) idx = 15;
    return colors[idx];
}

template<class Derived, class F>
bool MatrixBase<Derived, F>::save_svg(const char* filename,
    int cell_size,
    int margin,
    bool draw_grid) const {
    typedef typename F::T T;
    enum { Q = GFCardinality<F::p, F::r>::value };

    const int m = rows();
    const int n = cols();

    if (m < 0 || n < 0 || cell_size <= 0 || margin < 0) {
        return false;
    }

    FILE* f = std::fopen(filename, "w");
    if (!f) {
        return false;
    }

    const int width = 2 * margin + n * cell_size + 1;
    const int height = 2 * margin + m * cell_size + 1;

    std::fprintf(f,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        "version=\"1.1\" width=\"%d\" height=\"%d\" "
        "viewBox=\"0 0 %d %d\">\n",
        width, height, width, height);

    std::fprintf(f,
        "  <rect x=\"0\" y=\"0\" width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n",
        width, height);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            const int x = margin + j * cell_size;
            const int y = margin + i * cell_size;

            int idx = gf_svg_index<F>((*this)[i * n + j]);
            if (Q > 0) {
                idx %= Q;
                if (idx < 0) idx += Q;
            }

            const char* fill = gf_svg_color(idx);

            std::fprintf(f,
                "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" "
                "fill=\"%s\" stroke=\"none\"/>\n",
                x, y, cell_size, cell_size, fill);
        }
    }

    if (draw_grid) {
        const char* grid_color = "#B8B8B8";

        for (int i = 0; i <= m; ++i) {
            const int y = margin + i * cell_size;
            std::fprintf(f,
                "  <line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" "
                "stroke=\"%s\" stroke-width=\"1\"/>\n",
                margin, y, margin + n * cell_size, y, grid_color);
        }

        for (int j = 0; j <= n; ++j) {
            const int x = margin + j * cell_size;
            std::fprintf(f,
                "  <line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" "
                "stroke=\"%s\" stroke-width=\"1\"/>\n",
                x, margin, x, margin + m * cell_size, grid_color);
        }
    }

    std::fprintf(f,
        "  <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" "
        "fill=\"none\" stroke=\"#888888\" stroke-width=\"1.2\"/>\n",
        margin, margin, n * cell_size, m * cell_size);

    std::fprintf(f, "</svg>\n");
    std::fclose(f);

    return true;
}


template<int p, int r>
void pascal_pow_int(GF<p, r>* mat, int n, int power) {
    typedef GF<p, r> T;

    const T field_power = gf_from_int<p, r>(power);

    mat[0] = gf_from_int<p, r>(1);

    for (int i = 1; i < n; i++) {
        mat[i * n + 0] = gf_from_int<p, r>(0);

        mat[i] =
            gf_reduce<p, r>(
                mat[i - 1] * field_power
            );
    }

    for (int i = 1; i < n; i++) {
        for (int j = 1; j < n; j++) {
            mat[i * n + j] =
                gf_reduce<p, r>(
                    mat[i * n + j - 1] * field_power
                    + mat[(i - 1) * n + j - 1]
                );
        }
    }
}

template<int p, int r> 
void pascal_pow_field(GF<p, r>* mat, int n, GF<p, r> power) {
    typedef GF<p, r> T; mat[0] = T{ 1 }; 
    for (int i = 1; i < n; i++) { 
        mat[i * n + 0] = T{ 0 }; 
        mat[i] = gf_reduce<p, r>(mat[i - 1] * power); 
    } 
    for (int i = 1; i < n; i++) { 
        for (int j = 1; j < n; j++) { 
            mat[i * n + j] = gf_reduce<p, r>(mat[i * n + j - 1] * power + mat[(i - 1) * n + j - 1]); 
        } 
    } 
}

template<class F, int N, int Power=1>
struct PascalPowIntMatrix {
    typedef typename F::T T;

    static MatrixView<F> value() {
        static T buffer[N * N];
        static bool initialized = false;

        if (!initialized) {
            pascal_pow_int<F::p, F::r>(buffer, N, Power);
            initialized = true;
        }

        return MatrixView<F>{ buffer, N, N };
    }
};


template<class F, int N, int Power=1>
struct PascalTranslateFieldMatrix {
    typedef typename F::T T;

    static MatrixView<F> value() {
        static T buffer[N * N];
        static bool initialized = false;

        if (!initialized) {
            pascal_pow_field<F::p, F::r>(buffer, N, T{ Power });
            initialized = true;
        }

        return MatrixView<F>{ buffer, N, N };
    }
};



template<class F>
void fill_pascal_int(MatrixView<F> &out, int power_int) {
    assert(out.m == out.n);
    pascal_pow_int<F::p, F::r>(out.values, out.n, power_int);
}

template<class F>
void fill_pascal_field(MatrixView<F> &out, typename F::T power_field) {
    assert(out.m == out.n);
    pascal_pow_field<F::p, F::r>(out.values, out.n, power_field);
}



template<class F>
void tensor_product_full(const MatrixView<F>& A,
    const MatrixView<F>& B,
    MatrixView<F> result) {
    typedef typename F::T T;

    assert(result.m == A.m * B.m);
    assert(result.n == A.n * B.n);

    const T* Adata = A.values;
    const T* Bdata = B.values;
    T* Cdata = result.values;

    const int Am = A.m;
    const int An = A.n;
    const int Bm = B.m;
    const int Bn = B.n;
    const int Cn = result.n;

    for (int i = 0; i < Am; ++i) {
        const T* Arow = Adata + i * An;

        for (int i2 = 0; i2 < Bm; ++i2) {
            const T* Brow = Bdata + i2 * Bn;
            T* Crow = Cdata + (i * Bm + i2) * Cn;

            for (int j = 0; j < An; ++j) {
                const T a = Arow[j];
                T* Cblock = Crow + j * Bn;

                for (int j2 = 0; j2 < Bn; ++j2) {
                    Cblock[j2] = gf_reduce<F::p, F::r>(a * Brow[j2]);
                }
            }
        }
    }
}

template<class F>
void tensor_product_truncated(const MatrixView<F>& A,
    const MatrixView<F>& B,
    MatrixView<F> result,
    int truncate_row,
    int truncate_col) {
    typedef typename F::T T;

    const int full_m = A.m * B.m;
    const int full_n = A.n * B.n;

    assert(truncate_row >= 0 && truncate_row <= full_m);
    assert(truncate_col >= 0 && truncate_col <= full_n);

    assert(result.m >= truncate_row);
    assert(result.n >= truncate_col);

    const T* Adata = A.values;
    const T* Bdata = B.values;
    T* Cdata = result.values;

    const int Am = A.m;
    const int An = A.n;
    const int Bm = B.m;
    const int Bn = B.n;
    const int Cn = result.n;

    const int i_max = (truncate_row + Bm - 1) / Bm;
    const int j_max = (truncate_col + Bn - 1) / Bn;

    for (int i = 0; i < i_max; ++i) {
        const T* Arow = Adata + i * An;

        const int remaining_rows = truncate_row - i * Bm;
        const int i2_limit = std::min(Bm, remaining_rows);

        for (int i2 = 0; i2 < i2_limit; ++i2) {
            const T* Brow = Bdata + i2 * Bn;
            T* Crow = Cdata + (i * Bm + i2) * Cn;

            for (int j = 0; j < j_max; ++j) {
                const T a = Arow[j];
                T* Cblock = Crow + j * Bn;

                const int remaining_cols = truncate_col - j * Bn;
                const int j2_limit = std::min(Bn, remaining_cols);

                for (int j2 = 0; j2 < j2_limit; ++j2) {
                    Cblock[j2] = gf_reduce<F::p, F::r>(a * Brow[j2]);
                }
            }
        }
    }
}

template<class F>
void tensor_product(const MatrixView<F> A,
    const MatrixView<F> B,
    MatrixView<F> result,
    int truncate_row = -1,
    int truncate_col = -1) {
    if (truncate_row < 0 && truncate_col < 0) {
        tensor_product_full(A, B, result);
        return;
    }

    const int full_m = A.m * B.m;
    const int full_n = A.n * B.n;

    if (truncate_row < 0) truncate_row = full_m;
    if (truncate_col < 0) truncate_col = full_n;

    tensor_product_truncated(A, B, result, truncate_row, truncate_col);
}

template<class AExpr, class BExpr, class F>
Matrix<F> matmul(const MatrixBase<AExpr, F>& A_,
    const MatrixBase<BExpr, F>& B_) {
    const AExpr& A = A_.derived();
    const BExpr& B = B_.derived();

    Matrix<F> result(A.m, B.n);

    for (int i = 0; i < A.m; ++i) {
        for (int j = 0; j < B.n; ++j) {
            typename F::T acc = F::zero();

            for (int k = 0; k < A.n; ++k) {
                acc = acc + A.values[i * A.n + k]*B.values[k * B.n + j];
            }

            result.values[i * result.n + j] = gf_reduce<F::p, F::r>(acc);
        }
    }

    return result;
}

template<class AExpr, class BExpr, class F>
Matrix<F> operator*(const MatrixBase<AExpr, F>& A,
    const MatrixBase<BExpr, F>& B) {
    return matmul(A, B);
}

template<class AExpr, class BExpr, class F>
Matrix<F> operator+(const MatrixBase<AExpr, F>& A_,
    const MatrixBase<BExpr, F>& B_) {
    const AExpr& A = A_.derived();
    const BExpr& B = B_.derived();

    Matrix<F> result(A.m, A.n);

    for (int i = 0; i < A.m * A.n; ++i) {
        result.values[i] = A.values[i]+B.values[i];
    }

    return result;
}
template<class AExpr, class BExpr, class F>
Matrix<F> operator-(const MatrixBase<AExpr, F>& A_,
    const MatrixBase<BExpr, F>& B_) {
    const AExpr& A = A_.derived();
    const BExpr& B = B_.derived();

    Matrix<F> result(A.m, A.n);

    for (int i = 0; i < A.m * A.n; ++i) {
        result.values[i] = A.values[i]-B.values[i];
    }

    return result;
}

template<class Expr, class F>
std::ostream& operator<<(std::ostream& out,
    const MatrixBase<Expr, F>& M_) {
    const Expr& M = M_.derived();

    out << "{";
    for (int i = 0; i < M.m; ++i) {
        out << "{";
        for (int j = 0; j < M.n-1; ++j) {
            out << M.values[i * M.n + j] << ",\t";
        }
        out << M.values[i * M.n + M.n-1];
        if (i< M.m-1)
            out << "},\n";
        else
            out << "}\n";
    }
    out << "}";

    return out;
}



template<typename F>
void combine_matrices(const MatrixView<F>  *matrices, const int* row_choices, int num_matrices, int max_columns, MatrixView<F> result) {
    assert(result.n == max_columns);

    typedef typename F::T T;

    const int row_bytes = max_columns * sizeof(T);
    int row_num = 0;
    for (int i = 0; i < num_matrices; i++) {
        const MatrixView<F>& curMat = matrices[i];
        const T* __restrict src = curMat.data();
        T* __restrict dst = result.data() + row_num * max_columns;
        int r_i = row_choices[i];
        for (int row = 0; row < r_i; row++) {
            std::memcpy(dst, src, row_bytes);
            dst += max_columns;
            src += curMat.n;
        }
        row_num += r_i;
    }
}

template<class Derived, class F>
int MatrixBase<Derived,F>::checkDet(MatrixView<F> tmpMatT, MatrixView<F> tmpInv, MatrixView<F> tmpiAc, MatrixView<F> tmpriA) const {

    assert(rows() == cols());

	const MatrixBase<Derived, F>& mat = *this;
	int m = rows();
	if (m == 0) return -1;

	typedef typename F::T T;


	// a_ij aliases 
	T a00 = gf_reduce<F::p, F::r>(mat[0 * m + 0]);

	T det = a00;
	if (det == T{ 0 }) return 0;
	if (m == 1) return -1;

	T a01 = mat[0 * m + 1];
	T a10 = mat[1 * m + 0], a11 = mat[1 * m + 1];

	T det2x2 = gf_reduce<F::p, F::r>(a00 * a11 - a10 * a01);
	if (det2x2 == T{ 0 }) return 1;
	if (m == 2) return -1;

	T a02 = mat[0 * m + 2];
	T a12 = mat[1 * m + 2];
	T a20 = mat[2 * m + 0], a21 = mat[2 * m + 1], a22 = mat[2 * m + 2];

	T det2d_01 = (a10 * a21 - a20 * a11);                // rows(1,2), cols(0,1)
	T det2f_01 = (a00 * a21 - a20 * a01);                // rows(0,2), cols(0,1)
	T det3x3 = gf_reduce<F::p, F::r>(a02 * det2d_01 - a12 * det2f_01 + a22 * det2x2);
	if (det3x3 == T{ 0 }) return 2;
	if (m == 3) return -1;

	// extra 2x2 minors needed (still only rows 0..2)
	T a03 = mat[0 * m + 3];
	T a13 = mat[1 * m + 3];
	T a23 = mat[2 * m + 3];
	T a30 = mat[3 * m + 0], a31 = mat[3 * m + 1], a32 = mat[3 * m + 2], a33 = mat[3 * m + 3];

	T det2x2_02 = (a00 * a12 - a10 * a02);               // cols(0,2), rows(0,1)
	T det2d_02 = (a10 * a22 - a20 * a12);               // cols(0,2), rows(1,2)
	T det2f_02 = (a00 * a22 - a20 * a02);               // cols(0,2), rows(0,2)

	T det2x2_12 = (a01 * a12 - a11 * a02);               // cols(1,2), rows(0,1)
	T det2d_12 = (a11 * a22 - a21 * a12);               // cols(1,2), rows(1,2)
	T det2f_12 = (a01 * a22 - a21 * a02);               // cols(1,2), rows(0,2)

	// the three other 3x3 leading minors using column 3
	T D2 = (a03 * det2d_01 - a13 * det2f_01 + a23 * det2x2); // cols(0,1,3)
	T D1 = (a03 * det2d_02 - a13 * det2f_02 + a23 * det2x2_02); // cols(0,2,3)
	T D0 = (a03 * det2d_12 - a13 * det2f_12 + a23 * det2x2_12); // cols(1,2,3)
	T D3 = det3x3;                                                // cols(0,1,2)

	// 4x4 determinant by expansion along row 3
	T det4x4 = gf_reduce<F::p, F::r>(a30 * D0 - a31 * D1 + a32 * D2 - a33 * D3);
	if (det4x4 == T{ 0 }) return 3;
	if (m == 4) return -1;
	// else keep going...*/



// -----------------------------------------------------------------------------
// m >= 5, faster 5x5 test using Sylvester identity.
// Uses the already nonzero det3x3.
// -----------------------------------------------------------------------------

	T a04 = mat[0 * m + 4];
	T a14 = mat[1 * m + 4];
	T a24 = mat[2 * m + 4];
	T a34 = mat[3 * m + 4];

	T a40 = mat[4 * m + 0];
	T a41 = mat[4 * m + 1];
	T a42 = mat[4 * m + 2];
	T a43 = mat[4 * m + 3];
	T a44 = mat[4 * m + 4];


	// 3x3 minors involving column 4 and columns 0,1,2.
	// Rows are 0,1,2.
	T C014 = (a04 * det2d_01 - a14 * det2f_01 + a24 * det2x2);     // cols(0,1,4)
	T C024 = (a04 * det2d_02 - a14 * det2f_02 + a24 * det2x2_02);  // cols(0,2,4)
	T C124 = (a04 * det2d_12 - a14 * det2f_12 + a24 * det2x2_12);  // cols(1,2,4)

	// 4x4 bordered determinants around the leading 3x3 block.
	//
	// Same sign convention as your det4x4:
	// value = negative of the standard 4x4 determinant,
	// but the sign cancels in B33 * B44 - B34 * B43.

	T B33 = det4x4;                                                // rows 0,1,2,3 ; cols 0,1,2,3
	T B34 = (a30 * C124 - a31 * C024 + a32 * C014 - a34 * D3);     // rows 0,1,2,3 ; cols 0,1,2,4
	T B43 = (a40 * D0 - a41 * D1 + a42 * D2 - a43 * D3);     // rows 0,1,2,4 ; cols 0,1,2,3
	T B44 = (a40 * C124 - a41 * C024 + a42 * C014 - a44 * D3);     // rows 0,1,2,4 ; cols 0,1,2,4

	// Sylvester:
	//
	// det5x5 * det3x3 = B33 * B44 - B34 * B43
	//
	// Since det3x3 != 0 mod GF, det5x5 == 0 mod GF iff:
	// B33 * B44 - B34 * B43 == 0 mod GF.

	T det5_test = gf_reduce<F::p, F::r>(B33 * B44 - B34 * B43);

	if (det5_test == T{ 0 }) return 4;
	if (m == 5) return -1;


	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			tmpMatT[i * m + j] = mat[j * m + i];
		}
	}

	tmpInv[0] = F::div(T{ 1 }, det);

	for (int i = 1; i < m; i++) {

		const T* mat_row_i = mat.data() + i * m;
		const T* mat_col_i = tmpMatT.data() + i * m;
		// update the inverse of mat with sherman morrison
		T rac{ 0 };
		for (int u = 0; u < i; u++) {
			T ac{ 0 };
			T ra{ 0 };
			const T* invM_row_u = tmpInv.data() + u * m;
			for (int v = 0; v < i; v++) {
				ac += invM_row_u[v] * mat_col_i[v];
				ra += mat_row_i[v] * tmpInv[v * m + u];
			}
			tmpiAc[u] = gf_reduce<F::p, F::r>(ac);
			tmpriA[u] = gf_reduce<F::p, F::r>(ra);
			rac += mat_row_i[u] * tmpiAc[u];
		}
		T alpha = gf_reduce<F::p, F::r>(mat_row_i[i] - rac);

		if (alpha == T{ 0 }) return i;

		det *= alpha; // update the determinant with the determinant lemma

		T invalpha = F::div(T{ 1 }, alpha);

		for (int u = 0; u < i; u++) {
			T iAcuInvalpha = gf_reduce<F::p, F::r>(tmpiAc[u] * invalpha);
			for (int v = 0; v < i; v++) {
				tmpInv[u * m + v] += iAcuInvalpha * tmpriA[v];
			}
			tmpInv[u * m + i] = -iAcuInvalpha;
			tmpInv[i * m + u] = gf_reduce<F::p, F::r>(-(tmpriA[u] * invalpha));

		}
		tmpInv[i * m + i] = invalpha;

	}
	return -1;

}




template<class Derived, class F>
Matrix<F> MatrixBase<Derived, F>::get_inverse(MatrixView<F> tmpMatT, MatrixView<F> tmpiAc, MatrixView<F> tmpriA) const {

    assert(rows() == cols());

    const MatrixBase<Derived, F>& mat = *this;
    int m = rows();
   
    if (m == 0) {
        return Matrix<F>(0, 0);
    }
    
    Matrix<F> Inv(rows(), cols());

    typedef typename F::T T;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            tmpMatT[i * m + j] = mat[j * m + i];
        }
    }

    T a00 = gf_reduce<F::p, F::r>(mat[0 * m + 0]);
    Inv[0] = F::div(T{ 1 }, a00);


    for (int i = 1; i < m; i++) {

        const T* mat_row_i = mat.data() + i * m;
        const T* mat_col_i = tmpMatT.data() + i * m;
        // update the inverse of mat with sherman morrison
        T rac{ 0 };
        for (int u = 0; u < i; u++) {
            T ac{ 0 };
            T ra{ 0 };
            const T* invM_row_u = Inv.data() + u * m;
            for (int v = 0; v < i; v++) {
                ac += invM_row_u[v] * mat_col_i[v];
                ra += mat_row_i[v] * Inv[v * m + u];
            }
            tmpiAc[u] = gf_reduce<F::p, F::r>(ac);
            tmpriA[u] = gf_reduce<F::p, F::r>(ra);
            rac += mat_row_i[u] * tmpiAc[u];
        }
        T alpha = gf_reduce<F::p, F::r>(mat_row_i[i] - rac);

        if (alpha == T{ 0 }) {
            std::cout << "get_inverse on non invertible matrix, or would require pivoting" << std::endl;
            return Matrix<F>(0,0);
        }


        T invalpha = F::div(T{ 1 }, alpha);

        for (int u = 0; u < i; u++) {
            T iAcuInvalpha = gf_reduce<F::p, F::r>(tmpiAc[u] * invalpha);
            for (int v = 0; v < i; v++) {
                Inv[u * m + v] += iAcuInvalpha * tmpriA[v];
            }
            Inv[u * m + i] = -iAcuInvalpha;
            Inv[i * m + u] = gf_reduce<F::p, F::r>(-(tmpriA[u] * invalpha));

        }
        Inv[i * m + i] = invalpha;

    }
    return Inv;

}

template<class Derived, class F>
Matrix<F> MatrixBase<Derived, F>::get_inverse() const {

    int m = rows();
    typedef typename F::T T;
    T* values = new T[m * m  + m * 2];
    MatrixView<F> tmpMatT(values, m, m);
    MatrixView<F> tmpiAc(values + m * m, m, 1);
    MatrixView<F> tmpriA(values + m * m + m, m, 1);

    Matrix<F> inv = get_inverse(tmpMatT, tmpiAc, tmpriA);
    inv.reduce();
    delete[] values;
    return inv;
}


template<class Derived, class F>
int MatrixBase<Derived, F>::checkDet() const {

   int m = rows();
   typedef typename F::T T;
   T* values = new T[m * m * 2 + m * 2];
   MatrixView<F> tmpMatT(values, m, m);
   MatrixView<F> tmpInv(values+m*m, m, m);
   MatrixView<F> tmpiAc(values + 2*m * m, m, 1);
   MatrixView<F> tmpriA(values + 2*m * m+m, m, 1);

   int det = checkDet(tmpMatT, tmpInv, tmpiAc, tmpriA);
   delete[] values;
   return det;
}



// returns -1 if it is progressive
// otherwise returns the index of size m that failed determinant tests
template<class F>
int is_t0_progressive(const MatrixView<F>* all_matrices, int n_matrices,
    MatrixView<F> combination, // n x n 
    MatrixView<F> matT, // n x n 
    MatrixView<F> tmpInv, // n x n  (with n = max_size_to_test)
    MatrixView<F> tmpiAc, // n x 1 
    MatrixView<F> tmpriA // 1 x n 

    ) {  

    int m = all_matrices[0].cols();
    int s = n_matrices;

    int minMFail = m;
    bool succeed = true;

    std::vector<int>  rows_combinations(s, 0);
    int position_value[200];
    int position_value_size = 0;
    for (int i = 0; i <= m; i++) {
        position_value[position_value_size * 2] = 0;
        position_value[position_value_size * 2 + 1] = i;
        position_value_size++;
    }
    long long id = 0;
    while (position_value_size != 0) {
        id++;
        int pos = position_value[(position_value_size - 1) * 2];
        int value = position_value[(position_value_size - 1) * 2 + 1];
        position_value_size--;
        if (pos < s) {
            rows_combinations[pos] = value; 
        }
        if (pos == s - 1) {

            int sumk = 0;
            for (int i = 0; i <= pos; i++) {
                sumk += rows_combinations[i];
            }

            combination.m = sumk;
            combination.n = sumk;
            combine_matrices(all_matrices, &rows_combinations[0], s, sumk, combination);

            int mfail = combination.checkDet(matT, tmpInv, tmpiAc, tmpriA);

            if (mfail != -1) {
                return mfail;
            }

            /*if ( mfail != -1) {
                minMFail = std::min(minMFail, mfail);
                succeed = false;
            }*/

        }
        else {
            pos++;

            int sumk = 0;
            for (int i = 0; i < pos; i++) {
                sumk += rows_combinations[i];
            }
            for (int k = 0; k <= m - sumk; k++) {
                position_value[position_value_size * 2] = pos;
                position_value[position_value_size * 2 + 1] = k;
                position_value_size++;
            }
        }
    }

    return -1;
}

template<class F>
bool is_t0_progressive(const MatrixView<F>* all_matrices, int n_matrices) {

    typedef typename F::T T;
    int m = all_matrices[0].cols();
    T* values = new T[m * m * 3 + m * 2];
    //memset(values, 0, (m * m * 3 + m * 2) * sizeof(T));
    MatrixView<F> combination(values, m, m);
    MatrixView<F> matT(values + m * m, m, m);
    MatrixView<F> tmpInv(values + m*m*2, m, m); 
    MatrixView<F> tmpiAc(values+m*m*3, m, 1);
    MatrixView<F> tmpriA(values+m*m*3+m, 1, m);  
    
    

    int ret = is_t0_progressive(all_matrices, n_matrices, combination, matT, tmpInv, tmpiAc, tmpriA);
    delete[] values;

    return (ret==-1);

}


template<class Derived, class F>
int MatrixBase<Derived, F>::rank(MatrixView<F> tmpPivoted) const {
    typedef typename F::T T;

    const MatrixBase<Derived, F>& mat = *this;

    const int m = rows();
    const int n = cols();

    if (m == 0 || n == 0) {
        return 0;
    }

    assert(tmpPivoted.m >= m);
    assert(tmpPivoted.n >= n);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            tmpPivoted[i * tmpPivoted.n + j] =
                gf_reduce<F::p, F::r>(mat[i * n + j]);
        }
    }

    int rank = 0;

    for (int col = 0; col < n && rank < m; ++col) {
        int pivot = -1;

        // find non nul pivot in current column
        for (int row = rank; row < m; ++row) {
            if (tmpPivoted[row * tmpPivoted.n + col] != T{ 0 }) {
                pivot = row;
                break;
            }
        }

        if (pivot == -1) {
            continue;
        }

        // pivot on rank row
        if (pivot != rank) {
            for (int j = col; j < n; ++j) {
                std::swap(
                    tmpPivoted[rank * tmpPivoted.n + j],
                    tmpPivoted[pivot * tmpPivoted.n + j]
                );
            }
        }

        const T pivotValue = tmpPivoted[rank * tmpPivoted.n + col];

        // pivot elimination
        for (int row = rank + 1; row < m; ++row) {
            const T x = tmpPivoted[row * tmpPivoted.n + col];

            if (x == T{ 0 }) {
                continue;
            }

            const T coeff = F::div(x, pivotValue);

            // column col becomes 0
            tmpPivoted[row * tmpPivoted.n + col] = T{ 0 };

            for (int j = col + 1; j < n; ++j) {
                tmpPivoted[row * tmpPivoted.n + j] =
                    gf_reduce<F::p, F::r>(
                        tmpPivoted[row * tmpPivoted.n + j]
                        - coeff * tmpPivoted[rank * tmpPivoted.n + j]
                    );
            }
        }

        ++rank;
    }

    return rank;
}


template<class Derived, class F>
int MatrixBase<Derived, F>::rank() const {
    typedef typename F::T T;
    int m = rows();
    int n = cols();
    T* values = new T[m * n];
    MatrixView<F> tmpPivoted(values, m, n);
    int r = rank(tmpPivoted);
    delete[] values;

    return r;
}

template<class F>
std::vector<int> t_values_naive(const MatrixView<F>* all_matrices, int n_matrices)  {

    typedef typename F::T T;
    T* mem = new T[all_matrices[0].rows() * all_matrices[0].cols()*2];
    MatrixView<F> combination(mem, all_matrices[0].rows(), all_matrices[0].cols());
    MatrixView<F> tmpPivoted(mem+ all_matrices[0].rows() * all_matrices[0].cols(), all_matrices[0].rows(), all_matrices[0].cols());
    std::vector<int> t = t_values_naive(all_matrices, n_matrices, combination, tmpPivoted);
    delete[] mem;
    return t;
}

// returns a vector of t values (one value for each m)
template<class F>
std::vector<int> t_values_naive(const MatrixView<F>* all_matrices, int n_matrices,
    MatrixView<F> combination, // n x n 
    MatrixView<F> tmpPivoted // n x n
)  {
    // todo, implement https://arxiv.org/pdf/1910.02277

    typedef typename F::T T;

    int s = n_matrices;
    std::vector<int>  rows_combinations(s, 0);
    int position_value[200];
    int position_value_size = 0;
   
    std::vector<int> result;
    for (int m = 1; m <= all_matrices[0].cols(); m++) {

        int t_trial;
        for (t_trial = 0; t_trial < m; t_trial++) {

            position_value_size = 0;

            for (int i = 0; i <= m-t_trial; i++) {
                position_value[position_value_size * 2] = 0;
                position_value[position_value_size * 2 + 1] = i;
                position_value_size++;
            }
            long long id = 0;
            while (position_value_size != 0) {
                id++;
                int pos = position_value[(position_value_size - 1) * 2];
                int value = position_value[(position_value_size - 1) * 2 + 1];
                position_value_size--;
                if (pos < s) {
                    rows_combinations[pos] = value;
                }
                if (pos == s - 1) {

                    int sumk = 0;
                    for (int i = 0; i <= pos; i++) {
                        sumk += rows_combinations[i];
                    }

                    if (sumk == m - t_trial) {
                        combination.m = sumk ;
                        combination.n = sumk + t_trial;
                        combine_matrices(all_matrices, &rows_combinations[0], s, m, combination);


                        int rk = combination.rank(tmpPivoted);
                        if (rk < sumk) {
                            goto next_t_trial;
                        }
                    }


                }
                else {
                    pos++;

                    int sumk = 0;
                    for (int i = 0; i < pos; i++) {
                        sumk += rows_combinations[i];
                    }
                    for (int k = 0; k <= m - sumk-t_trial; k++) {
                        position_value[position_value_size * 2] = pos;
                        position_value[position_value_size * 2 + 1] = k;
                        position_value_size++;
                    }
                }
            }
            break;
        next_t_trial:
            continue;
        }
        result.push_back(t_trial);
    }

    return result;
}


// from here :
// https://www.sciencedirect.com/science/article/pii/S0377042719306740
// An algorithm to compute the t-value of a digital net and of its projections
// Pure ChatGPT implementation

template<class Callback>
void gray_compositions_rec(
    int parts,
    int total,
    bool rev,
    int pos,
    std::vector<int>& comp,
    Callback& cb
) {
    if (parts == 1) {
        comp[pos] = total;
        cb(comp);
        return;
    }

    if (!rev) {
        for (int v = 0; v <= total; ++v) {
            comp[pos] = v;
            bool subrev = (v & 1) ? !rev : rev;
            gray_compositions_rec(parts - 1, total - v, subrev, pos + 1, comp, cb);
        }
    }
    else {
        for (int v = total; v >= 0; --v) {
            comp[pos] = v;
            bool subrev = (v & 1) ? !rev : rev;
            gray_compositions_rec(parts - 1, total - v, subrev, pos + 1, comp, cb);
        }
    }
}

template<class Callback>
void gray_compositions(int parts, int total, Callback cb) {
    std::vector<int> comp(parts, 0);
    gray_compositions_rec(parts, total, false, 0, comp, cb);
}


template<class F>
struct RAREF_LT_State {
    typedef typename F::T T;

    int q;
    int k;

    Matrix<F> C;   // q x k
    Matrix<F> L;   // q x q
    Matrix<F> Tm;  // q x k, Tm = L*C

    std::vector<int> pivot_col;
    int rank_value;

    RAREF_LT_State(int q_, int k_)
        : q(q_),
        k(k_),
        C(q_, k_),
        L(q_, q_),
        Tm(q_, k_),
        pivot_col(q_, -1),
        rank_value(0) {
    }

    static bool is_zero(T x) {
        return x == T{ 0 };
    }

    void set_C_row_no_update(int row, const MatrixView<F>& M, int source_row) {
        for (int c = 0; c < k; ++c) {
            C[row * C.n + c] =
                gf_reduce<F::p, F::r>(M[source_row * M.n + c]);
        }
    }

    void set_identity_L() {
        for (int i = 0; i < q; ++i) {
            for (int j = 0; j < q; ++j) {
                L[i * L.n + j] = (i == j) ? T{ 1 } : T{ 0 };
            }
        }
    }

    void copy_C_to_T() {
        for (int i = 0; i < q; ++i) {
            for (int j = 0; j < k; ++j) {
                Tm[i * Tm.n + j] = C[i * C.n + j];
            }
        }
    }

    void swap_rows(Matrix<F>& M, int a, int b, int cols) {
        if (a == b) return;

        for (int c = 0; c < cols; ++c) {
            std::swap(M[a * M.n + c], M[b * M.n + c]);
        }
    }

    void row_scale(Matrix<F>& M, int row, T scale, int cols) {
        for (int c = 0; c < cols; ++c) {
            M[row * M.n + c] =
                gf_reduce<F::p, F::r>(M[row * M.n + c] * scale);
        }
    }

    void row_sub_scaled(Matrix<F>& M, int dst, int src, T scale, int cols) {
        if (is_zero(scale)) return;

        for (int c = 0; c < cols; ++c) {
            M[dst * M.n + c] =
                gf_reduce<F::p, F::r>(
                    M[dst * M.n + c] - scale * M[src * M.n + c]
                );
        }
    }

    void pivot_single_row(int row) {
        // Step 1 of Algorithm 1 for this row:
        // zero entries in already-pivot columns.
        for (int r = 0; r < q; ++r) {
            if (r == row) continue;

            const int pc = pivot_col[r];
            if (pc < 0) continue;

            const T coeff = Tm[row * Tm.n + pc];

            if (!is_zero(coeff)) {
                row_sub_scaled(Tm, row, r, coeff, k);
                row_sub_scaled(L, row, r, coeff, q);
            }
        }

        // Locate first nonzero coefficient.
        int pc = -1;
        for (int c = 0; c < k; ++c) {
            if (!is_zero(Tm[row * Tm.n + c])) {
                pc = c;
                break;
            }
        }

        if (pc < 0) {
            pivot_col[row] = -1;
            return;
        }

        // Normalize pivot to 1.
        const T inv_pivot = F::div(T{ 1 }, Tm[row * Tm.n + pc]);

        row_scale(Tm, row, inv_pivot, k);
        row_scale(L, row, inv_pivot, q);

        // Zero pivot column everywhere else.
        for (int r = 0; r < q; ++r) {
            if (r == row) continue;

            const T coeff = Tm[r * Tm.n + pc];

            if (!is_zero(coeff)) {
                row_sub_scaled(Tm, r, row, coeff, k);
                row_sub_scaled(L, r, row, coeff, q);
            }
        }

        pivot_col[row] = pc;
    }

    void compute_from_C() {
        set_identity_L();
        copy_C_to_T();

        std::fill(pivot_col.begin(), pivot_col.end(), -1);

        for (int row = 0; row < q; ++row) {
            pivot_single_row(row);
        }

        recompute_rank_value();
    }

    void recompute_rank_value() {
        rank_value = 0;

        for (int i = 0; i < q; ++i) {
            if (pivot_col[i] >= 0) {
                ++rank_value;
            }
        }
    }

    int rank() const {
        return rank_value;
    }

    bool full_rank() const {
        return rank_value == q;
    }

    int rightmost_pivot_plus_one() const {
        int r = 0;

        for (int i = 0; i < q; ++i) {
            if (pivot_col[i] >= 0 && pivot_col[i] + 1 > r) {
                r = pivot_col[i] + 1;
            }
        }

        return r;
    }

    void replace_row(int row, const MatrixView<F>& M, int source_row) {
        // Replace row `row` of C by row `source_row` of M,
        // and update L,Tm using Algorithm 2.

        // Algorithm 2, step 1:
        // find j such that L[j,row] != 0.
        int jrow = -1;
        for (int j = 0; j < q; ++j) {
            if (!is_zero(L[j * L.n + row])) {
                jrow = j;
                break;
            }
        }

        assert(jrow >= 0);

        if (jrow != row) {
            swap_rows(L, row, jrow, q);
            swap_rows(Tm, row, jrow, k);
            std::swap(pivot_col[row], pivot_col[jrow]);
        }

        const T a = L[row * L.n + row];
        assert(!is_zero(a));

        // Algorithm 2, step 2:
        // zero the rest of column `row` in L, applying same ops to Tm.
        for (int r = 0; r < q; ++r) {
            if (r == row) continue;

            const T coeff = L[r * L.n + row];

            if (!is_zero(coeff)) {
                const T lambda = F::div(coeff, a);

                row_sub_scaled(L, r, row, lambda, q);
                row_sub_scaled(Tm, r, row, lambda, k);
            }
        }

        // Algorithm 2, step 3:
        // isolate L[row,row] = a and replace C[row] and Tm[row].
        for (int c = 0; c < q; ++c) {
            L[row * L.n + c] = T{ 0 };
        }
        L[row * L.n + row] = a;

        for (int c = 0; c < k; ++c) {
            C[row * C.n + c] =
                gf_reduce<F::p, F::r>(M[source_row * M.n + c]);

            Tm[row * Tm.n + c] =
                gf_reduce<F::p, F::r>(a * C[row * C.n + c]);
        }

        // Algorithm 2, step 4:
        // pivot on this row.
        pivot_col[row] = -1;
        pivot_single_row(row);
        recompute_rank_value();
    }
};


template<class F>
std::vector<int> t_values(
    const MatrixView<F>* all_matrices,
    int n_matrices, int max_m = -1
) {
    assert(n_matrices > 0);

    const int s = n_matrices;
    const int k = (max_m==-1)?all_matrices[0].cols():max_m;

    for (int d = 0; d < s; ++d) {
        assert(all_matrices[d].rows() >= k);
        assert(all_matrices[d].cols() >= k);
    }

    std::vector<int> lq(k + 1, k + 1);

    for (int q = 1; q <= k; ++q) {
        RAREF_LT_State<F> state(q, k);

        std::vector<int> prev_comp(s, 0);
        std::vector<int> slot_dim(q, -1);
        std::vector<int> slot_row(q, -1);

        bool first = true;
        bool all_full_rank = true;
        int worst_needed_cols = 0;

        gray_compositions(s, q, [&](const std::vector<int>& comp) {
            if (!all_full_rank) {
                return;
            }

            if (first) {
                int out_row = 0;

                for (int d = 0; d < s; ++d) {
                    for (int r = 0; r < comp[d]; ++r) {
                        state.set_C_row_no_update(out_row, all_matrices[d], r);

                        slot_dim[out_row] = d;
                        slot_row[out_row] = r;

                        ++out_row;
                    }
                }

                assert(out_row == q);

                state.compute_from_C();

                first = false;
            }
            else {
                int from_dim = -1;
                int to_dim = -1;

                for (int d = 0; d < s; ++d) {
                    const int diff = comp[d] - prev_comp[d];

                    if (diff == -1) {
                        from_dim = d;
                    }
                    else if (diff == +1) {
                        to_dim = d;
                    }
                    else {
                        assert(diff == 0);
                    }
                }

                assert(from_dim >= 0);
                assert(to_dim >= 0);

                const int removed_row = prev_comp[from_dim] - 1;
                const int added_row = comp[to_dim] - 1;

                int slot = -1;

                for (int r = 0; r < q; ++r) {
                    if (slot_dim[r] == from_dim && slot_row[r] == removed_row) {
                        slot = r;
                        break;
                    }
                }

                assert(slot >= 0);

                state.replace_row(slot, all_matrices[to_dim], added_row);

                slot_dim[slot] = to_dim;
                slot_row[slot] = added_row;
            }

            if (!state.full_rank()) {
                all_full_rank = false;
                return;
            }

            const int needed_cols = state.rightmost_pivot_plus_one();

            if (needed_cols > worst_needed_cols) {
                worst_needed_cols = needed_cols;
            }

            prev_comp = comp;
            });

        if (!all_full_rank) {
            for (int qq = q; qq <= k; ++qq) {
                lq[qq] = k + 1;
            }
            break;
        }

        lq[q] = worst_needed_cols;
    }

    std::vector<int> result;
    result.reserve(k);

    for (int m = 1; m <= k; ++m) {
        int rho = 0;

        for (int q = 1; q <= m; ++q) {
            if (lq[q] <= m) {
                rho = q;
            }
        }

        result.push_back(m - rho);
    }

    return result;
}

/// end of ChatGPT implem of [Marion et al. 2020]


template<class F>
void get_points(
    const MatrixView<F>* all_matrices,
    int n_matrices,
    long long n_points,
    double* points
) {
 
    for (int k = 0; k < n_matrices; k++) {
        all_matrices[k].get_point_coordinates(n_points, points+k, n_matrices);
    }
    return;

}

template<class F>
std::vector<double> get_points(
    const MatrixView<F>* all_matrices,
    int n_matrices,
    long long n_points
) {
    std::vector<double> points(n_points * n_matrices);
    get_points(all_matrices, n_matrices, n_points, &points[0]);
    return points;
}

template<class F>
void fill_sobol(MatrixView<F> out,
    MatrixView<F> poly,
    MatrixView<F> V) {
    typedef typename F::T T;

    assert(out.m == out.n);
    assert(V.m == V.n);

    const int m = out.n;
    const int e = V.n;

    assert(e > 0);
    assert(poly.m * poly.n >= e);

    const int os = out.n;
    const int vs = V.n;

    if (m <= e) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < m; ++j) {
                out[i * os + j] = gf_reduce<F::p, F::r>(V[i * vs + j]);
            }
        }
        return;
    }


    for (int i = 0; i < e; ++i) {
        for (int j = 0; j < e; ++j) {
            out[i * os + j] = gf_reduce<F::p, F::r>(V[i * vs + j]);
        }
    }

    for (int ncur = e; ncur < m; ++ncur) {
        const int new_col = ncur;
        const int new_row = ncur;

        for (int j = 0; j <= new_col; ++j) {
            out[new_row * os + j] = T{ 0 };
        }

        for (int i = 0; i <= new_row; ++i) {
            out[i * os + new_col] = T{ 0 };
        }

        for (int k = 0; k < e; ++k) {
            const int src_col = ncur - e + k;
            const T coeff = poly[k];

            for (int i = 0; i <= ncur; ++i) {
                out[i * os + new_col] =
                    gf_reduce<F::p, F::r>(
                        out[i * os + new_col]
                        - coeff * out[i * os + src_col]
                    );
            }
        }

        const int feedback_col = ncur - e;

        for (int i = e; i <= ncur; ++i) {
            out[i * os + new_col] =
                gf_reduce<F::p, F::r>(
                    out[i * os + new_col]
                    + out[(i - e) * os + feedback_col]
                );
        }
    }
}



template<class F>
struct SobolMatrix : public Matrix<F> {
    typedef Matrix<F> Base;
    typedef typename F::T T;

    SobolMatrix(int n,
        MatrixView<F> poly,
        MatrixView<F> init)
        : Base(n, n) {
        fill_sobol<F>(this->view(), poly, init);
    }

};

template<class F>
struct SobolMatrixView : public MatrixView<F> {
    typedef MatrixView<F> Base;
    typedef typename F::T T;

    SobolMatrixView(MatrixView<F> out,
        MatrixView<F> poly,
        MatrixView<F> init)
        : Base(out) {
        fill_sobol<F>(
            static_cast<MatrixView<F>&>(*this),
            poly,
            init
        );
    }


    SobolMatrixView(T* buffer,
        int n,
        MatrixView<F> poly,
        MatrixView<F> init)
        : Base(buffer, n, n) {
        fill_sobol<F>(
            static_cast<MatrixView<F>&>(*this),
            poly,
            init
        );
    }


};

enum MatrixRangeOrder {
    MatrixRangeSequential,
    MatrixRangeRandomized
};

inline uint64_t gcd_u64(uint64_t a, uint64_t b) {
    while (b != 0) {
        const uint64_t r = a % b;
        a = b;
        b = r;
    }
    return a;
}

inline uint64_t splitmix64_next(uint64_t& x) {
    uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

template<class F>
struct MatrixRange {
    typedef typename F::T T;

    Matrix<F> minM;
    Matrix<F> maxM;
    Matrix<F> curM;

    std::vector<int> free_index;
    std::vector<int> free_min;
    std::vector<int> free_extent;

    uint64_t total;
    uint64_t start_rank;
    uint64_t stride;

    MatrixRangeOrder order;

    MatrixRange(MatrixView<F> min_,
        MatrixView<F> max_,
        MatrixRangeOrder order_ = MatrixRangeSequential,
        uint64_t seed = 0x123456789ABCDEF0ULL)
        : minM(min_),
        maxM(max_),
        curM(min_),
        total(1),
        start_rank(0),
        stride(1),
        order(order_) {
        assert(min_.m == max_.m);
        assert(min_.n == max_.n);

        enum { Q = GFCardinality<F::p, F::r>::value };

        const int N = minM.m * minM.n;

        for (int idx = 0; idx < N; ++idx) {
            const int a = gf_to_raw_index<F>(minM[idx]);
            const int b = gf_to_raw_index<F>(maxM[idx]);

            assert(0 <= a && a < Q);
            assert(0 <= b && b < Q);
            assert(a <= b);

            const int extent = b - a + 1;

            if (extent > 1) {
                free_index.push_back(idx);
                free_min.push_back(a);
                free_extent.push_back(extent);

                assert(total <=
                    std::numeric_limits<uint64_t>::max() / uint64_t(extent));

                total *= uint64_t(extent);
            }
        }

        if (order == MatrixRangeRandomized && total > 1) {
            uint64_t rng = seed;

            start_rank = splitmix64_next(rng) % total;

            // Pick a stride coprime with total.
            do {
                stride = splitmix64_next(rng) % total;
            } while (stride == 0 || gcd_u64(stride, total) != 1);
        }
        else {
            start_rank = 0;
            stride = 1;
        }
    }

    void decode_rank(uint64_t rank) {
        // Reset fixed and constrained entries.
        // This also handles entries with min == max.
        const int N = curM.m * curM.n;
        for (int i = 0; i < N; ++i) {
            curM[i] = minM[i];
        }

        // Mixed-radix decoding.
        // free_index[0] is the fastest-changing coordinate.
        for (size_t k = 0; k < free_index.size(); ++k) {
            const uint64_t base = uint64_t(free_extent[k]);
            const int digit = int(rank % base);
            rank /= base;

            const int raw_value = free_min[k] + digit;
            curM[free_index[k]] = T{ raw_value };
        }
    }

    struct iterator {
        MatrixRange<F>* range;
        uint64_t step;
        uint64_t rank;

        iterator(MatrixRange<F>* range_,
            uint64_t step_,
            uint64_t rank_)
            : range(range_),
            step(step_),
            rank(rank_) {
        }

        Matrix<F>& operator*() {
            range->decode_rank(rank);
            return range->curM;
        }

        Matrix<F>* operator->() {
            range->decode_rank(rank);
            return &range->curM;
        }

        iterator& operator++() {
            ++step;

            if (range->total > 1) {
                // rank = (rank + stride) % total, without overflow.
                const uint64_t remaining = range->total - rank;

                if (range->stride >= remaining) {
                    rank = range->stride - remaining;
                }
                else {
                    rank += range->stride;
                }
            }

            return *this;
        }

        bool operator!=(const iterator& other) const {
            return step != other.step;
        }
    };

    iterator begin() {
        return iterator(this, 0, start_rank);
    }

    iterator end() {
        return iterator(this, total, 0);
    }
};

enum DiscrepancyPlotMetric {
    DISCREPANCY_GENERALIZED_L2,
    DISCREPANCY_STAR
};

extern inline const char* discrepancy_metric_name(DiscrepancyPlotMetric metric);
extern inline double safe_log10(double x);
extern inline void svg_text(FILE* f, double x, double y, const char* text, int font_size = 12, const char* anchor = "middle", const char* fill = "#333333");
extern inline void svg_line(FILE* f, double x1, double y1, double x2, double y2, const char* stroke, double width = 1.0, const char* extra = "");
extern inline void svg_circle(FILE* f, double cx, double cy, double r, const char* fill, const char* stroke = "#333333");
extern inline long long integer_power_clamped(int base, int exponent);
extern inline void format_integer_label(long long v, char* out, int out_size);

template<class F>
bool plot_discrepancy_svg(const MatrixView<F>* all_matrices,
    int n_matrices,
    double m_min,
    double m_max,
    double octave_step,
    DiscrepancyPlotMetric metric,
    const char* filename,
    int width = 900,
    int height = 560,
    bool draw_iid_reference = true) {
    assert(all_matrices != 0);
    assert(n_matrices > 0);
    assert(m_max >= m_min);
    assert(octave_step > 0.0);
    assert(filename != 0);

    enum { Q = GFCardinality<F::p, F::r>::value };

    const int dim = n_matrices;

    struct Sample {
        double m_value;
        long long n_points;
        double discrepancy;
    };

    std::vector<Sample> samples;

    long long previous_n = -1;

    for (double m = m_min; m <= m_max + 1e-12; m += octave_step) {
        const double nf = std::pow(double(Q), m);
        long long n_points = static_cast<long long>(std::floor(nf + 0.5));

        if (n_points < 1) {
            n_points = 1;
        }

        // Avoid duplicate sample counts caused by rounding.
        if (n_points == previous_n) {
            continue;
        }

        previous_n = n_points;

        std::vector<double> points;
        points.resize(static_cast<size_t>(n_points) * static_cast<size_t>(dim));

        get_points<F>(
            all_matrices,
            n_matrices,
            n_points,
            &points[0]
        );

        double disc = 0.0;

        if (metric == DISCREPANCY_GENERALIZED_L2) {
            disc = generalized_l2_discrepancy(&points[0],
                static_cast<int>(n_points),
                dim);
        }
        else {
            disc = star_discrepancy(&points[0],
                static_cast<int>(n_points),
                dim);
        }

        Sample s;
        s.m_value = std::log(double(n_points)) / std::log(double(Q));
        s.n_points = n_points;
        s.discrepancy = disc;

        samples.push_back(s);
    }

    if (samples.empty()) {
        return false;
    }

    double xmin = samples.front().m_value;
    double xmax = samples.back().m_value;

    if (xmax <= xmin) {
        xmax = xmin + 1.0;
    }

    double ymin = samples[0].discrepancy;
    double ymax = samples[0].discrepancy;

    for (size_t i = 0; i < samples.size(); ++i) {
        ymin = std::min(ymin, samples[i].discrepancy);
        ymax = std::max(ymax, samples[i].discrepancy);
    }

    // Add an IID reference curve with slope N^{-1/2}.
    // The curve is normalized to the first discrepancy value for readability.
    std::vector<double> iid_reference(samples.size(), 0.0);

    if (draw_iid_reference && samples[0].discrepancy > 0.0) {
        const double n0 = double(samples[0].n_points);
        const double d0 = samples[0].discrepancy;

        for (size_t i = 0; i < samples.size(); ++i) {
            const double n = double(samples[i].n_points);
            iid_reference[i] = d0 * std::sqrt(n0 / n);

            ymin = std::min(ymin, iid_reference[i]);
            ymax = std::max(ymax, iid_reference[i]);
        }
    }

    const double min_positive = 1e-300;

    if (ymin <= 0.0) {
        ymin = min_positive;
    }

    if (ymax <= 0.0) {
        ymax = 1.0;
    }

    double log_ymin = safe_log10(ymin);
    double log_ymax = safe_log10(ymax);

    if (log_ymax <= log_ymin) {
        log_ymax = log_ymin + 1.0;
    }

    {
        const double pad = 0.08 * (log_ymax - log_ymin);
        log_ymin -= pad;
        log_ymax += pad;
    }

    const double left = 90.0;
    const double right = 30.0;
    const double top = 55.0;
    const double bottom = 80.0;

    const double plot_x0 = left;
    const double plot_y0 = top;
    const double plot_w = double(width) - left - right;
    const double plot_h = double(height) - top - bottom;
    const double plot_x1 = plot_x0 + plot_w;
    const double plot_y1 = plot_y0 + plot_h;

    auto map_x = [&](double x) -> double {
        return plot_x0 + (x - xmin) / (xmax - xmin) * plot_w;
        };

    auto map_y = [&](double y) -> double {
        const double ly = safe_log10(y);
        return plot_y1 - (ly - log_ymin) / (log_ymax - log_ymin) * plot_h;
        };

    FILE* f = std::fopen(filename, "w");
    if (!f) {
        return false;
    }

    std::fprintf(f,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        "version=\"1.1\" width=\"%d\" height=\"%d\" "
        "viewBox=\"0 0 %d %d\">\n",
        width, height, width, height);

    std::fprintf(f,
        "<rect x=\"0\" y=\"0\" width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n",
        width, height);

    std::fprintf(f,
        "<rect x=\"%.3f\" y=\"%.3f\" width=\"%.3f\" height=\"%.3f\" "
        "fill=\"#fbfbfb\" stroke=\"#cccccc\" stroke-width=\"1\"/>\n",
        plot_x0, plot_y0, plot_w, plot_h);

    {
        char title[512];
        std::snprintf(title, sizeof(title),
            "%s, dimension %d, GF(%d)",
            discrepancy_metric_name(metric),
            dim,
            Q);

        svg_text(f, width * 0.5, 28.0, title, 18, "middle", "#222222");
    }

    // X ticks are shown only at integer octaves, with labels equal to N = q^m.
    {
        const int m_tick_min = static_cast<int>(std::ceil(xmin - 1e-12));
        const int m_tick_max = static_cast<int>(std::floor(xmax + 1e-12));

        for (int mt = m_tick_min; mt <= m_tick_max; ++mt) {
            const double x = map_x(double(mt));

            svg_line(f, x, plot_y0, x, plot_y1, "#e0e0e0", 1.0);

            const long long n_tick = integer_power_clamped(Q, mt);

            char label[64];
            format_integer_label(n_tick, label, sizeof(label));

            svg_text(f, x, plot_y1 + 20.0, label, 12, "middle", "#333333");
        }
    }

    // Horizontal log-scale grid.
    {
        const int e_min = static_cast<int>(std::ceil(log_ymin));
        const int e_max = static_cast<int>(std::floor(log_ymax));

        for (int e = e_min; e <= e_max; ++e) {
            const double yv = std::pow(10.0, double(e));
            const double y = map_y(yv);

            svg_line(f, plot_x0, y, plot_x1, y, "#e0e0e0", 1.0);

            char label[64];
            std::snprintf(label, sizeof(label), "1e%d", e);
            svg_text(f, plot_x0 - 12.0, y + 4.0, label, 12, "end", "#333333");
        }
    }

    svg_line(f, plot_x0, plot_y1, plot_x1, plot_y1, "#333333", 1.4);
    svg_line(f, plot_x0, plot_y0, plot_x0, plot_y1, "#333333", 1.4);

    svg_text(f, plot_x0 + plot_w * 0.5, height - 28.0,
        "Number of points N", 14, "middle", "#222222");

    std::fprintf(f,
        "<text x=\"24\" y=\"%.3f\" "
        "font-family=\"Arial, Helvetica, sans-serif\" "
        "font-size=\"14\" text-anchor=\"middle\" fill=\"#222222\" "
        "transform=\"rotate(-90 24 %.3f)\">%s</text>\n",
        plot_y0 + plot_h * 0.5,
        plot_y0 + plot_h * 0.5,
        discrepancy_metric_name(metric));

    const char* curve_color = "#5B8DEF";
    const char* point_color = "#A7C7F9";
    const char* iid_color = "#D36C6C";

    // IID reference line.
    if (draw_iid_reference && samples.size() >= 2 && samples[0].discrepancy > 0.0) {
        std::fprintf(f,
            "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"2.0\" "
            "stroke-dasharray=\"7 5\" stroke-linejoin=\"round\" "
            "stroke-linecap=\"round\" points=\"",
            iid_color);

        for (size_t i = 0; i < samples.size(); ++i) {
            const double x = map_x(samples[i].m_value);
            const double y = map_y(iid_reference[i]);
            std::fprintf(f, "%.3f,%.3f ", x, y);
        }

        std::fprintf(f, "\"/>\n");
    }

    // Measured discrepancy curve.
    if (samples.size() >= 2) {
        std::fprintf(f,
            "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"2.2\" "
            "stroke-linejoin=\"round\" stroke-linecap=\"round\" points=\"",
            curve_color);

        for (size_t i = 0; i < samples.size(); ++i) {
            const double x = map_x(samples[i].m_value);
            const double y = map_y(samples[i].discrepancy);
            std::fprintf(f, "%.3f,%.3f ", x, y);
        }

        std::fprintf(f, "\"/>\n");
    }

    // Data points with native SVG tooltips.
    for (size_t i = 0; i < samples.size(); ++i) {
        const double x = map_x(samples[i].m_value);
        const double y = map_y(samples[i].discrepancy);

        std::fprintf(f,
            "<circle cx=\"%.3f\" cy=\"%.3f\" r=\"4.0\" "
            "fill=\"%s\" stroke=\"#2f5fa7\" stroke-width=\"1\">\n",
            x, y, point_color);

        std::fprintf(f,
            "<title>N=%lld, log_%d(N)=%.6g, discrepancy=%.12g</title>\n",
            samples[i].n_points,
            Q,
            samples[i].m_value,
            samples[i].discrepancy);

        std::fprintf(f, "</circle>\n");
    }

    // Compact legend.
    {
        const double lx = plot_x1 - 255.0;
        const double ly = plot_y0 + 16.0;

        std::fprintf(f,
            "<rect x=\"%.3f\" y=\"%.3f\" width=\"240\" height=\"55\" "
            "rx=\"6\" ry=\"6\" fill=\"#ffffff\" stroke=\"#d0d0d0\"/>\n",
            lx, ly);

        svg_circle(f, lx + 15.0, ly + 18.0, 4.0, point_color, "#2f5fa7");
        svg_line(f, lx + 27.0, ly + 18.0, lx + 57.0, ly + 18.0, curve_color, 2.2);

        svg_text(f, lx + 66.0, ly + 22.0,
            discrepancy_metric_name(metric),
            12, "start", "#333333");

        if (draw_iid_reference) {
            svg_line(f, lx + 10.0, ly + 39.0, lx + 57.0, ly + 39.0,
                iid_color, 2.0, "stroke-dasharray=\"7 5\"");

            svg_text(f, lx + 66.0, ly + 43.0,
                "IID uniform, slope N^{-1/2}",
                12, "start", "#333333");
        }
    }

    std::fprintf(f, "</svg>\n");
    std::fclose(f);

    return true;
}



template<int DIM>
inline double gl2_kernel_soa_dim(const double* Y,
    int npts,
    int i,
    int k) {
    // Y[d * npts + i] = 2 - x_{i,d}
    // Kernel = prod_d min(Y_i,d, Y_k,d)

    double prod = 1.0;

    for (int d = 0; d < DIM; ++d) {
        const double yi = Y[d * npts + i];
        const double yk = Y[d * npts + k];
        prod *= yi < yk ? yi : yk;
    }

    return prod;
}

extern inline double sanitize_unit_coord(double x);

template<int DIM>
double generalized_l2_discrepancy_squared_exact_dim(
    const double* points,
    int npts,
    int block_size = 64
) {
    assert(points != 0);
    assert(npts > 0);
    assert(DIM > 0);
    assert(block_size > 0);

    const double invN = 1.0 / double(npts);

    // Store Y in structure-of-arrays layout:
    // Y[d * npts + i] = 2 - points[i * DIM + d]
    std::vector<double> Y(size_t(DIM) * size_t(npts));

#pragma omp parallel for
    for (int i = 0; i < npts; ++i) {
        for (int d = 0; d < DIM; ++d) {
            const double x = points[i * DIM + d];
            Y[d * npts + i] = 2.0 - x;
        }
    }

    // term0 = (4/3)^DIM
    long double term0 = 1.0L;
    for (int d = 0; d < DIM; ++d) {
        term0 *= 4.0L / 3.0L;
    }

    // term1 = sum_i prod_d (3 - x_i,d^2) / 2
    long double term1_sum = 0.0L;

#pragma omp parallel for reduction(+:term1_sum)
    for (int i = 0; i < npts; ++i) {
        long double prod = 1.0L;

        for (int d = 0; d < DIM; ++d) {
            const double y = Y[d * npts + i];
            const double x = 2.0 - y;
            prod *= 0.5L * (3.0L - double(x) * double(x));
        }

        term1_sum += prod;
    }

    // Diagonal part of term2:
    // sum_i prod_d (2 - x_i,d)
    long double term2_diag = 0.0L;

#pragma omp parallel for reduction(+:term2_diag)
    for (int i = 0; i < npts; ++i) {
        long double prod = 1.0L;

        for (int d = 0; d < DIM; ++d) {
            prod *= Y[d * npts + i];
        }

        term2_diag += prod;
    }

    // Upper triangular part of term2.
    const int nb = (npts + block_size - 1) / block_size;

    long double term2_upper = 0.0L;

#pragma omp parallel for schedule(dynamic, 1) reduction(+:term2_upper)
    for (int bi = 0; bi < nb; ++bi) {
        const int i0 = bi * block_size;
        const int i1 = std::min(i0 + block_size, npts);

        long double local = 0.0L;

        for (int bj = bi; bj < nb; ++bj) {
            const int k0 = bj * block_size;
            const int k1 = std::min(k0 + block_size, npts);

            if (bi == bj) {
                // Same block: only i < k.
                for (int i = i0; i < i1; ++i) {
                    for (int k = i + 1; k < i1; ++k) {
                        local += gl2_kernel_soa_dim<DIM>(Y.data(), npts, i, k);
                    }
                }
            }
            else {
                // Different blocks: all pairs are upper-triangular.
                for (int i = i0; i < i1; ++i) {
                    for (int k = k0; k < k1; ++k) {
                        local += gl2_kernel_soa_dim<DIM>(Y.data(), npts, i, k);
                    }
                }
            }
        }

        term2_upper += local;
    }

    const double term2_sum =
        term2_diag + 2.0L * term2_upper;

    const double d2 =
        term0
        - 2.0L * invN * term1_sum
        + double(invN) * double(invN) * term2_sum;

    return (d2 < 0.0 && d2 > -1e-12) ? 0.0 : d2;
}



template<int DIM>
double generalized_l2_discrepancy_compile_time_dimension(const double* points,
    int npts,
    int block_size = 64) {
    return std::sqrt(
        generalized_l2_discrepancy_squared_exact_dim<DIM>(
            points,
            npts,
            block_size
        )
    );
}


typedef Field<2, 1> GF2;
typedef Field<3, 1> GF3;
typedef Field<2, 2> GF4;
typedef Field<5, 1> GF5;
typedef Field<7, 1> GF7;
typedef Field<2, 3> GF8;
typedef Field<3, 2> GF9;
typedef Field<11, 1> GF11;
typedef Field<13, 1> GF13;
typedef Field<2, 4> GF16;