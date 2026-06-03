#pragma once
#include <type_traits>
#include <iostream>
#include <ostream>
#include <vector>


#define MAX_GF 16

extern char add_non_prime[MAX_GF + 1][MAX_GF][MAX_GF];
extern char sub_non_prime[MAX_GF + 1][MAX_GF][MAX_GF];
extern char mul_non_prime[MAX_GF + 1][MAX_GF][MAX_GF];
extern char div_non_prime[MAX_GF + 1][MAX_GF][MAX_GF];
extern char neg_non_prime[MAX_GF + 1][MAX_GF];
extern char invGalois[MAX_GF + 1][MAX_GF];

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

    typename T operator[](int i) const {
        return derived().values[i];
    }
    typename T& operator[](int i) {
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
};

template<class F>
struct MatrixView : MatrixBase<MatrixView<F>, F> {
    typedef typename F::T T;

    T* values;
    int m, n;

    MatrixView() : values(0), m(0), n(0) {}

    MatrixView(T* values_, int m_, int n_)
        : values(values_), m(m_), n(n_) {
    }
    typename T operator[](int i) const {
        return values[i];
    }
    typename T& operator[](int i) {
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

    typename F::T operator[](int i) const {
        return values[i];
    }
    typename F::T& operator[](int i) {
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

template<class F, int N, int Power>
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


template<class F, int N, int Power>
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

#include <cassert>
#include <algorithm>

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
void tensor_product(const MatrixView<F>& A,
    const MatrixView<F>& B,
    MatrixView<F>& result,
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

    for (int i = 0; i < M.m; ++i) {
        for (int j = 0; j < M.n; ++j) {
            out << M.values[i * M.n + j] << '\t';
        }
        out << '\n';
    }

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

    const MatrixBase<Derived, F>& mat = *this;
    int m = rows();
    int n = cols();

    if (m == 0 || n == 0) {
        return 0;
    }

    typedef typename F::T T;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            tmpPivoted[i * m + j] = mat[i * m + j];  // first column unchanged
        }
    }
    for (int j = 0; j < n - 1; j++) {

        if (tmpPivoted[j * m + j] == T{ 0 }) {
            int pivSwap = -1;
            for (int jj = j + 1; jj < n; jj++) {
                if (tmpPivoted[j * m + jj] != T{ 0 }) {
                    pivSwap = jj;
                    break;
                }
            }
            if (pivSwap == -1) {
                break;
            }
            else {
                for (int i = 0; i < m; i++) {
                    std::swap(tmpPivoted[i * m + j], tmpPivoted[i * m + pivSwap]);
                }
            }
        }
        for (int jj = j + 1; jj < n; jj++) {
            T inverseCoeff = F::div(tmpPivoted[j * m + jj], tmpPivoted[j * m + j]);  ;
            for (int i = 0; i < m; i++) {
                tmpPivoted[i * m + jj] = gf_reduce<F::p, F::r>(tmpPivoted[i * m + jj] - inverseCoeff * tmpPivoted[i * m + j]);
            }
        }
    }

    int rank = 0;
    for (int i = 0; i < std::min(m, n); i++) {
        if (tmpPivoted[i * m + i] != T{ 0 })
            rank++;
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


typedef Field<2, 1> GF2;
typedef Field<3, 1> GF3;
typedef Field<2, 2> GF4;
typedef Field<5, 1> GF5;
typedef Field<2, 3> GF8;
typedef Field<3, 2> GF9;
typedef Field<2, 4> GF16;