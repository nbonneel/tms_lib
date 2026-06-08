#include "tms_lib.h"
#include <vector>
#include <chrono>
#include <list>

void test_arithmetic() {
	typedef GF5 F;
	typedef F::T T;

	T a{ 3 };
	T b{ 2 };
	T c{ 4 };
	auto d = a * F::div(a, b) + c;
	std::cout << d << std::endl;
}


void testPascal() {

	typedef GF4 F;
	typedef F::T T;

	const MatrixView<F>& P = PascalPowIntMatrix<F /*GF4*/, 9 /*size*/, 1 /*power*/>::value();
	const MatrixView<F>& P3I = PascalPowIntMatrix<F, 9, 3 >::value();
	const Matrix<F> P3F = PascalTranslateFieldMatrix<F, 9, 3 >::value();

	std::cout << "P = " << std::endl << P << std::endl;

	std::cout << "PascalPowIntMatrix(3) = " << std::endl << P3I << std::endl;   // the one that results in Pascal^3

	std::cout << "PascalTranslateFieldMatrix(3) = " << std::endl << P3F << std::endl; // the one used in Sobol' (amounts to translating polynomial by 3 " f(x+3))

	std::cout << "P^3 = " << std::endl << P*P*P << std::endl;
}

void testSobolCharacteristicMatrix() {

	typedef GF8 F;
	typedef F::T T;

	Matrix<F> poly1(1, 5, { 1, 4, 0, 0, 0 });  // x^5 + 4x + 1 (monic polynomial : x^5 is necessarily 1)
	Matrix<F> init1(5, 5, { 1, 1, 0, 2, 1,   0, 1, 1, 3, 2,   0, 0, 1, 1, 1,   0, 0, 0, 1, 0,  0, 0, 0, 0, 1 });

	SobolMatrix<F> S1(20, poly1.view(), init1.view());
	std::cout << "S1 : " << std::endl;
	std::cout << S1 << std::endl;


	Matrix<F> poly2(1, 5, { 2, 4, 0, 0, 0 });  // x^5 + 4x + 2 
	Matrix<F> init2(5, 5, { 1, 3, 2, 5, 1,   0, 1, 4, 1, 0,   0, 0, 1, 2, 1,   0, 0, 0, 1, 3,  0, 0, 0, 0, 1 });

	SobolMatrix<F> S2(20, poly2.view(), init2.view());
	std::cout << "S2 : " << std::endl;
	std::cout << S2 << std::endl;

	std::cout << "characteristic : " << std::endl;
	std::cout << S2 * S1.get_inverse() << std::endl;

}

void test_tensor_product() {

	typedef GF5 F;
	typedef F::T T;

	Matrix<F> triu(4, 4);
	triu.set_zero();
	for (int i = 0; i < 4; i++) {
		for (int j = i; j < 4; j++) {
			triu[i] = T{ 1 };
		}
	}

	const MatrixView<F>& P = PascalPowIntMatrix<F, 9>::value();

	Matrix<F> full_product(36, 36);
	tensor_product(P, triu.view(), full_product.view());

	Matrix<F> truncated_product(15, 15);
	tensor_product(P, triu.view(), truncated_product.view(), 15, 15);

	std::cout << "full_product = " << std::endl << full_product << std::endl;
	std::cout << "truncated_product = " << std::endl << truncated_product << std::endl;

}

void test_t0() {

	typedef GF5 F;
	typedef F::T T;

	Matrix<F> Id(9, 9);
	Id.set_id();

	const MatrixView<F>& P = PascalPowIntMatrix<F, 9>::value();
	const MatrixView<F>& P3F = PascalPowIntMatrix<F, 9, 3>::value();

	std::vector<MatrixView<F> > sequence(3);
	sequence[0] = Id.view();
	sequence[1] = P;
	sequence[2] = P3F;
	

	bool ret = is_t0_progressive(&sequence[0], 3);

	std::cout << "ret = " << (ret ? std::string("progressive") : std::string("not progressive")) << std::endl;

}

void test_rank() {

	typedef GF5 F;
	typedef F::T T;

	Matrix<F> Id(9, 9);
	Id.set_id();

	const MatrixView<F>& P = PascalPowIntMatrix<F, 9>::value();
	const MatrixView<F>& P3F = PascalPowIntMatrix<F, 9, 3>::value();

	std::vector<MatrixView<F> > matrices_to_combine(3);
	matrices_to_combine[0] = Id.view();
	matrices_to_combine[1] = P;
	matrices_to_combine[2] = P3F;
	
	int row_choice[3] = { 4,3,2 };

	Matrix<F> combination(9, 9);
	combine_matrices(&matrices_to_combine[0], row_choice, matrices_to_combine.size(), combination.cols(), combination.view());

	std::cout << "combination = " << std::endl << combination << std::endl;

	std::cout << "combination rank = " << combination.rank() << std::endl;

	Matrix<F> combination_lr = combination.view(); // destroy 2 ranks by 1 linear combination of rows, and 1 linear combination of columns
	for (int i = 0; i < 9; i++) {
		combination_lr[4 * 9 + i] = combination_lr[3 * 9 + i] + combination_lr[5 * 9 + i];
	}
	for (int i = 0; i < 9; i++) {
		combination_lr[i * 9 + 2] = combination_lr[i * 9 + 1] - T{ 2 }*combination_lr[i * 9 + 6];
	}

	std::cout << "combination_lr rank = " << combination_lr.rank() << std::endl;
}

void test_t_values() {

	typedef GF5 F;
	typedef F::T T;

	// 3 matrices from Victor. t={0, 1, 1, 1, 1, 2, 3, 2}
	/*T matrices[3][8][8] = {
		{{1, 1, 2, 2, 3, 1, 3, 2}, {0, 2, 3, 4, 1, 3, 2, 2}, {0, 0, 2, 0, 0, 4, 1, 4}, {0, 0, 0, 1, 3, 4, 1, 1},
		 {0, 0, 0, 0, 3, 2, 4, 0}, {0, 0, 0, 0, 0, 3, 0, 1}, {0, 0, 0, 0, 0, 0, 1, 1}, {0, 0, 0, 0, 0, 0, 0, 2}},
		{{4, 3, 2, 4, 1, 2, 0, 4}, {0, 4, 2, 1, 2, 1, 3, 1}, {0, 0, 1, 2, 4, 3, 3, 0}, {0, 0, 0, 4, 1, 1, 2, 4},
		 {0, 0, 0, 0, 1, 0, 3, 0}, {0, 0, 0, 0, 0, 1, 4, 0}, {0, 0, 0, 0, 0, 0, 4, 1}, {0, 0, 0, 0, 0, 0, 0, 4}},
		{{4, 3, 4, 4, 4, 2, 1, 4}, {0, 3, 3, 2, 1, 4, 2, 3}, {0, 0, 1, 4, 3, 2, 1, 2}, {0, 0, 0, 3, 3, 4, 1, 0},
		 {0, 0, 0, 0, 3, 4, 3, 1}, {0, 0, 0, 0, 0, 2, 0, 1}, {0, 0, 0, 0, 0, 0, 2, 0}, {0, 0, 0, 0, 0, 0, 0, 1}} };

	std::vector<MatrixView<F> > sequence(3);
	sequence[0] = MatrixView<F>(&matrices[0][0][0], 8, 8);
	sequence[1] = MatrixView<F>(&matrices[1][0][0], 8, 8);
	sequence[2] = MatrixView<F>(&matrices[2][0][0], 8, 8);*/


	Matrix<F> poly1(1, 5, { 4, 2, 4, 0, 2 });
	Matrix<F> init1(5, 5, { 1, 1, 0, 2, 1,   0, 1, 1, 3, 2,   0, 0, 1, 1, 1,   0, 0, 0, 1, 0,  0, 0, 0, 0, 1 });
	SobolMatrix<F> S1(40, poly1.view(), init1.view()); // 40x40 Sobol matrices

	Matrix<F> poly2(1, 5, { 1, 2, 3, 3, 4 });
	Matrix<F> init2(5, 5, { 1, 1, 3, 2, 1,   0, 1, 2, 3, 2,   0, 0, 1, 2, 1,   0, 0, 0, 1, 4,  0, 0, 0, 0, 1 });
	SobolMatrix<F> S2(40, poly2.view(), init2.view());

	Matrix<F> poly3(1, 3, { 3, 4, 0 });
	Matrix<F> init3(3, 3, { 1, 2, 2, 0, 4, 1, 0, 0, 1 });
	SobolMatrix<F> S3(40, poly3.view(), init3.view());

	Matrix<F> poly4(1, 4, { 3, 0, 3, 4 });
	Matrix<F> init4(4, 4, { 1, 2, 2, 3,   0, 4, 1, 2,   0, 0, 1, 1,  0, 0, 0, 1 });
	SobolMatrix<F> S4(40, poly4.view(), init4.view());

	std::vector<MatrixView<F> > sequence(4);
	sequence[0] = S1.view();
	sequence[1] = S2.view();
	sequence[2] = S3.view();
	sequence[3] = S4.view();

	auto time_point1 = std::chrono::high_resolution_clock::now();
	std::vector<int> all_t = t_values_naive(&sequence[0], sequence.size());
	auto time_point2 = std::chrono::high_resolution_clock::now();


	std::vector<int> all_t_v2 = t_values(&sequence[0], sequence.size());
	auto time_point3 = std::chrono::high_resolution_clock::now();

	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>((time_point2 - time_point1)).count() << " ms, t (naive) : " << std::endl;
	for (int i = 0; i < all_t.size(); i++)
		std::cout << all_t[i] << " ";
	std::cout << std::endl;


	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>((time_point3 - time_point2)).count() << " ms, t ([Marion et al. 2020]): " << std::endl;
	for (int i = 0; i < all_t_v2.size(); i++)
		std::cout << all_t_v2[i] << " ";
	std::cout << std::endl;

}

// lists all upper triangular matrices with exactly 1 in the diagonal, and at least 1 in the first row 
void test_list_matrices() {
	typedef GF3 F;
	typedef F::T T;

	const int mat_size = 3;
	Matrix<F> min_matrix(mat_size, mat_size);
	Matrix<F> max_matrix(mat_size, mat_size);
	min_matrix.set_id();
	for (int i = 0; i < mat_size; i++) {
		min_matrix[i] = 1;
	}

	max_matrix.set_id();
	for (int i = 0; i < mat_size; i++) {
		for (int j = i+1; j < mat_size; j++) {
			max_matrix[i * mat_size + j] = 2;
		}
	}

	MatrixRange<F> range_fully_randomized(min_matrix.view(), max_matrix.view(), MatrixRangeGlobalRandomized, 123456789);

	long long num_matrix = 0;
	for (auto it = range_fully_randomized.begin(); it != range_fully_randomized.end(); ++it) {
		MatrixView<F> A = it->view();
		std::cout << "Matrix " << num_matrix << std::endl;
		std::cout << A << std::endl;
		num_matrix++;
	}


	MatrixRange<F> range_DFS_randomized(
		min_matrix.view(),
		max_matrix.view(),
		MatrixRangeDFSValuePermuted,
		12345
	);

	num_matrix = 0;
	// parallel exploration of the matrices
#pragma omp parallel
	{
		Matrix<F> local(range_DFS_randomized.rows(), range_DFS_randomized.cols());

#pragma omp for schedule(dynamic, 1)
		for (long long s = 0; s < (long long)range_DFS_randomized.size(); ++s) {
			range_DFS_randomized.decode_step(uint64_t(s), local.view());

			MatrixView<F> A = local.view();
			int rank = A.rank();
#pragma omp critical
			{
				std::cout << "Matrix " << num_matrix << std::endl;
				std::cout << A << std::endl;
				std::cout << "rank = " << rank << std::endl; // ok, all matrices are full rank since triangular with positive diagonal
			}

#pragma atomic
			num_matrix++;
		}
	}
}

void test_draw_points() {

	typedef GF4 F;
	typedef F::T T;

	// theses Sobol' polynomials work nicely in GF3, but produce funny points in GF4
	Matrix<F> poly1(1, 3, { 1, 2, 0 });
	Matrix<F> init1(3, 3, { 1, 1, 1, 0, 1, 1, 0, 0, 1 });
	SobolMatrix<F> S1(10, poly1.view(), init1.view());

	Matrix<F> poly2(1, 3, { 2, 2, 0 });
	Matrix<F> init2(3, 3, { 1, 0, 1,  0, 1, 2, 0, 0, 1 });
	SobolMatrix<F> S2(10, poly2.view(), init2.view());

	std::vector<MatrixView<F> > m;
	m.push_back(S1.view());
	m.push_back(S2.view());
	
	draw_2D_points(S1.view(), S2.view(), std::pow(std::pow((int)F::p, (int)F::r), 7), "points.svg");

	S1.save_svg("matrix1.svg");
}

void test_discrepancy() {

	typedef GF5 F;
	typedef F::T T;

	const int n_pts = std::pow(std::pow((double)F::p, F::r), 4); //4.5 => 20 s

	T inits[4][5][5] = { 
	  { {1,     2,      3,      1,      4},
		{0,     1,      3,      4,      3},
		{0,     0,      1,      1,      3},
		{0,     0,      0,      1,      1},
		{0,     0,      0,      0,      1}
	},
	  { { 1,     4,		 2,		 3,		 4 },
		{ 0,     1,      1,      1,      4 },
		{ 0,     0,      1,      2,      2 },
		{ 0,     0,      0,      1,      2 },
		{ 0,     0,      0,      0,      1 }
	},
	  { { 1,	 1,		 2,		 2,		 4 },
		{ 0,     1,      4,      1,      1 },
		{ 0,     0,      1,      3,      2 },
		{ 0,     0,      0,      1,      3 },
		{ 0,     0,      0,      0,      1 }
	},
	  { { 1,	 3,		 3,		 4,		 4 },
		{ 0,     1,      2,      4,      2 },
		{ 0,     0,      1,      4,      3 },
		{ 0,     0,      0,      1,      4 },
		{ 0,     0,      0,      0,      1 }
	}
	};

	Matrix<F> poly(1, 5, { 1, 4, 0, 0, 0 } );
	MatrixView<F> init(&inits[0][0][0], 5, 5);

	std::vector< Matrix<F> > matrices(4);
	std::vector< MatrixView<F> > matrices_view(4);
	for (int i = 0; i < 4; i++) {
		poly[0] = i + 1;
		init.values = &inits[i][0][0];
		
		SobolMatrix<F> S(15, poly.view(), init);
		matrices[i] = S.view(); // copy
		matrices_view[i] = matrices[i].view();
	}


	std::vector<double> pts = get_points<F>(&matrices_view[0], matrices_view.size(), n_pts);
	double stardisc = star_discrepancy(&pts[0], n_pts, matrices_view.size());
	

	print_point_range(&pts[0], n_pts, matrices_view.size());
	std::cout << "n_pts: "<<n_pts<<", star discrepancy="<< stardisc << std::endl;

	double gl2disc = generalized_l2_discrepancy(&pts[0], n_pts, matrices_view.size());
	std::cout << "n_pts: " << n_pts << ", generalized l2 discrepancy=" << gl2disc << std::endl;

	plot_discrepancy_svg(&matrices_view[0], matrices_view.size(), 1, 6, 1 / 4., DISCREPANCY_GENERALIZED_L2, "plot_disc.svg");
}


void test_owen() {


	typedef GF5 F;
	typedef F::T T;

	const int m = 7;
	const int dim = 3;
	const int base = std::pow((double)F::p, F::r);
	const int n_pts = std::pow(base, m);

	// 3 matrices from Victor. t={0, 1, 1, 1, 1, 2, 3, 2}
	T matrices[3][8][8] = {
		{{1, 1, 2, 2, 3, 1, 3, 2}, {0, 2, 3, 4, 1, 3, 2, 2}, {0, 0, 2, 0, 0, 4, 1, 4}, {0, 0, 0, 1, 3, 4, 1, 1},
		 {0, 0, 0, 0, 3, 2, 4, 0}, {0, 0, 0, 0, 0, 3, 0, 1}, {0, 0, 0, 0, 0, 0, 1, 1}, {0, 0, 0, 0, 0, 0, 0, 2}},
		{{4, 3, 2, 4, 1, 2, 0, 4}, {0, 4, 2, 1, 2, 1, 3, 1}, {0, 0, 1, 2, 4, 3, 3, 0}, {0, 0, 0, 4, 1, 1, 2, 4},
		 {0, 0, 0, 0, 1, 0, 3, 0}, {0, 0, 0, 0, 0, 1, 4, 0}, {0, 0, 0, 0, 0, 0, 4, 1}, {0, 0, 0, 0, 0, 0, 0, 4}},
		{{4, 3, 4, 4, 4, 2, 1, 4}, {0, 3, 3, 2, 1, 4, 2, 3}, {0, 0, 1, 4, 3, 2, 1, 2}, {0, 0, 0, 3, 3, 4, 1, 0},
		 {0, 0, 0, 0, 3, 4, 3, 1}, {0, 0, 0, 0, 0, 2, 0, 1}, {0, 0, 0, 0, 0, 0, 2, 0}, {0, 0, 0, 0, 0, 0, 0, 1}} };

	std::vector<MatrixView<F> > sequence(3);
	sequence[0] = MatrixView<F>(&matrices[0][0][0], 8, 8);
	sequence[1] = MatrixView<F>(&matrices[1][0][0], 8, 8);
	sequence[2] = MatrixView<F>(&matrices[2][0][0], 8, 8);

	std::vector<int> t_vals = t_values(&sequence[0], dim);

	std::vector<double> pts = get_points<F>(&sequence[0], dim, n_pts);

	OwenTreeND tree = make_random_owen_tree_nd(dim, base, m+3, 12345);
	std::vector<double> scrambled(n_pts * dim);
	apply_owen_permutation_real(&pts[0], &scrambled[0], n_pts, dim, m+3, tree);

	int t_before = t_factor_pointset(&pts[0], n_pts, dim, base);
	int t_after = t_factor_pointset(&scrambled[0], n_pts, dim, base);
	int t_matrix = t_vals[m-1];
	std::cout << "t value according to matrix : " << t_matrix << std::endl;
	std::cout << "t value according to point set : " << t_before << std::endl;
	std::cout << "t value according to scrambled point set : " << t_after << std::endl;
	double gl2disc_before = generalized_l2_discrepancy(&pts[0], n_pts, dim);
	std::cout << "generalized l2 discrepancy before scrambling: " << gl2disc_before << std::endl;
	double gl2disc_after = generalized_l2_discrepancy(&scrambled[0], n_pts, dim);
	std::cout << "generalized l2 discrepancy after scrambling: " << gl2disc_after << std::endl;

}


int main() {

	test_arithmetic();

	testPascal();

	testSobolCharacteristicMatrix();

	test_tensor_product();

	test_t0();

	test_rank();

	test_t_values();

	test_list_matrices();

	test_draw_points();	

	test_discrepancy();	

	test_owen();

	return 0;
};