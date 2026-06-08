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

void test_t_values() {

	typedef GF5 F;
	typedef F::T T;

	// 3 matrices from Victor. t={0, 1, 1, 1, 1, 2, 3, 2}
	/*int matrices[3][8][8] = { 
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
	std::vector<int> all_t = t_values_naive(&sequence[0], 4);
	auto time_point2 = std::chrono::high_resolution_clock::now();
	

	std::vector<int> all_t_v2 = t_values(&sequence[0], 4);
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
	combine_matrices(&matrices_to_combine[0], row_choice, 3, 9, combination.view());

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

	MatrixRange<F> range(min_matrix.view(), max_matrix.view());

	long long num_matrix = 0;
	for (auto it = range.begin(); it != range.end(); ++it) {
		MatrixView<F> A = it->view();
		std::cout << "Matrix " << num_matrix << std::endl;
		std::cout << A << std::endl;
		num_matrix++;
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
	std::vector<int> vals = t_values(&m[0], 2, 5);
	draw_2D_points(S1.view(), S2.view(), std::pow(std::pow((int)F::p, (int)F::r), 7), "points.svg");

	S1.save_svg("matrix1.svg");
}

void test_discrepancy() {

	typedef GF5 F;
	typedef F::T T;

	const int n_pts = std::pow(std::pow((double)F::p, F::r), 2);

	T inits[4][5][5] = { { {1,    2,      3,      1,      4},
{0,     1,      3,      4,      3},
{0,     0,      1,      1,      3},
{0,     0,      0,      1,      1},
{0,     0,      0,      0,      1}
	},
	{ { 1, 4, 2, 3, 4 },
		{ 0,     1,      1,      1,      4 },
		{ 0,     0,      1,      2,      2 },
		{ 0,     0,      0,      1,      2 },
		{ 0,     0,      0,      0,      1 }
	},
	{ { 1, 1, 2, 2, 4 },
		{ 0,     1,      4,      1,      1 },
		{ 0,     0,      1,      3,      2 },
		{ 0,     0,      0,      1,      3 },
		{ 0,     0,      0,      0,      1 }
	},
	{ { 1, 3, 3, 4, 4 },
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

	plot_discrepancy_svg(&matrices_view[0], matrices_view.size(), 1, 5, 1 / 4., DISCREPANCY_STAR, "plot_disc.svg");
}



void generate_all_mat() {

	typedef GF5 F;
	typedef F::T T;
	enum { Q = GFCardinality<F::p, F::r>::value };

	const int test_size = 8;
	const int e = Q;
	const Matrix<F> P = PascalTranslateFieldMatrix<F, Q, 1 >::value();
	std::vector<MatrixView<F> > Ppow(12);
	Ppow[0] = PascalTranslateFieldMatrix<F, Q, 1 >::value();
	Ppow[1] = PascalTranslateFieldMatrix<F, Q, 2 >::value();
	Ppow[2] = PascalTranslateFieldMatrix<F, Q, 3 >::value();
	Ppow[3] = PascalTranslateFieldMatrix<F, Q, 4 >::value();
	Ppow[4] = PascalTranslateFieldMatrix<F, Q, 5 >::value();
	Ppow[5] = PascalTranslateFieldMatrix<F, Q, 6 >::value();
	Ppow[6] = PascalTranslateFieldMatrix<F, Q, 7 >::value();
	Ppow[7] = PascalTranslateFieldMatrix<F, Q, 8 >::value();
	Ppow[8] = PascalTranslateFieldMatrix<F, Q, 9 >::value();
	Ppow[9] = PascalTranslateFieldMatrix<F, Q, 10 >::value();
	Ppow[10] = PascalTranslateFieldMatrix<F, Q, 11 >::value();
	Ppow[11] = PascalTranslateFieldMatrix<F, Q, 12 >::value();


	Matrix<F> init(e, e);

	std::vector<MatrixView<F> > high_dim_sequence(Q * 2 - 1);

	high_dim_sequence[0] = PascalTranslateFieldMatrix<F, test_size, 0 >::value();
	high_dim_sequence[1] = PascalTranslateFieldMatrix<F, test_size, 1 >::value();
	high_dim_sequence[2] = PascalTranslateFieldMatrix<F, test_size, 2 >::value();
	high_dim_sequence[3] = PascalTranslateFieldMatrix<F, test_size, 3 >::value();
	high_dim_sequence[4] = PascalTranslateFieldMatrix<F, test_size, 4 >::value();
	/*high_dim_sequence[5] = PascalTranslateFieldMatrix<F, test_size, 5 >::value();
	high_dim_sequence[6] = PascalTranslateFieldMatrix<F, test_size, 6 >::value();
	high_dim_sequence[7] = PascalTranslateFieldMatrix<F, test_size, 7 >::value();
	high_dim_sequence[8] = PascalTranslateFieldMatrix<F, test_size, 8 >::value();
	high_dim_sequence[9] = PascalTranslateFieldMatrix<F, test_size, 9 >::value();
	high_dim_sequence[10] = PascalTranslateFieldMatrix<F, test_size, 10 >::value();*/
	for (int i = 0; i < Q - 1; i++) {
		Matrix<F>* m = new Matrix<F>(test_size, test_size);
		high_dim_sequence[Q + i] = m->view();
	}

	std::vector< Matrix<F> > polynomials(Q - 1);
	for (int i = 0; i < Q - 1; i++) {
		polynomials[i] = Matrix<F>(1, Q);
		polynomials[i].set_zero();
		polynomials[i][0] = T{ i + 1 };
		polynomials[i][1] = -T{ 1 };
	}

	int best_s = 10000;

	std::list<std::pair<int, int> > position_value;
	for (int i = 1; i < Q; i++)
		position_value.push_back(std::pair<int, int>(0, i)); // if a matrix is good, then k*matrix is good as well, so just let's fix m[0][0] to 1

	const int list_size = std::min(test_size - 1, Q - 1);
	int* mult = new int[list_size];

	long long numOptions = std::pow(Q - 1, list_size);
	long long curOption = 0;
	while (!position_value.empty()) {
		std::pair<int, int> top = position_value.back();
		position_value.pop_back();
		int pos = top.first;
		if (pos < list_size) {
			mult[pos] = top.second;
		}


		if (pos == list_size - 1) {
			curOption++;



			std::vector< Matrix<F> > all_pow_inits(Q - 1);
			for (int i = 0; i < Q - 1; i++) {

				Matrix<F> pow_init = Ppow[i];

				for (int j = 0; j < list_size; j++) {

					for (int k = 0; k < e; k++) {
						pow_init[(j + 1) * e + k] = gf_div<F::p, F::r>(pow_init[(j + 1) * e + k], T{ mult[j] });
					}
					for (int k = 0; k < e; k++) {
						pow_init[k * e + j + 1] *= T{ mult[j] };
					}
				}
				pow_init.reduce();
				all_pow_inits[i] = pow_init.view();

				fill_sobol(high_dim_sequence[Q + i], polynomials[i].view(), pow_init.view());
			}

			std::vector<int> ttest1 = t_values(&high_dim_sequence[0], Q);
			std::vector<int> ttestt2 = t_values(&high_dim_sequence[Q], Q - 1);

			std::vector<int> t = t_values(&high_dim_sequence[0], Q * 2 - 1);
			int cur_s = 0;
			for (int i = 0; i < t.size(); i++) {
				cur_s += t[i];
			}
			if (cur_s < best_s) {
				best_s = cur_s;
				std::cout << "t = " << std::endl;
				for (int i = 0; i < t.size(); i++) {
					std::cout << t[i] << " ";
				}
				std::cout << std::endl;

				std::cout << "inits = " << std::endl;
				for (int i = 0; i < Q - 1; i++) {
					std::cout << all_pow_inits[i];
					std::cout << std::endl;
				}

			}

			/*for (int i = 0; i < e; i++) {
				for (int j = 0; j < e; j++) {
					fprintf(f, "%u ", init[i * e + j]);
				}
				fprintf(f, "\n");
			}
			fprintf(f, "\n");*/

			if (curOption % 1000 == 0) {
				std::cout << curOption / (double)numOptions * 100 << " %" << std::endl;
			}
			/*std::cout << "matrix " << std::endl;
			printMat(&init[0], e);
			std::cout << "-------" << std::endl;*/


		}
		else {
			pos++;

			for (int k = 1; k < Q; k++) {
				position_value.push_back(std::pair<int, int>(pos, k));
			}

		}
	}


	delete[] mult;


	std::cout << "numOptions : " << numOptions << std::endl;
}



int main() {

	test_arithmetic();

	testPascal();

	test_tensor_product();

	test_t0();

	//test_rank();

	//test_t_values();

	//testSobolCharacteristicMatrix();

	test_draw_points();

	//test_list_matrices();

	test_discrepancy();

	//generate_all_mat();

	return 0;
};