#define TMS_USE_CUDA

#include "tms_lib.h"
#include <vector>
#include <chrono>
#include <list>
#include <sstream>

void test_arithmetic() {
	typedef GF5 F;
	typedef typename F::T T;

	T a{ 3 };
	T b{ 2 };
	T c{ 4 };
	T d = a * F::div(a, b) + c;
	std::cout << d << std::endl;
	d = gf_reduce<F::p, F::r>(d);
	std::cout << d << std::endl;
}


void testPascal() {

	typedef GF4 F;
	typedef typename F::T T;

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
	typedef typename F::T T;

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
	typedef typename F::T T;

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
	typedef typename F::T T;

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
	typedef typename F::T T;

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
	typedef typename F::T T;

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
	typedef typename F::T T;

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

#pragma omp atomic
			num_matrix++;
		}
	}
}

void test_draw_points() {

	typedef GF3 F;
	typedef typename F::T T;
	const int dim = 8;

	std::vector<SobolMatrix<F> > seq;
	for (int i = 0; i < dim; i++) {
		seq.push_back(SobolMatrix<F>(SOBOL_QUADQUAD_GF3, i, 12));
	}
	
	std::vector<MatrixView<F> > seq_view;
	for (int i = 0; i < dim; i++) {
		seq_view.push_back(seq[i].view());
	}

	
	// save two dimensions
	draw_2D_points(seq_view[0], seq_view[1], std::pow(std::pow((int)F::p, (int)F::r), 7), "points.svg");

	// save matrix
	seq_view[1].save_svg("matrix1.svg");



	// saves all projection pairs
	// first the set of projections we would like to highlight and which colors
	ProjectionHighlight highlights[] = {
		ProjectionHighlight(0, 1, "#F28E2B", 2.4),

		ProjectionHighlight(0, 2, "#2CA02C", 2.4),
		ProjectionHighlight(1, 2, "#2CA02C", 2.4),

		ProjectionHighlight(0, 3, "#2CA02C", 2.4),
		ProjectionHighlight(1, 3, "#2CA02C", 2.4),
		ProjectionHighlight(2, 3, "#2CA02C", 2.4),

		ProjectionHighlight(4, 5, "#F28E2B", 2.4),

		ProjectionHighlight(4, 6, "#2CA02C", 2.4),
		ProjectionHighlight(5, 6, "#2CA02C", 2.4),

		ProjectionHighlight(4, 7, "#2CA02C", 2.4),
		ProjectionHighlight(5, 7, "#2CA02C", 2.4),
		ProjectionHighlight(6, 7, "#2CA02C", 2.4)
	};
	std::vector<double> points = get_points(&seq_view[0], dim, std::pow(std::pow((int)F::p, (int)F::r), 5));
	draw_2d_projections_svg(&points[0], dim, points.size()/dim, "projections.svg", highlights, 12);
}



void test_owen() {


	typedef GF5 F;
	typedef typename F::T T;

	const int m = 8;
	const int dim = 3;
	const int base = std::pow((int)F::p, (int) F::r);
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

	int t_before = t_value_pointset(&pts[0], n_pts, dim, base);
	int t_after = t_value_pointset(&scrambled[0], n_pts, dim, base);
	int t_matrix = t_vals[m-1];
	std::cout << "t value according to matrix : " << t_matrix << std::endl;
	std::cout << "t value according to point set : " << t_before << std::endl;
	std::cout << "t value according to scrambled point set : " << t_after << std::endl;
	auto time_point1 = std::chrono::high_resolution_clock::now();
	double gl2disc_before = generalized_l2_discrepancy(&pts[0], n_pts, dim, 128, false);
	auto time_point2 = std::chrono::high_resolution_clock::now();
	double gl2disc_before_gpu = generalized_l2_discrepancy(&pts[0], n_pts, dim, 128, true);
	auto time_point3 = std::chrono::high_resolution_clock::now();

	
	std::cout << "generalized l2 discrepancy before scrambling (CPU): " << gl2disc_before << ", time taken: "<<std::chrono::duration_cast<std::chrono::milliseconds>((time_point2 - time_point1)).count() << " ms"<< std::endl;
	std::cout << "generalized l2 discrepancy before scrambling (GPU): " << gl2disc_before_gpu << ", time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>((time_point3 - time_point2)).count() << " ms" << std::endl;

	double gl2disc_after = generalized_l2_discrepancy(&scrambled[0], n_pts, dim);
	std::cout << "generalized l2 discrepancy after scrambling: " << gl2disc_after << std::endl;

	auto time_point4 = std::chrono::high_resolution_clock::now();
	double stardisc = star_discrepancy(&scrambled[0], std::min(n_pts,2500), dim, false);
	auto time_point5 = std::chrono::high_resolution_clock::now();
	std::cout << "star discrepancy in 3D ("<< std::min(n_pts, 2300)<<" points, CPU) : "<<stardisc << ", time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>((time_point5 - time_point4)).count() << " ms" << std::endl;

	auto time_point6 = std::chrono::high_resolution_clock::now();
	stardisc = star_discrepancy(&scrambled[0], std::min(n_pts, 2300), dim, true);
	auto time_point7 = std::chrono::high_resolution_clock::now();
	std::cout << "star discrepancy in 3D (" << std::min(n_pts, 2300) << " points, GPU) : " << stardisc << ", time taken: " << std::chrono::duration_cast<std::chrono::milliseconds>((time_point7 - time_point6)).count() << " ms" << std::endl;

}


void test_plot_discrepancy() {

	typedef GF2 F;
	typedef typename F::T T;
	enum { base = GFCardinality<F::p, F::r>::value };

	const int max_m = 9;
	const int dim = 5;
	const int max_n_pts = std::pow((double)base, max_m); //4.5 => 20 s
	const int nb_owen = 32;

	
	std::vector< Matrix<GF2> > matricesJoeKuo(dim);
	std::vector< MatrixView<GF2> > matricesJoeKuo_view(dim);
	for (int i = 0; i < dim; i++) {
		SobolMatrix<GF2> sobol(SOBOL_JOE_KUO_GF2, i, 25);
		matricesJoeKuo[i] = sobol.view(); // copy
		matricesJoeKuo_view[i] = matricesJoeKuo[i].view();
	}


	std::vector<double> pts = get_points<F>(&matricesJoeKuo_view[0], matricesJoeKuo_view.size(), max_n_pts);

	DiscrepancyCurve curve_star, curve_gl2;
	std::vector<double> discrepancies;
	for (double m = 1; m <= max_m; m += 1. / 8.) {

		const int n_pts = std::pow((double)base, m);

		double avg_star_disc = 0;
		double avg_gl2_disc = 0;
		
#pragma omp parallel for reduction(+:avg_gl2_disc)	reduction(+:avg_star_disc)	
		for (int i = 0; i < nb_owen; i++) {
			std::vector<double> scrambled(n_pts* dim);
			OwenTreeND tree = make_random_owen_tree_nd(dim, base, std::ceil(m), 12345+i*123);			
			apply_owen_permutation_real(&pts[0], &scrambled[0], n_pts, dim, std::ceil(m), tree);
			padd_least_significant_digits(&scrambled[0], n_pts, dim, base, std::ceil(m), 123456+i*456);

			// since we are computing discrepancies in parallel, let's do it on the CPU.
			double stardisc = star_discrepancy(&scrambled[0], n_pts, dim, false);
			double gl2disc = generalized_l2_discrepancy(&scrambled[0], n_pts, dim, 128, false);

			avg_star_disc += stardisc;
			avg_gl2_disc += gl2disc;
			
		}
		avg_star_disc /= nb_owen;
		avg_gl2_disc /= nb_owen;

		curve_star.values.push_back(avg_star_disc);
		curve_gl2.values.push_back(avg_gl2_disc);

		//curve.show_points = false;
		curve_star.n_points.push_back(n_pts);
		curve_gl2.n_points.push_back(n_pts);
	}
	curve_star.label = "Star discrepancy";
	curve_gl2.label = "Generalized l2 discrepancy";

	std::vector< DiscrepancyCurve> curves;
	curves.push_back(curve_star);
	curves.push_back(curve_gl2);

	std::vector< ReferenceCurveSpec > refs;
	ReferenceCurveSpec r1;
	r1.label = "N^{-1/2}";
	r1.exponent = 0.5;
	r1.dashed = true;
	refs.push_back(r1);

	ReferenceCurveSpec r2;
	r2.label = "N^{-1}";
	r2.exponent = 1.0;
	r2.dashed = true;
	r2.color = "#6CA87A";
	refs.push_back(r2);

	DiscrepancyPlotOptions opt;
	opt.title = "5D Star vs. Generalized L2 discrepancies (Joe &amp; Kuo Sobol)";
	opt.x_label = "Number of points N";
	opt.y_label = "Discrepancy value";
	opt.base = base;
	opt.x_tick_octave_step = 1.0;
	opt.y_scale = PLOT_Y_LOG10;

	plot_discrepancy_curves_svg(curves, refs, opt, "star_vs_l2_discrepancy.svg");
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

	test_owen();

	test_plot_discrepancy();	// about 5 minutes on my machine.

	


	return 0;
};