#include "tms_lib.h"
#include <vector>

int main() {

	typedef GF5 F;
	typedef F::T T;


	T a{ 1 };
	T b{ 2 };
	auto c = a* F::div(a, b);
	std::cout << c << std::endl;


	const MatrixView<F>& P = PascalPowIntMatrix<F /*GF4*/, 9 /*size*/, 1 /*power*/>::value();
	Matrix<F> inv = P.get_inverse();
	const MatrixView<F>& P3I = PascalPowIntMatrix<F, 9, 3 >::value();
	Matrix<F> P3F = PascalTranslateFieldMatrix<F, 9, 3 >::value();

	Matrix<F> triu(4, 4);
	triu.set_zero();
	for (int i = 0; i < 4; i++) {
		for (int j = i; j < 4; j++) {
			triu[i] = T{ 1 };
		}
	}

	Matrix<F> prod(15, 15);
	tensor_product_truncated(P, triu.view(), prod.view(), 15, 15);

	std::cout << "P = " << std::endl << P << std::endl;

	std::cout << "P3I = " << std::endl << P3I << std::endl;

	std::cout << "P3F = " << std::endl << P3F << std::endl;

	std::cout << "inv(P) = " << std::endl << inv << std::endl;

	Matrix<F> Id(9, 9);
	Id.set_id();

	std::vector<MatrixView<F> > tocombine(3);
	tocombine[1] = P;
	tocombine[2] = P3F.view();
	tocombine[0] = Id.view();
	int choice[3] = { 4,3,2 };

	Matrix<F> result(9, 9);
	combine_matrices(&tocombine[0], choice, 3, 9, result.view());

	std::cout << "mat = " <<std::endl << result << std::endl;

	bool ret = is_t0_progressive(&tocombine[0], 3);

	std::cout << "ret = " << (ret?std::string("progressive"): std::string("not progressive"))<< std::endl;
	std::cout << "rank = " << P.rank() << std::endl;

	Matrix<F> Plr = P; // destroy 2 ranks by 1 linear combination of rows, and 1 linear combination of columns
	for (int i = 0; i < 9; i++) {
		Plr[4 * 9 + i] = Plr[3 * 9 + i] + Plr[5 * 9 + i];
	}
	for (int i = 0; i < 9; i++) {
		Plr[i * 9 + 2] = Plr[i * 9 + 1] - T{ 2 }*Plr[i * 9 + 6];
	}

	std::cout << "rank7 = " << Plr.rank() << std::endl;

	// 3 matrices from Victor. t={0, 1, 1, 1, 1, 2, 3, 2}
	int debugs[3][8][8] = { {{1, 1, 2, 2, 3, 1, 3, 2}, {0, 2, 3, 4, 1, 3, 2, 2}, {0, 0, 2, 0, 0, 4, 1,
 4}, {0, 0, 0, 1, 3, 4, 1, 1}, {0, 0, 0, 0, 3, 2, 4, 0}, {0, 0, 0, 0, 0, 3,
 0, 1}, {0, 0, 0, 0, 0, 0, 1, 1}, {0, 0, 0, 0, 0, 0, 0, 2}}, {{4, 3, 2, 4, 1,
 2, 0, 4}, {0, 4, 2, 1, 2, 1, 3, 1}, {0, 0, 1, 2, 4, 3, 3, 0}, {0, 0, 0, 4,
 1, 1, 2, 4}, {0, 0, 0, 0, 1, 0, 3, 0}, {0, 0, 0, 0, 0, 1, 4, 0}, {0, 0, 0,
 0, 0, 0, 4, 1}, {0, 0, 0, 0, 0, 0, 0, 4}}, {{4, 3, 4, 4, 4, 2, 1, 4}, {0, 3,
 3, 2, 1, 4, 2, 3}, {0, 0, 1, 4, 3, 2, 1, 2}, {0, 0, 0, 3, 3, 4, 1, 0}, {0,
 0, 0, 0, 3, 4, 3, 1}, {0, 0, 0, 0, 0, 2, 0, 1}, {0, 0, 0, 0, 0, 0, 2, 0},
 {0, 0, 0, 0, 0, 0, 0, 1}} };

	std::vector<MatrixView<F> > tocombine2(3);
	tocombine2[0] = MatrixView<F>(&debugs[0][0][0], 8, 8);
	tocombine2[1] = MatrixView<F>(&debugs[1][0][0], 8, 8);
	tocombine2[2] = MatrixView<F>(&debugs[2][0][0], 8, 8);

	std::vector<int> all_t = t_values(&tocombine2[0], 3);


	return 0;
};