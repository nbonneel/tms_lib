#include "tmslib.h"
#include <vector>

int main() {

	typedef GF8 F;
	typedef F::T T;


	T a{ 1 };
	T b{ 2 };
	auto c = a* F::div(a, b);
	std::cout << c << std::endl;


	const MatrixView<F>& P = PascalPowIntMatrix<F /*GF4*/, 9 /*size*/, 1 /*power*/>::value();
	Matrix<F> inv = P.get_inverse();
	const MatrixView<F>& P3I = PascalPowIntMatrix<F, 9, 3 >::value();
	const MatrixView<F>& P3F = PascalTranslateFieldMatrix<F, 9, 3 >::value();

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
	tocombine[2] = P3F;
	tocombine[0] = Id.view();
	int choice[3] = { 4,3,2 };

	Matrix<F> result(9, 9);
	combine_matrices(&tocombine[0], choice, 3, 9, result.view());


	std::cout << "mat = " <<std::endl << result << std::endl;

	bool ret = is_t0_progressive(&tocombine[0], 3);

	std::cout << "ret = " << (ret?std::string("progressive"): std::string("not progressive"))<< std::endl;
	std::cout << "rank = " << P.rank() << std::endl;

	Matrix<F> Plr = P;
	for (int i = 0; i < 9; i++) {
		Plr[4 * 9 + i] = Plr[3 * 9 + i] + Plr[5 * 9 + i];
	}
	for (int i = 0; i < 9; i++) {
		Plr[i * 9 + 2] = Plr[i * 9 + 1] - T{ 2 }*Plr[i * 9 + 6];
	}

	std::cout << "rank7 = " << Plr.rank() << std::endl;
	return 0;
};