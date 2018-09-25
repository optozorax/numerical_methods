/* 4-й пункт тестирования. */
#include <sstream>
#include <cstdlib>
#include "matrix.h"
#include "sparse.h"

int main() {
	Matrix a_d;
	generateTestMatrix(10, 4, a_d);
	a_d.saveToFile("test1/a.txt");

	Vector x;
	x.generate(10);

	Matrix x1, f1;
	x.toDenseMatrix(x1, true);
	mul(a_d, x1, f1);
	f1.saveToFile("test1/f.txt");

	for (int k = 0; k < std::numeric_limits<real>::digits10; ++k) {
		MatrixProfileSymmetric a(a_d);
		a.getDiagonalElement(0) += pow(10.0, -k);

		Vector f;
		mul(a, x, f);

		solveSLAE_by_LDL(a, f);

		std::stringstream sout;
		#ifdef ALL_FLOAT
		sout << "f";
		#endif
		#ifdef ALL_DOUBLE
		sout << "d";
		#endif
		#ifdef ALL_FLOAT_WITH_DOUBLE
		sout << "fd";
		#endif
		sout << "_" << k+1 << ".txt";

		Vector sub;
		f.saveToFile("test1/result" + sout.str());
		f.negate();
		sum(x, f, sub);
		sub.saveToFile("test1/sub" + sout.str());
	}
}