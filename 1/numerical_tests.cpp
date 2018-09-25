#include <sstream>
#include <cstdlib>
#include "matrix.h"
#include "sparse.h"

//-----------------------------------------------------------------------------
void test1(void) {
	/* 4-й пункт тестирования. */
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

//-----------------------------------------------------------------------------
void test2(void) {
	/* 5-й пункт тестирования. */
	for (int k = 1; k < std::numeric_limits<real>::digits10+2; ++k) {
		Matrix a_d;
		generateGilbertMatrix(k, a_d);

		Vector x;
		x.generate(k);

		MatrixProfileSymmetric a(a_d);

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
		sout << "_" << k << ".txt";

		Vector sub;
		f.saveToFile("test2/result" + sout.str());
		f.negate();
		sum(x, f, sub);
		sub.saveToFile("test2/sub" + sout.str());
	}
}

//-----------------------------------------------------------------------------
void test3(void) {
	/* 7-й пункт тестирования. */
	for (int k = 1; k < std::numeric_limits<real>::digits10+2; ++k) {
		Matrix a;
		generateGilbertMatrix(k, a);

		Matrix x;
		generateVector(k, x);

		Matrix y;
		mul(a, x, y);

		Matrix x2;
		solveSLAE_byGaussMethod(a, y, x2);

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
		sout << "_" << k << ".txt";

		Vector sub, x2_s(x2), x_s(x);
		x2_s.saveToFile("test3/result" + sout.str());
		x2_s.negate();
		sum(x_s, x2_s, sub);
		sub.saveToFile("test3/sub" + sout.str());
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main() {
	test1();
	test2();
	test3();
}