/* 5-й пункт тестирования. */
#include <sstream>
#include <cstdlib>
#include "matrix.h"
#include "sparse.h"

int main() {
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