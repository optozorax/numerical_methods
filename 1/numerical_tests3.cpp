/* 7-й пункт тестирования. */
#include <sstream>
#include <cstdlib>
#include "matrix.h"
#include "sparse.h"

int main() {
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

		Matrix sub;
		x2.saveToFile("test3/result" + sout.str());
		x2.negate();
		sum(x, x2, sub);
		sub.saveToFile("test3/sub" + sout.str());
	}
}