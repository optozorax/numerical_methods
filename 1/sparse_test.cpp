#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "sparse.h"

//-----------------------------------------------------------------------------
bool testSolve(int size, int maxDistanceToDiagonal) {
	MatrixProfileSymmetric a;
	a.generate(size, 1, 20, maxDistanceToDiagonal);

	Vector x;
	x.generate(5, 1, 20);

	Vector y;
	mul(a, x, y);

	solveSLAE_by_LDL(a, y);
	y.negate();

	Vector sub;
	sum(x, y, sub);

	double sum = 0;
	for (int i = 0; i < sub.size(); ++i)
		sum += sub(i);

	return fabs(sum) < 0.0001;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
TEST_CASE("Random matrixes") {
	CHECK(testSolve(5, 2));
	CHECK(testSolve(50, 2));
	CHECK(testSolve(500, 3));
	CHECK(testSolve(5000, 50));
	CHECK(testSolve(50000, 50));
	//CHECK(testSolve(500000, 50));
	//CHECK(testSolve(5000000, 5));

	//system("pause");
}