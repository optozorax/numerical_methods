#define CATCH_CONFIG_RUNNER

#include "../1/catch.hpp"
#include "../1/matrix.h"
#include "../1/vector.h"
#include "diagonal.h"

typedef std::function<IterationsResult(const SolverSLAE_Iterative*, const MatrixDiagonal& a, const Vector& y, Vector& x)> method_function;

//-----------------------------------------------------------------------------
void testResidual(const MatrixDiagonal& a, const Vector& x, method_function method, real epsilon) {
	SolverSLAE_Iterative solver;
	solver.w = 1;
	solver.start = Vector(a.dimension(), 5);
	solver.epsilon = epsilon;
	solver.maxIterations = 1e5;

	Vector y;
	mul(a, x, y);

	Vector x1;
	auto result = method(&solver, a, y, x1);

	Vector y1;
	mul(a, x1, y1);

	y1.negate();
	sum(y, y1, y1);
	real relativeResidual = calcNorm(y1) / calcNorm(y);

	if (relativeResidual == relativeResidual) {
		CHECK(fabs(relativeResidual - result.relativeResidual) / relativeResidual < 0.01);

		if (result.iterations < solver.maxIterations) {
			CHECK(relativeResidual <= epsilon);
		}
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

TEST_CASE("Test residual of methods") {
	for (int i = 10; i < 100; ++i) {
		for (int j = 0; j < 3; ++j) {
			MatrixDiagonal a;
			auto format = generateRandomFormat(i, intRandom(10, Diagonal(i).calcDiagonalsCount()));
			generateDiagonallyDominantMatrix(i, format, intRandom(0, 10) % 2, a);

			Vector x;
			x.generate(i);

			testResidual(a, x, &SolverSLAE_Iterative::jacobi, 1e-10);
			testResidual(a, x, &SolverSLAE_Iterative::seidel, 1e-10);
		}
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main(int argc, char* const argv[]) {
	int result = Catch::Session().run(argc, argv);

	system("pause");
	return result;
}