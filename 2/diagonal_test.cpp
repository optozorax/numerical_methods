#define CATCH_CONFIG_RUNNER

#include "../1/catch.hpp"
#include "../1/matrix.h"
#include "../1/vector.h"
#include "diagonal.h"

/*
	Тестирование класса диагональной матрицы
		Создание по формату
		Создание из плотной матрицы
			Проверяются все методы
		Преобразование из плотной туда-сюда
	Тестирование итератора по строке
		Проходятся все элементы, проверяется, действительно ли они находятся на той позиции, что должны, все пройденные элементы отмечаются -1
		Если еще раз встречается на -1, то ошибка
		Если после прохода все элементы равны только -1 то все хорошо

		Это делается для различных методов
	Тестирование итерационных методов
		Тестирование с помощью итерационных методов в плотной матрице
			Для определенных матриц
			Для случайных диагональных матриц
			Для абсолютно плотных матриц
		Тестирование невязки
			Для определенных матриц
			Для случайных диагональных матриц
			Для абсолютно плотных матриц
*/

//-----------------------------------------------------------------------------
void generateTestMatrix(int n, std::vector<int> format, bool isNegative,  MatrixDiagonal& result) {
	result.resize(n, format);

	for (int i = 0; i < result.getDiagonalsCount(); ++i) {
		auto mit = result.posBegin(i);
		for (auto it = result.begin(i); it != result.end(i); ++it, ++mit) {
			if (isNegative)
				*it = -intRandom(0, 5);
			else
				*it = intRandom(0, 5);
		}
	}

	matrix_diagonal_line_iterator mit(n, format, false);
	for (; !mit.isEnd(); ++mit) {
		sumreal& sum = result.begin(0)[mit.i];
		sum = 0;
		for (; !mit.isLineEnd(); ++mit)
			if (mit.i != mit.j)
				sum += result.begin(mit.dn)[mit.di];
		sum = std::fabs(sum);
	}

	result.begin(0)[0] += 1;
}

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
			//CHECK(fabs(relativeResidual - epsilon) / relativeResidual < 0.05);
			CHECK(relativeResidual < epsilon);
		}
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

TEST_CASE("Test residual of methods") {
	for (int i = 10; i < 100; ++i) {
		for (int j = 0; j < 1; ++j) {
			MatrixDiagonal a;
			auto format = generateRandomFormat(i, intRandom(10, Diagonal(i).calcDiagonalsCount()));
			generateTestMatrix(i, format, true, a);

			Vector x;
			x.generate(i);

			testResidual(a, x, &SolverSLAE_Iterative::jacobi, 1e-6);
			testResidual(a, x, &SolverSLAE_Iterative::seidel, 1e-6);
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