#define CATCH_CONFIG_RUNNER

#include "catch.hpp"
#include "matrix.h"
#include "sparse.h"

//-----------------------------------------------------------------------------
void testFormatByDenseMatrix(const Matrix& a_d) {
	// Общий тест перевода матрицы из плотного формата в разреженный

	// Если плотная матрица не симметрична, то ее нельзя тестировать
	if (!a_d.isSymmetric())
		throw std::exception();

	// Переводим матрицу из плотного формата
	MatrixProfileSymmetric a_s(a_d);
	
	// Проверка размера
	CHECK(a_s.size() == a_d.height());

	// Проверка обратного преобразованя
	Matrix sub, a_d2;
	a_s.toDenseMatrix(a_d2);
	a_d2.negate();
	sum(a_d, a_d2, sub);
	CHECK(isNear(sumAllElementsAbs(sub), 0));

	//-------------------------------------------------------------------------
	// Проверка формата
	for (int i = 0; i < a_d.height(); ++i) {
		// Проверка диагональных элементов
		CHECK(isNear(a_s.getDiagonalElement(i), a_d(i, i)));

		// Определяем, с какого элемента начинается строка
		int j = 0;
		for (; j < i; ++j)
			if (a_d(i, j) != 0)
				break;

		// Проверка старта строки
		CHECK(a_s.getLineFirstElementPos(i) == j);

		// Проверка размера строки
		CHECK(a_s.getLineSize(i) == (i - j));

		// Проверка значений строки
		for (int k = j; k < i; ++k)
			CHECK(a_s.getLineElement(i, k - j) == a_d(i, k));
	}
}

//-----------------------------------------------------------------------------
void testVectorByDenseMatrix(const Matrix& x_d) {
	// Общий тест перевода вектора из плотного формата в разреженный

	// Если это не вектор, то проверить нельзя
	if (x_d.width() != 1)
		throw std::exception();

	// Переводим матрицу из плотного формата
	Vector x_s(x_d);
	
	// Проверка размера
	CHECK(x_s.size() == x_d.height());

	// Проверка обратного преобразованя
	Matrix sub, x_d2;
	x_s.toDenseMatrix(x_d2, true);
	x_d2.negate();
	sum(x_d, x_d2, sub);
	CHECK(isNear(sumAllElementsAbs(sub), 0));

	//-------------------------------------------------------------------------
	// Проверка формата
	for (int i = 0; i < x_d.height(); ++i)
		CHECK(isNear(x_s(i), x_d(i, 0)));
}

//-----------------------------------------------------------------------------
void testMultiplication(const Matrix& a_d, const Matrix& b_d) {
	// Общий тест перемножения матриц в профильном формате через перемножение матриц в плотном формате

	// Если матрица не является симметричной и ее нельзя перемножить на второй вектор, то тестировать это нельзя
	if (!a_d.isSymmetric() || b_d.width() != 1 || b_d.height() != a_d.height())
		throw std::exception();

	MatrixProfileSymmetric a_s(a_d);
	Vector b_s(b_d);

	// Перемножаем исходные матрицы в плотном формате
	Matrix c_d;
	mul(a_d, b_d, c_d);

	// Перемножаем исходные матрицы в профильном формате
	Vector c_s;
	mul(a_s, b_s, c_s);

	// Преобразуем полученную матрицу в плотный формат
	Matrix c_d2;
	c_s.toDenseMatrix(c_d2, true);

	// Находим разность между ними
	Matrix sub;
	c_d2.negate();
	sum(c_d, c_d2, sub);

	// Она должна быть близка нулю
	CHECK(isNear(sumAllElementsAbs(sub), 0));
}

//-----------------------------------------------------------------------------
void testLDL(const Matrix& l_d, const Matrix& d_d) {
	// Общий тест LDL^T разложения путем перемножения матриц L и D, и сравнения их с результатом разложения

	// Если матрицы L и D не являются теми, кем должны являться, то проверять нельзя
	if (!l_d.isLowerTriangular() || !l_d.isDiagonalIdentity() || !d_d.isDiagonal())
		throw std::exception();

	// Получаем транспонированную матрицу в плотном формате
	Matrix lt_d = l_d;
	transpose(lt_d);

	// a = l * d * l^T
	Matrix a_d, ld_d;
	mul(l_d, d_d, ld_d);
	mul(ld_d, lt_d, a_d);

	// Считаем LDL^T разложение в плотном формате
	MatrixProfileSymmetric a_s(a_d);
	calcLDL(a_s);

	// Проверяем полученные данные с имеющимися матрицами
	for (int i = 0; i < a_s.size(); ++i) {
		// Проверяем диагональные элементы
		CHECK(isNear(a_s.getDiagonalElement(i), d_d(i, i)));

		// Проверяем L
		for (int j = 0; j < a_s.getLineSize(i); ++j)
			CHECK(isNear(
				a_s.getLineElement(i, j), 
				l_d(i, j + a_s.getLineFirstElementPos(i))
			));
	}
}

//-----------------------------------------------------------------------------
void testSolve(const Matrix& a_d, const Matrix& x_d) {
	// Проверка решения полного цикла СЛАУ через умножение матриц
	if (!a_d.isSymmetric() || x_d.width() != 1 || x_d.height() != a_d.height())
		throw std::exception();

	// Генерируем матрицу и вектор заданного размера
	MatrixProfileSymmetric a_s(a_d);
	Vector x_s(x_d);

	// Получаем вектор: y = a * x
	Vector y_x;
	mul(a_s, x_s, y_x);

	// Решаем СЛАУ
	solveSLAE_by_LDL(a_s, y_x);
	y_x.negate(); // y_x = x*

	// Отнимаем полученное решение от изначально сгенерированного вектора
	Matrix sub, x2_d;
	y_x.toDenseMatrix(x2_d, true);
	sum(x2_d, x_d, sub); // sub = x - x*

	// Разность решения и изначального вектора должна быть близка к нулю
	CHECK(isNear(sumAllElementsAbs(sub), 0));
}

//-----------------------------------------------------------------------------
void testSolveOnBigData(int size, int maxDistanceToDiagonal) {
	// Проверка решения полного цикла СЛАУ через умножение случайно сгенерированных матриц

	// Генерируем матрицу и вектор заданного размера
	MatrixProfileSymmetric a;
	Vector x;
	a.generate(size, 1, 20, maxDistanceToDiagonal);
	x.generate(5, 1, 20);

	// Получаем вектор: y = a * x
	Vector y_x;
	mul(a, x, y_x);

	// Решаем СЛАУ
	solveSLAE_by_LDL(a, y_x);
	y_x.negate(); // y_x = x*

	// Отнимаем полученное решение от изначально сгенерированного вектора
	Matrix sub, x2_d, x_d;
	y_x.toDenseMatrix(x2_d, true);
	x.toDenseMatrix(x_d, true);
	sum(x2_d, x_d, sub); // sub = x - x*

	// Разность решения и изначального вектора должна быть близка к нулю
	CHECK(isNear(sumAllElementsAbs(sub), 0));
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
TEST_CASE("Dense converting test") {
	Matrix a;

	a.getFromVector(5, 5, {
		1,	0,	2,	0,	0,
		0,	3,	0,	0,	0,
		2,	0,	-1,	2,	4,
		0,	0,	2,	2.2,0,
		0,	0,	4,	0,	5,
	});
	testFormatByDenseMatrix(a);

	a.getFromVector(6, 6, {
		4,	1,	1,	0,	0,	0,
		1,	9,	6,	0,	0,	0,
		1,	6,	16,	0,	4,	0,
		0,	0,	0,	99,	-1,	0,
		0,	0,	4,	-1,	100,100,
		0,	0,	0,	0,	100,-99,
	});
	testFormatByDenseMatrix(a);

	generateSparseSymmetricMatrix(10, 0, 100, 0.5, a);
	testFormatByDenseMatrix(a);

	generateSparseSymmetricMatrix(100, 0, 100, 0.5, a);
	testFormatByDenseMatrix(a);

	generateSparseSymmetricMatrix(1000, 0, 100, 0.1, a);
	testFormatByDenseMatrix(a);

	//-------------------------------------------------------------------------
	Matrix x;

	x.getFromVector(1, 5, {
		0, 2, 3, 4, 5
	});
	testVectorByDenseMatrix(x);

	x.getFromVector(1, 6, {
		0, 2, 3, 4, 5, 0
	});
	testVectorByDenseMatrix(x);

	generateVector(10, 0, 100, x);
	testVectorByDenseMatrix(x);

	generateVector(100, 0, 100, x);
	testVectorByDenseMatrix(x);

	generateVector(100, 0, 100, x);
	testVectorByDenseMatrix(x);
}

//-----------------------------------------------------------------------------
TEST_CASE("Multiplication test") {
	Matrix a, x;

	a.getFromVector(5, 5, {
		1,	0,	2,	0,	0,
		0,	3,	0,	0,	0,
		2,	0,	-1,	2,	4,
		0,	0,	2,	2.2,0,
		0,	0,	4,	0,	5,
	});
	x.getFromVector(1, 5, {
		1, 2, 3, 4, 5
	});
	testMultiplication(a, x);

	a.getFromVector(6, 6, {
		4,	1,	1,	0,	0,	0,
		1,	9,	6,	0,	0,	0,
		1,	6,	16,	0,	4,	0,
		0,	0,	0,	99,	-1,	0,
		0,	0,	4,	-1,	100,100,
		0,	0,	0,	0,	100,-99,
	});
	x.getFromVector(1, 6, {
		1, 2, 3, 4, 5, 6
	});
	testMultiplication(a, x);

	generateSparseSymmetricMatrix(10, 0, 100, 0.5, a);
	generateVector(10, 0, 100, x);
	testMultiplication(a, x);

	generateSparseSymmetricMatrix(100, 0, 100, 0.5, a);
	generateVector(100, 0, 100, x);
	testMultiplication(a, x);

	generateSparseSymmetricMatrix(1000, 0, 100, 0.1, a);
	generateVector(1000, 0, 100, x);
	testMultiplication(a, x);
}

//-----------------------------------------------------------------------------
TEST_CASE("LDL test") {
	Matrix l, d;

	l.getFromVector(5, 5, {
		1,	0,	0,	0,	0,
		0,	1,	0,	0,	0,
		0,	2,	1,	0,	0,
		0,	3,	-1,	1,	0,
		0,	0,	0,	5,	1,
	});
	d.getFromVector(5, 5, {
		1,	0,	0,	0,	0,
		0,	6,	0,	0,	0,
		0,	0,	10,	0,	0,
		0,	0,	0,	1,	0,
		0,	0,	0,	0,	-1,
	});
	testLDL(l, d);

	l.getFromVector(6, 6, {
		1,	0,	0,	0,	0,	0,
		2,	1,	0,	0,	0,	0,
		1,	0,	1,	0,	0,	0,
		0,	0,	0,	1,	0,	0,
		0,	0,	-9,	5,	1,	0,
		1,	2,	3,	4,	5,	1,
	});
	d.getFromVector(6, 6, {
		10,	0,	0,	0,	0,	0,
		0,	9,	0,	0,	0,	0,
		0,	0,	8,	0,	0,	0,
		0,	0,	0,	7,	0,	0,
		0,	0,	0,	0,	6,	0,
		0,	0,	0,	0,	0,	5,
	});
	testLDL(l, d);

	generateLMatrix(10, 0, 100, 0.5, l);
	generateDiagonalMatrix(10, 1, 100, d);
	testLDL(l, d);

	generateLMatrix(100, 0, 100, 0.5, l);
	generateDiagonalMatrix(100, 1, 100, d);
	testLDL(l, d);

	generateLMatrix(1000, 0, 100, 0.1, l);
	generateDiagonalMatrix(1000, 1, 100, d);
	testLDL(l, d);
}

//-----------------------------------------------------------------------------
TEST_CASE("Random matrixes") {
	Matrix a, x;

	a.getFromVector(5, 5, {
		1,	0,	2,	0,	0,
		0,	3,	0,	0,	0,
		2,	0,	-1,	2,	4,
		0,	0,	2,	2.2,0,
		0,	0,	4,	0,	5,
	});
	x.getFromVector(1, 5, {
		1, 2, 3, 4, 5
	});
	testMultiplication(a, x);

	a.getFromVector(6, 6, {
		4,	1,	1,	0,	0,	0,
		1,	9,	6,	0,	0,	0,
		1,	6,	16,	0,	4,	0,
		0,	0,	0,	99,	-1,	0,
		0,	0,	4,	-1,	100,100,
		0,	0,	0,	0,	100,-99,
	});
	x.getFromVector(1, 6, {
		1, 2, 3, 4, 5, 6
	});
	testMultiplication(a, x);

	testSolveOnBigData(5000, 50);
	testSolveOnBigData(50000, 50);
	testSolveOnBigData(500000, 50);
	testSolveOnBigData(5000000, 5);
	testSolveOnBigData(10000000, 5);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main(int argc, char* const argv[]) {
	int result = Catch::Session().run(argc, argv);

	system("pause");
	return result;
}