#include <fstream>
#include <cmath>
#include <iomanip>
#include "diagonal.h"

typedef std::function<IterationsResult(const SolverSLAE_Iterative*, const MatrixDiagonal& a, const Vector& y, Vector& x)> method_function;

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

//-----------------------------------------------------------------------------
void makeTable(
	const MatrixDiagonal& a, 
	const Vector& x, 
	const Vector& y, 
	SolverSLAE_Iterative& solver, 
	std::string fileName, 
	method_function method
) {
	std::ofstream fout(fileName + ".txt");
	std::ofstream fout1("graph_" + fileName + ".txt");

	real wMin = 0;
	int itMin = solver.maxIterations;

	Vector x1(x.size());
	real xNorm = calcNorm(x);
	for (int i = 0; i < 200; ++i) {
		solver.w = i / 100.0;
		auto result = method(&solver, a, y, x1);
		x1.negate();
		sum(x1, x, x1);
		real xNorm1 = calcNorm(x1);

		// Если начинается ошибки после 100 итерации, то и потом ничего кроме них не будет, поэтому заканчиваем цикл
		if ((result.relativeResidual > solver.epsilon && i >= 100) || (result.relativeResidual != result.relativeResidual))
			break;

		// Считаем число обусловленности
		real va = xNorm1 / xNorm / result.relativeResidual;

		// Выводим таблицу в файл
		fout
			<< std::fixed << std::setprecision(2) << solver.w << "\t"
			<< std::scientific << std::setprecision(3) << xNorm1 << "\t"
			<< result.relativeResidual << "\t"
			<< va << "\t"
			<< result.iterations << std::endl;

		// Выводим другую таблицу в файл для построения графика зависимости скорости сходимости от веса
		fout1 << result.iterations << std::endl;

		// Находим минимум
		if (result.iterations < itMin) {
			wMin = solver.w;
			itMin = result.iterations;
		}
	}

	// Выводим вектор и сопутствующую таблицу для лучшего веса
	auto write_vector = [&] (real w) {
		Vector x2;
		solver.w = w;
		auto result = method(&solver, a, y, x1);
		x2 = x1;
		x2.negate();
		sum(x2, x, x2);
		real xNorm1 = calcNorm(x1);

		real va = xNorm1 / xNorm / result.relativeResidual;

		fout
			<< std::fixed << std::setprecision(2) << solver.w << "\t"
			<< std::scientific << std::setprecision(3) << xNorm1 << "\t"
			<< std::fixed << std::setprecision(16) << x1(0) << "\t"
			<< std::scientific << std::setprecision(3) << x2(0) << "\t"
			<< result.relativeResidual << "\t"
			<< va << "\t"
			<< result.iterations << std::endl;

		for (int i = 1; i < x2.size(); ++i) {
			fout
				<< "\t"
				<< "\t"
				<< std::fixed << std::setprecision(16) << x1(i) << "\t"
				<< std::scientific << std::setprecision(3) << x2(i) << "\t"
				<< "\t"
				<< "\t"
				<< std::endl;
		}
	};

	fout << std::endl;
	write_vector(wMin);

	fout.close();
	fout1.close();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main() {
	// Получаем необходимые матрицы
	MatrixDiagonal a, b;
	generateTestMatrix(10, makeSevenDiagonalFormat(10, 3, 2), true, a);
	generateTestMatrix(10, makeSevenDiagonalFormat(10, 2, 4), false, b);

	Vector x;
	x.generate(10);

	Vector y_a, y_b;
	mul(a, x, y_a);
	mul(b, x, y_b);

	// Сохраняем эти матрицы
	Matrix a_d, b_d;
	a.toDenseMatrix(a_d);
	b.toDenseMatrix(b_d);
	a_d.saveToFile("a.txt");
	b_d.saveToFile("b.txt");
	x.saveToFile("x.txt");
	y_a.saveToFile("y_a.txt");
	y_b.saveToFile("y_b.txt");

	// Начальные присвоения
	SolverSLAE_Iterative solver;
	solver.w = 0;
	solver.isLog = false;
	solver.start = Vector(10, 0);
	solver.epsilon = 1e-14;
	solver.maxIterations = 1e5;

	// Цикл с созданием таблицы
	makeTable(a, x, y_a, solver, "a_jacobi", &SolverSLAE_Iterative::jacobi);
	makeTable(a, x, y_a, solver, "a_seidel", &SolverSLAE_Iterative::seidel);
	makeTable(b, x, y_b, solver, "b_jacobi", &SolverSLAE_Iterative::jacobi);
	makeTable(b, x, y_b, solver, "b_seidel", &SolverSLAE_Iterative::seidel);

	//-------------------------------------------------------------------------
	Matrix c;
	c.getFromVector(2, 2, {
		1, 10,
		100, 1001
	});

	Vector x_c(2, 1);

	Vector y_c;
	mul(c, x_c, y_c);

	solver.start = Vector(2, 0);
	solver.epsilon = 1e-8;

	makeTable(c, x_c, y_c, solver, "c_jacobi", &SolverSLAE_Iterative::jacobi);
	makeTable(c, x_c, y_c, solver, "c_seidel", &SolverSLAE_Iterative::seidel);
}