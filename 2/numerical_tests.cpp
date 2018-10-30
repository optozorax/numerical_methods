#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "diagonal.h"

typedef std::function<IterationsResult(const SolverSLAE_Iterative*, const MatrixDiagonal& a, const Vector& y, Vector& x)> method_function;

//-----------------------------------------------------------------------------
void makeTable(
	const MatrixDiagonal& a, 
	const Vector& x_precise, 
	const Vector& y, 
	SolverSLAE_Iterative& solver, 
	std::string fileName
) {
	std::ofstream fout(fileName + ".tex");
	std::ofstream fout1(fileName + ".dat");

	// Выводим матрицу и остальное в виде формулы
	auto format = a.getFormat();
	fout << "$$ " << fileName << "=\\left(\\quad\\begin{matrix}\n";
	Matrix a_dense;
	a.toDenseMatrix(a_dense);
	for (int i = 0; i < a_dense.height(); i++) {
		for (int j = 0; j < a_dense.width(); j++) {
			int d = Diagonal(a_dense.height()).calcDiag_byLR(i, j);
			bool isOnFormat = std::find(format.begin(), format.end(), d) != format.end();
			if (isOnFormat)
				fout << "\\cellcolor{green!30}";
			fout << int(a_dense(i, j));
			if (j + 1 != a_dense.width())
				 fout << " & ";
		}
		if (i + 1 != a_dense.height())
			fout << " \\\\\n";
		else
			fout << " \n";
	}
	fout << "\\end{matrix}\\quad\\right), X=\\begin{pmatrix}";
	for (int i = 0; i < x_precise.size(); i++) {
		if (i + 1 != x_precise.size())
			fout << int(x_precise(i)) << " \\\\\n";
		else
			fout << int(x_precise(i)) << " \n";
	}
	fout << "\\end{pmatrix}, F=\\begin{pmatrix}";
	for (int i = 0; i < y.size(); i++) {
		if (i + 1 != y.size())
			fout << int(y(i)) << " \\\\\n";
		else
			fout << int(y(i)) << " \n";
	}
	fout << "\\end{pmatrix} $$\n\n";

	// Выводим параметры решателя
	int exponent = floor(log10(solver.epsilon));
	double number = solver.epsilon / pow(10.0, exponent);
	fout
		<< "$$ \\varepsilon = ";
	if (fabs(number - 1) >= 0.01)
		fout
			<< std::setprecision(2) << std::fixed << number
			<< " \\cdot ";
	fout 
		<< "10^{" << exponent << "}, \\quad iterations_{max} = "
		<< solver.maxIterations << ", \\quad start = \\begin{pmatrix} ";
	fout << std::defaultfloat;
	for (int i = 0; i < solver.start.size(); i++) {
		fout << solver.start(i);
		if (i + 1 != solver.start.size())
			fout << " & ";
		else
			fout << " ";
	}
	fout << "\\end{pmatrix}^T $$\n\n";

	std::vector<double> w1(200), w2(200);
	std::vector<Vector> x1(200), x2(200);
	std::vector<Vector> xsub1(200), xsub2(200);
	std::vector<double> rr1(200), rr2(200); // relativeResidual
	std::vector<double> va1(200), va2(200);
	std::vector<int> it1(200), it2(200);

	int min1, min2;
	int count1, count2;

	auto one_method = [&a, &x_precise, &y, &solver] (
		std::vector<double>& w,
		std::vector<Vector>& x, 
		std::vector<Vector>& xsub,
		std::vector<double>& rr, 
		std::vector<double>& va, 
		std::vector<int>& it,

		int& min, 
		int& count,

		method_function method
	) {
		min = 0;
		count = 200;
		Vector x_solve(x_precise.size());
		Vector x_sub(x_precise.size());
		real xNorm = calcNorm(x_precise);

		for (int i = 0; i < 200; ++i) {
			solver.w = i / 100.0;
			auto result = method(&solver, a, y, x_solve);

			// Если начинается ошибки после 100 итерации, то и потом ничего кроме них не будет, поэтому заканчиваем цикл
			if ((result.relativeResidual > solver.epsilon && i >= 100) || 
				(result.relativeResidual != result.relativeResidual)) {
				count = i;
				break;
			}

			// Вычисления разности точного и приближенного решений
			x_sub = x_solve;
			x_sub.negate();
			sum(x_sub, x_precise, x_sub);
			real x_subNorm = calcNorm(x_sub);

			w[i] = solver.w;
			x[i] = x_solve;
			xsub[i] = x_sub;
			rr[i] = result.relativeResidual;
			va[i] = x_subNorm / xNorm / result.relativeResidual;
			it[i] = result.iterations;

			// Находим минимум
			if (result.iterations < it[min])
				min = i;
		}
	};

	one_method(w1, x1, xsub1, rr1, va1, it1, min1, count1, &SolverSLAE_Iterative::jacobi);
	one_method(w2, x2, xsub2, rr2, va2, it2, min2, count2, &SolverSLAE_Iterative::seidel);

	// w, x, x-x*, относительная невязка, vA, итераций
	fout
		<< "\\setlength{\\tabcolsep}{2pt}\n"
		<< "\\tabulinesep=0.3mm\n"
		<< "\\noindent{\\scriptsize\\texttt{"
		<< "\\begin{longtabu}{\n"
		<< "|X[-1,c]||X[-1,c]|X[-1,c]|X[-1,c]|X[-1,c]|X[-1,c]|\n"
		<< "p{0.05cm}\n|X[-1,c]||X[-1,c]|X[-1,c]|X[-1,c]|X[-1,c]|X[-1,c]|}\n"
		<< "\\cline{1-6}\\cline{8-13}\n"
		<< "\\multicolumn{6}{|c|}{Метод Якоби} && \\multicolumn{6}{c|}{Метод Зейделя} \\\\\n"
		<< "\\cline{1-6}\\cline{8-13}\n"
		<< "$w$ & $x$ & $x-x^*$ & \\tcell{\\tiny Относительная\\\\\\tiny невязка} & {\\tiny $\\mathop{cond}(A) >$} & {\\tiny Итераций} && $w$ & $x$ & $x-x^*$ & \\tcell{\\tiny Относительная\\\\\\tiny невязка} & {\\tiny $\\mathop{cond}(A) >$} & {\\tiny Итераций} \\\\\n"
		<< "\\cline{1-6}\\cline{8-13}\n";

	auto write_vector = [&fout] (const Vector& a) {
		fout << "\\tcell{";
		for (int i = 0; i < a.size(); ++i)
			if (i + 1 != a.size())
				fout << a(i) << " \\\\ ";
			else
				fout << a(i) << "}";
	};

	auto write_line = [&] (int i, int colorNo) {
		int doublePrec = 16;
		if (i < count1) {
			if (colorNo == 1) {
				fout << std::fixed << std::setprecision(2) << "\\cellcolor{green!30}{" << w1[i] << "} & ";
				fout << std::fixed << std::setprecision(doublePrec) << "\\tiny{\\cellcolor{green!30}{";
				write_vector(x1[i]);
				fout << "}} & " << "\\tiny{\\cellcolor{green!30}{";
				fout << std::scientific << std::setprecision(1);
				write_vector(xsub1[i]);
				fout << "}} & ";
				fout << std::scientific << std::setprecision(1) << "\\cellcolor{green!30}{" << rr1[i] << "} & ";
				fout << std::fixed << std::setprecision(2) << "\\cellcolor{green!30}{" << va1[i] << "} & ";
				fout << "\\cellcolor{green!30}{" << it1[i] << "} & ";
			} else {
				fout << std::fixed << std::setprecision(2) << w1[i] << " & ";
				fout << std::fixed << std::setprecision(doublePrec) << "\\tiny{";
				write_vector(x1[i]);
				fout << "} & " << "\\tiny{";
				fout << std::scientific << std::setprecision(2);
				write_vector(xsub1[i]);
				fout << "} & ";
				fout << std::scientific << std::setprecision(2) << rr1[i] << " & ";
				fout << std::fixed << std::setprecision(2) << va1[i] << " & ";
				fout << it1[i] << " & ";
			}
		} else {
			fout << "& & & & & &";
		}

		fout << " & ";

		if (i < count2) {
			if (colorNo == 2) {
				fout << std::fixed << std::setprecision(2) << "\\cellcolor{green!30}{" << w2[i] << "} & ";
				fout << std::fixed << std::setprecision(doublePrec) << "\\tiny{\\cellcolor{green!30}{";
				write_vector(x2[i]);
				fout << "}} & " << "\\tiny{\\cellcolor{green!30}{";
				fout << std::scientific << std::setprecision(1);
				write_vector(xsub2[i]);
				fout << "}} & ";
				fout << std::scientific << std::setprecision(2) << "\\cellcolor{green!30}{" << rr2[i] << "} & ";
				fout << std::fixed << std::setprecision(2) << "\\cellcolor{green!30}{" << va2[i] << "} & ";
				fout << "\\cellcolor{green!30}{" << it2[i] << "} \\\\";
			} else {
				fout << std::fixed << std::setprecision(2) << w2[i] << " & ";
				fout << std::fixed << std::setprecision(doublePrec) << "\\tiny{";
				write_vector(x2[i]);
				fout << "} & " << "\\tiny{";
				fout << std::scientific << std::setprecision(1);
				write_vector(xsub2[i]);
				fout << "} & ";
				fout << std::scientific << std::setprecision(2) << rr2[i] << " & ";
				fout << std::fixed << std::setprecision(2) << va2[i] << " & ";
				fout << it2[i] << " \\\\";
			}
		} else {
			fout << "& & & & & \\\\";	
		}

		fout << "\n";
		fout << "\\cline{1-6}\\cline{8-13}\n";
	};

	for (int i = 0; i < std::max(count1, count2); i+=10) {
		if (min1 == i) {
			if (i != 0) write_line(i-1, 0);
			write_line(i, 1);
			write_line(i+1, 0);
		} else {
			if (min2 == i) {
				if (i != 0) write_line(i-1, 0);
				write_line(i, 2);
				write_line(i+1, 0);
			} else {
				write_line(i, 0);
			}
		}
		if (i + 10 > min1 && min1 > i) {	
			if (min1-1 != i) write_line(min1-1, 0);
			write_line(min1, 1);
			if (min1+1 != i+10) write_line(min1+1, 0);
		} else {
			if (i + 10 > min2 && min2 > i) {
				if (min2-1 != i) write_line(min2-1, 0);
				write_line(min2, 2);
				if (min2+1 != i+10) write_line(min2+1, 0);
			}
		}
	}

	fout1 << "w1\tit1\tw2\tit2" << std::endl;
	fout1 << std::fixed << std::setprecision(2);
	for (int i = 0; i < std::max(count1, count2); ++i) {
		if (i >= count1)
			fout1 << w1[count1-1] << "\t" << it1[count1-1] << "\t";
		else
			fout1 << w1[i] << "\t" << it1[i] << "\t";
		if (i >= count2)
			fout1 << w2[count2-1] << "\t" << it2[count2-1] << std::endl;
		else
			fout1 << w2[i] << "\t" << it2[i] << std::endl;
	}

	fout << "\\end{longtabu}}}\n\n";

	fout 
		<< "\\noindent\\begin{tikzpicture}\n"
		<< "\\begin{semilogyaxis}[xlabel=w,ylabel=Iterations,width=\\textwidth, height=6cm]\n"
		<< "\\addplot[red, no markers] table [y=it1, x=w1]{" << fileName << ".dat};\n"
		<< "\\addplot[blue, no markers] table [y=it2, x=w2]{" << fileName << ".dat};\n"
		<< "\\legend{Jacobi,Seidel}\n"
		<< "\\end{semilogyaxis}\n"
		<< "\\end{tikzpicture}";

	fout.close();
	fout1.close();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main() {
	// Получаем необходимые матрицы
	MatrixDiagonal a, b;
	generateDiagonallyDominantMatrix(10, makeSevenDiagonalFormat(10, 3, 2), true, a);
	generateDiagonallyDominantMatrix(10, makeSevenDiagonalFormat(10, 2, 4), false, b);

	Vector x;
	x.generate(10);

	Vector y_a, y_b;
	mul(a, x, y_a);
	mul(b, x, y_b);

	// Начальные присвоения
	SolverSLAE_Iterative solver;
	solver.w = 0;
	solver.isLog = false;
	solver.start = Vector(10, 0);
	solver.epsilon = 1e-14;
	solver.maxIterations = 1e5;

	// Создаем таблицы
	makeTable(a, x, y_a, solver, "A");
	makeTable(b, x, y_b, solver, "B");
}