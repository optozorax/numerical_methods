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

	//-------------------------------------------------------------------------
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

	//-------------------------------------------------------------------------
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

	//-------------------------------------------------------------------------
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

	// Заполняем массивы расчитанными данными
	one_method(w1, x1, xsub1, rr1, va1, it1, min1, count1, &SolverSLAE_Iterative::jacobi);
	one_method(w2, x2, xsub2, rr2, va2, it2, min2, count2, &SolverSLAE_Iterative::seidel);

	//-------------------------------------------------------------------------
	// Выводим преамбулу таблицы
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
		<< "\\hhline{*{6}{-}~*{6}{-}}\n"
		<< "$w$ & $x$ & $x-x^*$ & \\tcell{\\tiny Относительная\\\\\\tiny невязка} & {\\tiny $\\mathop{cond}(A) >$} & {\\tiny Итераций} && $w$ & $x$ & $x-x^*$ & \\tcell{\\tiny Относительная\\\\\\tiny невязка} & {\\tiny $\\mathop{cond}(A) >$} & {\\tiny Итераций} \\\\\n"
		<< "\\hhline{*{6}{-}~*{6}{-}}\n";

	//-------------------------------------------------------------------------
	auto write_vector = [&fout] (const Vector& a) {
		fout << "\\tcell{";
		for (int i = 0; i < a.size(); ++i)
			if (i + 1 != a.size())
				fout << a(i) << " \\\\ ";
			else
				fout << a(i) << "}";
	};

	auto write_line = [&] (int i, bool isWrite1, bool isWrite2, std::string color1, std::string color2) {
		int doublePrec = 16;
		if (isWrite1) {
			std::string before_cell = "\\cellcolor{" + color1 + "}{";
			std::string after_cell = "}";
			if (color1 == "white") {
				before_cell = "";
				after_cell = "";
			}

			fout << std::fixed << std::setprecision(2)
				 << before_cell << w1[i] << after_cell <<" & ";

			fout << std::fixed << std::setprecision(doublePrec) 
				 << "\\tiny{" << before_cell;
				 write_vector(x1[i]);
			fout << after_cell << "} & ";

			fout << std::scientific << std::setprecision(1)
 				 << "\\tiny{" << before_cell;
				 write_vector(xsub1[i]);
			fout << after_cell << "} & ";

			fout << std::scientific << std::setprecision(1) 
				 << before_cell << rr1[i] << after_cell << " & ";

			fout << std::fixed << std::setprecision(2) 
				 << before_cell << va1[i] << after_cell << " & ";

			fout << before_cell << it1[i] << after_cell << " & ";
		} else {
			fout << "& & & & & &";
		}

		fout << " & ";

		if (isWrite2) {
			std::string before_cell = "\\cellcolor{" + color2 + "}{";
			std::string after_cell = "}";
			if (color2 == "white") {
				before_cell = "";
				after_cell = "";
			}

			fout << std::fixed << std::setprecision(2)
				 << before_cell << w2[i] << after_cell <<" & ";

			fout << std::fixed << std::setprecision(doublePrec) 
				 << "\\tiny{" << before_cell;
				 write_vector(x2[i]);
			fout << after_cell << "} & ";

			fout << std::scientific << std::setprecision(1)
 				 << "\\tiny{" << before_cell;
				 write_vector(xsub2[i]);
			fout << after_cell << "} & ";

			fout << std::scientific << std::setprecision(1) 
				 << before_cell << rr2[i] << after_cell << " & ";

			fout << std::fixed << std::setprecision(2) 
				 << before_cell << va2[i] << after_cell << " & ";
				 
			fout << before_cell << it2[i] << after_cell << " \\\\";
		} else {
			fout << "& & & & & \\\\";	
		}

		fout << "\n";
		fout << "\\hhline{*{6}{-}~*{6}{-}}\n";
	};

	//-------------------------------------------------------------------------
	// Формируем массив тех значений, которые надо вывести
	int tableSize = std::max(count1, count2);
	std::vector<int> to_write;
	for (int i = 10; i < tableSize; i+=10)
		to_write.push_back(i);
	to_write.push_back(min1-1);
	to_write.push_back(min1);
	to_write.push_back(min1+1);
	to_write.push_back(min2-1);
	to_write.push_back(min2);
	to_write.push_back(min2+1);

	std::sort(to_write.begin(), to_write.end(), std::less<int>());

	// Удаляем дубликаты
	to_write.erase(std::unique(to_write.begin(), to_write.end()), to_write.end());

	// Удаляем отрицательные значения
	while (to_write.front() <= 0)
		to_write.erase(to_write.begin());

	// Удаляем значения за допустимыми пределами
	while (to_write.back() >= tableSize)
		to_write.pop_back();

	// Выводим каждую строку
	std::string orange = "orange!30";
	std::string green = "green!30";
	for (auto& i : to_write) {
		std::string color1 = "white", color2 = "white";
		if (i == min1)
			color1 = green;
		if (abs(i-min1) == 1)
			color1 = orange;
		if (i == min2)
			color2 = green;
		if (abs(i-min2) == 1)
			color2 = orange;
		write_line(i, i < count1, i < count2, color1, color2);
	}

	//-------------------------------------------------------------------------
	// Выводим данные для построения графика
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

	//-------------------------------------------------------------------------
	// Выводим график
	fout 
		<< "\\subsubsection{График зависимости числа итераций от параметра релаксации}\n\n"
		<< "\\noindent\\begin{tikzpicture}\n"
		<< "\\begin{semilogyaxis}[xlabel=w,ylabel=Итераций,width=\\textwidth, height=6cm]\n"
		<< "\\addplot[red, no markers] table [y=it1, x=w1]{" << fileName << ".dat};\n"
		<< "\\addplot[blue, no markers] table [y=it2, x=w2]{" << fileName << ".dat};\n"
		<< "\\legend{Метод Якоби,Метод Зейделя}\n"
		<< "\\end{semilogyaxis}\n"
		<< "\\end{tikzpicture}";

	fout.close();
	fout1.close();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main() {
	MatrixDiagonal a, b;
	Vector x;
	Vector y_a, y_b;
	SolverSLAE_Iterative solver;

	std::ifstream fin("num.txt");
	a.load(fin);
	b.load(fin);
	x.load(fin);
	solver.load(fin);
	fin.close();
		
	mul(a, x, y_a);
	mul(b, x, y_b);

	// Создаем таблицы
	makeTable(a, x, y_a, solver, "A");
	makeTable(b, x, y_b, solver, "B");

	Matrix a_dense, b_dense;
	a.toDenseMatrix(a_dense);
	b.toDenseMatrix(b_dense);

	std::ofstream fout("matrixes.txt");
	fout << std::defaultfloat;
	a_dense.save(fout);
	b_dense.save(fout);
	fout.close();

	/*// Сохраняем сгенерированные данные в файл
	std::ofstream fout("num.txt");

	generateDiagonallyDominantMatrix(10, makeSevenDiagonalFormat(10, 3, 2), true, a);
	b = a;
	for (int i = 0; i < b.getDiagonalsCount(); i++)
		for (auto j = b.begin(i); j != b.end(i); j++)
			(*j) = std::fabs(*j);
	x.generate(10);
	solver.w = 0;
	solver.isLog = false;
	solver.start = Vector(10, 0);
	solver.epsilon = 1e-14;
	solver.maxIterations = 1e5;

	a.save(fout);
	b.save(fout);
	x.save(fout);
	solver.save(fout);

	fout.close();*/
}