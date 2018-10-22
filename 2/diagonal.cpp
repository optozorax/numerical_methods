#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "diagonal.h"

//-----------------------------------------------------------------------------
MatrixDiagonal::MatrixDiagonal() : n(0) {
}

//-----------------------------------------------------------------------------
MatrixDiagonal::MatrixDiagonal(int n, std::vector<int> format) : n(n), fi(format) {
	if (format[0] != 0)
		throw std::exception();

	for (const auto& i : format)
		di.push_back(std::vector<real>(calcDiagonalSize(n, i), 0));
}

//-----------------------------------------------------------------------------
MatrixDiagonal::MatrixDiagonal(const Matrix& a) {
	if (a.width() != a.height())
		throw std::exception();

	n = a.width();

	// Определяем формат
	fi.clear();
	fi.push_back(0);
	for (int i = -n+1; i < n; ++i) if (i != 0) {
		auto mit = matrix_diagonal_iterator(n, i, false);
		auto mite = matrix_diagonal_iterator(n, i, true);
		for (; mit != mite; ++mit) {
			if (a(mit.i, mit.j) != 0) {
				fi.push_back(i);
				break;
			}
		}
	}

	// Создаем формат
	di.clear();
	for (const auto& i : fi)
		di.push_back(std::vector<real>(calcDiagonalSize(n, i), 0));

	// Обходим массив и записываем элементы
	for (int i = 0; i < getDiagonalsCount(); ++i) {
		auto mit = posBegin(i);
		for (auto it = begin(i); it != end(i); ++it, ++mit)
			*it = a(mit.i, mit.j);
	}
}

//-----------------------------------------------------------------------------
void MatrixDiagonal::toDenseMatrix(Matrix& dense) const {
	dense.resize(n, n, 0);

	// Обходим массив и записываем элементы
	for (int i = 0; i < getDiagonalsCount(); ++i) {
		auto mit = posBegin(i);
		for (auto it = begin(i); it != end(i); ++it, ++mit)
			dense(mit.i, mit.j) = *it;
	}	
}

//-----------------------------------------------------------------------------
int MatrixDiagonal::dimension(void) const {
	return n;
}

//-----------------------------------------------------------------------------
int MatrixDiagonal::getDiagonalsCount(void) const {
	return di.size();
}

//-----------------------------------------------------------------------------
int MatrixDiagonal::getDiagonalSize(int diagNo) const {
	return di[diagNo].size();
}

//-----------------------------------------------------------------------------
int MatrixDiagonal::getDiagonalPos(int diagNo) const {
	return fi[diagNo];
}

//-----------------------------------------------------------------------------
std::vector<int> MatrixDiagonal::getFormat(void) const {
	return fi;
}

//-----------------------------------------------------------------------------
matrix_diagonal_iterator MatrixDiagonal::posBegin(int diagNo) const {
	return matrix_diagonal_iterator(n, fi[diagNo], false);
} 

//-----------------------------------------------------------------------------
matrix_diagonal_iterator MatrixDiagonal::posEnd(int diagNo) const {
	return matrix_diagonal_iterator(n, fi[diagNo], true);
}

//-----------------------------------------------------------------------------
MatrixDiagonal::iterator MatrixDiagonal::begin(int diagNo) {
	return di[diagNo].begin();
}

//-----------------------------------------------------------------------------
MatrixDiagonal::const_iterator MatrixDiagonal::begin(int diagNo) const {
	return di[diagNo].begin();
}

//-----------------------------------------------------------------------------
MatrixDiagonal::iterator MatrixDiagonal::end(int diagNo) {
	return di[diagNo].end();
}

//-----------------------------------------------------------------------------
MatrixDiagonal::const_iterator MatrixDiagonal::end(int diagNo) const {
	return di[diagNo].end();
}

//-----------------------------------------------------------------------------
int MatrixDiagonal::calcDiagonalSize(int n, int i) {
	return n - std::abs(i);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
matrix_diagonal_iterator::matrix_diagonal_iterator(int n, int i1, bool isEnd) {
	int count = MatrixDiagonal::calcDiagonalSize(n, i1);
	if (!isEnd) {
		if (i1 < 0) {
			i = -i1;
			j = 0;
		} else {
			i = 0;
			j = i1;
		}
	} else {
		if (i1 < 0) {
			i = -i1 + count;
			j = 0 + count;
		} else {
			i = 0 + count;
			j = i1 + count;
		}
	}
}

//-----------------------------------------------------------------------------
matrix_diagonal_iterator& matrix_diagonal_iterator::operator++() {
	i++;
	j++;
	return *this;
}

//-----------------------------------------------------------------------------
matrix_diagonal_iterator matrix_diagonal_iterator::operator++(int) {
	i++;
	j++;
	return *this;
}

//-----------------------------------------------------------------------------
bool matrix_diagonal_iterator::operator==(const matrix_diagonal_iterator& b) const {
	return b.i == i && b.j == j;
}

//-----------------------------------------------------------------------------
bool matrix_diagonal_iterator::operator!=(const matrix_diagonal_iterator& b) const {
	return b.i != i || b.j != j;
}

//-----------------------------------------------------------------------------
matrix_diagonal_iterator& matrix_diagonal_iterator::operator+=(const ptrdiff_t& movement) {
	i += movement;
	j += movement;
	return *this;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
matrix_diagonal_line_iterator::matrix_diagonal_line_iterator(int n, std::vector<int> format, bool isOnlyLowTriangle) : n(n), m_isEnd(false), m_isLineEnd(false) {
	// Создаем обратное преобразование из формата диагонали в ее номер в формате
	for (int i = 0; i < format.size(); ++i)
		if ((isOnlyLowTriangle && format[i] < 0) || !isOnlyLowTriangle)
			m_map[format[i]] = i;

	// Создаем сортированный формат, чтобы по нему двигаться
	if (isOnlyLowTriangle) {
		for (int i = 0; i < format.size(); ++i)
			if (format[i] < 0)
				m_sorted_format.push_back(format[i]);
	} else 
		m_sorted_format = format;
	std::sort(m_sorted_format.begin(), m_sorted_format.end());

	line = 0;
	pos = 0;
	start = m_sorted_format.size() - 1;
	end = start;

	// Находим, с какой диагонали начинается текущая строка
	for (int i = 0; i < m_sorted_format.size(); ++i) {
		if (isIntersectDiagLine(n, line, m_sorted_format[i])) {
			start = i;
			break;
		}
	}

	// Находим на какой диагонали кончается текущая строка
	for (int i = 0; i < m_sorted_format.size(); ++i) {
		int j = m_sorted_format.size() - i - 1;
		if (isIntersectDiagLine(n, line, m_sorted_format[j])) {
			end = j;
			break;
		}
	}

	calcPos();
}

//-----------------------------------------------------------------------------
matrix_diagonal_line_iterator& matrix_diagonal_line_iterator::operator++() {
	if (!m_isEnd) {
		if (m_isLineEnd) {
			// Сдвигаемся на одну строку
			line++;

			// Определяем какие диагонали пересекают эту строку
			if (start != 0)
				if (isIntersectDiagLine(n, line, m_sorted_format[start-1]))
					start = start-1;

			if (end != 0)
				if (!isIntersectDiagLine(n, line, m_sorted_format[end]))
					if (start != end)
						end = end-1;

			m_isLineEnd = false;
			if (line == n)
				m_isEnd = true;

			pos = 0;
			calcPos();
		} else {
			// Сдвигаемся на один столбец
			pos++;
			calcPos();
		}
	}

	return *this;
}

//-----------------------------------------------------------------------------
matrix_diagonal_line_iterator matrix_diagonal_line_iterator::operator++(int) {
	return operator++();
}

//-----------------------------------------------------------------------------
bool matrix_diagonal_line_iterator::isLineEnd(void) const {
	return m_isLineEnd;
}

//-----------------------------------------------------------------------------
bool matrix_diagonal_line_iterator::isEnd(void) const {
	return m_isEnd;
}

//-----------------------------------------------------------------------------
void matrix_diagonal_line_iterator::calcPos(void) {
	// Вычисляет все текущие положения согласно переменным start, pos и формату
	if (!isIntersectDiagLine(n, line, m_sorted_format[end]) || (start + pos > end)) {
		m_isLineEnd = true;
		i = line;
		j = 0;
		d = 0;
		di = 0;
		dn = 0;
	} else {
		i = line;
		d = m_sorted_format[start + pos];
		dn = m_map[d];
		di = getStartIntersectDiagLine(n, i, d);
		j = getRowIntersectDiagLine(n, i, d);
	}
}

//-----------------------------------------------------------------------------
bool matrix_diagonal_line_iterator::isIntersectDiagLine(int n, int i, int d) {
	if (d == 0)
		return true;

	if (d < 0)
		return (i+d >= 0);

	if (d > 0)
		return (i < MatrixDiagonal::calcDiagonalSize(n, d));
}

//-----------------------------------------------------------------------------
int matrix_diagonal_line_iterator::getStartIntersectDiagLine(int n, int i, int d) {
	if (d == 0)
		return i+d;

	if (d < 0)
		return i+d;

	if (d > 0)
		return i;
}

//-----------------------------------------------------------------------------
int matrix_diagonal_line_iterator::getRowIntersectDiagLine(int n, int i, int d) {
	if (d == 0)
		return i+d;

	if (d < 0)
		return i+d;

	if (d > 0)
		return i+d;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
std::vector<int> makeSevenDiagonalFormat(int n, int m, int k) {
	std::vector<int> result;

	if (1+m+k >= n)
		throw std::exception();

	result.push_back(0);

	result.push_back(1);
	result.push_back(1+m);
	result.push_back(1+m+k);

	result.push_back(-1);
	result.push_back(-1-m);
	result.push_back(-1-m-k);

	return result;
}

//-----------------------------------------------------------------------------
std::vector<int> generateRandomFormat(int n, int diagonalsCount) {
	std::vector<int> result;
	result.push_back(0);

	// Создаем массив всех возможных диагоналей
	std::vector<int> diagonals;
	for (int i = -n+1; i < n; ++i)
		if (i != 0)
			diagonals.push_back(i);

	diagonalsCount = std::min<int>(diagonals.size(), diagonalsCount);

	// Заполняем результат случайными диагоналями из этого массива
	for (int i = 0; i < diagonalsCount; ++i) {
		int pos = intRandom(0, diagonals.size());
		result.push_back(diagonals[pos]);
		diagonals.erase(diagonals.begin() + pos);
	}

	return result;
}

//-----------------------------------------------------------------------------
void generateDiagonalMatrix(int n, int min, int max, std::vector<int> format, MatrixDiagonal& result) {
	result = MatrixDiagonal(n, format);
	for (int i = 0; i < result.getDiagonalsCount(); ++i) {
		auto mit = result.posBegin(i);
		for (auto it = result.begin(i); it != result.end(i); ++it, ++mit)
			(*it) = intRandom(min, max);
	}
}

//-----------------------------------------------------------------------------
bool mul(const MatrixDiagonal& a, const Vector& x, Vector& y) {
	if (x.size() != a.dimension())
		return false;

	y.resize(x.size());

	// Зануление результата
	y.zero();
	for (int i = 0; i < a.getDiagonalsCount(); ++i) {
		auto mit = a.posBegin(i);
		for (auto it = a.begin(i); it != a.end(i); ++it, ++mit)
			y(mit.i) += (*it) * x(mit.j);
	}

	return true;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
SolverSLAE_Iterative::SolverSLAE_Iterative() : 
	w(1), 
	isLog(false), 
	log(std::cout), 
	start(), 
	epsilon(0.00001), 
	maxIterations(100) {
}

//-----------------------------------------------------------------------------
int SolverSLAE_Iterative::jacobi(const MatrixDiagonal& a, const Vector& y, Vector& x) {
	if (a.dimension() != y.size() || start.size() != y.size())
		throw std::exception();

	Vector x1(y.size()); // x^(k+1)
	Vector f(y.size()); // A*x^(k)

	// Считаем норму матрицы: ее максимальный элемент по модулю
	real yNorm = y.getMax();

	// Цикл по итерациям
	x = start;
	real error = epsilon + 1;
	int i = 0;
	for (; i < maxIterations && error > epsilon; ++i) {
		// Умножем матрицу на решение
		//mul(a, x, x1);
		mulUpperTriangle(a, x, x1);
		for (int i = 0; i < a.getDiagonalsCount(); ++i) 
			if (a.getDiagonalPos(i) < 0) {
				auto mit = a.posBegin(i);
				for (auto it = a.begin(i); it != a.end(i); ++it, ++mit)
					x1(mit.i) += (*it) * x(mit.j);
			}

		// x^(k+1) = x^k + w/a(i, i) * x^(k+1)
		auto it = a.begin(0);
		for (int i = 0; i < x1.size(); ++i, ++it)
			x(i) += w / (*it) * (y(i)-x1(i));

		// Считаем невязку
		mul(a, x, f);

		// Считаем относительную погрешность
		real fNorm = f.getMax();
		error = fabs(yNorm - fNorm) / yNorm;

		// Выводим данные
		if (isLog)
			log << i << "\t" << std::scientific << std::setprecision(3) << error << std::endl;
	}

	return i;
}

//-----------------------------------------------------------------------------
int SolverSLAE_Iterative::seidel(const MatrixDiagonal& a, const Vector& y, Vector& x) {
	if (a.dimension() != y.size() || start.size() != y.size())
		throw std::exception();

	Vector x1(y.size()); // x^(k+1)
	Vector f(y.size()); // A*x^(k)

	// Считаем норму матрицы: ее максимальный элемент по модулю
	real yNorm = y.getMax();

	// Цикл по итерациям
	x = start;
	real error = epsilon + 1;

	int i = 0;
	for (; i < maxIterations && error > epsilon; ++i) {
		// Умножем матрицу на решение
		x1.zero();
		mulUpperTriangle(a, x, x1);

		// Проходим по нижнему треугольнику и считаем все параметры
		matrix_diagonal_line_iterator mit(a.dimension(), a.getFormat(), true);
		for (; !mit.isEnd(); ++mit) {
			for (; !mit.isLineEnd(); ++mit)
				x1(mit.i) += a.begin(mit.dn)[mit.di] * x(mit.j);
			x(mit.i) = x(mit.i) + w/a.begin(0)[mit.i] * (y(mit.i) - x1(mit.i));
		}

		// Считаем невязку
		mul(a, x, f);

		// Считаем относительную погрешность
		real fNorm = f.getMax();
		error = fabs(yNorm - fNorm) / yNorm;

		// Выводим данные
		if (isLog)
			log << i << "\t" << std::scientific << std::setprecision(3) << error << std::endl;
	}

	return i;
}

//-----------------------------------------------------------------------------
void SolverSLAE_Iterative::mulUpperTriangle(const MatrixDiagonal& a, const Vector& x, Vector& y) {
	if (x.size() != a.dimension())
		throw std::exception();

	y.resize(x.size());
	y.zero();
	for (int i = 0; i < a.getDiagonalsCount(); ++i) 
		if (a.getDiagonalPos(i) >= 0) {
			auto mit = a.posBegin(i);
			for (auto it = a.begin(i); it != a.end(i); ++it, ++mit)
				y(mit.i) += (*it) * x(mit.j);
		}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
SolverSLAE_Iterative_matrix::SolverSLAE_Iterative_matrix() :
	w(1),
	isLog(false),
	log(std::cout),
	start(),
	epsilon(0.00001),
	maxIterations(100) {
}

//-----------------------------------------------------------------------------
int SolverSLAE_Iterative_matrix::jacobi(const Matrix& a, const Vector& y, Vector& x) {
	if (a.width() != a.height() || a.width() != y.size() || start.size() != y.size())
		throw std::exception();

	Vector x1(y.size()); // x^(k+1)
	Vector f(y.size()); // A*x^(k)

						// Считаем норму матрицы: ее максимальный элемент по модулю
	real yNorm = y.getMax();

	// Цикл по итерациям
	x = start;
	real error = epsilon + 1;
	int iter = 0;
	for (; iter < maxIterations && error > epsilon; ++iter) {
		for (int i = 0; i < a.height(); ++i) {
			sumreal sum = 0;
			for (int j = 0; j < a.height(); ++j)
				sum += a(i, j) * x(j);
			x1(i) = x(i) + w / a(i, i) * (y(i) - sum);
		}

		x = x1;

		// Считаем невязку
		mul(a, x, f);

		// Считаем относительную погрешность
		real fNorm = f.getMax();
		error = fabs(yNorm - fNorm) / yNorm;

		// Выводим данные
		if (isLog)
			log << iter << "\t" << std::scientific << std::setprecision(3) << error << std::endl;
	}

	return iter;
}

//-----------------------------------------------------------------------------
int SolverSLAE_Iterative_matrix::seidel(const Matrix& a, const Vector& y, Vector& x) {
	if (a.width() != a.height() || a.width() != y.size() || start.size() != y.size())
		throw std::exception();

	Vector f(y.size()); // A*x^(k)

						// Считаем норму матрицы: ее максимальный элемент по модулю
	real yNorm = y.getMax();

	// Цикл по итерациям
	x = start;
	real error = epsilon + 1;
	int iter = 0;
	for (; iter < maxIterations && error > epsilon; ++iter) {
		for (int i = 0; i < a.height(); ++i) {
			sumreal sum = 0;
			for (int j = 0; j < a.height(); ++j)
				sum += a(i, j) * x(j);
			x(i) = x(i) + w / a(i, i) * (y(i) - sum);
		}

		// Считаем невязку
		mul(a, x, f);

		// Считаем относительную погрешность
		real fNorm = f.getMax();
		error = fabs(yNorm - fNorm) / yNorm;

		// Выводим данные
		if (isLog)
			log << iter << "\t" << std::scientific << std::setprecision(3) << error << std::endl;
	}

	return iter;
}