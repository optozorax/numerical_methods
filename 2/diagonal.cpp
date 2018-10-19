#include <cmath>
#include <iostream>
#include <iomanip>
#include "diagonal.h"

//-----------------------------------------------------------------------------
MatrixDiagonal::MatrixDiagonal() : n(0) {
}

//-----------------------------------------------------------------------------
MatrixDiagonal::MatrixDiagonal(int n, std::vector<int> format) : n(n) {
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
			if (a(mit->i, mit->j) != 0) {
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
			*it = a(mit->i, mit->j);
	}
}

//-----------------------------------------------------------------------------
void MatrixDiagonal::toDenseMatrix(Matrix& dense) const {
	dense.resize(n, n, 0);

	// Обходим массив и записываем элементы
	for (int i = 0; i < getDiagonalsCount(); ++i) {
		auto mit = posBegin(i);
		for (auto it = begin(i); it != end(i); ++it, ++mit)
			dense(mit->i, mit->j) = *it;
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

typedef matrix_diagonal_iterator mdi; // Для того, чтобы код был меньше

//-----------------------------------------------------------------------------
mdi::matrix_diagonal_iterator(int n, int i, bool isEnd) {
	int count = MatrixDiagonal::calcDiagonalSize(n, i);
	if (!isEnd) {
		if (i < 0)
			a = {-i, 0};
		else
			a = {0, i};
	} else {
		if (i < 0)
			a = {-i + count, 0 + count};
		else
			a = {0 + count, i + count};
	}
}

//-----------------------------------------------------------------------------
mdi& mdi::operator++() {
	a.i++;
	a.j++;
	return *this;
}

//-----------------------------------------------------------------------------
mdi  mdi::operator++(int) {
	a.i++;
	a.j++;
	return *this;
}

//-----------------------------------------------------------------------------
mdi& mdi::operator--() {
	a.i--;
	a.j--;
	return *this;
}

//-----------------------------------------------------------------------------
mdi  mdi::operator--(int) {
	a.i--;
	a.j--;
	return *this;
}

//-----------------------------------------------------------------------------
bool mdi::operator==(const mdi& b) const {
	return b.a.i == a.i && b.a.j == a.j;
}

//-----------------------------------------------------------------------------
bool mdi::operator!=(const mdi& b) const {
	return b.a.i != a.i || b.a.j != a.j;
}

//-----------------------------------------------------------------------------
mdi& mdi::operator+=(const ptrdiff_t& movement) {
	a.i += movement;
	a.j += movement;
	return *this;
}

//-----------------------------------------------------------------------------
mdi& mdi::operator-=(const ptrdiff_t& movement) {
	a.i -= movement;
	a.j -= movement;
	return *this;
}

//-----------------------------------------------------------------------------
mdi  mdi::operator+(const ptrdiff_t& movement) {
	mdi result = *this;
	result += movement;
	return result;
}

//-----------------------------------------------------------------------------
mdi  mdi::operator-(const ptrdiff_t& movement) {
	mdi result = *this;
	result -= movement;
	return result;
}

//-----------------------------------------------------------------------------
ptrdiff_t mdi::operator-(const mdi& b) {
	return b.a.i - a.i;
}

//-----------------------------------------------------------------------------
address& mdi::operator*() {
	return a;
}

//-----------------------------------------------------------------------------
const address& mdi::operator*() const {
	return a;
}

//-----------------------------------------------------------------------------
address* mdi::operator->() {
	return &a;
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
bool mul(const MatrixDiagonal& a, const Vector& x, Vector& y) {
	if (x.size() != a.dimension())
		return false;

	y.resize(x.size());

	// Зануление результата
	y.zero();
	for (int i = 0; i < a.getDiagonalsCount(); ++i) {
		auto mit = a.posBegin(i);
		for (auto it = a.begin(i); it != a.end(i); ++it, ++mit)
			y(mit->i) += (*it) * x(mit->j);
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
void SolverSLAE_Iterative::jacobi(const MatrixDiagonal& a, const Vector& y, Vector& x) {
	Vector x1(y.size()); // x^(k+1)
	Vector f(y.size()); // A*x^(k)

	// Считаем норму матрицы: ее максимальный элемент по модулю
	real yNorm = y.getMax();

	// Цикл по итерациям
	x = start;
	real error = epsilon + 1;
	for (int i = 0; i < maxIterations && error > epsilon; ++i) {
		// Умножем матрицу на решение
		mul(a, x, x1);

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
}

//-----------------------------------------------------------------------------
void SolverSLAE_Iterative::seidel(const MatrixDiagonal& a, const Vector& y, Vector& x) {

}

//-----------------------------------------------------------------------------
void SolverSLAE_Iterative::iterationStep(const MatrixDiagonal& a, Vector& x1, const Vector& x, const Vector& y) {
	
}