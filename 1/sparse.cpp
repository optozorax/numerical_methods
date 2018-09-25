#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "sparse.h"

//-----------------------------------------------------------------------------
MatrixProfileSymmetric::MatrixProfileSymmetric() {
}

//-----------------------------------------------------------------------------
MatrixProfileSymmetric::MatrixProfileSymmetric(const Matrix& a) {
	if (!a.isSymmetric())
		throw std::exception();

	di.resize(a.height(), 0);
	ai.resize(a.height() + 1, 0);
	for (int i = 0; i < a.height(); ++i) {
		di[i] = a(i, i);
		int zeroCount = 0;
		for (; zeroCount < i; ++zeroCount)
			if (a(i, zeroCount) != 0)
				break;
		int count = i - zeroCount;
		ai[i+1] = ai[i] + count;
		for (int j = 0; j < count; ++j)
			al.push_back(a(i, j + zeroCount));
	}
}

//-----------------------------------------------------------------------------
void MatrixProfileSymmetric::loadFromFile(std::string fileName) {
	std::ifstream fin(fileName);
	int n;

	fin >> n;
	di.resize(n, 0);
	for (int i = 0; i < n; ++i)
		fin >> di[i];

	fin >> n;
	al.resize(n, 0);
	for (int i = 0; i < n; ++i)
		fin >> al[i];

	fin >> n;
	ai.resize(n, 0);
	for (int i = 0; i < n; ++i)
		fin >> ai[i];

	fin.close();
}

//-----------------------------------------------------------------------------
void MatrixProfileSymmetric::saveToFile(std::string fileName) const {
	std::ofstream fout(fileName);

	fout.precision(std::numeric_limits<real>::max_digits10);
	int w = std::numeric_limits<real>::digits10 + 6;

	fout << di.size() << std::endl;
	for (const auto& i : di)
		fout << std::setw(w) << i;
	fout << std::endl;

	fout << ai.size() << std::endl;
	for (const auto& i : ai)
		fout << std::setw(w) << i;
	fout << std::endl;

	fout << al.size() << std::endl;
	for (const auto& i : al)
		fout << std::setw(w) << i;
	fout << std::endl;

	fout.close();
}

//-----------------------------------------------------------------------------
void MatrixProfileSymmetric::toDenseMatrix(Matrix& dense) const {
	if (dense.width() != size() && dense.height() != size())
		dense.resize(size(), size());

	for (int i = 0; i < size(); ++i) {
		int iLineStart = getLineFirstElementPos(i);
		int iLineSize = getLineSize(i);
		for (int j = 0; j < iLineSize; ++j) {
			dense(i, iLineStart + j) = getLineElement(i, j);
			dense(iLineStart + j, i) = getLineElement(i, j);
		}
		dense(i, i) = getDiagonalElement(i);
	}
}

//-----------------------------------------------------------------------------
void MatrixProfileSymmetric::negate(void) {
	for (auto& i : al)
		i = -i;
}

//-----------------------------------------------------------------------------
void MatrixProfileSymmetric::generate(int n, int min, int max, int maxDistanceToDiagonal) {
	di.resize(n, 0);
	ai.resize(n+1, 0);
	al.clear();
	al.reserve(maxDistanceToDiagonal * n);

	for (int i = 0; i < n; ++i) {
		di[i] = intRandom(min, max);
		int maxLineSize = i;
		int lineSize = intRandom(0, std::min(maxLineSize, maxDistanceToDiagonal));
		ai[i+1] = ai[i] + lineSize;
		for (int j = 0; j < lineSize; ++j)
			al.push_back(intRandom(min, max));
	}
}

//-----------------------------------------------------------------------------
int MatrixProfileSymmetric::size(void) const {
	return di.size();
}

//-----------------------------------------------------------------------------
real& MatrixProfileSymmetric::getDiagonalElement(int n) {
	return di[n];
}

//-----------------------------------------------------------------------------
const real& MatrixProfileSymmetric::getDiagonalElement(int n) const {
	return di[n];
}

//-----------------------------------------------------------------------------
int MatrixProfileSymmetric::getLineFirstElementPos(int lineNo) const {
	return lineNo - getLineSize(lineNo);
}

//-----------------------------------------------------------------------------
int MatrixProfileSymmetric::getLineSize(int lineNo) const {
	return ai[lineNo+1] - ai[lineNo];
}

//-----------------------------------------------------------------------------
real& MatrixProfileSymmetric::getLineElement(int lineNo, int elemNo) {
	return al[ai[lineNo] + elemNo];
}

//-----------------------------------------------------------------------------
const real& MatrixProfileSymmetric::getLineElement(int lineNo, int elemNo) const {
	return al[ai[lineNo] + elemNo];
}

//-----------------------------------------------------------------------------
std::vector<real>::iterator MatrixProfileSymmetric::getLineFirstElement(int lineNo) {
	return al.begin() + ai[lineNo];
}

//-----------------------------------------------------------------------------
std::vector<real>::const_iterator MatrixProfileSymmetric::getLineFirstElement(int lineNo) const {
	return al.begin() + ai[lineNo];
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
Vector::Vector() {
}

//-----------------------------------------------------------------------------
Vector::Vector(int size, real fill) : mas(size, 0) {
}

//-----------------------------------------------------------------------------
Vector::Vector(const Matrix& a) {
	if (a.width() == 1) {
		resize(a.height());
		for (int i = 0; i < a.height(); ++i)
			mas[i] = a(i, 0);
	} else if (a.height() == 1) {
		resize(a.width());
		for (int i = 0; i < a.width(); ++i)
			mas[i] = a(0, i);
	} else
		throw std::exception();
}

//-----------------------------------------------------------------------------
void Vector::toDenseMatrix(Matrix& dense, bool isVertical) const {
	if (isVertical) {
		dense.resize(1, mas.size());
		for (int i = 0; i < mas.size(); ++i)
			dense(i, 0) = mas[i];
	} else {
		dense.resize(mas.size(), 1);
		for (int i = 0; i < mas.size(); ++i)
			dense(0, i) = mas[i];
	}
}

//-----------------------------------------------------------------------------
void Vector::loadFromFile(std::string fileName) {
	std::ifstream fin(fileName);
	int n;
	fin >> n;
	resize(n);
	for (int i = 0; i < n; ++i)
		fin >> mas[i];
	fin.close();
}

//-----------------------------------------------------------------------------
void Vector::saveToFile(std::string fileName) const {
	std::ofstream fout(fileName);
	fout.precision(std::numeric_limits<real>::max_digits10);
	fout << mas.size() << std::endl;
	for (const auto& i : mas)
		fout << i << std::endl;
	fout.close();
}

//-----------------------------------------------------------------------------
void Vector::resize(int n, real fill) {
	if (mas.size() != n)
		mas.resize(n, fill);
}

//-----------------------------------------------------------------------------
void Vector::negate(void) {
	for (auto& i : mas)
		i = -i;
}

//-----------------------------------------------------------------------------
void Vector::generate(int n, int min, int max) {
	resize(n);
	for (int i = 0; i < n; ++i)
		operator()(i) = intRandom(min, max);
}

//-----------------------------------------------------------------------------
void Vector::generate(int n) {
	resize(n);
	for (int i = 0; i < n; ++i)
		operator()(i) = i+1;
}

//-----------------------------------------------------------------------------
int Vector::size(void) const {
	return mas.size();
}

//-----------------------------------------------------------------------------
real& Vector::operator()(int i) {
	return mas[i];
}

//-----------------------------------------------------------------------------
const real& Vector::operator()(int i) const {
	return mas[i];
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool mul(const MatrixProfileSymmetric& a, const Vector& x, Vector& y) {
	if (x.size() != a.size())
		return false;

	y.resize(x.size());

	// Зануление результата
	for (int i = 0; i < y.size(); ++i)
		y(i) = 0;

	// Умножение элементов из матрицы L
	for (int i = 0; i < a.size(); ++i) {
		int iLineStart = a.getLineFirstElementPos(i);
		int iLineSize = a.getLineSize(i);

		int j = iLineStart;
		for (int k = 0; k < iLineSize; ++j, ++k) {
			real elem = a.getLineElement(i, k);
			y(i) += elem * x(j);
			y(j) += elem * x(i);
		}
	}

	// Умножение диагональных элементов на вектор
	for (int i = 0; i < a.size(); ++i)
		y(i) += a.getDiagonalElement(i) * x(i);

	return true;
}

//-----------------------------------------------------------------------------
bool sum(const Vector& a, const Vector& b, Vector& result) {
	if (a.size() != b.size())
		return false;

	if (result.size() != a.size())
		result.resize(a.size());

	for (int i = 0; i < result.size(); ++i)
		result(i) = a(i) + b(i);

	return true;
}

//-----------------------------------------------------------------------------
sumreal sumAllElementsAbs(const Vector& a) {
	sumreal sum = 0;
	for (int i = 0; i < a.size(); ++i)
		sum += fabs(a(i));

	return sum;
}

//-----------------------------------------------------------------------------
void calcLDL(MatrixProfileSymmetric& a_l) {
	// Типо обращение к диагональному элементу
	auto d = [&a_l] (int i) -> real& {
		return a_l.getDiagonalElement(i);
	};

	// Перебираем все строки матрицы A
	for (int i = 0; i < a_l.size(); ++i) {
		// Сумма для диагонального элемента
		sumreal dsum = 0;

		int iLineStart = a_l.getLineFirstElementPos(i);
		int iLineSize = a_l.getLineSize(i);

		// Считаем элементы матрицы L
		for (int j = iLineStart; j < i; ++j) {
			real& currentElem = a_l.getLineElement(i, j - iLineStart);

			int jLineStart = a_l.getLineFirstElementPos(j);
			int jLineSize = a_l.getLineSize(j);

			int offset = std::max(iLineStart, jLineStart);
			int iLineOffset = offset - iLineStart;
			int jLineOffset = offset - jLineStart;

			int end = std::min(iLineStart + iLineSize, jLineStart + jLineSize);

			int count = end - offset;

			sumreal sum = 0;
			for (int k = 0; k < count; ++k)
				sum += d(offset + k) * 
					   a_l.getLineElement(j, k + jLineOffset) * 
					   a_l.getLineElement(i, k + iLineOffset);

			currentElem = (currentElem - sum) / d(j);

			dsum += d(j) * currentElem * currentElem;
		}

		// Считаем диагональный элемент
		d(i) = d(i) - dsum;
	}
}

//-----------------------------------------------------------------------------
void calcGaussianReverseOrder(const MatrixProfileSymmetric& a, Vector& y_x) {
	for (int i = y_x.size() - 1; i >= 0; --i) {
		sumreal sum = 0;
		int iLineStart = a.getLineFirstElementPos(i);
		int iLineSize = a.getLineSize(i);
		y_x(i) = y_x(i) - sum;
		for (int j = iLineStart; j < i; j++)
			y_x(j) -= y_x(i) * a.getLineElement(i, j - iLineStart);
	}
}

//-----------------------------------------------------------------------------
void calcGaussianFrontOrder(const MatrixProfileSymmetric& a, Vector& y_x) {
	for (int i = 0; i < y_x.size(); ++i) {
		int iLineStart = a.getLineFirstElementPos(i);
		int iLineSize = a.getLineSize(i);

		sumreal sum = 0;
		for (int j = 0; j < iLineSize; ++j)
			sum += a.getLineElement(i, j) * y_x(iLineStart + j);
		y_x(i) = y_x(i) - sum;
	}
}

//-----------------------------------------------------------------------------
void calcGaussianCentralOrder(const MatrixProfileSymmetric& d, Vector& y_x) {
	for (int i = 0; i < y_x.size(); ++i)
		y_x(i) = y_x(i) / d.getDiagonalElement(i);
		// TODO что если элемент равен нулю?
}

//-----------------------------------------------------------------------------
void solveSLAE_by_LDL(MatrixProfileSymmetric& a, Vector& y_x) {
	y_x.resize(a.size());

	calcLDL(a);
	calcGaussianFrontOrder(a, y_x);
	calcGaussianCentralOrder(a, y_x);
	calcGaussianReverseOrder(a, y_x);
}