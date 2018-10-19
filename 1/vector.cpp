#include <fstream>
#include <iomanip>
#include "vector.h"

//-----------------------------------------------------------------------------
Vector::Vector() {
}

//-----------------------------------------------------------------------------
Vector::Vector(int size, real fill) : mas(size, fill) {
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
	fout.precision(std::numeric_limits<real>::digits10);
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
void Vector::zero(void) {
	for (auto& i : mas)
		i = 0;
}

//-----------------------------------------------------------------------------
real Vector::getMax(void) const {
	real max = mas[0];
	for (auto& i : mas)
		if (i > max)
			max = i;
	return max;
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