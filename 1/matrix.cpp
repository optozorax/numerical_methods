#include <fstream>
#include <iomanip>
#include "matrix.h"

//-----------------------------------------------------------------------------
bool isNear(real a, real b) {
	if (a != 0) {
		if (fabs(a - b)/a > 0.0001)
			return false;
	} else {
		if (fabs(b) > 0.0001)
			return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
Matrix::Matrix(int n, int m, real fill) : m_matrix(m, std::vector<real>(n, fill)), m_n(n), m_m(m) {
}

//-----------------------------------------------------------------------------
void Matrix::loadFromFile(std::string fileName) {
	std::ifstream fin(fileName);

	m_matrix.clear();
	int n, m;
	fin >> n >> m;
	resize(n, m);
	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j < width(); ++j) {
			fin >> operator()(i, j);
		}
	}

	fin.close();
}

//-----------------------------------------------------------------------------
void Matrix::saveToFile(std::string fileName) const {
	std::ofstream fout(fileName);

	fout << m_n << '\t' << m_m << std::endl;

	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j < width(); ++j) {
			fout << std::setprecision(3) << std::fixed << std::setw(7) << operator()(i, j) << '\t';
		}
		fout << std::endl;
	}

	fout.close();
}

//-----------------------------------------------------------------------------
void Matrix::getFromVector(int n, int m, const std::vector<real>& data) {
	m_n = n;
	m_m = m;
	resize(n, m);
	for (int i = 0; i < data.size(); ++i)
		operator()(i % n, i / n) = data[i];
}

//-----------------------------------------------------------------------------
void Matrix::resize(int n, int m, real fill) {
	if (m_n != n || m_m != m) {
		m_n = n;
		m_m = m;
		m_matrix.resize(m_m, std::vector<real>(m_n, fill));
	}
}

//-----------------------------------------------------------------------------
void Matrix::negate(void) {
	for (auto& i : m_matrix) {
		for (auto& j : i) {
			j = -j;
		}
	}
}

//-----------------------------------------------------------------------------
bool Matrix::isSymmetric(void) const {
	if (height() != width())
		return false;

	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j <= i; ++j) {
			const real& a = operator()(i, j);
			const real& b = operator()(j, i);
			if (!isNear(a, b))
				return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool Matrix::isLowerTriangular(void) const {
	if (height() != width())
		return false;

	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j < i; ++j) {
			if (operator()(j, i) != 0)
				return false;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
bool Matrix::isUpperTriangular(void) const {
	if (height() != width())
		return false;

	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j < i; ++j) {
			if (operator()(i, j) != 0)
				return false;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
bool Matrix::isDiagonal(void) const {
	if (height() != width())
		return false;

	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j < i; ++j) {
			if (operator()(j, i) != 0 && operator()(i, j) != 0)
				return false;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
bool Matrix::isDiagonalIdentity(void) const {
	if (height() != width())
		return false;

	for (int i = 0; i < height(); ++i) {
		if (operator()(i, i) != 1)
			return false;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
real& Matrix::operator()(int i, int j) {
	return m_matrix[i][j];
}

//-----------------------------------------------------------------------------
const real& Matrix::operator()(int i, int j) const {
	return m_matrix[i][j];
}

//-----------------------------------------------------------------------------
int Matrix::width(void) const {
	return m_n;
}

//-----------------------------------------------------------------------------
int Matrix::height(void) const {
	return m_m;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void generateSparseSymmetricMatrix(int n, int min, int max, real percent, Matrix& result) {
	result.resize(n, n, 0);

	int count = percent * n * n;

	for (int k = 0; k < count; ++k) {
		int i = intRandom(0, n-1);
		int j = intRandom(0, i-1);
		result(i, j) = intRandom(min, max);
		result(j, i) = result(i, j);
	}

	for (int i = 0; i < n; i++)
		result(i, i) = intRandom(1, max - min);
}

//-----------------------------------------------------------------------------
void generateLMatrix(int n, int min, int max, real percent, Matrix& result) {
	result.resize(n, n, 0);

	int count = percent * n * n / 2;

	for (int k = 0; k < count; ++k) {
		int i = intRandom(0, n-1);
		int j = intRandom(0, i-1);
		result(i, j) = intRandom(min, max);
	}

	for (int i = 0; i < n; ++i) {
		result(i, i) = 1;
	}
}

//-----------------------------------------------------------------------------
void generateDiagonalMatrix(int n, int min, int max, Matrix& result) {
	result.resize(n, n, 0);

	for (int i = 0; i < n; ++i)
		result(i, i) = intRandom(min, max);
}

//-----------------------------------------------------------------------------
void generateVector(int n, int min, int max, Matrix& result) {
	result.resize(1, n, 0);

	for (int i = 0; i < n; ++i)
		result(i, 0) = intRandom(min, max);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool mul(const Matrix& a, const Matrix& b, Matrix& result) {
	// result = a * b
	if (a.width() != b.height())
		return false;

	result.resize(b.width(), a.height());

	for (int i = 0; i < b.width(); ++i) {
		for (int j = 0; j < a.height(); ++j) {
			real sum = 0;
			for (int k = 0; k < a.width(); ++k) {
				sum += a(j, k) * b(k, i);
			}
			result(j, i) = sum;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool sum(const Matrix& a, const Matrix& b, Matrix& result) {
	// result = a + b
	if (a.width() != b.width() || a.height() != b.height())
		return false;

	result.resize(a.width(), a.height());

	for (int i = 0; i < a.width(); ++i) {
		for (int j = 0; j < a.height(); ++j) {
			result(j, i) = a(j, i) + b(j, i);
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool transpose(Matrix& a) {
	// a = a^T
	if (a.height() != a.width())
		return false;

	for (int i = 0; i < a.height(); ++i) {
		for (int j = 0; j < i; ++j) {
			std::swap(a(j, i), a(i, j));
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool calcLDL(const Matrix& a, Matrix& l, Matrix& d) {
	// l * d * l^T = a
	if (!a.isSymmetric())
		return false;

	l.resize(a.width(), a.height(), 0);
	d.resize(a.width(), a.height(), 0);	

	for (int i = 0; i < a.height(); ++i) {
		// Считаем элементы матрицы L
		for (int j = 0; j < i; ++j) {
			real sum = 0;
			for (int k = 0; k < j; ++k)
				sum += d(k, k) * l(j, k) * l(i, k);

			if (fabs(d(j, j)) < 0.0001)
				l(i, j) = 0;
			else
				l(i, j) = (a(i, j) - sum) / d(j, j);
		}

		// Считаем диагональный элемент
		{
			real sum = 0;
			for (int j = 0; j < i; ++j)
				sum += d(j, j) * l(i, j) * l(i, j);
			d(i, i) = a(i, i) - sum;
		}
	}

	for (int i = 0; i < l.height(); i++)
		l(i, i) = 1;
	
	return true;
}

//-----------------------------------------------------------------------------
bool calcGaussianReverseOrder(const Matrix& l, const Matrix& y, Matrix& x) {
	// l * x = y, l - нижнетреугольная матрица
	if (!l.isLowerTriangular() || !l.isDiagonalIdentity())	
		return false;

	x.resize(1, y.height());

	for (int i = x.height() - 1; i >= 0; --i) {
		real sum = 0;
		for (int j = i; j < x.height(); ++j)
			sum += l(j, i) * x(j, 0);
		x(i, 0) = y(i, 0) - sum;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool calcGaussianFrontOrder(const Matrix& l, const Matrix& y, Matrix& x) {
	// l * x = y, l - верхнетреугольная матрица
	if (!l.isLowerTriangular() || !l.isDiagonalIdentity())	
		return false;

	x.resize(1, y.height());

	for (int i = 0; i < x.height(); ++i) {
		real sum = 0;
		for (int j = 0; j < i; ++j)
			sum += l(i, j) * x(j, 0);
		x(i, 0) = y(i, 0) - sum;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool calcGaussianCentralOrder(const Matrix& d, const Matrix& y, Matrix& x) {
	// d * x = y, d - диагональная матрица
	if (!d.isDiagonal())	
		return false;

	x.resize(1, y.height());

	for (int i = 0; i < x.height(); ++i)
		x(i, 0) = y(i, 0) / d(i, i);

	return true;
}

//-----------------------------------------------------------------------------
bool solveSLAE_by_LDL(const Matrix& a, const Matrix& y, Matrix& x) {
	// a * x = y, a - симметричная матрица
	Matrix l, d, z, w;

	if (!calcLDL(a, l, d))
		return false;

	if (!calcGaussianFrontOrder(l, y, z))
		return false;

	if (!calcGaussianCentralOrder(d, z, w))
		return false;

	if (!calcGaussianReverseOrder(l, w, x))
		return false;

	return true;
}