#pragma once

#include <string>
#include <vector>
#include "common.h"

typedef double real;

//-----------------------------------------------------------------------------
class Matrix
{
public:
	Matrix(int n = 0, int m = 0, real fill = 0); // n - количество столбцов, m - количество строк
	void loadFromFile(std::string fileName);
	void saveToFile(std::string fileName) const;

	void getFromVector(int n, int m, const std::vector<real>& data);

	void resize(int n, int m, real fill = 0);
	void negate(void);

	bool isSymmetric(void) const;
	bool isLowerTriangular(void) const;
	bool isUpperTriangular(void) const;
	bool isDiagonal(void) const;
	bool isDiagonalIdentity(void) const;
	bool isDegenerate(void) const;

	real& operator()(int i, int j);
	const real& operator()(int i, int j) const;

	int width(void) const;
	int height(void) const;

private:
	std::vector<std::vector<real>> m_matrix;
	int m_n, m_m;
};

//-----------------------------------------------------------------------------
void generateSparseSymmetricMatrix(
	int n,
	int min, int max, 
	real percent,
	Matrix& result
);

void generateLMatrix(
	int n,
	int min, int max,
	real percent,
	Matrix& result
);

void generateDiagonalMatrix(
	int n,
	int min, int max,
	Matrix& result
);

void generateVector(
	int n,
	int min, int max,
	Matrix& result
);

//-----------------------------------------------------------------------------
real random(void);
int intRandom(int min, int max);

bool mul(const Matrix& a, const Matrix& b, Matrix& result);
bool sum(const Matrix& a, const Matrix& b, Matrix& result);

bool transpose(Matrix& a);

bool calcLDL(const Matrix& a, Matrix& l, Matrix& d);
bool calcGaussianReverseOrder(const Matrix& l, const Matrix& y, Matrix& x);
bool calcGaussianFrontOrder(const Matrix& l, const Matrix& y, Matrix& x);
bool calcGaussianCentralOrder(const Matrix& d, const Matrix& y, Matrix& x);
bool solveSLAE_by_LDL(const Matrix& a, const Matrix& y, Matrix& x);

bool solevSLAE_byGaussMethod(const Matrix& a, const Matrix& y, Matrix& x);