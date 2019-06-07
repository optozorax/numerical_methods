#pragma once

#include <string>
#include <vector>
#include "matrix.h"
#include "common.h"
#include "vector.h"

//-----------------------------------------------------------------------------
/** Симметричная матрица в профильном формате. */
class MatrixProfileSymmetric
{
public:
	MatrixProfileSymmetric();
	MatrixProfileSymmetric(const Matrix& a);

	void loadFromFile(std::string fileName);
	void saveToFile(std::string fileName) const;

	void generate(int n, int min, int max, int maxDistanceToDiagonal);

	void toDenseMatrix(Matrix& dense) const;

	//-------------------------------------------------------------------------
	void negate(void);

	//-------------------------------------------------------------------------
	int size(void) const;

	myreal& getDiagonalElement(int n);
	const myreal& getDiagonalElement(int n) const;

	//-------------------------------------------------------------------------
	int getLineFirstElementPos(int lineNo) const;
	int getLineSize(int lineNo) const;

	// Нумерация начинается с 0, но дается элемент с индексом в плотной матрице как getLineFirstElementPos
	myreal& getLineElement(int lineNo, int elemNo);
	const myreal& getLineElement(int lineNo, int elemNo) const;

	std::vector<myreal>::iterator getLineFirstElement(int lineNo);
	std::vector<myreal>::const_iterator getLineFirstElement(int lineNo) const;

private:
	std::vector<myreal> di;
	std::vector<myreal> al;
	std::vector<int> ai;
};

//-----------------------------------------------------------------------------
bool mul(const MatrixProfileSymmetric& a, const Vector& x, Vector& y);

void calcLDL(MatrixProfileSymmetric& a_l);

void calcGaussianReverseOrder(const MatrixProfileSymmetric& l, Vector& y_x);
void calcGaussianFrontOrder(const MatrixProfileSymmetric& l, Vector& y_x);
void calcGaussianCentralOrder(const MatrixProfileSymmetric& d, Vector& y_x);

void solveSLAE_by_LDL(MatrixProfileSymmetric& a, Vector& y_x);