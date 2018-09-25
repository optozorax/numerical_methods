#pragma once

#include <string>
#include <vector>
#include "matrix.h"
#include "common.h"

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

	real& getDiagonalElement(int n);
	const real& getDiagonalElement(int n) const;

	//-------------------------------------------------------------------------
	int getLineFirstElementPos(int lineNo) const;
	int getLineSize(int lineNo) const;

	// Нумерация начинается с 0, но дается элемент с индексом в плотной матрице как getLineFirstElementPos
	real& getLineElement(int lineNo, int elemNo);
	const real& getLineElement(int lineNo, int elemNo) const;

	std::vector<real>::iterator getLineFirstElement(int lineNo);
	std::vector<real>::const_iterator getLineFirstElement(int lineNo) const;

private:
	std::vector<real> di;
	std::vector<real> al;
	std::vector<int> ai;
};

//-----------------------------------------------------------------------------
class Vector
{
public:
	Vector();
	Vector(int size, real fill = 0);
	Vector(const Matrix& a);

	void loadFromFile(std::string fileName);
	void saveToFile(std::string fileName) const;

	void generate(int n, int min, int max);
	void generate(int n);

	//-------------------------------------------------------------------------
	void resize(int n, real fill = 0);

	void negate(void);

	void toDenseMatrix(Matrix& dense, bool isVertical) const;

	//-------------------------------------------------------------------------
	int size(void) const;
	
	real& operator()(int i);
	const real& operator()(int i) const;

private:
	std::vector<real> mas;
};

//-----------------------------------------------------------------------------
bool mul(const MatrixProfileSymmetric& a, const Vector& x, Vector& y);
bool sum(const Vector& a, const Vector& b, Vector& result);

sumreal sumAllElementsAbs(const Vector& a);

void calcLDL(MatrixProfileSymmetric& a_l);

void calcGaussianReverseOrder(const MatrixProfileSymmetric& l, Vector& y_x);
void calcGaussianFrontOrder(const MatrixProfileSymmetric& l, Vector& y_x);
void calcGaussianCentralOrder(const MatrixProfileSymmetric& d, Vector& y_x);

void solveSLAE_by_LDL(MatrixProfileSymmetric& a, Vector& y_x);