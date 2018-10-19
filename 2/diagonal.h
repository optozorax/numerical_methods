#pragma once

#include <string>
#include <vector>
#include <iostream>
#include "../1/common.h"
#include "../1/vector.h"
#include "../1/matrix.h"

class MatrixDiagonal;
class matrix_diagonal_iterator;
class SolverSLAE_Iterative;

//-----------------------------------------------------------------------------
/** Матрица в диагональном формате. */
/** 0-я диагональ всегда главная диагональ. */
class MatrixDiagonal
{
public:
	typedef std::vector<real>::iterator iterator;
	typedef std::vector<real>::const_iterator const_iterator;

	//-------------------------------------------------------------------------
	MatrixDiagonal();
	MatrixDiagonal(int n, std::vector<int> format); // format[0] must be 0, because it's main diagonal
	MatrixDiagonal(const Matrix& a);

	//void loadFromFile(std::string fileName);
	//void saveToFile(std::string fileName) const;

	void toDenseMatrix(Matrix& dense) const;

	//-------------------------------------------------------------------------
	int dimension(void) const;
	int getDiagonalsCount(void) const;
	int getDiagonalSize(int diagNo) const;
	int getDiagonalPos(int diagNo) const;

	//-------------------------------------------------------------------------
	matrix_diagonal_iterator posBegin(int diagNo) const;
	matrix_diagonal_iterator posEnd(int diagNo) const;

	iterator begin(int diagNo);
	const_iterator begin(int diagNo) const;

	iterator end(int diagNo);
	const_iterator end(int diagNo) const;

	//-------------------------------------------------------------------------
	static int calcDiagonalSize(int n, int i);
private:
	std::vector<std::vector<real>> di;
	std::vector<int> fi;
	int n;
};

std::vector<int> makeSevenDiagonalFormat(int n, int m, int k);

bool mul(const MatrixDiagonal& a, const Vector& x, Vector& y);

//-----------------------------------------------------------------------------
/** Матричный "итератор" для движения по диагонали. */
struct address { int i, j; };
class matrix_diagonal_iterator : 
	std::iterator<std::random_access_iterator_tag, address>
{
public:
	matrix_diagonal_iterator(int n, int i, bool isEnd);

	matrix_diagonal_iterator& operator++();
	matrix_diagonal_iterator  operator++(int);
	matrix_diagonal_iterator& operator--();
	matrix_diagonal_iterator  operator--(int);

	bool operator==(const matrix_diagonal_iterator& b) const;
	bool operator!=(const matrix_diagonal_iterator& b) const;

	matrix_diagonal_iterator& operator+=(const ptrdiff_t& movement);
    matrix_diagonal_iterator& operator-=(const ptrdiff_t& movement);
    matrix_diagonal_iterator  operator+(const ptrdiff_t& movement);
    matrix_diagonal_iterator  operator-(const ptrdiff_t& movement);

    ptrdiff_t operator-(const matrix_diagonal_iterator& b);

	address& 		operator*();
    const address& 	operator*() const;
    address* 		operator->();
private:
	address a;

	address& operator[](int b) {}
};

//-----------------------------------------------------------------------------
class SolverSLAE_Iterative
{
public:
	SolverSLAE_Iterative();

	void jacobi(const MatrixDiagonal& a, const Vector& y, Vector& x);
	void seidel(const MatrixDiagonal& a, const Vector& y, Vector& x);

	void iterationStep(
		const MatrixDiagonal& a, 
		Vector& x1, 
		const Vector& x, 
		const Vector& y
	);

	double 			w;
	bool 			isLog;
	std::ostream& 	log;
	Vector 			start;
	double 			epsilon;
	int 			maxIterations;
};