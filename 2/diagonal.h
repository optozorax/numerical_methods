#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <map>
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

	void toDenseMatrix(Matrix& dense) const;

	//-------------------------------------------------------------------------
	int dimension(void) const;
	int getDiagonalsCount(void) const;
	int getDiagonalSize(int diagNo) const;
	int getDiagonalPos(int diagNo) const;

	std::vector<int> getFormat(void) const;

	//-------------------------------------------------------------------------
	matrix_diagonal_iterator posBegin(int diagNo) const;
	matrix_diagonal_iterator posEnd(int diagNo) const;

	iterator begin(int diagNo);
	const_iterator begin(int diagNo) const;

	iterator end(int diagNo);
	const_iterator end(int diagNo) const;

	//-------------------------------------------------------------------------
	static int calcDiagonalSize(int n, int d);
private:
	std::vector<std::vector<real>> di;
	std::vector<int> fi;
	int n;
};

std::vector<int> makeSevenDiagonalFormat(int n, int m, int k);
std::vector<int> generateRandomFormat(int n, int diagonalsCount);

void generateDiagonalMatrix(
	int n,
	int min, int max,
	std::vector<int> format, 
	MatrixDiagonal& result
);

bool mul(const MatrixDiagonal& a, const Vector& x, Vector& y);

//-----------------------------------------------------------------------------
/** Матричный "итератор" для движения по диагонали. */
class matrix_diagonal_iterator
{
public:
	matrix_diagonal_iterator(int n, int i, bool isEnd);

	matrix_diagonal_iterator& operator++();
	matrix_diagonal_iterator  operator++(int);

	bool operator==(const matrix_diagonal_iterator& b) const;
	bool operator!=(const matrix_diagonal_iterator& b) const;

	matrix_diagonal_iterator& operator+=(const ptrdiff_t& movement);

    int i, j;
};

/** Матричный "итератор" для движения по строке между различными диагоналями. Может обрабатывать как всю строку, так и только нижний треугольник. */
class matrix_diagonal_line_iterator
{
public:
	matrix_diagonal_line_iterator(int n, std::vector<int> format, bool isOnlyLowTriangle);

	matrix_diagonal_line_iterator& operator++();
	matrix_diagonal_line_iterator  operator++(int);

	bool isLineEnd(void) const;
	bool isEnd(void) const;

	int i, j; // i - текущая строка,  j - текущий столбец
	int d, dn, di; // d - формат текущей диагонали, d - номер текущей диагонали, di - номер текущего элемента в диагонали
private:
	std::map<int, int> m_map;
	std::vector<int> m_sorted_format;

	int line, start, end, pos;
	int n;

	bool m_isLineEnd;
	bool m_isEnd;

	static bool isIntersectDiagLine(int n, int i, int d);
	static int getStartIntersectDiagLine(int n, int i, int d);
	static int getRowIntersectDiagLine(int n, int i, int d);

	void calcPos(void);
};

//-----------------------------------------------------------------------------
/** Класс итеративного решателя СЛАУ для диагональной матрицы. */
class SolverSLAE_Iterative
{
public:
	SolverSLAE_Iterative();

	int jacobi(const MatrixDiagonal& a, const Vector& y, Vector& x);
	int seidel(const MatrixDiagonal& a, const Vector& y, Vector& x);

	double 			w;
	bool 			isLog;
	std::ostream& 	log;
	Vector 			start;
	double 			epsilon;
	int 			maxIterations;

private:
	void mulUpperTriangle(
		const MatrixDiagonal& a, 
		const Vector& x, 
		Vector& y
	);
};

//-----------------------------------------------------------------------------
/** Класс итеративного решателя СЛАУ для плотной матрицы. */
class SolverSLAE_Iterative_matrix
{
public:
	SolverSLAE_Iterative_matrix();

	int jacobi(const Matrix& a, const Vector& y, Vector& x);
	int seidel(const Matrix& a, const Vector& y, Vector& x);

	double 			w;
	bool 			isLog;
	std::ostream& 	log;
	Vector 			start;
	double 			epsilon;
	int 			maxIterations;
};