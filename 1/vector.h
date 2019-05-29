#pragma once

#include <string>
#include <vector>
#include "matrix.h"
#include "common.h"

//-----------------------------------------------------------------------------
class Vector
{
public:
	Vector();
	Vector(int size, myreal fill = 0);
	Vector(const Matrix& a);

	void loadFromFile(std::string fileName);
	void saveToFile(std::string fileName) const;

	void save(std::ostream& out) const;
	void load(std::istream& in);

	void generate(int n, int min, int max);
	void generate(int n);

	//-------------------------------------------------------------------------
	void resize(int n, myreal fill = 0);

	void negate(void);

	void zero(void);

	void toDenseMatrix(Matrix& dense, bool isVertical) const;

	//-------------------------------------------------------------------------
	int size(void) const;
	
	myreal& operator()(int i);
	const myreal& operator()(int i) const;

private:
	std::vector<myreal> mas;
};

//-----------------------------------------------------------------------------
sumreal sumAllElementsAbs(const Vector& a);
bool sum(const Vector& a, const Vector& b, Vector& result);
bool mul(const Matrix& a, const Vector& b, Vector& result);
myreal calcNorm(const Vector& a);