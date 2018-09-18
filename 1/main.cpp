#include <iostream>
#include <cstdlib>
#include "sparse.h"
#include "matrix.h"

// int countZeroElementsInLowerTriangle(const Matrix& a) {
// 	if (a.height() != a.width())
// 		return false;

// 	int count = 0;
// 	for (int i = 0; i < a.height(); ++i) {
// 		for (int j = 0; j < i; ++j) {
// 			if (fabs(a(i, j)) < 0.0001)
// 				count++;
// 			else
// 				break;
// 		}
// 	}
	
// 	return count;
// }

// int main1() {
// 	{
// 		int size = 15;

// 		Matrix l, d, a;
// 		generateLMatrix(size, 0, 10, 0.5, l);
// 		generateDiagonalMatrix(size, 0, 10, d);
// 		l.saveToFile("L.txt");
// 		d.saveToFile("D.txt");

// 		{
// 			Matrix ld, lt;
// 			mul(l, d, ld);
// 			lt = l;
// 			transpose(lt);
// 			mul(ld, lt, a);
// 			a.saveToFile("A.txt");
// 		}

// 		Matrix lNew, dNew;
// 		calcLDL(a, lNew, dNew);
// 		lNew.saveToFile("Lnew.txt");
// 		dNew.saveToFile("Dnew.txt");

// 		{
// 			Matrix ld, lt, aNew;
// 			mul(lNew, dNew, ld);
// 			lt = lNew;
// 			transpose(lt);
// 			mul(ld, lt, aNew);
// 			aNew.saveToFile("Anew.txt");

// 			Matrix sub;
// 			aNew.negate();
// 			sum(a, aNew, sub);
// 			sub.saveToFile("sub.txt");
// 		}
// 	}

// 	{
// 		int size2 = 100;
// 		Matrix a, l, d, aNew, ld, lt;
// 		generateSparseSymmetricMatrix(size2, 0, 20, 0.1, a);
// 		calcLDL(a, l, d);
// 		a.saveToFile("ABig.txt");
// 		l.saveToFile("LBig.txt");
// 		d.saveToFile("DBig.txt");

// 		mul(l, d, ld);
// 		lt = l;
// 		transpose(lt);
// 		mul(ld, lt, aNew);
// 		aNew.saveToFile("ABig2.txt");

// 		Matrix sub;
// 		aNew.negate();
// 		sum(a, aNew, sub);
// 		sub.saveToFile("SubBig.txt");

// 		std::cout << "Zeros in A: " << countZeroElementsInLowerTriangle(a) << std::endl;
// 		std::cout << "Zeros in L: " << countZeroElementsInLowerTriangle(l) << std::endl;
// 	}

// 	{
// 		int size = 100;
// 		Matrix a, x, y;
// 		generateSparseSymmetricMatrix(size, 0, 20, 0.1, a);
// 		generateVector(size, 0, 20, x);
// 		mul(a, x, y);

// 		/*Matrix l, d, z, w, y1, lt;
// 		calcLDL(a, l, d);
// 		lt = l; transpose(lt);
// 		mul(lt, x, w);
// 		mul(d, w, z);
// 		mul(l, z, y1);

// 		l.saveToFile("SLAE_L.txt");
// 		y.saveToFile("SLAE_y.txt");

// 		Matrix z1, w1, x1;
// 		calcGaussianFrontOrder(l, y, z1);
// 		calcGaussianCentralOrder(d, z, w1);
// 		calcGaussianReverseOrder(l, w, x1);*/

// 		Matrix x1;
// 		solveSLAE_by_LDL(a, y, x1);

// 		Matrix sub;
// 		x1.negate();
// 		sum(x, x1, sub);
// 		sub.saveToFile("SLAE_Sub.txt");
// 	}

// 	system("pause");
// }

int main() {
	// Test of multiplication
	Matrix denseA;
	generateSparseSymmetricMatrix(50, 1, 20, 0.1, denseA);
	//denseA.loadFromFile("sparse_matrix.txt");
	Matrix denseX;
	generateVector(50, 1, 20, denseX);
	//denseX.loadFromFile("sparse_vector.txt");

	Matrix denseY;
	mul(denseA, denseX, denseY);

	MatrixProfileSymmetric a(denseA);
	Vector x(denseX);
	Vector y(denseY);

	mul(a, x, y);

	Matrix denseY2;
	y.toDenseMatrix(denseY2, true);

	Matrix sub;
	denseY2.negate();
	sum(denseY, denseY2, sub);
	sub.saveToFile("sparse_sub.txt");

	Matrix denseL, denseD, denseLD;
	calcLDL(a);
	calcLDL(denseA, denseL, denseD);
	a.toDenseMatrix(denseLD);
	denseLD.saveToFile("sparse_ld.txt");
	denseL.saveToFile("sparse_l.txt");
	denseD.saveToFile("sparse_d.txt");

	Matrix denseX2;
	a = Matrix(denseA);
	solveSLAE_by_LDL(denseA, denseY, denseX2);
	solveSLAE_by_LDL(a, y);
	y.toDenseMatrix(denseY2, true);	
	denseY2.negate();
	sum(denseX2, denseY2, sub);
	sub.saveToFile("sparse_sub_SLAE.txt");
}