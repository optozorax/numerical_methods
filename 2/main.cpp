#include "../1/matrix.h"
#include "diagonal.h"

int main() {
	Matrix a_d;

	// Число обусловленности: 3.78802
	a_d.getFromVector(5, 5, {
		9, 0, 1, 0, 4, 
		0, 8, 0, 2, 0, 
		0, 0, 7, 0, 3, 
		1, 0, 0, 6, 0, 
		0, 2, 0, 0, 5, 
	});

	Matrix x_d;

	x_d.getFromVector(1, 5, {
		1, 2, 3, 4, 5
	});

	/*a_d.getFromVector(2, 2, {
		1, 10,
		100, 1001
	});

	Matrix x_d;

	x_d.getFromVector(1, 2, {
		1, 1
	});*/

	Matrix y_d;
	mul(a_d, x_d, y_d);

	Vector y(y_d);

	MatrixDiagonal a(a_d);

	SolverSLAE_Iterative solver;
	solver.w = 0.7;
	solver.isLog = true;
	solver.start = Vector(5, 0);
	/*solver.start = Vector(2, 0);
	solver.start(0) = 11.01;
	solver.start(1) = 0;
	solver.epsilon = 0.00001;
	solver.maxIterations = 100000;*/
	solver.epsilon = 0.000001;
	solver.maxIterations = 100;

	Vector x1;
	solver.jacobi(a, y, x1);

	system("pause");
}