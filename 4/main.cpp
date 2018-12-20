#include <iostream>
#include <iomanip>

#include "logic.h"
#include "objects.h"

int main() {
	solved_t solved;

	/*auto p = two_circles(circle(0, 0, 1), circle(2, 0, 1.5));
	xn_t x_0 = {-10, -10};*/

	auto p = three_circles(circle(0, 0, 1), circle(5, 0, 1), circle(10, 5, 1), false, true, false);
	xn_t x_0 = {0, 0, 50};

	cout << "Analytic jacobi matrix: " << endl;

	solved = solve(p.second, p.first, x_0, 100, 1e-10, true);

	cout << setprecision(8) << fixed;
	cout << "Point: " << solved.point << endl;
	cout << "Iterations: " << solved.iterations << endl;
	cout << setprecision(2) << scientific;
	cout << "Residual: " << solved.residual << endl;

	cout << endl;

	cout << "Numeric jacobi matrix: " << endl;

	solved = solve(calc_jacobi_matrix_numeric_functon(p.first), p.first, x_0, 100, 1e-10, true);

	cout << setprecision(8) << fixed;
	cout << "Point: " << solved.point << endl;
	cout << "Iterations: " << solved.iterations << endl;
	cout << setprecision(2) << scientific;
	cout << "Residual: " << solved.residual << endl;

	system("pause");
}