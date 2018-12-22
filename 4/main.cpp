#include <iostream>
#include <iomanip>

#include "logic.h"
#include "objects.h"

int main() {
	solved_t solved;

	/*auto p = two_circles(circle(0, 0, 1), circle(2, 0, 1.5));
	xn_t x_0 = {-10, -10};*/

	auto p = three_lines(
		{point(0, 0), point(1, 2)},
		{point(0, 0), point(2, 1)}, 
		{point(1, 2), point(2, 1)}
	);
	xn_t x_0 = {-3, -3};
	solved = solve(square_cast_3(get_sle_function(calc_jacobi_matrix_numeric_functon(p.second), p.second)), p.second, x_0, 100, 1e-10, true);

	/*auto p = two_circles_and_line(
		circle(0, 0, 1), circle(2, 0, 1.5),
		{point(1, 2), point(2, 1)}
	);
	xn_t x_0 = {-3, -3};
	solved = solve(square_cast_3(p.first), p.second, x_0, 100, 1e-10, true);*/

	/*auto p = three_circles(circle(0, 0, 1), circle(5, 0, 1), circle(10, 5, 1), false, true, false);
	xn_t x_0 = {0, 0, 50};*/

	cout << "Analytic jacobi matrix_t: " << endl;

	//solved = solve(p.second, p.first, x_0, 100, 1e-10, true);

	cout << setprecision(8) << fixed;
	cout << "Point: " << solved.point << endl;
	cout << "Iterations: " << solved.iterations << endl;
	cout << setprecision(2) << scientific;
	cout << "Residual: " << solved.residual << endl;

	cout << endl;

	// cout << "Numeric jacobi matrix_t: " << endl;

	// solved = solve(calc_jacobi_matrix_numeric_functon(p.first), p.first, x_0, 100, 1e-10, true);

	// cout << setprecision(8) << fixed;
	// cout << "Point: " << solved.point << endl;
	// cout << "Iterations: " << solved.iterations << endl;
	// cout << setprecision(2) << scientific;
	// cout << "Residual: " << solved.residual << endl;

	system("pause");
}