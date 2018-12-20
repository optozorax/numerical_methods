#include "logic.h"

//-----------------------------------------------------------------------------
double length(const xn_t& x) {
	double sum = 0;
	for (const auto& i : x)
		sum += i*i;
	return sqrt(sum);
}

//-----------------------------------------------------------------------------
double calc_partial_derivative_numeric(const fn1_t& f, const xn_t& x, int i) {
	assert(i < x.size());

	double step = x[i]*1e-9;

	auto xcopy = x;
	x[i] +=      0; double x0 = f(x); x[i] -=      0;
	x[i] +=   step; double x1 = f(x); x[i] -=   step;
	x[i] +=  -step; double x2 = f(x); x[i] -=  -step;
	x[i] += 2*step; double x3 = f(x); x[i] -= 2*step;

	double result = (-2*x2-3*x0+6*x1-x3)/(6.0*step);
	return result;
}

//-----------------------------------------------------------------------------
matrix calc_jacobi_matrix_numeric(const fnm_t& f, const xn_t& x) {
	matrix result(f.size(), xn_t(x.size()));
	for (int i = 0; i < f.size(); ++i) {
		for (int j = 0; j < x.size(); ++j) {
			result[i][j] = calc_partial_derivative_numeric(f[i], x, j);
		}
	}
	return result;
}

//-----------------------------------------------------------------------------
jnm_t calc_jacobi_matrix_numeric_functon(const fnm_t& f) {
	return bind(calc_jacobi_matrix_numeric, f, placeholders::_1);
}

//-----------------------------------------------------------------------------
solved_t solve(const jnm_t& j, const fnm_t& f, const xn_t& x_0, int maxiter, double eps, bool is_log) {
	
}

//-----------------------------------------------------------------------------
jnm_t square_cast_1(const jnm_t& j, const fnm_t& f) {

}

//-----------------------------------------------------------------------------
jnm_t square_cast_3(const jnm_t& j, const fnm_t& f) {

}