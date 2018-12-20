#include <iostream>
#include <iomanip>
#include <cassert>

#include "../1/matrix.h"
#include "logic.h"

//-----------------------------------------------------------------------------
double length(const xn_t& x) {
	double sum = 0;
	for (const auto& i : x)
		sum += i*i;
	return sqrt(sum);
}

//-----------------------------------------------------------------------------
void solve_gauss(const matrix& a, const xn_t& b, xn_t& dx) {
	Matrix a_m(a[0].size(), a.size()), b_m(1, b.size()), dx_m(1, b.size());
	for (int i = 0; i < a.size(); ++i) {
		for (int j = 0; j < a[i].size(); ++j)
			a_m(i, j) = a[i][j];
		b_m(i, 0) = b[i];
	}
	solveSLAE_byGaussMethod(a_m, b_m, dx_m);
	dx.clear();
	dx.resize(dx_m.height());
	for (int i = 0; i < dx_m.height(); ++i)
		dx[i] = dx_m(i, 0);
}

//-----------------------------------------------------------------------------
xn_t operator+(const xn_t& a, const xn_t& b) {
	assert(a.size() == b.size());
	xn_t result(a.size());
	for (int i = 0; i < a.size(); ++i)
		result[i] = a[i] + b[i];
	return result;
}

//-----------------------------------------------------------------------------
xn_t operator*(const xn_t& a, double b) {
	xn_t result(a.size());
	for (int i = 0; i < a.size(); ++i)
		result[i] = a[i] * b;
	return result;
}

//-----------------------------------------------------------------------------
xn_t operator*(double b, const xn_t& a) {
	return operator*(a, b);
}

//-----------------------------------------------------------------------------
ostream& operator<<(ostream& out, const xn_t& v) {
	out << "(";
	for (int i = 0; i < v.size()-1; ++i)
		out << v[i] << ", ";
	out << v.back() << ")";
	return out;
}

//-----------------------------------------------------------------------------
double calc_partial_derivative_numeric(const fn1_t& f, const xn_t& x_in, int i) {
	assert(i < x_in.size());

	double step = 1e-9;
	if (x_in[i] != 0)
		step = x_in[i]*step;

	auto x = x_in;
	x[i] = x[i] +      0; double x0 = f(x); x[i] = x[i] -      0;
	x[i] = x[i] +   step; double x1 = f(x); x[i] = x[i] -   step;
	x[i] = x[i] -   step; double x2 = f(x); x[i] = x[i] +   step;
	x[i] = x[i] + 2*step; double x3 = f(x); x[i] = x[i] - 2*step;

	double result = (-2*x2-3*x0+6*x1-x3)/(6.0*step);
	//double result = (x1 - x0) / step;
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
xn_t calc(const fnm_t& f, const xn_t& x) {
	xn_t result;
	for (const auto& i : f)
		result.push_back(i(x));
	return result;
}

//-----------------------------------------------------------------------------
solved_t solve(const jnm_t& j, const fnm_t& f, const xn_t& x_0, int maxiter, double eps, bool is_log) {
	vector<xn_t> process;
	process.push_back(x_0);
	xn_t x_k = x_0, x_kv, dx;
	int it = 0;
	while (it < maxiter && length(calc(f, x_k)) > eps) {
		auto A = j(x_k);
		auto b = calc(f, x_k);
		for (auto& i : b) i = -i; // b = -b
		solve_gauss(A, b, dx);

		double beta = 1;
		x_kv = x_k + beta*dx;
		while (length(calc(f, x_kv)) > length(calc(f, x_k))) {
			beta /= 2.0;
			x_kv = x_k + beta * dx;
		}

		x_k = x_kv;
		it++;

		process.push_back(x_k);

		if (is_log) {
			cout << "Iteration: " << setw(5) << it;
			cout << scientific << setprecision(2);
			cout << ", B: " << setw(8) << beta;
			cout << ", Residual: " << setw(8) << length(calc(f, x_k)) << endl;
		}
	}

	return {it, length(calc(f, x_k)), x_k, process};
}

//-----------------------------------------------------------------------------
jnm_t square_cast_1(const jnm_t& j, const fnm_t& f) {
	return {};
}

//-----------------------------------------------------------------------------
jnm_t square_cast_3(const jnm_t& j, const fnm_t& f) {
	return {};
}