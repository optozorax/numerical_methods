#pragma once

#include <functional>
#include <vector>

using namespace std;

typedef vector<double> xn_t; // x in R^n вектор
typedef vector<xn_t> matrix;
typedef function<double(const xn_t&)> fn1_t; // f : R^n -> R одна функция нелинейной системы
typedef vector<fn1_t> fnm_t; // f : R^n -> R^m нелинейная система
typedef function<matrix(const xn_t&)> jnm_t; // j : R^n -> R^(m*n) матрица якоби

double length(const xn_t& x);

double calc_partial_derivative_numeric(const fn1_t& f, const xn_t& x, int i);
matrix calc_jacobi_matrix_numeric(const fnn_t& f, const xn_t& x);
jnm_t calc_jacobi_matrix_numeric_functon(const fnn_t& f);

struct solved_t
{
	int iterations;
	double relative_residual;
	xn_t point;
};

solved_t solve(
	const jnm_t& j, 
	const fnm_t& f,
	const xn_t& x_0,
	int maxiter, 
	double eps,
	bool is_log
);

jnm_t square_cast_1(const jnm_t& j, const fnm_t& f);
jnm_t square_cast_3(const jnm_t& j, const fnm_t& f);

// bind(calc_jacobi_matrix_numeric, f, placeholders::_1);