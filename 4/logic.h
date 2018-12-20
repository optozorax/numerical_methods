#pragma once

#include <iostream>
#include <functional>
#include <vector>

using namespace std;

typedef vector<double> xn_t; // x in R^n вектор
typedef vector<xn_t> matrix;
typedef function<double(const xn_t&)> fn1_t; // f : R^n -> R одна функция нелинейной системы
typedef vector<fn1_t> fnm_t; // f : R^n -> R^m нелинейная система
typedef function<matrix(const xn_t&)> jnm_t; // j : R^n -> R^(m*n) матрица якоби

double length(const xn_t& x);

void solve_gauss(const matrix& a, const xn_t& b, xn_t& dx);

xn_t operator+(const xn_t& a, const xn_t& b);
xn_t operator*(const xn_t& a, double b);
xn_t operator*(double b, const xn_t& a);
ostream& operator<<(ostream& out, const xn_t& v);

double calc_partial_derivative_numeric(const fn1_t& f, const xn_t& x, int i);
matrix calc_jacobi_matrix_numeric(const fnm_t& f, const xn_t& x);
jnm_t calc_jacobi_matrix_numeric_functon(const fnm_t& f);

xn_t calc(const fnm_t& f, const xn_t& x);

struct solved_t
{
	int iterations;
	double residual;
	xn_t point;
	vector<xn_t> process;
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