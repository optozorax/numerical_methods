#pragma once

#include <iostream>
#include <functional>
#include <vector>

using namespace std;

typedef vector<double> xn_t; // x in R^n вектор
typedef vector<xn_t> matrix_t;
typedef function<double(const xn_t&)> fn_f; // f : R^n -> R одна функция нелинейной системы
typedef vector<fn_f> fnm_f; // f : R^n -> R^m нелинейная система
typedef function<matrix_t(const xn_t&)> jnm_f; // j : R^n -> R^(m*n) матрица якоби
typedef pair<matrix_t, xn_t> sle_t; // Тип СЛАУ
typedef function<sle_t(const xn_t&)> sle_f; // Функция, которая возвращает СЛАУ

double length(const xn_t& x);
xn_t operator+(const xn_t& a, const xn_t& b);
xn_t operator*(const xn_t& a, double b);
xn_t operator*(double b, const xn_t& a);
ostream& operator<<(ostream& out, const xn_t& v);

void solve_gauss(const matrix_t& a, const xn_t& b, xn_t& dx);

double calc_partial_derivative_numeric(const fn_f& f, const xn_t& x, int i);
matrix_t calc_jacobi_matrix_numeric(const fnm_f& f, const xn_t& x);
jnm_f calc_jacobi_matrix_numeric_functon(const fnm_f& f);

xn_t calc_vector_function(const fnm_f& f, const xn_t& x);
sle_f get_sle_function(const jnm_f& j, const fnm_f& f);
sle_f square_cast_1(const sle_f& s);
sle_f square_cast_2(const sle_f& s);
sle_f square_cast_3(const sle_f& s);

struct solved_t
{
	int iterations;
	double residual;
	xn_t point;
	vector<xn_t> process;
};

solved_t solve(
	const sle_f& s,
	const fnm_f& f,
	const xn_t& x_0,
	int maxiter, 
	double eps,
	bool is_log
);