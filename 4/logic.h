#pragma once

#include <iostream>
#include <functional>
#include <vector>
#include "../1/matrix.h"

#define myassert(A) {if (!(A)) { std::cerr << "ERROR: " << #A << endl; throw std::exception(); }}

using namespace std;

typedef vector<double> xn_t; // x in R^n вектор
typedef vector<xn_t> matrix_t;
typedef function<double(const xn_t&)> fn_f; // f : R^n -> R одна функция нелинейной системы
typedef vector<fn_f> fnmv_f; // vector<fn_f>
typedef function<xn_t(const xn_t&)> fnm_f; // f : R^n -> R^m нелинейная система
typedef function<matrix_t(const xn_t&)> jnm_f; // j : R^n -> R^(m*n) матрица якоби
typedef pair<matrix_t, xn_t> sle_t; // Тип СЛАУ
typedef function<sle_t(const xn_t&)> sle_f; // Функция, которая возвращает СЛАУ
typedef function<sle_f(const sle_f&)> sqr_f; // Функция, преобразующая СЛАУ к квадратному виду

Matrix to(const matrix_t& a);
Matrix to(const xn_t& a);
matrix_t to_matrix(const Matrix& a);
xn_t to_vec(const Matrix& a);

double length(const xn_t& x);
xn_t operator+(const xn_t& a, const xn_t& b);
xn_t operator-(const xn_t& a, const xn_t& b);
xn_t operator*(const xn_t& a, double b);
xn_t operator*(double b, const xn_t& a);
xn_t operator*(const xn_t& a, const xn_t& b);
ostream& operator<<(ostream& out, const xn_t& v);

void solve_gauss(const matrix_t& a, const xn_t& b, xn_t& dx);
void mul_t(const matrix_t& a, const matrix_t& b, matrix_t& result);
void mul_t(const matrix_t& a, const xn_t& b, xn_t& result);

double calc_partial_derivative_numeric(const fn_f& f, const xn_t& x, int i);
matrix_t calc_jacobi_matrix_numeric(const fnmv_f& f, const xn_t& x);
jnm_f calc_jacobi_matrix_numeric_functon(const fnmv_f& f);

xn_t calc_vector_function(const fnmv_f& f, const xn_t& x);
sle_f get_sle_function(const jnm_f& j, const fnmv_f& f);
sle_f square_cast_none(const sle_f& s);
sle_f square_cast_1(const sle_f& s);
sle_f square_cast_2(const sle_f& s);
sle_f square_cast_3(const sle_f& s);
sle_f square_cast_4(const sle_f& s);

fnm_f get_f(const sle_f& s);
//function<fnmv_f(const fnmv_f&)> function_mul(const xn_t& m);
sqr_f square_cast_mul(const xn_t& x);
sqr_f composition(const sqr_f& f, const sqr_f& g);

enum exit_type_t
{
	EXIT_ITER,
	EXIT_RESIDUAL,
	EXIT_BETA,
	EXIT_STEP,
	EXIT_ERROR
};

struct solved_t
{
	int iterations;
	double residual;
	xn_t point;
	vector<xn_t> x_process;
	vector<double> beta_process;
	vector<double> residual_process;
	exit_type_t exit_type;
};

solved_t solve(
	const sle_f& s,
	const xn_t& x_0,
	int maxiter, 
	double eps,
	bool is_log
);