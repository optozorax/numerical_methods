#include <iostream>
#include <iomanip>
#include <cassert>
#include <algorithm>

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
void solve_gauss(const matrix_t& a, const xn_t& b, xn_t& dx) {
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
double calc_partial_derivative_numeric(const fn_f& f, const xn_t& x_in, int i) {
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
matrix_t calc_jacobi_matrix_numeric(const fnm_f& f, const xn_t& x) {
	matrix_t result(f.size(), xn_t(x.size()));
	for (int i = 0; i < f.size(); ++i) {
		for (int j = 0; j < x.size(); ++j) {
			result[i][j] = calc_partial_derivative_numeric(f[i], x, j);
		}
	}
	return result;
}

//-----------------------------------------------------------------------------
jnm_f calc_jacobi_matrix_numeric_functon(const fnm_f& f) {
	return bind(calc_jacobi_matrix_numeric, f, placeholders::_1);
}

//-----------------------------------------------------------------------------
xn_t calc_vector_function(const fnm_f& f, const xn_t& x) {
	xn_t result;
	for (const auto& i : f)
		result.push_back(i(x));
	return result;
}

//-----------------------------------------------------------------------------
sle_f get_sle_function(const jnm_f& j, const fnm_f& f) {
	return [j, f] (const xn_t& x) -> sle_t {
		matrix_t a = j(x);
		xn_t b = calc_vector_function(f, x);
		return {a, b};
	};
}

//-----------------------------------------------------------------------------
sle_f square_cast_1(const sle_f& s) {
	return [s] (const xn_t& x) -> sle_t {
		// Получаем значение текущей СЛАУ
		sle_t res = s(x);
		auto& A = res.first;
		auto& b = res.second;

		int count = x.size() - A[0].size();

		// Находим номера элементов, для которых dF_i(x)/dx_j минимально при всех i
		vector<pair<double, int>> b_sorted;
		for (int j = 0; j < A[0].size(); ++j) {
			double max1 = fabs(A[0][j]);
			for (int i = 1; i < A.size(); ++i)
				max1 = max(max1, fabs(A[i][j]));
			b_sorted.push_back({max1, j});
		}

		sort(b_sorted.begin(), b_sorted.end(), [] (auto a, auto b) -> bool {
			return a.first < b.second;
		});

		vector<int> mins;
		for (int i = 0; i < count; ++i)
			mins.push_back(b_sorted[i].second);
		sort(mins.begin(), mins.end(), less<int>());

		int start = mins[0];

		// Добавляем к вектору нулевые элементы
		for (int i = 0; i < mins.size(); ++i)
			b.push_back(0);	

		auto make_vec = [] (int size, int where_one) -> vector<double> {
			vector<double> result(size, 0);
			result[where_one] = 1;
			return result;
		};

		// Добавляем к матрице новые строки
		for (auto& i : mins)
			A.push_back(make_vec(x.size(), i));

		return {A, b};
	};
}

//-----------------------------------------------------------------------------
sle_f square_cast_2(const sle_f& s) {
	return [s] (const xn_t& x) -> sle_t {
		// Получаем значение текущей СЛАУ
		sle_t res = s(x);
		auto& A = res.first;
		auto& b = res.second;

		int count = b.size() - x.size() + 1;

		// Находим номера элементов, для которых f_i(x) минимальны
		vector<pair<double, int>> b_sorted;
		for (int i = 0; i < b.size(); ++i)
			b_sorted.push_back({fabs(b[i]), i});

		sort(b_sorted.begin(), b_sorted.end(), [] (auto a, auto b) -> bool {
			return a.first < b.first;
		});

		vector<int> mins;
		for (int i = 0; i < count; ++i)
			mins.push_back(b_sorted[i].second);
		sort(mins.begin(), mins.end(), less<int>());

		int start = mins[0];

		// Удаляем лишние строки
		for (int i = mins.size()-1; i > 0; --i) {
			A.erase(A.begin() + mins[i]);
			b.erase(b.begin() + mins[i]);
		}

		return {A, b};
	};
}


//-----------------------------------------------------------------------------
sle_f square_cast_3(const sle_f& s) {
	return [s] (const xn_t& x) -> sle_t {
		// Получаем значение текущей СЛАУ
		sle_t res = s(x);
		auto& A = res.first;
		auto& b = res.second;

		int count = b.size() - x.size() + 1;

		// Находим номера элементов, для которых f_i(x) минимальны
		vector<pair<double, int>> b_sorted;
		for (int i = 0; i < b.size(); ++i)
			b_sorted.push_back({fabs(b[i]), i});

		sort(b_sorted.begin(), b_sorted.end(), [] (auto a, auto b) -> bool {
			return a.first < b.first;
		});

		vector<int> mins;
		for (int i = 0; i < count; ++i)
			mins.push_back(b_sorted[i].second);
		sort(mins.begin(), mins.end(), less<int>());

		int start = mins[0];

		// Строим новую матрицу Якоби
		for (int j = 0; j < A[start].size(); ++j)
			A[start][j] = 2*A[start][j] * b[start];

		for (int i = 1; i < mins.size(); ++i) {
			for (int j = 0; j < A[mins[i]].size(); ++j) {
				A[start][j] += 2 * A[mins[i]][j] * b[mins[i]];
			}
		}

		// Строим новый вектор правой части
		b[start] = b[start] * b[start]; 
		for (int i = 1; i < mins.size(); ++i)
			b[start] += b[mins[i]] * b[mins[i]];

		// Удаляем лишние строки
		for (int i = mins.size()-1; i > 0; --i) {
			A.erase(A.begin() + mins[i]);
			b.erase(b.begin() + mins[i]);
		}

		return {A, b};
	};
}

//-----------------------------------------------------------------------------
solved_t solve(const sle_f& s, const fnm_f& f, const xn_t& x_0, int maxiter, double eps, bool is_log) {
	vector<xn_t> process;
	process.push_back(x_0);
	xn_t x_k = x_0, x_kv, dx;
	int it = 0;
	while (it < maxiter && length(calc_vector_function(f, x_k)) > eps) {
		auto sle = s(x_k);
		auto& A = sle.first;
		auto& b = sle.second;
		for (auto& i : b) i = -i; // b = -b

		#ifdef _DEBUG
		int An = A.size();
		assert(An == b.size());
		for (auto& i : A) assert(An == i.size());
		#endif

		solve_gauss(A, b, dx);

		double beta = 1;
		x_kv = x_k + beta*dx;
		while (length(calc_vector_function(f, x_kv)) > length(calc_vector_function(f, x_k))) {
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
			cout << ", Residual: " << setw(8) << length(calc_vector_function(f, x_k)) << endl;
		}
	}

	return {it, length(calc_vector_function(f, x_k)), x_k, process};
}