#include <iostream>
#include <iomanip>
#include <algorithm>

#include "../1/matrix.h"
#include "logic.h"

//-----------------------------------------------------------------------------
Matrix to(const matrix_t& a) {
	Matrix result(a[0].size(), a.size());
	for (int i = 0; i < a.size(); ++i) {
		for (int j = 0; j < a[i].size(); ++j)
			result(i, j) = a[i][j];
	}
	return result;
}

//-----------------------------------------------------------------------------
Matrix to(const xn_t& a) {
	Matrix result(1, a.size());
	for (int i = 0; i < a.size(); ++i) {
		result(i, 0) = a[i];
	}
	return result;
}

//-----------------------------------------------------------------------------
matrix_t to_matrix(const Matrix& a) {
	matrix_t result(a.height(), vector<double>(a.width()));
	for (int i = 0; i < a.height(); ++i) {
		for (int j = 0; j < a.width(); ++j) {
			result[i][j] = a(i, j);
		}
	}
	return result;
}

//-----------------------------------------------------------------------------
xn_t to_vec(const Matrix& a) {
	xn_t result(a.height());
	for (int i = 0; i < a.height(); ++i)
		result[i] = a(i, 0);
	return result;
}

//-----------------------------------------------------------------------------
double length(const xn_t& x) {
	double sum = 0;
	for (const auto& i : x)
		sum += i*i;
	return sqrt(sum);
}

//-----------------------------------------------------------------------------
void solve_gauss(const matrix_t& a, const xn_t& b, xn_t& dx) {
	Matrix dx_m;
	solveSLAE_byGaussMethod(to(a), to(b), dx_m);
	dx = to_vec(dx_m);
}

//-----------------------------------------------------------------------------
void mul_t(const matrix_t& a, const matrix_t& b, matrix_t& result) {
	Matrix result_m;
	Matrix a_m = to(a);
	transpose(a_m);
	mul(a_m, to(b), result_m);
	result = to_matrix(result_m);
}

//-----------------------------------------------------------------------------
void mul_t(const matrix_t& a, const xn_t& b, xn_t& result) {
	Matrix result_m;
	Matrix a_m = to(a);
	transpose(a_m);
	mul(a_m, to(b), result_m);
	result = to_vec(result_m);
}

//-----------------------------------------------------------------------------
xn_t operator+(const xn_t& a, const xn_t& b) {
	myassert(a.size() == b.size());
	xn_t result(a.size());
	for (int i = 0; i < a.size(); ++i)
		result[i] = a[i] + b[i];
	return result;
}

//-----------------------------------------------------------------------------
xn_t operator-(const xn_t& a, const xn_t& b) {
	myassert(a.size() == b.size());
	xn_t result(a.size());
	for (int i = 0; i < a.size(); ++i)
		result[i] = a[i] - b[i];
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
	myassert(i < x_in.size());

	double step = 1e-9;
	if (x_in[i] != 0)
		step = x_in[i]*step;

	auto x = x_in;
	x[i] = x[i] +      0; double x0 = f(x); x[i] = x[i] -      0;
	x[i] = x[i] +   step; double x1 = f(x); x[i] = x[i] -   step;
	x[i] = x[i] -   step; double x2 = f(x); x[i] = x[i] +   step;
	x[i] = x[i] + 2*step; double x3 = f(x); x[i] = x[i] - 2*step;
	x[i] = x[i] - 2*step; double x4 = f(x); x[i] = x[i] + 2*step;

	double result = (-x3 + 8*x1 - 8*x2 + x4)/(12.0*step);
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

sle_f square_cast_none(const sle_f& s) {
	return [s] (const xn_t& x) -> sle_t {
		auto res = s(x);
		myassert(res.first[0].size() == res.first.size());
		myassert(res.first.size() == res.second.size());
		return res;
	};
}

//-----------------------------------------------------------------------------
sle_f square_cast_1(const sle_f& s) {
	return [s] (const xn_t& x) -> sle_t {
		// Получаем значение текущей СЛАУ
		sle_t res = s(x);
		auto& A = res.first;
		auto& b = res.second;

		myassert(A[0].size() > A.size());

		int count = x.size() - A.size();

		// Находим номера элементов, для которых dF_i(x)/dx_j минимально при всех i
		vector<pair<double, int>> b_sorted;
		for (int j = 0; j < A[0].size(); ++j) {
			double max1 = fabs(A[0][j]);
			for (int i = 1; i < A.size(); ++i)
				max1 = max(max1, fabs(A[i][j]));
			b_sorted.push_back({max1, j});
		}

		sort(b_sorted.begin(), b_sorted.end(), [] (auto a, auto b) -> bool {
			return a.first < b.first;
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

		myassert(A[0].size() < A.size());

		int count = b.size() - x.size();

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
		for (int i = mins.size()-1; i >= 0; --i) {
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

		myassert(A[0].size() < A.size());

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
sle_f square_cast_4(const sle_f& s) {
	return [s] (const xn_t& x) -> sle_t {
		// Получаем значение текущей СЛАУ
		sle_t res = s(x);
		auto& A = res.first;
		auto& b = res.second;

		myassert(A[0].size() < A.size());

		matrix_t AR;
		xn_t bR;
		mul_t(A, A, AR);
		mul_t(A, b, bR);

		return {AR, bR};
	};
}

//-----------------------------------------------------------------------------
solved_t solve(const sle_f& s, const fnm_f& f, const xn_t& x_0, int maxiter, double eps, bool is_log) {
	vector<xn_t> x_process;
	vector<double> beta_process;
	vector<double> residual_process;
	x_process.push_back(x_0);
	beta_process.push_back(0);
	residual_process.push_back(0);
	xn_t x_k = x_0, x_kv, dx;
	exit_type_t exit_type;

	double f_0 = length(calc_vector_function(f, x_k));
	int it = 0;
	while (true) {
		if (it > maxiter) {
			exit_type = EXIT_ITER;
			break;
		}

		if (length(calc_vector_function(f, x_k)) / f_0 < eps) {
			exit_type = EXIT_RESIDUAL;
			break;
		}
	
		auto sle = s(x_k);
		auto& A = sle.first;
		auto& b = sle.second;
		for (auto& i : b) i = -i; // b = -b

		#ifdef _DEBUG
		int An = A.size();
		myassert(An == b.size());
		for (auto& i : A) 
			myassert(An == i.size());
		#endif

		solve_gauss(A, b, dx);

		if (dx.size() == 0) {
			exit_type = EXIT_ERROR;
			break;
		}

		double beta = 1;
		x_kv = x_k + beta*dx;
		while (length(calc_vector_function(f, x_kv)) > length(calc_vector_function(f, x_k))) {
			beta /= 2.0;
			x_kv = x_k + beta * dx;
		}

		x_k = x_kv;
		it++;

		x_process.push_back(x_k);
		beta_process.push_back(beta);
		residual_process.push_back(length(calc_vector_function(f, x_k)) / f_0);

		if (is_log) {
			cout << "Iteration: " << setw(5) << it;
			cout << scientific << setprecision(2);
			cout << ", B: " << setw(8) << beta;
			cout << ", Residual: " << setw(8) << length(calc_vector_function(f, x_k)) << endl;
		}

		if (length(x_k - x_process[x_process.size() - 2]) < eps) {
			exit_type = EXIT_STEP;
			break;
		}

		if (fabs(beta) < eps) {
			exit_type = EXIT_BETA;
			break;
		}
	}

	return {it, length(calc_vector_function(f, x_k))/f_0, x_k, x_process, beta_process, residual_process, exit_type};
}