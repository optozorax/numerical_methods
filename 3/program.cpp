#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>

#include "../1/matrix.h"
#include "../1/vector.h"

using namespace std;

ostream& operator<<(ostream& out, const vector<double>& x) {
	out << "(";
	for (int i = 0; i < x.size()-1; ++i)
		out << x[i] << ", ";
	out << x.back() << ")";
	return out;
}

double operator*(const vector<double>& a, const vector<double>& b) {
	double sum = 0;
	for (int i = 0; i < a.size(); ++i)
		sum += a[i] * b[i];
	return sum;
}

struct matrix
{
	int n;
	vector<double> d, l, u;
	vector<int> i, j;

	void init(int n1) {
		n = n1;
		d.clear();
		l.clear();
		u.clear();
		i.clear();
		j.clear();
		d.resize(n);
		i.resize(n+1, 0);
	}

	void toDense(Matrix& m) const {
		m.resize(n, n, 0);
		for (int i = 0; i < n; ++i) {
			m(i, i) = d[i];
			for (int j = 0; j < lineElemCount(i); ++j) {
				m(i, lineElemRow(i, j)) = l[lineElemStart(i) + j];
				m(lineElemRow(i, j), i) = u[lineElemStart(i) + j];
			}
		}
	}

	int lineElemStart(int line) const { 
		return i[line]; 
	}
	int lineStart(int line) const { 
		return j[lineElemStart(line)]; 
	}
	int lineSize(int line) const { 
		return line - lineStart(line); 
	} 
	int lineElemRow(int line, int elem) const { 
		return j[lineElemStart(line) + elem]; 
	}
	int lineElemCount(int line) const { 
		return i[line+1]-i[line]; 
	}
};

struct matrix_iterator
{
	matrix_iterator(const matrix& m, int line) : m(m), line(line) {
		if (m.lineElemCount(line) != 0) {
			pos = m.lineStart(line);
			j = m.lineElemStart(line);
			is_on_elem = true;
			is_empty_line = false;
			line_size = m.lineElemCount(line);
		} else {
			pos = 0;
			j = 0;
			is_on_elem = false;
			is_empty_line = true;
			line_size = 0;
		}
	}

	matrix_iterator& operator++() {
		if (!is_empty_line && j+1 < m.lineElemStart(line) + line_size && pos+1 == m.lineElemRow(line, j+1-m.lineElemStart(line))) {
			is_on_elem = true;
			pos++;
			j++;
		} else {
			is_on_elem = false;
			pos++;
		}
		return *this;
	}

	int getpos(void) const { return pos; }
protected:
	bool is_on_elem, is_empty_line;
	const matrix& m;
	int line, line_size;
	int pos, j, end;
};

struct matrix_iterator_l : public matrix_iterator
{
	matrix_iterator_l(const matrix& m, int line) : matrix_iterator(m, line) {}
	double operator*() const {
		if (is_on_elem) return m.l[j];
		else return 0;
	}
};

struct matrix_iterator_u : public matrix_iterator
{
	matrix_iterator_u(const matrix& m, int row) : matrix_iterator(m, row) {}
	double operator*() const {
		if (is_on_elem) return m.u[j];
		else return 0;
	}
};

struct line_iterator
{
	line_iterator(const vector<pair<int, double>>& line) : line(line), i(0), pos(0), is_on_elem(line.size() != 0) {
		if (line.size() != 0) {
			pos = line[0].first;
			is_on_elem = true;
		}
	}

	line_iterator& operator++() {
		if (line.size() != 0 && pos+1 <= line.back().first && i + 1 < line.size() && pos+1 == line[i+1].first) {
			is_on_elem = true;
			pos++;
			i++;
		} else {
			is_on_elem = false;
			pos++;
		}
		return *this;
	}

	double operator*() const {
		if (is_on_elem) return line[i].second;
		else return 0;
	}

	int getpos(void) const { return pos; }
protected:
	const vector<pair<int, double>>& line;
	int i, pos;
	bool is_on_elem;
};

template<class T1, class T2>
void synchronize_iterators(T1& i, T2& j) {
	while (i.getpos() != j.getpos()) {
		if (i.getpos() < j.getpos())
			++i;
		else
			++j;
	}
}

void lu_decompose(const matrix& a, matrix& lu) {
	lu.init(a.n);
	for (int i = 0; i < a.n; ++i) {
		// Считаем элементы матрицы L
		vector<pair<int, double>> l_add, u_add;
		{
			matrix_iterator_l a_j(a, i);
			int iLineStart = a.lineStart(i);
			for (int j = iLineStart; j < i; ++j, ++a_j) {
				double sum = 0;
				if (j != iLineStart) {
					line_iterator l_k(l_add);
					matrix_iterator_u u_k(lu, j);
					synchronize_iterators(l_k, u_k);
					for (int k = l_k.getpos(); k < j; ++k, ++l_k, ++u_k)
						sum += (*l_k) * (*u_k);
				}

				double res = ((*a_j) - sum) / lu.d[j];
				if (res != 0)
					l_add.push_back({j, res});
			}
		}

		// Считаем элементы матрицы U
		{
			matrix_iterator_u a_j(a, i);
			int iRowStart = a.lineStart(i);
			for (int j = iRowStart; j < i; ++j, ++a_j) {
				double sum = 0;
				if (j != 0) {
					matrix_iterator_l l_k(lu, j);
					line_iterator u_k(u_add);
					synchronize_iterators(l_k, u_k);
					for (int k = l_k.getpos(); k < j; ++k, ++l_k, ++u_k)
						sum += (*l_k) * (*u_k);
				}

				double res = ((*a_j) - sum) / lu.d[j];
				if (res != 0) 
					u_add.push_back({j, res});
			}
		}

		// Находим такие элементы, которые имеются в одном массиве, но нет в другом
		{
			int j = 0;
			while (j < l_add.size() || j < u_add.size()) {
				if (j >= u_add.size()) {
					u_add.insert(u_add.begin() + j, {l_add[j].first, 0});
				} else 
				if (j >= l_add.size()) {
					l_add.insert(l_add.begin() + j, {u_add[j].first, 0});
				} else
				if (l_add[j].first != u_add[j].first) {
					if (u_add[j].first > l_add[j].first)
						u_add.insert(u_add.begin() + j, {l_add[j].first, 0});
					else
						l_add.insert(l_add.begin() + j, {u_add[j].first, 0});
				}
				j++;
			}
		}

		if (l_add.size() != u_add.size())
			throw std::exception();

		// Добавляем формат к обоим массивам
		lu.i[i+1] = lu.i[i] + l_add.size();
		for (int j = 0; j < l_add.size(); ++j) {
			lu.j.push_back(l_add[j].first);
			lu.l.push_back(l_add[j].second);
			lu.u.push_back(u_add[j].second);
		}

		// Считаем элементы диагонали
		{
			double sum = 0;
			if (i != 0) {
				for (int k = 0; k < l_add.size(); ++k)
					sum += l_add[k].second * u_add[k].second;
			}

			double res = sqrt((a.d[i]) - sum);
			lu.d[i] = res;
		}
	}
}

void mul(const matrix& a, vector<double>& x_y) {
	vector<double> result(a.n, 0);

	for (int i = 0; i < a.n; ++i) {
		int start = a.lineElemStart(i);
		int size = a.lineElemCount(i);
		for (int j = 0; j < size; j++) {
			result[i] += a.l[start + j] * x_y[a.lineElemRow(i, j)];
			result[a.lineElemRow(i, j)] += a.u[start + j] * x_y[i];
		}
	}

	// Умножение диагональных элементов на вектор
	for (int i = 0; i < a.n; ++i)
		result[i] += a.d[i] * x_y[i];

	x_y = result;
}

void mul_t(const matrix& a, vector<double>& x_y) {
	vector<double> result(a.n, 0);

	for (int i = 0; i < a.n; ++i) {
		int start = a.lineElemStart(i);
		int size = a.lineElemCount(i);
		for (int j = 0; j < size; j++) {
			result[i] += a.u[start + j] * x_y[a.lineElemRow(i, j)];
			result[a.lineElemRow(i, j)] += a.l[start + j] * x_y[i];
		}
	}

	// Умножение диагональных элементов на вектор
	for (int i = 0; i < a.n; ++i)
		result[i] += a.d[i] * x_y[i];

	x_y = result;
}

void mul_l_invert_t(const matrix& l, vector<double>& y_x) {
	for (int i = l.n - 1; i >= 0; i--) {
		int start = l.lineElemStart(i);
		int size = l.lineElemCount(i);

		y_x[i] /= l.d[i];
		for (int j = 0; j < size; ++j)
			y_x[l.lineElemRow(i, j)] -= y_x[i] * l.l[start + j];
	}
}

void mul_u_invert_t(const matrix& u, vector<double>& y_x) {
	for (int i = 0; i < u.n; ++i) {
		int start = u.lineElemStart(i);
		int size = u.lineElemCount(i);

		sumreal sum = 0;
		for (int j = 0; j < size; ++j)
			sum += u.u[start + j] * y_x[u.lineElemRow(i, j)];
		y_x[i] = (y_x[i] - sum) / u.d[i];
	}
}

void mul_l_invert(const matrix& l, vector<double>& y_x) {
	for (int i = 0; i < l.n; ++i) {
		int start = l.lineElemStart(i);
		int size = l.lineElemCount(i);

		sumreal sum = 0;
		for (int j = 0; j < size; ++j)
			sum += l.l[start + j] * y_x[l.lineElemRow(i, j)];
		y_x[i] = (y_x[i] - sum) / l.d[i];
	}
}

void mul_u_invert(const matrix& u, vector<double>& y_x) {
	for (int i = u.n-1; i >= 0; i--) {
		int start = u.lineElemStart(i);
		int size = u.lineElemCount(i);

		y_x[i] /= u.d[i];
		for (int j = 0; j < size; ++j)
			y_x[u.lineElemRow(i, j)] -= y_x[i] * u.u[start + j];
	}
}

void mul_u(const matrix& u, vector<double>& x_y) {
	vector<double> result(u.n, 0);

	for (int i = 0; i < u.n; ++i) {
		int start = u.lineElemStart(i);
		int size = u.lineElemCount(i);
		for (int j = 0; j < size; j++) {
			result[u.lineElemRow(i, j)] += u.u[start + j] * x_y[i];
		}
	}

	// Умножение диагональных элементов на вектор
	for (int i = 0; i < u.n; ++i)
		result[i] += u.d[i] * x_y[i];

	x_y = result;
}

void mul(const vector<double>& d, vector<double>& x_y) {
	for (int i = 0; i < d.size(); i++)
		x_y[i] /= d[i];
}

void mul_invert(const vector<double>& d, vector<double>& x_y) {
	for (int i = 0; i < d.size(); i++)
		x_y[i] *= d[i];
}

double length(const vector<double>& mas) {
	double sum = 0;
	for (auto& i : mas)
		sum += mas * mas;
	return sqrt(sum);
}

class SLAU
{
public:

void read(string dir) {
	ifstream fin;

	fin.open(dir + "/kuslau.txt");
	fin >> n >> maxiter >> eps;
	fin.close();

	a.n = n;

	a.d.resize(n);
	fin.open(dir + "/di.txt");
	for (auto& i : a.d) fin >> i;
	fin.close();

	f.resize(n);
	fin.open(dir + "/pr.txt");
	for (auto& i : f) fin >> i;
	fin.close();

	a.i.resize(n+1);
	fin.open(dir + "/ig.txt");
	for (auto& i : a.i) { fin >> i; i--; }
	fin.close();

	a.j.resize(a.i.back());
	fin.open(dir + "/jg.txt");
	for (auto& i : a.j) { fin >> i; i--; }
	fin.close();

	a.l.resize(a.i.back());
	fin.open(dir + "/ggl.txt");
	for (auto& i : a.l) { fin >> i; }
	fin.close();

	a.u.resize(a.i.back());
	fin.open(dir + "/ggu.txt");
	for (auto& i : a.u) { fin >> i; }
	fin.close();

	x.resize(n);
	t1.resize(n);
	t2.resize(n);
}

void msg1() {
	/*
	как?
	*/
}

void msg2() {
	lu_decompose(a, lu);

	x.resize(n, 0); // x_0 = (0, 0, ...)

	r = x;
	mul(a, r); // r = A*x_0
	for (int i = 0; i < n; ++i)
		r[i] = f[i]-r[i];
	mul_l_invert(lu, r);
	mul_l_invert_t(lu, r);
	mul_t(a, r);
	mul_u_invert_t(lu, r);

	mul_u(lu, x);

	z = r;
	double rr = r*r;
	double ff = f*f;

	int i = 0;
	while (true) {
		t1 = z;
		mul_u_invert(lu, t1);
		mul(a, t1);
		mul_l_invert(lu, t1);
		mul_l_invert_t(lu, t1);
		mul_t(a, t1);
		mul_u_invert(lu, t1);
		double alpha = (rr) / (t1*z);
		for (int i = 0; i < n; ++i) {
			x[i] += alpha * z[i];
			r[i] -= alpha * t1[i];
		}
		double rr2 = r*r;
		double beta = rr2/rr;
		rr = rr2;
		for (int i = 0; i < n; ++i)
			z[i] = r[i] + beta * z[i];
		double residual = length(r) / length(f);
		i++;

		cout << "Iter: " << i << "\tResidual: " << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}

	mul_u_invert(lu, x);
}

void msg3() {
	lu_decompose(a, lu);

	x.resize(n, 0); // x_0 = (0, 0, ...)

	r = x;
	mul(a, r); // r = A*x_0
	for (int i = 0; i < n; ++i)
		r[i] = f[i]-r[i];
	mul(a.d, r);
	mul(a.d, r);
	mul_t(a, r);
	mul(a.d, r);

	mul_invert(a.d, x);

	z = r;
	double rr = r*r;
	double ff = f*f;

	int i = 0;
	while (true) {
		t1 = z;
		mul(a.d, r);
		mul(a, t1);
		mul(a.d, r);
		mul(a.d, r);
		mul_t(a, t1);
		mul(a.d, r);
		double alpha = (rr) / (t1*z);
		for (int i = 0; i < n; ++i) {
			x[i] += alpha * z[i];
			r[i] -= alpha * t1[i];
		}
		double rr2 = r*r;
		double beta = rr2/rr;
		rr = rr2;
		for (int i = 0; i < n; ++i)
			z[i] = r[i] + beta * z[i];
		double residual = rr/ff;
		i++;

		residual = length(r) / length(f);

		cout << "Iter: " << i << "\tResidual: " << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}

	mul_u_invert(lu, x);
}

void los1() {
	x.resize(n, 0); // x_0 = (0, 0, ...)

	r = x;
	mul(a, r); // r = A*x_0
	for (int i = 0; i < n; i++) // r = f - A*x_0
		r[i] = f[i] - r[i]; 

	z = r;

	p = z;
	mul(a, p); // p = A z

	double residual = r*r;

	int i = 0;
	while (true) {
		double pp = p*p;
		double alpha = (p*r) / pp;
		for (int i = 0; i < n; ++i) {
			x[i] += alpha * z[i];
			r[i] -= alpha * p[i];
		}
		t1 = r;
		mul(a, t1);
		double beta = -(p*t1) / pp;
		for (int i = 0; i < n; ++i) {
			z[i] = r[i] + beta * z[i];
			p[i] = t1[i] + beta * p[i];
		}
		residual -= alpha * alpha * pp;
		i++;

		residual = length(r) / length(f);

		cout << "Iter: " << i << "\tResidual: " << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}
}

void los2() {
	lu_decompose(a, lu);
	x.resize(n, 0); // x_0 = (0, 0, ...)

	r = x;
	mul(a, r); // r = A*x_0
	for (int i = 0; i < n; i++) // r = f - A*x_0
		r[i] = f[i] - r[i];
	mul_l_invert(lu, r); // r = L^-1(f - A*x_0)

	z = r;
	mul_u_invert(lu, z); // z = U^-1 r

	p = z;
	mul(a, p); // p = A z
	mul_l_invert(lu, p); // p = L^-1 A z

	double residual = r*r;

	int i = 0;
	while (true) {
		double pp = p*p;
		double alpha = (p*r) / pp;
		for (int i = 0; i < n; ++i) {
			x[i] += alpha * z[i];
			r[i] -= alpha * p[i];
		}
		t1 = r;
		mul_u_invert(lu, t1);
		t2 = t1;
		mul(a, t2);
		mul_l_invert(lu, t2);
		double beta = -(p*t2) / pp;
		for (int i = 0; i < n; ++i) {
			z[i] = t1[i] + beta * z[i];
			p[i] = t2[i] + beta * p[i];
		}
		residual -= alpha * alpha * pp;
		i++;

		residual = length(r) / length(f);

		cout << "Iter: " << setw(4) << i << "\tResidual: " << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}
}

void los3() {
	x.resize(n, 0); // x_0 = (0, 0, ...)

	r = x;
	mul(a, r); // r = A*x_0
	for (int i = 0; i < n; i++) // r = f - A*x_0
		r[i] = f[i] - r[i];
	mul(a.d, r); // r = L^-1(f - A*x_0)

	z = r;
	mul(a.d, z); // z = U^-1 r

	p = z;
	mul(a, p); // p = A z
	mul(a.d, p); // p = L^-1 A z

	double residual = r*r;

	int i = 0;
	while (true) {
		double pp = p*p;
		double alpha = (p*r) / pp;
		for (int i = 0; i < n; ++i) {
			x[i] += alpha * z[i];
			r[i] -= alpha * p[i];
		}
		t1 = r;
		mul(a.d, t1);
		t2 = t1;
		mul(a, t2);
		mul(a.d, t2);
		double beta = -(p*t2) / pp;
		for (int i = 0; i < n; ++i) {
			z[i] = t1[i] + beta * z[i];
			p[i] = t2[i] + beta * p[i];
		}
		residual -= alpha * alpha * pp;
		i++;

		residual = length(r) / length(f);

		cout << "Iter: " << i << "\tResidual: " << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}
}

int n, maxiter;
double eps;

matrix a, lu;

vector<double> f;

vector<double> r, z, p;

vector<double> x, t1, t2;

};

int main() {
	string dir;
	cout << "Enter dir: ";
	//cin >> dir;

	SLAU s;
	s.read(dir);

	ofstream fout(dir + "/matrix.txt");
	vector<double> f1 = s.f;
	Matrix m, l, u, a, sub;
	Vector pr(s.f.size()), pr1, f1_m(f1.size());
	for (int i = 0; i < s.f.size(); i++) pr(i) = f1[i];
	s.a.toDense(m);

	lu_decompose(s.a, s.lu);
	matrix l_s = s.lu, u_s = s.lu;
	l_s.u.clear(); l_s.u.resize(l_s.l.size(), 0);
	u_s.l.clear(); u_s.l.resize(u_s.u.size(), 0);
	l_s.toDense(l);
	u_s.toDense(u);
	mul(l, u, a);
	a.negate();
	sum(m, a, sub);

	//mul(m, pr, pr1);
	//mul(s.a, f1);

	//transpose(m);
	//mul(m, pr, pr1);
	//mul_t(s.a, f1);

	//mul(u, pr, pr1);
	//mul_u(s.lu, f1);

	/*pr(3) = 115;
	mul(l, pr, pr1);
	for (int i = 0; i < pr1.size(); i++) f1[i] = pr1(i);
	pr1 = pr;
	mul_l_invert(s.lu, f1);*/

	/*pr(3) = 115;
	transpose(u);
	mul(u, pr, pr1);
	for (int i = 0; i < pr1.size(); i++) f1[i] = pr1(i);
	pr1 = pr;
	mul_u_invert_t(s.lu, f1);*/

	/*pr(3) = 115;
	mul(u, pr, pr1);
	for (int i = 0; i < pr1.size(); i++) f1[i] = pr1(i);
	pr1 = pr;
	mul_u_invert(s.lu, f1);*/

	/*pr(3) = 115;
	transpose(l);
	mul(l, pr, pr1);
	for (int i = 0; i < pr1.size(); i++) f1[i] = pr1(i);
	pr1 = pr;
	mul_l_invert_t(s.lu, f1);*/

	pr1 = pr;
	for (int i = 0; i < s.f.size(); i++) pr1(i) = s.f[i];
	//mul(s.a, s.f);
	s.los2();
	f1 = s.x;

	for (int i = 0; i < f1.size(); i++) f1_m(i) = f1[i];

	for (int i = 0; i < sub.height(); i++) {
		for (size_t j = 0; j < sub.width(); j++) {
			if (fabs(sub(i, j)) < 0.000001)
				sub(i, j) = 0;
		}
	}

	fout << "m:" << endl;
	m.save(fout);
	fout << endl;

	fout << "l:" << endl;
	l.save(fout);
	fout << endl;

	fout << "u:" << endl;
	u.save(fout);
	fout << endl;

	fout << "sub:" << endl;
	sub.save(fout);
	fout << endl;

	fout << "pr1_m:" << endl;
	pr1.save(fout);
	fout << endl;

	fout << "pr1:" << endl;
	f1_m.save(fout);
	fout << endl;

	fout.close();


	/*chrono::high_resolution_clock::time_point t1, t2;

	cout << "LOS 1:" << endl;
	t1 = chrono::high_resolution_clock::now();
	s.los1();
	t2 = chrono::high_resolution_clock::now();
	cout << "Time: " << chrono::high_resolution_clock::duration_cast<chrono::high_resolution_clock::seconds>(t2-t1).count() << endl;
	cout << s.x;

	cout << "LOS 2:" << endl;
	t1 = now();
	s.los2();
	t2 = now();
	cout << "Time: " << duration_cast<seconds>(t2-t1).count() << endl;
	cout << s.x;

	cout << "LOS 3:" << endl;
	t1 = now();
	s.los3();
	t2 = now();
	cout << "Time: " << duration_cast<seconds>(t2-t1).count() << endl;
	cout << s.x;*/

	//system("pause");
}