#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>

#include "../1/matrix.h"

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

	/*void addElem(int i1, int j1, double elem) {
		if (i[i1]-i[i1-1] == 0) {
			i[i1] = j.size();
			i[i1+1] = j.size();
		}
		j.push_back(j1);
		i[i1+1]++;
		if (j1 < i1)
			l.push_back(elem);
		else 
			u.push_back(elem);
	}*/

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
		if (!is_empty_line && j+1 < line_size && pos+1 == m.lineElemRow(line, j+1)) {
			is_on_elem = true;
			pos++;
			j++;
		} else {
			is_on_elem = false;
			pos++;
		}
		return *this;
	}

	int getRow(void) const { return pos; }
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
	line_iterator(const vector<pair<int, double>>& line, int pos) : line(line), i(0), pos(pos), is_on_elem(line.size() != 0) {}

	line_iterator& operator++() {
		if (line.size() != 0 && pos <= line.back().first && pos+1 == line[i+1].first) {
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
protected:
	const vector<pair<int, double>>& line;
	int i, pos;
	bool is_on_elem;
};

void mul_l_invert_t(const matrix& l, vector<double>& y_x) {

}

void mul_u_invert_t(const matrix& u, vector<double>& y_x) {

}

void mul_l_invert(const matrix& l, vector<double>& y_x) {

}

void mul_u_invert(const matrix& u, vector<double>& y_x) {

}

void mul_u(const matrix& u, vector<double>& y_x) {

}

void lu_decompose(const matrix& a, matrix& l, matrix& u) {
	l.init(a.n);
	u.init(a.n);
	for (int i = 0; i < a.n; ++i) {
		// Считаем элементы матрицы L
		vector<pair<int, double>> l_add, u_add;
		{
			matrix_iterator_l a_j(a, i);
			int iLineStart = a.lineStart(i);
			for (int j = iLineStart; j < i; ++j, ++a_j) {
				double sum = 0;
				if (j != 0) {
					line_iterator l_k(l_add, a_j.getRow());
					matrix_iterator_u u_k(u, j);
					for (int k = 0; k < j; ++k, ++l_k, ++u_k)
						sum += (*l_k) * (*u_k);
				}

				double res = ((*a_j) - sum) / l.d[j];
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
					matrix_iterator_l l_k(l, j);
					line_iterator u_k(u_add, a_j.getRow());
					for (int k = 0; k < j; ++k, ++l_k, ++u_k)
						sum += (*l_k) * (*u_k);
				}

				double res = ((*a_j) - sum) / u.d[j];
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
		l.i[i+1] = l.i[i] + l_add.size();
		for (int j = 0; j < l_add.size(); ++j) {
			l.j.push_back(l_add[j].first);
			l.l.push_back(l_add[j].second);
		}

		u.i[i+1] = u.i[i] + u_add.size();
		for (int j = 0; j < u_add.size(); ++j) {
			u.j.push_back(u_add[j].first);
			u.u.push_back(u_add[j].second);
		}

		// Считаем элементы диагонали
		{
			double sum = 0;
			if (i != 0) {
				matrix_iterator_l l_k(l, i);
				matrix_iterator_u u_k(u, i);
				for (int k = 0; k < i; ++k, ++l_k, ++u_k)
					sum += (*l_k) * (*u_k);
			}

			double res = sqrt((a.d[i]) - sum);
			l.d[i] = res;
			u.d[i] = res;
		}
	}

	l.u.resize(l.l.size(), 0);
	u.l.resize(u.u.size(), 0);
}

void mul(const matrix& a, vector<double>& x_y) {

}

void mul_t(const matrix& a, vector<double>& x_y) {

}

void mul(const vector<double>& d, vector<double>& x_y) {

}

void mul_invert(const vector<double>& d, vector<double>& x_y) {

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
	lu_decompose(a, l, u);

	x.resize(n, 0); // x_0 = (0, 0, ...)

	r = x;
	mul(a, r); // r = A*x_0
	for (int i = 0; i < n; ++i)
		r[i] = f[i]-r[i];
	mul_l_invert(l, r);
	mul_l_invert_t(l, r);
	mul_t(a, r);
	mul_u_invert_t(u, r);

	mul_u(u, x);

	z = r;
	double rr = r*r;
	double ff = f*f;

	int i = 0;
	while (true) {
		t1 = z;
		mul_u_invert(u, t1);
		mul(a, t1);
		mul_l_invert(l, t1);
		mul_l_invert_t(l, t1);
		mul_t(a, t1);
		mul_u_invert(u, t1);
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

		cout << "Iter: " << i << "\tResidual: " << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}
}

void msg3() {
	lu_decompose(a, l, u);

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

		cout << "Iter: " << i << "\tResidual: " << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}
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

		cout << "Iter: " << i << "\tResidual: " << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}
}

void los2() {
	lu_decompose(a, l, u);
	x.resize(n, 0); // x_0 = (0, 0, ...)

	r = x;
	mul(a, r); // r = A*x_0
	for (int i = 0; i < n; i++) // r = f - A*x_0
		r[i] = f[i] - r[i];
	mul_l_invert(l, r); // r = L^-1(f - A*x_0)

	z = r;
	mul_u_invert(u, z); // z = U^-1 r

	p = z;
	mul(a, p); // p = A z
	mul_l_invert(l, p); // p = L^-1 A z

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
		mul_u_invert(u, t1);
		t2 = t1;
		mul(a, t2);
		mul_l_invert(l, t2);
		double beta = -(p*t2) / pp;
		for (int i = 0; i < n; ++i) {
			z[i] = t1[i] + beta * z[i];
			p[i] = t2[i] + beta * p[i];
		}
		residual -= alpha * alpha * pp;
		i++;

		cout << "Iter: " << i << "\tResidual: " << residual << endl;
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

		cout << "Iter: " << i << "\tResidual: " << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}
}

int n, maxiter;
double eps;

matrix a, l, u;

vector<double> f;

vector<double> r, z, p;

vector<double> x, t1, t2;

};

int main() {
	string dir;
	cout << "Enter dir: ";
	//cin >> dir;
	dir = "test1";

	SLAU s;
	s.read(dir);

	ofstream fout(dir + "/matrix.txt");
	Matrix m, l, u, a, sub;
	s.a.toDense(m);

	lu_decompose(s.a, s.l, s.u);
	s.l.toDense(l);
	s.u.toDense(u);
	mul(l, u, a);
	a.negate();
	sum(m, a, sub);

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

	system("pause");
}