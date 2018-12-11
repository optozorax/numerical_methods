#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;

ostream& operator<<(ostream& out, const vector<double>& x) {
	out << "(";
	for (int i = 0; i < x.size()-1; ++i)
		out << x[i] << ", ";
	out << x.back() << ")";
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

void lu_decompose(const matrix& a, matrix& to) {
	for (int i = 0; i < a.n; ++i) {
		double lsum = 0, usum = 0;
		// L, D, U
	}

	for (int i = 0; i < n; i++) {
		double sum_di = 0;
		int i0 = ia[i];
		int i1 = ia[i+1];
		for (int k = i0; k < i1; k++) {
			int j = ja[k];
			int j0 = ia[j];
			int j1 = ia[j + 1];
			double sum_u = 0, sum_l = 0;
			int ki = i0;
			int kj = j0;
			while (ki < k && kj < j1) {
				int j_ki = ja[ki];
				int j_kj = ja[kj];
				if (j_ki == j_kj) {
					sum_l += b_al[ki] * b_au[kj];
					sum_u += b_au[ki] * b_al[kj];
					ki++;
					kj++;
				} else if (j_ki > j_kj) 
					kj++;
				else 
					ki++;
			}
			b_al[k] = a_al[k] - sum_l;
			b_au[k] = (a_au[k] - sum_u) / b_di[j];
			sum_di += b_al[k] * b_au[k];
		}
		b_di[i] = a_di[i] - sum_di;
	}
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

	fin.open(dir + "kuslau.txt");
	fin >> n >> maxiter >> eps;
	fin.close();

	a.n = n;

	a.d.resize(n);
	fin.open(dir + "di.txt");
	for (auto& i : a.d) fin >> i;
	fin.close();

	f.resize(n);
	fin.open(dir + "pr.txt");
	for (auto& i : f) fin >> i;
	fin.close();

	a.i.resize(n+1);
	fin.open(dir + "ig.txt");
	for (auto& i : a.i) { fin >> i; i--; }
	fin.close();

	a.j.resize(a.i.back());
	fin.open(dir + "jg.txt");
	for (auto& i : a.j) { fin >> i; i--; }
	fin.close();

	a.l.resize(a.i.back());
	fin.open(dir + "ggl.txt");
	for (auto& i : a.l) { fin >> i; i--; }
	fin.close();

	a.u.resize(a.i.back());
	fin.open(dir + "ggu.txt");
	for (auto& i : a.u) { fin >> i; i--; }
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
		residual = rr/ff;
		i++;

		cout << "Iter: " << i << "\tResidual: " << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}
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
		residual = rr/ff;
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
	r = f - r; // r = f - A*x_0

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
	lu_decompose(a, lu);
	x.resize(n, 0); // x_0 = (0, 0, ...)

	r = x;
	mul(a, r); // r = A*x_0
	r = f - r; // r = f - A*x_0
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

		cout << "Iter: " << i << "\tResidual: " << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}
}

void los3() {
	x.resize(n, 0); // x_0 = (0, 0, ...)

	r = x;
	mul(a, r); // r = A*x_0
	r = f - r; // r = f - A*x_0
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

matrix a, lu;

vector<double> f;

vector<double> r, z, p;

vector<double> x, t1, t2;

};

int main() {
	using namespace chrono;
	using namespace chrono::high_resolution_clock;

	string dir;
	cout << "Enter dir: ";
	cin >> dir;

	SLAU s;
	s.read(dir);

	time_point t1, t2;

	cout << "LOS 1:" << endl;
	t1 = now();
	s.los1();
	t2 = now();
	cout << "Time: " << duration_cast<seconds>(t2-t1).count() << endl;
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
	cout << s.x;

	system("pause");
}