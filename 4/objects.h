#pragma once

#include "logic.h"

//-----------------------------------------------------------------------------
struct point 
{ 
	point(double x, double y) : x(x), y(y) {} 
	point(const xn_t& x) : x(x[0]), y(x[1]) {
		myassert(x.size() == 2);
	}
	double x, y; 
};

struct circle {
	circle(double x, double y, double r) : c(x, y), r(r) {}
	circle(point c, double r) : c(c), r(r) {}
	circle(const xn_t& x) : c(x[0], x[1]), r(x[2]) {
		myassert(x.size() == 3);
	}
	point c; double r;
};

struct line { point a, b; };

//-----------------------------------------------------------------------------
pair<sle_f, fnmv_f> one_circle(circle a);
pair<sle_f, fnmv_f> two_circles(circle a, circle b);
pair<sle_f, fnmv_f> two_circles_and_line(circle a, circle b, line l);
pair<sle_f, fnmv_f> three_lines(line a, line b, line c);

pair<sle_f, fnmv_f> sin_and_line(point b); // a is {0, 0}

pair<sle_f, fnmv_f> three_circles(
	circle a, circle b, circle c,
	bool a_in, bool b_in, bool c_in
);

//-----------------------------------------------------------------------------
inline double sqr(double a) {return a*a;}
inline double sign(double a) { if (a == 0) return 0; else return (a<0)?-1:1; }

double dist(const point& a, const point& b);

double line_f(const line& l, const point& x);
double line_d_x(const line& l, const point& x);
double line_d_y(const line& l, const point& x);

double circle_f(const circle& a, const point& x);
double circle_d_x(const circle& a, const point& x);
double circle_d_y(const circle& a, const point& x);

double circles_f(const circle& a, const circle& x, bool in);
double circles_d_x(const circle& a, const circle& x, bool in);
double circles_d_y(const circle& a, const circle& x, bool in);
double circles_d_r(const circle& a, const circle& x, bool in);