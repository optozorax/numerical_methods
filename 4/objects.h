#pragma once

#include "main.h"

//-----------------------------------------------------------------------------
struct point { double x, y; };
struct circle { point c; double r; };
struct line { point a, b; };

point to(const xn_t& x);
circle to(const xn_t& x);

//-----------------------------------------------------------------------------
pair<fnm_t, jnm_t> two_circles(circle a, circle b);
pair<fnm_t, jnm_t> two_circles_and_line(circle a, circle b, line l);
pair<fnm_t, jnm_t> three_lines(line a, line b, line c);

pair<fnm_t, jnm_t> sin_and_line(point b); // a is {0, 0}

pair<fnm_t, jnm_t> three_circles(
	circle a, circle b, circle c,
	bool a_in, bool b_in, bool c_in
);

// pair<fnm_t, jnm_t> three_circles_2(
// 	circle a, circle b, circle c, 
// 	bool a_in, bool b_in, bool c_in
// );

//-----------------------------------------------------------------------------
inline double sqr(double a) {return a*a;}

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