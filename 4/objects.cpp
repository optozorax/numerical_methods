#include <cassert>

#include "objects.h"

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
pair<sle_f, fnm_f> two_circles(circle a, circle b) {
	//-------------------------------------------------------------------------
	fn_f f1 = [a] (const xn_t& x) -> double { return circle_f(a, x); };
	fn_f f2 = [b] (const xn_t& x) -> double { return circle_f(b, x); };

	fnm_f f = {f1, f2};

	//-------------------------------------------------------------------------
	jnm_f j = [a, b] (const xn_t& x) -> matrix_t {
		matrix_t result(2, xn_t(2));
		result[0][0] = circle_d_x(a, x);
		result[0][1] = circle_d_y(a, x);

		result[1][0] = circle_d_x(b, x);
		result[1][1] = circle_d_y(b, x);

		return result;
	};

	return {get_sle_function(j, f), f};
}

//-----------------------------------------------------------------------------
pair<sle_f, fnm_f> two_circles_and_line(circle a, circle b, line l) {
	//-------------------------------------------------------------------------
	fn_f f1 = [a] (const xn_t& x) -> double { return circle_f(a, x); };
	fn_f f2 = [b] (const xn_t& x) -> double { return circle_f(b, x); };
	fn_f f3 = [l] (const xn_t& x) -> double { return line_f(l, x); };

	fnm_f f = {f1, f2, f3};

	//-------------------------------------------------------------------------
	jnm_f j = [a, b, l] (const xn_t& x) -> matrix_t {
		matrix_t result(3, xn_t(2));
		result[0][0] = circle_d_x(a, x);
		result[0][1] = circle_d_y(a, x);

		result[1][0] = circle_d_x(b, x);
		result[1][1] = circle_d_y(b, x);

		result[2][0] = line_d_x(l, x);
		result[2][1] = line_d_y(l, x);

		return result;
	};

	return {get_sle_function(j, f), f};
}

//-----------------------------------------------------------------------------
pair<sle_f, fnm_f> three_lines(line a, line b, line c) {
	//-------------------------------------------------------------------------
	fn_f f1 = [a] (const xn_t& x) -> double { return line_f(a, x); };
	fn_f f2 = [b] (const xn_t& x) -> double { return line_f(b, x); };
	fn_f f3 = [c] (const xn_t& x) -> double { return line_f(c, x); };

	fnm_f f = {f1, f2, f3};

	//-------------------------------------------------------------------------
	jnm_f j = [a, b, c] (const xn_t& x) -> matrix_t {
		matrix_t result(3, xn_t(2));
		result[0][0] = line_d_x(a, x);
		result[0][1] = line_d_y(a, x);

		result[1][0] = line_d_x(b, x);
		result[1][1] = line_d_y(b, x);

		result[2][0] = line_d_x(c, x);
		result[2][1] = line_d_y(c, x);

		return result;
	};

	return {get_sle_function(j, f), f};
}

//-----------------------------------------------------------------------------
pair<sle_f, fnm_f> sin_and_line(point b) {
	//-------------------------------------------------------------------------
	line l = {point(0, 0), b};
	fn_f f1 = [l] (const xn_t& x) -> double {return line_f(l, x); };
	fn_f f2 = [] (const xn_t& x) -> double {
		assert(x.size() == 2);
		return x[1]-sin(x[0]);
	};

	fnm_f f = {f1, f2};

	//-------------------------------------------------------------------------
	jnm_f j = [l] (const xn_t& x) -> matrix_t {
		assert(x.size() == 2);
		matrix_t result(2, xn_t(2));
		result[0][0] = line_d_x(l, x);
		result[0][1] = line_d_y(l, x);

		result[1][0] = -cos(x[0]);
		result[1][1] = 1;

		return result;
	};

	return {get_sle_function(j, f), f};
}

//-----------------------------------------------------------------------------
pair<sle_f, fnm_f> three_circles(circle a, circle b, circle c, bool a_in, bool b_in, bool c_in) {
	//-------------------------------------------------------------------------
	fn_f f1 = [a, a_in] (const xn_t& x) -> double { return circles_f(a, x, a_in); };
	fn_f f2 = [b, b_in] (const xn_t& x) -> double { return circles_f(b, x, b_in); };
	fn_f f3 = [c, c_in] (const xn_t& x) -> double { return circles_f(c, x, c_in); };

	fnm_f f = {f1, f2, f3};

	//-------------------------------------------------------------------------
	jnm_f j = [a, b, c, a_in, b_in, c_in] (const xn_t& x) -> matrix_t {
		matrix_t result(3, xn_t(3));
		result[0][0] = circles_d_x(a, x, a_in);
		result[0][1] = circles_d_y(a, x, a_in);
		result[0][2] = circles_d_r(a, x, a_in);

		result[1][0] = circles_d_x(b, x, b_in);
		result[1][1] = circles_d_y(b, x, b_in);
		result[1][2] = circles_d_r(b, x, b_in);

		result[2][0] = circles_d_x(c, x, c_in);
		result[2][1] = circles_d_y(c, x, c_in);
		result[2][2] = circles_d_r(c, x, c_in);

		return result;
	};

	return {get_sle_function(j, f), f};
}

//-----------------------------------------------------------------------------
// sle_f three_circles_2(circle a, circle b, circle c, bool a_in, bool b_in, bool c_in) {
// }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
double dist(const point& a, const point& b) {
	return sqrt(sqr(a.x-b.x) + sqr(a.y-b.y));
}

//-----------------------------------------------------------------------------
double line_f(const line& l, const point& x) {
	const double relative = 1e-10;
	double length = dist(l.a, l.b);

	bool first_zero = fabs(l.b.x-l.a.x)/length < relative;
	bool second_zero = fabs(l.b.y-l.a.y)/length < relative;

	assert(!(first_zero && second_zero));

	if (first_zero) 
		return x.x-l.a.x;
	else 
	if (second_zero) 
		return x.y-l.a.y;
	else 
		return (x.x-l.a.x)/(l.b.x-l.a.x) - (x.y-l.a.y)/(l.b.y-l.a.y);
}

//-----------------------------------------------------------------------------
double line_d_x(const line& l, const point& x) {
	const double relative = 1e-10;
	double length = dist(l.a, l.b);

	bool first_zero = fabs(l.b.x-l.a.x)/length < relative;
	bool second_zero = fabs(l.b.y-l.a.y)/length < relative;

	assert(!(first_zero && second_zero));

	if (first_zero) 
		return 1;
	else 
	if (second_zero) 
		return 0;
	else 
		return 1.0/(l.b.x-l.a.x);
}

//-----------------------------------------------------------------------------
double line_d_y(const line& l, const point& x) {
	const double relative = 1e-10;
	double length = dist(l.a, l.b);

	bool first_zero = fabs(l.b.x-l.a.x)/length < relative;
	bool second_zero = fabs(l.b.y-l.a.y)/length < relative;

	assert(!(first_zero && second_zero));

	if (first_zero) 
		return 0;
	else 
	if (second_zero) 
		return 1;
	else 
		return 1.0/(l.b.y-l.a.y);
}

//-----------------------------------------------------------------------------
double circle_f(const circle& a, const point& x) {
	return sqr(a.c.x-x.x) + sqr(a.c.y-x.y) - sqr(a.r);
}

//-----------------------------------------------------------------------------
double circle_d_x(const circle& a, const point& x) {
	return -2*(a.c.x-x.x);
}

//-----------------------------------------------------------------------------
double circle_d_y(const circle& a, const point& x) {
	return -2*(a.c.y-x.y);
}

//-----------------------------------------------------------------------------
double circles_f(const circle& a, const circle& x, bool in) {
	double result = sqr(a.c.x-x.c.x) + sqr(a.c.y-x.c.y);
	if (in)
		return result - sqr(a.r-fabs(x.r));
	else
		return result - sqr(a.r+fabs(x.r));
}

//-----------------------------------------------------------------------------
double circles_d_x(const circle& a, const circle& x, bool in) {
	return -2*(a.c.x-x.c.x);
}

//-----------------------------------------------------------------------------
double circles_d_y(const circle& a, const circle& x, bool in) {
	return -2*(a.c.y-x.c.y);
}

//-----------------------------------------------------------------------------
double circles_d_r(const circle& a, const circle& x, bool in) {
	if (in)
		return -2*(x.r-a.r*sign(x.r));
	else
		return -2*(x.r+a.r*sign(x.r));
}