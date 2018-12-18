#include "objects.h"

//-----------------------------------------------------------------------------
pair<fnm_t, jnm_t> two_circles(circle a, circle b) {
	 fn1_t f1 = [=a] (const xn_t& x) -> double {
	 	#ifdef _DEBUG
		if (x.size() != 2)
			throw exception();
		#endif
	 	return (x[0]-a.x)*(x[0]-a.x) + (x[1]-a.y)*(x[1]-a.y) - a.r*a.r;
	 };

	 fn1_t f2 = [=b] (const xn_t& x) -> double {
	 	#ifdef _DEBUG
		if (x.size() != 2)
			throw exception();
		#endif
	 	return (x[0]-b.x)*(x[0]-b.x) + (x[1]-b.y)*(x[1]-b.y) - b.r*b.r;
	 };

	 fnm_t f = {f1, f2};

	 jnm_t j = [=a, =b] (const xn_t& x) -> matrix {
	 	matrix result(2, xn_t(2));
	 	result[0][0] = 2*(x[0]-a.x);
	 	result[0][1] = 2*(x[1]-a.y);
	 	result[1][0] = 2*(x[0]-b.x);
	 	result[1][1] = 2*(x[1]-b.y);
	 };
}

//-----------------------------------------------------------------------------
pair<fnm_t, jnm_t> two_circles_and_line(circle a, circle b, line l) {

}

//-----------------------------------------------------------------------------
pair<fnm_t, jnm_t> three_lines(line a, line b, line c) {
	// (y-y_1)/(y_2-y_1) = (x-x_1)/(x_2-x_1)

	fn1_t f1 = [=a] (const xn_t& x) -> double {
	 	#ifdef _DEBUG
		if (x.size() != 2)
			throw exception();
		#endif
	 	return (x[0]-a.x)*(x[0]-a.x) + (x[1]-a.y)*(x[1]-a.y) - a.r*a.r;
	 };

	 fn1_t f2 = [=b] (const xn_t& x) -> double {
	 	#ifdef _DEBUG
		if (x.size() != 2)
			throw exception();
		#endif
	 	return (x[0]-b.x)*(x[0]-b.x) + (x[1]-b.y)*(x[1]-b.y) - b.r*b.r;
	 };

	 fnm_t f = {f1, f2};

	 jnm_t j = [=a, =b] (const xn_t& x) -> matrix {
	 	matrix result(2, xn_t(2));
	 	result[0][0] = 2*(x[0]-a.x);
	 	result[0][1] = 2*(x[1]-a.y);
	 	result[1][0] = 2*(x[0]-b.x);
	 	result[1][1] = 2*(x[1]-b.y);
	 };
}

//-----------------------------------------------------------------------------
pair<fnm_t, jnm_t> sin_and_line(point b) {

}