#pragma once

//-------------------------------------------------------------------------
class vec2
{
public:
	double x, y;

	vec2();
	explicit vec2(double a);
	vec2(double x, double y);

	vec2& negate(void);
	vec2& normalize(void);

	double length(void) const;
	double lengthSqr(void) const;

	//---------------------------------------------------------------------
	vec2& operator+=(const vec2& a);
	vec2& operator-=(const vec2& a);
	vec2& operator*=(double a);
	vec2& operator/=(double a);
};

//-------------------------------------------------------------------------
vec2 operator-(const vec2& a);
vec2 operator+(const vec2& a, const vec2& b);
vec2 operator-(const vec2& a, const vec2& b);
vec2 operator*(const vec2& a, double k);
vec2 operator/(const vec2& a, double k);
vec2 operator*(double k, const vec2& a);
vec2 operator/(double k, const vec2& a);

//=============================================================================
//=============================================================================
//=============================================================================

//-----------------------------------------------------------------------------
inline vec2::vec2() : vec2(0) {}

//-----------------------------------------------------------------------------
inline vec2::vec2(double a) : x(a), y(a) {}

//-----------------------------------------------------------------------------
inline vec2::vec2(double x, double y) : x(x), y(y) {}

//-----------------------------------------------------------------------------
inline vec2& vec2::negate(void) {
	x = -x;
	y = -y;
	return *this;
}

//-----------------------------------------------------------------------------
inline vec2& vec2::normalize(void) {
	double length1 = length();
	if (length1 != 0) {
		x /= length1;
		y /= length1;
	}
	return *this;
}

//-----------------------------------------------------------------------------
inline double vec2::length(void) const {
	return std::sqrt(x*x + y*y);
}

//-----------------------------------------------------------------------------
inline double vec2::lengthSqr(void) const {
	return x*x + y*y;
}

//-----------------------------------------------------------------------------
inline vec2& vec2::operator+=(const vec2& a) {
	x += a.x;
	y += a.y;
	return *this;
}

//-----------------------------------------------------------------------------
inline vec2& vec2::operator-=(const vec2& a) {
	x -= a.x;
	y -= a.y;
	return *this;
}

//-----------------------------------------------------------------------------
inline vec2& vec2::operator*=(double a) {
	x *= a;
	y *= a;
	return *this;
}

//-----------------------------------------------------------------------------
inline vec2& vec2::operator/=(double a) {
	#ifdef _DEBUG
	if (a == 0)
		throw std::exception("Divide by zero");
	#endif
	x /= a;
	y /= a;
	return *this;
}

//-----------------------------------------------------------------------------
inline vec2 operator-(const vec2& a) {
	return vec2(-a.x, -a.y);
}

//-----------------------------------------------------------------------------
inline vec2 operator+(const vec2& a, const vec2& b) {
	return vec2(a.x + b.x, a.y + b.y);
}

//-----------------------------------------------------------------------------
inline vec2 operator-(const vec2& a, const vec2& b) {
	return vec2(a.x - b.x, a.y - b.y);
}

//-----------------------------------------------------------------------------
inline vec2 operator*(const vec2& a, double k) {
	return vec2(a.x * k, a.y * k);
}

//-----------------------------------------------------------------------------
inline vec2 operator/(const vec2& a, double k) {
	#ifdef _DEBUG
	if (k == 0)
		throw std::exception("Divide by zero");
	#endif
	return vec2(a.x / k, a.y / k);
}

//-----------------------------------------------------------------------------
inline vec2 operator*(double k, const vec2& a) {
	return vec2(a.x * k, a.y * k);
}

//-----------------------------------------------------------------------------
inline vec2 operator/(double k, const vec2& a) {
	#ifdef _DEBUG
	if (a.x == 0 || a.y == 0)
		throw std::exception("Divide by zero");
	#endif
	return vec2(k / a.x, k/ a.y);
}