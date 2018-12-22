#include "find_borders.h"

//-----------------------------------------------------------------------------
FindBorders::FindBorders(int maxSize, int border, bool isRevert) : m_maxSize(maxSize), m_border(border), m_isRevert(isRevert) {
	double inf = std::numeric_limits<double>::infinity();
	min = vec2(inf, inf);
	max = vec2(-inf, -inf);
}

//-----------------------------------------------------------------------------
void FindBorders::process(vec2 point) {
	if (point.x > max.x) max.x = point.x;
	if (point.y > max.y) max.y = point.y;
	if (point.x < min.x) min.x = point.x;
	if (point.y < min.y) min.y = point.y;
}

//-----------------------------------------------------------------------------
void FindBorders::process(circle cr) {
	vec2 c(cr.c.x, cr.c.y);
	vec2 sz(cr.r, cr.r);
	process(c + sz);
	process(c - sz);
}

//-----------------------------------------------------------------------------
void FindBorders::process(const std::vector<vec2>& points) {
	for (const auto& i : points)
		process(i);
}

//-----------------------------------------------------------------------------
void FindBorders::finish(void) {
	if (max.x - min.x > max.y - min.y) {
		m_size.x = m_maxSize;
		m_size.y = (max.y - min.y)/(max.x - min.x) * m_size.x;

		m_scale = m_maxSize/(max.x - min.x);
	} else {
		m_size.y = m_maxSize;
		m_size.x = (max.x - min.x)/(max.y - min.y) * m_size.y;

		m_scale = m_maxSize/(max.y - min.y);
	}
	m_offset = -min;
	m_size += vec2(m_border) * 2;
}

//-----------------------------------------------------------------------------
vec2 FindBorders::toImg(vec2 point) const {
	vec2 result = (point + m_offset) * m_scale + vec2(m_border, m_border);
	if (m_isRevert)
		result.y = m_size.y - result.y - 1;
	return result;
}

//-----------------------------------------------------------------------------
double FindBorders::toImg(double a) const {
	return m_scale * a;
}

//-----------------------------------------------------------------------------
vec2 FindBorders::fromImg(vec2 point) const {
	if (m_isRevert)
		point.y = m_size.y - point.y - 1;
	return (point - vec2(m_border, m_border)) / m_scale - m_offset;
}

//-----------------------------------------------------------------------------
vec2 FindBorders::getCalculatedSize(void) const {
	return m_size;
}

//-----------------------------------------------------------------------------
bool FindBorders::isInside(vec2 pos) const {
	return (pos.x > min.x) && (pos.y > min.y) && (pos.x < max.x) && (pos.y < max.y);
}