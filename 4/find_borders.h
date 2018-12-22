#pragma once

#include <vector>
#include "vector2.h"
#include "logic.h"
#include "objects.h"

class FindBorders
{
public:
	FindBorders(int maxSize, int border, bool isRevert);
	void process(vec2 point);
	void process(circle c);
	void process(const std::vector<vec2>& points);
	void finish(void);
	vec2 toImg(vec2 point) const;
	double toImg(double a) const;
	vec2 fromImg(vec2 point) const;
	vec2 getCalculatedSize(void) const;
	bool isInside(vec2 pos) const;
	vec2 min, max;
private:
	int m_maxSize, m_border;
	double m_scale;
	vec2 m_offset;
	vec2 m_size;
	bool m_isRevert;
};