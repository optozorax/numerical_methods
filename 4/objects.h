#pragma once

#include "main.h"

struct point { double x, y; };
struct circle { point c; double r; };
struct line { point a, b; };

pair<fnm_t, jnm_t> two_circles(circle a, circle b);
pair<fnm_t, jnm_t> two_circles_and_line(circle a, circle b, line l);
pair<fnm_t, jnm_t> three_lines(line a, line b, line c);

pair<fnm_t, jnm_t> sin_and_line(point b); // a is {0, 0}