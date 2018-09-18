#include <cmath>
#include "common.h"

//-----------------------------------------------------------------------------
double random(void) {
	return std::rand() / double(RAND_MAX);
}

//-----------------------------------------------------------------------------
int intRandom(int min, int max) {
	return min + random() * (max - min);
}