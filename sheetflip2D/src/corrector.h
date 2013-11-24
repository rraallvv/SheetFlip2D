/*
 *  corrector.h
 *  flip2d
 */

#include "common.h"
#include "sorter.h"
#include <vector>

namespace corrector {
	void resample( sorter *sort, FLOAT p[2], FLOAT u[2], FLOAT re );
	void correct( sorter *sort, std::vector<particle *> &particle, FLOAT dt, FLOAT re, bool anisotropic=false );
};