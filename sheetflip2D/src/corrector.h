/*
 *  corrector.h
 *  flip2d
 */

#include "common.h"
#include "sorter.h"
#include <vector>

namespace corrector {
	void resample( Sorter *sorter, FLOAT p[2], FLOAT u[2], FLOAT re );
	void correct( Sorter *sorter, std::vector<particle *> &particle, FLOAT dt, FLOAT re, bool anisotropic=false );
};