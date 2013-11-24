/*
 *  tension.h
 *  flip2d
 */

#include "common.h"
#include "sorter.h"
#include <vector>

namespace tension {
	void add_surface_tension( sorter *sort, std::vector<particle *> &particles, FLOAT dt, FLOAT min_dens, FLOAT max_dens, FLOAT tension );
}