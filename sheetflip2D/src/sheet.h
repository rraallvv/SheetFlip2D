/*
 *  sheet.h
 *  flip2d
 */

#include <vector>
#include "sorter.h"
#include "common.h"

namespace sheet {
	void keepThinSheet( char **A, sorter *sort, std::vector<particle *> &particles, FLOAT obb_rate, FLOAT density, FLOAT wall_thick );
	void collapseThinSheet( char **A, sorter *sort, std::vector<particle *> &particles, FLOAT obb_rate, FLOAT density );
}