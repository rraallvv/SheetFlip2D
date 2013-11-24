/*
 *  implicit.h
 *  flip2d
 */

#include "common.h"
#include "sorter.h"
#include <vector>

namespace implicit {
    void stretchPosition( OBB& obb, FLOAT p[2], FLOAT center[2], FLOAT minimum=0.1, FLOAT maximum=10.0 );
	double implicit_func( sorter *sort, char **DeepZone, FLOAT p[2], FLOAT density );
};