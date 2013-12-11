//
//  adaptive.h
//  flip2d
//

#include "common.h"
#include "sorter.h"

namespace adaptive {
	void initial_merge( Sorter *sorter, char **A, std::vector<particle *> &particles, FLOAT density );
	void computeDeepZone( char **A, char **DeepZone, int gn );
    char** resample( Sorter *sorter, char **A, std::vector<particle *> &particles, FLOAT density, bool onlysplit=false );
};