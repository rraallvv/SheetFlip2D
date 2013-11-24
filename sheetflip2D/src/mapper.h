/*
 *  mapper.h
 *  flip2d
 */

#include "common.h"
#include "sorter.h"
#include <vector>

class sorter;
namespace mapper {
	void mapP2G( sorter *sort, std::vector<particle *> &particles, FLOAT ***grid, int gn );
	void mapG2P( std::vector<particle *> &particles, FLOAT ***grid, int gn );
	void fetchVelocity_RK2( char **A, FLOAT p[2], FLOAT u[2], FLOAT ***grid, int gn, FLOAT dt );
	void fetchVelocity( FLOAT p[2], FLOAT u[2], FLOAT ***grid, int gn );
}