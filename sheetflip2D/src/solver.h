/*
 *  solver.h
 *  smoke
 *
 */

#include "common.h"

namespace solver {
	// Solve Ax = b
	void setSubcell( int value );
	void solve( char **A, Vector2 **F[], FLOAT **p, FLOAT **b, int n );
}
