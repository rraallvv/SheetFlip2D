/*
 *  implicit.cpp
 *  flip2d
 */

#include "implicit.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include "utility.h"
using namespace std;

void implicit::stretchPosition( OBB& obb, FLOAT p[2], FLOAT center[2], FLOAT minimum, FLOAT maximum ) {
	// Compute Relative Pos
	FLOAT rp[2] = { p[0]-center[0], p[1]-center[1] };
	
	// Compute Dot Product
	FLOAT dot[3];
	for( int k=0; k<2; k++ ) {
		dot[k] = rp[0]*obb.u[0][k] + rp[1]*obb.u[1][k];
	}
	
	// Scale By Eigen Value
	for( int k=0; k<2; k++ ) {
		FLOAT r = obb.c[k];
		r = fmax(minimum,fmin(maximum,r));
		dot[k] /= r;
	}
	
	// Compute Final Position
	for( int k=0; k<2; k++ ) {
		p[k] = center[k] + dot[0]*(obb.u[k][0]) + dot[1]*obb.u[k][1];
	}
}

static double implicit_func( vector<particle *> &neighbors, char **DeepZone, FLOAT p[2], FLOAT density, int gn) {
	double phi = 8.0*density/gn;
	int i = min(gn-1,max(0,gn*p[0]));
	int j = min(gn-1,max(0,gn*p[1]));
	if( DeepZone[i][j] ) {
		phi = 0.2*hypotf((i+0.5)/gn-p[0],(j+0.5)/gn-p[1]);
	}
	
	for( int m=0; m<neighbors.size(); m++ ) {
		particle &np = *neighbors[m];
		if( np.type == WALL ) {
			continue;
		}
		FLOAT stretch_pos[2] = { p[0], p[1] };
        implicit::stretchPosition( np.obb, stretch_pos, np.p, 0.25 );
		double d = (1.0/comp_rad(np.level,1.0))*hypotf(np.p[0]-stretch_pos[0],np.p[1]-stretch_pos[1]);
		if( d < phi ) {
			phi = d;
		}
	}
	return phi - 1.0/gn;
}

double implicit::implicit_func( Sorter *sorter, char **DeepZone, FLOAT p[2], FLOAT density ) {
	int gn = sorter->getCellSize();
	// Find Neighbors
#if 0 // ?
	return hypotf(p[0]-0.5,p[1]-0.5)-0.2;
#else
	vector<particle *> neighbors = sorter->getNeigboringParticles_cell(fmax(0,fmin(gn-1,gn*p[0])),
																	 fmax(0,fmin(gn-1,gn*p[1])),
																	 1,1);
	return implicit_func( neighbors, DeepZone, p, density, gn );
#endif
}
