/*
 *  mapper.cpp
 *  flip2d
 */

#include "mapper.h"
#include "utility.h"
#include "kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

#define RE			1.4
#define FOR_EVERY_PARTICLE for( int n=0; n<particles.size(); n++ ) { particle *p = particles[n];

void mapper::mapP2G( sorter *sort, vector<particle *> &particles, FLOAT ***grid, int gn ) {
	
	// Compute Mapping
	OPENMP_FOR FOR_EVERY_CELL(gn+1) {
		
		// Variales for Particle Sorter
		vector<particle *> neighbors;
		
		// Map X Grids
		if( j < gn ) {
			FLOAT px[2] = { (FLOAT)i, (FLOAT)(j+0.5) };
			FLOAT sumw = 0.0;
			FLOAT sumx = 0.0;
			neighbors = sort->getNeigboringParticles_wall(i,j,1,2);
			for( int n=0; n<neighbors.size(); n++ ) {
				particle *p = neighbors[n];
				if( p->type == FLUID ) {
					FLOAT x = fmax(0,fmin(gn,gn*p->p[0]));
					FLOAT y = fmax(0,fmin(gn,gn*p->p[1]));
					FLOAT w = p->m * kernel::sharp_kernel(hypot2(x-px[0],y-px[1]),comp_rad(p->level,RE));
					sumx += w*p->u[0];
					sumw += w;
				}
			}
			grid[0][i][j] = sumw ? sumx/sumw : 0.0;
		}
		
		// Map Y Grids
		if( i < gn ) {
			FLOAT py[2] = { (FLOAT)(i+0.5), (FLOAT)j };
			FLOAT sumw = 0.0;
			FLOAT sumy = 0.0;
			neighbors = sort->getNeigboringParticles_wall(i,j,2,1);
			for( int n=0; n<neighbors.size(); n++ ) {
				particle *p = neighbors[n];
				if( p->type == FLUID ) {
					FLOAT x = fmax(0,fmin(gn,gn*p->p[0]));
					FLOAT y = fmax(0,fmin(gn,gn*p->p[1]));
					FLOAT w = p->m * kernel::sharp_kernel(hypot2(x-py[0],y-py[1]),comp_rad(p->level,RE));
					sumy += w*p->u[1];
					sumw += w;
				}
			}
			grid[1][i][j] = sumw ? sumy/sumw : 0.0;
		}
	} END_FOR;
}

void mapper::mapG2P( vector<particle *> &particles, FLOAT ***grid, int gn ) {
	OPENMP_FOR FOR_EVERY_PARTICLE {
		fetchVelocity( p->p, p->u, grid, gn );
	} END_FOR;
}

void mapper::fetchVelocity_RK2( char **A, FLOAT p[2], FLOAT u[2], FLOAT ***grid, int gn, FLOAT dt ) {
	FLOAT p1[2];
	FLOAT u1[2];
	fetchVelocity( p, u, grid, gn );
	for( int n=0; n<2; n++ ) p1[n] = p[n]+dt*u[n];
	
	int i = fmax(0,fmin(gn-1,gn*p1[0]));
	int j = fmax(0,fmin(gn-1,gn*p1[1]));
	if( A[i][j] == FLUID ) {
		fetchVelocity( p1, u1, grid, gn );
		for( int n=0; n<2; n++ ) u[n] = 0.5*(u[n]+u1[n]);
	}
}

static FLOAT linear ( FLOAT **d, FLOAT x, FLOAT y, int w, int h ) {
	x = fmax(0.0,fmin(w,x));
	y = fmax(0.0,fmin(h,y));
	int i = min(x,w-2);
	int j = min(y,h-2);
	
	return ((i+1-x)*d[i][j]+(x-i)*d[i+1][j])*(j+1-y) + ((i+1-x)*d[i][j+1]+(x-i)*d[i+1][j+1])*(y-j);
}

void mapper::fetchVelocity( FLOAT p[2], FLOAT u[2], FLOAT ***grid, int gn ) {
#if 1
	u[0] = linear( grid[0], gn*p[0], gn*p[1]-0.5, gn+1, gn );
	u[1] = linear( grid[1], gn*p[0]-0.5, gn*p[1], gn, gn+1 );
#else
	FLOAT x = fmax(0,fmin(gn-1,gn*p[0]));
	FLOAT y = fmax(0,fmin(gn-1,gn*p[1]));
	int i = x;
	int j = y;
	u[0] = (1.0-(x-i))*grid[0][i][j] + (x-i)*grid[0][i+1][j];
	u[1] = (1.0-(y-j))*grid[1][i][j] + (y-j)*grid[1][i][j+1];
#endif
}
















