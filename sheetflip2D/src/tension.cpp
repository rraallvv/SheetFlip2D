/*
 *  tension.cpp
 *  flip2d
 */

#include "tension.h"
#include "utility.h"
#include <math.h>
using namespace std;

void tension::add_surface_tension( sorter *sort, vector<particle *> &particles, FLOAT dt, FLOAT min_dens, FLOAT max_dens, FLOAT tension ) {	
	int gn = sort->getCellSize();
	static FLOAT **color = NULL;
	static FLOAT **color_div = NULL;
	static FLOAT ***color_grad = NULL;
	FLOAT h = 1.0/gn;
	
	// Allocate Memory For Color Grid
	if( ! color ) {
		color = alloc2D<FLOAT>(gn);
		color_div = alloc2D<FLOAT>(gn);
		color_grad = new FLOAT **[3];
		color_grad[0] = alloc2D<FLOAT>(gn);
		color_grad[1] = alloc2D<FLOAT>(gn);
	}
	
	// First Compute Color On Grid
	OPENMP_FOR FOR_EVERY_CELL(gn) {
		FLOAT p[2] = { (FLOAT)((i+0.5)*h), (FLOAT)((j+0.5)*h) };
		vector<particle *> neighbors = sort->getNeigboringParticles_cell(i,j,1,1);
		FLOAT wsum = 0.0;
		for( int m=0; m<neighbors.size(); m++ ) {
			particle np = *neighbors[m];
			if( np.type == WALL ) continue;
			FLOAT d = gn*hypotf( np.p[0]-p[0], np.p[1]-p[1] );
			wsum += expf(-d);
		}
		color[i][j] = wsum;
	} END_FOR
	
	// Next Compute Gradient On Grid
	OPENMP_FOR FOR_EVERY_CELL(gn) {		
		color_grad[0][i][j] = (color[min(gn-1,i+1)][j] - color[max(0,i-1)][j]) / (2.0*h);
		color_grad[1][i][j] = (color[i][min(gn-1,j+1)] - color[i][max(0,j-1)]) / (2.0*h);
	} END_FOR
	
	// Map Gradient Onto Particles
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		particle &p = *particles[n];
		if( p.type == WALL ) continue;
		if( p.dens > max_dens || p.dens < min_dens ) {
			p.n[0] = 0.0;
			p.n[1] = 0.0;
		} else {
			FLOAT x = fmax(0.5,fmin(gn-1.5,gn*p.p[0]))-0.5;
			FLOAT y = fmax(0.5,fmin(gn-1.5,gn*p.p[1]))-0.5;
			int i = x;
			int j = y;		
			for( int m=0; m<2; m++ ) {
				p.n[m] = ( (x-i)*color_grad[m][i+1][j+1]+(1-x+i)*color_grad[m][i][j+1] ) * (y-j) +
				( (x-i)*color_grad[m][i+1][j]+(1-x+i)*color_grad[m][i][j] ) * (1-y+j);
			}
		}
	}
	
	// Compute Divergence On Grid
	OPENMP_FOR FOR_EVERY_CELL(gn) {		
		color_div[i][j] =	color[min(gn-1,i+1)][j] +
							color[max(0,i-1)][j] +
							color[i][min(gn-1,j+1)] +
							color[i][max(0,j-1)] - 4.0*color[i][j];
		color_div[i][j] /= h*h;
	} END_FOR
	
	// Map Divergence Onto Particles
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		particle &p = *particles[n];
		
		vector<particle *> neighbors = sort->getNeigboringParticles_cell(fmax(0,fmin(gn-1,gn*p.p[0])),
																		 fmax(0,fmin(gn-1,gn*p.p[1])),1,1);
		int nump = 0;
		for( int i=0; i<neighbors.size(); i++ ) nump += neighbors[i]->type == FLUID;
			
		p.curv = 0.0;
		if( p.type == WALL ) continue;
		if( p.dens > max_dens || p.dens < min_dens || nump < 2 ) {
		} else {
			FLOAT x = fmax(0.5,fmin(gn-1.5,gn*p.p[0]))-0.5;
			FLOAT y = fmax(0.5,fmin(gn-1.5,gn*p.p[1]))-0.5;
			int i = x;
			int j = y;		
			FLOAT nlen = hypotf(p.n[0],p.n[1]);
			FLOAT c =	((x-i)*color_div[i+1][j+1]+(1-x+i)*color_div[i][j+1] ) * (y-j) +
						((x-i)*color_div[i+1][j]+(1-x+i)*color_div[i][j] ) * (1-y+j);
			FLOAT curv = 0.0;
			if( nlen > 20.0 ) {
				curv = -c/nlen;
				p.n[0] /= nlen;
				p.n[1] /= nlen;
			}
			
			// Add Surace Tension
			p.u[0] += dt*curv*tension*p.n[0];
			p.u[1] += dt*curv*tension*p.n[1];
			p.curv = curv;
		}
	}
}

