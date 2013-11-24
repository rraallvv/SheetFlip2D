/*
 *  corrector.cpp
 *  flip2d
 */

#include "corrector.h"
#include "kernel.h"
#include "utility.h"
#include "implicit.h"
#include <math.h>

using namespace std;

#define SPRING		50

void corrector::resample( sorter *sort, FLOAT p[2], FLOAT u[2], FLOAT re ) {
	// Variables for Neighboring Particles
	std::vector<particle *> neighbors;
	int cell_size = sort->getCellSize();
	FLOAT wsum = 0.0;
	FLOAT save[2] = { u[0], u[1] };
	u[0] = u[1] = 0.0;
	
	// Gather Neighboring Particles
	neighbors = sort->getNeigboringParticles_cell(max(0,min(cell_size-1,cell_size*p[0])),
												  max(0,min(cell_size-1,cell_size*p[1])),1,1);
	for( int n=0; n<neighbors.size(); n++ ) {
		particle *np = neighbors[n];
		if( np->type == FLUID ) {
			FLOAT dist2 = hypot2(p[0]-np->p[0],p[1]-np->p[1]);
			FLOAT w = np->m*kernel::sharp_kernel(dist2,re);
			u[0] += w * np->u[0];
			u[1] += w * np->u[1];
			wsum += w;
		}
	}
	if( wsum ) {
		u[0] /= wsum;
		u[1] /= wsum;
	} else {
		u[0] = save[0];
		u[1] = save[1];
	}
}

void corrector::correct( sorter *sort, std::vector<particle *> &particles, FLOAT dt, FLOAT re, bool anisotropic ) {
	// Variables for Neighboring Particles
	sort->sort(particles);
	int cell_size = sort->getCellSize();
	
	// Compute Pseudo Moved Point
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {		
		if( particles[n]->type == FLUID ) {
			particle *p = particles[n];
			FLOAT spring[2] = { 0.0, 0.0 };
			FLOAT x = max(0,min(cell_size,cell_size*p->p[0]));
			FLOAT y = max(0,min(cell_size,cell_size*p->p[1]));
			std::vector<particle *> neighbors = sort->getNeigboringParticles_cell(x,y,1,1);
			for( int n=0; n<neighbors.size(); n++ ) {
				particle *np = neighbors[n];
				if( p != np ) {
              
					FLOAT n0pos[2] = { np->p[0], np->p[1] };
					FLOAT n1pos[2] = { np->p[0], np->p[1] };
					if( anisotropic ) {
						implicit::stretchPosition(p->obb, n0pos, p->p, 0.1, 0.99);
						implicit::stretchPosition(np->obb,n1pos, np->p, 0.1, 0.99);
					}
					
					FLOAT pos[2] = {(FLOAT)((n0pos[0]+n1pos[0])*0.5),(FLOAT)((n0pos[1]+n1pos[1])*0.5)};
					FLOAT dist = hypotf(p->p[0]-pos[0],p->p[1]-pos[1]);
					FLOAT are = (comp_rad(np->level,re)+comp_rad(p->level,re))*0.5;
					FLOAT w = SPRING * kernel::smooth_kernel(dist*dist,are);
					
					if( dist > 0.1*re ) {
						spring[0] += w * (p->p[0]-pos[0]) / dist * are;
						spring[1] += w * (p->p[1]-pos[1]) / dist * are;
					} else {
						if( np->type == FLUID ) {
							spring[0] += 0.01*re/dt*(rand()%101)/100.0;
							spring[1] += 0.01*re/dt*(rand()%101)/100.0;
						} else {
							spring[0] += 0.05*re/dt*np->n[0];
							spring[1] += 0.05*re/dt*np->n[1];
						}
					}
				}
			}
			p->tmp[0][0] = p->p[0] + dt*spring[0];
			p->tmp[0][1] = p->p[1] + dt*spring[1];
		}
	}
	
	// Resample New Velocity
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		if( particles[n]->type == FLUID ) {
			particle *p = particles[n];
			p->tmp[1][0] = p->u[0];
			p->tmp[1][1] = p->u[1];
			resample( sort, p->tmp[0], p->tmp[1], a1*comp_rad(p->level,re) );
		}
	}
	
	// Update
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) { 
		if( particles[n]->type == FLUID ) {
			particle *p = particles[n];
			p->p[0] = p->tmp[0][0];
			p->p[1] = p->tmp[0][1];
			p->u[0] = p->tmp[1][0];
			p->u[1] = p->tmp[1][1];
		}
	}
}