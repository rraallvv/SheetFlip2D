/*
 *  OBB.cpp
 *  flip2d
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include <Accelerate/Accelerate.h>

#include "OBB.h"
#include "kernel.h"
#include "utility.h"

using namespace std;

OBB isotropicOBB() {
	OBB obb;
	for( int k=0; k<2; k++ ) {
		obb.u[k][0] = k==0;
		obb.u[k][1] = k==1;
        obb.c[k] = 1.0;
	}
    obb.ratio = 1.0;
	return obb;
}

void buildOBB( char **A, char **DeepZone, Sorter *sorter, std::vector<particle *> &particles ) {
	int gn = sorter->getCellSize();
	FLOAT min_r = 1.0/gn;
	
	// Build OBB
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		particle &p = *particles[n];
		
		int i = fmin(gn-1,fmax(0,p.p[0]*gn));
		int j = fmin(gn-1,fmax(0,p.p[1]*gn));
		
		if( A[i][j]==WALL || DeepZone[i][j] || p.level > 1 || p.dens > a7 ) p.obb = isotropicOBB();
		else if( p.type == FLUID ) p.obb = buildNearbyOBB(sorter,&p,2.0*min_r);
	}
}

OBB buildNearbyOBB( Sorter *sorter, particle *p, FLOAT re ) {
	std::vector<particle *> neighbors;
	int cell_size = sorter->getCellSize();

	// Gather Neighboring Particles
	neighbors = sorter->getNeigboringParticles_cell(max(0,min(cell_size-1,cell_size*p->p[0])),
												  max(0,min(cell_size-1,cell_size*p->p[1])),2,2);
	
	std::vector<particle *> final_neighbors;
	for( int n=0; n<neighbors.size(); n++ ) {
		particle &np = *neighbors[n];
		if( np.type == FLUID && hypotf(np.p[0]-p->p[0],np.p[1]-p->p[1]) < re ) {
			final_neighbors.push_back(&np);
		}
	}
	
	return buildOBB(final_neighbors,p->p,re);
}

OBB buildOBB( vector<particle *> particles, FLOAT cp[2], FLOAT re ) {
	
    int pn = particles.size();
    
	// Compute Center Position
    FLOAT c[2] = { 0.0, 0.0 };
    {
        FLOAT wsum = 0.0;
        for( int n=0; n<pn; n++ ) {
            particle &p = *particles[n];
            FLOAT w = kernel::smooth_kernel(hypot2(p.p[0]-cp[0],p.p[1]-cp[1]),re);
            c[0] += w*p.p[0];
            c[1] += w*p.p[1];
            wsum += w;
        }
        c[0] /= wsum;
        c[1] /= wsum;
    }
    
	// Allocate 2x2 Matrix
	gsl_matrix *M = gsl_matrix_alloc(2,2);
	
	// Compute Variance-covariance Matrix
	FLOAT A[2][2] = { {0.0, 0.0}, {0.0, 0.0}};
	FOR_EVERY_CELL(2) {
        FLOAT wsum = 0.0;
		for( int n=0; n<pn; n++ ) {
			particle &p = *particles[n];
            FLOAT w = kernel::smooth_kernel(hypot2(p.p[0]-c[0],p.p[1]-c[1]),re);
			A[i][j] += w*(p.p[i]-c[i])*(p.p[j]-c[j]);
            wsum += w;
		}
        if( wsum ) A[i][j] /= wsum;
	} END_FOR;
	
	// Fill In
	FOR_EVERY_CELL(2) {
		gsl_matrix_set(M, i, j, A[i][j] );
	} END_FOR;
	
	// Compute Eigen Vectors
	gsl_vector *r = gsl_vector_alloc(2);
	gsl_vector *w = gsl_vector_alloc(2);
	gsl_matrix *v = gsl_matrix_alloc(2,2);
	gsl_linalg_SV_decomp(M, v, r, w);
		
	// Extract Eigen Vectors
	OBB obb;
	if (gsl_vector_get(r,0) > gsl_vector_get(r,1)) {
		obb.u[0][0] = gsl_matrix_get(v, 0, 0);
		obb.u[0][1] = gsl_matrix_get(v, 0, 1);
		obb.u[1][0] = gsl_matrix_get(v, 1, 0);
		obb.u[1][1] = gsl_matrix_get(v, 1, 1);
	} else {
		obb.u[0][0] = gsl_matrix_get(v, 1, 0);
		obb.u[0][1] = gsl_matrix_get(v, 1, 1);
		obb.u[1][0] = gsl_matrix_get(v, 0, 0);
		obb.u[1][1] = gsl_matrix_get(v, 0, 1);
	}

	// Put Eigen Value
	for( int i=0; i<2; i++ ) {
		obb.c[i] = gsl_vector_get( r, i ) / (0.1*re*re);
	}
    
	FLOAT a = gsl_vector_get( r, 0 );
	obb.ratio = a ? obb.ratio = gsl_vector_get( r, 1 ) / a : 0.0;
	
	// Deallocation
	gsl_matrix_free(M);
	gsl_matrix_free(v);
	gsl_vector_free(r);
	gsl_vector_free(w);
	// gsl_eigen_symmv_free(w);
	
	return obb;
}