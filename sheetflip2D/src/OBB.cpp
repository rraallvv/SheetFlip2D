/*
 *  OBB.cpp
 *  flip2d
 */

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

void SingularValueDecomposition(FLOAT mIn[2][2], FLOAT vOut[2], FLOAT mOut[2][2]) // pointer to top-left corner
{
	double A[2][2];
	double U[2][2];
	double V[2][2];
    double S[2];
	
	for (int j=0; j<2; j++) {
		for (int i=0; i<2; i++) {
			A[j][i] = mIn[i][j];
		}
	}
	
    char JOBU='A';
    char JOBVT='A';
    int LWORK=-1;
    double test;
    int INFO;
	
	int N = 2;
	
    // Allocate memory
    dgesvd_(&JOBU, &JOBVT, &N, &N,
			(double*)A, &N,
			(double*)S,
			(double*)U, &N,
			(double*)V, &N,
			&test, &LWORK, &INFO);
    LWORK=test;
    int size=int(test);
    double WORK[size];
	
    // Compute SVD
    dgesvd_(&JOBU, &JOBVT, &N, &N,
			(double*)A, &N,
			(double*)S,
			(double*)U, &N,
			(double*)V, &N,
			&WORK[0], &LWORK, &INFO);
	
    // Output as doubles
	for (int i=0; i<2; i++) {
		vOut[i] = S[i];
	}
	
	for (int j=0; j<2; j++) {
		for (int i=0; i<2; i++) {
			mOut[i][j] = V[i][j];
		}
	}
}

OBB buildOBB( vector<particle *> particles, FLOAT cp[2], FLOAT re ) {
	
    int pn = (int)particles.size();
    
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
	
	// Compute Eigen Vectors
	FLOAT r[2];
	FLOAT v[2][2];
	
	SingularValueDecomposition(A, r, v);
	
	// Extract Eigen Vectors
	OBB obb;
	obb.u[0][0] = v[0][0];
	obb.u[0][1] = v[1][0];
	obb.u[1][0] = v[0][1];
	obb.u[1][1] = v[1][1];
	
	// Put Eigen Value
	for( int i=0; i<2; i++ ) {
		obb.c[i] = r[i] / (0.1*re*re);
	}
    
	FLOAT a = r[0];
	obb.ratio = a ? obb.ratio = r[1] / a : 0.0;
	
	return obb;
}
