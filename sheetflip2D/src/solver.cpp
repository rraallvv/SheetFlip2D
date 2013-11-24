/*
 *  solver.cpp
 *  smoke
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solver.h"
#include "utility.h"

#undef a1
#undef a2

#define USE_PRECOND   1
static char subcell = 0;

// Clamped Fetch
inline static FLOAT ref( FLOAT **x, int pos[2], int n ) {
	if( pos[0] < 0 || pos[0] > n-1 || pos[1] < 0 || pos[1] > n-1 ) return 0.0;
	return x[pos[0]][pos[1]];
}

// Ans = Ax
static void compute_Ax( FLOAT **M[], FLOAT **x, FLOAT **ans, int n ) {
	FOR_EVERY_CELL(n) {
		int q[][2] = { {i,j}, {i-1,j}, {i,j-1}, {i+1,j}, {i,j+1} };
		ans[i][j] = 0.0;
		if( M[0][i][j] ) {
			for( int m=0; m<5; m++ ) ans[i][j] += M[m][i][j]*ref(x,q[m],n);
		}
	} END_FOR
}

// ans = x^T * x
static double product( FLOAT **M[], FLOAT **x, FLOAT **y, int n ) {
	double ans = 0.0;
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			if( M[0][i][j] ) ans += x[i][j]*y[i][j];
		}
	}
	return ans;
}

// x = 0
static void clear( FLOAT **x, int n ) {
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			x[i][j] = 0.0;
		}
	}
}

static void flip( FLOAT **x, int n ) {
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			x[i][j] = -x[i][j];
		}
	}
}

// x <= y
static void copy( FLOAT **x, FLOAT **y, int n ) {
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			x[i][j] = y[i][j];
		}
	}
}
				 
// Ans = x + a*y
static void op( FLOAT **M[], FLOAT **x, FLOAT **y, FLOAT **ans, FLOAT a, int n ) {
	static FLOAT **tmp = alloc2D<FLOAT>(n);
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			if( M[0][i][j] ) tmp[i][j] = x[i][j]+a*y[i][j];
            else tmp[i][j] = 0.0;
		}
	}
	copy(ans,tmp,n);
}

// r = b - Ax
static void residual( FLOAT **M[], FLOAT **x, FLOAT **b, FLOAT **r, int n ) {
	compute_Ax(M,x,r,n);
	op( M, b, r, r, -1.0, n );
}

static inline FLOAT square( FLOAT a ) {
	return a*a;
}

static void buildPreconditioner( FLOAT **P, FLOAT **M[], int n ) {
	clear(P,n);
	FLOAT a = 0.25;
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			int q[][2] = { {i,j}, {i-1,j}, {i,j-1}, {i+1,j}, {i,j+1} };
			P[i][j] = 0.0;
			if( M[0][i][j] ) {
				double diag = M[0][i][j];
				double e = diag;
				for( int m=1; m<3; m++ ) e -= square(M[m][i][j]*ref(P,q[m],n));
				if( e < a*diag ) e = diag;
				P[i][j] = 1.0/sqrt(e);
			}
		}
	}
}

static void applyPreconditioner( FLOAT **z, FLOAT **r, FLOAT **P, FLOAT **M[], int n ) {
#if USE_PRECOND
	static FLOAT **qq = alloc2D<FLOAT>(n);
	clear(qq,n);
	
	// Lq = r
	for( int i=0; i<n; i++ ) {
		for( int j=0; j<n; j++ ) {
			int q[][2] = { {i,j}, {i-1,j}, {i,j-1}, {i+1,j}, {i,j+1} };
			if( M[0][i][j] ) {
				double down_neibors = 0.0;
				for( int m=1; m<3; m++ ) down_neibors += M[m][i][j]*ref(P,q[m],n)*ref(qq,q[m],n);
				double t = r[i][j] - down_neibors;
				qq[i][j] = t*P[i][j];
			} else {
				qq[i][j] = 0.0;
			}
		}
	}
	
	// L^T z = q
	for( int i=n-1; i>=0; i-- ) {
		for( int j=n-1; j>=0; j-- ) {
			int q[][2] = { {i,j}, {i-1,j}, {i,j-1}, {i+1,j}, {i,j+1} };
			if( M[0][i][j] ) {
				double top_neibors = 0.0;
				for( int m=3; m<5; m++ ) top_neibors += M[m][i][j]*ref(P,q[0],n)*ref(z,q[m],n);
				double t = qq[i][j] - top_neibors;
				z[i][j] = t*P[i][j];
			} else {
				z[i][j] = 0.0;
			}
		}
	}
#else
	copy(z,r,n);
#endif
}

static void conjGrad( FLOAT **M[], FLOAT **P, FLOAT **p, FLOAT **b, int n ) {
	// Pre-allocate Memory
	static FLOAT **r = alloc2D<FLOAT>(n);
	static FLOAT **z = alloc2D<FLOAT>(n);
	static FLOAT **s = alloc2D<FLOAT>(n);
	
	compute_Ax( M, p, z, n );					// z = applyA(x)
	op( M, b, z, r, -1.0, n );                  // r = b-Ax
	applyPreconditioner(z,r,P,M,n);				// Apply Conditioner z = f(r)
	copy(s,z,n);								// s = z
	
	double a = product( M, z, r, n );			// a = z . r
    double error2;
	for( int k=0; k<n*n; k++ ) {
		compute_Ax( M, s, z, n );				// z = applyA(s)
		double alpha = a/product( M, z, s, n );	// alpha = a/(z . s)
		op( M, p, s, p, alpha, n );				// p = p + alpha*s
		op( M, r, z, r, -alpha, n );			// r = r - alpha*z;
		error2 = product( M, r, r, n );         // error2 = r . r
        
		if( error2/(n*n) < 1.0e-6 ) {
            // printf( "Converged. error2 = %e. k/(n*n) = %f\n", error2/(n*n), k/(double)(n*n) );
            return;
        }
		applyPreconditioner(z,r,P,M,n);			// Apply Conditioner z = f(r)
		double a2 = product( M, z, r, n );		// a2 = z . r
		double beta = a2/a;
		op( M, z, s, s, beta, n );				// s = z + beta*s
		a = a2;
	}
    
    printf( "Not Converged. error2 = %e\n", error2/(n*n) );
    exit(0);
}

static void buildMatrix( char **A, Vector2 **F[], FLOAT **M[], int n ) {
	FLOAT scale = n*n;
	FOR_EVERY_CELL(n) {
		if( A[i][j]==FLUID ) {
			int q[][2] = { {i,j}, {i-1,j}, {i,j-1}, {i+1,j}, {i,j+1} };
			M[0][i][j] = 4.0*scale;		// Diag
			for( int m=1; m<5; m++ ) {	// Neighbors
				int ni =  q[m][0];
				int nj =  q[m][1];
				M[m][i][j] = 0.0;
				bool outRegion = ( ni<0 || ni>=n || nj<0 || nj>=n );
				if( outRegion ) {
					M[0][i][j] -= scale;
				} else if( A[ni][nj]==FLUID ) {
					M[m][i][j] = -scale;
				} else {
					if( subcell ) {
						Vector2 Fv[] = { Vector2() ,F[0][i][j], F[1][i][j], F[0][i+1][j], F[1][i][j+1] };
						M[0][i][j] -= scale*Fv[m][0];
					} else {
						M[0][i][j] -= scale*(A[ni][nj]==WALL);
					}
				}
			}
		} else {
			for( int m=0; m<5; m++ ) M[m][i][j] = 0.0;
		}
	} END_FOR;
}

void solver::setSubcell( int value ) {
	subcell = value;
}

void solver::solve( char **A, Vector2 **F[], FLOAT **p, FLOAT **b, int n ) {
	static FLOAT ***M = NULL;
#if USE_PRECOND
	static FLOAT **P = alloc2D<FLOAT>(n);
#else
	static FLOAT **P = NULL;
#endif
	
	// Allocate Matrix Entry
	if( ! M ) {
		M = new FLOAT **[5];
		for( int i=0; i<5; i++ ) M[i] = alloc2D<FLOAT>(n);
	}
	
	// Build Matrix
	buildMatrix( A, F, M, n );

	// Flip Divergence
	flip(b,n);
	
#if USE_PRECOND
	// Build Modified Incomplete Cholesky Precondioner Matrix
	buildPreconditioner(P,M,n);
#endif
	
	// Conjugate Gradient Method
	conjGrad(M,P,p,b,n);
}