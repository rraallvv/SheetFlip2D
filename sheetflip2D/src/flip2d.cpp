/*
 *  flip2d.cpp
 *  smoke
 *
 */

#include "flip2d.h"
#include "solver.h"
#include "utility.h"
#include "mapper.h"
#include "sorter.h"
#include "corrector.h"
#include "tension.h"
#include "sheet.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "write_bmp.h"
#include "kernel.h"
//#include "SPH2D.h"
#include "implicit.h"
#include "OBB.h"
#include "adaptive.h"

#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <sys/time.h>
#elif defined(WIN32)
#include "glut.h"
#include <windows.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#include <sys/time.h>
#endif

using namespace std;

#define		N			25
#define		DT			1.2e-2
#define		FLIP		1
#define     TEST        0
#define		ALPHA		0.99
#define		GRAVITY		9.8
#define		DENSITY		0.5
#define		OBB_RATE	0.1
#define		RECORDING	0
#define		WALL_THICK	(1.0/N)
#define		SPH			0
#define     DIM         2

static FLOAT ***u = NULL;		// Access Bracket u[DIM][X][Y] ( Staggered Grid )
static FLOAT ***u_save = NULL;	// Access Bracket u[DIM][X][Y] ( Staggered Grid )
static FLOAT **p = NULL;		// Equivalent to p[N][N]
static FLOAT **d = NULL;		// Equivalent to d[N][N]
static FLOAT **L = NULL;		// Equivalent to L[N][N]
static FLOAT **SL = NULL;		// Equivalent to SL[N][N] (Surface Levelset)
static Vector2 ***F = NULL;     // Access Bracket u[DIM][X][Y][DIM] ( Staggered Grid )
static char **A = NULL;
static FLOAT ***wall_normal = NULL;
static Sorter *sorter = NULL;
static FLOAT max_dens = 0.0;
static FLOAT volume0 = 0.0;
static FLOAT y_volume0 = 0.0;

static bool place_sphere = 1;
static FLOAT sphere_r = 0.3;
static FLOAT sphere_pos[2] = { 1.0, 0.0 };

static vector<particle *> particles;
static int timeStep = 0;
static int doInit = 0;
static int init_num = 0;

// Variables
static int correction = 1;
static int splitting = 0;
static int subcell = 1;
static int drawOval = 1;
static int anisotropic_spring = 1;
static int adaptive_sampling = 1;
static int correct_volume = 1;
static int hide_help = 0;
static char **DeepZone = NULL;
static char **DeepArea = NULL;

static unsigned long getMicroseconds() {
#if defined(_WIN32)
	LARGE_INTEGER nFreq, Time;
	QueryPerformanceFrequency(&nFreq);
	QueryPerformanceCounter(&Time);
	return (double)Time.QuadPart / nFreq.QuadPart * 1000000;
#else
	struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec*1000000 + tv.tv_usec;
#endif
}

static double dumptime() {
	static unsigned prevMicroSec = getMicroseconds();
	unsigned curMicroSec = getMicroseconds();
	double res = (curMicroSec - prevMicroSec)/1000000.0;
	prevMicroSec = curMicroSec;
	return res;
}

static void compute_wall_normal() {
	// Sort Particles
	sorter->sort(particles);
	
	// Compute Wall Normal
	for( int n=0; n<particles.size(); n++ ) {
		particle *p = particles[n];
		p->n[0] = p->n[1] = 0.0;
		int i = fmin(N-1,fmax(0,p->p[0]*N));
		int j = fmin(N-1,fmax(0,p->p[1]*N));
		if( p->type == WALL ) {
			if( p->p[0] <= 1.1*WALL_THICK ) {
				p->n[0] = 1.0;
			} 
			if( p->p[0] >= 1.0-1.1*WALL_THICK ) {
				p->n[0] = -1.0;
			} 
			if( p->p[1] <= 1.1*WALL_THICK ) {
				p->n[1] = 1.0;
			} 
			if( p->p[1] >= 1.0-1.1*WALL_THICK ) {
				p->n[1] = -1.0;
			} 
			
			if( p->n[0] == 0.0 && p->n[1] == 0.0 ) {
				vector<particle *> neighbors = sorter->getNeigboringParticles_cell(i,j,1,1);
				for( int n=0; n<neighbors.size(); n++ ) {
					particle *np = neighbors[n];
					if( p!=np && np->type == WALL ) {
						FLOAT d = hypotf(p->p[0]-np->p[0],p->p[1]-np->p[1]);
						FLOAT w = 1.0/d;
						p->n[0] += w*(p->p[0]-np->p[0])/d;
						p->n[1] += w*(p->p[1]-np->p[1])/d;
					}
				}
			}
			
			if( place_sphere ) {
				if( i>0 && i<N-1 && j>0 && j<N-1 ) {
					if( hypot(p->p[0]-sphere_pos[0],p->p[1]-sphere_pos[1])<sphere_r ) {
						p->n[0] = p->p[0]-sphere_pos[0];
						p->n[1] = p->p[1]-sphere_pos[1];
					}
				}
			}
		}
		FLOAT d = hypotf(p->n[0],p->n[1]);
		if( d ) {
			p->n[0] /= d;
			p->n[1] /= d;
		}
	}
	
	FOR_EVERY_CELL(N) {
		if( A[i][j] == WALL ) {
			vector<particle *> array = sorter->getNeigboringParticles_cell(i,j,0,0);
			for( int m=0; m<array.size(); m++ ) {
				if( array[m]->type == WALL ) {
					wall_normal[0][i][j] = array[m]->n[0];
					wall_normal[1][i][j] = array[m]->n[1];
				}
			}
		}
	} END_FOR;
	
	// Extrapolate Wall Normal
	for( int cnt=0; cnt<10; cnt++ ) {
		OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
			particle *p = particles[n];
			if( p->type == WALL && p->n[0] == 0.0 && p->n[1] == 0.0 ) {
				int i = fmin(N-1,fmax(0,p->p[0]*N));
				int j = fmin(N-1,fmax(0,p->p[1]*N));
				vector<particle *> neighbors = sorter->getNeigboringParticles_cell(i,j,1,1);
				p->tmp[0][0] = 0.0;
				p->tmp[0][1] = 0.0;
				for( int n=0; n<neighbors.size(); n++ ) {
					particle *np = neighbors[n];
					if( p!=np && np->type == WALL ) {
						p->tmp[0][0] += np->n[0];
						p->tmp[0][1] += np->n[1];
					}
				}
			}
		}
		
		OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
			particle *p = particles[n];
			if( p->type == WALL && p->n[0] == 0.0 && p->n[1] == 0.0 ) {
				FLOAT d = hypotf(p->tmp[0][0],p->tmp[0][1]);
				if( d ) {
					int i = fmin(N-1,fmax(0,p->p[0]*N));
					int j = fmin(N-1,fmax(0,p->p[1]*N));
					p->n[0] = p->tmp[0][0]/d;
					p->n[1] = p->tmp[0][1]/d;
					wall_normal[0][i][j] = p->n[0];
					wall_normal[1][i][j] = p->n[1];
				}
			}
		}
	}
	
	sorter->sort(particles);
	sorter->markWater(A,DENSITY);
}

static void computeDensity() {
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		// Find Neighbors
		int gn = sorter->getCellSize();
		if( particles[n]->type == WALL ) {
			particles[n]->dens = 1.0;
			continue;
		}
		
		FLOAT *p = particles[n]->p;
		
		vector<particle *> neighbors = sorter->getNeigboringParticles_cell(fmax(0,fmin(gn-1,gn*p[0])),
																		 fmax(0,fmin(gn-1,gn*p[1])),1,1);
		FLOAT wsum = 0.0;
		for( int m=0; m<neighbors.size(); m++ ) {
			particle np = *neighbors[m];
			if( np.type == WALL ) continue;
			FLOAT d2 = hypot2( np.p[0]-p[0], np.p[1]-p[1] );
			FLOAT w = np.m*kernel::smooth_kernel(d2, a2*DENSITY/N);
			wsum += w;
		}
		particles[n]->dens = wsum / max_dens;
	}
}

static void pushParticle( FLOAT pos[3], FLOAT u[3], char type ) {
	particle *p = new particle;
	//FLOAT xj = (type==FLUID)*0.001*((rand()%101)/50.0-1.0)/N;
	//FLOAT yj = (type==FLUID)*0.001*((rand()%101)/50.0-1.0)/N;
	p->p[0] = pos[0];//+xj;
	p->p[1] = pos[1];//+yj;
	p->u[0] = u[0];
	p->u[1] = u[1];
	p->n[0] = 0.0;
	p->n[1] = 0.0;
	p->f[0] = 0.0;
	p->f[1] = 0.0;
	p->dens = 1.0;
	p->curv = 0.0;
	p->type = type;
	p->mark = 0;
	p->pref = NULL;
	p->m = 1.0;
	p->level = 1;
	p->split_cnt = 0;
	mass aMass = { 1.0, 0, p };
	p->hosts.push_back(aMass);
	particles.push_back(p);
}

static int waterDrop( FLOAT x, FLOAT y ) {
	// place_sphere = 1;
	// if( hypot(x-sphere_pos[0],y-sphere_pos[1]) < sphere_r ) {
	//	return WALL;
	// }
	
	if( (x < 1.0 && y < 0.75) /*|| hypot(x-0.5,y-0.75) < 0.12*/ ) return FLUID;
	else return AIR;
}

static int sliderDambreak( FLOAT x, FLOAT y ) {
	if( x < 0.4 && y < 0.55 ) return FLUID;
	if( y-x < -0.65 ) return WALL;
	return AIR;
}

static int towerDrop( FLOAT x, FLOAT y ) {
	if( x > 0.4 && x < 0.6 && y < 0.4 ) return WALL;
	if( hypot(x-0.5,y-0.8) < 0.1 ) return FLUID;
	else return AIR;
}

static int particleSetting( FLOAT x, FLOAT y ) {
	place_sphere = 0;
	if( init_num % 2 == 0 ) {
		return sliderDambreak(x,y);
	} else {
		return waterDrop(x,y);
	}
	return WALL;
}

void flip2d::init() {
	timeStep = 0;
	volume0 = 0.0;
	
	// Allocate Variables
	if( ! p ) p = alloc2D<FLOAT>(N);	
	if( ! d ) d = alloc2D<FLOAT>(N);
	if( ! A ) A = alloc2D<char>(N);
	if( ! L ) L = alloc2D<FLOAT>(N);
	if( ! SL ) SL = alloc2D<FLOAT>(N);
	if( ! u ) {
		u = new FLOAT **[3];
		u[0] = alloc2D<FLOAT>(N+1);
		u[1] = alloc2D<FLOAT>(N+1);
		
		u_save = new FLOAT **[3];
		u_save[0] = alloc2D<FLOAT>(N+1);
		u_save[1] = alloc2D<FLOAT>(N+1);
	}
    if( ! F ) {
        F = new Vector2 **[2];
        F[0] = alloc2D<Vector2>(N+1);
        F[1] = alloc2D<Vector2>(N+1);
    }
	if( ! wall_normal ) {
		wall_normal = new FLOAT **[2];
		wall_normal[0] = alloc2D<FLOAT>(N);
		wall_normal[1] = alloc2D<FLOAT>(N);
	}
	
	// Clear Variables
	FOR_EVERY_X_FLOW(N) {
		u_save[0][i][j] = u[0][i][j] = 0.0;
	} END_FOR
	
	FOR_EVERY_Y_FLOW(N) {
		u_save[1][i][j] = u[1][i][j] = 0.0;
	} END_FOR
	
	// Allocate And Place Particles
	for( int n=0; n<particles.size(); n++ ) {
		delete particles[n];
	}
	particles.clear();
	
	// Allocate Sorter
	if( ! sorter ) sorter = new Sorter(N);
	
	// This Is A Test Part. We Generate Pseudo Particles To Measure Maximum Particle Density
	FLOAT h = DENSITY/N;
	FOR_EVERY_CELL(10) {
		particle *p = new particle;
		p->p[0] = (i+0.5)*h;
		p->p[1] = (j+0.5)*h;
		p->type = FLUID;
		p->m = 1.0;
		particles.push_back(p);
	} END_FOR
	sorter->sort(particles);

	max_dens = 1.0;
	computeDensity();
	max_dens = 0.0;
	for( int n=0; n<particles.size(); n++ ) {
		particle *p = particles[n];
		max_dens = fmax(max_dens,p->dens);
		delete p;
	}
	// printf( "max_dens = %f\n", max_dens );
	particles.clear();
	
	FLOAT zeros[2] = { 0.0, 0.0 };
#if TEST
    for( FLOAT x=0.0; x < 1.0; x += DENSITY/N ) {
		for( FLOAT y=0.0; y < 1.0; y += DENSITY/N ) {
            if( hypot(x-0.25,y-0.5) < 0.15 ) {
				particle *p = new particle;
				FLOAT pos[2] = { x, y };
				pushParticle( pos, zeros, FLUID );
			}
        }
    }
#else

	double w = DENSITY*WALL_THICK;
	for( int i=0; i < N/DENSITY; i++ ) {
		for( int j=0; j < N/DENSITY; j++ ) {
			double x = i*w+w/2.0;
			double y = j*w+w/2.0;
			
			if( x > WALL_THICK && y > WALL_THICK ) {
				bool inside = particleSetting(x,y) == FLUID;
				if( inside ) {
					FLOAT pos[2] = { (FLOAT)x, (FLOAT)y };
					pushParticle( pos, zeros, FLUID );
				}
			}
		}
	}
	
	w = 1.0/N;
#if ! SPH
	for( int i=0; i < N; i++ ) {
		for( int j=0; j < N; j++ ) {
			double x = i*w+w/2.0;
			double y = j*w+w/2.0;
#else
	for( int i=0; i < N/DENSITY; i++ ) {
		for( int j=0; j < N/DENSITY; j++ ) {
			double x = (i*w+w/2.0)*DENSITY;
			double y = (j*w+w/2.0)*DENSITY;
#endif
			if( x < WALL_THICK || x > 1.0-WALL_THICK ||
               y < WALL_THICK || y > CEILING-WALL_THICK || particleSetting(x,y) == WALL ) {
				FLOAT pos[2] = { (FLOAT)x, (FLOAT)y };
				pushParticle( pos, zeros, WALL );
			}
		}
	}

#endif
    
	// Remove Particles That Stuck On Wal Cells
	sorter->sort(particles);
	sorter->markWater(A,DENSITY);
	
	for( vector<particle *>::iterator iter=particles.begin(); iter!=particles.end(); ) {
		particle &p = **iter;
		if( p.type == WALL ) {
			iter++;
			continue;
		}
		int i = fmin(N-1,fmax(0,p.p[0]*N));
		int j = fmin(N-1,fmax(0,p.p[1]*N));
		if( A[i][j] == WALL ) {
			delete *iter;
			iter = particles.erase(iter);
		} else {
			iter ++;
		}
	}
	
	// Compute Normal for Walls
	compute_wall_normal();
	
	// Initial Merge
	if( adaptive_sampling ) adaptive::initial_merge( sorter, A, particles, DENSITY );
	
	// Turn On Blending
	glEnable(GL_BLEND);
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	
	system( "rm -rf screenshot/*" );
}

static void pourWater( int limit ) {
    if( timeStep > limit ) return;
    FLOAT numr[2] = { 0.3, 0.45 };
	double w = DENSITY/N;
    for( double x=w+w/2.0; x < 1.0-w/2.0; x += w ) {
		if( x >= numr[0] && x <= numr[1] ) {
            FLOAT pos[2] = { (FLOAT)x, 1.0 - WALL_THICK - 2.5*DENSITY/N };
			FLOAT u[2] = { 0.0, -0.5*DENSITY/N/DT };
			pushParticle( pos, u, FLUID );
        }
    }
}

void flip2d::reshape( int w, int h ) {
	FLOAT margin = 0.5/N;
	glViewport(0, 0, w, h);
	glLoadIdentity();
	glOrtho(-margin,1.0+margin,-margin,1.0+margin,-1.0,1.0);
}

void raw_drawBitmapString( const char *string) {
	while (*string) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *string++);
}

static void calcMarchingPoints( int i, int j, int gn, FLOAT **L, FLOAT p[8][2], int &pnum ) {
	pnum = 0;
	FLOAT w = 1.0/gn;
	int quads[][2] = { {i, j}, {i+1, j}, {i+1, j+1}, {i, j+1} };
	for( int n=0; n<4; n++ ) {
		// Inside Liquid
		if( L[quads[n][0]][quads[n][1]] < 0.0 ) {
			p[pnum][0] = w*(0.5+quads[n][0]);
			p[pnum][1] = w*(0.5+quads[n][1]);
			pnum ++;
		}
		// If Line Crossed
		if( L[quads[n][0]][quads[n][1]] * L[quads[(n+1)%4][0]][quads[(n+1)%4][1]] < 0 ) {
			// Calculate Cross Position
			FLOAT y0 = L[quads[n][0]][quads[n][1]];
			FLOAT y1 = L[quads[(n+1)%4][0]][quads[(n+1)%4][1]];
			FLOAT a = y0/(y0-y1);
			FLOAT p0[2] = { (FLOAT)(w*(0.5+quads[n][0])), (FLOAT)(w*(0.5+quads[n][1])) };
			FLOAT p1[2] = { (FLOAT)(w*(0.5+quads[(n+1)%4][0])), (FLOAT)(w*(0.5+quads[(n+1)%4][1])) };
			p[pnum][0] = (1.0-a)*p0[0]+a*p1[0];
			p[pnum][1] = (1.0-a)*p0[1]+a*p1[1];
			pnum ++;
		}
	}	
}

static void drawMarchingCube() {	
	int gn = N;
	
	// Paint Distance Field
	glColor4d(0.1,0.2,0.3,1.0);
	FOR_EVERY_CELL(gn-1) {
		FLOAT p[8][2];
		int pnum;
		calcMarchingPoints( i, j, gn, SL, p, pnum );
		FLOAT v[pnum][2];
		for( int m=0; m<pnum; m++ ) {
			v[m][0] = p[m][0];
			v[m][1] = p[m][1];
		}
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(2, GL_FLOAT, 0, v);
		glDrawArrays(GL_TRIANGLE_FAN, 0, pnum);
		glDisableClientState(GL_VERTEX_ARRAY);
	} END_FOR;
}

void drawBitmapString( const char *string, void *font=GLUT_BITMAP_HELVETICA_12) {
	while (*string) glutBitmapCharacter(font, *string++);
}

static void render() {
#if 1 && ! TEST // Paint Fluid Domain
	FOR_EVERY_CELL(N) {
		if( A[i][j] != AIR ) {
			double h = 1.0/N;
			double pos[2] = {i*h,j*h};
			if( A[i][j] == FLUID ) glColor4d(0.4,0.6,1.0,0.0);
			else if( A[i][j] == WALL ) glColor4d(1.0,0.7,0.4,0.3);
			FLOAT v[4][2];
			unsigned indx[6] = { 0, 1, 3, 1, 2, 3 };
			v[0][0] = pos[0];
			v[0][1] = pos[1];
			v[1][0] = pos[0]+h;
			v[1][1] = pos[1];
			v[2][0] = pos[0]+h;
			v[2][1] = pos[1]+h;
			v[3][0] = pos[0];
			v[3][1] = pos[1]+h;
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(2, GL_FLOAT, 0, v);
			glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, indx);
			glDisableClientState(GL_VERTEX_ARRAY);
		}
	} END_FOR;
#endif
	
#if 0 // Paint Color
	FOR_EVERY_CELL(N) {
		double h = 1.0/N;
		double p[2] = {i*h,j*h};
		double c = -color_div[i][j];
		glColor4d(1.0,1.0,0.0,c*0.001);
		glBegin(GL_QUADS);
		glVertex2d(p[0],p[1]);
		glVertex2d(p[0]+h,p[1]);
		glVertex2d(p[0]+h,p[1]+h);
		glVertex2d(p[0],p[1]+h);
		glEnd();
	} END_FOR;
#endif
	
#if 1 // Paint What?
	drawMarchingCube();
#endif
	
#if 0 // Paint Filled Cells
	FOR_EVERY_CELL(N) {
		if( DeepArea[i][j] ) {
			double h = 1.0/N;
			double pos[2] = {i*h,j*h};
			glColor4d(0.8,0.3,0.3,0.8);
			glBegin(GL_QUADS);
			glVertex2d(pos[0],pos[1]);
			glVertex2d(pos[0]+h,pos[1]);
			glVertex2d(pos[0]+h,pos[1]+h);
			glVertex2d(pos[0],pos[1]+h);
			glEnd();
		}
	} END_FOR;
#endif
	
#if 0 // Paint What?
	FOR_EVERY_CELL(N) {
		if( DeepZone[i][j] ) {
			double h = 1.0/N;
			double pos[2] = {i*h,j*h};
			glColor4d(0.8,0.6,0.3,0.8);
			glBegin(GL_QUADS);
			glVertex2d(pos[0],pos[1]);
			glVertex2d(pos[0]+h,pos[1]);
			glVertex2d(pos[0]+h,pos[1]+h);
			glVertex2d(pos[0],pos[1]+h);
			glEnd();
		}
	} END_FOR;
#endif
	
#if 0 // Paint Liquid Boundary Dots
    glPointSize(2.0);
    glBegin(GL_POINTS);
	glColor4d(1.0,0.0,1.0,1.0);
    // Draw X Boundary Condition
	FOR_EVERY_X_FLOW(N) {
        if( ! F[0][i][j].empty() ) {
            double h = 1.0/N;
            double p[2] = {i*h,j*h+h/2.0};
            glVertex2d(p[0],p[1]);
        }
	} END_FOR
    // Draw Y Boundary Condition
	FOR_EVERY_Y_FLOW(N) {
        if( ! F[1][i][j].empty() ) {
            double h = 1.0/N;
            double p[2] = {i*h+h/2.0,j*h};
            glVertex2d(p[0],p[1]);
        }
	} END_FOR
    glEnd();
    glPointSize(1.0);
#endif
    
#if 0 // Paint Normals
	glBegin(GL_LINES);	
	for( int n=0; n<particles.size(); n++ ) {
		if( particles[n]->type == FLUID ) {
			//FLOAT k = particles[n]->obb.ratio;
			FLOAT s = -0.001*particles[n]->curv;
			glVertex2d(particles[n]->p[0],particles[n]->p[1]);
			glVertex2d(particles[n]->p[0]+s*particles[n]->n[0],
					   particles[n]->p[1]+s*particles[n]->n[1]);
		}
	}
	glEnd();
#endif
	
	FLOAT s = 10.0;
#if 0 // Paint Velocity  Projected Inside The Walls
	glColor4d(1.0,1.0,0.0,0.8);
	FOR_EVERY_CELL(N) {
		if( A[i][j] == WALL ) {
			FLOAT h = 1.0/N;
			FLOAT p[2] = {(FLOAT)(i*h+h/2.0),(FLOAT)(j*h+h/2.0)};
			FLOAT v[2] = {(FLOAT)(0.5*u[0][i][j]+0.5*u[0][i+1][j]),(FLOAT)(0.5*u[1][i][j]+0.5*u[1][i][j+1])};
			glBegin(GL_LINES);
			glVertex2d(p[0],p[1]);
			glVertex2d(p[0]+s*DT*v[0],p[1]+s*DT*v[1]);
			glEnd();
		}
	} END_FOR
#endif
	
#if 0 // Draw X flow
	glColor4d(0.0,0.0,1.0,1.0);
	FOR_EVERY_X_FLOW(N) {
		double h = 1.0/N;
		double p[2] = {i*h,j*h+h/2.0};
		glBegin(GL_LINES);
		glVertex2d(p[0],p[1]);
		glVertex2d(p[0]+s*DT*u[0][i][j],p[1]);
		glEnd();
	} END_FOR
	
	// Draw Y Flow
	glColor4d(1.0,0.0,0.0,1.0);
	FOR_EVERY_Y_FLOW(N) {
		double h = 1.0/N;
		double p[2] = {i*h+h/2.0,j*h};
		glBegin(GL_LINES);
		glVertex2d(p[0],p[1]);
		glVertex2d(p[0],p[1]+s*DT*u[1][i][j]);
		glEnd();
	} END_FOR
#endif
	
	// Paint Particles
#if 1 // Position
	glPointSize(2.0);
	
	FLOAT v[particles.size()][2];
	FLOAT c[particles.size()][4];

	for( int n=0; n<particles.size(); n++ ) {
		FLOAT a = particles[n]->dens;
		a = 0.8;
		if( particles[n]->type == FLUID ) {
			if( particles[n]->mark == 0 ) {
				if( particles[n]->level > 1 ) { c[n][0] = 0.0; c[n][1] = 1.0; c[n][2] = 0.3; c[n][3] = a; }
				else { c[n][0] = 0.5; c[n][1] = 0.7; c[n][2] = 1.0; c[n][3] = a; }
			} else {
				c[n][0] = 1.0; c[n][1] = 0.7; c[n][2] = 0.5; c[n][3] = a;
			}
		} else {
			c[n][0] = 0.7; c[n][1] = 0.7; c[n][2] = 0.0; c[n][3] = 1.0;
		}
		v[n][0] = particles[n]->p[0];
		v[n][1] = particles[n]->p[1];
	}
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	glVertexPointer(2, GL_FLOAT, 0, v);
	glColorPointer(4, GL_FLOAT, 0, c);
	glDrawArrays(GL_POINTS, 0, particles.size());
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	glPointSize(1.0);
	glColor4d(1.0,1.0,1.0,1.0);
	
#if 0 // Velocity
	glBegin(GL_LINES);	
	for( int n=0; n<particles.size(); n++ ) {
		glVertex2d(particles[n]->p[0],particles[n]->p[1]);
		glVertex2d(particles[n]->p[0]+DT*s*particles[n]->u[0],
				   particles[n]->p[1]+DT*s*particles[n]->u[1]);
	}
	glEnd();
#endif
	
#if 1 // Normal
	/*
	glBegin(GL_LINES);	
	for( int n=0; n<particles.size(); n++ ) {
		glVertex2d(particles[n]->p[0],particles[n]->p[1]);
		glVertex2d(particles[n]->p[0]+0.01*particles[n]->n[0],
				   particles[n]->p[1]+0.01*particles[n]->n[1]);
	}
	glEnd();
	*/
	glColor4d(1.0,1.0,1.0,0.8);
	FLOAT v1[N][N][2][2];
	int num_lines = 0;
	FOR_EVERY_CELL(N) {
		if(A[i][j]==WALL) {
			FLOAT h = 1.0/N;
			FLOAT p[2] = { (FLOAT)(i*h+h/2.0),(FLOAT)(j*h+h/2.0)};
			FLOAT v[2] = { (FLOAT)(wall_normal[0][i][j]), (FLOAT)(wall_normal[1][i][j]) };
			int line_i = num_lines/N;
			int line_j = num_lines%N;
			v1[line_i][line_j][0][0] = p[0];
			v1[line_i][line_j][0][1] = p[1];
			v1[line_i][line_j][1][0] = p[0]+0.5*h*v[0];
			v1[line_i][line_j][1][1] = p[1]+0.5*h*v[1];
			num_lines++;
		}
	} END_FOR
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(2, GL_FLOAT, 0, v1);
	glDrawArrays(GL_LINES, 0, num_lines*2);
	glDisableClientState(GL_VERTEX_ARRAY);

#endif
	
    if( drawOval && ! SPH ) {
        // OBB
        for( int n=0; n<particles.size(); n++ ) {
            if( particles[n]->type == FLUID ) {
                FLOAT a = particles[n]->level > 1 ? 0.8 : 0.4;
                if( particles[n]->mark == 0 ) {
					if( particles[n]->level > 1 ) glColor4d(0.0,1.0,0.3,a);
					else glColor4d(1.0,1.0,1.0,a);
				} else glColor4d(1.0,0.0,0.0,a);
                FLOAT dtheta = 20;
				FLOAT v[(int)dtheta][2];
				
                for( int i=0; i<dtheta; i++ ) {
                    FLOAT theta = 2*PI*i/dtheta;
                    OBB obb = particles[n]->obb;
                    for ( int i=0; i<2; i++ ) obb.c[i] = 2.0/obb.c[i];
					int level = particles[n]->level;
					FLOAT r = comp_rad(level,DENSITY/N);
                    FLOAT pos[] = { (FLOAT)(cos(theta)*r+particles[n]->p[0]), (FLOAT)(sin(theta)*r+particles[n]->p[1]) };
                    implicit::stretchPosition(obb,pos,particles[n]->p);
					v[i][0] = pos[0];
					v[i][1] = pos[1];
                }

				glEnableClientState(GL_VERTEX_ARRAY);
				glVertexPointer(2, GL_FLOAT, 0, v);
				glDrawArrays(GL_LINE_LOOP, 0, dtheta);
				glDisableClientState(GL_VERTEX_ARRAY);
            }
        }
    }
    glColor4d(1.0,1.0,1.0,1.0);
	
#endif
    
#if TEST
    glRasterPos2d(0.04, 0.95);
    char tmp[64];
    sprintf( tmp, "TimeStep=%d", timeStep );
	raw_drawBitmapString(tmp);
#endif
	
	// Display Usage
	if( ! hide_help ) {
		int winsize[2] = { glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT) };
		double dw = 1.0/(double)winsize[0];
		double dh = 1.0/(double)winsize[1];
		int peny = 20;
		glColor4d( 1.0, 1.0, 1.0, 1.0 );
		for( int j=0; j<10; j++ ) {
			glRasterPos2d(25*dw, 1.0-(j+1)*(peny)*dh-10*dh);
			switch(j) {
				case 0:
					drawBitmapString("Push \"c\" to toggle the spring force correction");
					if( correction ) drawBitmapString( " ( current: Enabled ).");
					else drawBitmapString( " ( current: Disabled ).");
					break;
				case 1:
					drawBitmapString("Push \"s\" to toggle the particle split");
					if( splitting ) drawBitmapString( " ( current: Enabled ).");
					else drawBitmapString( " ( current: Disabled ).");
					break;
				case 2:
					drawBitmapString("Push \"v\" to toggle order of boundary accuracy.");
					if( subcell ) drawBitmapString( " ( current: 2nd order ).");
					else drawBitmapString( " ( current: 1st order ).");
					break;
				case 3:
					drawBitmapString("Push \"o\" to toggle anisotropic view.");
					if( drawOval ) drawBitmapString( " ( current: Shown ).");
					else drawBitmapString( " ( current: Hiden ).");
					break;
				case 4:
					drawBitmapString("Push \"z\" to toggle adaptive sampling.");
					if( adaptive_sampling ) drawBitmapString( " ( current: Enabled ).");
					else drawBitmapString( " ( current: Disabled ).");
					break;
				case 5:
					drawBitmapString("Push \"x\" to toggle volume correction.");
					if( correct_volume ) drawBitmapString( " ( current: Enabled ).");
					else drawBitmapString( " ( current: Disabled ).");
					break;
				case 6:
					drawBitmapString("Push \"h\" to toggle help display.");
					break;
				case 7:
					drawBitmapString("Push \"r\" to reset.");
					break;
			}
		}
	} else {
		glRasterPos2d(0.05,0.93);
		char timeStr[64];
		static float fps = 0;
		float frame = 1.0f/dumptime();
		if (!isinf(frame))
			fps = 0.1f*frame + 0.9f*fps;
		sprintf( timeStr, "FPS: %.2f", fps );
		drawBitmapString(timeStr,GLUT_BITMAP_HELVETICA_18);
	}
}

static void write_svg() {	
	int width = glutGet(GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);
	
	char path[256];
	static int counter = 0;
	
	sprintf( path, "screenshot/frame_%d.svg", counter++ );
	
	FILE *fp;
	if( ! (fp = fopen( path, "w" ))) exit(1);
	
	// Write SVG Header
	fprintf( fp, "<?xml version=\"1.0\" standalone=\"no\"?>\n" );
	fprintf( fp, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n" );
	fprintf( fp, "<svg width=\"30cm\" height=\"30cm\" viewBox=\"0 0 %d %d\"\n", (int)width, (int)width );
	fprintf( fp, "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n" );
	
	const char *colors[] = { "mediumslateblue", "orange", "orangered", "mediumseagreen" };
	
	// Write Particles
	for( int t=0; t<2; t++ ) {
		char type = t == 0 ? WALL : FLUID;
		for( int n=0; n<particles.size(); n++ ) {	
			particle &p = *particles[n];
			if( p.type != type ) continue;
			
			if( p.p[1] > CEILING ) continue;
			double r = p.type == FLUID ? 0.7*0.5*DENSITY*width/(double)N : 0.5*width/(double)N;
			
			// Compute Rotation
			double theta = acos(p.obb.u[0][0])/PI*180.0;
			if( p.obb.u[0][1] < 0.0 ) theta += 2.0*(180.0 - theta);
			double rad = comp_rad(p.level,1.0);
			
			if( p.type == FLUID ) {
				fprintf( fp, "<ellipse transform=\"translate(%lf %lf) rotate(%lf)\" rx=\"%lf\" ry=\"%lf\" "
						"fill=\"%s\" stroke=\"none\" />\n", 
						width*p.p[0], width*(1.0 - p.p[1]), -theta, rad*r*fmax(0.5,p.obb.c[0]), rad*r*fmax(0.5,p.obb.c[1]),
						(p.mark > 0 ? colors[2] : (p.level>1 ? colors[3] : colors[0]) ));
			} else {
				fprintf( fp,  "<rect x=\"%lf\" y=\"%lf\" width=\"%lf\" height=\"%lf\" fill=\"%s\" stroke=\"black\" stroke-width=\"2\" />\n",
							width*p.p[0]-r, width*(1.0-p.p[1])-r, 2*r, 2*r, colors[1] );
			}
		}
	}
	
	// Write SVG Footer
	fprintf( fp, "</svg>\n" );
	fclose(fp);
}

static bool write_frame() {
	int width = glutGet(GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);
	unsigned char *buffer = new unsigned char[width*height*4];
	
	glFlush();
	glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
	glReadPixels( 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, buffer );
	
	char name[256];
	static int counter = 0;
    
    // ffmpeg -qscale 4 -r 60 -b 9600 -i frame_%d.bmp out.mp4
	sprintf( name, "screenshot/frame_%d.bmp", counter++ );
	write_bmp( name, buffer, width, height*CEILING, true );
	
	delete [] buffer;
	return true;
}

static void save_grid() {
	FOR_EVERY_X_FLOW(N) {
		u_save[0][i][j] = u[0][i][j];
	} END_FOR
	
	FOR_EVERY_Y_FLOW(N) {
		u_save[1][i][j] = u[1][i][j];
	} END_FOR
}

static void subtract_grid() {
	FOR_EVERY_X_FLOW(N) {
		u_save[0][i][j] = u[0][i][j] - u_save[0][i][j];
	} END_FOR
	
	FOR_EVERY_Y_FLOW(N) {
		u_save[1][i][j] = u[1][i][j] - u_save[1][i][j];
	} END_FOR
}

static void comp_divergence() {
	FLOAT h = 1.0/N;
	FOR_EVERY_CELL(N) {
		FLOAT div = A[i][j] == FLUID ? (u[0][i+1][j]-u[0][i][j]) + (u[1][i][j+1]-u[1][i][j]) : 0.0;
		d[i][j] = div/h;
	} END_FOR
}

static inline char wall_chk( char a ) {
	return a == WALL ? 1.0 : -1.0;
}

static void enforce_boundary() {
	FOR_EVERY_X_FLOW(N) {
		if( i==0 || i==N ) u[0][i][j] = 0.0;
		if( i<N && i>0 && wall_chk(A[i][j])*wall_chk(A[i-1][j]) < 0 ) {
			u[0][i][j] = 0.0;
		}
	} END_FOR
	
	FOR_EVERY_Y_FLOW(N) {
		if( j==0 || j==N ) u[1][i][j] = 0.0;
		if( j<N && j>0 && wall_chk(A[i][j])*wall_chk(A[i][j-1]) < 0 ) {
			u[1][i][j] = 0.0;
		}
	} END_FOR
}

static void buildBoundaryConditions ( Vector2 ***F, FLOAT **d ) {
    FOR_EVERY_X_FLOW(N) {
        F[0][i][j][0] = F[0][i][j][1] = 0.0;
    } END_FOR
    
    FOR_EVERY_Y_FLOW(N) {
        F[1][i][j][0] = F[1][i][j][1] = 0.0;
    } END_FOR
    
    FLOAT h = 1.0/N;
    FOR_EVERY_CELL(N) {
        if( i>0 && i<N-1 && j>0 && j<N-1 ) {
            if( A[i][j] == FLUID ) {
                for( int dim=0; dim<DIM; dim++ ) {      // Dimension
                    for( int dir=-1; dir<=1; dir+=2 ) { // Direction
                        
                        // Cell Center Position
                        int ni = i+(dim==0)*dir;
                        int nj = j+(dim==1)*dir;
                        
                        // Staggered Position
                        int si = i+(dim==0)*(dir==1);
                        int sj = j+(dim==1)*(dir==1);
                        
                        // If Neighbor Is An Air Cell
                        if( A[ni][nj] == AIR ) {
                            F[dim][si][sj][0] = L[ni][nj]*L[i][j] < 0.0 ? L[ni][nj]/fmin(1.0e-6,L[i][j]) : 0.0;
                        }
                        
                        // If Neighbor Is A Wall Cell
                        if( A[ni][nj] == WALL ) {
                            FLOAT dotProduct = 0.0;
                            FLOAT cu[DIM];
                            for( int dm=0; dm<DIM; dm++ ) {
                                cu[dm] = dim==dm ? u[dm][si][sj] : 0.5*(u[dm][i][j]+u[dm][i+(dm==0)][j+(dm==1)]);
                            }
                            
                            FLOAT wn[2] = { wall_normal[0][ni][nj], wall_normal[1][ni][nj] };
                            //FLOAT wn[2] = { -dir*(dim==0), -dir*(dim==1) };
                            for( int dm=0; dm<DIM; dm++ ) {
                                dotProduct += wn[dm]*cu[dm];
                            }
                            F[dim][si][sj][0] = 1.0;
                            F[dim][si][sj][1] = -fabs(wn[dim])*dotProduct*h;
                            d[i][j] += -F[dim][si][sj][1]/(h*h);
                        }
                    }
                }
            }
        }
    } END_FOR
}

static FLOAT computeAccurateVolume( FLOAT **L ) {
	FLOAT volume = 0.0;
	FOR_EVERY_CELL(N-1) {
		FLOAT p[8][2];
		int pnum;
		calcMarchingPoints( i, j, N, L, p, pnum );
		for( int m=0; m<pnum; m++ ) {
			volume += p[m][0]*p[(m+1)%pnum][1]-p[m][1]*p[(m+1)%pnum][0];
		}
	} END_FOR;
	
	return volume*0.5;
}

static FLOAT computeVolume( FLOAT **L ) {
	FLOAT volume = 0.0;
	FLOAT w = 1.0/N;
	FOR_EVERY_CELL(N-1) {
		volume += (A[i][j]==FLUID)*w*w;
	} END_FOR;
	return volume;
}

static FLOAT computeVolumeError( FLOAT volume0, FLOAT **L ) {
	FLOAT vdiv = 0.0;
	FLOAT curVolume = computeVolume(L);
	
	FLOAT x = (curVolume - volume0)/volume0;
	y_volume0 += x*DT;
	
	FLOAT kp = 2.3 / (25.0 * DT);
	FLOAT ki = kp*kp/16.0;
	vdiv = -(kp * x + ki * y_volume0) / (x + 1.0);
	return vdiv;
}

static void compute_pressure() {
	// Clear Pressure
	FOR_EVERY_CELL(N) {
		p[i][j] = 0.0;
	} END_FOR
	
	// Calculate LevelSet
	FOR_EVERY_CELL(N) {
		L[i][j] = sorter->levelset(i,j,DENSITY);
	} END_FOR;
    
	// Compute Initial Volume if Neccessary
	if( ! volume0 ) volume0 = computeVolume(L);
#if 0 // Save Data
	FILE *fp = fopen("/Users/mscp/Desktop/volume.dat","a");
	fprintf( fp, "%d %f\n", timeStep, computeVolume(L)/volume0 );
	fclose(fp);
#endif
	if( correct_volume ) {
		FLOAT vdiv = computeVolumeError( volume0, L );
		FOR_EVERY_CELL(N) {
			if( A[i][j] == FLUID ) d[i][j] -= vdiv;
		} END_FOR;
	}
	
    // Compute Boundary Condition
    buildBoundaryConditions(F,d);
	
	// Solve Ap = d ( p = Pressure, d = Divergence )
	solver::setSubcell(subcell);
	solver::solve( A, F, p, d, N );
	
	// Instant Wall Separating Condition
    FOR_EVERY_CELL(N) {
        if( j >= N-2 ) {
            if( p[i][j] < p[i][j-1] ) p[i][j] = p[i][j-1];
        }
    } END_FOR;
	
	FOR_EVERY_CELL(N) {
		if( j >= N-2 ) p[i][j] = p[i][j+1];
	} END_FOR
}

static void extrapolate_velocity() {
	// Mark Fluid Cell Face
	static char **mark[2] = { alloc2D<char>(N+1), alloc2D<char>(N+1) };
	static char **wall_mark[2] = { alloc2D<char>(N+1), alloc2D<char>(N+1) };
	static FLOAT **tmp_u[2] = { alloc2D<FLOAT>(N+1), alloc2D<FLOAT>(N+1) };
	
	OPENMP_FOR FOR_EVERY_X_FLOW(N) {
		mark[0][i][j] = (i>0 && A[i-1][j]==FLUID) || (i<N && A[i][j]==FLUID);
		wall_mark[0][i][j] = (i<=0 || A[i-1][j]==WALL) && (i>=N || A[i][j]==WALL);
	} END_FOR
	
	OPENMP_FOR FOR_EVERY_Y_FLOW(N) {
		mark[1][i][j] = (j>0 && A[i][j-1]==FLUID) || (j<N && A[i][j]==FLUID);
		wall_mark[1][i][j] = (j<=0 || A[i][j-1]==WALL) && (j>=N || A[i][j]==WALL);
	} END_FOR
	
	// Now Extrapolate
	FOR_EVERY_CELL(N+1) {
		for( int n=0; n<DIM; n++ ) {
			if( n!=0 && i>N-1 ) continue;
			if( n!=1 && j>N-1 ) continue;
			
			tmp_u[n][i][j] = u[n][i][j];
			if( ! mark[n][i][j] && wall_mark[n][i][j] ) {
				int wsum = 0;
				FLOAT sum = 0.0;
				int q[][2] = { {i-1,j}, {i+1,j}, {i,j-1}, {i,j+1} };
				for( int qk=0; qk<4; qk++ ) {
					if( q[qk][0]>=0 && q[qk][0]<N+(n==0) && q[qk][1]>=0 && q[qk][1]<N+(n==1) ) {
						if( mark[n][q[qk][0]][q[qk][1]] ) {
							wsum ++;
							sum += u[n][q[qk][0]][q[qk][1]];
						}
					}
				}
				if( wsum ) tmp_u[n][i][j] = sum/wsum;
			}
		}
	} END_FOR;
	
	OPENMP_FOR FOR_EVERY_CELL(N+1) {
		for( int n=0; n<DIM; n++ ) {
			if( n!=0 && i>N-1 ) continue;
			if( n!=1 && j>N-1 ) continue;
			u[n][i][j] = tmp_u[n][i][j];
		}
	} END_FOR;
}

static void subtract_pressure() {
	FLOAT h = 1.0/N;
	FOR_EVERY_X_FLOW(N) {
		if( i>0 && i<N ) {
			FLOAT pf = p[i][j];
			FLOAT pb = p[i-1][j];
			if( subcell && !F[0][i][j].empty() ) {
				pf = A[i-1][j] == FLUID ? F[0][i][j][0]*pb+F[0][i][j][1] : pf;
				pb = A[i][j] == FLUID ? F[0][i][j][0]*pf+F[0][i][j][1] : pb;
			} else {
				if( A[i-1][j]==WALL ) pb = pf;
				if( A[i][j]==WALL ) pf = pb;
			}
			u[0][i][j] -= (pf-pb)/h;
		}
	} END_FOR;
	
	FOR_EVERY_Y_FLOW(N) {
		if( j>0 && j<N ) {
			FLOAT pf = p[i][j];
			FLOAT pb = p[i][j-1];
			if( subcell && !F[1][i][j].empty() ) {
				pf = A[i][j-1] == FLUID ? F[1][i][j][0]*pb+F[1][i][j][1] : pf;
				pb = A[i][j] == FLUID ? F[1][i][j][0]*pf+F[1][i][j][1] : pb;
			} else {
				if( A[i][j-1]==WALL ) pb = pf;
				if( A[i][j]==WALL ) pf = pb;
			}
			u[1][i][j] -= (pf-pb)/h;
		}
	} END_FOR;
}

static void add_ExtForce() {
	for( int n=0; n<particles.size(); n++ ) {
		// Add Gravity
		//particles[n]->f[0] += -DT*GRAVITY;
		particles[n]->f[1] += -DT*GRAVITY;
	}
	
	for( int n=0; n<particles.size(); n++ ) {
		// Add External Force
		if( particles[n]->type == FLUID ) {
			particles[n]->u[0] += particles[n]->f[0];
			particles[n]->u[1] += particles[n]->f[1];
			particles[n]->f[0] = particles[n]->f[1] = 0.0;
		} else {
			particles[n]->u[0] = particles[n]->u[1] = 0.0;
		}
	}
}

static void advect_particle() {
	
	// Advect Particle by RK2
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		if( particles[n]->type == FLUID ) {
			FLOAT p1[2];
			FLOAT vel[2];
			mapper::fetchVelocity( particles[n]->p, vel, u, N );
#if 0 // ?
			for( int k=0; k<2; k++ ) p1[k] = particles[n]->p[k]+DT*vel[k];	
			int i = fmax(0,fmin(N-1,N*p1[0]));
			int j = fmax(0,fmin(N-1,N*p1[1]));
			if( A[i][j] == FLUID ) {
				FLOAT u1[2];
				mapper::fetchVelocity( p1, u1, u, N );
				for( int k=0; k<2; k++ ) vel[k] = 0.5*(vel[k]+u1[k]);
			}
#endif
			for( int k=0; k<2; k++ ) {
				particles[n]->p[k] += DT*vel[k];
			}
		}
	}
	
	// Sort
	sorter->sort(particles);
	
	// Constraint Outer Wall
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		FLOAT r = 1.01*WALL_THICK;
		for( int k=0; k<2; k++ ) {
			if( particles[n]->type == FLUID ) {
				particles[n]->p[k] = fmax(r,fmin(1.0-r,particles[n]->p[k]));
			}
		}
#if 1 // What?
		particle *p = particles[n];
		if( p->type == FLUID ) {
			int i = fmin(N-1,fmax(0,p->p[0]*N));
			int j = fmin(N-1,fmax(0,p->p[1]*N));
            
			vector<particle *> neighbors = sorter->getNeigboringParticles_cell(i,j,1,1);
			for( int n=0; n<neighbors.size(); n++ ) {
				particle *np = neighbors[n];
				double re = 0.5/N + 0.5*comp_rad(np->level,DENSITY/N);
				if( np->type == WALL ) {
					double dist = hypotf(p->p[0]-np->p[0],p->p[1]-np->p[1]);					
					if( dist < re ) {	
						
						FLOAT normal[2] = { np->n[0], np->n[1] };
						if( normal[0] == 0.0 && normal[1] == 0.0 && dist ) {
							for( int c=0; c<2; c++ ) normal[c] = (p->p[c]-np->p[c])/dist;
						}
						p->p[0] += (re-dist)*normal[0];
						p->p[1] += (re-dist)*normal[1];
						FLOAT dot = p->u[0] * normal[0] + p->u[1] * normal[1];
						p->u[0] -= dot*normal[0];
						p->u[1] -= dot*normal[1];
					}
				}
			}
			
			// Force Push...
			int max_exit = 500;
			while(1) {
				int i = fmin(N-1,fmax(0,p->p[0]*N));
				int j = fmin(N-1,fmax(0,p->p[1]*N));
				if( A[i][j] != WALL ) {
					break;
				} else {
					FLOAT normal[2] = { wall_normal[0][i][j], wall_normal[1][i][j] };
					p->p[0] += 0.1/N*normal[0];
					p->p[1] += 0.1/N*normal[1];
					max_exit --;
				}
				if( ! max_exit ) {
					printf( "Faild to push out particles out of walls...\n" );
					exit(0);
				}
			}
		}
#endif
	}
}

static FLOAT streamfunc( FLOAT x, FLOAT y ) {
    return 1.0/PI*pow(sin(PI*x),2)*pow(sin(PI*y),2);
}

static void streamflow( int inverse, FLOAT p[2], FLOAT out[2] ) {
    FLOAT h = 0.0001;
    FLOAT sgn = inverse ? -1.0 : 1.0;
    out[0] = sgn*(streamfunc(p[0],p[1]+h)-streamfunc(p[0],p[1]-h))/(2.0*h);
    out[1] = -sgn*(streamfunc(p[0]+h,p[1])-streamfunc(p[0]-h,p[1]))/(2.0*h);
}

static void rungekutta4_streamflow( FLOAT p[2], FLOAT out[2], FLOAT dt=DT, int inverse=0 ) {
    const int dim = 2;
    FLOAT k[4][dim];
    FLOAT tmp[dim];
    streamflow(inverse,p,k[0]);    // k[0]
    
    for( int i=0; i<dim; i++ ) tmp[i] = p[i]+0.5*dt*k[0][i];
    streamflow(inverse,tmp,k[1]);  // k[1]
    
    for( int i=0; i<dim; i++ ) tmp[i] = p[i]+0.5*dt*k[1][i];
    streamflow(inverse,tmp,k[2]);  // k[2]
    
    for( int i=0; i<dim; i++ ) tmp[i] = p[i]+dt*k[2][i];
    streamflow(inverse,tmp,k[3]);  // k[3]
    
    for( int i=0; i<dim; i++ ) out[i] = (k[0][i]+2.0*k[1][i]+2.0*k[2][i]+k[3][i])/6.0;
}

static void solve_picflip() {
    
    // Map Particles Onto Grid
	sorter->sort(particles);
	mapper::mapP2G(sorter,particles,u,N);
	sorter->markWater(A,DENSITY);
	
#if 1 // Compute Surface Tension Force
	tension::add_surface_tension( sorter, particles, DT, 0.0, 0.4, 0.03 );
#endif
    
    // Solve Velocity On Grid
	save_grid();
	if( ! subcell ) enforce_boundary();
	comp_divergence();
	compute_pressure();
	subtract_pressure();
	if( ! subcell ) enforce_boundary();
	extrapolate_velocity();
	subtract_grid();

#if FLIP
	// Copy Current Velocity
	for( int n=0; n<particles.size(); n++ ) {
		particles[n]->tmp[5][0] = particles[n]->u[0];
		particles[n]->tmp[5][1] = particles[n]->u[1];
	}
	
	// Map Changes Back To Particles
	mapper::mapG2P(particles,u_save,N);
    
	// Set Tmp As FLIP Velocity
	for( int n=0; n<particles.size(); n++ ) {
		particles[n]->tmp[5][0] = particles[n]->u[0] + particles[n]->tmp[5][0];
		particles[n]->tmp[5][1] = particles[n]->u[1] + particles[n]->tmp[5][1];
	}
	
	// Set u As PIC Velocity
	mapper::mapG2P(particles,u,N);
	
	// Interpolate
	for( int n=0; n<particles.size(); n++ ) {
		particles[n]->u[0] = (1.0-ALPHA)*particles[n]->u[0] + ALPHA*particles[n]->tmp[5][0];
		particles[n]->u[1] = (1.0-ALPHA)*particles[n]->u[1] + ALPHA*particles[n]->tmp[5][1];
	}
	
#else
	// Map Changes Back To Particles
	mapper::mapG2P(particles,u,N);
#endif
}

static void simulate() {
	// Reset Timer
	dumptime();
	
    // Pour Water
	pourWater(-1);
    
	// Sort Particles
	sorter->sort(particles);
	
	// Mark Cell
	sorter->markWater(A,DENSITY);
	
	// Adaptive Resampling
	DeepZone = adaptive::resample( sorter, A, particles, DENSITY, !adaptive_sampling );
	
	// Compute Average Position And Density
	computeDensity();
	
	// Insert Thin Particles
	if( splitting ) {
		sheet::keepThinSheet(A,sorter,particles,OBB_RATE,DENSITY,WALL_THICK);
	}
	sheet::collapseThinSheet(A,sorter,particles,OBB_RATE,DENSITY);
	buildOBB(A,DeepZone,sorter,particles);
    
#if TEST
    const int turn_cnt = 3000;
    for( int n=0; n<particles.size(); n++ ) {
        if( timeStep <= 2*turn_cnt+1 ) {
            rungekutta4_streamflow( particles[n]->p, particles[n]->u, DT, timeStep > turn_cnt);
			for( int k=0; k<2; k++ ) particles[n]->p[k] += DT*particles[n]->u[k];
        } else {
            particles[n]->u[0] = particles[n]->u[1] = 0.0;
        }
    }
#else
	
	// Add Gravity
	add_ExtForce();
	
	// Solve Fluid
	solve_picflip();
	
#endif // TEST
	
	// Generate Deep Area
	if( ! DeepArea ) DeepArea = alloc2D<char>(N);
	FOR_EVERY_CELL(N) {
		DeepArea[i][j] = 0;
		if( A[i][j] != WALL ) {
			if( DeepZone[i][j] ) {
				DeepArea[i][j] = 1;
			} else {
				vector<particle *> cell_particles = sorter->getNeigboringParticles_cell(i,j,1,1);
				int counts[2] = { 0, 0 };
				for( int n=0; n<cell_particles.size(); n++ ) {
					if( cell_particles[n]->level > 1 ) counts[1] ++;
					else counts[0]++;
				}
				DeepArea[i][j] = (counts[1]>0.5*counts[0]);
				if( ! DeepArea[i][j] ) {
					int q[][2] = { {i-1,j}, {i+1,j}, {i,j-1}, {i,j+1} };
					int fcnt = 0;
					for( int m=0; m<4; m++ ) {
						int qi = q[m][0];
						int qj = q[m][1];
						if( qi<0 || qi>N-1 || qj<0 || qj>N-1 || A[qi][qj]==WALL || A[qi][qj]==FLUID ) fcnt++;
					}
					if( fcnt == 4 ) DeepArea[i][j] = 1;
				}
			}
		}
	} END_FOR;
	
#if ! TEST
	// Advect Particle
	advect_particle();
	
	// Correct Position
	if( correction ) corrector::correct(sorter,particles,DT,DENSITY/N,anisotropic_spring);
#endif
		
	// Generate Surface LevelSet
	FOR_EVERY_CELL(N) {
		FLOAT p[] = { (FLOAT)((i+0.5)/N), (FLOAT)((j+0.5)/N) };
		SL[i][j] = implicit::implicit_func(sorter,DeepArea,p,DENSITY);
		if( A[i][j] == WALL ) SL[i][j] = 1.0;
	} END_FOR;
}

void flip2d::display() {
	
	// Reset If Requested
	if( doInit ) {
		flip2d::init();
	}
	
#if SPH
	sorter->sort(particles);
	sorter->markWater(A,DENSITY);
	add_ExtForce();
	
	static char **tmpDeep = alloc2D<char>(N);
	DeepZone = tmpDeep;
	FOR_EVERY_CELL(N) {
		DeepZone[i][j] = 0;
	} END_FOR;
	
	SPH2D::solve_SPH(sorter, particles, DENSITY/N, DT, doInit || timeStep==0 );
	
	sorter->sort(particles);
	sorter->markWater(A,DENSITY);
	
	// Insert Thin Particles
	if( splitting ) {
		computeDensity();
		sheet::keepThinSheet(A,sorter,particles,OBB_RATE,DENSITY,WALL_THICK);
	}
	sheet::collapseThinSheet(A,sorter,particles,OBB_RATE,DENSITY);
	buildOBB(A,DeepZone,sorter,particles);
#else
	simulate();
#endif
	// Rernder Result
	render();
	
	// Write Image
#if RECORDING
	if( timeStep++ % 3 == 0 ) {
		write_svg();
		write_frame();
	}
    //if( timeStep > 1600 ) exit(0);
#endif
	
	// Increment Time Step
    timeStep++;
	doInit = 0;
}

void flip2d::keyDown( unsigned char key ) {
	switch(key) {
		case 'r':
			doInit = 1;
			init_num ++;
			break;
		case 'c':
			correction = ! correction;
			break;
		case 's':
			splitting = ! splitting;
			break;
        case 'v':
            subcell = ! subcell;
            break;
		case 'x':
			correct_volume = ! correct_volume;
			break;
        case 'o':
            drawOval = ! drawOval;
            break;
		case 'z':
			adaptive_sampling = ! adaptive_sampling;
			break;
		case 'h':
			hide_help = ! hide_help;
			break;
		case '\e':
			exit(0);
			break;
	}
}

void flip2d::mouse( FLOAT x, FLOAT y, int state ) {
}

void flip2d::motion( FLOAT x, FLOAT y, FLOAT dx, FLOAT dy ) {
	int i = fmin(N-1,fmax(0,x*N));
	int j = fmin(N-1,fmax(0,y*N));
	FLOAT s = 0.125/DT;
	vector<particle *> neighbors = sorter->getNeigboringParticles_cell(i,j,1,1);
	for( int n=0; n<neighbors.size(); n++ ) {
		particle *p = neighbors[n];
		if( hypotf(x-p->p[0],y-p->p[1]) < 0.5 ) {
			p->f[0] = s*dx;
			p->f[1] = s*dy;
		}
	}
}












