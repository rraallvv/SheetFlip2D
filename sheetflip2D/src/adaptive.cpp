//
//  adaptive.cpp
//  flip2d
//

#include "adaptive.h"
#include "corrector.h"
#include "utility.h"

#define MAX_SPLIT	4

using namespace std;
static const int depeth = 1;
static char **Z = NULL;

static bool split( Sorter *sorter, char **A, char **Z, std::vector<particle *> &particles, FLOAT density, int gn ) {
	bool didSplit = false;
	vector<particle *> new_particles;
	
	for( int n=0; n<particles.size(); n++ ) {	
		particles[n]->remove = 0;
		if( particles[n]->type == FLUID ) {
			particle &p = *particles[n];
			if( p.level==1 ) continue;
			if( p.mark ) {
				printf( "This should not happen !\n" );
				exit(0);
				continue;
			}
			int i = fmax(0,fmin(gn-1,gn*p.p[0]));
			int j = fmax(0,fmin(gn-1,gn*p.p[1]));
			// If It Lies In Shallow Zone
			if( ! Z[i][j] && A[i][j] != WALL ) {
				FLOAT w = 0.25*density/gn;
				FLOAT theta = 2.0*PI*(rand()%1001)/1000.0;
				for( int m=0; m<p.level; m++ ) {
					particle *newp = new particle;
					newp->p[0] = p.p[0]+w*cos(theta+m*2*PI/p.level);
					newp->p[1] = p.p[1]+w*sin(theta+m*2*PI/p.level);
					newp->u[0] = p.u[0];
					newp->u[1] = p.u[1];
					newp->n[0] = p.n[0];
					newp->n[1] = p.n[1];
					newp->f[0] = p.f[0];
					newp->f[1] = p.f[1];
					newp->dens = p.dens;
					newp->curv = p.curv;
					newp->type = FLUID;
					newp->mark = 0;
					newp->m = 1.0;
					newp->level = 1;
					mass aMass = { newp->m, 0, newp };
					newp->hosts.push_back(aMass);
					new_particles.push_back(newp);
					newp->split_cnt = 0;
				}
				p.remove = 1;
				didSplit = true;
			}
		}
	}
	
	// Remove
	cleanParticles(sorter,particles);
	
	// Insert New Particles
	particles.insert( particles.end(), new_particles.begin(), new_particles.end() );
	
	if( new_particles.size() ) sorter->sort(particles);
	return didSplit;
}

static bool merge( Sorter *sorter, char **A, char **Z, std::vector<particle *> &particles, FLOAT density, int gn) {
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		particles[n]->pref = NULL;
		particles[n]->remove = 0;
		if( particles[n]->type == FLUID ) {
			particle &p = *particles[n];
			if( p.level>=MAX_SPLIT ) continue;
			if( p.mark ) continue;
			if( p.split_cnt > 1 ) continue;
			
			int i = fmax(0,fmin(gn-1,gn*p.p[0]));
			int j = fmax(0,fmin(gn-1,gn*p.p[1]));
			// If It Lies In Deep Zone
			if( Z[i][j] || A[i][j] == WALL ) {
				std::vector<particle *> neighbors = sorter->getNeigboringParticles_cell(i,j,1,1);
				FLOAT dist2 = 999.0;
				for( int m=0; m<neighbors.size(); m++ ) {
					particle &np = *neighbors[m];
					if( np.type == WALL ) continue;
					if( np.mark ) continue;
					if( np.split_cnt > 1 ) continue;
					if( np.level+p.level > MAX_SPLIT ) continue;
					
					int ni = fmax(0,fmin(gn-1,gn*np.p[0]));
					int nj = fmax(0,fmin(gn-1,gn*np.p[1]));
					if( ! Z[ni][nj] ) continue;
					
					FLOAT d2 = length2(np.p,p.p);
					if( np.type == FLUID && particles[n] != neighbors[m] && d2<1.0/(gn*gn) ) {
						if( d2 < dist2 ) {
							dist2 = d2;
							particles[n]->pref = neighbors[m];
						}
					}
				}
			}
		}
	}
	
	// Go Actually Merge
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		bool doMerge = (
						particles[n]->pref && 
						particles[n]->pref->pref == particles[n] && 
						particles[n] > particles[n]->pref );
		
		if( doMerge ) {
			particle &p = *particles[n];
			particle &pref = *p.pref;
			pref.remove = 1;
			// Slightly Move
			p.p[0] = 0.5*(p.p[0]+pref.p[0]);
			p.p[1] = 0.5*(p.p[1]+pref.p[1]);
			p.u[0] = 0.5*(p.u[0]+pref.u[0]);
			p.u[1] = 0.5*(p.u[1]+pref.u[1]);
			p.m = p.m+pref.m;
			p.hosts[0].m = p.m;
			p.level += pref.level;
		}
	}
	
	// Remove
	return cleanParticles(sorter,particles);;
}

// Deep Zone Will Be Marked With "1"
static void findDeepZone( char **A, char **Z, int depth, int gn ) {
	static char **Z_tmp = alloc2D<char>(gn);
	OPENMP_FOR FOR_EVERY_CELL(gn) {
		Z[i][j] = 1;
		if( A[i][j] == AIR ) {
			int cnt=0;
			for( int nx=-1; nx<=1; nx++ ) {
				for( int ny=-1; ny<=1; ny++ ) {
					int fi = i+nx;
					int fj = j+ny;
					if( fi>=0 && fi<gn && fj>=0 && fj<gn ) {
						if( A[fi][fj] == AIR ) cnt ++;
					}
				}
			}
			if( cnt > 3 ) Z[i][j] = 0;
		}
		Z_tmp[i][j] = Z[i][j];
	} END_FOR;
	
	for( int d=0; d<depth; d++ ) {
		OPENMP_FOR FOR_EVERY_CELL(gn) {
			for( int nx=-1; nx<=1; nx++ ) {
				for( int ny=-1; ny<=1; ny++ ) {
					int fi = i+nx;
					int fj = j+ny;
					if( fi>=0 && fi<gn && fj>=0 && fj<gn ) {
						if( ! Z[fi][fj] ) {
							Z_tmp[i][j] = 0;
						}
					}
				}
			}
		} END_FOR;
		
		OPENMP_FOR FOR_EVERY_CELL(gn) {
			Z[i][j] = Z_tmp[i][j];
		} END_FOR;
	}
	
	OPENMP_FOR FOR_EVERY_CELL(gn) {
		if( A[i][j] == WALL ) Z[i][j] = 0;
	} END_FOR;
}

void adaptive::computeDeepZone( char **A, char **DeepZone, int gn ) {
	findDeepZone( A, DeepZone, depeth+1, gn );
}

void adaptive::initial_merge( Sorter *sorter, char **A, std::vector<particle *> &particles, FLOAT density ) {
	int gn = sorter->getCellSize();
	if( density != 0.5 ) {
		printf( "Density must be 0.5. Exiting...\n" );
		exit(0);
	}
	
	if( ! Z ) Z = alloc2D<char>(gn);
	
	// First Find Deep Cell Zone
	findDeepZone( A, Z, depeth+1, gn );
	
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		particles[n]->pref = NULL;
		particles[n]->remove = 0;
	}
	
	OPENMP_FOR FOR_EVERY_CELL(gn) {
		if( Z[i][j] ) {
			vector<particle *> cell_particles = sorter->getNeigboringParticles_cell( i, j, 0, 0 );
			if( cell_particles.size() == MAX_SPLIT ) {
				particle &p = *cell_particles[0];
				FLOAT pos[2] = { (FLOAT)((i+0.5)/gn), (FLOAT)((j+0.5)/gn) };
				// Slightly Move
				p.p[0] = pos[0];
				p.p[1] = pos[1];
				p.u[0] = p.u[0];
				p.u[1] = p.u[1];
				p.m = 4.0;
				p.hosts[0].m = p.m;
				p.level = MAX_SPLIT;
				for( int n=1;n<MAX_SPLIT;n++ ) cell_particles[n]->remove = 1;
			}
		}
	} END_FOR;
	
	// Remove
	cleanParticles(sorter,particles);
}

char** adaptive::resample( Sorter *sorter, char **A, std::vector<particle *> &particles, FLOAT density, bool onlysplit ) {
	int gn = sorter->getCellSize();
	if( ! Z ) Z = alloc2D<char>(gn);
	
	if( ! onlysplit ) {
		// First Find Deep Cell Zone
		findDeepZone( A, Z, depeth+1, gn );
		
		// Second Merge Particles
		while( merge( sorter, A, Z, particles, density, gn )) {};
	}
	
	if( ! onlysplit ) {
		// First Find Deep Cell Zone
		findDeepZone( A, Z, depeth, gn );
	} else {
		FOR_EVERY_CELL(gn) {
			Z[i][j] = 0;
		} END_FOR;
	}
	
	// Third Split Particles
	while( split( sorter, A, Z, particles, density, gn )) {};
	
	sorter->markWater(A,density);
	findDeepZone( A, Z, depeth, gn );
	
	return Z;
}