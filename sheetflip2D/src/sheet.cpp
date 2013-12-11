/*
 *  sheet.cpp
 *  flip2d
 */

#include "sheet.h"
#include "utility.h"
#include "kernel.h"
#include <list>
#include "OBB.h"
#include <math.h>
using namespace std;

static void recalculate_mass( std::vector<particle *> particles ) {
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		particles[n]->split_cnt = 0;
	}
	
	for( int n=0; n<particles.size(); n++ ) {
		particle &p = *particles[n];
		if( p.type == WALL ) continue;
		
		FLOAT sum = 0.0;
		for( int m=0; m<p.hosts.size(); m++ ) {
			mass &aMass = p.hosts[m];
			sum += aMass.m;
			p.hosts[m].ref->split_cnt ++;
		}
		p.m = sum;
	}
}

static void compute_density_among_candidates( std::vector<particle *> &candidates, Sorter *cand_sort ) {
	int gn = cand_sort->getCellSize();
	cand_sort->sort(candidates);
	
	OPENMP_FOR for( int n=0; n<candidates.size(); n++ ) {
		particle &p = *candidates[n];
		// Find Neighbors
		vector<particle *> neighbors = cand_sort->getNeigboringParticles_cell(fmax(0,fmin(gn-1,gn*p.p[0])),
																			  fmax(0,fmin(gn-1,gn*p.p[1])),1,1);
		// Find A Nearby Particle
		FLOAT wsum = 0.0;
		for( int m=0; m<neighbors.size(); m++ ) {
			particle &np = *neighbors[m];
			if( np.type == WALL ) continue;
			wsum += kernel::smooth_kernel( hypot2(np.p[0]-p.p[0],np.p[1]-p.p[1]), 1.5/gn );
		}
		p.dens = wsum;
	}
}

static particle * search( FLOAT from[2], std::vector<particle *> &candidates, FLOAT d0, Sorter *cand_sort ) {
	if( from[0] > 0.0 ) {
		// Pop out every nearby particles
		int gn = cand_sort->getCellSize();
		vector<particle *> neighbors = cand_sort->getNeigboringParticles_cell(fmax(0,fmin(gn-1,gn*from[0])),
																			  fmax(0,fmin(gn-1,gn*from[1])),1,1);
		
		// Mark Nearby Particles
		particle *moderate = NULL;
		FLOAT mind = 99999.0;
		
		for( int m=0; m<neighbors.size(); m++ ) {
			particle &np = *neighbors[m];
			if( np.type == WALL ) continue;
			FLOAT len2 = hypot2( np.p[0]-from[0], np.p[1]-from[1] );
			if( len2 < (a9*d0)*(a9*d0) ) {
				np.remove = 1;
			} else {
				if( len2 < (a11*d0)*(a11*d0) ) {
					if( len2 < mind ) {
						mind = len2;
						moderate = neighbors[m];
					}
				}				
			}
		}
		
		// Pop Out From List
		if( moderate ) moderate->remove = 1;
		for( vector<particle *>::iterator iter=candidates.begin(); iter!=candidates.end(); ) {
			if( (*iter)->remove ) {
				if( *iter != moderate ) delete *iter;
				iter = candidates.erase(iter);
			} else iter++;
		}
		cand_sort->sort(candidates);
		
		if( moderate ) {
			return moderate;
		} else {
			FLOAT dummy[2] = { -1.0, -1.0 };
			return search( dummy, candidates, d0, cand_sort );
		}
	} else {
		FLOAT mindens = 9999999.0;
		vector<particle *>::iterator pick;
		bool found = false;
		for( vector<particle *>::iterator iter=candidates.begin(); iter!=candidates.end(); iter++) {
			if( (*iter)->dens < mindens ) {
				mindens = (*iter)->dens;
				pick = iter;
				found = true;
			}
		}
		
		if( found ) {
			particle *res = *pick;
			candidates.erase(pick);
			return res;
		}
	}
	return NULL;
}

void sheet::keepThinSheet( char **A, Sorter *sorter, std::vector<particle *> &particles, FLOAT obb_rate, FLOAT density, FLOAT wall_thick ) {
	vector<particle *> new_particles;
	int gn = sorter->getCellSize();
	FLOAT d0 = density/gn;
	
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		particle &p = *particles[n];
		particle_iterator *last = NULL;
		p.tmps = NULL;
        
		//if( p.m < 0.01 ) continue;
		if( p.dens < a6 ) continue;
		if( p.level > 1 ) continue;
		//if( p.p[1] > 0.9 ) continue;
		
        // Skip If Adjacent To Up-Bottom Walls
        if(0) {
            int i = fmin(gn-2.0*gn*wall_thick,fmax(2.0*gn*wall_thick,p.p[0]*gn));
            int j = fmin(gn-2.0*gn*wall_thick,fmax(2.0*gn*wall_thick,p.p[1]*gn));
            if( A[i][max(0,j-1)] == WALL || A[i][min(gn-1,j+1)] == WALL ) continue;
        }
        
		if( p.type == FLUID && p.obb.ratio < obb_rate ) {
			// Find Neighbors
			vector<particle *> neighbors = sorter->getNeigboringParticles_cell(fmax(0,fmin(gn-1,gn*p.p[0])),
																			 fmax(0,fmin(gn-1,gn*p.p[1])),1,1);
				
			// Find A Nearby Particle
			for( int m=0; m<neighbors.size(); m++ ) {
				particle &np = *neighbors[m];
				if( np.type == WALL ) continue;
				if( np.level > 1 ) continue;
				if( np.dens < a6 ) continue;
				
				// Check The Sparcity
				// if( 0.5*(np.dens+p.dens) > a12 ) continue;
				
				if( p.p[0] < np.p[0] ) {
                    
                    // Skip If Adjacent To Up-Bottom Walls
                    if(0) {
                        int i = fmin(gn-2.0*gn*wall_thick,fmax(2.0*gn*wall_thick,np.p[0]*gn));
                        int j = fmin(gn-2.0*gn*wall_thick,fmax(2.0*gn*wall_thick,np.p[1]*gn));
                        if( A[i][max(0,j-1)] == WALL || A[i][min(gn-1,j+1)] == WALL ) continue;
                    }
                    
					// Check If The Distance Is Too Far
					FLOAT d2 = hypot2( np.p[0]-p.p[0], np.p[1]-p.p[1] );
					if( d2 > a9*d0*a9*d0 && d2 < a10*a10*d0*d0) {
						
						// Check This Pair Is Separating
						if( (np.p[0]-p.p[0])*(np.u[0]-p.u[0]) + (np.p[1]-p.p[1])*(np.u[1]-p.u[1]) > 0 ) {
							
							// Found Candidate. Check The Position Is Sparse Around There
							FLOAT pos[2] = { (FLOAT)(0.5*(np.p[0]+p.p[0])), (FLOAT)(0.5*(np.p[1]+p.p[1])) };
                            vector<particle *> neighbors_cand = sorter->getNeigboringParticles_cell(fmax(0,fmin(gn-1,gn*pos[0])),
                                                                                                  fmax(0,fmin(gn-1,gn*pos[1])),1,1);
                            
							bool sparse = true;
							for( int q=0; q<neighbors_cand.size(); q++ ) {
								particle qp = *neighbors_cand[q];
								if( qp.type == WALL ) continue;
								if( hypot2( qp.p[0]-pos[0], qp.p[1]-pos[1] ) < (a9*d0)*(a9*d0) ) {
									sparse = false;
									break;
								}
							}
							
							// If Sparse, Insert A New Particle
							if( sparse ) {
								particle *newp = new particle;
								newp->p[0] = pos[0];
								newp->p[1] = pos[1];			
								newp->u[0] = 0.5*(np.u[0]+p.u[0]);
								newp->u[1] = 0.5*(np.u[1]+p.u[1]);
								newp->n[0] = 0.5*(np.n[0]+p.n[0]);
								newp->n[1] = 0.5*(np.n[1]+p.n[1]);
								newp->f[0] = 0.0;
								newp->f[1] = 0.0;
								newp->type = FLUID;
								newp->pref = NULL;
								newp->curv = 0.5*(np.curv+p.curv);
								newp->mark = p.mark+1;
								
								// Build Pre-Hosts
								mass aMass = { 0.0, 0, particles[n] };
								mass bMass = { 0.0, 0, neighbors[m] };
								newp->hosts.push_back(aMass);
								newp->hosts.push_back(bMass);
								
								newp->remove = 0;
								newp->tmps = NULL;
								newp->dens = 0.5*(np.dens+p.dens);
								newp->tmp[0][0] = newp->dens;
                                newp->obb = p.obb;
								
								newp->level = p.level;
								
                                // SPH
                                newp->SPH_dens0 = 0.5*(np.SPH_dens0+p.SPH_dens0);
                                
								particle_iterator *iter = new particle_iterator;
								iter->p = newp;
								iter->next = NULL;
								if( ! last ) {
									p.tmps = iter;
									last = iter;
								} else {
									last->next = iter;
									last = iter;
								}
							}
						}
					}
				}
			}
		}
	}
	
	// Collect Candidates
    vector<particle *> generated_particles;
    vector<particle *> inserted_particles;
    
	bool refresh = false;
	// Pick Up Generated Particles
    int b4_size = particles.size();
	for( int n=0; n<b4_size; n++ ) {
		particle &p = *particles[n];
		while( p.tmps ) {
			generated_particles.push_back( p.tmps->p );
			particle_iterator *old = p.tmps;
			p.tmps = p.tmps->next;
			delete old;
			refresh = true;
		}
	}
	
	// Sort Particle
	static Sorter *cand_sort = new Sorter(gn);
	cand_sort->sort( generated_particles );
	
	// Compute Density Among Candidates
	compute_density_among_candidates( generated_particles, cand_sort );
	
	// Insert Particles
	FLOAT from[2] = { -1.0, -1.0 };
	while( generated_particles.size() ) {
		particle *newone = search( from, generated_particles, d0, cand_sort );
		if( ! newone ) break;
		
		// Split Mass
		particle *parents[2] = { newone->hosts[0].ref, newone->hosts[1].ref };
		newone->hosts.clear();
		for( int side=0; side<2; side++ ) {
			for( int n=0; n<parents[side]->hosts.size(); n++ ) {
				mass &aMass = parents[side]->hosts[n];
				mass newMass = { (FLOAT)(aMass.m/3.0), 0, aMass.ref };
				
				aMass.m *= (2.0/3.0);
				newone->hosts.push_back(newMass);
			}
		}
		
		newone->dens = newone->tmp[0][0];
		particles.push_back( newone );
		from[0] = newone->p[0];
		from[1] = newone->p[1];
		refresh = true;
	}
	
	// Resort If Neccessary
	if( refresh ) sorter->sort(particles);
	
	// Recalculate Mass
	recalculate_mass(particles);
}

void sheet::collapseThinSheet( char **A, Sorter *sorter, std::vector<particle *> &particles, FLOAT obb_rate, FLOAT density ) {
	int gn = sorter->getCellSize();
	FLOAT d0 = density/gn;
	
	// Remove Non-Thin Particles
#if 1 // What?
	OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		particle &p = *particles[n];
		p.remove = 0;
		if( p.type == FLUID ) {
			if( p.mark > 0 && (p.dens > a12 && p.obb.ratio > obb_rate && rand()%2==0 ) ) {
				p.remove = 1;
			}
		}
	}
    
    // Check If The Inserted Particle Is Too Close To Some Other Particle
    OPENMP_FOR for( int n=0; n<particles.size(); n++ ) {
		particle &p = *particles[n];
        if( p.mark ) {
            vector<particle *> neighbors = sorter->getNeigboringParticles_cell(fmax(0,fmin(gn-1,gn*p.p[0])),
                                                                             fmax(0,fmin(gn-1,gn*p.p[1])),1,1);
            bool sparse = true;
            for( int q=0; q<neighbors.size(); q++ ) {
                particle &qp = *neighbors[q];
                if( qp.type == WALL ) continue;
                if( particles[n] > neighbors[q] ) {
                    if( hypot2( qp.p[0]-p.p[0], qp.p[1]-p.p[1] ) < (a13*d0)*(a13*d0) ) {
                        sparse = false;
                        break;
                    }
                }
            }
            if( ! sparse && rand()%2==0 ) p.remove = 1;
        }
    }
	
	// Remove
	cleanParticles(sorter,particles);
#endif
	
	
	// Recalculate Mass
	recalculate_mass(particles);
}
