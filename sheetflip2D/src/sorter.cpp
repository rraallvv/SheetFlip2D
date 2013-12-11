/*
 *  sorter.cpp
 *  flip2d
 */

#include "sorter.h"
#include "utility.h"
#include "adaptive.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

Sorter::Sorter( int gn ) {
	cells = alloc2D<vector<particle *> >(gn);
	this->gn = gn;
}

Sorter::~Sorter() {
}

void Sorter::sort( std::vector<particle *> &particles ) {
	// Clear All Cells
	FOR_EVERY_CELL(gn) {
		cells[i][j].clear();
	} END_FOR
	
	// Store Into The Cells
	for( int n=0; n<particles.size(); n++ ) { 
		particle *p = particles[n];
		FLOAT pos[2];
		for( int k=0; k<2; k++ ) {
			pos[k] = p->p[k];
		}
		int i = fmax(0,fmin(gn-1,gn*pos[0]));
		int j = fmax(0,fmin(gn-1,gn*pos[1]));
		cells[i][j].push_back(p);
	}
}

std::vector<particle *> Sorter::getNeigboringParticles_wall( int i, int j, int w, int h ) {
	std::vector<particle *> res;
	int sum = 0;
	for( int si=i-w; si<=i+w-1; si++ ) for( int sj=j-h; sj<=j+h-1; sj++ ) {
		if( si < 0 || si > gn-1 || sj < 0 || sj > gn-1 ) continue;
		sum += cells[si][sj].size();
	}
	res.resize(sum);
	
	int index = 0;
	for( int si=i-w; si<=i+w-1; si++ ) for( int sj=j-h; sj<=j+h-1; sj++ ) {
		if( si < 0 || si > gn-1 || sj < 0 || sj > gn-1 ) continue;
		for( int k=0; k<cells[si][sj].size(); k++ ) { 
			particle *p = cells[si][sj][k];
			res[index++] = p;
		}
	}
	return res;
}

std::vector<particle *> Sorter::getNeigboringParticles_cell( int i, int j, int w, int h ) {
	std::vector<particle *> res;
	int sum = 0;
	for( int si=i-w; si<=i+w; si++ ) for( int sj=j-h; sj<=j+h; sj++ ) {
		if( si < 0 || si > gn-1 || sj < 0 || sj > gn-1 ) continue;
		sum += cells[si][sj].size();
	}
	res.resize(sum);
	
	int index = 0;
	for( int si=i-w; si<=i+w; si++ ) for( int sj=j-h; sj<=j+h; sj++ ) {
		if( si < 0 || si > gn-1 || sj < 0 || sj > gn-1 ) continue;
		for( int k=0; k<cells[si][sj].size(); k++ ) { 
			particle *p = cells[si][sj][k];
			res[index++] = p;
		}
	}
	return res;
}

FLOAT Sorter::levelset( int i, int j, FLOAT density ) {
	FLOAT accm = 0.0;
	for( int k=0; k<cells[i][j].size(); k++ ) { 
		if( cells[i][j][k]->type == FLUID ) {
			accm += cells[i][j][k]->m;
		} else {
			return 1.0;
		}
	}
	FLOAT n0 = 1.0/(density*density);
	return 0.2*n0-accm;
}

void Sorter::markWater( char **A, FLOAT density ) {
	FOR_EVERY_CELL(gn) {
		bool Wall = false;
		for( int k=0; k<cells[i][j].size(); k++ ) { 
			if( cells[i][j][k]->type == WALL ) {
				Wall = true;
				break;
			}
		}
		if( Wall ) A[i][j] = WALL;
		else A[i][j] = levelset( i, j, density ) < 0.0 ? FLUID : AIR;
	} END_FOR;
	
	// Search For Merged Particle
	FOR_EVERY_CELL(gn) {
		if( A[i][j] == WALL ) continue;
		int found_cnt = 0;
		for( int nx=-1; nx<=1; nx++ ) {
			for( int ny=-1; ny<=1; ny++ ) {
				int fi = i+nx;
				int fj = j+ny;
				if( fi>=0 && fi<gn && fj>=0 && fj<gn ) {
					int found[2] = { 0, 0 };
					for( int k=0; k<cells[fi][fj].size(); k++ ) { 
						if( cells[fi][fj][k]->type == FLUID ) {
							if( cells[fi][fj][k]->level == 1 ) {
								found[0]++;
							} else {
								found[1]++;
							}
						}
					}
					if( found[0] < found[1] ) found_cnt++;
				}
				if( cells[fi][fj].empty() ) found_cnt--;
			}
		}
		if( found_cnt > 0 ) A[i][j] = FLUID;
	} END_FOR;
}