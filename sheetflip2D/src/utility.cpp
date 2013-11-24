/*
 *  utility.cpp
 *  flip2d
 */

#include "utility.h"
#undef min
#undef max
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void my_rand_shuffle( std::vector<ipos> &waters ) {
	random_shuffle( waters.begin(), waters.end() );
}

int cleanParticles( sorter *sort, std::vector<particle *> &particles ) {
	int cnt = 0;
	int size = 0;
	for( int n=0; n<particles.size(); n++ ) {
		particle &p = *particles[n];
		if( p.remove ) {
			if( p.mark > 0 ) {
				for( int n=0; n<p.hosts.size(); n++ ) {
					mass &aMass = p.hosts[n];
					particle *ref = aMass.ref;
					ref->hosts[0].m += aMass.m;
				}
			}
		} else {
			size++;
		}
	}
	if( particles.size() == size ) return 0;
	
	std::vector<particle *> new_particles;
	new_particles.resize(size);
	int index=0;
	for( int n=0; n<particles.size(); n++ ) {
		if( particles[n]->remove ) {
			delete particles[n];
			cnt ++;
		} else {
			new_particles[index++] = particles[n];
		}
	}
	particles = new_particles;
	sort->sort(particles);
	return cnt;
}