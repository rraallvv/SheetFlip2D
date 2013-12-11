/*
 *  sorter.h
 *  flip2d
 */

#include "common.h"
#include <vector>
#ifndef _SORTER_H
#define _SORTER_H

class Sorter {
public:
	Sorter( int gn );
	~Sorter();
	
	void sort( std::vector<particle *> &particles );
	std::vector<particle *> getNeigboringParticles_wall( int i, int j, int w, int h );
	std::vector<particle *> getNeigboringParticles_cell( int i, int j, int w, int h );
	FLOAT levelset( int i, int j, FLOAT density );
	
	int	 getCellSize(){ return gn; }
	void markWater( char **A, FLOAT density );
	
protected:
	std::vector<particle *> **cells;
	int gn;
};

#endif