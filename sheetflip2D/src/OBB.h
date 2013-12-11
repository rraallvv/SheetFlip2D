/*
 *  OBB.h
 *  flip2d
 */

#include "common.h"
#include "sorter.h"
#include <vector>

#ifndef _OBB_H
#define _OBB_H

OBB isotropicOBB();
void buildOBB( char **A, char **DeepZone, Sorter *sorter, std::vector<particle *> &particles );
OBB buildNearbyOBB( Sorter *sorter, particle *p, FLOAT re );
OBB buildOBB( std::vector<particle *> particles, FLOAT cp[2], FLOAT re );

#endif