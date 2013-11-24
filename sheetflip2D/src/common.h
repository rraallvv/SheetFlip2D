/*
 *  common.h
 *  flip2d
 */

#include "math.h"
#include <vector>

#ifndef _COMMON_H
#define _COMMON_H

#define FLOAT		float

#define AIR			0
#define FLUID		1
#define WALL		(-1)

#define	CEILING		1.0
#define PI          3.14159265

///// Parameters

#define a1		1.0
#define	a2		4.0
#define a6      0.001
#define a7		0.7
#define a9		0.8
#define a10		3.5
#define a11		2.0
#define a12		0.3
#define a13		0.2

// Vector Type.
typedef struct _Vector2 {
	FLOAT x, y;
	FLOAT& operator[]( int i ) {
		if( i==0 ) return x;
		else return y;
	}
	_Vector2 operator+(_Vector2 v) {
		_Vector2 r = { v.x+x,v.y+y };
		return r;
	}
	_Vector2 operator-(_Vector2 v) {
		_Vector2 r = { x-v.x,y-v.y };
		return r;
	}
	void operator+=(_Vector2 v) {
		x+=v.x;
		y+=v.y;
	}
	_Vector2 operator*(double s) {
		_Vector2 r = { (FLOAT)s*x,(FLOAT)s*y };
		return r;
	}
	double operator*(_Vector2 v) {
		return x*v.x+y*v.y;
	}
	_Vector2 operator/(double s) {
		_Vector2 r = { (FLOAT)(x/s),(FLOAT)(y/s) };
		return r;
	}
	friend _Vector2 operator*( double s, _Vector2 vec )
	{
		return vec*s;
	}
	bool operator==(_Vector2 v) {
		return x==v.x && y==v.y;
	}
	bool empty(){ return (x==0.0 && y==0.0); }
	double length(){ return hypot(x,y); }
	double length2(){ return x*x+y*y; }
	_Vector2 normalize() {
		return (*this)/length();
	}
} Vector2;
const static Vector2 zeroVec = { 0.0, 0.0 };

////////////////////////////////////////

struct OBB {
	FLOAT u[2][2];
	FLOAT c[2];
	FLOAT ratio;
};

struct _particle;
typedef struct _mass {
	FLOAT m;
	int index;
	_particle *ref;
} mass;

struct _particle_iterator;
typedef struct _particle {
	FLOAT p[2];
	FLOAT u[2];
	FLOAT f[2];
	FLOAT n[2];
	char type;
	char remove;
	FLOAT dens;
	FLOAT curv;
	FLOAT m;
	_particle *pref;
	char level; // Adaptive Resampling Level
	std::vector<mass> hosts;
	int split_cnt;
	FLOAT tmp[6][2];
	OBB obb;
	_particle_iterator *tmps;
	char mark;
    //// Variable For SPH
    FLOAT SPH_dens0;
    FLOAT SPH_dens;
    FLOAT SPH_press;
    Vector2 F_press;
	Vector2 F_visc;
} particle;

typedef struct _particle_iterator {
	particle *p;
	_particle_iterator *next;
} particle_iterator;

typedef struct _ipos {
	int i; int j;
} ipos;

#endif