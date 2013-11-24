/*
 *  flip2d.h
 *  smoke
 *
 */

#include "common.h"

namespace flip2d {
	void init();
	void reshape( int w, int h );
	void display();
	void mouse( FLOAT x, FLOAT y, int state );
	void motion( FLOAT x, FLOAT y, FLOAT dx, FLOAT dy );
	void keyDown( unsigned char key );
}