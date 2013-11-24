#include "kernel.h"
#include <math.h>


FLOAT kernel::smooth_kernel( FLOAT r2, FLOAT h ) {
#if 0
    r2 = r2/(h*h);
    return fmax( 1.0-r2, 0.0 );
#else
    r2 = r2/(h*h);
    FLOAT hh = 0.1;
    FLOAT hh2 = (1.0+hh)*(1.0+hh);
    return r2 < 1+hh ? (1.0-r2)*(hh2-r2) : 0.0;
#endif
}

FLOAT kernel::sharp_kernel( FLOAT r2, FLOAT h ) {
    return fmax( h*h/fmax(r2,1.0e-5) - 1.0, 0.0 );
}
