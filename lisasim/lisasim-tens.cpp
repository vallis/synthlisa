/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

#include <math.h>

#include "lisasim-tens.h"

// Vector member functions in the header file (all inlined)

// Tensor member functions

Tensor& Tensor::setproduct(const Tensor &fac1, const Tensor &fac2) {
    c[0] = fac1.c[0]*fac2.c[0] + fac1.c[1]*fac2.c[3] + fac1.c[2]*fac2.c[6];
    c[3] = fac1.c[3]*fac2.c[0] + fac1.c[4]*fac2.c[3] + fac1.c[5]*fac2.c[6];
    c[6] = fac1.c[6]*fac2.c[0] + fac1.c[7]*fac2.c[3] + fac1.c[8]*fac2.c[6];

    c[1] = fac1.c[0]*fac2.c[1] + fac1.c[1]*fac2.c[4] + fac1.c[2]*fac2.c[7];
    c[4] = fac1.c[3]*fac2.c[1] + fac1.c[4]*fac2.c[4] + fac1.c[5]*fac2.c[7];
    c[7] = fac1.c[6]*fac2.c[1] + fac1.c[7]*fac2.c[4] + fac1.c[8]*fac2.c[7];

    c[2] = fac1.c[0]*fac2.c[2] + fac1.c[1]*fac2.c[5] + fac1.c[2]*fac2.c[8];
    c[5] = fac1.c[3]*fac2.c[2] + fac1.c[4]*fac2.c[5] + fac1.c[5]*fac2.c[8];
    c[8] = fac1.c[6]*fac2.c[2] + fac1.c[7]*fac2.c[5] + fac1.c[8]*fac2.c[8];

    return *this;
}

Tensor& Tensor::seteuler(double b, double l, double p) {
    c[0] =  cos(p)*sin(l) - cos(l)*sin(b)*sin(p);
    c[1] = -cos(l)*cos(p)*sin(b) - sin(l)*sin(p);
    c[2] = -cos(l)*cos(b);

    c[3] = -cos(l)*cos(p) - sin(l)*sin(b)*sin(p);
    c[4] = -cos(p)*sin(l)*sin(b) + cos(l)*sin(p);
    c[5] = -cos(b)*sin(l);
    
    c[6] =  cos(b)*sin(p);
    c[7] =  cos(b)*cos(p);
    c[8] = -sin(b);

    return *this;
}

Tensor& Tensor::settranspose() {
    double tmp;
    
    tmp=c[1]; c[1]=c[3]; c[3]=tmp;
    tmp=c[2]; c[2]=c[6]; c[6]=tmp;
    tmp=c[5]; c[5]=c[7]; c[7]=tmp;

    return *this;
}

Tensor& Tensor::settranspose(const Tensor& tens) {
    c[0] = tens.c[0]; c[1] = tens.c[3]; c[3] = tens.c[1];
    c[4] = tens.c[4]; c[2] = tens.c[6]; c[6] = tens.c[2];
    c[8] = tens.c[8]; c[5] = tens.c[7]; c[7] = tens.c[5];

    return *this;
}
