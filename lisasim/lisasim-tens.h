#ifndef _LISASIM_TENS_H_
#define _LISASIM_TENS_H_

#include "nr.h"
#include "nrutil.h"

class Tensor : public Mat_IO_DP {
    public:

    Tensor();
    
    void seteuler(double d, double a, double p);
    
    void setproduct(Tensor &fac1, Tensor &fac2);
                                
    void settranspose();
};

class Vector : public Vec_IO_DP {
    public:
    
    Vector();
    
    void setproduct(Tensor &mat, Vector &vec);

    double dotproduct(Vector &vec);
};

#endif /* _LISASIM_TENS_H_ */
