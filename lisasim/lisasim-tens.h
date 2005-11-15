/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */
 
#ifndef _LISASIM_TENS_H_
#define _LISASIM_TENS_H_

class Tensor;

class Vector {
 private:
    friend class Tensor;
    double c[3];

 public:
    Vector() {};

    Vector(double a) {
	c[0] = a; c[1] = a; c[2] = a;
    }

    double& operator[](int i) {
	return c[i];
    };

    Vector& operator=(const Vector &vec) {
	c[0] = vec.c[0]; c[1] = vec.c[1]; c[2] = vec.c[2];
	return *this;
    };

    Vector& setsum(const Vector &add1,const Vector &add2);
    Vector& setdifference(const Vector &add1,const Vector &add2);

    Vector& setproduct(const double fac);
    Vector& setproduct(const Tensor &mat, const Vector &vec);

    double dotproduct();
    double dotproduct(const Vector &vec);
    friend double dotproduct(const Vector &vec1,const Vector &vec2);

    Vector& setnormalized();
};

class Tensor {
 private:
    friend class Vector;
    double c[9];
    
 public:
    Tensor() {};

    Tensor(double a) {
	for(int i=0;i<9;i++)
	    c[i] = a;
    }
    
    double* operator[](int i) {
	return c + i*3;
    };

    Tensor& operator=(const Tensor &tens) {
	c[0] = tens.c[0]; c[1] = tens.c[1]; c[2] = tens.c[2];
	c[3] = tens.c[3]; c[4] = tens.c[4]; c[5] = tens.c[5];
	c[6] = tens.c[6]; c[7] = tens.c[7]; c[8] = tens.c[8];
	return *this;
    };

    Tensor& setproduct(const Tensor &fac1, const Tensor &fac2);

    Tensor& settranspose();
    Tensor& settranspose(const Tensor& tens);

    Tensor& seteuler(double b, double l, double p);
};

#endif /* _LISASIM_TENS_H_ */
