#include "lisasim-tens.h"

Tensor::Tensor() : Mat_IO_DP(0.0,3,3) {}

void Tensor::setproduct(Tensor &fac1, Tensor &fac2) {
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++) {
            double pr = 0.0;
                
            for(int k=0;k<3;k++)
                pr += fac1[i][k] * fac2[k][j];
                
            (*this)[i][j] = pr;                    
        }
}

void Tensor::seteuler(double d, double a, double p) {
    (*this)[0][0] =  cos(p)*sin(a) - cos(a)*sin(d)*sin(p);
    (*this)[0][1] = -cos(a)*cos(p)*sin(d) - sin(a)*sin(p);
    (*this)[0][2] = -cos(a)*cos(d);

    (*this)[1][0] = -cos(a)*cos(p) - sin(a)*sin(d)*sin(p);
    (*this)[1][1] = -cos(p)*sin(a)*sin(d) + cos(a)*sin(p);
    (*this)[1][2] = -cos(d)*sin(a);
    
    (*this)[2][0] =  cos(d)*sin(p);
    (*this)[2][1] =  cos(d)*cos(p);
    (*this)[2][2] = -sin(d);
}

void Tensor::settranspose() {
    double tmp;
    
    tmp=(*this)[0][1]; (*this)[0][1]=(*this)[1][0]; (*this)[1][0]=tmp;
    tmp=(*this)[0][2]; (*this)[0][2]=(*this)[2][0]; (*this)[2][0]=tmp;
    tmp=(*this)[1][2]; (*this)[1][2]=(*this)[2][1]; (*this)[2][1]=tmp;
}

Vector::Vector() : Vec_IO_DP(0.0,3) {}

void Vector::setproduct(Tensor &mat, Vector &vec) {
    for(int i=0;i<3;i++) {
        double pr=0.0;
        
        for(int k=0;k<3;k++)
            pr += mat[i][k] * vec[k];

        (*this)[i] = pr;
    }
}

double Vector::dotproduct(Vector &vec) {
    double pr=0.0;
    
    for(int k=0;k<3;k++)
        pr += (*this)[k] * vec[k];

    return pr;
}
