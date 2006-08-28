/* $Id$
 * $Date$
 * $Author$
 * $Revision$
 */

#ifndef _LISASIM_TDISIGNAL_H_
#define _LISASIM_TDISIGNAL_H_

#include "lisasim-tdi.h"
#include "lisasim-lisa.h"
#include "lisasim-wave.h"

class TDIsignal : public TDI {
 private:
    LISA *lisa, *phlisa;
    WaveObject *wave;
    
    double psi(Wave *nwave, Vector &lisan, double t);
    
 public:
    TDIsignal(LISA *mylisa, WaveObject *mywave);

    // change the physical LISA

    void setphlisa(LISA *mylisa);

    // won't reset wave objects

    void reset();

    // defined here only for comparison with the LISA simulator

    double M(double t);
    double N(double t);
    double O(double t);

    double y(int send, int link, int recv, int ret1, int ret2, int ret3, double t);
    double y(int send, int link, int recv, int ret1, int ret2, int ret3, int ret4, int ret5, int ret6, int ret7, double t);

    double Phi(int slink,double t);
};

#endif /* _LISASIM_TDISIGNAL_H_ */
