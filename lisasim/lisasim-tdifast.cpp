#include "lisasim.h"




#define Check(value) (CheckError((value), (__FILE__), (__LINE__)))
//#define Check(value) 1

void CheckError(bool value, char *file, int linenum) {
  if (!value) {
    cout << "Value not cached in "<< file << " line "<< linenum <<"."<< endl;
    exit(0);
  }
  //return(1);
}

// double TDI::psi(int arm, double t) {
double TDIfast::psi(double lisanTemp[], double twave) {
  Vector lisan;
    // lisa->putn(lisan,arm,t);  // Cache?
  /*    Check(Checkputn[arm]);
    lisan = Storedputn[arm][0];  // Incorrect Testing
  */
  lisan[0] = lisanTemp[0];
  lisan[1] = lisanTemp[1];
  lisan[2] = lisanTemp[2];

    Tensor cwave;
    wave->putwave(cwave,twave);
    
    Vector tmp;
    tmp.setproduct(cwave,lisan);

    return(0.5 * lisan.dotproduct(tmp));
}
/*
double TDIfast::retard(int craft, double t) {
    Vector lisap;
    lisa->putp(lisap,craft,t);
    cout << "Should not be used" << endl;
    return(-lisap.dotproduct(wave->k));
} */

double TDIfast::retardfast(int craft, int tIndex) {
    Vector lisap;
#ifdef DEBUG    
    Check(Checkputp[craft]);
#endif
    lisap = Storedputp[craft][tIndex];
    
    return(-lisap.dotproduct(wave->k));
}



TDIfast::TDIfast(LISA *mylisa, Wave *mywave, double mysrate, long mysamples):TDI(mylisa,mywave) {
  lisa = mylisa;
  wave = mywave;
  srate = mysrate * 3.17098E-8; // convert from seconds to years
  samples = mysamples;

  cout << "sample rate in years is" << srate << endl;

  Storedarmlength = new double*[4];
  Checkarmlength = new bool[4];
  for (int i = 0; i < 4; i++) {
    Checkarmlength[i] = 0;
  }

  Storedputn = new Vector*[4];
  Checkputn = new bool[4];
  for (int i = 0; i < 4; i++) {
    Checkputn[i] = 0;
  }  

  Storedputp = new Vector*[4];
  Checkputp = new bool[4];
  for (int i = 0; i < 4; i++) {
    Checkputp[i] = 0;
  }  
  
  StoredRetardedTime = new double*[4*4*4];
  CheckRetardedTime = new bool**[4];
  //  StoredputnRet = new Vector*****[4];
  StoredputnRet = new double**[2*4*4*4*4];
  CheckputnRet = new bool***[4];
  
  

  // i = ret1, j = ret2, k = ret3, l = link, (In all cases 0=subtract nothing, 1...3 = subtract armlength(1...3, tIndex))
  for (int i = 0; i < 4; i++) {
    //    StoredRetardedTime[i] = new double**[4];
    CheckRetardedTime[i] = new bool*[4];
    //StoredputnRet[i] = new Vector****[4];
    CheckputnRet[i] = new bool**[4];
    for (int j = 0; j < 4; j++) {
      //StoredRetardedTime[i][j] = new double*[4];
      CheckRetardedTime[i][j] = new bool[4];
      //StoredputnRet[i][j] = new Vector***[4];
      CheckputnRet[i][j] = new bool*[4];
      for (int k = 0; k < 4; k++) {
	CheckRetardedTime[i][j][k] = 0;
	//StoredputnRet[i][j][k] = new Vector**[4];
	CheckputnRet[i][j][k] = new bool[4];
	for (int l = 0; l < 4; l++){
	  //StoredputnRet[i][j][k][l] = new Vector*[2];
	  //	  CheckputnRet[i][j][k][l] = new bool;
	  CheckputnRet[i][j][k][l] = 0;
	  
	}
      }
    }
  }
  
}

inline int SRTI(int i, int j, int k) {
  return(k + 4*j + 16*i);
}

inline int SPNRI(int i, int j, int k, int l, int d) {
  return(d + 2*l+ 2*4*k + 2*4*4*j + 2*4*4*4*i);
}



TDIfast::~TDIfast() {
  
  
  for (int i = 0; i < 4; i++) {
    if (Checkarmlength[i]) {
      delete [] Storedarmlength[i];
    }
  }
  delete [] Storedarmlength;
  delete [] Checkarmlength;

  for (int i = 0; i < 4; i++) {
    if (Checkputn[i]) {
      delete [] Storedputn[i];
    }
  }
  delete [] Storedputn;
  delete [] Checkputn;

  for (int i = 0; i < 4; i++) {
    if (Checkputp[i]) {
      delete [] Storedputp[i];
    }
  }
  delete [] Storedputp;
  delete [] Checkputp;

  // i = ret1, j = ret2, k = ret3, l = link, (In all cases 0=subtract nothing, 1...3 = subtract armlength(1...3, tIndex))
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
	for (int l = 0; l < 4; l++){
	  if (CheckputnRet[i][j][k][l]) {
	    for (int tIndex = 0; tIndex <= samples; tIndex++){
	      delete [] StoredputnRet[SPNRI(i,j,k,l,0)][tIndex];
	      delete [] StoredputnRet[SPNRI(i,j,k,l,1)][tIndex];
	    }
	    delete [] StoredputnRet[SPNRI(i,j,k,l,0)];
	    delete [] StoredputnRet[SPNRI(i,j,k,l,1)];
	  }
	  //delete []  StoredputnRet[i][j][k][l];
	}
	//delete [] StoredputnRet[i][j][k];
	delete [] CheckputnRet[i][j][k];
	if (CheckRetardedTime[i][j][k]) {
	  delete [] StoredRetardedTime[SRTI(i,j,k)];
	}
      }
      //delete [] StoredputnRet[i][j];
      delete [] CheckputnRet[i][j];
      //      delete [] StoredRetardedTime[i][j];
      delete [] CheckRetardedTime[i][j];
      
    }
    //delete [] StoredputnRet[i];
    delete [] CheckputnRet[i];
    //delete [] StoredRetardedTime[i];
    delete [] CheckRetardedTime[i];
  }
  delete [] StoredputnRet;
  delete [] CheckputnRet;
  delete [] StoredRetardedTime;
  delete [] CheckRetardedTime;

}


int TDIfast::yc(int send, int link, int recv, int ret1, int ret2, int ret3) {

  
  
  int armArray[4];
  int putpArray[3];

  armArray[0] = ret1;
  armArray[1] = ret2;
  armArray[2] = ret3;
  armArray[3] = link;

  putpArray[0] = send;
  putpArray[1] = link;
  putpArray[2] = recv;

  for (int i=0; i< 4; i++) {
    if(armArray[i]) {
      if (!Checkarmlength[armArray[i]]) {
	Storedarmlength[armArray[i]] = new double[samples];
	for (int j = 0; j < samples; j++) {
	  Storedarmlength[armArray[i]][j] = lisa->armlength(armArray[i],j*srate);
	}
	Checkarmlength[armArray[i]] = 1;
      }
    }
  }

  if (!CheckRetardedTime[ret1][ret2][ret3]) {
    StoredRetardedTime[SRTI(ret1,ret2,ret3)] = new double[samples];
  }

  if (!CheckputnRet[ret1][ret2][ret3][link]) {
    StoredputnRet[SPNRI(ret1,ret2,ret3,link,0)] = new double*[samples];
    StoredputnRet[SPNRI(ret1,ret2,ret3,link,1)] = new double*[samples];
  }
    
  if ((!CheckRetardedTime[ret1][ret2][ret3]) || (!CheckputnRet[ret1][ret2][ret3][link])) {
    for (int tIndex = 0; tIndex < samples; tIndex++) {

      if (!CheckRetardedTime[ret1][ret2][ret3]) {
	double retardedtime = tIndex*srate;  
	
	if(ret1) retardedtime -= Storedarmlength[ret1][tIndex];
	if(ret2) retardedtime -= Storedarmlength[ret2][tIndex];
	if(ret3) retardedtime -= Storedarmlength[ret3][tIndex];
	
	StoredRetardedTime[SRTI(ret1,ret2,ret3)][tIndex] = retardedtime;
      }
      
      if (!CheckputnRet[ret1][ret2][ret3][link]) {
	Vector temp;
	
	StoredputnRet[SPNRI(ret1,ret2,ret3,link,0)][tIndex] = new double[3];
	StoredputnRet[SPNRI(ret1,ret2,ret3,link,1)][tIndex] = new double[3];
	
	lisa->putn(StoredputnRet[SPNRI(ret1,ret2,ret3,link,0)][tIndex], link, StoredRetardedTime[SRTI(ret1,ret2,ret3)][tIndex]);
	lisa->putn(StoredputnRet[SPNRI(ret1,ret2,ret3,link,1)][tIndex], link, StoredRetardedTime[SRTI(ret1,ret2,ret3)][tIndex] - Storedarmlength[link][tIndex]);
      }

    }
    CheckRetardedTime[ret1][ret2][ret3] = 1;
    CheckputnRet[ret1][ret2][ret3][link] = 1;  
  }
  
  if (!Checkputn[link]) {
    Storedputn[link] = new Vector[samples];
    for (int j = 0; j < samples; j++) {
      lisa->putn(Storedputn[link][j],link,j*srate);
    }
    Checkputn[link] = 1;
  }
  
  
  for (int i=0; i< 3; i++) {
    if (!Checkputp[putpArray[i]]) {
      Storedputp[putpArray[i]] = new Vector[samples];
      for (int j = 0; j < samples; j++) {
	lisa->putp(Storedputp[putpArray[i]][j],putpArray[i],j*srate);
      }
      Checkputp[putpArray[i]] = 1;
    }
  }
  
  
  // STORE retard(send,t) retard(recv,t)

  
  //    return (( psi(link, retardedtime + retard(send, t) - lisa->armlength(link,t)) -
  //              psi(link, retardedtime + retard(recv, t)) ) / denom );
  return (1);
  //  return (( psi(link, retardedtime - lisa->armlength(link,t), retardedtime + retard(send, t) - lisa->armlength(link,t)) -
  //	    psi(link, retardedtime, retardedtime + retard(recv, t)) ) / denom );
} 

double TDIfast::y(int send, int link, int recv, int ret1, int ret2, int ret3, int tIndex) {
  double retardedtime; // = tIndex*srate;
#ifdef DEBUG    
  Check(CheckRetardedTime[ret1][ret2][ret3]);
#endif 
  retardedtime = StoredRetardedTime[SRTI(ret1,ret2,ret3)][tIndex];
  
  /*    if(ret1) Check(Checkarmlength[ret1]);
	if(ret2) Check(Checkarmlength[ret2]);
	if(ret3) Check(Checkarmlength[ret3]);
	
	if(ret1) retardedtime -= Storedarmlength[ret1][tIndex];
	if(ret2) retardedtime -= Storedarmlength[ret2][tIndex];
	if(ret3) retardedtime -= Storedarmlength[ret3][tIndex];
  */
#ifdef DEBUG    
  Check(Checkputn[link]);
#endif
  Vector linkn;
  linkn = Storedputn[link][tIndex];
  //    lisa->putn(linkn,link,t);  
    
  double denom = linkn.dotproduct(wave->k);
  if( (link == 3 && recv == 1) || (link == 2 && recv == 3) || (link == 1 && recv == 2))
    denom = 1.0 - denom;
  else
    denom = 1.0 + denom;

//    return (( psi(link, retardedtime + retard(send, t) - lisa->armlength(link,t)) -
//              psi(link, retardedtime + retard(recv, t)) ) / denom );
    
#ifdef DEBUG    
    Check(Checkarmlength[link]);
#endif
    //    return (( psi(link, retardedtime - Storedarmlength[link][tIndex], retardedtime + retard(send, tIndex) - Storedarmlength[link][tIndex]) -
    //        psi(link, retardedtime, retardedtime + retard(recv, tIndex)) ) / denom );
    //

    return (( psi( StoredputnRet[SPNRI(ret1,ret2,ret3,link,1)][tIndex], retardedtime + retardfast(send, tIndex) - Storedarmlength[link][tIndex]) -
	      psi( StoredputnRet[SPNRI(ret1,ret2,ret3,link,0)][tIndex], retardedtime + retardfast(recv, tIndex)) ) / denom );
    
}


int TDIfast::CacheX(void) {
  
  return( yc(1, 3, 2, 3, 2, 2) -
	  yc(1, 2, 3, 2, 3, 3) +
	  yc(2, 3, 1, 2, 2, 0) -
	  yc(3, 2, 1, 3, 3, 0) +
	  yc(1, 2, 3, 2, 0, 0) -
	  yc(1, 3, 2, 3, 0, 0) + 
	  yc(3, 2, 1, 0, 0, 0) -
	  yc(2, 3, 1, 0, 0, 0) );
}

int TDIfast::CacheY(void) {
        
  return( yc(2, 1, 3, 1, 3, 3) -
	  yc(2, 3, 1, 3, 1, 1) +
	  yc(3, 1, 2, 3, 3, 0) -
	  yc(1, 3, 2, 1, 1, 0) +
	  yc(2, 3, 1, 3, 0, 0) -
	  yc(2, 1, 3, 1, 0, 0) + 
	  yc(1, 3, 2, 0, 0, 0) -
	  yc(3, 1, 2, 0, 0, 0) );
}

int TDIfast::CacheZ(void) {
  
  return( yc(3, 2, 1, 2, 1, 1) -
	  yc(3, 1, 2, 1, 2, 2) +
	  yc(1, 2, 3, 1, 1, 0) -
	  yc(2, 1, 3, 2, 2, 0) +
	  yc(3, 1, 2, 1, 0, 0) -
	  yc(3, 2, 1, 2, 0, 0) + 
	  yc(2, 1, 3, 0, 0, 0) -
	  yc(1, 2, 3, 0, 0, 0) );
}


/*
double TDIfast::M(double time) {
    double t = time * 3.17098E-8;
    
    return( y(1, 2, 3, 2, 0, 0, t) -
            y(1, 3, 2, 3, 0, 0, t) + 
            y(3, 2, 1, 0, 0, 0, t) -
            y(2, 3, 1, 0, 0, 0, t) );
}

double TDIfast::N(double time) {
    double t = time * 3.17098E-8;
        
    return( y(2, 3, 1, 3, 0, 0, t) -
            y(2, 1, 3, 1, 0, 0, t) + 
            y(1, 3, 2, 0, 0, 0, t) -
            y(3, 1, 2, 0, 0, 0, t) );
}

double TDIfast::O(double time) {
    double t = time * 3.17098E-8;

    return( y(3, 1, 2, 1, 0, 0, t) -
            y(3, 2, 1, 2, 0, 0, t) + 
            y(2, 1, 3, 0, 0, 0, t) -
            y(1, 2, 3, 0, 0, 0, t) );
}
*/

double TDIfast::Xfast(int tIndex) {
  return( y(1, 3, 2, 3, 2, 2, tIndex) -
	  y(1, 2, 3, 2, 3, 3, tIndex) +
	  y(2, 3, 1, 2, 2, 0, tIndex) -
	  y(3, 2, 1, 3, 3, 0, tIndex) +
	  y(1, 2, 3, 2, 0, 0, tIndex) -
	  y(1, 3, 2, 3, 0, 0, tIndex) + 
	  y(3, 2, 1, 0, 0, 0, tIndex) -
	  y(2, 3, 1, 0, 0, 0, tIndex) );
}

double TDIfast::Yfast(int tIndex) {
  return( y(2, 1, 3, 1, 3, 3, tIndex) -
	  y(2, 3, 1, 3, 1, 1, tIndex) +
	  y(3, 1, 2, 3, 3, 0, tIndex) -
	  y(1, 3, 2, 1, 1, 0, tIndex) +
	  y(2, 3, 1, 3, 0, 0, tIndex) -
	  y(2, 1, 3, 1, 0, 0, tIndex) + 
	  y(1, 3, 2, 0, 0, 0, tIndex) -
	  y(3, 1, 2, 0, 0, 0, tIndex) );
}

double TDIfast::Zfast(int tIndex) {
  return( y(3, 2, 1, 2, 1, 1, tIndex) -
	  y(3, 1, 2, 1, 2, 2, tIndex) +
	  y(1, 2, 3, 1, 1, 0, tIndex) -
	  y(2, 1, 3, 2, 2, 0, tIndex) +
	  y(3, 1, 2, 1, 0, 0, tIndex) -
	  y(3, 2, 1, 2, 0, 0, tIndex) + 
	  y(2, 1, 3, 0, 0, 0, tIndex) -
	  y(1, 2, 3, 0, 0, 0, tIndex) );
}

/*
double TDIfast::alpha(double time) {
    double t = time * 3.17098E-8;
        
    return( y(3, 2, 1, 0, 0, 0, t) -
            y(2, 3, 1, 0, 0, 0, t) +
            y(2, 1, 3, 2, 0, 0, t) -
            y(3, 1, 2, 3, 0, 0, t) +
            y(1, 3, 2, 1, 2, 0, t) -
            y(1, 2, 3, 1, 3, 0, t)  );
}

double TDIfast::beta(double time) {
    double t = time * 3.17098E-8;
        
    return( y(1, 3, 2, 0, 0, 0, t) -
            y(3, 1, 2, 0, 0, 0, t) +
            y(3, 2, 1, 3, 0, 0, t) -
            y(1, 2, 3, 1, 0, 0, t) +
            y(2, 1, 3, 2, 3, 0, t) -
            y(2, 3, 1, 2, 1, 0, t)  );
}

double TDIfast::gamma(double time) {
    double t = time * 3.17098E-8;
        
    return( y(2, 1, 3, 0, 0, 0, t) -
            y(1, 2, 3, 0, 0, 0, t) +
            y(1, 3, 2, 1, 0, 0, t) -
            y(2, 3, 1, 2, 0, 0, t) +
            y(3, 2, 1, 3, 1, 0, t) -
            y(3, 1, 2, 3, 2, 0, t)  );
}

double TDIfast::zeta(double time) {
    double t = time * 3.17098E-8;
        
    return( y(1, 3, 2, 2, 0, 0, t) -
            y(1, 2, 3, 3, 0, 0, t) +
            y(2, 1, 3, 3, 0, 0, t) -
            y(2, 3, 1, 1, 0, 0, t) +
            y(3, 2, 1, 1, 0, 0, t) -
            y(3, 1, 2, 2, 0, 0, t)  );
}
*/
