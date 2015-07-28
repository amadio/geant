#ifndef ROOT_TNudyAlias
#define ROOT_TNudyAlias

#include <TObject.h>
#include <TRandom.h>

#include "TNudyCore.h"

//--------------------------------------------------
// Provides random number generation using a discrete
// probability distribution using the Alias Method
//--------------------------------------------------

class TNudyAlias : public TObject {
private:
  int fLen; // Length of data
  // Alias table
  double *fP;    //[fLen]
  double *fX;    //[fLen]
  double *fA;    //[fLen]
  double *fR;    //[fLen]
  TRandom *fRnd; // Uniform random number generation
public:
  TNudyAlias();
  virtual ~TNudyAlias();
  TNudyAlias(double *p, double *x, const int len, unsigned int seed = 65539);
  void DumpTable();
  double Random();
  double GetP(int i) {
    if (i < fLen)
      return fP[i];
    else
      return -1;
  };
  double GetX(int i) {
    if (i < fLen)
      return fX[i];
    else
      return -1;
  };
  int GetLen() { return fLen; }
#ifdef TNUDYALIAS_MULTITHREAD
  class TNudyComStruct {
    // Class to communicate between thread handler and objects
  public:
    TNudyAlias *fAl;
    int fI;
    TNudyComStruct(TNudyAlias *a, int j) {
      fAl = a;
      fI = j;
    }
    virtual ~TNudyComStruct() {}
  };
  static void *ThreadHandle(void *ptr);
  double *fMult; //! When generating using a multithreaded approach
  int fMultLen;  //! Number of random values to be generated using the multi threaded approach
  double *Randoms(int n);
#endif

  ClassDef(TNudyAlias, 1)
};

#endif
