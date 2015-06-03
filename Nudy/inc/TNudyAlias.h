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
    Int_t fLen; //Length of data
    //Alias table
    Double_t *fP; //[fLen]
    Double_t *fX; //[fLen]
    Double_t *fA; //[fLen]
    Double_t *fR; //[fLen]
    TRandom *fRnd; //Uniform random number generation
  public:
    TNudyAlias();
    virtual ~TNudyAlias();
    TNudyAlias(Double_t *p,Double_t *x,const Int_t len,UInt_t seed=65539);
    void DumpTable();
    Double_t Random();
    Double_t GetP(Int_t i) {if(i<fLen) return fP[i]; else return -1;};
    Double_t GetX(Int_t i) {if(i<fLen) return fX[i]; else return -1;};
    Int_t GetLen() {return fLen;}
#ifdef TNUDYALIAS_MULTITHREAD
    class TNudyComStruct {
      //Class to communicate between thread handler and objects
      public:
        TNudyAlias *fAl;
        Int_t fI;
        TNudyComStruct(TNudyAlias *a,Int_t j) {fAl=a;fI=j;}
        virtual ~TNudyComStruct() {}
    };
    static void* ThreadHandle(void *ptr);
    Double_t *fMult; //! When generating using a multithreaded approach
    Int_t fMultLen;//! Number of random values to be generated using the multi threaded approach
    Double_t* Randoms(Int_t n);
#endif

  ClassDef(TNudyAlias,1)
};

#endif
