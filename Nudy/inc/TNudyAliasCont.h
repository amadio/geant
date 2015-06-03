#ifndef ROOT_TNudyAliasCont
#define ROOT_TNudyAliasCont

#include <TObject.h>
#include <TRandom.h>
#include <TArrayD.h>

#include "TNudyAlias.h"
#include "TNudyCore.h"

//--------------------------------------------------
// This class takes a continuous piecewise linear probability
// distribution and generates a random number according to that
// probability distributions using the Alias method
//-------------------------------------------------- 

class TNudyAliasCont : public TObject {
  private:
    Int_t fLen; // Length of data
    TNudyAlias *fChooseBin; // Use alias method to choose bin
    //Probability distribution
    TArrayD *fP; // Probability distribution
    TArrayD *fX; // Energy distribution
    TArrayD *fInterX; // Built Energy distribution
    TArrayD *fInterP; // Built Probability distribution
    TArrayD *fTx; //! Temporary Energy distribution
    TArrayD *fTp; //! Temporary Probability distribution
    Double_t fTa; //! Temporary Alpha value
    Double_t fInterAlpha; // Alpha for the next distribution
    TRandom *fRan; // To generate unifrom random numbers
    TRandom *fRnd; // To generate unifrom random numbers
    Double_t fAlpha;// Stores the number which identifies the distribution
  public:
    TNudyAliasCont();
    TNudyAliasCont(Double_t *p,Double_t *x,const Int_t len,Double_t alpha,UInt_t seed=65539);
    TNudyAliasCont(TArrayD *data,Double_t alpha,UInt_t seed);
    void Initialize(Double_t *p,Double_t *x,const Int_t len,Double_t alpha,UInt_t seed=65539);
    virtual ~TNudyAliasCont();
    void DumpTable();
    Double_t GetP(Int_t i) {if(i>=0 && i<fLen) return fP->GetAt(i); else return -1;}
    Double_t GetX(Int_t i) {if(i>=0 && i<fLen) return fX->GetAt(i); else return -1;}
    Double_t GetAlpha() {return fAlpha;}
    Double_t Uniform(Double_t x=1) {if(fRnd) return fRnd->Uniform(x); else return -1;}
    Double_t RandomBin() {if(fChooseBin) return fChooseBin->Random(); else return -1;}
    Double_t Random(IntScheme_t iScheme=kLinear,Int_t binNo=-1,AliasDist_t arr=kOriginal,Double_t r1 = -1,Double_t r2 = -1);
    Double_t ImprovedInterpolation(Double_t alpha);
    Double_t SelectBin(Double_t ein, AliasDist_t distribution=kOriginal);
    static void BuildIntermediate(TNudyAliasCont *x,const Int_t len);
    static Double_t StatisticalInterpolation(TNudyAliasCont *dist1, TNudyAliasCont *dist2,Double_t alpha);
#ifdef TNUDYALIAS_MULTITHREAD
    Double_t *fMult; //! Random numbers to be generated in parallel
    Int_t fMultLen; //! Number if random numbres to be generated in parallel
    Double_t* Randoms(Int_t n,IntScheme_t iScheme=kLinear);
#endif
    
    ClassDef(TNudyAliasCont,1)
};

#endif
