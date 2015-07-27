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
    int fLen; // Length of data
    TNudyAlias *fChooseBin; // Use alias method to choose bin
    //Probability distribution
    TArrayD *fP; // Probability distribution
    TArrayD *fX; // Energy distribution
    TArrayD *fInterX; // Built Energy distribution
    TArrayD *fInterP; // Built Probability distribution
    TArrayD *fTx; //! Temporary Energy distribution
    TArrayD *fTp; //! Temporary Probability distribution
    double fTa; //! Temporary Alpha value
    double fInterAlpha; // Alpha for the next distribution
    TRandom *fRan; // To generate unifrom random numbers
    TRandom *fRnd; // To generate unifrom random numbers
    double fAlpha;// Stores the number which identifies the distribution
  public:
    TNudyAliasCont();
    TNudyAliasCont(double *p,double *x,const int len,double alpha,unsigned int seed=65539);
    TNudyAliasCont(TArrayD *data,double alpha,unsigned int seed);
    void Initialize(double *p,double *x,const int len,double alpha,unsigned int seed=65539);
    virtual ~TNudyAliasCont();
    void DumpTable();
    double GetP(int i) {if(i>=0 && i<fLen) return fP->GetAt(i); else return -1;}
    double GetX(int i) {if(i>=0 && i<fLen) return fX->GetAt(i); else return -1;}
    double GetAlpha() {return fAlpha;}
    double Uniform(double x=1) {if(fRnd) return fRnd->Uniform(x); else return -1;}
    double RandomBin() {if(fChooseBin) return fChooseBin->Random(); else return -1;}
    double Random(IntScheme_t iScheme=kLinear,int binNo=-1,AliasDist_t arr=kOriginal,double r1 = -1,double r2 = -1);
    double ImprovedInterpolation(double alpha);
    double SelectBin(double ein, AliasDist_t distribution=kOriginal);
    static void BuildIntermediate(TNudyAliasCont *x,const int len);
    static double StatisticalInterpolation(TNudyAliasCont *dist1, TNudyAliasCont *dist2,double alpha);
#ifdef TNUDYALIAS_MULTITHREAD
    double *fMult; //! Random numbers to be generated in parallel
    int fMultLen; //! Number if random numbres to be generated in parallel
    double* Randoms(int n,IntScheme_t iScheme=kLinear);
#endif
    
    ClassDef(TNudyAliasCont,1)
};

#endif
