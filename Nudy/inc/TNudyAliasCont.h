#ifndef ROOT_TNudyAliasCont
#define ROOT_TNudyAliasCont

#include <TRandom.h>

#include "TNudyTypes.h"
#include "TNudyAlias.h"

//--------------------------------------------------
// This class takes a continuous piecewise linear probability
// distribution and generates a random number according to that
// probability distributions using the Alias method
//--------------------------------------------------

class TNudyAliasCont {
private:
  int fLen;               // Length of data
  TNudyAlias *fChooseBin; // Use alias method to choose bin
  // Probability distribution
  double *fP;         //[fLen] Probability distribution
  double *fX;         //[fLen] Energy distribution
  double *fInterX;    //[fLen] Built Energy distribution
  double *fInterP;    //[fLen] Built Probability distribution
  double *fTx;        //! Temporary Energy distribution
  double *fTp;        //! Temporary Probability distribution
  double fTa;         //! Temporary Alpha value
  double fInterAlpha; // Alpha for the next distribution
  TRandom *fRnd;      // To generate unifrom random numbers
  double fAlpha;      // Stores the number which identifies the distribution
public:
  TNudyAliasCont();
  TNudyAliasCont(double *p, double *x, int len, double alpha, unsigned int seed = 65539);
  TNudyAliasCont(double *data, int len, double alpha, unsigned int seed);
  void Initialize(double *p, double *x, const int len, double alpha, unsigned int seed = 65539);
  virtual ~TNudyAliasCont();
  void DumpTable();
  double GetP(int i) const {
    if (i >= 0 && i < fLen)
       return fP[i];
    else
      return -1;
  }
  double GetX(int i) const {
    if (i >= 0 && i < fLen)
       return fX[i];
    else
      return -1;
  }
  double GetAlpha() const { return fAlpha; }
  double Uniform(double x = 1) {
    if (fRnd)
      return fRnd->Uniform(x);
    else
      return -1;
  }
  double RandomBin() {
    if (fChooseBin)
      return fChooseBin->Random();
    else
      return -1;
  }
  double Random(IntScheme_t iScheme = kLinear, int binNo = -1, AliasDist_t arr = kOriginal, double r1 = -1,
                double r2 = -1);
  double ImprovedInterpolation(double alpha);
  double SelectBin(double ein, AliasDist_t distribution = kOriginal);
  static void BuildIntermediate(TNudyAliasCont *x, const int len);
  static double StatisticalInterpolation(TNudyAliasCont *dist1, TNudyAliasCont *dist2, double alpha);
#ifdef TNUDYALIAS_MULTITHREAD
  double *fMult; //! Random numbers to be generated in parallel
  int fMultLen;  //! Number if random numbres to be generated in parallel
  double *Randoms(int n, IntScheme_t iScheme = kLinear);
#endif

#ifdef USE_ROOT
  ClassDef(TNudyAliasCont, 1)
#endif
};

#endif
