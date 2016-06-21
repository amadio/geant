#ifndef TNudyEndfAng_H
#define TNudyEndfAng_H

#include "TNudyEndfRecoPoint.h"

class  TNudyEndfAng : public TNudyEndfRecoPoint {

public: 
  TNudyEndfAng ();
  TNudyEndfAng (TNudyEndfFile *file);
  virtual ~TNudyEndfAng ();
  std::vector<std::vector<std::vector<double> > >cosPdf4OfMts;        // cosine and pdf from file 4 for each reaction
  std::vector<std::vector<std::vector<double> > >cosCdf4OfMts;        // cosine and cdf from file 4 for each reaction
  std::vector<std::vector<double> > energy4OfMts;       // incident energy in file 4 for each reaction
  std::vector<std::vector<int> > Mt4Values;             // MT values for which angular distributions are given in file 4
  std::vector<std::vector<int> > Mt4Lct;                // CM and Lab flag for angular distributions as given in file 4
  std::vector<int> MtNumbers;				// MT numbers
  std::vector<int> MtLct;				// LCT numbers
private:
  double recursionLinearLeg(int i, double x1, double x2, double pdf1, double pdf2);
  double recursionLinearProb(double x1, double x2, double pdf1, double pdf2);
  double A, AWR, ABN, QX;                         // standard ENDF parameters
  std::vector<double>ein,cdf,pdf,lCoef1;
  std::vector<std::vector<double> >pdf2d,cdf2d,lCoef;
  std::vector<double> cosFile4;
  std::vector<double> cosPdfFile4;
  std::vector<double> cosCdfFile4;
  std::vector<int> nbt1,int1;
  int nr, np;                         // standard ENDF parameters
#ifdef USE_ROOT
  ClassDef(TNudyEndfAng, 1) // class for an ENDF reconstruction
#endif
};
#endif
