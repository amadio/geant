#ifndef TNudyEndfAng_H
#define TNudyEndfAng_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iterator>
#include <map>
#include <algorithm>
#include <iomanip>
#include <dirent.h>
#include "TRandom3.h"
#include "Math/SpecFuncMathMore.h"
#include <TCollection.h>
#include <TIterator.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TClass.h>
#include <TObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <Rtypes.h>
#include <TTreeReader.h>
#include <TKey.h>

#include "TNudyCore.h"
#include "TNudyDB.h"
#include "TNudyElementTable.h"
#include "TNudyENDF.h"
#include "TNudyEndfCont.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfINTG.h"
#include "TNudyEndfList.h"
#include "TNudyEndfMat.h"
#include "TNudyEndfTape.h"
#include "TNudyEndfRecord.h"
#include "TNudyEndfSec.h"
#include "TNudyEndfTab1.h"
#include "TNudyEndfTab2.h"
#include "TNudyLibrary.h"
#include "TNudyManager.h"
#include "TNudySubLibrary.h"
#include "TVNudyModel.h"

#include <Riostream.h>
#include "TNudyEndfRecoPoint.h"
#include <RConfig.h>

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
  ClassDef(TNudyEndfAng, 1) // class for an ENDF reconstruction
};
#endif