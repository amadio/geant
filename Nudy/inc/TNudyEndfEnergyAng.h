#ifndef TNudyEndfEnergyAng_H
#define TNudyEndfEnergyAng_H

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

#define PI 2.0 * asin(1.0)

class  TNudyEndfEnergyAng : public TNudyEndfRecoPoint {

public: 
  TNudyEndfEnergyAng ();
  TNudyEndfEnergyAng (TNudyEndfFile *file);
  virtual ~TNudyEndfEnergyAng ();
  std::vector<std::vector<std::vector<double> > >energyPdf5OfMts;        // cosine, energy and pdf from file 5 for each reaction
  std::vector<std::vector<std::vector<double> > >energyCdf5OfMts;        // cosine, energy and cdf from file 5 for each reaction
  std::vector<std::vector<double> > energy5OfMts;       // incident energy in file 5 for each reaction
  std::vector<std::vector<int> > Mt4Values;             // MT values for which angular distributions are given in file 4
  std::vector<int> MtNumbers;				// MT numbers
  std::vector<std::vector<int> > Mt4Lct;                // CM and Lab flag for angular distributions as given in file 4
  std::vector<int> MtLct;				// LCT numbers
private:
  double recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2);
  double recursionLinearFile4(int i, double x1, double x2, double pdf1, double pdf2);
  double recursionLinearProb(double x1, double x2, double pdf1, double pdf2);

  double A, AWR, ABN, QX;                         // standard ENDF parameters
  int NR, NP;                         // standard ENDF parameters for range and interpolation
  int NR2, NE;                            // standard ENDF parameters for range and interpolation
  int NR3, NE2;				// standard ENDF parameters for range and interpolation
//  double *INorm, QValue[999];				// ENDF parameter and Q values from file 3
  std::vector<double> cosFile4, energyFile5;
  std::vector<double> cosPdfFile4, energyPdfFile5;
  std::vector<double> cosCdfFile4, energyCdfFile5;
  std::vector<double> edes6,f06,r6,a6;
  std::vector<double> fE1, fP1, fE2, fP2, fE3, fP3, INorm;
  std::vector<int> nbt1,int1;
  int nr1, np1;                         // standard ENDF parameters
  std::vector<int> nbt2,int2;
  int nr2, np2;                         // standard ENDF parameters
  std::vector<int> nbt3,int3;
  int nr3, np3;                         // standard ENDF parameters
  std::vector<double>ein,cdf,pdf,lCoef1;
  std::vector<std::vector<double> >pdf2d,cdf2d,lCoef;
  ClassDef(TNudyEndfEnergyAng, 1) // class for an ENDF reconstruction
};
#endif