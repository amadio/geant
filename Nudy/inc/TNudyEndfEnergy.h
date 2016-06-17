#ifndef TNudyEndfEnergy_H
#define TNudyEndfEnergy_H

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

class  TNudyEndfEnergy : public TNudyEndfRecoPoint {

public: 
  TNudyEndfEnergy ();
  TNudyEndfEnergy (TNudyEndfFile *file);
  virtual ~TNudyEndfEnergy ();
  double cmToLabElasticE(double inE, double cmCos, double awr);
  double cmToLabElasticCosT(double cmCos, double awr);
  double cmToLabInelasticE(double cmEOut, double inE, double cmCos, double awr);
  double cmToLabInelasticCosT(double labEOut, double cmEOut, double inE, double cmCos, double awr);
  std::vector<std::vector<std::vector<double> > >cosPdf4OfMts;        // cosine and pdf from file 4 for each reaction
  std::vector<std::vector<std::vector<double> > >cosCdf4OfMts;        // cosine and cdf from file 4 for each reaction
  std::vector<std::vector<double> > energy4OfMts;       // incident energy in file 4 for each reaction
  std::vector<std::vector<int> > Mt4Values;             // MT values for which angular distributions are given in file 4
  std::vector<int> MtNumbers;				// MT numbers
private:
  double recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2);
  double recursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy);
  double recursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy);
  double recursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy);
  void cdfGenerateE(std::vector<double> &x1,std::vector<double> &x2);
  void Sort(std::vector<double> &x1,std::vector<double> &x2);
  double Thinning(std::vector<double> &x1, std::vector<double> &x2);
  double ThinningDuplicate(std::vector<double> &x1);
  double ThinningDuplicate(std::vector<double> &x1,std::vector<double> &x2);

  double A, AWR, ABN, QX;                         // standard ENDF parameters
  int *NBT,*NBT1,*INT1, NR, NP;                         // standard ENDF parameters for range and interpolation
  int *NBT2, *INT2, NR2, NE;                            // standard ENDF parameters for range and interpolation
  int *NBT3, *INT3, NR3, NE2;				// standard ENDF parameters for range and interpolation
  double *INorm, QValue[999];				// ENDF parameter and Q values from file 3
  double sigDiff;					// precision for cross-section reconstruction
  double *fE1_file5,*fp91_file5,*fE2_file5,*fp92_file5,
         *fE3_file5,*fp93_file5; 			// file5 energy and probabilities
  TArrayD *lCoef, xengr;
  std::vector<double> energyFile5;
  std::vector<double> energyPdfFile5;
  std::vector<double> energyCdfFile5;
  std::vector<double> energy, sigma, sigmaT; 
  std::vector<double> eneTemp,sigTemp;			// temporary vectors to store energy and sigma
  ClassDef(TNudyEndfEnergy, 1) // class for an ENDF reconstruction
};
#endif