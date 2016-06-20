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
  std::vector<std::vector<std::vector<double> > >energyPdf5OfMts;        // energy and pdf from file 5 for each reaction
  std::vector<std::vector<std::vector<double> > >energyCdf5OfMts;        // energy and cdf from file 5 for each reaction
  std::vector<std::vector<double> > energy5OfMts;       // incident energy in file 5 for each reaction
  std::vector<std::vector<int> > Mt5Values;             // MT values for which energy distributions are given in file 5
  std::vector<int> MtNumbers;				// MT numbers temporary
private:
  double recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2);
  double recursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy);
  double recursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy);
  double recursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy);

  double A, AWR, ABN, QX;                         // standard ENDF parameters
  int NR, NP;                         // standard ENDF parameters for range and interpolation
  double QValue[999];				// ENDF parameter and Q values from file 3
  double sigDiff;					// precision for cross-section reconstruction
  std::vector<double> fE1, fP1, fE2, fP2, fE3, fP3, INorm;
  std::vector<int> nbt1,int1;
  int nr1, np1;                         // standard ENDF parameters
  std::vector<int> nbt2,int2;
  int nr2, np2;                         // standard ENDF parameters
  std::vector<int> nbt3,int3;
  int nr3, np3;                         // standard ENDF parameters
  std::vector<double> energyFile5;
  std::vector<double> energyPdfFile5;
  std::vector<double> energyCdfFile5;
  std::vector<double>ein, cdf, pdf;
  std::vector<std::vector<double> >cdf2d, pdf2d;
  ClassDef(TNudyEndfEnergy, 1) // class for an ENDF energy reconstruction
};
#endif