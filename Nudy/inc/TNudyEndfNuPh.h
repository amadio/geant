#ifndef TNudyEndfNuPh_H
#define TNudyEndfNuPh_H

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

class  TNudyEndfNuPh : public TNudyEndfRecoPoint {

public: 
  TNudyEndfNuPh ();
  TNudyEndfNuPh (TNudyEndfFile *file);
  virtual ~TNudyEndfNuPh ();
private:
  double recursionLinearNuPh(double x1, double x2, double sig1, double sig2, std::vector<double> x, std::vector<double> sig);
  double A, AWR, ABN, QX;                         // standard ENDF parameters
  int NR, NP;                         // standard ENDF parameters for range and interpolation
  std::vector<double> eintFile1,nutFile1,einFile1,nuFile1;
  std::vector<double> eindFile1,nudFile1,einphFile1,phFile1;
  std::vector<double> einfFile1,heatFile1;
  std::vector<double> cnc;
  std::vector<int> nbt1,int1;
  ClassDef(TNudyEndfNuPh, 1) // class for an ENDF reconstruction
};
#endif