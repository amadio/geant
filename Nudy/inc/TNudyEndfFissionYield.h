#ifndef TNudyEndfFissionYield_H
#define TNudyEndfFissionYield_H

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

class  TNudyEndfFissionYield : public TNudyEndfRecoPoint {

public: 
  TNudyEndfFissionYield ();
  TNudyEndfFissionYield (TNudyEndfFile *file);
  virtual ~TNudyEndfFissionYield ();
private:

  double A, AWR, ABN, QX;                                               // standard ENDF parameters
  std::vector<double> ein, einc;				        // incident energy
  std::vector<std::vector<double> >zafp,fps,zafpc,fpsc,yi, dyi, yc, dyc;// charge, mass, yield (independent and cummulative)
  std::vector<double> zafp1,fps1,zafpc1,fpsc1,yi1, dyi1, yc1, dyc1;     // charge, mass, yield (independent and cummulative)
  ClassDef(TNudyEndfFissionYield, 1)                                    // class for an ENDF fission yield reconstruction
};
#endif