//
// S.Y. Jun & J. Apostolaskis, April 2014
//
#ifndef TabulatedDataManager_HH
#define TabulatedDataManager_HH 1

#include "globals.hh"
#include "GXTrack.hh"
#include "TFile.h"

class TabulatedDataManager
{
public:

  TabulatedDataManager();	
  virtual ~TabulatedDataManager();	

  void PrepareTable(const char* pnam, const char* reac);
private:

  //tabulated physics data of VP 
  TFile *fxsec;
  TFile *ffsta;

  //this should array for multiple materials
  G4double* mxsec;
  GXTrack* fstrack_h;
};

#endif
