#ifndef SteppingAnalysis_H
#define SteppingAnalysis_H 1

#ifdef USE_ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TBranch.h"
#endif

class G4Step;
class HistogramManager;

class SteppingAnalysis 
{
public:
  SteppingAnalysis();
  ~SteppingAnalysis();

public:
  void DoIt(const G4Step* aStep);

public:

  void FillCrossSections( const G4Step* aStep );

//root
private:
  SteppingAnalysis(const SteppingAnalysis&); // Not implemented
  SteppingAnalysis& operator=(const SteppingAnalysis&); // Not implemented

  //histogram Manager
  HistogramManager* theHisto;
 
};

#endif


