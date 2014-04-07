//
// S.Y. Jun & J. Apostolaskis, April 2014
//
#ifndef TabulatedProcess_HH
#define TabulatedProcess_HH 1

#include "G4WrapperProcess.hh"
#include "G4ProcessType.hh"

//tabulated physics 
#include "TFile.h"
#include "GXTrack.hh"

class G4VParticleChange;

class TabulatedProcess : public G4WrapperProcess {
  
public:

  TabulatedProcess(G4String processName, G4ProcessType processType);	
  virtual ~TabulatedProcess();	

  // Override PostStepGetPhysicalInteractionLength method
  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                G4double previousStepSize,
                                                G4ForceCondition* condition);  

  // Override PostStepDoIt method
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);

  void Print(const G4Step& astep);
  void PrepareTable_t(const char* pnam, const char* reac);
private:
  G4VParticleChange* particleChange;

  //tabulated physics data of VP 
  TFile *fxsec;
  TFile *ffsta;

  //this should array for multiple materials
  G4double* mxsec;
  GXTrack* fstrack_h;
};

#endif
