//
// S.Y. Jun & J. Apostolaskis, April 2014
//
#ifndef TabulatedProcess_HH
#define TabulatedProcess_HH 1

#include "G4WrapperProcess.hh"
#include "G4ProcessType.hh"
#include "GXTrack.hh"

#include "TFile.h"

//class TabulatedDataManager;
class G4VParticleChange;

class TabulatedProcess : public G4WrapperProcess {
  
public:

  TabulatedProcess(G4String processName, G4ProcessType processType);	
  virtual ~TabulatedProcess();	

  // Override PostStepGetPhysicalInteractionLength method
  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                G4double previousStepSize,
                                                G4ForceCondition* condition);  
  G4double MeanFreePath(const G4Track& track);

  // Override PostStepDoIt method
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);

  void Print(const G4Step& astep);
private:
  void PrepareTable(const char* pnam, const char* reac);

private:
  G4VParticleChange* particleChange;

  //tabulated physics data of VP - move to PhysicsList
  //  TabulatedDataManager* dataManager;
  static TFile *fxsec;
  static TFile *ffsta;

  G4double* fXsecTable;  
  GXTrack* fFinalState;

};

#endif
