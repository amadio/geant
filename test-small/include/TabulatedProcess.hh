//
// S.Y. Jun & J. Apostolaskis, April 2014
//
#ifndef TabulatedProcess_HH
#define TabulatedProcess_HH 1

#include "G4WrapperProcess.hh"
#include "G4ProcessType.hh"
#include "GXTrack.hh"

#include "TPartIndex.h"

class TabulatedDataManager;
class G4VParticleChange;
class G4DynamicParticle;

class TabulatedProcess : public G4WrapperProcess {
  
public:

  TabulatedProcess(G4String procName, G4ProcessType procType);	
  TabulatedProcess(G4String procName, G4ProcessType procType, 
		   G5proc indexReac);	
  virtual ~TabulatedProcess();	

  // Override PostStepGetPhysicalInteractionLength method
  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                G4double previousStepSize,
                                                G4ForceCondition* condition);  
  G4double SetupForMaterial(const G4Track& track);
  G4double MeanFreePath(const G4Track& track);

  // Override PostStepDoIt method
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);

  void ResetSecondaryDefinition(G4int pid);
  void Print(const G4Step& astep);

private:
  G4int fMaterialIndex;
  G5proc fReaction;
  G4VParticleChange* fParticleChange;
  G4ParticleDefinition* fSecDefinition;
  TabulatedDataManager* theDataManager;
  std::vector<GXTrack*>  fSecParticles;
};

#endif
