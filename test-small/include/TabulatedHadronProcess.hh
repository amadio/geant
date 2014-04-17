//
// S.Y. Jun & J. Apostolaskis, April 2014
//
#ifndef TabulatedHadronProcess_HH
#define TabulatedHadronProcess_HH 1

#include "G4WrapperProcess.hh"
#include "G4ProcessType.hh"

class TabulatedDataManager;
class G4VParticleChange;

class TabulatedHadronProcess : public G4WrapperProcess {
  
public:

  TabulatedHadronProcess(G4String processName, G4ProcessType processType);	
  virtual ~TabulatedHadronProcess();	

  // Override PostStepGetPhysicalInteractionLength method
  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                G4double previousStepSize,
                                                G4ForceCondition* condition);  

  // Override PostStepDoIt method
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);

  void Print(const G4Step& astep);
private:
  G4VParticleChange* particleChange;
  TabulatedDataManager* theDataManager;
};

#endif
