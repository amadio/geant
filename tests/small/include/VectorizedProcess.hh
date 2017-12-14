//
// S.Y. Jun & J. Apostolaskis, April 2014
//
#ifndef VectorizedProcess_HH
#define VectorizedProcess_HH 1

#include "G4WrapperProcess.hh"
#include "G4ProcessType.hh"

class G4VParticleChange;

class VectorizedProcess : public G4WrapperProcess {
  
public:

  VectorizedProcess(G4String processName, G4ProcessType processType);	
  virtual ~VectorizedProcess();	
  
  // Override PostStepDoIt  method
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);

  void Print(const G4Step& astep);

private:
  VectorizedProcess(const VectorizedProcess&); // Not implemented
  VectorizedProcess& operator=(const VectorizedProcess&); // Not implemented
  G4VParticleChange* particleChange;
};

#endif
