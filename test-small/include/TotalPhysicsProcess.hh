//
// Authors: M. Nowak, S.Y. Jun & J. Apostolakis, April 2014
//
// Purpose:
//    Single process to replace all (G4) physics processes using tabulations 
//      of x-section and sampled reactions.
// 
// This must be the *only* discrete physics process registered for a particle 
//    (it uses the total cross-section, summed over the G4/real physics processes)

#ifndef TotalPhysicsProcess_HH
#define TotalPhysicsProcess_HH 1

#include "G4VDiscreteProcess.hh"
#include "G4ProcessType.hh"
#include "GXTrack.hh"

class TabulatedDataManager;
class TTabPhysMgr; // Alternative to using TabulatedDataManager

class G4VParticleChange;

class TotalPhysicsProcess : public G4VDiscreteProcess {
  
public:

  TotalPhysicsProcess(G4String processName);
  virtual ~TotalPhysicsProcess();	

  // Override PostStepGetPhysicalInteractionLength method
 /* G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                G4double previousStepSize,
                                                G4ForceCondition* condition);  */
  G4double GetMeanFreePath(const G4Track& track,
                           G4double previousStepSize,
                           G4ForceCondition* condition);

  int      SetupForMaterial(const G4Track& track);
   // For material from the track, find the corresponding Root material
   // Return (and store) its index.

   // Return the products from all processes
  G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);

  void Print(const G4Step& astep);
  
  const G4ParticleDefinition* ParticleDefinition(G4int ipdg);

private:
  G4int                 fMaterialIndex;
  G4int                 fParticleId; 
  G4VParticleChange*    fParticleChange;
  // G4ParticleDefinition* fSecDefinition;
  TabulatedDataManager* theDataManager;
  // TTabPhysMgr *pTTabPhysMgr;

  std::vector<GXTrack*> fSecParticles;
};
#endif
