#ifndef GPVEmProcess_h
#define GPVEmProcess_h 1

#include "GPMaterial.h"
#include "GXTrack.h"
//#include "G4ParticleChangeForGamma.hh"

class GPVEmProcess // : public G4VDiscreteProcess
{
public:

  FQUALIFIER GPVEmProcess(curandState* devStates,
			  int threadId);

  FQUALIFIER ~GPVEmProcess();

  FQUALIFIER void Clear();

  FQUALIFIER void InitialiseProcess( /* G4VEmModel* model */ );

  FQUALIFIER void StartTracking(GXTrack* track);

  FQUALIFIER G4double PostStepGetPhysicalInteractionLength(
                             GXTrack* track,
                             GPMaterial* material,
                             G4double   previousStepSize,
                             GPForceCondition* condition);

  FQUALIFIER GPVParticleChange* PostStepDoIt(GXTrack* track, 
					     GPMaterial* material);

  FQUALIFIER G4double CrossSectionPerVolume(G4double kineticEnergy,
					    GPMaterial* material);

  FQUALIFIER G4double GetMeanFreePath(GPMaterial* material,
				      G4double kineticEnergy,
				      GPForceCondition* condition);

  FQUALIFIER G4double MeanFreePath(GPMaterial* material, G4double kineticEnergy);

  FQUALIFIER G4double ComputeCrossSectionPerAtom(G4double kineticEnergy, 
						 G4double Z, G4double A=0., 
						 G4double cut=0.0);

  FQUALIFIER GPElement* GetCurrentElement();
 
  FQUALIFIER void SetMinKinEnergy(G4double e);

  FQUALIFIER void SetMaxKinEnergy(G4double e);

  FQUALIFIER int GPlrint(double ad);

  FQUALIFIER void DefineMaterial(GPMaterial* material);

  FQUALIFIER G4double GetLambdaFromTable(G4double kinEnergy);

  FQUALIFIER G4double ComputeCurrentLambda(G4double kinEnergy);

  FQUALIFIER G4double GetCurrentLambda(G4double kinEnergy);

  FQUALIFIER void ComputeIntegralLambda(G4double kinEnergy);

  FQUALIFIER G4double GetLambdaFromTablePrim(G4double kinEnergy);

  FQUALIFIER G4double GetLambda(G4double kinEnergy, 
				GPMaterial* material);

  FQUALIFIER G4double MinKinEnergy();

  FQUALIFIER G4double MaxKinEnergy();

  FQUALIFIER void SetPolarAngleLimit(G4double a);

  FQUALIFIER G4double PolarAngleLimit();

  FQUALIFIER void SetLambdaFactor(G4double val);

  //G4VProcess

  FQUALIFIER void ResetNumberOfInteractionLengthLeft();

  FQUALIFIER void EndTracking();

private:
  //curand                  
  int fThreadId;
  curandState* fDevStates;

  GPVEmModel*                  currentModel;
  GPPhysicsTable*              theLambdaTable;
  G4double*                    theDensityFactor;
  G4int*                       theDensityIdx;

  G4double                     minKinEnergy;
  G4double                     maxKinEnergy;
  G4int                        nLambdaBins;
  G4double                     lambdaFactor;
  G4double                     polarAngleLimit;

  GPVParticleChange            fParticleChange;
  GPMaterial*                  currentMaterial;

  G4double                     mfpKinEnergy;
  G4double                     preStepKinEnergy;
  G4double                     preStepLambda;
  G4double                     fFactor;

  G4bool isInitialised;

  //G4VProcess
  G4double theNumberOfInteractionLengthLeft;
  G4double currentInteractionLength;
  G4double theInitialNumberOfInteractionLength;

};

#endif
