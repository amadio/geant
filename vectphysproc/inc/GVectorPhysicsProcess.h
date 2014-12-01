#ifndef GVECTORPHYSICSPROCESS_H
#define GVECTORPHYSICSPROCESS_H

#include "Geant/Config.h"

#ifndef GEANT_PHYSICSPROCESS
#include "PhysicsProcess.h"
#endif

// Geant-V related
class TGeoMaterial;
class GeantTrack_v;

// Interface to vector physics models
class GVComptonProcess;

//______________________________________________________________________________
class GVectorPhysicsProcess : public PhysicsProcess
{
private:
  GVComptonProcess    *fVComptonProcess; 
  double               fEnergyLimit;       // tracking cut in kinetic energy [GeV]

public:
  GVectorPhysicsProcess();
  GVectorPhysicsProcess(double energyLimit);
  virtual ~GVectorPhysicsProcess();

  virtual void Initialize();
  virtual void PostStepFinalStateSampling( TGeoMaterial *mat,
                                           Int_t ntracks, 
                                           GeantTrack_v &tracks,
                                           Int_t &nout, 
                                           Int_t tid);

  // these are not active !!! 
  //
  virtual void ComputeIntLen(TGeoMaterial * /*mat*/,
                             Int_t /*ntracks*/, 
                             GeantTrack_v & /*tracks*/,
                             Double_t * /*lengths*/, 
                             Int_t /*tid*/)                            {}
  
  virtual void PostStep(     TGeoMaterial * /*mat*/,
                             Int_t /*ntracks*/,
                             GeantTrack_v &/*tracks*/, 
                             Int_t & /*nout*/, 
                             Int_t /*tid*/)                            {}         

  virtual void PostStepTypeOfIntrActSampling(     TGeoMaterial * /*mat*/,
                                                  Int_t /*ntracks*/,
                                                  GeantTrack_v & /*tracks*/, 
                                                  Int_t /*tid*/)       {} 

  virtual void AtRest(       Int_t /*ntracks*/,
                             GeantTrack_v &/*tracks*/, 
                             Int_t &/*nout*/, 
                             Int_t /*tid*/)                             {}

  virtual void Eloss(        TGeoMaterial */*mat*/,
                             Int_t /*ntracks*/,
                             GeantTrack_v &/*tracks*/,
                             Int_t &/*nout*/,
                             Int_t /*tid*/)                             {}

  virtual void ApplyMsc(     TGeoMaterial */*mat*/,
                             Int_t /*ntracks*/,
                             GeantTrack_v &/*tracks*/,
                             Int_t /*tid*/)                             {}

private:
   GVectorPhysicsProcess (const GVectorPhysicsProcess  &);//no imp.	
   GVectorPhysicsProcess & operator=(const GVectorPhysicsProcess  &);//no imp.

ClassDef(GVectorPhysicsProcess,1)
};

#endif

