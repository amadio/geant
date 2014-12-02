#ifndef GVECTORPHYSICSPROCESS_H
#define GVECTORPHYSICSPROCESS_H

#include "Geant/Config.h"

#ifndef GEANT_PHYSICSPROCESS
#include "PhysicsProcess.h"
#endif

// GeantV related
class TGeoMaterial; // ROOT: we won't use this
class GeantTrack_v;

// Interface to vector physics models
class GVComptonProcess;

class GVectorPhysicsProcess : public PhysicsProcess
{
private:
  GVComptonProcess    **fVComptonProcess; 
  double                fEnergyLimit;       // tracking cut in kinetic energy [GeV]
  int                   fNumThreads;        // number of working threads

public:
  GVectorPhysicsProcess();
  GVectorPhysicsProcess(double energyLimit, int numThreads);
  virtual ~GVectorPhysicsProcess();

  virtual void Initialize();
  virtual void PostStepFinalStateSampling( TGeoMaterial* /*mat*/,
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

   // we need this while vecprot_v2/inc/PhysicsProcess is derived from TNamed
   ClassDef(GVectorPhysicsProcess,1)
};

#endif

