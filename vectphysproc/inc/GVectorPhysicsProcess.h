#ifndef GVECTORPHYSICSPROCESS_H
#define GVECTORPHYSICSPROCESS_H

#include "Geant/Config.h"
#include "Geant/Typedefs.h"

#ifndef GEANT_PHYSICSPROCESS
#include "PhysicsProcess.h"
#endif

#include "GeantFwd.h"

// Interface to vector physics models
class GVComptonProcess;

class GVectorPhysicsProcess : public PhysicsProcess
{
private:
  GVComptonProcess    **fVComptonProcess; 
  double                fEnergyLimit;       // tracking cut in kinetic energy [GeV]
  int                   fNumThreads;        // number of working threads

public:
  using GeantTrack_v = Geant::GeantTrack_v;

  GVectorPhysicsProcess();
  GVectorPhysicsProcess(double energyLimit, int numThreads);
  virtual ~GVectorPhysicsProcess();

  virtual void Initialize();
  virtual void PostStepFinalStateSampling( Material_t* /*mat*/,
                                           Int_t ntracks, 
                                           GeantTrack_v &tracks,
                                           Int_t &nout, 
                                           Int_t tid);

  // these are not active !!! 
  //
  virtual void ComputeIntLen(Material_t * /*mat*/,
                             Int_t /*ntracks*/, 
                             GeantTrack_v & /*tracks*/,
                             Double_t * /*lengths*/, 
                             Int_t /*tid*/)                            {}
  
  virtual void PostStep(     Material_t * /*mat*/,
                             Int_t /*ntracks*/,
                             GeantTrack_v &/*tracks*/, 
                             Int_t & /*nout*/, 
                             Int_t /*tid*/)                            {}         

  virtual void PostStepTypeOfIntrActSampling(     Material_t * /*mat*/,
                                                  Int_t /*ntracks*/,
                                                  GeantTrack_v & /*tracks*/, 
                                                  Int_t /*tid*/)       {} 

  virtual void AtRest(       Int_t /*ntracks*/,
                             GeantTrack_v &/*tracks*/, 
                             Int_t &/*nout*/, 
                             Int_t /*tid*/)                             {}

  virtual void Eloss(        Material_t */*mat*/,
                             Int_t /*ntracks*/,
                             GeantTrack_v &/*tracks*/,
                             Int_t &/*nout*/,
                             Int_t /*tid*/)                             {}

  virtual void ApplyMsc(     Material_t */*mat*/,
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

