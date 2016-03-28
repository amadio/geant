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
  using GeantTaskData = Geant::GeantTaskData;

  GVectorPhysicsProcess();
  GVectorPhysicsProcess(double energyLimit, int numThreads);
  virtual ~GVectorPhysicsProcess();

  virtual void Initialize();
  virtual void PostStepFinalStateSampling( Material_t* /*mat*/,
                                           int ntracks, 
                                           GeantTrack_v &tracks,
                                           int &nout, 
                                           GeantTaskData *td);

  // these are not active !!! 
  //
  virtual void ComputeIntLen(Material_t * /*mat*/,
                             int /*ntracks*/, 
                             GeantTrack_v & /*tracks*/,
                             double * /*lengths*/, 
                             GeantTaskData * /*tid*/)                            {}
  
  virtual void PostStep(     Material_t * /*mat*/,
                             int /*ntracks*/,
                             GeantTrack_v &/*tracks*/, 
                             int & /*nout*/, 
                             GeantTaskData * /*tid*/)                            {}         

  virtual void PostStepTypeOfIntrActSampling(     Material_t * /*mat*/,
                                                  int /*ntracks*/,
                                                  GeantTrack_v & /*tracks*/, 
                                                  GeantTaskData * /*tid*/)       {} 

  virtual void AtRest(       int /*ntracks*/,
                             GeantTrack_v &/*tracks*/, 
                             int &/*nout*/, 
                             int /*tid*/)                             {}

  virtual void Eloss(        Material_t */*mat*/,
                             int /*ntracks*/,
                             GeantTrack_v &/*tracks*/,
                             int &/*nout*/,
                             int /*tid*/)                             {}

  virtual void ApplyMsc(     Material_t */*mat*/,
                             int /*ntracks*/,
                             GeantTrack_v &/*tracks*/,
                             int /*tid*/)                             {}

private:
   GVectorPhysicsProcess (const GVectorPhysicsProcess  &);//no imp.	
   GVectorPhysicsProcess & operator=(const GVectorPhysicsProcess  &);//no imp.

   // we need this while vecprot_v2/inc/PhysicsProcess is derived from TNamed
   ClassDef(GVectorPhysicsProcess,1)
};

#endif

