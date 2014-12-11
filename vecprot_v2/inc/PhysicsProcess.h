#ifndef GEANT_PHYSICSPROCESS
#define GEANT_PHYSICSPROCESS
// Toy physics processes for our propagator prototype. Currently including:
// - single scattering as a discrete process
// - energy loss as continuous process
// - generic interaction as discrete process, producing secondaries

#include "Geant/Config.h"

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TGeoMaterial;
class GeantTrack;
class GeantTrack_v;
class TGenPhaseSpace;

#include "TMutex.h"

//______________________________________________________________________________
class PhysicsProcess : public TNamed
{
public:
enum EProcessType {
  kDiscrete   = BIT(14),
  kContinuous = BIT(15)
};

public:
  PhysicsProcess() : TNamed() {}
  PhysicsProcess(const char *name) : TNamed(name,"") {}
  virtual ~PhysicsProcess() {}
  
  Bool_t       IsType(EProcessType type) {return TObject::TestBit(type);}
  virtual void Initialize() {}
  virtual void ComputeIntLen(TGeoMaterial *mat,
                             Int_t ntracks, 
                             GeantTrack_v &tracks,
                             Double_t *lengths, 
                             Int_t tid)                             = 0;
  
  virtual void PostStep(     TGeoMaterial *mat,
                             Int_t ntracks,
                             GeantTrack_v &tracks, 
                             Int_t &nout, 
                             Int_t tid)                             = 0;

  // # smapling: -target atom and type of the interaction for each primary tracks 
  //             -all inf. regarding sampling output is stored in the tracks
  virtual void PostStepTypeOfIntrActSampling(     TGeoMaterial *mat,
                                                  Int_t ntracks,
                                                  GeantTrack_v &tracks, 
                                                  Int_t tid)        = 0;

  // # sampling final states for each primary tracks based on target atom and
  //    interaction type sampled by PostStepTypeOfIntrActSampling;
  // # upadting primary track properties and inserting secondary tracks;
  // # number of inserted secondary tracks will be stored in nout at termination 
  virtual void PostStepFinalStateSampling(        TGeoMaterial *mat, 
                                                  Int_t ntracks, 
                                                  GeantTrack_v &tracks,
                                                  Int_t &nout, 
                                                  Int_t tid)        = 0;

  virtual void AtRest(       Int_t /*ntracks*/,
                             GeantTrack_v &/*tracks*/, 
                             Int_t &/*nout*/, 
                             Int_t /*tid*/)                             {}
  GEANT_CUDA_DEVICE_CODE
  virtual void Eloss(        TGeoMaterial */*mat*/,
                             Int_t /*ntracks*/,
                             GeantTrack_v &/*tracks*/,
                             Int_t &/*nout*/,
                             Int_t /*tid*/)                                 {}
  GEANT_CUDA_DEVICE_CODE
  virtual void ApplyMsc(     TGeoMaterial */*mat*/,
                             Int_t /*ntracks*/,
                             GeantTrack_v &/*tracks*/,
                             Int_t /*tid*/)                             {}

  ClassDef(PhysicsProcess,1)    // Physics process base class
};

#endif
