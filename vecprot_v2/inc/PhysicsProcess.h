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

//______________________________________________________________________________
class ScatteringProcess : public PhysicsProcess
{
public:
  ScatteringProcess() : PhysicsProcess() {TObject::SetBit(kDiscrete);}
  ScatteringProcess(const char *name) : PhysicsProcess(name) {TObject::SetBit(kDiscrete);}
  virtual ~ScatteringProcess() {}
  
  virtual void ComputeIntLen(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Double_t *lengths, Int_t tid);
  virtual void PostStep(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Int_t &nout, Int_t tid);
  // dummy methods just keep these ancient Scattering class;
  // however, we shouldn't need them anymore !!! so we should remove this class 
  virtual void PostStepTypeOfIntrActSampling( TGeoMaterial* /* mat */, 
                                              Int_t /* ntracks*/, 
                                              GeantTrack_v &/* tracks*/,
                                              Int_t /*tid*/) {}
  virtual void PostStepFinalStateSampling(    TGeoMaterial* /* mat */, 
                                              Int_t /* ntracks*/, 
                                              GeantTrack_v &/* tracks*/,
                                              Int_t & /*nout*/,
                                              Int_t /*tid*/) {}

  ClassDef(ScatteringProcess,1)    // Single scattering process
};

//______________________________________________________________________________
class ElossProcess : public PhysicsProcess
{
public:
  ElossProcess() : PhysicsProcess() {TObject::SetBit(kContinuous);}
  ElossProcess(const char *name) : PhysicsProcess(name) {TObject::SetBit(kContinuous);}
  virtual ~ElossProcess() {}
  
//  static Double_t     Bbf1(Double_t *x, Double_t *par);
  static Double_t     BetheBloch(const GeantTrack_v &tracks, Int_t itrack, Double_t tz, Double_t ta, Double_t rho);
//  void                PlotBB(Double_t z, Double_t a, Double_t rho, Double_t bgmin=1e-2, Double_t bgmax=1e6);
  virtual void ComputeIntLen(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Double_t *lengths, Int_t tid);
  virtual void PostStep(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Int_t &nout, Int_t tid);
  // dummy methods just keep these ancient ElossProcess class;
  // however, we shouldn't need them anymore !!! so we should remove this class 
  virtual void PostStepTypeOfIntrActSampling( TGeoMaterial* /* mat */, 
                                              Int_t /* ntracks*/, 
                                              GeantTrack_v &/* tracks*/,
                                              Int_t /*tid*/) {}
  virtual void PostStepFinalStateSampling(    TGeoMaterial* /* mat */, 
                                              Int_t /* ntracks*/, 
                                              GeantTrack_v &/* tracks*/,
                                              Int_t & /*nout*/,
                                              Int_t /*tid*/) {}


  ClassDef(ElossProcess,1)    // Energy loss process
};

//______________________________________________________________________________
class InteractionProcess : public PhysicsProcess
{
private:
  Int_t                fNthreads; // Number of threads
  TGenPhaseSpace      *fGen; //[fNthreads] Phase space generator
  TMutex               fMutex; //! mutex
  
  InteractionProcess(const InteractionProcess&); // not implemented
  InteractionProcess &operator=(const InteractionProcess&); // not implemented

public:
  InteractionProcess() : PhysicsProcess(), fNthreads(0), fGen(0), fMutex() {TObject::SetBit(kDiscrete);}
  InteractionProcess(const char *name);
  virtual ~InteractionProcess();
  
  virtual void ComputeIntLen(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Double_t *lengths, Int_t tid);
  virtual void PostStep(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Int_t &nout, Int_t tid);
  // dummy methods just keep these ancient InteractionProcess class;
  // however, we shouldn't need them anymore !!! so we should remove this class 
  virtual void PostStepTypeOfIntrActSampling( TGeoMaterial* /* mat */, 
                                              Int_t /* ntracks*/, 
                                              GeantTrack_v &/* tracks*/,
                                              Int_t /*tid*/) {}
  virtual void PostStepFinalStateSampling(    TGeoMaterial* /* mat */, 
                                              Int_t /* ntracks*/, 
                                              GeantTrack_v &/* tracks*/,
                                              Int_t & /*nout*/,
                                              Int_t /*tid*/) {}


  ClassDef(InteractionProcess,1)    // Single scattering process
};

#endif
