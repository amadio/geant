#ifndef GEANT_PHYSICSPROCESS
#define GEANT_PHYSICSPROCESS
// Toy physics processes for our propagator prototype. Currently including:
// - single scattering as a discrete process
// - energy loss as continuous process
// - generic interaction as discrete process, producing secondaries

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TGeoVolume;
class GeantTrack;

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
  
  bool       IsType(EProcessType type) {return TObject::TestBit(type);}
  virtual void ComputeIntLen(TGeoVolume *vol,
                             int ntracks, 
                             int *trackin, 
                             double *lengths)                             = 0;
  virtual void PostStep(     TGeoVolume *vol,
                             int ntracks,
                             int *trackin, 
                             int &nout, 
                             int* trackout)                             = 0;
  static void StepManager(int iproc, int npart, int */*particles*/, int nout, int */*partnext*/);
  ClassDef(PhysicsProcess,1)    // Physics process base class
};

//______________________________________________________________________________
class ScatteringProcess : public PhysicsProcess
{
public:
  ScatteringProcess() : PhysicsProcess() {TObject::SetBit(kDiscrete);}
  ScatteringProcess(const char *name) : PhysicsProcess(name) {TObject::SetBit(kDiscrete);}
  virtual ~ScatteringProcess() {}
  
  virtual void ComputeIntLen(TGeoVolume *vol,
                             int ntracks, 
                             int *trackin, 
                             double *lengths);
  virtual void PostStep(     TGeoVolume *vol,
                             int ntracks, 
                             int *trackin, 
                             int &nout, 
                             int* trackout);
  ClassDef(ScatteringProcess,1)    // Single scattering process
};

//______________________________________________________________________________
class ElossProcess : public PhysicsProcess
{
public:
  ElossProcess() : PhysicsProcess() {TObject::SetBit(kContinuous);}
  ElossProcess(const char *name) : PhysicsProcess(name) {TObject::SetBit(kContinuous);}
  virtual ~ElossProcess() {}
  
  static double     Bbf1(double *x, double *par);
  static double     BetheBloch(GeantTrack* track, double tz, double ta, double rho);
  void                PlotBB(double z, double a, double rho, double bgmin=1e-2, double bgmax=1e6);
  virtual void ComputeIntLen(TGeoVolume *vol,
                             int ntracks, 
                             int *trackin, 
                             double *lengths);
  virtual void PostStep(     TGeoVolume *vol,
                             int ntracks, 
                             int *trackin, 
                             int &nout, 
                             int* trackout);
  ClassDef(ElossProcess,1)    // Energy loss process
};

//______________________________________________________________________________
class InteractionProcess : public PhysicsProcess
{
public:
  InteractionProcess() : PhysicsProcess() {TObject::SetBit(kDiscrete);}
  InteractionProcess(const char *name) : PhysicsProcess(name) {TObject::SetBit(kDiscrete);}
  virtual ~InteractionProcess() {}
  
  virtual void ComputeIntLen(TGeoVolume *vol,
                             int ntracks, 
                             int *trackin, 
                             double *lengths);
  virtual void PostStep(     TGeoVolume *vol,
                             int ntracks, 
                             int *trackin, 
                             int &nout, 
                             int* trackout);
  ClassDef(InteractionProcess,1)    // Single scattering process
};

#endif
