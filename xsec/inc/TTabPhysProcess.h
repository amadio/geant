#ifndef GEANT_TTABPHYSPROCESS
#define GEANT_TTABPHYSPROCESS
//______________________________________________________________________________
// Generic process for tabulated physics
//______________________________________________________________________________

#ifndef GEANT_PHYSICSPROCESS
#include "PhysicsProcess.h"
#endif

class TGeoMaterial;
class GeantTrack;
class GeantTrack_v;
class TTabPhysMgr;

//______________________________________________________________________________
class TTabPhysProcess : public PhysicsProcess
{
private:
  TTabPhysMgr           *fMgr;            //! Tabulated physics manager
  TString                fXsecFileName;   // Name of Xsec file
  TString                fFinalSFileName; // Name of final states file
public:
  TTabPhysProcess();
  TTabPhysProcess(const char *name, const char *fxsec, const char *ffstate);
  virtual ~TTabPhysProcess() {}

  virtual void Initialize();
  virtual void ComputeIntLen(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Double_t *lengths, Int_t tid);
  virtual void PostStep(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Int_t &nout, Int_t tid);
  virtual void AtRest(Int_t ntracks, GeantTrack_v &tracks, Int_t &nout, Int_t tid);
  virtual void Eloss(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Int_t &nout, Int_t tid);
  virtual void ApplyMsc(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid);

private:
   TTabPhysProcess(const TTabPhysProcess &);//no imp.	
   TTabPhysProcess& operator=(const TTabPhysProcess &);//no imp.

  ClassDef(TTabPhysProcess,1)    // Generic tabulated physics process
};   

#endif
