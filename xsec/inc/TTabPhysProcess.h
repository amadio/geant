#ifndef GEANT_TTABPHYSPROCESS
#define GEANT_TTABPHYSPROCESS
//______________________________________________________________________________
// Generic process for tabulated physics
//______________________________________________________________________________

#include "Geant/Config.h"

#ifndef GEANT_PHYSICSPROCESS
#include "PhysicsProcess.h"
#endif

class TGeoMaterial;
class GeantTrack;
class GeantTrack_v;
class TTabPhysMgr;
class GeantTaskData;

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
  virtual void ComputeIntLen(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Double_t *lengths, GeantTaskData *td);

  // # dummy method: PostStep has been splitted up into two parts (see below)     
  virtual void PostStep( TGeoMaterial* /*mat*/ , 
                         Int_t /*ntracks*/, 
                         GeantTrack_v& /*tracks */,
                         Int_t& /*nout*/,
                         GeantTaskData */*td*/) {}

  // # smapling: -target atom and type of the interaction for each primary tracks 
  //             -all inf. regarding sampling output is stored in the tracks
  virtual void PostStepTypeOfIntrActSampling( TGeoMaterial *mat,
                                              Int_t ntracks, 
                                              GeantTrack_v &tracks, 
                                              GeantTaskData *td);

  // # sampling final states for each primary tracks based on target atom and
  //    interaction type sampled by PostStepTypeOfIntrActSampling;
  // # upadting primary track properties and inserting secondary tracks;
  // # number of inserted secondary tracks will be stored in nout at termination 
  virtual void PostStepFinalStateSampling( TGeoMaterial *mat,
                                           Int_t ntracks, 
                                           GeantTrack_v &tracks,
                                           Int_t &nout, 
                                           GeantTaskData *td);


  virtual void AtRest(Int_t ntracks, GeantTrack_v &tracks, Int_t &nout, GeantTaskData *td);
  GEANT_CUDA_DEVICE_CODE
  virtual void Eloss(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Int_t &nout, GeantTaskData *td);
  GEANT_CUDA_DEVICE_CODE
  virtual void ApplyMsc(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, GeantTaskData *td);

private:
   TTabPhysProcess(const TTabPhysProcess &);//no imp.	
   TTabPhysProcess& operator=(const TTabPhysProcess &);//no imp.

  ClassDef(TTabPhysProcess,1)    // Generic tabulated physics process
};   

#endif
