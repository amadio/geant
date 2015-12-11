#ifndef GEANT_TTABPHYSPROCESS
#define GEANT_TTABPHYSPROCESS
//______________________________________________________________________________
// Generic process for tabulated physics
//______________________________________________________________________________

#include "Geant/Config.h"

#ifndef GEANT_PHYSICSPROCESS
#include "PhysicsProcess.h"
#endif

#include "base/Global.h"
#include "Geant/Typedefs.h"

#include "GeantFwd.h"
class TTabPhysMgr;

//______________________________________________________________________________
class TTabPhysProcess : public PhysicsProcess {
private:
  TTabPhysMgr *fMgr;       //! Tabulated physics manager
  TString fXsecFileName;   // Name of Xsec file
  TString fFinalSFileName; // Name of final states file
public:
  TTabPhysProcess();
  TTabPhysProcess(const char *name, const char *fxsec, const char *ffstate);
  virtual ~TTabPhysProcess() {}

  virtual void Initialize();
  virtual void ComputeIntLen(Material_t *mat, int ntracks, GeantTrack_v &tracks, double *lengths,
                             GeantTaskData *td);

  // # dummy method: PostStep has been splitted up into two parts (see below)
  virtual void PostStep(Material_t * /*mat*/, int /*ntracks*/, GeantTrack_v & /*tracks */, int & /*nout*/,
                        GeantTaskData * /*td*/) {}

  // # smapling: -target atom and type of the interaction for each primary tracks
  //             -all inf. regarding sampling output is stored in the tracks
  virtual void PostStepTypeOfIntrActSampling(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td);

  // # sampling final states for each primary tracks based on target atom and
  //    interaction type sampled by PostStepTypeOfIntrActSampling;
  // # upadting primary track properties and inserting secondary tracks;
  // # number of inserted secondary tracks will be stored in nout at termination
  virtual void PostStepFinalStateSampling(Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout,
                                          GeantTaskData *td);

  virtual void AtRest(int ntracks, GeantTrack_v &tracks, int &nout, GeantTaskData *td);
  GEANT_CUDA_DEVICE_CODE
  virtual void Eloss(Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout, GeantTaskData *td);
  GEANT_CUDA_DEVICE_CODE
  virtual void ApplyMsc(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td);

private:
  TTabPhysProcess(const TTabPhysProcess &);            // no imp.
  TTabPhysProcess &operator=(const TTabPhysProcess &); // no imp.

  ClassDef(TTabPhysProcess, 1) // Generic tabulated physics process
};

#endif
