#ifndef GEANT_TTABPHYSPROCESS
#define GEANT_TTABPHYSPROCESS
//______________________________________________________________________________
// Generic process for tabulated physics
//______________________________________________________________________________

#include "Geant/Config.h"

#include "PhysicsProcessOld.h"

#include "base/Global.h"
#include "Geant/Typedefs.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class TTabPhysMgr;

//______________________________________________________________________________
class TTabPhysProcess : public PhysicsProcessOld {

private:
  TTabPhysMgr *fMgr;       //! Tabulated physics manager
  std::string fXsecFileName;   // Name of Xsec file
  std::string fFinalSFileName; // Name of final states file
public:
  TTabPhysProcess();
  TTabPhysProcess(const char *name, const char *fxsec, const char *ffstate);
  virtual ~TTabPhysProcess() {}

  virtual void Initialize();
  virtual void ComputeIntLen(Material_t *mat, int ntracks, GeantTrack_v &tracks,
                             GeantTaskData *td);

  // # smapling: -target atom and type of the interaction for each primary tracks
  //             -all inf. regarding sampling output is stored in the tracks
  VECCORE_ATT_DEVICE
  virtual void PostStepTypeOfIntrActSampling(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td);

  // # sampling final states for each primary tracks based on target atom and
  //    interaction type sampled by PostStepTypeOfIntrActSampling;
  // # upadting primary track properties and inserting secondary tracks;
  // # number of inserted secondary tracks will be stored in nout at termination
  VECCORE_ATT_DEVICE
  virtual void PostStepFinalStateSampling(Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout,
                                          GeantTaskData *td);

  virtual void AtRest(int ntracks, GeantTrack_v &tracks, int &nout, GeantTaskData *td);
  VECCORE_ATT_DEVICE
  virtual void Eloss(Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout, GeantTaskData *td);
  VECCORE_ATT_DEVICE
  virtual void ApplyMsc(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td);

//=== N E W   I N T E R F A C E S ===//
  VECCORE_ATT_HOST_DEVICE
  virtual void ComputeIntLen(GeantTrack *track, GeantTaskData *td);
  VECCORE_ATT_HOST_DEVICE
  virtual void ComputeIntLen(TrackVec_t &tracks, GeantTaskData *td);
  // # smapling: -target atom and type of the interaction for each primary tracks
  //             -all inf. regarding sampling output is stored in the tracks

  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepTypeOfIntrActSampling(GeantTrack *track, GeantTaskData *td);
  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepTypeOfIntrActSampling(TrackVec_t &tracks, GeantTaskData *td);

  // # sampling final states for each primary tracks based on target atom and
  //    interaction type sampled by PostStepTypeOfIntrActSampling;
  // # upadting primary track properties and inserting secondary tracks;
  // # number of inserted secondary tracks will be stored in nout at termination
  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepFinalStateSampling(GeantTrack *track, int &nout, TrackVec_t &output, GeantTaskData *td);
  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepFinalStateSampling(TrackVec_t &tracks, int &nout, TrackVec_t &output, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  virtual void AtRest(TrackVec_t &tracks, int &nout, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  virtual void Eloss(GeantTrack *track, int &nout, TrackVec_t &output, GeantTaskData *td);
  VECCORE_ATT_HOST_DEVICE
  virtual void Eloss(TrackVec_t &tracks, int &nout, TrackVec_t &output, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  virtual void ApplyMsc(Material_t *mat, TrackVec_t &tracks, GeantTaskData *td);

//===================================//

private:
  TTabPhysProcess(const TTabPhysProcess &);            // no imp.
  TTabPhysProcess &operator=(const TTabPhysProcess &); // no imp.

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
