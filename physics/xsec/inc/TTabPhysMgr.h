#ifndef TTabPhysMgr_H
#define TTabPhysMgr_H

#include "Geant/Config.h"
#include "GeantTrack.h"

#define MAXNELEMENTS 20 // max number of elements in one material(TMXsec)

// Singleton that handles tabulated physics data loaded for the
// materials present in the geometry.

class TEXsec;
class TMXsec;
class TEFstate;
class TPDecay;

#include "Geant/Typedefs.h"
#ifdef USE_VECGEOM_NAVIGATOR
#include "base/Global.h"
#else
class TGeoManager;
class TGeoMaterial;
#endif

#include "GeantFwd.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class TTabPhysMgr {
public:
  using GeantTrack = Geant::GeantTrack;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;

private:
  int fNelements;         // Total number of elements in the geometry
  int fNmaterials;        // Total number of materials in the geometry
  TEXsec **fElemXsec;     // Array of x-section pointers per element
  TEFstate **fElemFstate; // Array of final state pointers per element
  TMXsec **fMatXsec;      // Array of x-section pointers per material
  TPDecay *fDecay;        // Decay tables for each particles
#ifndef USE_VECGEOM_NAVIGATOR
  TGeoManager *fGeom; // Pointer to the geometry manager
#endif
  bool *fHasNCaptureAtRest; // do the particle have nCapture at rest?

  static TTabPhysMgr *fgInstance; // Singleton instance

public:
  TTabPhysMgr();
  TTabPhysMgr(const char *xsecfilename, const char *finalsfilename);
  virtual ~TTabPhysMgr();
  static TTabPhysMgr *Instance(const char *xsecfilename = 0, const char *finalsfilename = 0);
  static const char *ClassName() { return "TTabPhysMgr"; }
  // Rotation+boost utility
  void TransformLF(int indref, GeantTrack_v &tracks, int nproducts, int indprod,
                   GeantTrack_v &output); // not. imp. but done
  // API used by particle transport
  VECCORE_ATT_HOST_DEVICE
  void ApplyMsc(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td);
  VECCORE_ATT_HOST_DEVICE
  int Eloss(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td);
  void ProposeStep(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td);
  int SampleDecay(int ntracks, GeantTrack_v &tracksin, GeantTrack_v &tracksout); // not. imp.

  // # sampling target, type of interaction, final states;
  // # updating primary track properties and inserting secondary tracks;
  // int SampleInt(int imat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td);

  // # smapling: target atom and type of the interaction for each primary tracks
  void SampleTypeOfInteractions(int imat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td);

  // # sampling final states for each primary tracks based on target atom and
  //    interaction type sampled in SampleTypeOfInteractionsInt;
  // # upadting primary track properties and inserting secondary tracks;
  // # return: number of inserted secondary tracks
  int SampleFinalStates(int imat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  void GetRestFinStates(int partindex, TMXsec *mxs, double energyLimit, GeantTrack_v &tracks, int iintrack,
                        int &nTotSecPart, GeantTaskData *td);
  void SampleDecayInFlight(int partindex, TMXsec *mxs, double energyLimit, GeantTrack_v &tracks, int iintrack,
                           int &nTotSecPart, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  bool HasRestProcess(int gvindex);

  void RotateNewTrack(double oldXdir, double oldYdir, double oldZdir, GeantTrack &track);
  VECCORE_ATT_HOST_DEVICE
  void RotateNewTrack(double oldXdir, double oldYdir, double oldZdir, GeantTrack_v &tracks, int itrack);
  void RotateTrack(GeantTrack &track, double theta, double phi);
  VECCORE_ATT_HOST_DEVICE
  void RotateTrack(GeantTrack_v &tracks, int itrack, double theta, double phi);

  // get current version number
  int VersionMajor() const { return fgVersion / 1000 / 1000; }
  int VersionMinor() const { return fgVersion / 1000 - VersionMajor() * 1000; }
  int VersionSub() const { return fgVersion - VersionMajor() * 1000000 - VersionMinor() * 1000; }
  const char *GetVersion() const;

  TPDecay* GetDecayTable() {return fDecay;}

//=== N E W   I N T E R F A C E S ===//
  // Rotation+boost utility
  VECCORE_ATT_HOST_DEVICE
  void TransformLF(int indref, TrackVec_t &tracks, int nproducts, int indprod,
                   TrackVec_t &output); // not. imp. but done

  VECCORE_ATT_HOST_DEVICE
  void ApplyMsc(Material_t *mat, TrackVec_t &tracks, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  int Eloss(GeantTrack *track, TrackVec_t &output, GeantTaskData *td);
  VECCORE_ATT_HOST_DEVICE
  int Eloss(TrackVec_t &tracks, TrackVec_t &output, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  void ProposeStep(GeantTrack *track, GeantTaskData *td);
  VECCORE_ATT_HOST_DEVICE
  void ProposeStep(TrackVec_t &tracks, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  int SampleDecay(TrackVec_t &tracks, TrackVec_t &output); // not. imp.
  // # sampling target, type of interaction, final states;
  // # updating primary track properties and inserting secondary tracks;
  // int SampleInt(int imat, TrackVec_t &tracks, GeantTaskData *td);
 // VECCORE_ATT_HOST_DEVICE
 // int SampleInt(int imat, TrackVec_t &tracks, GeantTaskData *td);

  // # smapling: target atom and type of the interaction for each primary tracks
  VECCORE_ATT_HOST_DEVICE
  void SampleTypeOfInteractions(GeantTrack *track, GeantTaskData *td);
  VECCORE_ATT_HOST_DEVICE
  void SampleTypeOfInteractions(TrackVec_t &tracks, GeantTaskData *td);

  // # sampling final states for each primary tracks based on target atom and
  //    interaction type sampled in SampleTypeOfInteractionsInt;
  // # upadting primary track properties and inserting secondary tracks;
  // # return: number of inserted secondary tracks
  VECCORE_ATT_HOST_DEVICE
  int SampleFinalStates(GeantTrack *track, TrackVec_t &output, GeantTaskData *td);
  VECCORE_ATT_HOST_DEVICE
  int SampleFinalStates(TrackVec_t &tracks, TrackVec_t &output, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  void GetRestFinStates(int partindex, TMXsec *mxs, double energyLimit, GeantTrack *track,
                        int &nTotSecPart, TrackVec_t &output, GeantTaskData *td);
  VECCORE_ATT_HOST_DEVICE
  void SampleDecayInFlight(int partindex, TMXsec *mxs, double energyLimit, GeantTrack *track,
                           int &nTotSecPart, TrackVec_t &output, GeantTaskData *td);

//===================================//

private:
  TTabPhysMgr(const TTabPhysMgr &);            // no imp.
  TTabPhysMgr &operator=(const TTabPhysMgr &); // no imp.

  // current version number
  static const int fgVersion = 1000002;

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif // TTabPhysMgr_H
