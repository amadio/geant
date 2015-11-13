#ifndef TTabPhysMgr_H
#define TTabPhysMgr_H

#include "Geant/Config.h"

#ifdef USE_ROOT
#include "Rtypes.h"
#endif

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
  GEANT_CUDA_DEVICE_CODE
  void ApplyMsc(Material_t *mat, int ntracks, GeantTrack_v &tracks, GeantTaskData *td);
  GEANT_CUDA_DEVICE_CODE
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

  GEANT_CUDA_DEVICE_CODE
  void GetRestFinStates(int partindex, TMXsec *mxs, double energyLimit, GeantTrack_v &tracks, int iintrack,
                        int &nTotSecPart, GeantTaskData *td);
  void SampleDecayInFlight(int partindex, TMXsec *mxs, double energyLimit, GeantTrack_v &tracks, int iintrack,
                           int &nTotSecPart, GeantTaskData *td);

  GEANT_CUDA_DEVICE_CODE
  bool HasRestProcess(int gvindex);

  void RotateNewTrack(double oldXdir, double oldYdir, double oldZdir, GeantTrack &track);
  GEANT_CUDA_DEVICE_CODE
  void RotateNewTrack(double oldXdir, double oldYdir, double oldZdir, GeantTrack_v &tracks, int itrack);
  void RotateTrack(GeantTrack &track, double theta, double phi);
  GEANT_CUDA_DEVICE_CODE
  void RotateTrack(GeantTrack_v &tracks, int itrack, double theta, double phi);

  // get current version number
  int VersionMajor() const { return fgVersion / 1000 / 1000; }
  int VersionMinor() const { return fgVersion / 1000 - VersionMajor() * 1000; }
  int VersionSub() const { return fgVersion - VersionMajor() * 1000000 - VersionMinor() * 1000; }
  const char *GetVersion() const;

private:
  TTabPhysMgr(const TTabPhysMgr &);            // no imp.
  TTabPhysMgr &operator=(const TTabPhysMgr &); // no imp.

  // current version number
  static const int fgVersion = 1000002;

#ifdef USE_ROOT
  ClassDefNV(TTabPhysMgr, 2)
#endif
};

#endif // TTabPhysMgr_H
