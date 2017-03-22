//===--- PhysicsProcess.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file PhysicsProcessOld.h
 * @brief Definition of physical processes used in Geant-V prototype
 * @details Toy physics processes for our propagator prototype.
 * Currently including:
 * 1. Single scattering as a discrete process;
 * 2. Energy loss as continuous process;
 * 3. Generic interaction as discrete process, producing secondaries.
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_PHYSICSPROCESS
#define GEANT_PHYSICSPROCESS
#include "Geant/Config.h"
#include "Geant/Typedefs.h"
#include "GeantTrack.h"

#include "base/Global.h"
//#include "GeantFwd.h"

#include <string>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Class describing physics processes
 */
class GeantTrack_v;
class GeantTaskData;

class PhysicsProcessOld {

public:

  /**
   * @enum EProcessType
   * @brief Enum ProcessType definition
   */
  enum EProcessType : char { kDiscrete = 1, kContinuous = 2 };

protected:
  std::string fName;

private:
  EProcessType fType;

public:
  /**
   * @brief PhysicsProcess constructor
   */
  PhysicsProcessOld() {}
  /**
   * @brief PhysicsProcess parametrized constructor
   *
   * @param name Name of physics process
   */
  PhysicsProcessOld(const char *name) : fName(name) {}

  /** @brief PhysicsProcess destructor */
  virtual ~PhysicsProcessOld();


  /** @brief Record type information */
  void SetType(EProcessType type) { fType = type; }

  /**
   * @brief Function that check type
   *
   * @param type EProcessType type
   * @return Boolean value
   */
  bool IsType(EProcessType type) { return fType == type; }

  /** @brief Function of initialization */
  virtual void Initialize() {}

  /**
   * @brief Function that computes interaction length
   *
   * @param mat Material_t material
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks_v
   * @param lengths Partial process lengths
   * @param td Thread data
   */
  virtual void ComputeIntLen(Material_t *mat, int ntracks, GeantTrack_v &tracks,
                             GeantTaskData *td) = 0;
  /**
   * @brief Post step type of intraction sampling function
   * @details Sampling:
   * 1. Target atom and type of the interaction for each primary tracks
   * 2. All inf. regarding sampling output is stored in the tracks
   *
   * @param mat Material_t material
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks_v
   * @param td  Thread data
   */
  VECCORE_ATT_DEVICE
  virtual void PostStepTypeOfIntrActSampling(Material_t *mat, int ntracks, GeantTrack_v &tracks,
                                             GeantTaskData *td) = 0;

  /**
   * @brief Post step final state sampling function
   * @details Sampling final states for each primary tracks based on target atom and
   * interaction type sampled by PostStepTypeOfIntrActSampling;
   * updating primary track properties and inserting secondary tracks;
   * number of inserted secondary tracks will be stored in nout at termination;
   *
   * @param mat Material_t material
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks_v
   * @param nout Number of tracks in the output
   * @param td Thread data
   */
  VECCORE_ATT_DEVICE
  virtual void PostStepFinalStateSampling(Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout,
                                          GeantTaskData *td) = 0;
 /**
   * @todo  Need to be implemented
   */
  virtual void AtRest(int /*ntracks*/, GeantTrack_v & /*tracks*/, int & /*nout*/, GeantTaskData * /*td*/) {}

  /**
   * @todo Need to be implemented
   */
  VECCORE_ATT_DEVICE
  virtual void Eloss(Material_t * /*mat*/, int /*ntracks*/, GeantTrack_v & /*tracks*/, int & /*nout*/,
                     GeantTaskData * /*td*/) {}

  /**
   * @todo Need to be implemented
   */
  VECCORE_ATT_DEVICE
  virtual void ApplyMsc(Material_t * /*mat*/, int /*ntracks*/, GeantTrack_v & /*tracks*/, GeantTaskData * /*td*/) {}

//=== N E W   I N T E R F A C E S ===//
  /**
   * @brief Function that computes interaction length
   *
   * @param mat Material_t material
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks_v
   * @param lengths Partial process lengths
   * @param td Thread data
   */
   VECCORE_ATT_HOST_DEVICE
   virtual void ComputeIntLen(GeantTrack *track, GeantTaskData *td) = 0;

   VECCORE_ATT_HOST_DEVICE
   virtual void ComputeIntLen(TrackVec_t &tracks, GeantTaskData *td) = 0;
  /**
   * @brief Post step type of intraction sampling function
   * @details Sampling:
   * 1. Target atom and type of the interaction for each primary tracks
   * 2. All inf. regarding sampling output is stored in the tracks
   *
   * @param mat Material_t material
   * @param tracks Vector of tracks_v
   * @param td  Thread data
   */
  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepTypeOfIntrActSampling(GeantTrack *track, GeantTaskData *td) = 0;
  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepTypeOfIntrActSampling(TrackVec_t &tracks, GeantTaskData *td) = 0;

  /**
   * @brief Post step final state sampling function
   * @details Sampling final states for each primary tracks based on target atom and
   * interaction type sampled by PostStepTypeOfIntrActSampling;
   * updating primary track properties and inserting secondary tracks;
   * number of inserted secondary tracks will be stored in nout at termination;
   *
   * @param mat Material_t material
   * @param ntracks Number of tracks
   * @param tracks Vector of tracks_v
   * @param nout Number of tracks in the output
   * @param td Thread data
   */
  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepFinalStateSampling(GeantTrack *track, int &nout, TrackVec_t &output, GeantTaskData *td) = 0;
  VECCORE_ATT_HOST_DEVICE
  virtual void PostStepFinalStateSampling(TrackVec_t &tracks, int &nout, TrackVec_t &output, GeantTaskData *td) = 0;

 /**
   * @todo  Need to be implemented
   */
  virtual void AtRest(TrackVec_t & /*tracks*/, int & /*nout*/, GeantTaskData * /*td*/) {}

  VECCORE_ATT_HOST_DEVICE
  virtual void Eloss(GeantTrack */*track*/, int &/*nout*/, TrackVec_t &/*output*/,
                     GeantTaskData * /*td*/) {}

  VECCORE_ATT_HOST_DEVICE
  virtual void Eloss(TrackVec_t &/*tracks*/, int &/*nout*/, TrackVec_t &/*output*/,
                     GeantTaskData * /*td*/) {}

  /**
   * @todo Need to be implemented
   */
  VECCORE_ATT_HOST_DEVICE
  virtual void ApplyMsc(Material_t * /*mat*/, TrackVec_t & /*tracks*/, GeantTaskData * /*td*/) {}

//===================================//

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
