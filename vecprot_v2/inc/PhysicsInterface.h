//===--- PhysicsInterface.h - Geant-V -----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file PhysicsInterface.h
 * @brief Physics interface to real physics simulation.
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_PHYSICSINTERFACE
#define GEANT_PHYSICSINTERFACE
#include "Geant/Config.h"
#include "Geant/Typedefs.h"

#include "base/Global.h"
#include "GeantFwd.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {
  class GeantPropagator;
  class SimulationStage;
  class TrackDataMgr;
}
}

/**
 * @brief Class describing physics interface
 */
class PhysicsInterface {

public:
  using GeantTrack_v  = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;
  using TrackDataMgr  = Geant::TrackDataMgr;

public:
  /**
   * @brief PhysicsInterface constructor
   */
  PhysicsInterface() {}

  /** @brief PhysicsInterface destructor */
  virtual ~PhysicsInterface();

  /** @brief Attach task data if needed */
  virtual void AttachUserData(GeantTaskData *) {}

  /** @brief Function of initialization */
  virtual void Initialize() {}

  // Interface methods to obtain physics realted symulation stages when V3 is used.
  // These methods are called from the Geant::GeantPropagator::CreateSimulationStages
  // methods (when real-physics is used) to obtain the pointers to the physics
  // simulation stages defined in the real-physics library.
  /** @brief Obtain/create physics step limit computation stage.
    *
    * @param[in,out] prop  Pointer to the propagator object that requires the simulation stage.
    * @return     Pointer to a created ComputeIntLen real-physics simulation stage object.
    */
  virtual  Geant::SimulationStage* CreateComputeIntLStage(Geant::GeantPropagator *prop) = 0;
  virtual  Geant::SimulationStage* CreatePrePropagationStage(Geant::GeantPropagator *prop) = 0;
  virtual  Geant::SimulationStage* CreatePostPropagationStage(Geant::GeantPropagator *prop) = 0;

  /** @brief Obtain/create along step action (continuous part) computation stage.
    *
    * @param[in,out] prop  Pointer to the propagator object that requires the simulation stage.
    * @return     Pointer to a created AlongStepAction real-physics simulation stage object.
    */
  virtual  Geant::SimulationStage* CreateAlongStepActionStage(Geant::GeantPropagator *prop) = 0;
  /** @brief Obtain/create post step action (discrete part) computation stage.
    *
    * @param[in,out] prop  Pointer to the propagator object that requires the simulation stage.
    * @return     Pointer to a created PostStepAction real-physics simulation stage object.
    */
  virtual  Geant::SimulationStage* CreatePostStepActionStage(Geant::GeantPropagator *prop) = 0;


  /**
   * @brief Method that computes the physics step limit.
   *
   * Both continuous and discrete physics step limits are included. The shorter will be provided.
   *
   * @param[in,out] mat Material_t material
   * @param[in,out] ntracks Number of tracks
   * @param[in,out] tracks Vector of tracks_v
   * @param[in,out] lengths Partial process lengths
   * @param[in,out] td Thread data
   */
  virtual void ComputeIntLen(Material_t * /*mat*/, int /*ntracks*/, GeantTrack_v & /*tracks*/, double * /*lengths*/,
                             GeantTaskData * /*td*/) {}//= 0;


  /**
   * @brief Method that performs the along-the-step, i.e. the continuous part, action of the in-flight interactions if
   *        any.
   *
   * @param[in,out] mat Material_t material
   * @param[in,out] ntracks Number of input tracks
   * @param[in,out] tracks Vector of input tracks_v
   * @param[in,out] nout  Number of secondary tracks created
   * @param[in,out] td Thread data
   */
  virtual void AlongStepAction(Material_t * /*mat*/, int /*ntracks*/, GeantTrack_v & /*tracks*/, int & /*nout*/,
                               GeantTaskData * /*td*/) {}

  /**
   * @brief Method that performs the post step action i.e. discrete interaction.
   *
   * @param[in,out] mat Material_t material
   * @param[in,out] ntracks Number of input tracks
   * @param[in,out] tracks Vector of input tracks_v
   * @param[in,out] nout  Number of secondary tracks created
   * @param[in,out] td Thread data
   */
  virtual void PostStepAction(Material_t * /*mat*/, int /*ntracks*/, GeantTrack_v & /*tracks*/, int & /*nout*/,
                              GeantTaskData * /*td*/) {} //= 0;

  /**
   * @brief Method that performs the at rest action i.e. discrete interaction.
   *
   * @param[in,out] mat Material_t material
   * @param[in,out] ntracks Number of input tracks
   * @param[in,out] tracks Vector of input tracks_v
   * @param[in,out] nout  Number of secondary tracks created
   * @param[in,out] td Thread data
   */
  virtual void AtRestAction(Material_t * /*mat*/, int /*ntracks*/, GeantTrack_v & /*tracks*/, int & /*nout*/,
                            GeantTaskData * /*td*/) {};


};

#endif
