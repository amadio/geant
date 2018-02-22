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

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
  class Propagator;
  class SimulationStage;
  class TrackDataMgr;
}
}

/**
 * @brief Class describing physics interface
 */
class PhysicsInterface {

public:
  using GeantTaskData = geant::GeantTaskData;
  using TrackDataMgr  = geant::TrackDataMgr;

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
  // These methods are called from the geant::Propagator::CreateSimulationStages
  // methods (when real-physics is used) to obtain the pointers to the physics
  // simulation stages defined in the real-physics library.
  /** @brief Obtain/create physics step limit computation stage.
    *
    * @param[in,out] prop  Pointer to the propagator object that requires the simulation stage.
    * @return     Pointer to a created ComputeIntLen real-physics simulation stage object.
    */
  virtual  geant::SimulationStage* CreateComputeIntLStage(geant::Propagator *prop) = 0;
  virtual  geant::SimulationStage* CreatePrePropagationStage(geant::Propagator *prop) = 0;
  virtual  geant::SimulationStage* CreatePostPropagationStage(geant::Propagator *prop) = 0;

  /** @brief Obtain/create along step action (continuous part) computation stage.
    *
    * @param[in,out] prop  Pointer to the propagator object that requires the simulation stage.
    * @return     Pointer to a created AlongStepAction real-physics simulation stage object.
    */
  virtual  geant::SimulationStage* CreateAlongStepActionStage(geant::Propagator *prop) = 0;
  /** @brief Obtain/create post step action (discrete part) computation stage.
    *
    * @param[in,out] prop  Pointer to the propagator object that requires the simulation stage.
    * @return     Pointer to a created PostStepAction real-physics simulation stage object.
    */
  virtual  geant::SimulationStage* CreatePostStepActionStage(geant::Propagator *prop) = 0;

  virtual  geant::SimulationStage* CreateAtRestActionStage(geant::Propagator *prop) = 0;

};

#endif
