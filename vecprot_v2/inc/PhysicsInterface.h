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

/**
 * @brief Class describing physics interface
 */
class PhysicsInterface {

public:
  using GeantTrack_v  = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;

public:
  /**
   * @brief PhysicsInterface constructor
   */
  PhysicsInterface() {}

  /** @brief PhysicsInterface destructor */
  virtual ~PhysicsInterface();

  /** @brief Function of initialization */
  virtual void Initialize() {}

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
