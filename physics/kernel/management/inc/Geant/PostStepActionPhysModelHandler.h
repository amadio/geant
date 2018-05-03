//===--- PostStepActionPhysModelHandler.h - Geant-V -------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file PostStepActionPhysModelHandler.h
 * @brief Implementation of post step action handler for vectorized EM transport
 * @author Vitalii Drohan
 */
//===----------------------------------------------------------------------===//
#ifndef POSTSTEPACTIONPHYSPROCESSHANDLER_H
#define POSTSTEPACTIONPHYSPROCESSHANDLER_H

// from geantV
#include "Geant/Handler.h"
// from geantV
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {
class Propagator;
class Track;
class TaskData;
class Basket;
}
}

namespace geantphysics {
class EMModel;
/**
 * @brief   Handler for post-step actions specific for physics process (discrete interactions).
 * @class   PostStepActionPhysProcessHandler
 *
 * This handler is used only for specific EM physics model in order to use vectorization over basket of particles which
 * are undergoing same model sampling.
 * (based on PostStepActionHandler by M. Novak)
 */

class PostStepActionPhysModelHandler : public geant::Handler {
public:
  /** @brief Default constructor */
  PostStepActionPhysModelHandler() : geant::Handler() {}

  /**
   * @brief Default constructor
   * @param propagator Propagator working with this handler
   */
  PostStepActionPhysModelHandler(int threshold, geant::Propagator *propagator, int modelIdx);

  /** @brief dtr */
  virtual ~PostStepActionPhysModelHandler();

  /** @brief Scalar DoIt interface */
  virtual void DoIt(geant::Track *track, geant::Basket &output, geant::TaskData *td);

  /** @brief Vector DoIt interface. Base class implements it as a loop. */
  virtual void DoIt(geant::Basket &input, geant::Basket &output, geant::TaskData *td);

private:
  PostStepActionPhysModelHandler(const PostStepActionPhysModelHandler &) = delete;
  PostStepActionPhysModelHandler &operator=(const PostStepActionPhysModelHandler &) = delete;

  void DoItVector(geant::Track **gtracks, int N, geant::Basket &output, geant::TaskData *td);
  void DoItScalar(geant::Track **gtracks, int N, geant::Basket &output, geant::TaskData *td);
  EMModel *fModel;
};

} // namespace geantphysics

#endif // POSTSTEPACTIONHANDLER_H
