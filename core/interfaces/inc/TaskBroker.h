//===--- TaskBroker.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file TaskBroker.h
 * @brief Implementation of task broker in Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TASKBROKER
#define GEANT_TASKBROKER

#include "Geant/Config.h"

/**
 * @defgroup TGEANT_TASKBROKER GeantV TaskBroker
 *
 * @{
 */

namespace vecgeom {
#ifndef VECCORE_CUDA
  inline namespace cxx {
#else
  namespace cxx {
#endif
    class VPlacedVolume;
  }
}

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;
class GeantBasket;
class GeantConfig;
class GeantPropagator;

/**
 * @brief Class TaskBroker
 */
class TaskBroker {
protected:

  /** @struct TaskData */
  struct TaskData {};

public:

  /** @brief TaskBroker destructor */
  virtual ~TaskBroker(){};

  /** @brief Virtual function that check validity */
  virtual bool IsValid() = 0;
  typedef TaskData *Stream;

  /** @brief Return true if the Task Broker can only process a subset of the particular */
  virtual bool IsSelective() const = 0;

  /** @brief Create the baskets for each stream */
  virtual void CreateBaskets(GeantPropagator *config) = 0;

  /** @brief Virtual function that get next stream */
  virtual Stream GetNextStream() = 0;

  /** @brief Virtual function that adds a track from a basket to be handled by the broker
   *
   * @param itr Track index
   * @param basket Reference to basket to copy from (output array)
   * @return Task broker took ownership of the track or not
   */
  virtual bool addTrack(int itr, GeantBasket &basket) = 0;

  /**
   * @brief Virtual function that provides run task
   *
   * @param threadid Thread ID
   * @param basket GeantBasket object
   */
  virtual void runTask(geant::GeantTaskData &td, GeantBasket &basket) = 0;

  /**
   * @brief Virtual function that launch tasks
   *
   * @param wait Wait parameter (by default false)
   */
  virtual Stream launchTask(bool wait = false) = 0;

  /** @brief Virtual function that checks wait tasks */
  virtual void waitForTasks() = 0;

  /** @brief Virtual function that return total work done */
  virtual long GetTotalWork() = 0;

  /** @brief Virtual function that return number of stream */
  virtual unsigned int GetNstream() = 0;

  /** @brief If the coprocessor has outstanding work, return it */
  virtual GeantBasket *GetBasketForTransport(geant::GeantTaskData &td) = 0;

  /** @brief Tell the tasks which priotizer to use */
  virtual int SetPrioritizer() = 0;

  /** @brief Prepare the geometry for the device and upload it to the device's memory */
  virtual bool UploadGeometry(vecgeom::VPlacedVolume const *const volume = nullptr) = 0;

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif // GEANT_TASKBROKER
/** @} */
