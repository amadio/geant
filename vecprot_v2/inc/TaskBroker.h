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

/**
 * @defgroup TGEANT_TASKBROKER GeantV TaskBroker
 *
 * @{
 */
class GeantBasket;
namespace vecgeom {
   inline namespace cxx {
      class VPlacedVolume;
   }
}


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
  virtual void runTask(int threadid, GeantBasket &basket) = 0;

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
  virtual int GetNstream() = 0;

   /** @brief Tell the tasks which priotizer to use */
   virtual int SetPrioritizer() = 0;

   /** @brief Prepare the geometry for the device and upload it to the device's memory */
   virtual bool UploadGeometry(vecgeom::VPlacedVolume const *const volume = nullptr) = 0;

};

#endif // GEANT_TASKBROKER
/** @} */
