
#ifndef CMSFULLAPP_H
#define CMSFULLAPP_H

#ifndef GEANT_VAPPLICATION
#include "Geant/UserApplication.h"
#endif

#include "Geant/Typedefs.h"
#include "Geant/Fwd.h"
#include "Geant/TaskData.h"

namespace GEANT_IMPL_NAMESPACE {
namespace geant {
class RunManager;
class TaskDataHandle;
class Event;
class Track;
}
}

#include "CMSData.h"

#include <mutex>
#include <vector>

namespace cmsapp {

class CMSParticleGun;

class CMSFullApp : public geant::UserApplication {
public:
  CMSFullApp(geant::RunManager *runmgr, CMSParticleGun *gun);
  virtual ~CMSFullApp();

  /** @brief Interface method to allow registration of user defined thread local data. */
  virtual void AttachUserData(geant::TaskData *td);

  /** @brief Interface method to initialize the application. */
  virtual bool Initialize();

  /** @brief Interace method that is called at the end of each simulation step. */
  virtual void SteppingActions(geant::Track &track, geant::TaskData *td);

  /** @brief Interace method that is called when the transportation of an event (including all primary and their
    *        secondary particles) is completed .*/
  virtual void FinishEvent(geant::Event *event);

  /** @brief Interface method that is called at the end of the simulation (when the transportation of all events are
    *        are completed). */
  virtual void FinishRun();

  void SetPerformanceMode(bool val) { fIsPerformance = val; }

private:
  /** @brief Copy constructor TestEm5 (deleted) */
  CMSFullApp(const CMSFullApp &) = delete;
  /** @brief Operator= for TestEm5 (deleted) */
  CMSFullApp &operator=(const CMSFullApp &) = delete;

private:
  bool fIsPerformance;
  bool fInitialized;
  int fNumPrimaryPerEvent;
  int fNumBufferedEvents;
  // user defined thread local data structure handlers to obtain the thread local data structures (defined and
  // registered by the user) during the simulation (in the SteppingActions(i.e. at the end of each simulation step),
  // Digitization(i.e. at the end of an event) and FinishRun (i.e. at the end of the simulation):
  //
  // 1. merged from all working threads when transportation of an event (with all corresponding primary and secondary
  //    particles) are completed
  geant::TaskDataHandle<CMSThreadDataEvents> *fDataHandlerEvents;
  // a unique, run-global user defined data structure to store cumulated quantities per primary particle type
  // during the simulation
  CMSData *fData;

  // mutex to prevent multiple threads writing into the unique, run-global, unique CMSData object (in the FinishEvent
  // after the merge of the user defined per-event data distributed among the threads)
  std::mutex fMutex;
  //
  CMSParticleGun *fGun;
};

} // namespace cmsapp

#endif // CMSFULLAPP_H
