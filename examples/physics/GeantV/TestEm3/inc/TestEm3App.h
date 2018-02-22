
#ifndef TESTEM3APP_H
#define TESTEM3APP_H

#ifndef GEANT_VAPPLICATION
#include "GeantVApplication.h"
#endif

#include "Geant/Typedefs.h"
#include "GeantFwd.h"
#include "GeantTaskData.h"

namespace GEANT_IMPL_NAMESPACE {
  namespace geant {
    class RunManager;
    class GeantTaskDataHandle;
    class Event;
    class Track;
  }
}


#include "TestEm3Data.h"

#include <mutex>
#include <vector>


namespace userapplication {

/**
 * @brief GeantV implementation of the Geant4 TestEm3 application.
 *
 * The application simulates passage of particles (type and energy configurable) through a simple calorimeter (material,
 * thickness, number of absorbers/layers, secondary production cuts are configurable). The main purpose of the simulation
 * is to generate results for energy deposition and number/type of generated secondaries. However, several other quantities
 * (like energy leakage[both primary and secondary] and energy balance, mean number of charged and neutral
 * steps in the target, mean number of secondary partciles per particle type, transmitting/backscattering coefficients,
 * etc.) will also be collected and reported at the simulation.
 *
 * @class   TestEm3App
 * @author  M. Novak
 * @date    November 2017
 */
//
class TestEm3DetectorConstruction;
class TestEm3PrimaryGenerator;
//
/** @brief TestEm3App user application */
class TestEm3App : public geant::GeantVApplication {
public:

  /** @brief Constructor TestEm3App */
  TestEm3App(geant::RunManager *runmgr, TestEm3DetectorConstruction *det, TestEm3PrimaryGenerator *gun);

  /** @brief Destructor TestEm3App*/
  virtual ~TestEm3App();

  /** @brief Interface method to allow registration of user defined thread local data. */
  virtual void AttachUserData(geant::GeantTaskData *td);

  /** @brief Interface method to initialize the application. */
  virtual bool Initialize();

  /** @brief Interace method that is called at the end of each simulation step. */
  virtual void SteppingActions(geant::Track &track, geant::GeantTaskData *td);

  /** @brief Interace method that is called when the transportation of an event (including all primary and their
    *        secondary particles) is completed .*/
  virtual void FinishEvent(geant::Event *event);

  /** @brief Interface method that is called at the end of the simulation (when the transportation of all events are
    *        are completed). */
  virtual void FinishRun();


  void         SetPerformanceMode(bool val) { fIsPerformance = val; }

private:
  /** @brief Copy constructor TestEm3App (deleted) */
  TestEm3App(const TestEm3App &) = delete;
  /** @brief Operator= for TestEm3App (deleted) */
  TestEm3App &operator=(const TestEm3App &) = delete;


private:
  bool             fIsPerformance;
  bool             fInitialized;
  int	             fNumAbsorbers;
  std::vector<int> fAbsorberLogicalVolumeID;
  // some data regarding the number of primaries per event and number of buffered events (i.e. number of event-slots)
  // these data will be obtained from the RunManager::GeantConfig object at initialization
  int              fNumPrimaryPerEvent;
  int              fNumBufferedEvents;
  //
  double           fPrimaryParticleCharge;
  // user defined thread local data structure handlers to obtain the thread local data structures (defined and
  // registered by the user) during the simulation (in the SteppingActions(i.e. at the end of each simulation step),
  // Digitization(i.e. at the end of an event) and FinishRun (i.e. at the end of the simulation):
  //
  // 1. merged from all working threads when transportation of an event (with all corresponding primary and secondary
  //    particles) are completed
  geant::TaskDataHandle<TestEm3ThreadDataEvents>  *fDataHandlerEvents;
  // a unique, run-global user defined data structure to store cumulated quantities per primary particle during the simulation
  TestEm3Data  *fData;
  //
  TestEm3DetectorConstruction      *fDetector;
  TestEm3PrimaryGenerator          *fPrimaryGun;
  // mutex to prevent multiple threads writing into the unique, run-global TestEm5Data object (in the Digitization after
  // the merge of the user defined per-event data distributed among the threads)
  std::mutex fMutex;
};

}      // namespace userapplication

#endif // TESTEM3_H
