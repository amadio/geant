
#ifndef TESTEM3DATA_H
#define TESTEM3DATA_H

#include <vector>


namespace userapplication {

/**
 * @brief   User defined data structures for the TestEm3App GeantV application to store/handle scoring data.
 * @author  M Novak
 * @date    November 2017
 *
 * GeantV takes multiple events(number of event-slots=number of buffered events), multiple primary partciles per event
 * and transport them on the same time by multiple working threads. It means that tracks (primary, secondary) that
 * belong to more than one events are distributed among multiple threads. GeantV associates a thread local data storage
 * (Geant::GeantTaskData) to each working threads and provides the possibility to the user(user application) to register
 * application specific thread local (data) objects in this thread local storage. When the transportation of an event,
 * that currently occupies one of the event-slots (one out of the "number of buffered events" possible places) is
 * completed (i.e. all the primary particles that belongs to the corresponding event and their secondaries are fully
 * transported) an application interface method [Digitize(Geant::GeantEvent*)] is invoked. GeantV provides the
 * possibility to merge the thread local user defined data (filled in the SteppingAction() interface method after each
 * simulation step in a thread local way) related to the finished event and distributed among multiple working threads.
 * This merge will return with a pointer to one of the user defined thread local data that contains the result of the
 * merge. At this point, the user can obtain all information that was defined and filled in the SteppingAction related
 * to the transported event and insert into a global data structure that cummulates information during the simulation.
 *
 * This file contains data structure to describe/handle data:
 * - per-primary particle: many simulated quantities will be collected and their mean values per primary particle will
 *                         be computed during the simulation
 * - per-event           : a GeantV event can contain more than one primary particles that are transported on the same
 *                         time. So the per-event data contains as many per-primary data structures as number of primary
 *                         particles per event.
 * - per-thread          : GeantV takes "number of buffered events" events and transport them on the same time. So the
 *                         per-thread data contains as many per-event data structures as number of events are taken and
 *                         transported on the same time.
 * Some of these per-thread data might need to be merged (from the threads) after an event is transported (after all
 * the primaries building up the event are transported like energy deposit per primary partcile if sima is computed).
 * These data will be described in the ThreadDataEvents object in this application and follows the structure described
 * above. Other data do not need to be merged from the threads after an event is completed but enough to be merged at
 * the very end of the simulation (like the angular distribution histogram). These data will be described by the
 * ThreadDataRun object.
 * An additional, global data structure, TestEm3Data is defined in this application to accumulate data that belongs to
 * an event/primary particle and available after the event is transported and the corresponding thread local data are
 * merged. Only one object is created from this global data structure and writing into this global data structure is
 * protected (more than one threads can finish an (different)event and try to write into this global data structure).
 */


// Data structure per-primary particle for TestEm3App
class TestEm3DataPerPrimary {
public:
  TestEm3DataPerPrimary(int numabs);
 ~TestEm3DataPerPrimary();


  void   AddChargedStep()        { fNumChargedSteps += 1.;  }
  double GetChargedSteps() const { return fNumChargedSteps; }

  void   AddNeutralStep()        { fNumNeutralSteps += 1.;  }
  double GetNeutralSteps() const { return fNumNeutralSteps; }

  void   AddChargedTrackL(double val, int absorber) { fChargedTrackL[absorber] += val; }
  double GetChargedTrackL(int absorber) const       { return fChargedTrackL[absorber]; }

  void   AddNeutralTrackL(double val, int absorber) { fNeutralTrackL[absorber] += val; }
  double GetNeutralTrackL(int absorber) const       { return fNeutralTrackL[absorber]; }

  void   AddGamma()           { fNumGammas += 1.;     }
  double GetGammas()    const { return fNumGammas;    }

  void   AddElectron()        { fNumElectrons += 1.;  }
  double GetElectrons() const { return fNumElectrons; }

  void   AddPositron()        { fNumPositrons += 1.;  }
  double GetPositrons() const { return fNumPositrons; }

  void   AddEdepInAbsorber(double val, int absorber)     { fEdepInAbsorber[absorber] += val; }
  double GetEdepInAbsorber(int absorber)           const { return fEdepInAbsorber[absorber]; }

  void   Clear();
  TestEm3DataPerPrimary& operator+=(const TestEm3DataPerPrimary& other);

private:
  int     fNumAbsorbers;        // number of absorbers
  double  fNumChargedSteps;     // mean number of charged steps per primary in calorimeter
  double  fNumNeutralSteps;     // mean number of neutral steps per primary in calorimeter
  double  fNumGammas;           // mean number of secondary gamma particles per primary
  double  fNumElectrons;        // mean number of secondary electron particles per primary
  double  fNumPositrons;        // mean number of secondary positron particles per primary
  std::vector<double>  fEdepInAbsorber;      // mean energy deposit per primary in the absorber
  std::vector<double>  fChargedTrackL;       // mean number of charged track length per primary in absorber
  std::vector<double>  fNeutralTrackL;       // mean number of neutral track length per primary in absorber
};

// Global data structure to accumulate per-primary data during the simulation. The only one object from this class
// will be updated each time an event(with the corresponding primaries) is transported.
class TestEm3Data {
public:
  TestEm3Data(int numabs);
 ~TestEm3Data();

 void   AddChargedSteps(double val) { fNumChargedSteps += val; fNumChargedSteps2 += val*val; }
 double GetChargedSteps()  const    { return fNumChargedSteps;  }
 double GetChargedSteps2() const    { return fNumChargedSteps2; }

 void   AddNeutralSteps(double val) { fNumNeutralSteps += val; fNumNeutralSteps2 += val*val; }
 double GetNeutralSteps() const     { return fNumNeutralSteps;  }
 double GetNeutralSteps2() const    { return fNumNeutralSteps2; }

 void   AddChargedTrackL(double val, int absorber) { fChargedTrackL[absorber] += val; fChargedTrackL2[absorber] += val*val; }
 double GetChargedTrackL(int absorber)  const    { return fChargedTrackL[absorber];  }
 double GetChargedTrackL2(int absorber) const    { return fChargedTrackL2[absorber]; }

 void   AddNeutralTrackL(double val, int absorber) { fNeutralTrackL[absorber] += val; fNeutralTrackL2[absorber] += val*val; }
 double GetNeutralTrackL(int absorber)  const    { return fNeutralTrackL[absorber];  }
 double GetNeutralTrackL2(int absorber) const    { return fNeutralTrackL2[absorber]; }

 void   AddGammas(double val)    { fNumGammas += val;    }
 double GetGammas()  const       { return fNumGammas;    }

 void   AddElectrons(double val) { fNumElectrons += val; }
 double GetElectrons()  const    { return fNumElectrons; }

 void   AddPositrons(double val) { fNumPositrons += val; }
 double GetPositrons()  const    { return fNumPositrons; }

 void   AddEdepInAbsorber(double val, int absorber)       { fEdepInAbsorber[absorber] += val; fEdepInAbsorber2[absorber] += val*val; }
 double GetEdepInAbsorber(int absorber)           const { return fEdepInAbsorber[absorber];  }
 double GetEdepInAbsorber2(int absorber)          const { return fEdepInAbsorber2[absorber]; }

 void   Clear();
 // add data after one primary particle finished tracking
 void   AddDataPerPrimary(TestEm3DataPerPrimary& data);

private:
  int     fNumAbsorbers;       // number of absorbers
  double  fNumChargedSteps;    // mean number of charged steps per primary in target
  double  fNumChargedSteps2;   // mean number of charged steps per primary square in target
  double  fNumNeutralSteps;    // mean number of neutral steps per primary in target
  double  fNumNeutralSteps2;   // mean number of neutral steps per primary square in target

  double  fNumGammas;          // mean number of secondary gamma particles per primary
  double  fNumElectrons;       // mean number of secondary electron particles per primary
  double  fNumPositrons;       // mean number of secondary positron particles per primary

  std::vector<double>  fEdepInAbsorber;     // mean energy deposit per primary in the target
  std::vector<double>  fEdepInAbsorber2;    // mean energy deposit per primary in the target square
  std::vector<double>  fChargedTrackL;      // mean number of charged track length per primary in target
  std::vector<double>  fChargedTrackL2;     // mean number of charged track length  per primary square in target
  std::vector<double>  fNeutralTrackL;      // mean number of neutral track length  per primary in target
  std::vector<double>  fNeutralTrackL2;     // mean number of neutral track length  per primary square in target
};



// Data structure per-event for TestEm3(contain as many per-primary data structures as number of primaries in one event)
class TestEm3DataPerEvent {
public:
  TestEm3DataPerEvent(int nprimperevent, int numabs);
 ~TestEm3DataPerEvent() {/*nothing to do*/}

  int GetNumberOfPrimaryPerEvent() const { return fNumPrimaryPerEvent; }
  void Clear();
  TestEm3DataPerPrimary& GetDataPerPrimary(int primindx) { return fPerPrimaryData[primindx]; }
  TestEm3DataPerEvent& operator+=(const TestEm3DataPerEvent& other);

private:
  int fNumPrimaryPerEvent;
  std::vector<TestEm3DataPerPrimary>  fPerPrimaryData; // as many as primary particle in an event
};



// Thread local data structure for TestEm3 to collecet/handle thread local multiple per-event data structures (as many
// per-event data structures as number of events are transported on the same time). Each of the currently transported
// events occupies one possible event-slot and per-event data can be identified by the index of this event-slot. This
// user defined thread local data needs to implement both the Merge and Clear methods for a given event-slot index:
// these methods are called when a data per-event,that corresponds to the completed event(with a given event-slot index),
// is merged from all threads.
class TestEm3ThreadDataEvents {
public:
  TestEm3ThreadDataEvents(int nevtbuffered, int nprimperevent, int numabs);
 ~TestEm3ThreadDataEvents() {/*nothing to do*/}

  void  Clear(int evtslotindx) { fPerEventData[evtslotindx].Clear(); }
//  void  Clear() {
//    for (int i=0; i<fNumBufferedEvents; ++i) Clear(i);
//  }
  bool  Merge(int evtslotindx, const TestEm3ThreadDataEvents &other);

  TestEm3DataPerEvent& GetDataPerEvent(int evtslotindx) { return fPerEventData[evtslotindx]; }
  const TestEm3DataPerEvent& GetDataPerEvent(int evtslotindx) const { return fPerEventData[evtslotindx]; }

private:
  int                                fNumBufferedEvents;
  std::vector<TestEm3DataPerEvent>   fPerEventData;
};


}       // namespace userapplication

#endif  // TESTEM3DATA_H
