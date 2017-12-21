
#ifndef TESTEM5DATA_H
#define TESTEM5DATA_H

#include <vector>

namespace userapplication {

/**
 * @brief   User defined data structures for the TestEm5 GeantV application to store/handle scoring data.
 * @author  M Novak
 * @date    July 2017
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
 * An additional, global data structure, TestEm5Data is defined in this application to accumulate data that belongs to
 * an event/primary particle and available after the event is transported and the corresponding thread local data are
 * merged. Only one object is created from this global data structure and writing into this global data structure is
 * protected (more than one threads can finish an (different)event and try to write into this global data structure).
 */


class Hist;

// Data structure per-primary particle for TestEm5
class TestEm5DataPerPrimary {
public:
  TestEm5DataPerPrimary();
 ~TestEm5DataPerPrimary() { /*nothing to do*/}

  void   AddChargedStep()        { fNumChargedSteps += 1.;  }
  double GetChargedSteps() const { return fNumChargedSteps; }

  void   AddNeutralStep()        { fNumNeutralSteps += 1.;  }
  double GetNeutralSteps() const { return fNumNeutralSteps; }

  void   AddChargedTrackL(double val) { fChargedTrackL += val; }
  double GetChargedTrackL() const     { return fChargedTrackL; }

  void   AddNeutralTrackL(double val) { fNeutralTrackL += val; }
  double GetNeutralTrackL() const     { return fNeutralTrackL; }

  void   AddGamma()           { fNumGammas += 1.;     }
  double GetGammas()    const { return fNumGammas;    }

  void   AddElectron()        { fNumElectrons += 1.;  }
  double GetElectrons() const { return fNumElectrons; }

  void   AddPositron()        { fNumPositrons += 1.;  }
  double GetPositrons() const { return fNumPositrons; }

  void   AddPrimaryTransmitted()       { fNumPrimaryTrans += 1.;  }
  double GetPrimaryTransmitted() const { return fNumPrimaryTrans; }

  void   AddPrimaryReflected()         { fNumPrimaryRefl += 1.;  }
  double GetPrimaryReflected()   const { return fNumPrimaryRefl; }

  void   SetOneTransmitted(double val) { if (!fNumOneTrans) fNumOneTrans = val;   }
  double GetOneTransmitted() const     { return fNumOneTrans; }

  void   SetOneReflected(double val)   { if (!fNumOneRefl) fNumOneRefl = val;   }
  double GetOneReflected() const       { return fNumOneRefl;  }

  void   AddEdepInTarget(double val)       { fEdepInTarget += val; }
  double GetEdepInTarget()           const { return fEdepInTarget; }

  void   AddELeakPrimary(double val)         { fELeakPrimary += val;   }
  double GetELeakPrimary()             const { return fELeakPrimary;   }

  void   AddELeakSecondary(double val)       { fELeakSecondary += val; }
  double GetELeakSecondary()           const { return fELeakSecondary; }

  void   Clear();
  TestEm5DataPerPrimary& operator+=(const TestEm5DataPerPrimary& other);

private:
  double  fNumChargedSteps;     // mean number of charged steps per primary in target
  double  fNumNeutralSteps;     // mean number of neutral steps per primary in target
  double  fChargedTrackL;       // mean number of charged track length per primary in target
  double  fNeutralTrackL;       // mean number of neutral track length per primary in target
  double  fNumGammas;           // mean number of secondary gamma particles per primary
  double  fNumElectrons;        // mean number of secondary electron particles per primary
  double  fNumPositrons;        // mean number of secondary positron particles per primary
  double  fNumPrimaryTrans;     // mean number of primary particles transmitted
  double  fNumPrimaryRefl;      // mean number of primary particles reflected
  double  fNumOneTrans;         // mean number of events with at last one particle (same charge as primary) transmitted
  double  fNumOneRefl;          // mean number of events with at last one particle (same charge as primary) reflected
  double  fEdepInTarget;        // mean energy deposit per primary in the target
  double  fELeakPrimary;        // mean primary particle energy leakage per primary particles
  double  fELeakSecondary;      // mean secondary particle energy leakage per primary particles

};

// Global data structure to accumulate per-primary data during the simulation. The only one object from this class
// will be updated each time an event(with the corresponding primaries) is transported.
class TestEm5Data {
public:
  TestEm5Data();
 ~TestEm5Data() { /*nothing to do*/ }

 void   AddChargedSteps(double val) { fNumChargedSteps += val; fNumChargedSteps2 += val*val; }
 double GetChargedSteps()  const    { return fNumChargedSteps;  }
 double GetChargedSteps2() const    { return fNumChargedSteps2; }

 void   AddNeutralSteps(double val) { fNumNeutralSteps += val; fNumNeutralSteps2 += val*val; }
 double GetNeutralSteps() const     { return fNumNeutralSteps;  }
 double GetNeutralSteps2() const    { return fNumNeutralSteps2; }

 void   AddChargedTrackL(double val) { fChargedTrackL += val; fChargedTrackL2 += val*val; }
 double GetChargedTrackL()  const    { return fChargedTrackL;  }
 double GetChargedTrackL2() const    { return fChargedTrackL2; }

 void   AddNeutralTrackL(double val) { fNeutralTrackL += val; fNeutralTrackL2 += val*val; }
 double GetNeutralTrackL()  const    { return fNeutralTrackL;  }
 double GetNeutralTrackL2() const    { return fNeutralTrackL2; }

 void   AddGammas(double val)    { fNumGammas += val;    }
 double GetGammas()  const       { return fNumGammas;    }

 void   AddElectrons(double val) { fNumElectrons += val; }
 double GetElectrons()  const    { return fNumElectrons; }

 void   AddPositrons(double val) { fNumPositrons += val; }
 double GetPositrons()  const    { return fNumPositrons; }

 void   AddPrimaryTransmitted(double val)       { fNumPrimaryTrans += val; }
 double GetPrimaryTransmitted()           const { return fNumPrimaryTrans; }

 void   AddPrimaryReflected(double val)         { fNumPrimaryRefl += val; }
 double GetPrimaryReflecred()             const { return fNumPrimaryRefl; }

 void   AddOneTransmitted(double val)       { if (val>0.) fNumOneTrans += 1; }
 double GetOneTransmitted()           const { return fNumOneTrans; }

 void   AddOneReflected(double val)         { if (val>0.) fNumOneRefl += 1.;  }
 double GetOneReflected()             const { return fNumOneRefl;  }

 void   AddEdepInTarget(double val)       { fEdepInTarget += val; fEdepInTarget2 += val*val; }
 double GetEdepInTarget()           const { return fEdepInTarget;  }
 double GetEdepInTarget2()          const { return fEdepInTarget2; }

 void   AddELeakPrimary(double val)         { fELeakPrimary += val; fELeakPrimary2 += val*val; }
 double GetELeakPrimary()             const { return fELeakPrimary;  }
 double GetELeakPrimary2()            const { return fELeakPrimary2; }

 void   AddELeakSecondary(double val)       { fELeakSecondary += val; fELeakSecondary2 += val*val; }
 double GetELeakSecondary()           const { return fELeakSecondary;  }
 double GetELeakSecondary2()          const { return fELeakSecondary2; }

 void   Clear();
 // add data after one primary particle finished tracking
 void   AddDataPerPrimary(TestEm5DataPerPrimary& data);

private:
  double  fNumChargedSteps;    // mean number of charged steps per primary in target
  double  fNumChargedSteps2;   // mean number of charged steps per primary square in target
  double  fNumNeutralSteps;    // mean number of neutral steps per primary in target
  double  fNumNeutralSteps2;   // mean number of neutral steps per primary square in target

  double  fChargedTrackL;      // mean number of charged track length per primary in target
  double  fChargedTrackL2;     // mean number of charged track length  per primary square in target
  double  fNeutralTrackL;      // mean number of neutral track length  per primary in target
  double  fNeutralTrackL2;     // mean number of neutral track length  per primary square in target

  double  fNumGammas;          // mean number of secondary gamma particles per primary
  double  fNumElectrons;       // mean number of secondary electron particles per primary
  double  fNumPositrons;       // mean number of secondary positron particles per primary

  double  fNumPrimaryTrans;    // mean number of parimary particle transported
  double  fNumPrimaryRefl;     // mean number of parimary particle reflected
  double  fNumOneTrans;        // mean number of events with at last one particle (same charge as primary) transmitted
  double  fNumOneRefl;         // mean number of events with at last one particle (same charge as primary) reflected

  double  fEdepInTarget;       // mean energy deposit per primary in the target
  double  fEdepInTarget2;      // mean energy deposit per primary in the target square
  double  fELeakPrimary;       // mean primary particle energy leakage per primary particles
  double  fELeakPrimary2;      // mean primary particle energy leakage per primary particles square
  double  fELeakSecondary;     // mean secondary particle energy leakage per primary particles
  double  fELeakSecondary2;    // mean secondary particle energy leakage per primary particles square
};



// Data structure per-event for TeatEm5(contain as many per-primary data structures as number of primaries in one event)
class TestEm5DataPerEvent {
public:
  TestEm5DataPerEvent(int nprimperevent);
 ~TestEm5DataPerEvent() {/*nothing to do*/}

  int GetNumberOfPrimaryPerEvent() const { return fNumPrimaryPerEvent; }
  void Clear();
  TestEm5DataPerPrimary& GetDataPerPrimary(int primindx) { return fPerPrimaryData[primindx]; }
  TestEm5DataPerEvent& operator+=(const TestEm5DataPerEvent& other);

private:
  int fNumPrimaryPerEvent;
  std::vector<TestEm5DataPerPrimary>  fPerPrimaryData; // as many as primary particle in an event
};



// Thread local data structure for TestEm5 to collecet/handle thread local multiple per-event data structures (as many
// per-event data structures as number of events are transported on the same time). Each of the currently transported
// events occupies one possible event-slot and per-event data can be identified by the index of this event-slot. This
// user defined thread local data needs to implement both the Merge and Clear methods for a given event-slot index:
// these methods are called when a data per-event,that corresponds to the completed event(with a given event-slot index),
// is merged from all threads.
class TestEm5ThreadDataEvents {
public:
  TestEm5ThreadDataEvents(int nevtbuffered, int nprimperevent);
 ~TestEm5ThreadDataEvents() {/*nothing to do*/}

  void  Clear(int evtslotindx) { fPerEventData[evtslotindx].Clear(); }
//  void  Clear() {
//    for (int i=0; i<fNumBufferedEvents; ++i) Clear(i);
//  }
  bool  Merge(int evtslotindx, const TestEm5ThreadDataEvents &other);

  TestEm5DataPerEvent& GetDataPerEvent(int evtslotindx) { return fPerEventData[evtslotindx]; }
  const TestEm5DataPerEvent& GetDataPerEvent(int evtslotindx) const { return fPerEventData[evtslotindx]; }

private:
  int                                fNumBufferedEvents;
  std::vector<TestEm5DataPerEvent>   fPerEventData;
};


// Thread local data structure for TestEm5 to collecet/handle thread local run-global data structures. The user defined
// needs to implement both the Merge and Clear methods: these methods are called when the simulation is completed and
// these thread local run-global data are merged from the working threads.
class TestEm5ThreadDataRun {
public:
  TestEm5ThreadDataRun();
 ~TestEm5ThreadDataRun();

 void   CreateHisto1(int nbins, double min, double max);
 Hist*  GetHisto1() const { return fHisto1; }
 // nothing to clear: per-thread histograms will be merged at the end of the run
 void   Clear(int /*evtslotindx*/) {}
 bool   Merge(int /*evtslotindx*/, const TestEm5ThreadDataRun& other);

private:
  // simple histogram to store angular distributions per working-threads (they will be merged at the end of run)
  Hist    *fHisto1;
};


}       // namespace userapplication

#endif  // TESTEM5DATA_H
