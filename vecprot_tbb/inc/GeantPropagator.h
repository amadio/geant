#ifndef GEANT_PROPAGATOR
#define GEANT_PROPAGATOR

#include "Rtypes.h"

#ifndef ROOT_TMutex
#include "TMutex.h"
#endif

#include <vector>

#include "tbb/atomic.h"
#include "tbb/concurrent_queue.h"
#include "tbb/enumerable_thread_specific.h"

#include "GeantThreadData.h"

class TRandom;
class TArrayI;
class TF1;
class TH1I;
class TTree;
class TFile;
class TStopwatch;
class TGeoHMatrix;
class TGeoRotation;
class TGeoVolume;
class TGeoHelix;
class PhysicsProcess;
class GeantTrack;
struct GeantEvent;
class GeantBasket;
class GeantOutput;
class GeantVolumeBasket;
class WorkloadManager;
class GeantTrackCollection;

typedef tbb::enumerable_thread_specific<GeantThreadData> PerThread;

#include "tbb/task_scheduler_observer.h"
class myObserver : public tbb::task_scheduler_observer
{
public:
	void on_scheduler_entry (bool is_worker);
	void on_scheduler_exit (bool is_worker);
};

class GeantPropagator
{
public:
// data members to be made private
   int       fNthreads;    // Number of threads
   int       fNevents;     // Number of buffered
   int       fNtotal;      // Total number of events
   int      *fNtracks;     //[fNevents] Number of tracks {array of [fNevents]}

   tbb::atomic<Long64_t>    fNadded;
   tbb::atomic<Long64_t>    fNtransported; // Number of transported tracks

   Long64_t    fNsafeSteps;  // Number of fast steps within safety
   Long64_t    fNsnextSteps; // Number of steps where full snext computation is needed
   int       fNprocesses;  // Number of active processes
   int       fElossInd;    // Index of eloss process
   int       fNstart;      // Cumulated initial number of tracks
   int       fMaxTracks;   // Maximum number of tracks per event
   int       fMaxThreads;  // Maximum number of threads
   int       fDebugTrk;    // Track to debug
   int       fMaxSteps;    // Maximum number of steps per track

   int       fNperBasket;  // Number of tracks per basket
   int       fMaxPerBasket; // Maximum number of tracks per basket
   int       fMaxPerEvent; // Maximum number of tracks per event

   double    fNaverage;    // Average number of tracks per event
   double    fVertex[3];   // Vertex position
   double    fEmin;        // Min energy threshold
   double    fEmax;        // Max energy threshold
   double    fBmag;        // Mag field

   Bool_t      fUsePhysics;  // Enable/disable physics
   Bool_t      fUseDebug;    // Use debug mode
   Bool_t      fUseGraphics; // graphics mode
   Bool_t      fSingleTrack; // Use single track transport mode
   Bool_t      fFillTree;    // Enable I/O
   TMutex      fTracksLock;  // Mutex for adding tracks

   WorkloadManager *fWMgr;   // Workload manager
   GeantOutput* fOutput;      // Output object

   TF1             *fKineTF1;   //
   TTree           *fOutTree;   // Output tree
   TFile           *fOutFile;   // Output file
   TStopwatch      *fTimer;     // Timer

   PhysicsProcess **fProcesses; //![fNprocesses] Array of processes
   GeantTrack     **fTracks;    //![fMaxTracks]  Array of tracks
   std::vector<GeantTrack> fTracksStorage; // buffer for tracks allocation

   GeantEvent     **fEvents;    //![fNevents]    Array of events

   PerThread fTBBthreadData;

   static GeantPropagator *fgInstance;

public:
   Bool_t fGarbageCollMode;

   int fMinFeeder;
   int fNevToPrioritize;
   int fDispThr;
   int fDispThrPriority;

   tbb::atomic<int>* fEventsStatus;
   tbb::atomic<int> fNimportedEvents;
   tbb::atomic<int> fPriorityRange[2];

   tbb::atomic<int> fCollsWaiting;
   tbb::atomic<int> fTracksWaiting;
   TMutex fPropTaskLock;
   TMutex fDispTaskLock;

   void SetPriorityRange (int min, int max) { fPriorityRange[0]=min; fPriorityRange[1]=max; }

   static void* GlobalObserver (void* arg);
   tbb::concurrent_queue<int> observerSigsQueue;
	myObserver propObserver;

   tbb::atomic<int> pnTasksTotal;
   tbb::atomic<int> ppTasksTotal;
   tbb::atomic<int> dTasksTotal;

	tbb::atomic<int> pnTasksRunning;
	tbb::atomic<int> ppTasksRunning;
	tbb::atomic<int> dTasksRunning;
   tbb::atomic<int> niter;
   tbb::atomic<int> niter2;
   tbb::atomic<int> niter3;
   tbb::atomic<int> niter4;

   // In time
   TH1I* numOfTracksTransportedInTime;
	TH1I* numOfPNtasks;
	TH1I* numOfPPtasks;
	TH1I* numOfDtasks;
	TH1I* sizeOfFQ;
	TH1I* sizeOfPFQ;
	TH1I* sizeOfCQ;

   // In iter
   TH1I *hnb;
   TH1I *htracks;
   TH1I* numOfTracksTransportedInIter;
   /*TH1I* numOfnPropTasks;
   TH1I* numOfpPropTasks;
   TH1I* numOfDispTasks;*/

   // Statistics
	TH1I* numOfCollsPerTask;
	TH1I* numOfTracksInBasket;
	TH1I* numOfTracksInPriorBasket;
	TH1I* numOfTracksInColl;

public:
   GeantPropagator();
   virtual ~GeantPropagator();

   int            AddTrack(GeantTrack *track);
   GeantTrack      *AddTrack(int evslot);
   void             StopTrack(GeantTrack *track);
   int            GetElossInd() const {return fElossInd;}
   Bool_t           LoadGeometry(const char *filename="geometry.root");
   int            ImportTracks(int nevents, double average, int startevent=0, int startslot=0);
   void             Initialize();
   void             InjectCollection(GeantTrackCollection* inColl);
   static           GeantPropagator *Instance();
   void             PhysicsSelect(int ntracks, int *trackin);
   void             PrintParticles(int *trackin, int ntracks);
   PhysicsProcess  *Process(int iproc) const {return fProcesses[iproc];}
   void             PropagatorGeom(const char *geomfile="geometry.root",
                                   Bool_t graphics=kFALSE,
                                   Bool_t single=kFALSE);
   void             SelectTracksForProcess(int iproc, int ntotransport, int *particles, int &ntodo, int *parttodo);

private:
   GeantPropagator(const GeantPropagator&); // Not implemented
   GeantPropagator& operator=(const GeantPropagator&); // Not implemented
};
#endif
