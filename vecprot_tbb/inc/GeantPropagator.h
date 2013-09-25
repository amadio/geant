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
class TH1F;
class TH1S;
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

class GeantPropagator
{
public:
// data members to be made private
   Int_t       fNthreads;    // Number of threads
   Int_t       fNevents;     // Number of buffered
   Int_t       fNtotal;      // Total number of events
   Int_t      *fNtracks;     //[fNevents] Number of tracks {array of [fNevents]}

   tbb::atomic<Long64_t>    fNadded;
   tbb::atomic<Long64_t>    fNtransported; // Number of transported tracks

   Long64_t    fNsafeSteps;  // Number of fast steps within safety
   Long64_t    fNsnextSteps; // Number of steps where full snext computation is needed
   Int_t       fNprocesses;  // Number of active processes
   Int_t       fElossInd;    // Index of eloss process
   Int_t       fNstart;      // Cumulated initial number of tracks
   Int_t       fMaxTracks;   // Maximum number of tracks per event
   Int_t       fMaxThreads;  // Maximum number of threads
   Int_t       fDebugTrk;    // Track to debug
   Int_t       fMaxSteps;    // Maximum number of steps per track

   Int_t       fNperBasket;  // Number of tracks per basket
   Int_t       fMaxPerBasket; // Maximum number of tracks per basket
   Int_t       fMaxPerEvent; // Maximum number of tracks per event

   Double_t    fNaverage;    // Average number of tracks per event
   Double_t    fVertex[3];   // Vertex position
   Double_t    fEmin;        // Min energy threshold
   Double_t    fEmax;        // Max energy threshold
   Double_t    fBmag;        // Mag field

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
   tbb::atomic<Bool_t> fPrioritize;
   Bool_t fGarbageCollMode;

   Int_t fMinFeeder;
   Int_t fNevToPrioritize;
   Int_t fDispThr;
   Int_t fDispThrPriority;

   tbb::atomic<Int_t>* fEventsStatus;
   tbb::atomic<Int_t> fNimportedEvents;
   tbb::atomic<Int_t> fPriorityRange[2];

   tbb::atomic<Int_t> fCollsWaiting;
   tbb::atomic<Int_t> fTracksWaiting;
   TMutex fPropTaskLock;
   TMutex fDispTaskLock;

   tbb::atomic<Int_t> niter;
   TH1F *hnb;
   TH1F *htracks;

   void SetPriorityRange (Int_t min, Int_t max) { fPriorityRange[0]=min; fPriorityRange[1]=max; }

   static void* GlobalObserver (void* arg);
   tbb::concurrent_queue<Int_t> observerSigsQueue;
	tbb::atomic<Int_t> pnTasksRunning;
	tbb::atomic<Int_t> ppTasksRunning;
	tbb::atomic<Int_t> dTasksRunning;
	TH1S* numOfPNtasks;
	TH1S* numOfPPtasks;
	TH1S* numOfDtasks;
	TH1S* sizeOfFQ;
	TH1S* sizeOfPFQ;
	TH1S* sizeOfCQ;
	TH1S* numOfCollsPerTask;
	TH1S* numOfTracksInBasket;
	TH1S* numOfTracksInPriorBasket;
	TH1S* numOfTracksInColl;

   TH1I* numOfTracksTransportedInTime;
   TH1I* numOfTracksTransportedInIter;

public:
   GeantPropagator();
   virtual ~GeantPropagator();

   Int_t            AddTrack(GeantTrack *track);
   GeantTrack      *AddTrack(Int_t evslot);
   void             StopTrack(GeantTrack *track);
   Int_t            GetElossInd() const {return fElossInd;}
   Bool_t           LoadGeometry(const char *filename="geometry.root");
   Int_t            ImportTracks(Int_t nevents, Double_t average, Int_t startevent=0, Int_t startslot=0);
   void             Initialize();
   void             InjectCollection(GeantTrackCollection* inColl);
   static           GeantPropagator *Instance();
   void             PhysicsSelect(Int_t ntracks, Int_t *trackin);
   void             PrintParticles(Int_t *trackin, Int_t ntracks);
   PhysicsProcess  *Process(Int_t iproc) const {return fProcesses[iproc];}
   void             PropagatorGeom(const char *geomfile="geometry.root",
                                   Bool_t graphics=kFALSE,
                                   Bool_t single=kFALSE);
   void             SelectTracksForProcess(Int_t iproc, Int_t ntotransport, Int_t *particles, Int_t &ntodo, Int_t *parttodo);

private:
   GeantPropagator(const GeantPropagator&); // Not implemented
   GeantPropagator& operator=(const GeantPropagator&); // Not implemented
};
#endif
