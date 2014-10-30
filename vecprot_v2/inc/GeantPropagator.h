#ifndef GEANT_PROPAGATOR
#define GEANT_PROPAGATOR

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef GEANT_TRACK
#include "GeantTrack.h"
#endif

#ifndef ROOT_TMutex
#include "TMutex.h"
#endif

#include <vector>

#if __cplusplus >= 201103L
#include <atomic>
#endif

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
class GeantBasketMgr;
class WorkloadManager;
class GeantThreadData;
class GeantVApplication;

class GeantPropagator : public TObject
{
public:
// data members to be made private
   Int_t       fNthreads;    // Number of threads
   Int_t       fNevents;     // Number of buffered
   Int_t       fNtotal;      // Total number of events
#if __cplusplus >= 201103L
   std::atomic<Long64_t>    fNtransported; // Number of transported tracks
#else
   Long64_t    fNtransported; // Number of transported tracks   
#endif   
   Long64_t    fNprimaries;  // Number of primary tracks
   Long64_t    fNsafeSteps;  // Number of fast steps within safety
   Long64_t    fNsnextSteps; // Number of steps where full snext computation is needed
   Long64_t    fNphysSteps;  // Number of steps to physics process
   Int_t       fNprocesses;  // Number of active processes
   Int_t       fElossInd;    // Index of eloss process
   Int_t       fNstart;      // Cumulated initial number of tracks
   Int_t       fMaxTracks;   // Maximum number of tracks per event
   Int_t       fMaxThreads;  // Maximum number of threads
   Int_t       fNminThreshold; // Threshold for starting transporting a basket
   Int_t       fDebugTrk;    // Track to debug
   Int_t       fMaxSteps;    // Maximum number of steps per track
   Int_t       fNperBasket;  // Number of tracks per basket
   Int_t       fMaxPerBasket; // Maximum number of tracks per basket
   Int_t       fMaxPerEvent; // Maximum number of tracks per event
   Int_t       fMaxDepth;    // Maximum geometry depth
   
   Double_t    fNaverage;    // Average number of tracks per event
   Double_t    fVertex[3];   // Vertex position
   Double_t    fEmin;        // Min energy threshold
   Double_t    fEmax;        // Max energy threshold
   Double_t    fBmag;        // Mag field
   
   Bool_t      fUsePhysics;  // Enable/disable physics
   Bool_t      fUseDebug;    // Use debug mode
   Bool_t      fUseGraphics; // graphics mode
   Bool_t      fTransportOngoing; // Flag for ongoing transport
   Bool_t      fSingleTrack; // Use single track transport mode
   Bool_t      fFillTree;    // Enable I/O
   Bool_t      fUseMonitoring; // Monitoring thread
   TMutex      fTracksLock;  // Mutex for adding tracks
   
   WorkloadManager *fWMgr;   // Workload manager
   GeantVApplication *fApplication; // User application
   GeantOutput      *fOutput;      // Output object
   
   TTree           *fOutTree;   // Output tree
   TFile           *fOutFile;   // Output file
   TStopwatch      *fTimer;     // Timer
   
   PhysicsProcess  *fProcess;   // For now the only generic process pointing to the tabulated physics
//   PhysicsProcess **fProcesses; //![fNprocesses] Array of processes
   GeantTrack_v    *fStoredTracks;    //! Stored array of tracks (history?)

   // Data per event
   Int_t           *fNtracks;   //[fNevents] Number of tracks {array of [fNevents]}
   GeantEvent     **fEvents;    //![fNevents]    Array of events

   GeantThreadData **fThreadData; //![fNthreads]
   
   static GeantPropagator *fgInstance;
public:
   GeantPropagator();
   virtual ~GeantPropagator();
   
   // Temporary track for the current caller thread
#if __cplusplus >= 201103L
   Long64_t         GetNtransported() const {return fNtransported.load();}
#else
   Long64_t         GetNtransported() const {return fNtransported;}
#endif   
   GeantTrack      &GetTempTrack(Int_t tid=-1);
   Int_t            AddTrack(GeantTrack &track);
   Int_t            DispatchTrack(GeantTrack &track);
   void             StopTrack(const GeantTrack_v &tracks, Int_t itr);
   Int_t            GetElossInd() const {return fElossInd;}
   Bool_t           LoadGeometry(const char *filename="geometry.root");
#if USE_VECGEOM_NAVIGATOR == 1
   Bool_t           LoadVecGeomGeometry();
#endif
   Int_t            ImportTracks(Int_t nevents, Double_t average, Int_t startevent=0, Int_t startslot=0);
   void             Initialize();
   static 
   GeantPropagator *Instance(Int_t ntotal=0, Int_t nbuffered=0);
   void             ProposeStep(Int_t ntracks, GeantTrack_v &tracks, Int_t tid);
   void             ApplyMsc(Int_t ntracks, GeantTrack_v &tracks, Int_t tid);
//   PhysicsProcess  *Process(Int_t iproc) const {return fProcesses[iproc];}
   PhysicsProcess  *Process() const {return fProcess;}
   void             PropagatorGeom(const char *geomfile="geometry.root",
                                   Int_t nthreads=4,
                                   Bool_t graphics=kFALSE,
                                   Bool_t single=kFALSE);

private:
   GeantPropagator(const GeantPropagator&); // Not implemented
   GeantPropagator& operator=(const GeantPropagator&); // Not implemented
   ClassDef(GeantPropagator, 1)
};
#endif
