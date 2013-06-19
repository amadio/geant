#ifndef GEANT_PROPAGATOR
#define GEANT_PROPAGATOR

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TMutex
#include "TMutex.h"
#endif

#include <vector>

class TRandom;
class TArrayI;
class TF1;
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
class GeantThreadData;

class GeantPropagator : public TObject
{
public:
// data members to be made private
   Int_t       fNthreads;    // Number of threads
   Int_t       fNevents;     // Number of buffered
   Int_t       fNtotal;      // Total number of events
   Int_t      *fNtracks;     //[fNevents] Number of tracks {array of [fNevents]}
   Long64_t    fNtransported; // Number of transported tracks
   Long64_t    fNsafeSteps;  // Number of fast steps within safety
   Long64_t    fNsnextSteps; // Number of steps where full snext computation is needed
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
   TMutex      fTracksLock;  // Mutex for adding tracks
   
   WorkloadManager *fWMgr;   // Workload manager
   GeantOutput *
               fOutput;      // Output object
   
   TF1             *fKineTF1;   //
   TTree           *fOutTree;   // Output tree
   TFile           *fOutFile;   // Output file
   TStopwatch      *fTimer;     // Timer
   
   PhysicsProcess **fProcesses; //![fNprocesses] Array of processes
   GeantTrack     **fTracks;    //![fMaxTracks]  Array of tracks
   std::vector<GeantTrack> fTracksStorage; // buffer for tracks allocation

   GeantEvent     **fEvents;    //![fNevents]    Array of events
   Int_t           *fTracksPerBasket; //![fNthreads]
   GeantTrackCollection **fCollections; //![fNthreads]  Track collections
   UInt_t          *fWaiting;           //![fNthreads] Threads in waiting flag
   GeantThreadData **fThreadData; //![fNthreads]
   
   static GeantPropagator *fgInstance;
public:
   GeantPropagator();
   virtual ~GeantPropagator();
   
   Int_t            AddTrack(GeantTrack *track);
   GeantTrack      *AddTrack(Int_t evslot);
   void             StopTrack(GeantTrack *track);
   Int_t            GetElossInd() const {return fElossInd;}
   UInt_t           GetNwaiting() const;
   Bool_t           LoadGeometry(const char *filename="geometry.root");
   GeantVolumeBasket *
                    ImportTracks(Int_t nevents, Double_t average, Int_t startevent=0, Int_t startslot=0);
   void             Initialize();
   void             InjectCollection(Int_t tid);
   GeantBasket     *InjectBasket(GeantBasket *basket);
   static 
   GeantPropagator *Instance();
   void             PhysicsSelect(Int_t ntracks, Int_t *trackin, Int_t tid);
   void             PrintParticles(Int_t *trackin, Int_t ntracks, Int_t tid);
   PhysicsProcess  *Process(Int_t iproc) const {return fProcesses[iproc];}
   void             PropagatorGeom(const char *geomfile="geometry.root",
                                   Int_t nthreads=4,
                                   Bool_t graphics=kFALSE,
                                   Bool_t single=kFALSE);
   void             SelectTracksForProcess(Int_t iproc, Int_t ntotransport, Int_t *particles, Int_t &ntodo, Int_t *parttodo);
   ClassDef(GeantPropagator, 1)
};
#endif
