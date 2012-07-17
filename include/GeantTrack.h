#ifndef GEANT_TRACK
#define GEANT_TRACK

#include "globals.h"
#include "TMath.h"
#include "TMutex.h"

class TGeoBranchArray;
class GeantVolumeBasket;

const Double_t kB2C = -0.299792458e-3;
enum TrackStatus_t {kAlive, kKilled, kBoundary};

//______________________________________________________________________________
struct GeantTrack {
   Int_t    event;     // event number
   Int_t    evslot;    // event slot
   Int_t    particle;  // index of corresponding particle
   Int_t    pdg;       // particle pdg code
   Species_t species;  // particle species
   TrackStatus_t status; // track status
   Int_t    charge;    // particle charge
   Double_t mass;      // particle mass
   Int_t    process;   // current process
   Double_t xpos;      // position
   Double_t ypos;
   Double_t zpos;
   Double_t px;        // momentum
   Double_t py;
   Double_t pz;
   Double_t e;         // energy
   Double_t pstep;     // selected physical step
   Double_t step;      // current step
   Double_t snext;     // straight distance to next boundary
   Double_t safety;    // safe distance to any boundary
   Bool_t   frombdr;   // true if starting from boundary
   Int_t    izero;     // number of small steps used to catch errors
   Int_t    nsteps;    // number of steps made
   TGeoBranchArray *path; // path for this particle in the geometry
   TGeoBranchArray *nextpath; // path for next volume
   Bool_t   pending;
   
   GeantTrack() : event(-1),evslot(-1), particle(-1),pdg(0),species(kHadron),status(kAlive),
                  charge(0),mass(0),process(-1),xpos(0),ypos(0),zpos(0),
                  px(0),py(0),pz(0),e(0), pstep(1.E20), step(0), snext(0), 
                  safety(0), frombdr(false), izero(0), nsteps(0), path(0) ,nextpath(0), pending(false){}
   GeantTrack(Int_t ipdg);
   ~GeantTrack();
   Double_t           Curvature() {return (charge)?TMath::Abs(kB2C*gPropagator->fBmag/Pt()):0.;}
   void               Direction(Double_t dir[3]);
   Bool_t             IsAlive() const {return (status != kKilled);}
   Bool_t             IsOnBoundary() const {return (status == kBoundary);}
   void               Kill()        {status = kKilled;}
   void               Print(Int_t trackindex=0) const;
   Bool_t             PropagateInFieldSingle(Double_t step, Bool_t checkcross, Int_t itr);
   GeantVolumeBasket *PropagateInField(Double_t step, Bool_t checkcross, Int_t itr);
   GeantVolumeBasket *PropagateStraight(Double_t step, Int_t itrack);
   Double_t           Pt()    const {return TMath::Sqrt(px*px+py*py);}
   Double_t           P()     const {return TMath::Sqrt(px*px+py*py+pz*pz);}
   Double_t           Gamma() const {return mass?e/mass:TMath::Limits<double>::Max();}
   Double_t           Beta()  const {return P()/e;}
};

//______________________________________________________________________________
struct GeantEvent {
   Int_t  event;        // event number
   Int_t  slot;         // fixed slot number
   Int_t  ntracks;      // number of tracks
   Int_t  ndone;        // number of done tracks
   TMutex the_mutex;    // mutex for this

   GeantEvent() : event(0), slot(0), ntracks(0), ndone(0), the_mutex() {}
   GeantEvent(Int_t ntr) : event(0), slot(0), ntracks(ntr), ndone(0), the_mutex() {}
   ~GeantEvent() {}
   
   void               AddTrack();
   void               Reset() {ntracks = ndone = 0;}
   void               StopTrack();
   Bool_t             Transported() const {return ((ntracks>0) && (ntracks==ndone));}
};

#endif
