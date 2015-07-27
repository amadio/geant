#ifndef GEANT_TRACK
#define GEANT_TRACK

#include "globals.h"
#include "TMutex.h"

using std::numeric_limits;

class TGeoBranchArray;
class GeantVolumeBasket;

const double kB2C = -0.299792458e-3;
enum TrackStatus_t { kAlive, kKilled, kBoundary };

//______________________________________________________________________________
class GeantTrack {
public:
  Int_t event;          // event number
  Int_t evslot;         // event slot
  Int_t particle;       // index of corresponding particle
  Int_t pdg;            // particle pdg code
  Int_t fGVcode;        // GV particle code
  Species_t species;    // particle species
  TrackStatus_t status; // track status
  Int_t charge;         // particle charge
  double mass;        // particle mass
  Int_t process;        // current process
  double xpos;        // position
  double ypos;
  double zpos;
  double px; // momentum
  double py;
  double pz;
  double e;                // energy
  double pstep;            // selected physical step
  double step;             // current step
  double snext;            // straight distance to next boundary
  double safety;           // safe distance to any boundary
  Bool_t frombdr;            // true if starting from boundary
  Int_t izero;               // number of small steps used to catch errors
  Int_t nsteps;              // number of steps made
  TGeoBranchArray *path;     // path for this particle in the geometry
  TGeoBranchArray *nextpath; // path for next volume
  Bool_t pending;

  GeantTrack()
      : event(-1), evslot(-1), particle(-1), pdg(0), fGVcode(0), species(kHadron), status(kAlive), charge(0), mass(0),
        process(-1), xpos(0), ypos(0), zpos(0), px(0), py(0), pz(0), e(0), pstep(1.E20), step(0), snext(0), safety(0),
        frombdr(false), izero(0), nsteps(0), path(0), nextpath(0), pending(false) {}
  GeantTrack(const GeantTrack &other);
  GeantTrack &operator=(const GeantTrack &other);
  GeantTrack(Int_t ipdg);
  virtual ~GeantTrack();
  double Curvature() const;
  void Direction(double dir[3]);
  Bool_t IsAlive() const { return (status != kKilled); }
  Bool_t IsOnBoundary() const { return (status == kBoundary); }
  void Kill() { status = kKilled; }
  void Print(Int_t trackindex = 0) const;
  GeantVolumeBasket *PropagateInField(double step, Bool_t checkcross, Int_t itr);
  GeantVolumeBasket *PropagateStraight(double step, Int_t itrack);
  double Pt() const { return sqrt(px * px + py * py); }
  double P() const { return sqrt(px * px + py * py + pz * pz); }
  double Gamma() const { return mass ? e / mass : numeric_limits<double>.max(); }
  double Beta() const { return P() / e; }

  void Reset();
};

#endif // GEANT_TRACK
