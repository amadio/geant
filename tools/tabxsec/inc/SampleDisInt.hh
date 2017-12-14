// #include "G4Material.hh"
class G4Material;
class G4ParticleDefinition;
class G4VProcess;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class TFinState;
#include "G4LorentzVector.hh"

int SampDisInt(
               G4Material* material,
               G4ThreeVector *pos,
               G4DynamicParticle *dpart,
               G4double de,
               G4VProcess* proc,
               G4int    nevt,
               G4int    verbose,
               TFinState& fs);

G4double GetNuclearMass( G4int, G4int, G4int ); // G4Material* material );
const G4MaterialCutsCouple* FindMaterialCutsCouple( G4Material* mat );
void  checkBalance(const G4LorentzVector &porig, const G4LorentzVector &pmom, G4int bnum, G4DynamicParticle *secs, G4int n, G4LorentzVector &ptest, G4double &perr, G4int &berr);
G4bool rescaleEnergy(const G4LorentzVector &porig, G4DynamicParticle *secs, G4int n, G4double eleft, G4double etot);


// -------- Simple structure to hold one final states
struct Finstat_t {
  G4int npart;
  G4float weight;
  G4float kerma;
  G4float en;
  G4bool survived;
  G4int *pid;
  G4float *mom;
};

G4int SampleOne(G4Material* material,
                G4ThreeVector *pos,
                G4DynamicParticle *dpart,
                G4VProcess* proc,
                G4int    verbose,
                Finstat_t& fs);

constexpr char const* tStatus[6] = {"Alive: Continue the tracking",
  "StopButAlive: Invoke active rest physics processes and kill the current track afterward",
  "StopAndKill: Kill the current track",
  "KillTrackAndSecondaries: Kill the current track and also associated secondaries",
  "Suspend: Suspend the current track",
  "PostponeToNextEvent: Postpones the tracking of thecurrent track to the next event"};

