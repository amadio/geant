// #include "G4Material.hh"
class G4Material;
class G4ParticleDefinition;
class G4VProcess;
class G4DynamicParticle;
class G4MaterialCutsCouple;
class TFinState;

int SampDisInt(
               G4Material* material,
               G4ThreeVector *pos,
               G4DynamicParticle *dpart,
               G4VProcess* proc,
               G4int    nevt,
               G4int    verbose,
               TFinState& fs);

G4double GetNuclearMass( G4int, G4int, G4int ); // G4Material* material );
const G4MaterialCutsCouple* FindMaterialCutsCouple( G4Material* mat );

// -------- Simple structure to hold one final states
struct Finstat_t {
  G4bool survived;
  G4int npart;
  G4float kerma;
  G4float weight;
  G4int *pid;
  G4float (*mom)[3];
};

G4int SampleOne(G4Material* material,
                G4ThreeVector *pos,
                G4DynamicParticle *dpart,
                G4VProcess* proc,
                G4int    verbose,
                Finstat_t& fs);

static const char* tStatus[6] = {"Alive: Continue the tracking",
  "StopButAlive: Invoke active rest physics processes and kill the current track afterward",
  "StopAndKill: Kill the current track",
  "KillTrackAndSecondaries: Kill the current track and also associated secondaries",
  "Suspend: Suspend the current track",
  "PostponeToNextEvent: Postpones the tracking of thecurrent track to the next event"};
