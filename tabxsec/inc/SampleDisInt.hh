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

