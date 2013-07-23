// #include "G4Material.hh"
class G4Material;
class G4ParticleDefinition;
class G4VProcess;
class G4DynamicParticle;
class TFinState;

    int SampDisInt(
		   G4Material* material, 
		   G4ThreeVector *pos,
		   G4DynamicParticle *dpart,
		   G4VProcess* proc,
		   G4int    nevt,
		   G4int    verbose,
		   TFinState& fs);
     G4double GetNuclearMass( G4int Z, G4int N );
