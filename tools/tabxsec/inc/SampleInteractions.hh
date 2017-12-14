// #include "G4Material.hh"
class G4Material;
class G4ParticleDefinition;
class G4VProcess;
class G4DynamicParticle;
class TFinState;

    int SampleInteractions( G4Material* material,
			    const G4ThreeVector *pos,
                            const G4ParticleDefinition* part,
                            G4VProcess* proc,
                            G4double energy,
                            G4double sigmae,
                            G4double step,
                            G4int    noInteractions,
                            G4int    verbose
       );

    int SampleDiscreteInteractions(  G4Material* material,
		       const G4ThreeVector *pos,
                       const G4ParticleDefinition* part,
                       G4VProcess* proc,
                       G4double energy,
                       G4double sigmae,
                       G4double step,
                       G4int    noInteractions,
                       G4int    verbose
                       );

