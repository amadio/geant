//---- ROOT stuff
#include <TH1F.h>
#include <TMath.h>
#include <TFile.h>
#include <TEXsec.h>
#include <TEFstate.h>
#include <TRandom.h>
#include <TPartIndex.h>
#include <TClass.h>
// #include <TVirtualStreamerInfo.h>
#include <TDatabasePDG.h>
#include <TMap.h>
#include <TObjString.h>
#include <TFinState.h>
#include "TPDecay.h"

#include <unistd.h>

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "TParticlePDG.h"
#include "TPartIndex.h"

void RegisterG4Particles()
{
  // Store the particle table into PartIndex
  G4ParticleDefinition* particleDef;


  G4ParticleTable* theParticleTable= G4ParticleTable::GetParticleTable();
  G4int np=theParticleTable->size();
  G4int *pPDG = new G4int[np];

  int   npGood=0;
  for(G4int i=0; i<np; ++i) {
        // G4cout << " i = " << i << G4endl;
     particleDef= theParticleTable->GetParticle(i);
        // G4cout << " Particle[ " << i << " ]: address " << particle ; 
           // <<  G4endl;
     // if( (particle != 0) && (particle->GetBaryonNumber()<5) ){
     TParticlePDG *AddParticleToPdgDatabase(const G4String& name,
                                            G4ParticleDefinition* particleDefinition);

     TParticlePDG* pdpdg_i = AddParticleToPdgDatabase(particleDef->GetParticleName(),particleDef);

     if( pdpdg_i ){
       pPDG[npGood]=pdpdg_i->PdgCode();
       npGood++;
     }else{
       std::cerr << "Refused to Add Particle to Root Pdg Database "
           << " for particle " << particleDef->GetParticleName() << std::endl;
     }
  }
  TPartIndex::I()->SetPartTable(pPDG,np);
}

static TMap partDict;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//_____________________________________________________________________________
TParticlePDG *AddParticleToPdgDatabase(const G4String& name,
                                       G4ParticleDefinition* particleDefinition)
{
  // Add the particle definition in TDatabasePDG
  // Code by Ivana Hrivnacova
  
  // Return if particle was already added
  G4int pdgEncoding = particleDefinition->GetPDGEncoding();
  TParticlePDG *particlePDG=0;
  if(pdgEncoding)
    particlePDG
    = TDatabasePDG::Instance()->GetParticle(pdgEncoding);
  else {
    TObjString *conversion = (TObjString*) partDict.GetValue((const char *)particleDefinition->GetParticleName());
    if( conversion )
      particlePDG
	    = TDatabasePDG::Instance()->GetParticle(conversion->String().Data());
    else
      particlePDG
	    = TDatabasePDG::Instance()->GetParticle((const char*) name);
  }
  if ( particlePDG )  return particlePDG;
  
  // Get particle data
  G4String g4Name = particleDefinition->GetParticleName();
  G4int pdgQ = G4int(particleDefinition->GetPDGCharge()/CLHEP::eplus);
  // !! here we do not save dynamic charge but the static one
  G4String g4Type = particleDefinition->GetParticleType();
  G4String rootType = g4Type;
  if ( g4Type == "nucleus" ||  g4Type == "anti_nucleus") rootType="Ion";
  
  G4bool verbose = FALSE;
  if (verbose) {
    G4cout << "Adding particle to TDatabasePDG " << G4endl;
    G4cout << "   name:   " << g4Name << G4endl;
    G4cout << "   g4name: " << name << G4endl;
    G4cout << "   PDG:    " << pdgEncoding << G4endl;
    G4cout << "   pdgQ:   " << pdgQ << G4endl;
    G4cout << "   type:   " << rootType << G4endl;
  }
  
  // Add particle to TDatabasePDG
  return TDatabasePDG::Instance()
  ->AddParticle(name, g4Name,
                particleDefinition->GetPDGMass()/CLHEP::GeV,
                particleDefinition->GetPDGStable(),
                particleDefinition->GetPDGWidth()/CLHEP::GeV,
                pdgQ*3, rootType, pdgEncoding);
}

