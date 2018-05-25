
#include "MyPrimaryGeneratorAction.hh"

#include "MyPrimaryGeneratorMessenger.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
// just for the final printout
#include "MyDetectorConstruction.hh"

const G4int      MyPrimaryGeneratorAction::gInfNumPrimaryPerEvt =  1;
const G4int      MyPrimaryGeneratorAction::gSupNumPrimaryPerEvt = 10;

const G4double   MyPrimaryGeneratorAction::gInfBeamEnergy       =   1.*GeV;
const G4double   MyPrimaryGeneratorAction::gSupBeamEnergy       = 100.*GeV;

// these static variables stores the gun configuration just for the Print()
G4int            MyPrimaryGeneratorAction::gNumPrimaryPerEvt(-1);
G4double         MyPrimaryGeneratorAction::gPrimaryEnergy(-1.);
std::string      MyPrimaryGeneratorAction::gPrimaryType("");
G4ThreeVector    MyPrimaryGeneratorAction::gPrimaryDir(0.,0.,0.);


// These are the particle types that can be used as primary beam particle, on a event-by-event based.
const G4int MyPrimaryGeneratorAction::gNumberCandidateParticles = 3;
const G4String MyPrimaryGeneratorAction::gNameParticlesVector[MyPrimaryGeneratorAction::gNumberCandidateParticles] = {
  "e-",
  "e+",
  "gamma"
};
const std::map<G4String,G4int> MyPrimaryGeneratorAction::gPrimaryNameToIndexMap = {
  {"e-",0},
  {"e+",1},
  {"gamma",2}
};


MyPrimaryGeneratorAction::MyPrimaryGeneratorAction() {
  fIsUserNumPrimaryPerEvt     = false;
  fIsUserPrimaryType          = false;
  fIsUserPrimaryDir           = false;
  fIsUserPrimaryEnergy        = false;
  fNumPrimaryPerEvt           = 1;
  fPrimaryParticleName        = "e-";
  fParticleGun                = new G4ParticleGun(1);
  fParticleTable              = G4ParticleTable::GetParticleTable();
  fPrimaryParticleEnergy      =  10.*GeV;
  fPrimaryParticlePosition    = G4ThreeVector(0.0,0.0,0.0);
  //
  fGunMessenger  = new MyPrimaryGeneratorMessenger(this);
}


MyPrimaryGeneratorAction::~MyPrimaryGeneratorAction() {
  delete fParticleGun;
  delete fGunMessenger;
}


void MyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  if (!fIsUserNumPrimaryPerEvt) {
    G4int idum        = static_cast<G4int>(G4UniformRand()*(gSupNumPrimaryPerEvt-gInfNumPrimaryPerEvt+1));
    fNumPrimaryPerEvt = std::max(gInfNumPrimaryPerEvt + idum,1);
  }
  for (G4int i=0; i<fNumPrimaryPerEvt; ++i) {
    // Select randomly the primary particle if it was not set by the user
    if (!fIsUserPrimaryType) {
      G4int caseBeamParticle = static_cast<G4int>(G4UniformRand()*MyPrimaryGeneratorAction::gNumberCandidateParticles);
      fPrimaryParticleName = gNameParticlesVector[caseBeamParticle];
    }
    // Select randomly the beam energy if it was not set by the user.
    if (!fIsUserPrimaryEnergy) {
      fPrimaryParticleEnergy = gInfBeamEnergy+G4UniformRand()*(gSupBeamEnergy-gInfBeamEnergy);
    }
    // Select random direction if it was not set by the user
    if (!fIsUserPrimaryDir) {
      fPrimaryParticleDirection = G4RandomDirection();
    }
    // Beam position: always the origin i.e. [0,0,0].
    //
    // Set the particle gun
    G4ParticleDefinition* pDef = fParticleTable->FindParticle(fPrimaryParticleName);
    if (!pDef) {
      G4cerr<< " *** ERROR in MyPrimaryGeneratorAction::GeneratePrimaries() " << G4endl
            << "       UNKNOWN PRIMARY PARTICLE WITH NAME = "
            << fPrimaryParticleName << G4endl
            << G4endl;
      exit(-1);
    }
    fParticleGun->SetParticleDefinition       (pDef                     );
    fParticleGun->SetParticleEnergy           (fPrimaryParticleEnergy   );
    fParticleGun->SetParticlePosition         (fPrimaryParticlePosition );
    fParticleGun->SetParticleMomentumDirection(fPrimaryParticleDirection);
    //
    fParticleGun->GeneratePrimaryVertex(anEvent);
 }
}

void  MyPrimaryGeneratorAction::SetNumPrimaryPerEvt(G4int pperevt) {
  fNumPrimaryPerEvt         = pperevt;
  gNumPrimaryPerEvt         = fNumPrimaryPerEvt;
  fIsUserNumPrimaryPerEvt   = true;
}

void  MyPrimaryGeneratorAction::SetPrimaryEnergy(G4double ekin) {
  fPrimaryParticleEnergy    = ekin;
  gPrimaryEnergy            = fPrimaryParticleEnergy;
  fIsUserPrimaryEnergy      = true;
}

void  MyPrimaryGeneratorAction::SetPrimaryName(const G4String& pname) {

  // Check the name provided first 
  G4ParticleDefinition* pDef = fParticleTable->FindParticle(pname);
  if (!pDef) {   
     G4cerr<< " *** ERROR> Proposed particle name " << pname << "  is not known.  It is IGNORED."
           << G4endl << G4endl;
     // exit(1);  // Good for batch mode - bad for interactive ... 
  } else { 
     fPrimaryParticleName      = pname;
     gPrimaryType              = fPrimaryParticleName;
     fIsUserPrimaryType        = true;
  }
}

void  MyPrimaryGeneratorAction::SetPrimaryDirection(const G4ThreeVector &pdir) {
  fPrimaryParticleDirection = pdir;
  fPrimaryParticleDirection.setMag(1.);
  gPrimaryDir = fPrimaryParticleDirection;
  fIsUserPrimaryDir         = true;
}

G4int MyPrimaryGeneratorAction::GetPrimaryTypeIndex(const G4String& pname) {
  G4int indx = gPrimaryNameToIndexMap.find(pname)->second;
  return indx;
}

// will give proper results only at the end of the run
void  MyPrimaryGeneratorAction::Print() {
  G4String str = "";
  G4double magFValue = MyDetectorConstruction::GetFieldValue();
  if( MyDetectorConstruction::IsFieldUniform() ){
     if (magFValue>0.) {
        str += "  Magnetic field           : constant " + std::to_string(magFValue/tesla) + " [ Tesla] field \n";
     } else {
        str += "  Magnetic field           : without magnetic field \n";
     }
  } else {
     str += "  Magnetic field           :  simplified RZ field from map ; central value = " +
        std::to_string(magFValue/tesla) + " [ Tesla] field \n";
  }
  str   += "  Primaries were generated : \n";
  if (gNumPrimaryPerEvt<0) {
    str += "     Primaries per event   : random for ecah event on ["
           + std::to_string(gInfNumPrimaryPerEvt) + ", " + std::to_string(gSupNumPrimaryPerEvt) + "]\n";
  } else {
    str += "     Primaries per event   : " + std::to_string(gNumPrimaryPerEvt) + "\n";
  }
  if (gPrimaryEnergy<0.) {
    str += "     Primary energy        : random for each primary on ["
           + std::to_string(gInfBeamEnergy/GeV) + " GeV, " +std::to_string(gSupBeamEnergy/GeV) + " GeV] \n";
  } else {
    str += "     Primary energy        : " + std::to_string(gPrimaryEnergy/GeV) + " [GeV] \n";
  }
  if (gPrimaryDir.mag()==0.) {
    str += "     Primary direction     : isotropic for each primary \n";
  } else {
    G4String sdir= "[";
    sdir += std::to_string(gPrimaryDir.x()) + ", "
          + std::to_string(gPrimaryDir.y()) + ", "
          + std::to_string(gPrimaryDir.x()) + "]\n";
    str += "     Primary direction     : " + sdir;
  }
  if (gPrimaryType=="") {
    str += "     Primary type          : randomly selected for each primary from \n";
    for (G4int i=0; i<gNumberCandidateParticles; i++) {
      str += "       type index: " + std::to_string(i) + ",  name: " + gNameParticlesVector[i] + "\n";
    }
  } else {
    str += "     Primary type          : " + gPrimaryType + "\n";
  }
  G4cout << " \n  ======= Info On Run Conditions ======================================================== \n"
         << str
         << "  --------------------------------------------------------------------------------------- "
         << G4endl;
}
