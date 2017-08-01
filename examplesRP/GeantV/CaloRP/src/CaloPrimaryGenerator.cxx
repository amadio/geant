
#include "CaloPrimaryGenerator.h"
#include "CaloDetectorConstruction.h"
#include "Particle.h"

// from geantV
#include "GeantTrack.h"

#include <cmath>

namespace userapplication {

CaloPrimaryGenerator::CaloPrimaryGenerator(CaloDetectorConstruction *det) : fDetector(det) {
  fPrimaryParticleName = "e-";
  fPrimaryPerEvent     = 1;
  fParticle            = nullptr;
  //
  fPDG                 = 0;
  fGVPartIndex         = 0;
  //
  fPrimaryEkin         = 15.7*geant::MeV;
  //
  fXPos                = 0.;
  fYPos                = 0.;
  fZPos                = 0.;
  //
  fXDir                = 1.;
  fYDir                = 0.;
  fZDir                = 0.;
  //
  fMass                = 0.;
  fCharge              = 0.;
  fETotal              = 0.;
  fPTotal              = 0.;
}

CaloPrimaryGenerator::~CaloPrimaryGenerator() {}

void CaloPrimaryGenerator::InitPrimaryGenerator() {
  fParticle            = geantphysics::Particle::GetParticleByName(fPrimaryParticleName);
  if (!fParticle) {
    std::cerr<< "   \n *** ERROR::UserPrimaryGenerator::InitPrimaryGenerator()    \n"
             << "          unknown particle name = " << fPrimaryParticleName << " \n"
             << std::endl;
    exit(-1);
  }
  fPDG                 = fParticle->GetPDGCode();
  fGVPartIndex         = fParticle->GetInternalCode();
  fMass                = fParticle->GetPDGMass();
  fCharge              = fParticle->GetPDGCharge();
  fETotal              = fPrimaryEkin + fMass;
  fPTotal              = std::sqrt((fETotal - fMass) * (fETotal + fMass));
  //
  fXPos                = -0.25*(fDetector->GetWorldX()+fDetector->GetDetectorX());
  fYPos                = 0.;
  fZPos                = 0.;
  //
  fXDir                = 1.;
  fYDir                = 0.;
  fZDir                = 0.;
}


Geant::GeantEventInfo CaloPrimaryGenerator::NextEvent() {
  Geant::GeantEventInfo current;
  current.ntracks = fPrimaryPerEvent;
  current.xvert   = fXPos;
  current.yvert   = fYPos;
  current.zvert   = fZPos;
  return current;
}


void CaloPrimaryGenerator::GetTrack(int /*n*/, Geant::GeantTrack &gtrack) {
  gtrack.SetPDG(fPDG);
  gtrack.SetGVcode(fGVPartIndex);
  gtrack.fXpos = fXPos;
  gtrack.fYpos = fYPos;
  gtrack.fZpos = fZPos;
  gtrack.fXdir = fXDir;
  gtrack.fYdir = fYDir;
  gtrack.fZdir = fZDir;
  //
  gtrack.SetCharge(fCharge);
  gtrack.SetMass(fMass);
  gtrack.fE = fETotal;
  gtrack.SetP(fPTotal);
}


}  // namespace userapplication
