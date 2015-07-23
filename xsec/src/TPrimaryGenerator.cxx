
#include "TPrimaryGenerator.h"

#include "TMath.h"
#ifdef USE_VECGEOM_NAVIGATOR
#include "volumes/Particle.h"
#else
#include "TDatabasePDG.h"
#endif
#include "GeantTrack.h"

ClassImp(TPrimaryGenerator)

    //______________________________________________________________________________
    TPrimaryGenerator::TPrimaryGenerator()
    : fPDG(11),                        // PDG code of the primary: 11 -> e-
      fPartEkin(0.03),                 // kinetic energy of the primary [GeV] : 30 MeV
      fXPos(0.),                       // (x,y,z) position of the primary particles: (0,0,0)
      fYPos(0.), fZPos(0.), fXDir(0.), // direction vector of the primary particles: (0,0,1)
      fYDir(0.), fZDir(1.), fGVPartIndex(-1), fPartPDG(0), fMass(0), fCharge(0), fPTotal(0), fETotal(0) {
  // init all remaining members
  InitPrimaryGenerator();
}

//______________________________________________________________________________
TPrimaryGenerator::TPrimaryGenerator(Int_t partpdg, Double_t partekin, Double_t xpos, Double_t ypos, Double_t zpos,
                                     Double_t xdir, Double_t ydir, Double_t zdir)
    : fPDG(partpdg),                         // PDG code of the primary particle
      fPartEkin(partekin),                   // kinetic energy of the primary [GeV]
      fXPos(xpos),                           // (x,y,z) position of the primary particles
      fYPos(ypos), fZPos(zpos), fXDir(xdir), // direction vector of the primary particles
      fYDir(ydir), fZDir(zdir), fGVPartIndex(-1), fPartPDG(0), fMass(0), fCharge(0), fPTotal(0), fETotal(0) {
  // ensure normality of the direction vector
  Double_t norm = TMath::Sqrt(fXDir * fXDir + fYDir * fYDir + fZDir * fZDir);
  fXDir /= norm;
  fYDir /= norm;
  fZDir /= norm;
  // init all remaining members
  InitPrimaryGenerator();
}

//______________________________________________________________________________
TPrimaryGenerator::~TPrimaryGenerator() {}

//______________________________________________________________________________
void TPrimaryGenerator::SetParticleByPDGCode(Int_t pdgcode) {
  fPDG = pdgcode;
  // update all remaining members
  InitPrimaryGenerator();
}

//______________________________________________________________________________
void TPrimaryGenerator::InitPrimaryGenerator() {
#ifdef USE_VECGEOM_NAVIGATOR
  Particle::CreateParticles();
#endif
  // set GV particle index
  fGVPartIndex = TPartIndex::I()->PartIndex(fPDG);
// set TDatabasePDG ptr
#ifdef USE_VECGEOM_NAVIGATOR
  fPartPDG = const_cast<Particle *>(&Particle::GetParticle(fPDG));
#else
  fPartPDG = TDatabasePDG::Instance()->GetParticle(fPDG);
#endif
  // set rest mass [GeV]
  fMass = fPartPDG->Mass();
  // set charge
  fCharge = fPartPDG->Charge() / 3.;
  // set total energy [GeV]
  fETotal = fPartEkin + fMass;
  // set total momentum [GeV]
  fPTotal = TMath::Sqrt((fETotal - fMass) * (fETotal + fMass));
}

//______________________________________________________________________________
void TPrimaryGenerator::SetParticleXYZDir(Double_t xdir, Double_t ydir, Double_t zdir) {
  // ensure normality of the direction vector
  Double_t norm = TMath::Sqrt(xdir * xdir + ydir * ydir + zdir * zdir);
  fXDir = xdir / norm;
  fYDir = ydir / norm;
  fZDir = zdir / norm;
}

//______________________________________________________________________________
void TPrimaryGenerator::InitPrimaryTrack(Geant::GeantTrack &gtrack) {
  gtrack.SetPDG(fPDG);
  gtrack.SetGVcode(fGVPartIndex);
  gtrack.fXpos = fXPos;
  gtrack.fYpos = fYPos;
  gtrack.fZpos = fZPos;
  gtrack.fXdir = fXDir;
  gtrack.fYdir = fYDir;
  gtrack.fZdir = fZDir;

  gtrack.SetCharge(fCharge);
  gtrack.SetMass(fMass);
  gtrack.fE = fETotal;
  gtrack.SetP(fPTotal);
}
