#include "Geant/TPrimaryGenerator.h"

#include "Geant/Error.h"

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
TPrimaryGenerator::TPrimaryGenerator(int partpdg, double partekin, double xpos, double ypos, double zpos, double xdir,
                                     double ydir, double zdir)
    : fPDG(partpdg),                         // PDG code of the primary particle
      fPartEkin(partekin),                   // kinetic energy of the primary [GeV]
      fXPos(xpos),                           // (x,y,z) position of the primary particles
      fYPos(ypos), fZPos(zpos), fXDir(xdir), // direction vector of the primary particles
      fYDir(ydir), fZDir(zdir), fGVPartIndex(-1), fPartPDG(0), fMass(0), fCharge(0), fPTotal(0), fETotal(0) {
  // ensure normality of the direction vector
  double norm = sqrt(fXDir * fXDir + fYDir * fYDir + fZDir * fZDir);
  fXDir /= norm;
  fYDir /= norm;
  fZDir /= norm;
  // init all remaining members
  InitPrimaryGenerator();
}

//______________________________________________________________________________
TPrimaryGenerator::~TPrimaryGenerator() {}

//______________________________________________________________________________
void TPrimaryGenerator::SetParticleByPDGCode(int pdgcode) {
  fPDG = pdgcode;
  // update all remaining members
  InitPrimaryGenerator();
}

//______________________________________________________________________________
void TPrimaryGenerator::InitPrimaryGenerator() {
  Particle_t::CreateParticles();
  // set GV particle index
  fGVPartIndex = TPartIndex::I()->PartIndex(fPDG);
// set TDatabasePDG ptr
  fPartPDG = const_cast<Particle_t *>(&Particle_t::GetParticle(fPDG));
  // set rest mass [GeV]
  fMass = fPartPDG->Mass();
  // set charge
  fCharge = fPartPDG->Charge();
  if ((int)fCharge != fCharge)
     geant::Error("TPrimaryGenerator::InitPrimaryGenerator()","Unsupported charge: %f\n",fCharge);


  // set total energy [GeV]
  fETotal = fPartEkin + fMass;
  // set total momentum [GeV]
  fPTotal = sqrt((fETotal - fMass) * (fETotal + fMass));
}

//______________________________________________________________________________
void TPrimaryGenerator::SetParticleXYZDir(double xdir, double ydir, double zdir) {
  // ensure normality of the direction vector
  double norm = sqrt(xdir * xdir + ydir * ydir + zdir * zdir);
  fXDir = xdir / norm;
  fYDir = ydir / norm;
  fZDir = zdir / norm;
}

//______________________________________________________________________________
void TPrimaryGenerator::InitPrimaryTrack(geant::Track &gtrack) {
  gtrack.SetPDG(fPDG);
  gtrack.SetGVcode(fGVPartIndex);
  gtrack.SetPosition(fXPos, fYPos, fZPos);
  gtrack.SetDirection(fXDir, fYDir, fZDir);
  gtrack.SetCharge(fCharge);
  gtrack.SetMass(fMass);
  gtrack.SetE(fETotal);
  gtrack.SetP(fPTotal);
}
