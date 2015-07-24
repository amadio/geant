#include "GunGenerator.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "management/GeoManager.h"
typedef vecgeom::GeoManager TGeoManager;
#define gGeoManager &vecgeom::GeoManager::Instance()
#else
#include "TGeoManager.h"
#endif

#include "TRandom.h"
#ifdef USE_VECGEOM_NAVIGATOR
#include "volumes/Particle.h"
using vecgeom::Particle;
#else
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#endif
#include "GeantTrack.h"

ClassImp(GunGenerator)

    //______________________________________________________________________________
    GunGenerator::GunGenerator()
    : average(0), fPDG(11),            // PDG code of the primary: 11 -> e-
      fPartEkin(0.03),                 // kinetic energy of the primary [GeV] : 30 MeV
      fXPos(0.),                       // (x,y,z) position of the primary particles: (0,0,0)
      fYPos(0.), fZPos(0.), fXDir(0.), // direction vector of the primary particles: (0,0,1)
      fYDir(0.), fZDir(1.), fGVPartIndex(-1), fPartPDG(0), fMass(0), fCharge(0), fPTotal(0), fETotal(0),
      numberoftracks(0), rndgen(0) {}

GunGenerator::GunGenerator(int aver, int partpdg, double partekin, double xpos, double ypos, double zpos, double xdir,
                           double ydir, double zdir)
    : average(aver), fPDG(partpdg),          // PDG code of the primary particle
      fPartEkin(partekin),                   // kinetic energy of the primary [GeV]
      fXPos(xpos),                           // (x,y,z) position of the primary particles
      fYPos(ypos), fZPos(zpos), fXDir(xdir), // direction vector of the primary particles
      fYDir(ydir), fZDir(zdir), fGVPartIndex(-1), fPartPDG(0), fMass(0), fCharge(0), fPTotal(0), fETotal(0),
      numberoftracks(0), rndgen(0) {
  // ensure normality of the direction vector
  double norm = sqrt(fXDir * fXDir + fYDir * fYDir + fZDir * fZDir);
  fXDir /= norm;
  fYDir /= norm;
  fZDir /= norm;

  rndgen = new TRandom();
}

//______________________________________________________________________________
GunGenerator::~GunGenerator() { delete rndgen; }

//______________________________________________________________________________
void GunGenerator::InitPrimaryGenerator() {
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
  fPTotal = sqrt((fETotal - fMass) * (fETotal + fMass));
}

//______________________________________________________________________________
int GunGenerator::NextEvent() {
  //
  if (average == 1)
    numberoftracks = 1;
  else
    numberoftracks = rndgen->Poisson(average);
  // here are generate an event with ntracks

  for (int nn = 1; nn <= numberoftracks; nn++) {
    // here I would normally push back the generated particles to some vector
    // no need to do it in this specific case, because all the particles are the same
  }

  return numberoftracks;
}

//______________________________________________________________________________
void GunGenerator::GetTrack(int /*n*/, Geant::GeantTrack &gtrack) {
  // here I get the n-th generated track and copy it to gtrack
  // they are all the same here, so no dependence on n

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
