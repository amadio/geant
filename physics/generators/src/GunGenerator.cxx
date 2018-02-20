#include "GunGenerator.h"
#include "GeantTaskData.h"
#include "Geant/Typedefs.h"
#include "Geant/Error.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GunGenerator::GunGenerator()
    : fAverage(0), fPDG(11),           // PDG code of the primary: 11 -> e-
      fPartEkin(0.03),                 // kinetic energy of the primary [GeV] : 30 MeV
      fXPos(0.),                       // (x,y,z) position of the primary particles: (0,0,0)
      fYPos(0.), fZPos(0.), fXDir(0.), // direction vector of the primary particles: (0,0,1)
      fYDir(0.), fZDir(1.), fGVPartIndex(-1), fPartPDG(0), fMass(0), fCharge(0), fPTotal(0), fETotal(0),
      fRndgen(0) {
  // Dummy constructor
}

//______________________________________________________________________________
GunGenerator::GunGenerator(int aver, int partpdg, double partekin, double xpos, double ypos, double zpos, double xdir,
                           double ydir, double zdir)
    : fAverage(aver), fPDG(partpdg),
      fPartEkin(partekin),
      fXPos(xpos), fYPos(ypos), fZPos(zpos),
      fXDir(xdir), fYDir(ydir), fZDir(zdir),
      fGVPartIndex(-1), fPartPDG(0), fMass(0), fCharge(0),
      fPTotal(0), fETotal(0),
      fRndgen(0) {
  // Constructor
  // ensure normality of the direction vector
  double norm = sqrt(fXDir * fXDir + fYDir * fYDir + fZDir * fZDir);
  fXDir /= norm;
  fYDir /= norm;
  fZDir /= norm;
#ifdef USE_ROOT
  fRndgen = new TRandom();
#else
  fRndgen = &RNG::Instance();
#endif
}

//______________________________________________________________________________
GunGenerator::~GunGenerator() {
  // Destructor
  delete fRndgen;
}

//______________________________________________________________________________
void GunGenerator::InitPrimaryGenerator() {
#ifdef USE_VECGEOM_NAVIGATOR
  Particle_t::CreateParticles();
#endif
  // set GV particle index
  fGVPartIndex = TPartIndex::I()->PartIndex(fPDG);
// set TDatabasePDG ptr
#ifdef USE_VECGEOM_NAVIGATOR
  fPartPDG = const_cast<Particle_t *>(&Particle_t::GetParticle(fPDG));
#else
  fPartPDG = TDatabasePDG::Instance()->GetParticle(fPDG);
#endif
  // set rest mass [GeV]
  fMass = fPartPDG->Mass();
  // set charge
  fCharge = fPartPDG->Charge();
#ifndef USE_VECGEOM_NAVIGATOR
  fCharge /= 3.;
#endif
  if ((int)fCharge != fCharge)
     Geant::Error("TPrimaryGenerator::InitPrimaryGenerator()","Unsupported charge: %f\n",fCharge);

  // set total energy [GeV]
  fETotal = fPartEkin + fMass;
  // set total momentum [GeV]
  fPTotal = sqrt((fETotal - fMass) * (fETotal + fMass));
}

//______________________________________________________________________________
GeantEventInfo GunGenerator::NextEvent(Geant::GeantTaskData* /*td*/) {
  //
  int ntracks = 1;
  if (fAverage > 1)
    ntracks = fRndgen->Poisson(fAverage);
  // here are generate an event with ntracks

  for (int nn = 1; nn <= ntracks; nn++) {
    // here I would normally push back the generated particles to some vector
    // no need to do it in this specific case, because all the particles are the same
  }

  GeantEventInfo current;
  current.ntracks = ntracks;
  current.xvert = fXPos;
  current.yvert = fYPos;
  current.zvert = fZPos;
  return current;
}

//______________________________________________________________________________
void GunGenerator::GetTrack(int /*n*/, Geant::GeantTrack &gtrack, Geant::GeantTaskData* /*td*/) {
  // here I get the n-th generated track and copy it to gtrack
  // they are all the same here, so no dependence on n

  gtrack.SetPDG(fPDG);
  gtrack.SetGVcode(fGVPartIndex);
  gtrack.SetPosition(fXPos, fYPos, fZPos);
  gtrack.SetDirection(fXDir, fYDir, fZDir);
  gtrack.SetCharge(fCharge);
  gtrack.SetMass(fMass);
  gtrack.SetE(fETotal);
  gtrack.SetP(fPTotal);
}

} // GEANT_IMPL_NAMESPACE
} // Geant
