#include "Geant/ParticleGun.h"
#include "Geant/Typedefs.h"
#include "Geant/Error.h"
#include "Geant/Track.h"

#include "Geant/Particle.h"

namespace geantphysics {

//______________________________________________________________________________
ParticleGun::ParticleGun()
    : fAverage(0), fPartEkin(0.03),    // kinetic energy of the primary [GeV] : 30 MeV
      fXPos(0.),                       // (x,y,z) position of the primary particles: (0,0,0)
      fYPos(0.), fZPos(0.), fXDir(0.), // direction vector of the primary particles: (0,0,1)
      fYDir(0.), fZDir(1.), fGVPartIndex(-1), fMass(0), fCharge(0), fPTotal(0), fETotal(0), fRndgen(0)
{
  // Dummy constructor
}

//______________________________________________________________________________
ParticleGun::ParticleGun(int aver, int gvcode, double partekin, double xpos, double ypos, double zpos, double xdir,
                         double ydir, double zdir)
    : fAverage(aver), fPartEkin(partekin), fXPos(xpos), fYPos(ypos), fZPos(zpos), fXDir(xdir), fYDir(ydir), fZDir(zdir),
      fGVPartIndex(gvcode), fMass(0), fCharge(0), fPTotal(0), fETotal(0), fRndgen(0)
{
  // Constructor
  // ensure normality of the direction vector
  double norm = sqrt(fXDir * fXDir + fYDir * fYDir + fZDir * fZDir);
  fXDir /= norm;
  fYDir /= norm;
  fZDir /= norm;
  //#ifdef USE_ROOT
  //  fRndgen = new TRandom();
  //#else
  fRndgen = &vecgeom::RNG::Instance();
  //#endif
}

//______________________________________________________________________________
ParticleGun::~ParticleGun()
{
  // Destructor
  delete fRndgen;
}

//______________________________________________________________________________
void ParticleGun::InitPrimaryGenerator()
{

  const Particle *part = Particle::GetParticleByInternalCode(fGVPartIndex);
  if (part) {
    fMass = part->GetPDGMass();
    // set charge
    fCharge = part->GetPDGCharge();
    // set total energy [GeV]
    fETotal = fPartEkin + fMass;
    // set total momentum [GeV]
    fPTotal = sqrt((fETotal - fMass) * (fETotal + fMass));
  }
}

//______________________________________________________________________________
EventInfo ParticleGun::NextEvent(geant::TaskData * /*td*/)
{
  //
  int ntracks = 1;
  ntracks     = fAverage;
  //  if (fAverage>1)
  //    ntracks = fRndgen->Poisson(fAverage);
  // here are generate an event with ntracks

  // for (int nn = 1; nn <= ntracks; nn++) {
  // here I would normally push back the generated particles to some vector
  // no need to do it in this specific case, because all the particles are the same
  //}

  EventInfo current;
  current.ntracks = ntracks;
  current.xvert   = fXPos;
  current.yvert   = fYPos;
  current.zvert   = fZPos;
  return current;
}

//______________________________________________________________________________
void ParticleGun::GetTrack(int /*n*/, geant::Track &gtrack, geant::TaskData * /*td*/)
{
  // here I get the n-th generated track and copy it to gtrack
  // they are all the same here, so no dependence on n
  gtrack.SetGVcode(fGVPartIndex);
  gtrack.SetPosition(fXPos, fYPos, fZPos);
  gtrack.SetDirection(fXDir, fYDir, fZDir);

  gtrack.SetCharge(fCharge);
  gtrack.SetMass(fMass);
  gtrack.SetE(fETotal);
  gtrack.SetP(fPTotal);
}
} // namespace geantphysics
