
#include "LHCbParticleGun.h"

#include "Particle.h"
#include "PhysicalConstants.h"

// from geantV
#include "GeantTrack.h"
#include "GeantTaskData.h"

#include <cmath>

namespace lhcbapp {


const int      LHCbParticleGun::gInfNumPrimaryPerEvt =  1;
const int      LHCbParticleGun::gSupNumPrimaryPerEvt = 10;

const double   LHCbParticleGun::gInfBeamEnergy       =   1.*geant::GeV;
const double   LHCbParticleGun::gSupBeamEnergy       = 100.*geant::GeV;

// these static variables stores the gun configuration just for the Print()
int          LHCbParticleGun::gNumPrimaryPerEvt(-1);
double       LHCbParticleGun::gPrimaryEnergy(-1.);
std::string  LHCbParticleGun::gPrimaryType("");
double       LHCbParticleGun::gPrimaryDir[3] = {0.,0.,0.};

// These are the particle types that can be used as primary beam particle, on a event-by-event based.
const int LHCbParticleGun::gNumberCandidateParticles = 3;
const std::string LHCbParticleGun::gNameParticlesVector[LHCbParticleGun::gNumberCandidateParticles] = {
  "e-",
  "e+",
  "gamma"
};
const std::map<std::string,int> LHCbParticleGun::gPrimaryNameToIndexMap = {
  {"e-",   0},
  {"e+",   1},
  {"gamma",2}
};


LHCbParticleGun::LHCbParticleGun()  {
  fIsUserNumPrimaryPerEvt     = false;
  fIsUserPrimaryType          = false;
  fIsUserPrimaryDir           = false;
  fIsUserPrimaryEnergy        = false;
  fNumPrimaryPerEvt           =  1;
  fPrimaryParticleName        = "e-";
  fPrimaryParticleEnergy      = 10.*geant::GeV;
  // position will be fixed to (0.,0.,0.)
  fXPos                       = 0.;
  fYPos                       = 0.;
  fZPos                       = 0.;
  fXDir                       = 0.;
  fYDir                       = 1.;
  fZDir                       = 0.;
}


LHCbParticleGun::~LHCbParticleGun() {}


// set number of primary vertex and the fixed vertex position and give it back to the EventServer
Geant::GeantEventInfo LHCbParticleGun::NextEvent(Geant::GeantTaskData* td) {
  Geant::GeantEventInfo current;
  // check if user set fix number of primaries per event and generate randomly if not
  int numPrimPerEvt = fNumPrimaryPerEvt;
  if (!fIsUserNumPrimaryPerEvt) {
    int idum        = static_cast<int>(td->fRndm->uniform()*(gSupNumPrimaryPerEvt-gInfNumPrimaryPerEvt+1));
    numPrimPerEvt   = std::max(gInfNumPrimaryPerEvt + idum,1);
  }
  current.ntracks = numPrimPerEvt;
  // primary vertex is fixed: always at (0,0,0)
  current.xvert   = fXPos;
  current.yvert   = fYPos;
  current.zvert   = fZPos;
  return current;
}


void LHCbParticleGun::GetTrack(int /*n*/, Geant::GeantTrack &gtrack, Geant::GeantTaskData* td) {
  // Select randomly the primary particle if it was not set by the user
  std::string &pParticleName = fPrimaryParticleName;
  if (!fIsUserPrimaryType) {
    int caseBeamParticle = static_cast<int>(td->fRndm->uniform()*LHCbParticleGun::gNumberCandidateParticles);
    pParticleName        = gNameParticlesVector[caseBeamParticle];
  }
  // Select randomly the beam energy if it was not set by the user.
  double pParticleEkin = fPrimaryParticleEnergy;
  if (!fIsUserPrimaryEnergy) {
    pParticleEkin = gInfBeamEnergy+td->fRndm->uniform()*(gSupBeamEnergy-gInfBeamEnergy);
  }
  // Select random direction if it was not set by the user
  double pParticleDir[3] = {fXDir, fYDir, fZDir};
  if (!fIsUserPrimaryDir) {
    double *rndArray = td->fDblArray;
    td->fRndm->uniform_array(2, rndArray);
    const double cost = 1.-2.*rndArray[0];
    double sint2      = 1.-cost*cost;
    if (sint2<0.) {
      sint2 = 0.;
    }
    const double sint = std::sqrt(sint2);
    const double phi  = rndArray[1]*geant::kTwoPi;
    pParticleDir[0]   = sint*std::cos(phi);
    pParticleDir[1]   = sint*std::sin(phi);
    pParticleDir[2]   = cost;
  }
  // get the particle
  const geantphysics::Particle* pParticle = geantphysics::Particle::GetParticleByName(pParticleName);
  if (!pParticle) {
    std::cerr<< "   \n *** ERROR::LHCbParticleGun::GetTrack()                      \n"
             << "          unknown particle name = " << pParticle << " \n"
             << std::endl;
    exit(-1);
  }
  const int    pPartPDG   = pParticle->GetPDGCode();
  const int    pPartGVc   = pParticle->GetInternalCode();
  const double pPartMass  = pParticle->GetPDGMass();
  const double pPartChrg  = pParticle->GetPDGCharge();
  const double pPartEt    = pParticleEkin + pPartMass;
  const double pPartPt    = std::sqrt((pPartEt-pPartMass)*(pPartEt+pPartMass));
  //
  // set the primary track
  gtrack.SetPDG(pPartPDG);
  gtrack.SetGVcode(pPartGVc);
  gtrack.SetPosition(fXPos, fYPos, fZPos);
  gtrack.SetDirection(pParticleDir[0],pParticleDir[1],pParticleDir[2]);
  gtrack.SetCharge(pPartChrg);
  gtrack.SetMass(pPartMass);
  gtrack.SetE(pPartEt);
  gtrack.SetP(pPartPt);
}


void  LHCbParticleGun::SetNumPrimaryPerEvt(const int pperevt) {
  fNumPrimaryPerEvt         = pperevt;
  gNumPrimaryPerEvt         = fNumPrimaryPerEvt;
  fIsUserNumPrimaryPerEvt   = true;
}

void  LHCbParticleGun::SetPrimaryEnergy(const double ekin) {
  fPrimaryParticleEnergy    = ekin;
  gPrimaryEnergy            = fPrimaryParticleEnergy;
  fIsUserPrimaryEnergy      = true;
}

void  LHCbParticleGun::SetPrimaryName(const std::string& pname) {
  fPrimaryParticleName      = pname;
  gPrimaryType              = fPrimaryParticleName;
  fIsUserPrimaryType        = true;
}

void  LHCbParticleGun::SetPrimaryDirection(const double pdir[3]) {
  double norm = 1./std::sqrt(pdir[0]*pdir[0] + pdir[1]*pdir[1] + pdir[2]*pdir[2]);
  fXDir             = pdir[0]*norm;
  fYDir             = pdir[1]*norm;
  fZDir             = pdir[2]*norm;
  gPrimaryDir[0]    = fXDir;
  gPrimaryDir[1]    = fYDir;
  gPrimaryDir[2]    = fZDir;
  fIsUserPrimaryDir = true;
}

int LHCbParticleGun::GetPrimaryTypeIndex(const std::string& pname) {
  int indx = gPrimaryNameToIndexMap.find(pname)->second;
  return indx;
}


// will give proper results only at the end of the run
void  LHCbParticleGun::Print() {
  using geant::GeV;
  std::string str = "  Primaries were generated : \n";
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
  if (!(gPrimaryDir[0] && gPrimaryDir[1] && gPrimaryDir[2])) {
    str += "     Primary direction     : isotropic for each primary \n";
  } else {
    std::string sdir= "[";
    sdir += std::to_string(gPrimaryDir[0]) + ", "
          + std::to_string(gPrimaryDir[1]) + ", "
          + std::to_string(gPrimaryDir[2]) + "]\n";
    str += "     Primary direction     : " + sdir;
  }
  if (gPrimaryType=="") {
    str += "     Primary type          : randomly selected for each primary from \n";
    for (int i=0; i<gNumberCandidateParticles; i++) {
      str += "       type index: " + std::to_string(i) + ",  name: " + gNameParticlesVector[i] + "\n";
    }
  } else {
    str += "     Primary type          : " + gPrimaryType + "\n";
  }
  std::cout << " \n  ======= Info On Run Conditions ======================================================== \n"
            << str
            << "  --------------------------------------------------------------------------------------- "
            << std::endl;
}



}  // namespace lhcbapp
