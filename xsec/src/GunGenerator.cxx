#include "GunGenerator.h"

#include "TGeoManager.h"

#include "TMath.h"
#include "TRandom.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "GeantTrack.h"
#include "GeantThreadData.h"


ClassImp(GunGenerator)

//______________________________________________________________________________
GunGenerator::GunGenerator():
  fPDG(11),             // PDG code of the primary: 11 -> e-
  fPartEkin(0.03),      // kinetic energy of the primary [GeV] : 30 MeV
  fXPos(0.),            // (x,y,z) position of the primary particles: (0,0,0) 
  fYPos(0.),
  fZPos(0.),
  fXDir(0.),            // direction vector of the primary particles: (0,0,1)
  fYDir(0.),
  fZDir(1.),
  fGVPartIndex(-1),
  fPartPDG(0),
  fMass(0),
  fCharge(0),
  fPTotal(0),
  fETotal(0)
{}

GunGenerator::GunGenerator(Int_t partpdg, Double_t partekin, 
                           Double_t xpos, Double_t ypos, Double_t zpos,
                           Double_t xdir, Double_t ydir, Double_t zdir):
  fPDG(partpdg),                 // PDG code of the primary particle
  fPartEkin(partekin),           // kinetic energy of the primary [GeV]
  fXPos(xpos),                   // (x,y,z) position of the primary particles
  fYPos(ypos),
  fZPos(zpos),
  fXDir(xdir),                   // direction vector of the primary particles
  fYDir(ydir),
  fZDir(zdir),
  fGVPartIndex(-1),
  fPartPDG(0),
  fMass(0),
  fCharge(0),
  fPTotal(0),
  fETotal(0)
{
  // ensure normality of the direction vector
  Double_t norm = TMath::Sqrt(fXDir*fXDir+fYDir*fYDir+fZDir*fZDir);
  fXDir /=norm;
  fYDir /=norm;
  fZDir /=norm;

}

//______________________________________________________________________________
GunGenerator::~GunGenerator()
{
}


//______________________________________________________________________________
void GunGenerator::InitPrimaryGenerator(){
  // set GV particle index
  fGVPartIndex= TPartIndex::I()->PartIndex(fPDG);
  // set TDatabasePDG ptr
  fPartPDG = TDatabasePDG::Instance()->GetParticle(fPDG);
  // set rest mass [GeV]
  fMass    = fPartPDG->Mass();
  //set charge 
  fCharge  = fPartPDG->Charge()/3.;
  // set total energy [GeV]
  fETotal  =  fPartEkin + fMass;
  // set total momentum [GeV]
  fPTotal  =  TMath::Sqrt((fETotal-fMass)*(fETotal+fMass));  
}

//______________________________________________________________________________
Int_t GunGenerator::NextEvent(){
    //
    numberoftracks = TRandom().Poisson(average);
    // here are generate an event with ntracks

    for( Int_t nn=1; nn<=numberoftracks; nn++ ) {
      // I push back the generated particles to some vector

      
    }
    
    return numberoftracks;
}

//______________________________________________________________________________
void GunGenerator::GetTrack(Int_t n, GeantTrack &gtrack){
  // here I get the n-th generated track and copy it to gtrack


     gtrack.SetPDG(fPDG);
     gtrack.SetG5code(fGVPartIndex); 
     gtrack.fXpos = fXPos;
     gtrack.fYpos = fYPos;
     gtrack.fZpos = fZPos;
     gtrack.fXdir = fXDir;
     gtrack.fYdir = fYDir;
     gtrack.fZdir = fZDir;

     gtrack.SetCharge(fCharge);
     gtrack.SetMass(fMass);
     gtrack.fE    = fETotal;
     gtrack.SetP(fPTotal);
}
   


