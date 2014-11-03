
#ifndef GunGenerator_h
#define GunGenerator_h

#include "TPartIndex.h"
#include "PrimaryGenerator.h"
#include "TRandom.h"

class TParticlePDG;
class GeantTrack;

class GunGenerator: public PrimaryGenerator{
 private:
  Int_t average;

  Int_t     fPDG;             // PDG code of parimary particles
  Double_t  fPartEkin;        // kinetic energy of the primary [GeV]
  Double_t  fXPos;            // (x,y,z) position of the primary particles 
  Double_t  fYPos;
  Double_t  fZPos;
  Double_t  fXDir;            // direction vector of the primary particles
  Double_t  fYDir;
  Double_t  fZDir;
  // additional members
  Int_t     fGVPartIndex;     // GV particle index of the primary
  TParticlePDG  *fPartPDG;
  Double_t  fMass;            // rest mass of the primary [GeV]
  Double_t  fCharge;          // charge of the primary 
  Double_t  fPTotal;          // total momentum of the primary [GeV]
  Double_t  fETotal;          // total energy of the primary [GeV]

  Int_t numberoftracks;
    
  TRandom* rndgen;
 public:
    GunGenerator();
    GunGenerator(Int_t aver, Int_t partpdg, Double_t partekin, 
                 Double_t xpos, Double_t ypos, Double_t zpos,
                 Double_t xdir, Double_t ydir, Double_t zdir);

    ~GunGenerator();

  // set one GeantTrack primary track properties
    virtual void InitPrimaryGenerator();
    virtual Int_t NextEvent();
    virtual void GetTrack(Int_t n, GeantTrack &gtrack);

 private:

   GunGenerator(const GunGenerator &);//no imp.	
   GunGenerator& operator=(const GunGenerator &);//no imp.

   ClassDef(GunGenerator,1)
};

#endif
