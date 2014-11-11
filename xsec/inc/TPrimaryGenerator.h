#ifndef ROOT_TPrimaryGenerator
#define ROOT_TPrimaryGenerator

#include "TPartIndex.h"

class TParticlePDG;
class GeantTrack;


class TPrimaryGenerator{
 private:
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

 public:
  TPrimaryGenerator();
  TPrimaryGenerator(Int_t partpdg, Double_t partekin, Double_t xpos, 
                    Double_t ypos, Double_t zpos, Double_t xdir, Double_t ydir,
                    Double_t zdir);
  virtual ~TPrimaryGenerator();
  
  void SetEkin(Double_t partekin){ fPartEkin = partekin; }
  Double_t GetEkin()const {return fPartEkin;}
  void  SetParticleByPDGCode(Int_t pdgcode);
  Int_t GetParticlePDGCode()const{return fPDG;} 
  void  SetParticleXYZPosition(Double_t x, Double_t y, Double_t z){fXPos=x; fYPos=y; fZPos=z;} 
  Double_t GetParticleXPos()const{return fXPos;}
  Double_t GetParticleYPos()const{return fYPos;}
  Double_t GetParticleZPos()const{return fZPos;}
  void  SetParticleXYZDir(Double_t xdir, Double_t ydir, Double_t zdir);
  Double_t GetParticleXDir()const{return fXDir;}
  Double_t GetParticleYDir()const{return fYDir;}
  Double_t GetParticleZDir()const{return fZDir;}

  // --
  Int_t GetParticleGVIndex()    const{return fGVPartIndex;} 
  TParticlePDG* GetParticlePDG()const{return fPartPDG;}
  Double_t GetParticleMass()    const{return fMass;}
  Double_t GetParticleCharge()  const{return fCharge;}
  Double_t GetparticlePTotal()  const{return fPTotal;}
  Double_t GetparticleETotal()  const{return fETotal;}
  

  // set one GeantTrack primary track properties
  void InitPrimaryTrack(GeantTrack &gtrack);

 private:
   void InitPrimaryGenerator();

   TPrimaryGenerator(const TPrimaryGenerator &);//no imp.	
   TPrimaryGenerator& operator=(const TPrimaryGenerator &);//no imp.

   ClassDef(TPrimaryGenerator,1)
};

#endif
