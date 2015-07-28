#ifndef ROOT_TNudyCore
#define ROOT_TNudyCore

#include <TObject.h>
#include <TGeoElement.h>
#include <TGeoManager.h>
#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TList.h>
#include <TSystem.h>

#include "TNudyTypes.h"

//Set to turn on multi threaded generation of random numbers
#define TNUDYALIAS_MULTITHREAD

class TNudyCore : public TNamed {
 protected:
  static TNudyCore* fgInstance;
  TGeoManager *fGeom;
  TGeoElementTable *fTable;
  TDatabasePDG *fPdgDB;
  TList *fListOfObjects;
  TNudyCore();  
  TNudyCore(const TNudyCore& core): TNamed(core){}
  TNudyCore& operator=(const TNudyCore& core)
    {if(this!=&core){TNamed::operator=(core);}return *this;}
  
 public:

  virtual ~TNudyCore();//Public Destructor
  static TNudyCore* Instance();//Returns Instance of TNudyManager
  //Calculation functions
  Double_t LinearInterpolation(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x);//Linear Interpolation
  Double_t BilinearInterploation(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t z11,Double_t z12,Double_t z21,Double_t z22,Double_t x,Double_t y);//Biliniear Interpolation
  void TrapezoidalIntegral(Double_t *xpts, Double_t *ypts,const Int_t npts,Double_t *out);//Calculates integral of discrete points
  void CumulativeIntegral(Double_t *x, Double_t *y, Double_t *q, Int_t len);
  Int_t BinarySearch(Double_t *array,Int_t len, Double_t val);
  Double_t InterpolateScale(Double_t x[2],Double_t y[2],Int_t law, Double_t xx);
  Double_t Interpolate(Int_t *nbt, Int_t *interp, Int_t nr, Double_t *x, Double_t*y, Int_t np, Double_t xx);
  char * ExpandReaction(Reaction_t reac);
  //Model Key checking/generation functions
  Int_t IsMaterial(const TGeoElementRN* endf, const char* key);
  Int_t IsTemperature(const ULong_t temp,const char *key);
  Int_t IsReaction(const Reaction_t r, const char *key);
  char * GetKey(const TGeoElementRN *mat,Reaction_t reac,ULong_t temp);

  //Uniform access functions
  TGeoElementTable* GetElementTable(){return fTable;}
  TDatabasePDG* GetDatabasePDG(){return fPdgDB;}
  TGeoElementRN* GetMaterial(Int_t ENDFcode) { return fTable->GetElementRN(ENDFcode);}
  TGeoElementRN* GetMaterial(Int_t a, Int_t z, Int_t iso = 0) { return fTable->GetElementRN(a,z,iso);}
  const THashList* GetParticleList(){return fPdgDB->ParticleList();}
  TParticlePDG* GetParticlePDG(Int_t pdgCode);
  TParticlePDG* GetParticlePDG(const char* name);
  TParticle* GetParticle(Int_t pdgCode);
  TParticle* GetParticle(const char* name);			 
  void MemProfile();
  
  ClassDef(TNudyCore,1)
};

#endif
