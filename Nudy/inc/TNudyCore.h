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
  double LinearInterpolation(double x1,double y1,double x2,double y2,double x);//Linear Interpolation
  double BilinearInterploation(double x1,double y1,double x2,double y2,double z11,double z12,double z21,double z22,double x,double y);//Biliniear Interpolation
  void TrapezoidalIntegral(double *xpts, double *ypts,const Int_t npts,double *out);//Calculates integral of discrete points
  void CumulativeIntegral(double *x, double *y, double *q, Int_t len);
  Int_t BinarySearch(double *array,Int_t len, double val);
  double InterpolateScale(double x[2],double y[2],Int_t law, double xx);
  double Interpolate(Int_t *nbt, Int_t *interp, Int_t nr, double *x, double*y, Int_t np, double xx);
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
