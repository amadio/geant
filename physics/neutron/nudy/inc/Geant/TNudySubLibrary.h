#ifndef ROOT_TNudySubLibrary
#define ROOT_TNudySubLibrary

#include "TNamed.h"
class TVNudyModel;
class TNudyEndfMat;
#include "Geant/TNudyTypes.h"
class TBtree;
class TParticlePDG;
class TGeoElementRN;

class TNudySubLibrary : public TNamed {
public:
  TNudySubLibrary();
  TNudySubLibrary(TParticlePDG *projectile);
  virtual ~TNudySubLibrary();
  void AddModel(TVNudyModel *model);
  void ReadMat(TNudyEndfMat *material);
  void SetProjectile(TParticlePDG *particle) { fProjectile = particle; }
  void ListModels();
  TParticlePDG *GetProjectile() { return fProjectile; }
  TVNudyModel *GetModel(const TGeoElementRN *mat, const Reaction_t reac, const unsigned long temp);
  TBtree *GetAllModels(const TGeoElementRN *mat = NULL, const Reaction_t reac = kNoReaction,
                       const unsigned long temp = 0);

private:
  TBtree *fIndex;            // Btree storing all Models
  TBtree *fBuffer;           // Buffer storing results of last query
  TParticlePDG *fProjectile; // Projectile particle
#ifdef USE_ROOT
  ClassDef(TNudySubLibrary, 1)
#endif
};
#endif
