#ifndef ROOT_TNudySubLibrary
#define ROOT_TNudySubLibrary

#include <TBtree.h>
#include <TList.h>
#include "TNudyCore.h"
#include "TNudyEndfMat.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfSec.h"
#include "TVNudyModel.h"

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
  ClassDef(TNudySubLibrary, 1)
};
#endif
