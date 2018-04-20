#ifndef ROOT_TNudyCore
#define ROOT_TNudyCore

class TGeoManager;
class TGeoElementRN;
class THashList;
class TParticle;
#include "TGeoElement.h"
#include "TDatabasePDG.h"
#include "TNudyTypes.h"

namespace Nudy {

class TNudyCore {
protected:
  static TNudyCore *fgInstance;
  TGeoManager *fGeom;
  TGeoElementTable *fTable;
  TDatabasePDG *fPdgDB;
  TNudyCore();
  TNudyCore(const TNudyCore &core) = delete;            // not implemented
  TNudyCore &operator=(const TNudyCore &core) = delete; // not implemented

public:
  virtual ~TNudyCore();         // Public Destructor
  static TNudyCore *Instance(); // Returns Instance of TNudyManager
  // Calculation functions
  double ThinningDuplicate(std::vector<double> &x1);
  double ThinningDuplicate(std::vector<double> &x1, std::vector<double> &x2);
  int BinarySearch(std::vector<double> array, int len, double val);
  double Interpolate(std::vector<int> nbt1, std::vector<int> int1, int nr, std::vector<double> x, std::vector<double> y,
                     int np, double xx);
  void cdfGenerateT(std::vector<double> &x1, std::vector<double> &x2, std::vector<double> &x3);
  double CmToLabElasticE(double inE, double cmCos, double awr);
  double CmToLabElasticCosT(double cmCos, double awr);
  double CmToLabInelasticE(double cmEOut, double inE, double cmCos, double awr);
  double CmToLabInelasticCosT(double labEOut, double cmEOut, double inE, double cmCos, double awr);
  void Sort(std::vector<double> &x1, std::vector<double> &x2);
  double LinearInterpolation(double x1, double y1, double x2, double y2, double x); // Linear Interpolation
  double BilinearInterploation(double x1, double y1, double x2, double y2, double z11, double z12, double z21,
                               double z22, double x, double y); // Biliniear Interpolation
  void TrapezoidalIntegral(double *xpts, double *ypts, const int npts,
                           double *&out); // Calculates integral of discrete points
  void CumulativeIntegral(double *x, double *y, double *q, int len);
  int BinarySearch(double *array, int len, double val);
  double InterpolateScale(double x[2], double y[2], int law, double xx);
  double Interpolate(int *nbt, int *interp, int nr, double *x, double *y, int np, double xx);
  char *ExpandReaction(Reaction_t reac);
  // Model Key checking/generation functions
  int IsMaterial(const TGeoElementRN *endf, const char *key);
  int IsTemperature(const unsigned long temp, const char *key);
  int IsReaction(const Reaction_t r, const char *key);
  char *GetKey(const TGeoElementRN *mat, Reaction_t reac, unsigned long temp);

  // Uniform access functions
  TGeoElementTable *GetElementTable() { return fTable; }
  TDatabasePDG *GetDatabasePDG() { return fPdgDB; }
  TGeoElementRN *GetMaterial(int ENDFcode) { return fTable->GetElementRN(ENDFcode); }
  TGeoElementRN *GetMaterial(int a, int z, int iso = 0) { return fTable->GetElementRN(a, z, iso); }
  const THashList *GetParticleList() { return fPdgDB->ParticleList(); }
  TParticlePDG *GetParticlePDG(int pdgCode);
  TParticlePDG *GetParticlePDG(const char *name);
  TParticle *GetParticle(int pdgCode);
  TParticle *GetParticle(const char *name);
  void MemProfile();

#ifdef USE_ROOT
  ClassDef(TNudyCore, 1)
#endif
};

} // namespace
#endif
