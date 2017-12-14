#ifndef ROOT_TVNudyModel
#define ROOT_TVNudyModel

#include <TNamed.h>
#include "TArrayD.h"
#include "TArrayI.h"
#include "TRandom3.h"
class TNudyAlias;
class TNudyAliasCont;
class TNudyEndfFile;
class TNudyEndfSec;
class TNudyEndfTab1;
class TNudyEndfTab2;
class TGeoElementRN;
class TParticlePDG;
class TNudyEndfMat;
#include "TNudyTypes.h"

class TVNudyModel : public TNamed {
public:
  TVNudyModel();
  TVNudyModel(TGeoElementRN *mat, Reaction_t reac, unsigned long temp, TParticlePDG *projectile,
              TNudyEndfMat *material);
  virtual ~TVNudyModel();

  void ReadFile(TNudyEndfMat *material);
  void DisplayData(FileData_t file);
  void DumpData(FileData_t file);
  TGeoElementRN *GetMaterial();
  const char *GetMaterialName();
  int GetZ();
  int GetA();
  int GetISO();
  Reaction_t GetReaction();
  unsigned long GetTemp();
  virtual double GetXSect(double e);
  virtual double GetEo(double ein);
  virtual double GetAo(double ein);
  TArrayD *GetFile5Data(double ein);
  TArrayD *GetFile5ProcessedData(double ein);

private:
  int fMAT;                  // Material number
  unsigned long fTemp;       // Temperature for evaluation of data
  int fEndf;                 // Endf code of the material
  int fPdg;                  // Pdgcode of the projectile
  TGeoElementRN *fMaterial;  //! Material
  Reaction_t fReaction;      // Reaction
  TParticlePDG *fProjectile; //! Projectile for reaction

  // Data from RENDF file
  int fEXSect_length; // Length of two arrays of energy and cross-section from file 3
  // Array to store data about energy
  double *fE_file3; //[fEXSect_length]
  // Array to store data about XSects f:(fE_file3)->(fXSect_file3)
  double *fXSect_file3; //[fEXSect_length]

  // File 4
  int f4nens;
  TArrayD f4eins;
  TNudyAliasCont *fAPAlias; //[f4nens]
  double f4Tein;            //!
  double f4Tel;             //!

  // File 5
  TRandom3 fRnd; //!
  TArrayD xengr; //
  int nens;      //
  int nperc;     //
  int maxpop;    //!
  TArrayI nEout; //
  //  double EPtable[4000][400]; //!
  TArrayD *fEPtable; //[nens]
  // double EPtable[1100][200]; //!
  TArrayD *fPerc;           //[nens]
  double f5Tein;            //!
  double f5Tel;             //!
  TNudyAliasCont *fEPAlias; //[nens]

  // Functions to read data from RENDF format
  void ReadFile3(TNudyEndfFile *file);
  void ReadFile4(TNudyEndfFile *file);
  void ReadFile5(TNudyEndfFile *file);
  void File5_Pass1(TNudyEndfSec *sec);
  void File5_Pass2(TNudyEndfSec *sec);
  void File5_Pass3();

  // Function to check consistency
  int CheckLinear(TNudyEndfTab1 *tab);
  int CheckLinear(TNudyEndfTab2 *tab);
  void Linearize(TNudyEndfTab1 *tab);
  // void Compare(TNudyEndfTab1* oldtab, TNudyEndfTab1 *newtab);
  void PopulateGrid(int index);
  void FillGrid(double u, int nep, TNudyEndfTab1 *tab, TNudyEndfTab1 *pe);
  int EoExists(int index, double ef);

#ifdef USE_ROOT
  ClassDef(TVNudyModel, 1)
#endif
};

#endif
