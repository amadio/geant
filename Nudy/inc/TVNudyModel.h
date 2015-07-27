#ifndef ROOT_TVNudyModel
#define ROOT_TVNudyModel

#include <TObjArray.h>
#include <TArrayD.h>
#include <TList.h>
#include <TArrayI.h>
#include <TNamed.h>
#include "TNudyCore.h"
#include "TNudyEndfMat.h"
#include "TNudyEndfSec.h"
#include "TNudyEndfTab1.h"
#include "TNudyEndfTab2.h"
#include "TNudyEndfList.h"
#include "TNudyEndfINTG.h"
#include "TNudyAliasCont.h"
#include "TRandom3.h"

class TVNudyModel : public TNamed {
 public:
  TVNudyModel();
  TVNudyModel(TGeoElementRN *mat, Reaction_t reac, ULong_t temp,TParticlePDG* projectile ,TNudyEndfMat *material);
  virtual ~TVNudyModel();

  void ReadFile(TNudyEndfMat *material);
  void DisplayData(FileData_t file);
  void DumpData(FileData_t file);
  TGeoElementRN* GetMaterial();
  const char * GetMaterialName();
  Int_t GetZ();
  Int_t GetA();
  Int_t GetISO();
  Reaction_t GetReaction();
  ULong_t GetTemp();
  virtual double GetXSect(double e);
  virtual double GetEo(double ein);
  virtual double GetAo(double ein);
  TArrayD* GetFile5Data(double ein);
  TArrayD* GetFile5ProcessedData(double ein);

 private:
  Int_t fMAT; //Material number
  ULong_t fTemp; //Temperature for evaluation of data
  Int_t fEndf; // Endf code of the material
  Int_t fPdg;  // Pdgcode of the projectile
  TGeoElementRN *fMaterial; //! Material
  Reaction_t fReaction; // Reaction
  TParticlePDG* fProjectile; //! Projectile for reaction
  

  //Data from RENDF file
  Int_t fEXSect_length;//Length of two arrays of energy and cross-section from file 3
  //Array to store data about energy
  double *fE_file3;//[fEXSect_length]
  //Array to store data about XSects f:(fE_file3)->(fXSect_file3)
  double *fXSect_file3;//[fEXSect_length]
  
  //File 4
  Int_t f4nens; 
  TArrayD f4eins;
  TNudyAliasCont *fAPAlias; //[f4nens]
  double f4Tein; //!
  double f4Tel; //!

  //File 5  
  TRandom3 fRnd; //!
  TArrayD xengr; //
  Int_t nens; //
  Int_t nperc; //
  Int_t maxpop; //!
  TArrayI nEout; //
  //  double EPtable[4000][400]; //!
  TArrayD *fEPtable; //[nens]
  //double EPtable[1100][200]; //!
  TArrayD *fPerc; //[nens]
  double f5Tein; //!
  double f5Tel; //!
  TNudyAliasCont *fEPAlias;  //[nens]

  //Functions to read data from RENDF format
  void ReadFile3(TNudyEndfFile *file);
  void ReadFile4(TNudyEndfFile *file);
  void ReadFile5(TNudyEndfFile *file);
  void File5_Pass1(TNudyEndfSec* sec);
  void File5_Pass2(TNudyEndfSec* sec);
  void File5_Pass3();


  //Function to check consistency
  Int_t CheckLinear(TNudyEndfTab1 *tab);
  Int_t CheckLinear(TNudyEndfTab2 *tab);
  void Linearize(TNudyEndfTab1 *tab);
  //void Compare(TNudyEndfTab1* oldtab, TNudyEndfTab1 *newtab);
  void PopulateGrid(Int_t index);
  void FillGrid(double u, Int_t nep, TNudyEndfTab1 *tab, TNudyEndfTab1 *pe);
  Int_t EoExists(Int_t index, double ef);
  
  
  ClassDef(TVNudyModel,1)
};

#endif
