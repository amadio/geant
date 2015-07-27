#ifndef ROOT_TNudyElementTable
#define ROOT_TNudyElementTable

#include "RConfigure.h"
#include "Riostream.h"
#include "TGeoElement.h"
#include "TGeoManager.h"
#include "TNudyElementRN.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TBox.h"
#include "TList.h"
#include "TButton.h"
const int kENABLED = 1;
const int kDISABLED = 0;

class TNudyElementTable: public TObject{

 private:
  int fState;
  float fOx,fOy;
  TGeoElementTable* fTable;
  TBox fWindow;
  TList fEleBox;
  float fLOD;
  TCanvas* fRNTable;
  TGeoManager* fGeom;
  void DrawUI();
  TList fControls;
  void InitializeControls();
 public:
  TNudyElementTable();
  virtual ~TNudyElementTable();
  void Draw(Option_t *option="");
  void ZoomIn();
  void ZoomOut();
  void MoveUp();
  void MoveDown();
  void MoveLeft();
  void MoveRight();
  void Update();
  ClassDef(TNudyElementTable,1) //Table of RadioNucleides
};
#endif
 
