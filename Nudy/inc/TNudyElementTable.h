#ifndef ROOT_TNudyElementTable
#define ROOT_TNudyElementTable

class TCanvas;
class TGeoElementTable;
class TGeoManager;
#include "TBox.h"
#include "TList.h"


const int kENABLED = 1;
const int kDISABLED = 0;

class TNudyElementTable {

private:
  int fState;
  float fOx, fOy;
  TGeoElementTable *fTable;
  TBox fWindow;
  TList fEleBox;
  float fLOD;
  TCanvas *fRNTable;
  TGeoManager *fGeom;
  void DrawUI();
  TList fControls;
  void InitializeControls();

public:
  TNudyElementTable();
  virtual ~TNudyElementTable();
  void Draw(const char *option = "");
  void ZoomIn();
  void ZoomOut();
  void MoveUp();
  void MoveDown();
  void MoveLeft();
  void MoveRight();
  void Update();
#ifdef USE_ROOT
  ClassDef(TNudyElementTable, 1) // Table of RadioNucleides
#endif
};
#endif
