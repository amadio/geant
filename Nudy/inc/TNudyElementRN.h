#ifndef ROOT_TNudyElementRN
#define ROOT_TNudyElementRN

#include "RConfigure.h"
#include "Riostream.h"
#include "TGeoElement.h"
#include "TBox.h"
#include "TObject.h"
#include "TPaveText.h"
#include "TList.h"
#include "TColor.h"
#include "Rtypes.h"
#include "TColor.h"
#include "TRint.h"
#include <string.h>
#include "TStyle.h"
#include "TROOT.h"

class TNudyElementRN: public TObject{

 private:
  int fCoSize;
  Float_t fX,fY;
  TGeoElementRN* fEle;
  TBox* fBox;
  TPaveText fInfo;
  Color_t GetColor(double halfLife);
 public:
  static double fCCodeRange[26];
  static Int_t fCCodeColor[26][3];

  Float_t fSize;
  Float_t fScale;
  Float_t fPadding;
  
  TColor fColors[26];
  TNudyElementRN();
  TNudyElementRN(TGeoElementRN* elem,Float_t fX, Float_t fY);
  virtual ~TNudyElementRN() {};
  void Draw(Option_t* option="");
  void Move(Float_t x, Float_t y);
  void SetColorCode(TList* cCodeRange, TList* cCodeColor);
  Int_t GetA() {return fEle->AtomicNo();}
  Int_t GetZ() {return fEle->MassNo();}
  
  ClassDef(TNudyElementRN,1) //Radio Nucleide Element
};
#endif
