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

class TNudyElementRN : public TObject {

private:
  int fCoSize;
  float fX, fY;
  TGeoElementRN *fEle;
  TBox *fBox;
  TPaveText fInfo;
  Color_t GetColor(double halfLife);

public:
  static double fCCodeRange[26];
  static int fCCodeColor[26][3];

  float fSize;
  float fScale;
  float fPadding;

  TColor fColors[26];
  TNudyElementRN();
  TNudyElementRN(TGeoElementRN *elem, float fX, float fY);
  virtual ~TNudyElementRN(){delete fBox;}
  void Draw(const char *option = "");
  void Move(float x, float y);
  void SetColorCode(TList *cCodeRange, TList *cCodeColor);
  int GetA() { return fEle->AtomicNo(); }
  int GetZ() { return fEle->MassNo(); }

#ifdef USE_ROOT
  ClassDef(TNudyElementRN, 1) // Radio Nucleide Element
#endif
};
#endif
