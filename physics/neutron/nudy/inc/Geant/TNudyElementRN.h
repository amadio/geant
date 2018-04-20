#ifndef ROOT_TNudyElementRN
#define ROOT_TNudyElementRN

/*
#include "RConfigure.h"
#include "Riostream.h"
#include "TBox.h"
#include "TObject.h"
#include "TList.h"
#include "TColor.h"
#include "Rtypes.h"
#include "TColor.h"
#include "TRint.h"
#include <string.h>
#include "TStyle.h"
#include "TROOT.h"
*/

#include "TGeoElement.h"
class TBox;
#include "TPaveText.h"
#include "TColor.h"

namespace Nudy {
class TNudyElementRN : public TObject {

private:
  int fCoSize;
  float fX, fY;
  TGeoElementRN *fEle;
  TBox *fBox;
  TPaveText fInfo;
  Color_t GetColor(double halfLife);

  static double fCCodeRange[26];
  static int fCCodeColor[26][3];

  float fSize;
  float fScale;
  float fPadding;

  TColor fColors[26];

public:
  TNudyElementRN();
  TNudyElementRN(TGeoElementRN *elem, float fX, float fY);
  virtual ~TNudyElementRN() { delete fBox; }
  void Draw(const char *option = "");
  void Move(float x, float y);
  void SetColorCode(TList *cCodeRange, TList *cCodeColor);
  int GetA() { return fEle->AtomicNo(); }
  int GetZ() { return fEle->MassNo(); }
  void SetScale(float scale) { fScale = scale; }

#ifdef USE_ROOT
  ClassDef(TNudyElementRN, 1) // Radio Nucleide Element
#endif
};

} // namespace
#endif
