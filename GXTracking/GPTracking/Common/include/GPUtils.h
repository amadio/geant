#ifndef GPUtils_hh
#define GPUtils_hh 1

#include "GPTypeDef.h"

extern "C" {

FQUALIFIER
int GPimax(int x, int y);

FQUALIFIER
G4double GPfmax(G4double x, G4double y);

FQUALIFIER
G4double GPfmin(G4double x, G4double y);

FQUALIFIER
G4double GPsqr(G4double x);
}

#endif
