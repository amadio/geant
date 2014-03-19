#include "GPUtils.h"

FQUALIFIER
int GPimax(int x, int y)
{
  if(x>y) return x;
  else return y;
}

FQUALIFIER
G4double GPfmax(G4double x, G4double y)
{
  if(x>y) return x;
  else return y;
}

FQUALIFIER
G4double GPfmin(G4double x, G4double y)
{
  if(x<y) return x;
  else return y;
}

FQUALIFIER
G4double GPsqr(G4double x)
{
  return x*x;
}
