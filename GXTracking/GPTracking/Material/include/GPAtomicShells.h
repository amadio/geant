#ifndef GPATOMICSHELLS_H
#define GPATOMICSHELLS_H

#include "GPTypeDef.h"

namespace GPAtomicShells 
{

   extern VARTYPE G4int fNumberOfShells[101];

// The total shell number is:
// 1 + G4AtomicShells::TotalNumberOfShells(100) = 1 + 1539 = 1540 

   extern VARTYPE G4int fIndexOfShells[101];

   extern VARTYPE G4double fBindingEnergies[1540];

   extern VARTYPE G4int fNumberOfElectrons[1540];

}; //end of name space

extern "C" {

  FQUALIFIER
  G4int    GPAtomicShells_GetNumberOfShells(G4int Z);

  FQUALIFIER
  G4int    GPAtomicShells_GetNumberOfElectrons(G4int Z, G4int SubshellNb);

  FQUALIFIER
  G4double GPAtomicShells_GetBindingEnergy (G4int Z, G4int SubshellNb);

  FQUALIFIER
  G4double GPAtomicShells_GetTotalBindingEnergy (G4int Z);
}

#endif   // end of G4AtomicShells.hh 
