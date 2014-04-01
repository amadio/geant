// ClassName:   G4Element

#ifndef GPELEMENT_HH
#define GPELEMENT_HH 1

#include "GPTypeDef.h"

struct GPElement
{
  G4double fZeff;              // Effective atomic number
  G4double fNeff;              // Effective number of nucleons
  G4double fAeff;              // Effective mass of a mole

  // Derived data members (computed from the basic data members)
  G4double fCoulomb;           // Coulomb correction factor
  G4double fRadTsai;           // Tsai formula for the radiation length
};

extern "C" {

  FQUALIFIER
  void GPElement_Constructor(GPElement *This, G4double zeff, G4double aeff);

  FQUALIFIER  
  void GPElement_ComputeDerivedQuantities(GPElement *This);

  FQUALIFIER  
  void GPElement_ComputeCoulombFactor(GPElement *This);

  FQUALIFIER  
  void GPElement_ComputeLradTsaiFactor(GPElement *This);

  FQUALIFIER
  G4double GPElement_GetZ(GPElement *This); 

  FQUALIFIER
  G4double GPElement_GetN(GPElement *This); 

  FQUALIFIER
  G4double GPElement_GetA(GPElement *This); 

  FQUALIFIER
  G4double GPElement_GetfCoulomb(GPElement *This); 

  FQUALIFIER
  G4double GPElement_GetfRadTsai(GPElement *This);
}

#endif
