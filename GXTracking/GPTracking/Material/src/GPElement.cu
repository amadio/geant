#include "GPElement.h"
#include "GPPhysicalConstants.h"
#include "GPSystemOfUnits.h"

// Constructor to Generate an element from scratch

FQUALIFIER
void GPElement_Constructor(GPElement *This, G4double zeff, G4double aeff)
{
  This->fZeff   = zeff;
  This->fNeff   = aeff/(g/mole);
  This->fAeff   = aeff;

  if(This->fNeff < 1.0) This->fNeff = 1.0;

  GPElement_ComputeDerivedQuantities(This);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FQUALIFIER
void GPElement_ComputeDerivedQuantities(GPElement *This)
{
  // some basic functions of the atomic number
  // Radiation Length
  GPElement_ComputeCoulombFactor(This);
  GPElement_ComputeLradTsaiFactor(This); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FQUALIFIER
void GPElement_ComputeCoulombFactor(GPElement *This)
{
  //
  //  Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)

  const G4double k1 = 0.0083 , k2 = 0.20206 ,k3 = 0.0020 , k4 = 0.0369 ;

  G4double az2 = 
    (fine_structure_const*This->fZeff)*(fine_structure_const*This->fZeff);
  G4double az4 = az2 * az2;

  This->fCoulomb = (k1*az4 + k2 + 1./(1.+az2))*az2 - (k3*az4 + k4)*az4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FQUALIFIER
void GPElement_ComputeLradTsaiFactor(GPElement *This)
{
  //  Compute Tsai's Expression for the Radiation Length
  //  (Phys Rev. D50 3-1 (1994) page 1254)

  const G4double Lrad_light[]  = {5.31  , 4.79  , 4.74 ,  4.71} ;
  const G4double Lprad_light[] = {6.144 , 5.621 , 5.805 , 5.924} ;
  
  const G4double logZ3 = log(This->fZeff)/3.;

  G4double Lrad, Lprad;
  G4int iz = (G4int)(This->fZeff+0.5) - 1 ;

  if (iz <= 3) { Lrad = Lrad_light[iz] ;  Lprad = Lprad_light[iz] ; }
  else { Lrad = log(184.15) - logZ3 ; Lprad = log(1194.) - 2*logZ3;}

  This->fRadTsai = 
    4*alpha_rcl2*This->fZeff*(This->fZeff*(Lrad-This->fCoulomb) + Lprad); 
}

FQUALIFIER
G4double GPElement_GetZ(GPElement *This) 
{
  return This->fZeff;
}

FQUALIFIER
G4double GPElement_GetN(GPElement *This) 
{
  return This->fNeff;
}

FQUALIFIER
G4double GPElement_GetA(GPElement *This) 
{
  return This->fAeff;
}

FQUALIFIER
G4double GPElement_GetfCoulomb(GPElement *This) 
{
  return This->fCoulomb;
}

FQUALIFIER
G4double GPElement_GetfRadTsai(GPElement *This) 
{
  return This->fRadTsai;
}

