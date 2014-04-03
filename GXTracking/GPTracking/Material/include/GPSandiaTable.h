#ifndef GPSANDIATABLE_HH
#define GPSANDIATABLE_HH

#include "GPTypeDef.h"
#include "GPMaterial.h"

struct GPSandiaTable
{
  GPMaterial* fMaterial;
  G4int       fCumulInterval[101];
  G4double    fSandiaCofPerAtom[4];
  G4double    fnulcof[4];

  // G4int           fMatNbOfIntervals;
  // G4OrderedTable* fMatSandiaMatrix;
};

extern "C" {

  FQUALIFIER
  void GPSandiaTable_Constructor(GPSandiaTable* This);

  FQUALIFIER
  G4double* GPSandiaTable_GetSandiaCofPerAtom(GPSandiaTable*This, 
					      G4int Z, G4double energy);

  FQUALIFIER
  G4double* GPSandiaTable_GetSandiaCofForMaterial(GPSandiaTable*This,
						  G4double energy);

  FQUALIFIER
  G4double GPSandiaTable_GetZtoA(G4int Z);

  FQUALIFIER
  G4int GPSandiaTable_GetNbOfIntervals(G4int Z);

  FQUALIFIER
  G4double GPSandiaTable_GetIonizationPot(G4int Z);
}

#endif 
