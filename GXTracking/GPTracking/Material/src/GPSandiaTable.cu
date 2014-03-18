#include "GPSandiaTable.h"
#include "GPStaticSandiaData.h"
#include "GPMaterial.h"
#include "GPPhysicalConstants.h"

#include <assert.h>

using namespace GPSandiaData;
 
FQUALIFIER
void GPSandiaTable_Constructor(GPSandiaTable* This /* GPMaterial *material */)
{
  This->fMaterial = 0; //material;

  //build the CumulInterval array

  This->fCumulInterval[0] = 1;

  for (G4int Z=1; Z<101; ++Z) {
    This->fCumulInterval[Z] = This->fCumulInterval[Z-1] + fNbOfIntervals[Z];
  }
  
  //initialisation of fnulcof
  for (G4int i=0; i<4; ++i)  This->fnulcof[i] = 0.0;

  //compute macroscopic Sandia coefs for a material   
  //  fMatNbOfIntervals   = 0;
  //  fMatSandiaMatrix    = 0; 
  //  ComputeMatSandiaMatrix(); // mma

}
							
FQUALIFIER
G4double* GPSandiaTable_GetSandiaCofPerAtom(GPSandiaTable*This, 
					    G4int Z, G4double energy)
{

  G4double Emin  = fSandiaTable[This->fCumulInterval[Z-1]][0]*keV;
  G4double Iopot = fIonizationPotentials[Z]*eV;
  if (Iopot > Emin) Emin = Iopot;
   
  G4int interval = fNbOfIntervals[Z] - 1;
  G4int row = This->fCumulInterval[Z-1] + interval;
  while ((interval>0) && (energy<fSandiaTable[row][0]*keV)) {
    --interval;
    row = This->fCumulInterval[Z-1] + interval;
  }
  if (energy >= Emin)
    {        
      G4double AoverAvo = Z*amu/fZtoAratio[Z];
         
      This->fSandiaCofPerAtom[0]=AoverAvo*funitc[1]*fSandiaTable[row][1];     
      This->fSandiaCofPerAtom[1]=AoverAvo*funitc[2]*fSandiaTable[row][2];     
      This->fSandiaCofPerAtom[2]=AoverAvo*funitc[3]*fSandiaTable[row][3];     
      This->fSandiaCofPerAtom[3]=AoverAvo*funitc[4]*fSandiaTable[row][4];
    }
  else 
    {
      This->fSandiaCofPerAtom[0] = This->fSandiaCofPerAtom[1] 
	= This->fSandiaCofPerAtom[2] = This->fSandiaCofPerAtom[3] = 0.;
    }                
  return This->fSandiaCofPerAtom;     
}

FQUALIFIER
G4double* GPSandiaTable_GetSandiaCofForMaterial(GPSandiaTable*This,
						 G4double energy)
{
  G4double* x = &(This->fnulcof[0]);
  //@@@ this is dummy implementation and should be extended if necessary
  /*
  if (energy >= (*(*(This->fMatSandiaMatrix))[0])[0]) {
   
    G4int interval = This->fMatNbOfIntervals - 1;
    while ((interval>0)&&(energy<(*(*(This->fMatSandiaMatrix))[interval])[0])) 
      {interval--;} 
    x = &((*(*(This->fMatSandiaMatrix))[interval])[1]);
  }
  */
  return x;
}

FQUALIFIER
G4double GPSandiaTable_GetZtoA(G4int Z)
{
  assert (Z>0 && Z<101);
  return fZtoAratio[Z];
}

FQUALIFIER
G4int GPSandiaTable_GetNbOfIntervals(G4int Z)
{
  assert (Z>0 && Z<101);  
  return fNbOfIntervals[Z];
}

FQUALIFIER
G4double GPSandiaTable_GetIonizationPot(G4int Z)
{
  assert (Z>0 && Z<101);
  return fIonizationPotentials[Z]*eV;
}

