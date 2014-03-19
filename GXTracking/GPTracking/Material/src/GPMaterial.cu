#include "GPMaterial.h"
#include <assert.h>
#ifdef GPUDEBUG
  #ifndef __CUDA_ARCH__
    #include <iostream>
  #endif
#endif

FQUALIFIER
void GPMaterial_Constructor(GPMaterial *This,
			    G4double density, 
			    G4double a, G4double z)
{
  This->fIndex = -1;
  This->fNumberOfElements = 1;

  This->fDensity  = density;
  This->fA  = a;
  This->fZ  = z;
  This->fRadlen = 0;

  This->TotNbOfElectPerVolume = 0.;
  This->TotNbOfElectPerVolume = 0.;

  GPElement aElement;
  GPElement_Constructor(&aElement,z,a);

  This->fElementVector[0] = aElement;
  This->fMassFractionVector[0] = 1.0;

  GPMaterial_ComputeDerivedQuantities(This);

}

FQUALIFIER
void GPMaterial_Constructor_ByElement(GPMaterial *This, G4double density)
{
  This->fIndex = -1;
  This->fNumberOfElements = 0;

  This->fDensity  = density;
  This->fA  = 0;
  This->fZ  = 0;
  This->fRadlen = 0;

  This->TotNbOfElectPerVolume = 0.;
  This->TotNbOfElectPerVolume = 0.;

  for (G4int i = 0 ; i < kElementMax ; ++i) {
    This->fMassFractionVector[i] = 0.0;
    This->VecNbOfAtomsPerVolume[i] = 0.0;
  }

}

FQUALIFIER
G4int GPMaterial_GetIndex(GPMaterial *This)
{
  return This->fIndex;
}

FQUALIFIER
void GPMaterial_SetIndex(GPMaterial *This, G4int index)
{
  This->fIndex = index;
}

FQUALIFIER
size_t GPMaterial_GetNumberOfElements(GPMaterial *This)
{
  return This->fNumberOfElements;
}

FQUALIFIER
G4double GPMaterial_GetDensity(GPMaterial *This)
{
  return This->fDensity;
}

FQUALIFIER
G4double GPMaterial_GetZ(GPMaterial *This)
{
  return This->fZ;
}

FQUALIFIER
G4double GPMaterial_GetA(GPMaterial *This)
{
  return This->fA;
}

FQUALIFIER
G4double GPMaterial_GetRadlen(GPMaterial *This)
{
  return This->fRadlen;
}

FQUALIFIER
void GPMaterial_ComputeDerivedQuantities(GPMaterial *This)
{
  // Header routine to compute various properties of material.

  // Number of atoms per volume (per element), total nb of electrons per volume
  G4double Zi, Ai;

  This->TotNbOfElectPerVolume = 0.;
  This->TotNbOfElectPerVolume = 0.;

  for (size_t i=0; i<This->fNumberOfElements; ++i) {
    Zi = GPElement_GetZ(&(This->fElementVector[i]));
    Ai = GPElement_GetA(&(This->fElementVector[i]));
    This->VecNbOfAtomsPerVolume[i] = 
      Avogadro*This->fDensity*This->fMassFractionVector[i]/Ai;
    This->TotNbOfAtomsPerVolume += This->VecNbOfAtomsPerVolume[i];
    This->TotNbOfElectPerVolume += This->VecNbOfAtomsPerVolume[i]*Zi;
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
    std::cout<<"GPMaterial: CompDerivedQuantities: "<< Zi <<' '<< Ai <<' '<< This->VecNbOfAtomsPerVolume[i]
	     <<' '<< This->TotNbOfAtomsPerVolume <<' '<< This->TotNbOfElectPerVolume << std::endl;
#endif
#endif
  }
        
  GPMaterial_ComputeRadiationLength(This);
}

FQUALIFIER
void GPMaterial_ComputeRadiationLength(GPMaterial *This)
{
  G4double radinv = 0.0 ;
  for (size_t i=0;i<This->fNumberOfElements;++i) {
    radinv += This->VecNbOfAtomsPerVolume[i]*GPElement_GetfRadTsai(&(This->fElementVector[i]));
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
    std::cout<<"GPMaterial::ComputeRadiationLength: i="<< i <<" radinv="<< radinv
	     <<", NbAtosPerVol="<< This->VecNbOfAtomsPerVolume[i]
	     <<", Elem fRadTsai="<< GPElement_GetfRadTsai(&(This->fElementVector[i])) << std::endl;
#endif
#endif
  }
  This->fRadlen = (radinv <= 0.0 ? DBL_MAX : 1./radinv);
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
  std::cout<<"GPMaterial::ComputeRadiationLength: final fRadLen="<< This->fRadlen <<" radinv="<< radinv << std::endl;
#endif
#endif
}

FQUALIFIER
void GPMaterial_AddElement(GPMaterial *This,
			   GPElement &element, G4double fraction)
{
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
  std::cout<<"GPMaterial_AddElement: fNumberOfElements="<< This->fNumberOfElements << std::endl;
#endif
#endif
  if(fraction < 0.0 || fraction > 1.0) {
    assert(0); 
  }
  else {
    //check the total fraction
    G4double totalFraction = 0.0;
    for(int i = 0 ; i < This->fNumberOfElements ; i++) {
      totalFraction += This->fMassFractionVector[This->fNumberOfElements];     
    }
    if((totalFraction+fraction)>1.0) {
      assert(0); 
    }
  }

  G4double Zi, Ai;

  if (This->fNumberOfElements < kElementMax ) {

    This->fElementVector[This->fNumberOfElements] = element;
    This->fMassFractionVector[This->fNumberOfElements] = fraction;

    Zi = GPElement_GetZ(&element);
    Ai = GPElement_GetA(&element);

    This->fZ += Zi*fraction;
    This->fA += Ai*fraction;

    This->VecNbOfAtomsPerVolume[This->fNumberOfElements] = 
      Avogadro*This->fDensity*This->fMassFractionVector[This->fNumberOfElements]/Ai;
    This->TotNbOfAtomsPerVolume += 
      This->VecNbOfAtomsPerVolume[This->fNumberOfElements];
    This->TotNbOfElectPerVolume += 
      This->VecNbOfAtomsPerVolume[This->fNumberOfElements]*Zi;

    ++(This->fNumberOfElements);

#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
    std::cout<<"GPMaterial: AddElement(elem,fraction): "<< Zi <<' '<< Ai <<' '<< This->VecNbOfAtomsPerVolume[This->fNumberOfElements]
	     <<' '<< This->TotNbOfAtomsPerVolume <<' '<< This->TotNbOfElectPerVolume << std::endl;
#endif
#endif
  }
  else {
    assert(0);
  }    

  //update radiation length
  GPMaterial_ComputeRadiationLength(This);
}

FQUALIFIER
GPElement* GPMaterial_GetElementVector(GPMaterial *This) 
{
  return &(This->fElementVector[0]);
}

FQUALIFIER
G4double* GPMaterial_GetVecNbOfAtomsPerVolume(GPMaterial *This) 
{
  return &(This->VecNbOfAtomsPerVolume[0]);
}

FQUALIFIER
G4double* GPMaterial_GetMassFractionVector(GPMaterial *This) 
{
  return &(This->fMassFractionVector[0]);
}

FQUALIFIER
G4double GPMaterial_GetTotNbOfAtomsPerVolume(GPMaterial *This)
{
  return This->TotNbOfAtomsPerVolume;
}

FQUALIFIER
G4double GPMaterial_GetTotNbOfElectPerVolume(GPMaterial *This)
{
  return This->TotNbOfElectPerVolume;
}

FQUALIFIER
G4double GPMaterial_GetElectronDensity(GPMaterial *This)
{
  return This->TotNbOfElectPerVolume;
}
