#include "GXMagneticField.h"
#include "GPConstants.h"
#include "GXMagneticFieldTexture.h"

FQUALIFIER 
void GXMagneticField_Constructor(GXMagneticField *This,
				 GXFieldMap *fieldMapArray)
{
  This->fieldMap = fieldMapArray;
}

FQUALIFIER 
void GXMagneticField_GetFieldValue( GXMagneticField *This,
				    const  G4double point[3],
				    G4double *bField ) 
{
  //default values
  bField[0] = bField[1] = bField[2] = 0.;

  //magnetic field unit
  const G4double tesla = 0.001;

  G4double x        = point[0];        
  G4double y        = point[1];    

  //conversion: [mm] to [cm]
  G4double rho      = 0.1*sqrt(x*x+y*y);
  G4double z        = 0.1*point[2]; 

  if(abs(z)>1600.0 || rho > 900.0) {
    ; // leave bField at default values
  }
  else {
    int ir = int(rho);
    int iz = int(z + noffZ);

#ifdef __CUDA_ARCH__
    bField[2] = tesla*tex2D(texZ, z+0.5+noffZ, ir);
#else
    G4double lbv, ubv;
    G4double dr = rho - ir;
    G4double dz = z + noffZ - iz;

    int index = iz+ir*nbinZ;

    lbv = This->fieldMap[index].Bz;
    ubv = This->fieldMap[index+nbinZ].Bz;
    
    bField[2] = tesla * (lbv + (ubv-lbv)*dz);
#endif

    if(rho>0) {
      G4double tmp;
#ifdef  __CUDA_ARCH__
      G4double lookup_rho = rho;
      if (iz > noffZ ) {
        iz -= noffZ;
        lookup_rho += nbinR;
      }
      tmp = tesla * tex2D(texR, lookup_rho+0.5, iz) / rho;
#else
      lbv = This->fieldMap[index].Br;
      ubv = This->fieldMap[index+1].Br;
      tmp = tesla * (lbv + (ubv-lbv)*dr) / rho;
#endif
      bField[0] = tmp*x;
      bField[1] = tmp*y;
    }
  }
}

