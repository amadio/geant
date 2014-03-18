#include "GXMagneticField.h"

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

  //  double rho = std::sqrt(point[0]*point[0]+point[1]*point[1]); 

  //magnetic field map array 
  //  const int nbinR = 900+1;
  const int nbinZ = 2*1600+1;
  const int noffZ = 1600;

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
    G4double dr = rho - ir;
  
    int iz = int(z + noffZ);
    G4double dz = z + noffZ - iz;
    
    int index = iz+ir*nbinZ;
    
    G4double Bz_lb = This->fieldMap[index].Bz;
    G4double Bz_ub = This->fieldMap[index+nbinZ].Bz;
    
    bField[2] = tesla*(Bz_lb + (Bz_ub-Bz_lb)*dz);  
    
    if(rho>0) {
      G4double Br_lb = This->fieldMap[index].Br;
      G4double Br_ub = This->fieldMap[index+1].Br;
      G4double tmp = tesla*(Br_lb + (Br_ub-Br_lb)*dr)/rho;
      bField[0] = tmp*x;
      bField[1] = tmp*y;
    }
  }
  //use cuda texture memory for fieldMap

}

