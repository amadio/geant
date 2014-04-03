#ifndef GPROTATIONMATRIX_HH
#define GPROTATIONMATRIX_HH

#include "GPTypeDef.h"
#include "GPThreeVector.h"

struct GPRotationMatrix
{
  G4double
  rxx, rxy, rxz, 
  ryx, ryy, ryz, 
  rzx, rzy, rzz;
};

extern "C" {

FQUALIFIER
GPRotationMatrix GPRotationMatrix_create(G4double xx, G4double xy, G4double xz,
                                         G4double yx, G4double yy, G4double yz,
                                         G4double zx, G4double zy, G4double zz);
FQUALIFIER
GPThreeVector GPRotationMatrix_apply(const GPRotationMatrix *This, 
                                     GPThreeVector p);

FQUALIFIER
GPRotationMatrix GPRotationMatrix_mult(const GPRotationMatrix *This, 
                                       const GPRotationMatrix *m);

FQUALIFIER
GPRotationMatrix GPRotationMatrix_transform(GPRotationMatrix *This, 
                                            const GPRotationMatrix *m);

FQUALIFIER
GPRotationMatrix GPRotationMatrix_inverse(const GPRotationMatrix *This);

FQUALIFIER
GPRotationMatrix GPRotationMatrix_invert(GPRotationMatrix *This);

FQUALIFIER
GPRotationMatrix GPRotationMatrix_rotateAxes(GPRotationMatrix *This,
                                             GPThreeVector newX,
                                             GPThreeVector newY,
                                             GPThreeVector newZ);
FQUALIFIER
void GPRotationMatrix_rotateX(GPRotationMatrix *This, G4double a);

FQUALIFIER
void GPRotationMatrix_rotateY(GPRotationMatrix *This, G4double a);

FQUALIFIER
void GPRotationMatrix_rotateZ(GPRotationMatrix *This, G4double a);

}

#endif
