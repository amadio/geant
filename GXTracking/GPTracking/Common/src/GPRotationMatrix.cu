#include "GPRotationMatrix.h"

FQUALIFIER
GPRotationMatrix GPRotationMatrix_create(G4double mxx, G4double mxy, G4double mxz,
					 G4double myx, G4double myy, G4double myz,
					 G4double mzx, G4double mzy, G4double mzz)
{
  GPRotationMatrix r = { mxx, mxy, mxz, 
			 myx, myy, myz, 
			 mzx, mzy, mzz};
  return r;
}

FQUALIFIER
GPThreeVector GPRotationMatrix_apply(const GPRotationMatrix *This, 
				     GPThreeVector p)
{
  return GPThreeVector_create(This->rxx*p.x + This->rxy*p.y + This->rxz*p.z,
			      This->ryx*p.x + This->ryy*p.y + This->ryz*p.z,
			      This->rzx*p.x + This->rzy*p.y + This->rzz*p.z);
}

FQUALIFIER
GPRotationMatrix GPRotationMatrix_mult(const GPRotationMatrix *This, 
				       const GPRotationMatrix *m)
{
  return GPRotationMatrix_create(This->rxx*m->rxx + This->rxy*m->ryx + This->rxz*m->rzx,
				 This->rxx*m->rxy + This->rxy*m->ryy + This->rxz*m->rzy,
				 This->rxx*m->rxz + This->rxy*m->ryz + This->rxz*m->rzz,
				 This->ryx*m->rxx + This->ryy*m->ryx + This->ryz*m->rzx,
				 This->ryx*m->rxy + This->ryy*m->ryy + This->ryz*m->rzy,
				 This->ryx*m->rxz + This->ryy*m->ryz + This->ryz*m->rzz,
				 This->rzx*m->rxx + This->rzy*m->ryx + This->rzz*m->rzx,
				 This->rzx*m->rxy + This->rzy*m->ryy + This->rzz*m->rzy,
				 This->rzx*m->rxz + This->rzy*m->ryz + This->rzz*m->rzz);
}

FQUALIFIER
GPRotationMatrix GPRotationMatrix_transform(GPRotationMatrix *This, 
					    const GPRotationMatrix *m)
{
  *This = GPRotationMatrix_mult(m,This);
  return *This;
}

FQUALIFIER
GPRotationMatrix GPRotationMatrix_inverse(const GPRotationMatrix *This)
{
  return GPRotationMatrix_create(This->rxx, This->ryx, This->rzx, 
				 This->rxy, This->ryy, This->rzy, 
				 This->rxz, This->ryz, This->rzz );
}

FQUALIFIER
GPRotationMatrix GPRotationMatrix_invert(GPRotationMatrix *This)
{
  return *This = GPRotationMatrix_inverse(This);
}

FQUALIFIER
GPRotationMatrix GPRotationMatrix_rotateAxes(GPRotationMatrix *This,
					     GPThreeVector newX,
					     GPThreeVector newY,
					     GPThreeVector newZ)
{
  GPRotationMatrix m = GPRotationMatrix_create(newX.x, newY.x, newZ.x,
					       newX.y, newY.y, newZ.y,
					       newX.z, newY.z, newZ.z);
  return GPRotationMatrix_transform(This,&m);
}

FQUALIFIER
void GPRotationMatrix_rotateX(GPRotationMatrix *This, 
			      G4double a)
{
  G4double c = cos(a);
  G4double s = sin(a);
  G4double x = This->ryx; 
  G4double y = This->ryy;
  G4double z = This->ryz; 

  This->ryx = c*x - s*This->rzx;
  This->ryy = c*y - s*This->rzy;
  This->ryz = c*z - s*This->rzz;
  This->rzx = s*x + c*This->rzx;
  This->rzy = s*y + c*This->rzy;
  This->rzz = s*z + c*This->rzz;
}

FQUALIFIER
void GPRotationMatrix_rotateY(GPRotationMatrix *This, 
			      G4double a)
{
  G4double c = cos(a);
  G4double s = sin(a);
  G4double x = This->rzx;
  G4double y = This->rzy;
  G4double z = This->rzz; 

  This->rzx = c*x - s*This->rxx;
  This->rzy = c*y - s*This->rxy;
  This->rzz = c*z - s*This->rxz;
  This->rxx = s*x + c*This->rxx;
  This->rxy = s*y + c*This->rxy;
  This->rxz = s*z + c*This->rxz;
}

FQUALIFIER
void GPRotationMatrix_rotateZ(GPRotationMatrix *This, G4double a)
{
  G4double c = cos(a);
  G4double s = sin(a);
  G4double x = This->rxx;
  G4double y = This->rxy;
  G4double z = This->rxz; 

  This->rxx = c*x - s*This->ryx;
  This->rxy = c*y - s*This->ryy;
  This->rxz = c*z - s*This->ryz;
  This->ryx = s*x + c*This->ryx;
  This->ryy = s*y + c*This->ryy;
  This->ryz = s*z + c*This->ryz;
}
