#ifndef ROOT_TGeoMatrix_PLAIN
#define ROOT_TGeoMatrix_PLAIN
/*************************************************************************
 * Geometrical transformations. TGeoMatrix - base class, TGeoTranslation *
 * TGeoRotation, TGeoScale, TGeoCombiTrans, TGeoGenTrans .               *
 *************************************************************************/

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

#include <cmath>

//--- globals 
const double kNullVector[3]       =       {0.0,  0.0,  0.0};

const double kIdentityMatrix[3*3] =       {1.0,  0.0,  0.0,
                                             0.0,  1.0,  0.0,
                                             0.0,  0.0,  1.0};

const double kUnitScale[3]        =       {1.0,  1.0,  1.0};

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoMatrix - base class for geometrical transformations.               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

class TGeoMatrix : public TNamed
{
public:
enum EGeoTransfTypes {
   kGeoIdentity  = 0,
   kGeoTranslation  = BIT(17),
   kGeoRotation     = BIT(18),
   kGeoScale        = BIT(19),
   kGeoReflection   = BIT(20),
   kGeoRegistered   = BIT(21),
   kGeoSavePrimitive = BIT(22),
   kGeoMatrixOwned   = BIT(23),
   kGeoCombiTrans   = kGeoTranslation | kGeoRotation,
   kGeoGenTrans     = kGeoTranslation | kGeoRotation | kGeoScale
};

protected:
   TGeoMatrix(const TGeoMatrix &other);

public :
   TGeoMatrix();

   // a pure translation
   TGeoMatrix(double dx, double dy, double dz)
     {
       fTranslation[0]=dx;
       fTranslation[1]=dy;
       fTranslation[2]=dz;
       fRotationMatrix[0] = 1.;
       fRotationMatrix[1] = 0.;
       fRotationMatrix[2] = 0.; 
       fRotationMatrix[3] = 0.;
       fRotationMatrix[4] = 1.;
       fRotationMatrix[5] = 0.;
       fRotationMatrix[6] = 0.;
       fRotationMatrix[7] = 0.;
       fRotationMatrix[8] = 1.;
     }
   TGeoMatrix(double dx, double dy, double dz, double phi, double theta, double psi)
     {
       fTranslation[0]=dx;
       fTranslation[1]=dy;
       fTranslation[2]=dz;
       SetAngles(phi,theta,psi);
     }
   // translation and Euler angles
   ~TGeoMatrix();

   TGeoMatrix& operator=(const TGeoMatrix &matrix);
   bool      operator ==(const TGeoMatrix &other) const;
   
   bool               IsIdentity()    const {return !TestBit(kGeoGenTrans);}
   bool               IsTranslation() const {return TestBit(kGeoTranslation);}
   bool               IsRotation()    const {return TestBit(kGeoRotation);}
   bool               IsReflection()  const {return TestBit(kGeoReflection);}
   bool               IsScale()       const {return TestBit(kGeoScale);}
   bool               IsCombi()       const {return (TestBit(kGeoTranslation) 
                                               && TestBit(kGeoRotation));}
   bool               IsGeneral()     const {return (TestBit(kGeoTranslation) 
                            && TestBit(kGeoRotation) && TestBit(kGeoScale));}
   bool               IsRegistered()  const {return TestBit(kGeoRegistered);}
   bool               IsRotAboutZ()   const;
   void                 GetHomogenousMatrix(double *hmat) const;
   char                *GetPointerName() const;

   int              GetByteCount() const;
   double const    *GetTranslation()    const { return fTranslation; }
   double const    *GetRotationMatrix() const { return fRotationMatrix; }
   inline void MasterToLocal(const double *master, double *local) const
   {
     // convert a point by multiplying its column vector (x, y, z, 1) to matrix
     // if idendity?? 

     // do the translation
     double const *tr  = GetTranslation();
     double mt0  = master[0]-tr[0];
     double mt1  = master[1]-tr[1];
     double mt2  = master[2]-tr[2];

     // do the rotation
     double const *rot = GetRotationMatrix();
     local[0] = mt0*rot[0] + mt1*rot[3] + mt2*rot[6];
     local[1] = mt0*rot[1] + mt1*rot[4] + mt2*rot[7];
     local[2] = mt0*rot[2] + mt1*rot[5] + mt2*rot[8];
   };

 private:
   // data members
   double             fTranslation[3];  // translation vector
   double             fRotationMatrix[3*3];   // rotation matrix

//_____________________________________________________________________________
   void SetAngles(double phi, double theta, double psi)
   {
     // Set matrix elements according to Euler angles
     double degrad = M_PI/180.;
     double sinphi = sin(degrad*phi);
     double cosphi = cos(degrad*phi);
     double sinthe = sin(degrad*theta);
     double costhe = cos(degrad*theta);
     double sinpsi = sin(degrad*psi);
     double cospsi = cos(degrad*psi);
     
     fRotationMatrix[0] =  cospsi*cosphi - costhe*sinphi*sinpsi;
     fRotationMatrix[1] = -sinpsi*cosphi - costhe*sinphi*cospsi;
     fRotationMatrix[2] =  sinthe*sinphi;
     fRotationMatrix[3] =  cospsi*sinphi + costhe*cosphi*sinpsi;
     fRotationMatrix[4] = -sinpsi*sinphi + costhe*cosphi*cospsi;
     fRotationMatrix[5] = -sinthe*cosphi;
     fRotationMatrix[6] =  sinpsi*sinthe;
     fRotationMatrix[7] =  cospsi*sinthe;
     fRotationMatrix[8] =  costhe;
   }
};

#endif
