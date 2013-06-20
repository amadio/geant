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
const Double_t kNullVector[3]       =       {0.0,  0.0,  0.0};

const Double_t kIdentityMatrix[3*3] =       {1.0,  0.0,  0.0,
                                             0.0,  1.0,  0.0,
                                             0.0,  0.0,  1.0};

const Double_t kUnitScale[3]        =       {1.0,  1.0,  1.0};

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
   TGeoMatrix(Double_t dx, Double_t dy, Double_t dz)
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
   TGeoMatrix(Double_t dx, Double_t dy, Double_t dz, Double_t phi, Double_t theta, Double_t psi)
     {
       fTranslation[0]=dx;
       fTranslation[1]=dy;
       fTranslation[2]=dz;
       SetAngles(phi,theta,psi);
     }
   // translation and Euler angles
   ~TGeoMatrix();

   TGeoMatrix& operator=(const TGeoMatrix &matrix);
   Bool_t      operator ==(const TGeoMatrix &other) const;
   
   Bool_t               IsIdentity()    const {return !TestBit(kGeoGenTrans);}
   Bool_t               IsTranslation() const {return TestBit(kGeoTranslation);}
   Bool_t               IsRotation()    const {return TestBit(kGeoRotation);}
   Bool_t               IsReflection()  const {return TestBit(kGeoReflection);}
   Bool_t               IsScale()       const {return TestBit(kGeoScale);}
   Bool_t               IsCombi()       const {return (TestBit(kGeoTranslation) 
                                               && TestBit(kGeoRotation));}
   Bool_t               IsGeneral()     const {return (TestBit(kGeoTranslation) 
                            && TestBit(kGeoRotation) && TestBit(kGeoScale));}
   Bool_t               IsRegistered()  const {return TestBit(kGeoRegistered);}
   Bool_t               IsRotAboutZ()   const;
   void                 GetHomogenousMatrix(Double_t *hmat) const;
   char                *GetPointerName() const;

   Int_t              GetByteCount() const;
   Double_t const    *GetTranslation()    const { return fTranslation; }
   Double_t const    *GetRotationMatrix() const { return fRotationMatrix; }
   inline void MasterToLocal(const Double_t *master, Double_t *local) const
   {
     // convert a point by multiplying its column vector (x, y, z, 1) to matrix
     // if idendity?? 

     // do the translation
     Double_t const *tr  = GetTranslation();
     Double_t mt0  = master[0]-tr[0];
     Double_t mt1  = master[1]-tr[1];
     Double_t mt2  = master[2]-tr[2];

     // do the rotation
     Double_t const *rot = GetRotationMatrix();
     local[0] = mt0*rot[0] + mt1*rot[3] + mt2*rot[6];
     local[1] = mt0*rot[1] + mt1*rot[4] + mt2*rot[7];
     local[2] = mt0*rot[2] + mt1*rot[5] + mt2*rot[8];
   };

 private:
   // data members
   Double_t             fTranslation[3];  // translation vector
   Double_t             fRotationMatrix[3*3];   // rotation matrix

//_____________________________________________________________________________
   void SetAngles(Double_t phi, Double_t theta, Double_t psi)
   {
     // Set matrix elements according to Euler angles
     Double_t degrad = M_PI/180.;
     Double_t sinphi = sin(degrad*phi);
     Double_t cosphi = cos(degrad*phi);
     Double_t sinthe = sin(degrad*theta);
     Double_t costhe = cos(degrad*theta);
     Double_t sinpsi = sin(degrad*psi);
     Double_t cospsi = cos(degrad*psi);
     
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
