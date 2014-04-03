//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4AffineTransform.hh,v 1.6 2006-06-29 18:30:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4AffineTransform
//
// Class description:
//
// A class for geometric affine transformations [see, eg. Foley & Van Dam]
// Supports efficient arbitrary rotation & transformation of vectors and the
// computation of compound & inverse transformations. A `rotation flag' is
// maintained internally for greater computational efficiency for transforms
// that do not involve rotation.
//
// Interfaces to the CLHEP classes G4ThreeVector & G4RotationMatrix
//
// For member function descriptions, see comments by declarations. For
// additional clarification, also check the `const' declarations for
// functions & their parameters.
//
// Member data:
//
//      G4double rxx,rxy,rxz; 
//      G4double ryx,ryy,ryz;  A 3x3 rotation matrix - net rotation
//      G4double rzx,rzy,rzz;
//      G4double tx,ty,tz;     Net translation 

// History:
// Paul R C Kent 6 Aug 1996 - initial version
//
// 19.09.96 E.Chernyaev:
// - direct access to the protected members of the G4RotationMatrix class
//   replaced by access via public access functions            
// - conversion of the rotation matrix to angle & axis used to get
//   a possibility to remove "friend" from the G4RotationMatrix class
// --------------------------------------------------------------------
#ifndef GPAFFINETRANSFORM_HH
#define GPAFFINETRANSFORM_HH

#include "GPTypeDef.h"
#include "GPThreeVector.h"
#include "GPRotationMatrix.h"

struct GPAffineTransform
{
  G4double rxx,rxy,rxz;
  G4double ryx,ryy,ryz;
  G4double rzx,rzy,rzz;
  G4double tx,ty,tz;
};

extern "C" {

FQUALIFIER 
void GPAffineTransform_Constructor( GPAffineTransform *This );

FQUALIFIER
GPAffineTransform GPAffineTransform_create_id(void);

FQUALIFIER 
void GPAffineTransform_Constructor2( GPAffineTransform *This,
                                     const GPRotationMatrix rot,
                                     const GPThreeVector tlate );
FQUALIFIER 
void GPAffineTransform_Constructor3( GPAffineTransform *This, 
				     const GPRotationMatrix *rot,
				     const GPThreeVector tlate);

FQUALIFIER 
GPAffineTransform GPAffineTransform_create_id(void);

FQUALIFIER 
void GPAffineTransform_Set(GPAffineTransform *This,
                           GPAffineTransform rhs);

FQUALIFIER 
void GPAffineTransform_Vector(GPAffineTransform *This,
                              const GPThreeVector tlate);

FQUALIFIER 
void GPAffineTransform_Matrix(GPAffineTransform *This,
                              const GPRotationMatrix rot);

FQUALIFIER 
void GPAffineTransform_Elements( GPAffineTransform *This, 
                  const G4double prxx,const G4double prxy,const G4double prxz,
                  const G4double pryx,const G4double pryy,const G4double pryz,
                  const G4double przx,const G4double przy,const G4double przz,
		  const G4double ptx,const G4double pty,const G4double ptz);

FQUALIFIER 
GPAffineTransform GPAffineTransform_Product(GPAffineTransform *This,
                                            const GPAffineTransform *tf1,
                                            const GPAffineTransform *tf2);

FQUALIFIER 
GPAffineTransform GPAffineTransform_InverseProduct( GPAffineTransform *This,
                                                    GPAffineTransform *tf1,
                                                    GPAffineTransform *tf2);

FQUALIFIER
GPThreeVector GPAffineTransform_TransformPoint(const GPAffineTransform *This, 
                                               const GPThreeVector vec);

FQUALIFIER
GPThreeVector GPAffineTransform_TransformAxis(const GPAffineTransform *This,
                                              const GPThreeVector axis);

FQUALIFIER 
GPThreeVector GPAffineTransform_ApplyPointTransform(const GPAffineTransform *This,
                                                    GPThreeVector vec);

FQUALIFIER
GPThreeVector GPAffineTransform_ApplyAxisTransform(const GPAffineTransform *This,
                                                   GPThreeVector axis);

FQUALIFIER
GPAffineTransform GPAffineTransform_Inverse(GPAffineTransform *This);

FQUALIFIER
GPAffineTransform GPAffineTransform_Invert(GPAffineTransform *This);

FQUALIFIER
G4bool GPAffineTransform_IsRotated( GPAffineTransform *This );

FQUALIFIER 
G4bool GPAffineTransform_IsTranslated(GPAffineTransform *This);

FQUALIFIER 
GPRotationMatrix GPAffineTransform_NetRotation( GPAffineTransform *This );

FQUALIFIER
GPThreeVector GPAffineTransform_NetTranslation(GPAffineTransform *This);

FQUALIFIER 
void GPAffineTransform_SetNetRotation(GPAffineTransform *This,
                                      const GPRotationMatrix rot);
FQUALIFIER
void GPAffineTransform_SetNetTranslation(GPAffineTransform *This,
                                         const GPThreeVector tlate);


}

#endif
