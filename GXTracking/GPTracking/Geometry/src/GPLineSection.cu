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
// $Id: GPLineSection.cc,v 1.10 2006-06-29 18:24:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------

#include "GPLineSection.h" 

FQUALIFIER
void GPLineSection_Constructor( GPLineSection *This,
				const GPThreeVector& PntA, 
				const GPThreeVector& PntB )
{ 
  This->EndpointA = GPThreeVector_mult(PntA,1.0);
  This->VecAtoB = GPThreeVector_sub(PntB,PntA);
  This->fABdistanceSq = GPThreeVector_mag2(This->VecAtoB) ;  
}

FQUALIFIER
G4double GPLineSection_Dist(GPLineSection *This,
			    GPThreeVector OtherPnt )
{
  G4double       dist_sq;  
  GPThreeVector  VecAZ;
  G4double sq_VecAZ, inner_prod, unit_projection ; 

  VecAZ= GPThreeVector_sub(OtherPnt,This->EndpointA);
  sq_VecAZ = GPThreeVector_mag2(VecAZ);

  inner_prod= GPThreeVector_dot(This->VecAtoB,VecAZ );
   
  //  Determine  Projection(AZ on AB) / Length(AB) 
  //
  if( This->fABdistanceSq != 0.0 )
  {
    //  unit_projection= inner_prod * InvsqDistAB();
    unit_projection = inner_prod/This->fABdistanceSq;

    if( (0. <= unit_projection ) && (unit_projection <= 1.0 ) )
    {
      dist_sq= sq_VecAZ -  unit_projection * inner_prod ;
    }
    else
    {
     //  The perpendicular from the point to the line AB meets the line
     //   in a point outside the line segment!
     
      if( unit_projection < 0. ) // A is the closest point
      {
        dist_sq= sq_VecAZ;  
      }
      else                       // B is the closest point
      {
	GPThreeVector   EndpointB = GPThreeVector_add(This->EndpointA,
						      This->VecAtoB);
        GPThreeVector   VecBZ =     GPThreeVector_sub(OtherPnt,
						      EndpointB);
        dist_sq =  GPThreeVector_mag2(VecBZ);
      }
    }
  }
  else
  {
    dist_sq = GPThreeVector_mag2(GPThreeVector_sub(OtherPnt,
						   This->EndpointA)) ;   
  }  
  if( dist_sq < 0.0 ) dist_sq = 0.0 ;

  return sqrt(dist_sq) ;  
}

FQUALIFIER
G4double GPLineSection_Distline( const GPThreeVector& OtherPnt, 
				 const GPThreeVector& LinePntA, 
				 const GPThreeVector& LinePntB )
{
  //GPLineSection( LinePntA, LinePntB );  // Line from A to B
  GPLineSection LineAB;
  GPLineSection_Constructor(&LineAB, LinePntA, LinePntB);

  //  return LineAB.Dist( OtherPnt );
  return GPLineSection_Dist(&LineAB, OtherPnt );
}
