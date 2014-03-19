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
// $Id: GPVoxelLimits.cc,v 1.11 2006-06-29 18:34:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class GPVoxelLimits
// 
// Implementation
//
// History:
//
// 14.03.02 V. Grichine, cosmetics
// 13.07.95 P.Kent Initial version
// --------------------------------------------------------------------

#include "GPVoxelLimits.h"

//#include "G4ios.hh"

///////////////////////////////////////////////////////////////////////////
//
// Empty constructor and destructor
//

FQUALIFIER
void GPVoxelLimits_Constructor(GPVoxelLimits *This)
{
  This->fxAxisMin = -kInfinity;
  This->fxAxisMax =  kInfinity;
  This->fyAxisMin = -kInfinity;
  This->fyAxisMax =  kInfinity;
  This->fzAxisMin = -kInfinity;
  This->fzAxisMax =  kInfinity;
}

/*
GPVoxelLimits_~GPVoxelLimits()
{
}
*/

///////////////////////////////////////////////////////////////////////////
//
// Further restrict limits
// No checks for illegal restrictions
//

FQUALIFIER
void GPVoxelLimits_AddLimit( GPVoxelLimits *This,
			     const EAxis pAxis, 
			     const G4double pMin,
			     const G4double pMax )
{
  if ( pAxis == kXAxis )
  {
    if ( pMin > This->fxAxisMin ) This->fxAxisMin = pMin ;    
    if ( pMax < This->fxAxisMax ) This->fxAxisMax = pMax ;    
  }
  else if ( pAxis == kYAxis )
  {
    if ( pMin > This->fyAxisMin ) This->fyAxisMin = pMin ;    
    if ( pMax < This->fyAxisMax ) This->fyAxisMax = pMax ;
  }
  else
  { 
    //    assert( pAxis == kZAxis ) ;
    if ( pMin > This->fzAxisMin ) This->fzAxisMin = pMin ;
    if ( pMax < This->fzAxisMax ) This->fzAxisMax = pMax ;
  }
}

///////////////////////////////////////////////////////////////////////////
//
// ClipToLimits
//
// Clip the line segment pStart->pEnd to the volume described by the
// current limits. Return true if the line remains after clipping,
// else false, and leave the vectors in an undefined state.
//
// Process:
//
// Use Cohen-Sutherland clipping in 3D
// [Fundamentals of Interactive Computer Graphics,Foley & Van Dam]
//

FQUALIFIER
G4bool GPVoxelLimits_ClipToLimits( GPVoxelLimits *This, 
				   GPThreeVector* pStart,
				   GPThreeVector* pEnd      )
{
  G4int sCode, eCode ;
  G4bool remainsAfterClip ;
    
  // Determine if line is trivially inside (both outcodes==0) or outside
  // (logical AND of outcodes !=0)

  sCode = GPVoxelLimits_OutCode(This,*pStart) ;
  eCode = GPVoxelLimits_OutCode(This,*pEnd)   ;

  if ( sCode & eCode )
  {
    // Trivially outside, no intersection with region

    remainsAfterClip = false;
  }
  else if ( sCode == 0 && eCode == 0 )
  {
    // Trivially inside, no intersections

    remainsAfterClip = true ;
  }
  else
  {
    // Line segment *may* cut volume boundaries
    // At most, one end point is inside

    G4double x1, y1, z1, x2, y2, z2 ;

    /*
    x1 = GPThreeVector_x(pStart) ;
    y1 = GPThreeVector_y(pStart) ;
    z1 = GPThreeVector_z(pStart) ;

    x2 = GPThreeVector_x(pEnd) ;
    y2 = GPThreeVector_y(pEnd) ;
    z2 = GPThreeVector_z(pEnd) ;
    */
    x1 = pStart->x ;
    y1 = pStart->y ;
    z1 = pStart->z ;

    x2 = pEnd->x ;
    y2 = pEnd->y ;
    z2 = pEnd->z ;

    /*
    if( std::abs(x1-x2) < kCarTolerance*kCarTolerance)
    {
      G4cout<<"x1 = "<<x1<<"\t"<<"x2 = "<<x2<<G4endl; 
    }   
    if( std::abs(y1-y2) < kCarTolerance*kCarTolerance)
    {
      G4cout<<"y1 = "<<y1<<"\t"<<"y2 = "<<y2<<G4endl; 
    }   
    if( std::abs(z1-z2) < kCarTolerance*kCarTolerance)
    {
      G4cout<<"z1 = "<<z1<<"\t"<<"z2 = "<<z2<<G4endl; 
    } 
    */  
    while ( sCode != eCode )
    {
      // Copy vectors to work variables x1-z1,x2-z2
      // Ensure x1-z1 lies outside volume, swapping vectors and outcodes
      // if necessary

      if ( sCode )
      {
        if ( sCode & 0x01 )  // Clip against This->fxAxisMin
        {
          z1 += (This->fxAxisMin-x1)*(z2-z1)/(x2-x1);
          y1 += (This->fxAxisMin-x1)*(y2-y1)/(x2-x1);
          x1  = This->fxAxisMin;
        }
        else if ( sCode & 0x02 ) // Clip against This->fxAxisMax
        {
          z1 += (This->fxAxisMax-x1)*(z2-z1)/(x2-x1);
          y1 += (This->fxAxisMax-x1)*(y2-y1)/(x2-x1);
          x1  = This->fxAxisMax ;
        }
        else if ( sCode & 0x04 )  // Clip against This->fyAxisMin
        {
          x1 += (This->fyAxisMin-y1)*(x2-x1)/(y2-y1);
          z1 += (This->fyAxisMin-y1)*(z2-z1)/(y2-y1);
          y1  = This->fyAxisMin;
        }
        else if ( sCode & 0x08 )  // Clip against This->fyAxisMax
        {
          x1 += (This->fyAxisMax-y1)*(x2-x1)/(y2-y1);
          z1 += (This->fyAxisMax-y1)*(z2-z1)/(y2-y1);
          y1  = This->fyAxisMax;
        }
        else if ( sCode & 0x10 )  // Clip against This->fzAxisMin
        {
          x1 += (This->fzAxisMin-z1)*(x2-x1)/(z2-z1);
          y1 += (This->fzAxisMin-z1)*(y2-y1)/(z2-z1);
          z1  = This->fzAxisMin;
        }
        else if ( sCode & 0x20 )  // Clip against This->fzAxisMax
        {
          x1 += (This->fzAxisMax-z1)*(x2-x1)/(z2-z1);
          y1 += (This->fzAxisMax-z1)*(y2-y1)/(z2-z1);
          z1  = This->fzAxisMax;
        }
      }
      if ( eCode )  // Clip 2nd end: repeat of 1st, but 1<>2
      {
        if ( eCode & 0x01 )  // Clip against This->fxAxisMin
        {
          z2 += (This->fxAxisMin-x2)*(z1-z2)/(x1-x2);
          y2 += (This->fxAxisMin-x2)*(y1-y2)/(x1-x2);
          x2  = This->fxAxisMin;
        }
        else if ( eCode & 0x02 )  // Clip against This->fxAxisMax
        {
          z2 += (This->fxAxisMax-x2)*(z1-z2)/(x1-x2);
          y2 += (This->fxAxisMax-x2)*(y1-y2)/(x1-x2);
          x2  = This->fxAxisMax;
        }
        else if ( eCode & 0x04 )  // Clip against This->fyAxisMin
        {
          x2 += (This->fyAxisMin-y2)*(x1-x2)/(y1-y2);
          z2 += (This->fyAxisMin-y2)*(z1-z2)/(y1-y2);
          y2  = This->fyAxisMin;
        }
        else if (eCode&0x08)  // Clip against This->fyAxisMax
        {
          x2 += (This->fyAxisMax-y2)*(x1-x2)/(y1-y2);
          z2 += (This->fyAxisMax-y2)*(z1-z2)/(y1-y2);
          y2  = This->fyAxisMax;
        }
        else if ( eCode & 0x10 )  // Clip against This->fzAxisMin
        {
          x2 += (This->fzAxisMin-z2)*(x1-x2)/(z1-z2);
          y2 += (This->fzAxisMin-z2)*(y1-y2)/(z1-z2);
          z2  = This->fzAxisMin;
        }
        else if ( eCode & 0x20 )  // Clip against This->fzAxisMax
        {
          x2 += (This->fzAxisMax-z2)*(x1-x2)/(z1-z2);
          y2 += (This->fzAxisMax-z2)*(y1-y2)/(z1-z2);
          z2  = This->fzAxisMax;
        }
      }
      //  G4endl; G4cout<<"x1 = "<<x1<<"\t"<<"x2 = "<<x2<<G4endl<<G4endl;
      *pStart = GPThreeVector_create(x1,y1,z1);
      *pEnd   = GPThreeVector_create(x2,y2,z2);

      sCode  = GPVoxelLimits_OutCode(This,*pStart);
      eCode  = GPVoxelLimits_OutCode(This,*pEnd);
    }
    if ( sCode == 0 && eCode == 0 ) remainsAfterClip = true;
    else                            remainsAfterClip = false;
  }
  return remainsAfterClip;
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate the `outcode' for the specified vector:
// The following bits are set:
//   0      pVec.x()<This->fxAxisMin && IsXLimited()
//   1      pVec.x()>This->fxAxisMax && IsXLimited()
//   2      pVec.y()<This->fyAxisMin && IsYLimited()
//   3      pVec.y()>This->fyAxisMax && IsYLimited()
//   4      pVec.z()<This->fzAxisMin && IsZLimited()
//   5      pVec.z()>This->fzAxisMax && IsZLimited()
//

FQUALIFIER
G4int GPVoxelLimits_OutCode( GPVoxelLimits *This,
			     GPThreeVector pVec )
{
  G4int code = 0 ;                // The outcode

  if ( GPVoxelLimits_IsXLimited(This) )
  {
    if ( GPThreeVector_x(pVec) < This->fxAxisMin ) code |= 0x01 ;
    if ( GPThreeVector_x(pVec) > This->fxAxisMax ) code |= 0x02 ;
  }
  if ( GPVoxelLimits_IsYLimited(This) )
  {
    if ( GPThreeVector_y(pVec) < This->fyAxisMin ) code |= 0x04 ;
    if ( GPThreeVector_y(pVec) > This->fyAxisMax ) code |= 0x08 ;
  }
  if ( GPVoxelLimits_IsZLimited(This) )
  {
    if ( GPThreeVector_z(pVec) < This->fzAxisMin ) code |= 0x10 ;
    if ( GPThreeVector_z(pVec) > This->fzAxisMax ) code |= 0x20 ;
  }
  return code;
}

///////////////////////////////////////////////////////////////////////////////
/*
std::ostream& operator << (std::ostream& os, const GPVoxelLimits& pLim)
{
    os << "{";
    if (pLim.IsXLimited())
        {
            os << "(" << pLim.GetMinXExtent() 
               << "," << pLim.GetMaxXExtent() << ") ";
        }
    else
        {
            os << "(-,-) ";
        }
    if (pLim.IsYLimited())
        {
            os << "(" << pLim.GetMinYExtent() 
               << "," << pLim.GetMaxYExtent() << ") ";
        }
    else
        {
            os << "(-,-) ";
        }
    if (pLim.IsZLimited())
        {
            os << "(" << pLim.GetMinZExtent()
               << "," << pLim.GetMaxZExtent() << ")";
        }
    else
        {
            os << "(-,-)";
        }
    os << "}";
    return os;
}
*/

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
// $Id: GPVoxelLimits.icc,v 1.4 2006-06-29 18:33:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GPVoxelLimits __Host__ __Device__ implementation
//
// --------------------------------------------------------------------

FQUALIFIER
G4double GPVoxelLimits_GetMaxXExtent(GPVoxelLimits *This )
{
  return This->fxAxisMax;
}

FQUALIFIER
G4double GPVoxelLimits_GetMaxYExtent(GPVoxelLimits *This )
{
  return This->fyAxisMax;
}

FQUALIFIER
G4double GPVoxelLimits_GetMaxZExtent(GPVoxelLimits *This )
{
  return This->fzAxisMax;
}

FQUALIFIER
G4double GPVoxelLimits_GetMinXExtent(GPVoxelLimits *This )
{
  return This->fxAxisMin;
}

FQUALIFIER
G4double GPVoxelLimits_GetMinYExtent(GPVoxelLimits *This )
{
  return This->fyAxisMin;
}

FQUALIFIER
G4double GPVoxelLimits_GetMinZExtent(GPVoxelLimits *This )
{
  return This->fzAxisMin;
}

FQUALIFIER
G4double GPVoxelLimits_GetMaxExtent(GPVoxelLimits *This ,
				     const EAxis pAxis)
{
  if (pAxis==kXAxis)
  {
    return GPVoxelLimits_GetMaxXExtent(This);
  }
  else if (pAxis==kYAxis)
  {
    return GPVoxelLimits_GetMaxYExtent(This);
  }
  else 
  {
    //    assert(pAxis==kZAxis);
    return GPVoxelLimits_GetMaxZExtent(This);
  }
}

FQUALIFIER
G4double GPVoxelLimits_GetMinExtent( GPVoxelLimits *This ,
				     const EAxis pAxis)
{
  if (pAxis==kXAxis)
  {
    return GPVoxelLimits_GetMinXExtent(This);
  }
  else if (pAxis==kYAxis)
  {
    return GPVoxelLimits_GetMinYExtent(This);
  }
  else 
  {
    //    assert(pAxis==kZAxis);
    return GPVoxelLimits_GetMinZExtent(This);
  }
}

FQUALIFIER
G4bool GPVoxelLimits_IsXLimited(GPVoxelLimits *This )
{
  return (This->fxAxisMin==-kInfinity&&This->fxAxisMax==kInfinity) ? false : true;
}

FQUALIFIER
G4bool GPVoxelLimits_IsYLimited(GPVoxelLimits *This )
{
  return (This->fyAxisMin==-kInfinity&&This->fyAxisMax==kInfinity) ? false : true;
}

FQUALIFIER
G4bool GPVoxelLimits_IsZLimited(GPVoxelLimits *This )
{
  return (This->fzAxisMin==-kInfinity&&This->fzAxisMax==kInfinity) ? false : true;
}

FQUALIFIER
G4bool GPVoxelLimits_IsLimited(GPVoxelLimits *This )
{
  return (GPVoxelLimits_IsXLimited(This)||
	  GPVoxelLimits_IsYLimited(This)||
	  GPVoxelLimits_IsZLimited(This));
}

FQUALIFIER
G4bool GPVoxelLimits_IsLimited2(GPVoxelLimits *This ,
				const EAxis pAxis)
{
  if (pAxis==kXAxis)
  {
    return GPVoxelLimits_IsXLimited(This);
  }
  else if (pAxis==kYAxis)
  {
    return GPVoxelLimits_IsYLimited(This);
  }
  else 
  {
    //    assert(pAxis==kZAxis);
    return GPVoxelLimits_IsZLimited(This);
  }
}

FQUALIFIER
G4bool GPVoxelLimits_Inside( GPVoxelLimits *This , 
			     GPThreeVector pVec)
{
  return ((GPVoxelLimits_GetMinXExtent(This)<=GPThreeVector_x(pVec)) &&
	  (GPVoxelLimits_GetMaxXExtent(This)>=GPThreeVector_x(pVec)) &&
	  (GPVoxelLimits_GetMinYExtent(This)<=GPThreeVector_y(pVec)) &&
	  (GPVoxelLimits_GetMaxYExtent(This)>=GPThreeVector_y(pVec)) &&
	  (GPVoxelLimits_GetMinZExtent(This)<=GPThreeVector_z(pVec)) &&
	  (GPVoxelLimits_GetMaxZExtent(This)>=GPThreeVector_z(pVec)) ) 
    ? true : false;
}
