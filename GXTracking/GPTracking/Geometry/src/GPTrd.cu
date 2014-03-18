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
// $Id: GPTrd.cc,v 1.38 2010-10-19 15:42:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Implementation for GPTrd class
//
// History:
//
// 28.04.05 V.Grichine: new SurfaceNormal according to J. Apostolakis proposal 
// 26.04.05, V.Grichine, new SurfaceNoramal is default
// 07.12.04, V.Grichine, SurfaceNoramal with edges/vertices.
// 07.05.00, V.Grichine, in d = DistanceToIn(p,v), if d<0.5*kCarTolerance, d=0
//    ~1996, V.Grichine, 1st implementation based on old code of P.Kent
//
//////////////////////////////////////////////////////////////////////////////

#include "GPTrd.h"
//#include "GPVSolid.h"

//#include "G4VPVParameterisation.hh"
//#include "GPVoxelLimits.h"
//#include "GPAffineTransform.h"
//#include "Randomize.hh"

//#include "G4VGraphicsScene.hh"
//#include "G4Polyhedron.hh"
//#include "G4NURBS.hh"
//#include "G4NURBSbox.hh"

//#include "GPThreeVectorList.h"
#include "GPUtils.h"

//using namespace CLHEP;

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths

FQUALIFIER
void GPTrd_Constructor(GPTrd *This,
		       G4double pdx1,  G4double pdx2,
		       G4double pdy1,  G4double pdy2,
		       G4double pdz )
{
  This->fSolid.fType = kTrd;
  GPTrd_CheckAndSetAllParameters (This,pdx1, pdx2, pdy1, pdy2, pdz);
}

/////////////////////////////////////////////////////////////////////////
//
// Set and check (coplanarity) of trd parameters

FQUALIFIER
void GPTrd_CheckAndSetAllParameters ( GPTrd *This,
				      G4double pdx1,  G4double pdx2,
				      G4double pdy1,  G4double pdy2,
				      G4double pdz ) 
{
  if ( pdx1>0&&pdx2>0&&pdy1>0&&pdy2>0&&pdz>0 )
  {
    This->fDx1=pdx1; This->fDx2=pdx2;
    This->fDy1=pdy1; This->fDy2=pdy2;
    This->fDz=pdz;
  }
  else
  {
    if ( pdx1>=0 && pdx2>=0 && pdy1>=0 && pdy2>=0 && pdz>=0 )
    {
      // G4double  Minimum_length= (1+per_thousand) * kCarTolerance/2.;
      // FIX-ME : temporary solution for ZERO or very-small parameters
      //
      G4double  Minimum_length= kCarTolerance/2.;
      This->fDx1=GPfmax(pdx1,Minimum_length); 
      This->fDx2=GPfmax(pdx2,Minimum_length); 
      This->fDy1=GPfmax(pdy1,Minimum_length); 
      This->fDy2=GPfmax(pdy2,Minimum_length); 
      This->fDz=GPfmax(pdz,Minimum_length);
    }
    else
    {
      ;
      /*
      std::ostringstream message;
      message << "Invalid negative dimensions for Solid: " << GetName()
              << G4endl
              << "          X - " << pdx1 << ", " << pdx2 << G4endl
              << "          Y - " << pdy1 << ", " << pdy2 << G4endl
              << "          Z - " << pdz;
      G4Exception("GPTrd_CheckAndSetAllParameters()",
                  "GeomSolids0002", FatalException, message);
      */ 
   }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
//GPTrd_~GPTrd()
//{
//}
//
//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
//GPTrd_GPTrd(const GPTrd& rhs)
//  : G4CSGSolid(rhs), This->fDx1(rhs.This->fDx1), This->fDx2(rhs.This->fDx2),
//    This->fDy1(rhs.This->fDy1), This->fDy2(rhs.This->fDy2), This->fDz(rhs.This->fDz)
//{
//}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
/*
GPTrd& GPTrd_operator = (const GPTrd& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   This->fDx1 = rhs.This->fDx1; This->fDx2 = rhs.This->fDx2;
   This->fDy1 = rhs.This->fDy1; This->fDy2 = rhs.This->fDy2;
   This->fDz = rhs.This->fDz;

   return *this;
}
*/
////////////////////////////////////////////////////////////////////////////
//
//

FQUALIFIER
void GPTrd_SetAllParameters (GPTrd *This, 
			     G4double pdx1, G4double pdx2, G4double pdy1, 
			     G4double pdy2, G4double pdz ) 
{
  GPTrd_CheckAndSetAllParameters (This, pdx1, pdx2, pdy1, pdy2, pdz);
}


/////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
/*
void GPTrd_ComputeDimensions(       G4VPVParameterisation* p,
                               const G4int n,
                               const GPVPhysicalVolume* pRep )
{
  p->ComputeDimensions(*this,n,pRep);
}
*/

///////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

FQUALIFIER
G4bool GPTrd_CalculateExtent(GPTrd *This, 
			     const EAxis pAxis,
			     GPVoxelLimits pVoxelLimit,
			     GPAffineTransform pTransform,
			     G4double *pMin, G4double *pMax )
{
  if (!GPAffineTransform_IsRotated(&pTransform))
  {
    // Special case handling for unrotated solids
    // Compute x/y/z mins and maxs respecting limits, with early returns
    // if outside limits. Then switch() on pAxis

    G4double xoffset,xMin,xMax;
    G4double yoffset,yMin,yMax;
    G4double zoffset,zMin,zMax;

    zoffset= GPAffineTransform_NetTranslation(&pTransform).z;
    zMin=zoffset-This->fDz;
    zMax=zoffset+This->fDz;
    if (GPVoxelLimits_IsZLimited(&pVoxelLimit))
    {
      if ( (zMin>GPVoxelLimits_GetMaxZExtent(&pVoxelLimit)+kCarTolerance)
        || (zMax<GPVoxelLimits_GetMinZExtent(&pVoxelLimit)-kCarTolerance) )
      {
        return false;
      }
        else
      {
        if (zMin<GPVoxelLimits_GetMinZExtent(&pVoxelLimit))
        {
          zMin=GPVoxelLimits_GetMinZExtent(&pVoxelLimit);
        }
        if (zMax>GPVoxelLimits_GetMaxZExtent(&pVoxelLimit))
        {
          zMax=GPVoxelLimits_GetMaxZExtent(&pVoxelLimit);
        }
      }
    }
    xoffset=  GPAffineTransform_NetTranslation(&pTransform).x;
    if (This->fDx2 >= This->fDx1)
    { 
      xMax =  xoffset+(This->fDx1+This->fDx2)/2+(zMax-zoffset)*(This->fDx2-This->fDx1)/(2*This->fDz) ;
      xMin = 2*xoffset - xMax ;
    }
    else
    {
      xMax =  xoffset+(This->fDx1+This->fDx2)/2+(zMin-zoffset)*(This->fDx2-This->fDx1)/(2*This->fDz) ;
      xMin =  2*xoffset - xMax ;
    }   
    if (GPVoxelLimits_IsXLimited(&pVoxelLimit))
    {
      if ( (xMin>GPVoxelLimits_GetMaxXExtent(&pVoxelLimit)+kCarTolerance)
        || (xMax<GPVoxelLimits_GetMinXExtent(&pVoxelLimit)-kCarTolerance) )
      {
        return false;
      }
      else
      {
        if (xMin<GPVoxelLimits_GetMinXExtent(&pVoxelLimit))
        {
          xMin=GPVoxelLimits_GetMinXExtent(&pVoxelLimit);
        }
        if (xMax>GPVoxelLimits_GetMaxXExtent(&pVoxelLimit))
        {
          xMax=GPVoxelLimits_GetMaxXExtent(&pVoxelLimit);
        }
      }
    }
    yoffset=  GPAffineTransform_NetTranslation(&pTransform).y;
    if(This->fDy2 >= This->fDy1)
    {
      yMax = yoffset+(This->fDy2+This->fDy1)/2+(zMax-zoffset)*(This->fDy2-This->fDy1)/(2*This->fDz) ;
      yMin = 2*yoffset - yMax ;
    }
    else
    {
      yMax = yoffset+(This->fDy2+This->fDy1)/2+(zMin-zoffset)*(This->fDy2-This->fDy1)/(2*This->fDz) ;
      yMin = 2*yoffset - yMax ;  
    }  
    if (GPVoxelLimits_IsYLimited(&pVoxelLimit))
    {
      if ( (yMin>GPVoxelLimits_GetMaxYExtent(&pVoxelLimit)+kCarTolerance)
        || (yMax<GPVoxelLimits_GetMinYExtent(&pVoxelLimit)-kCarTolerance) )
      {
        return false;
      }
      else
      {
        if (yMin<GPVoxelLimits_GetMinYExtent(&pVoxelLimit))
        {
          yMin=GPVoxelLimits_GetMinYExtent(&pVoxelLimit);
        }
        if (yMax>GPVoxelLimits_GetMaxYExtent(&pVoxelLimit))
        {
          yMax=GPVoxelLimits_GetMaxYExtent(&pVoxelLimit);
        }
      }
    }

    switch (pAxis)
    {
      case kXAxis:
        *pMin=xMin;
        *pMax=xMax;
        break;
      case kYAxis:
        *pMin=yMin;
        *pMax=yMax;
        break;
      case kZAxis:
        *pMin=zMin;
        *pMax=zMax;
        break;
      default:
        break;
    }

    // Add 2*Tolerance to avoid precision troubles ?
    //
    *pMin-=kCarTolerance;
    *pMax+=kCarTolerance;

    return true;
  }
  else
  {
    // General rotated case - create and clip mesh to boundaries

    G4bool existsAfterClip=false;
    GPThreeVectorList vertices = GPTrd_CreateRotatedVertices(This,pTransform);

    *pMin=+kInfinity;
    *pMax=-kInfinity;

    // Calculate rotated vertex coordinates
    //
    GPVSolid_ClipCrossSection(&vertices,0,&pVoxelLimit,pAxis,pMin,pMax);
    GPVSolid_ClipCrossSection(&vertices,4,&pVoxelLimit,pAxis,pMin,pMax);
    GPVSolid_ClipBetweenSections(&vertices,0,&pVoxelLimit,pAxis,pMin,pMax);
      
    if (*pMin!=kInfinity||*pMax!=-kInfinity)
    {
      existsAfterClip=true;

      // Add 2*tolerance to avoid precision troubles
      //
      *pMin-=kCarTolerance;
      *pMax+=kCarTolerance;
        
    }
    else
    {
      // Check for case where completely enveloping clipping volume
      // If point inside then we are confident that the solid completely
      // envelopes the clipping volume. Hence set min/max extents according
      // to clipping volume extents along the specified axis.

      GPThreeVector clipCentre =  GPThreeVector_create(
      (GPVoxelLimits_GetMinXExtent(&pVoxelLimit)+GPVoxelLimits_GetMaxXExtent(&pVoxelLimit))*0.5,
      (GPVoxelLimits_GetMinYExtent(&pVoxelLimit)+GPVoxelLimits_GetMaxYExtent(&pVoxelLimit))*0.5,
      (GPVoxelLimits_GetMinZExtent(&pVoxelLimit)+GPVoxelLimits_GetMaxZExtent(&pVoxelLimit))*0.5);
        
      GPAffineTransform invT = GPAffineTransform_Inverse( &pTransform );
      GPThreeVector aT = GPAffineTransform_TransformPoint(&invT, clipCentre);

      if (GPTrd_Inside(This,aT) !=kOutside)
      {
        existsAfterClip=true;
        *pMin=GPVoxelLimits_GetMinExtent(&pVoxelLimit,pAxis);
        *pMax=GPVoxelLimits_GetMaxExtent(&pVoxelLimit,pAxis);
      }
    }
    //    delete vertices;
    return existsAfterClip;
  }
}

///////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

FQUALIFIER
EInside GPTrd_Inside( GEOMETRYLOC GPTrd *This,
		      GPThreeVector p )
{  
  EInside in=kOutside;
  G4double x,y,zbase1,zbase2;
    
  if (fabs(p.z)<=This->fDz-kCarTolerance/2)
  {
    zbase1=p.z+This->fDz;  // Dist from -ve z plane
    zbase2=This->fDz-p.z;  // Dist from +ve z plane

    // Check whether inside x tolerance
    //
    x=0.5*(This->fDx2*zbase1+This->fDx1*zbase2)/This->fDz - kCarTolerance/2;
    if (fabs(p.x)<=x)
    {
      y=0.5*((This->fDy2*zbase1+This->fDy1*zbase2))/This->fDz - kCarTolerance/2;
      if (fabs(p.y)<=y)
      {
        in=kInside;
      }
      else if (fabs(p.y)<=y+kCarTolerance)
      {
        in=kSurface;
      }
    }
    else if (fabs(p.x)<=x+kCarTolerance)
    {
      // y = y half width of shape at z of point + tolerant boundary
      //
      y=0.5*((This->fDy2*zbase1+This->fDy1*zbase2))/This->fDz + kCarTolerance/2;
      if (fabs(p.y)<=y)
      {
        in=kSurface;
      }
    }
  }
  else if (fabs(p.z)<=This->fDz+kCarTolerance/2)
  {
    // Only need to check outer tolerant boundaries
    //
    zbase1=p.z+This->fDz;  // Dist from -ve z plane
    zbase2=This->fDz-p.z;   // Dist from +ve z plane

    // x = x half width of shape at z of point plus tolerance
    //
    x=0.5*(This->fDx2*zbase1+This->fDx1*zbase2)/This->fDz + kCarTolerance/2;
    if (fabs(p.x)<=x)
    {
      // y = y half width of shape at z of point
      //
      y=0.5*((This->fDy2*zbase1+This->fDy1*zbase2))/This->fDz + kCarTolerance/2;
      if (fabs(p.y)<=y) in=kSurface;
    }
  }  
  return in;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If two sides are equidistant, normal of first side (x/y/z) 
// encountered returned

FQUALIFIER
GPThreeVector GPTrd_SurfaceNormal( GEOMETRYLOC GPTrd *This,
				   GPThreeVector p )
{
  GPThreeVector norm = GPThreeVector_create(0.,0.,0.);
  GPThreeVector sumnorm = GPThreeVector_create(0.,0.,0.);
  G4int noSurfaces = 0; 
  G4double z = 2.0*This->fDz, tanx, secx, newpx, widx;
  G4double tany, secy, newpy, widy;
  G4double distx, disty, distz, fcos;
  G4double delta = 0.5*kCarTolerance;

  tanx  = (This->fDx2 - This->fDx1)/z;
  secx  = sqrt(1.0+tanx*tanx);
  newpx = fabs(p.x)-p.z*tanx;
  widx  = This->fDx2 - This->fDz*tanx;

  tany  = (This->fDy2 - This->fDy1)/z;
  secy  = sqrt(1.0+tany*tany);
  newpy = fabs(p.y)-p.z*tany;
  widy  = This->fDy2 - This->fDz*tany;

  distx = fabs(newpx-widx)/secx;       // perp. distance to x side
  disty = fabs(newpy-widy)/secy;       //                to y side
  distz = fabs(fabs(p.z)-This->fDz);  //                to z side

  fcos              = 1.0/secx;
  GPThreeVector nX  = GPThreeVector_create( fcos,0,-tanx*fcos);
  GPThreeVector nmX = GPThreeVector_create(-fcos,0,-tanx*fcos);

  fcos              = 1.0/secy;
  GPThreeVector nY  = GPThreeVector_create(0, fcos,-tany*fcos);
  GPThreeVector nmY = GPThreeVector_create(0,-fcos,-tany*fcos);
  GPThreeVector nZ  = GPThreeVector_create( 0, 0,  1.0);
 
  if (distx <= delta)      
  {
    noSurfaces ++;
    if ( p.x >= 0.) sumnorm = GPThreeVector_add(sumnorm,nX);//sumnorm += nX;
    else            sumnorm = GPThreeVector_add(sumnorm,nmX);//sumnorm += nmX;
  }
  if (disty <= delta)
  {
    noSurfaces ++;
    if ( p.y >= 0.) sumnorm = GPThreeVector_add(sumnorm,nY);//sumnorm += nY;
    else            sumnorm = GPThreeVector_add(sumnorm,nmY);//sumnorm += nmY;
  }
  if (distz <= delta)  
  {
    noSurfaces ++;
    if ( p.z >= 0.) sumnorm = GPThreeVector_add(sumnorm,nZ);//sumnorm += nZ;
    else            sumnorm = GPThreeVector_sub(sumnorm,nZ);//sumnorm -= nZ; 
  }
  if ( noSurfaces == 0 )
  {
    norm = GPTrd_ApproxSurfaceNormal(This,p);
  }
  else if ( noSurfaces == 1 ) norm = sumnorm;
  else                        norm = GPThreeVector_unit(sumnorm);
  return norm;   
}


/////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

FQUALIFIER
GPThreeVector GPTrd_ApproxSurfaceNormal( GEOMETRYLOC GPTrd *This,
					 GPThreeVector p )
{
  GPThreeVector norm = GPThreeVector_create(0,0,0);
  G4double z,tanx,secx,newpx,widx;
  G4double tany,secy,newpy,widy;
  G4double distx,disty,distz,fcos;

  z=2.0*This->fDz;

  tanx=(This->fDx2-This->fDx1)/z;
  secx=sqrt(1.0+tanx*tanx);
  newpx=fabs(p.x)-p.z*tanx;
  widx=This->fDx2-This->fDz*tanx;

  tany=(This->fDy2-This->fDy1)/z;
  secy=sqrt(1.0+tany*tany);
  newpy=fabs(p.y)-p.z*tany;
  widy=This->fDy2-This->fDz*tany;

  distx=fabs(newpx-widx)/secx;  // perpendicular distance to x side
  disty=fabs(newpy-widy)/secy;  //                        to y side
  distz=fabs(fabs(p.z)-This->fDz);  //                        to z side

  // find closest side
  //
  if (distx<=disty)
  { 
    if (distx<=distz) 
    {
      // Closest to X
      //
      fcos=1.0/secx;
      // normal=(+/-std::cos(ang),0,-std::sin(ang))
      if (p.x>=0)
        norm=GPThreeVector_create(fcos,0,-tanx*fcos);
      else
        norm=GPThreeVector_create(-fcos,0,-tanx*fcos);
    }
    else
    {
      // Closest to Z
      //
      if (p.z>=0)
        norm=GPThreeVector_create(0,0,1);
      else
        norm=GPThreeVector_create(0,0,-1);
    }
  }
  else
  {  
    if (disty<=distz)
    {
      // Closest to Y
      //
      fcos=1.0/secy;
      if (p.y>=0)
        norm=GPThreeVector_create(0,fcos,-tany*fcos);
      else
        norm=GPThreeVector_create(0,-fcos,-tany*fcos);
    }
    else 
    {
      // Closest to Z
      //
      if (p.z>=0)
        norm=GPThreeVector_create(0,0,1);
      else
        norm=GPThreeVector_create(0,0,-1);
    }
  }
  return norm;   
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside
// - return kInfinity if no intersection
//
// ALGORITHM:
// For each component, calculate pair of minimum and maximum intersection
// values for which the particle is in the extent of the shape
// - The smallest (MAX minimum) allowed distance of the pairs is intersect
// - Z plane intersectin uses tolerance
// - XZ YZ planes use logic & *SLIGHTLY INCORRECT* tolerance
//   (this saves at least 1 sqrt, 1 multiply and 1 divide... in applicable
//    cases)
// - Note: XZ and YZ planes each divide space into four regions,
//   characterised by ss1 ss2
// NOTE:
//
// `Inside' safe - meaningful answers given if point is inside the exact
// shape.

FQUALIFIER
G4double GPTrd_DistanceToIn2( GEOMETRYLOC GPTrd *This, 
			      GPThreeVector p,
			      GPThreeVector v )
{  
  G4double snxt = kInfinity ;    // snxt = default return value
  G4double smin,smax;
  G4double s1,s2,tanxz,tanyz,ds1,ds2;
  G4double ss1,ss2,sn1=0.,sn2=0.,Dist;

  if ( v.z )  // Calculate valid z intersect range
  {
    if ( v.z > 0 )   // Calculate smax: must be +ve or no intersection.
    {
      Dist = This->fDz - p.z ;  // to plane at +dz

      if (Dist >= 0.5*kCarTolerance)
      {
        smax = Dist/v.z ;
        smin = -(This->fDz + p.z)/v.z ;
      }
      else  return snxt ;
    }
    else // v.z <0
    {
      Dist=This->fDz+p.z;  // plane at -dz

      if ( Dist >= 0.5*kCarTolerance )
      {
        smax = -Dist/v.z ;
        smin = (This->fDz - p.z)/v.z ;
      }
      else return snxt ; 
    }
    if (smin < 0 ) smin = 0 ;
  }
  else // v.z=0
  {
    if (fabs(p.z) >= This->fDz ) return snxt ;     // Outside & no intersect
    else
    {
      smin = 0 ;    // Always inside z range
      smax = kInfinity;
    }
  }

  // Calculate x intersection range
  //
  // Calc half width at p.z, and components towards planes

  tanxz = (This->fDx2 - This->fDx1)*0.5/This->fDz ;
  s1    = 0.5*(This->fDx1+This->fDx2) + tanxz*p.z ;  // x half width at p.z
  ds1   = v.x - tanxz*v.z ;       // Components of v towards faces at +-x
  ds2   = v.x + tanxz*v.z ;
  ss1   = s1 - p.x ;         // -delta x to +ve plane
                               // -ve when outside
  ss2   = -s1 - p.x ;        // -delta x to -ve plane
                               // +ve when outside

  if (ss1 < 0 && ss2 <= 0 )
  {
    if (ds1 < 0)   // In +ve coord Area
    {
      sn1 = ss1/ds1 ;

      if ( ds2 < 0 ) sn2 = ss2/ds2 ;           
      else           sn2 = kInfinity ;
    }
    else return snxt ;
  }
  else if ( ss1 >= 0 && ss2 > 0 )
  {
    if ( ds2 > 0 )  // In -ve coord Area
    {
      sn1 = ss2/ds2 ;

      if (ds1 > 0)  sn2 = ss1/ds1 ;      
      else          sn2 = kInfinity;      
        
    }
    else   return snxt ;
  }
  else if (ss1 >= 0 && ss2 <= 0 )
  {
    // Inside Area - calculate leaving distance
    // *Don't* use exact distance to side for tolerance
    //                                             = ss1*std::cos(ang xz)
    //                                             = ss1/sqrt(1.0+tanxz*tanxz)
    sn1 = 0 ;

    if ( ds1 > 0 )
    {
      if (ss1 > 0.5*kCarTolerance) sn2 = ss1/ds1 ; // Leave +ve side extent
      else                         return snxt ;   // Leave immediately by +ve 
    }
    else  sn2 = kInfinity ;
      
    if ( ds2 < 0 )
    {
      if ( ss2 < -0.5*kCarTolerance )
      {
        Dist = ss2/ds2 ;            // Leave -ve side extent
        if ( Dist < sn2 ) sn2 = Dist ;
      }
      else  return snxt ;
    }    

  }
  else if (ss1 < 0 && ss2 > 0 )
  {
    // Within +/- plane cross-over areas (not on boundaries ss1||ss2==0)

    if ( ds1 >= 0 || ds2 <= 0 )
    {   
      return snxt ;
    }
    else  // Will intersect & stay inside
    {
      sn1  = ss1/ds1 ;
      Dist = ss2/ds2 ;
      if (Dist > sn1 ) sn1 = Dist ;
      sn2 = kInfinity ;
    }
  }

  // Reduce allowed range of distances as appropriate

  if ( sn1 > smin ) smin = sn1 ;
  if ( sn2 < smax ) smax = sn2 ;

  // Check for incompatible ranges (eg z intersects between 50 ->100 and x
  // only 10-40 -> no intersection)

  if ( smax < smin ) return snxt ;

  // Calculate valid y intersection range 
  // (repeat of x intersection code)

  tanyz = (This->fDy2-This->fDy1)*0.5/This->fDz ;
  s2    = 0.5*(This->fDy1+This->fDy2) + tanyz*p.z ;  // y half width at p.z
  ds1   = v.y - tanyz*v.z ;       // Components of v towards faces at +-y
  ds2   = v.y + tanyz*v.z ;
  ss1   = s2 - p.y ;         // -delta y to +ve plane
  ss2   = -s2 - p.y ;        // -delta y to -ve plane

  if ( ss1 < 0 && ss2 <= 0 )
  {
    if (ds1 < 0 ) // In +ve coord Area
    {
      sn1 = ss1/ds1 ;
      if ( ds2 < 0 )  sn2 = ss2/ds2 ;
      else            sn2 = kInfinity ;
    }
    else   return snxt ;
  }
  else if ( ss1 >= 0 && ss2 > 0 )
  {
    if ( ds2 > 0 )  // In -ve coord Area
    {
      sn1 = ss2/ds2 ;
      if ( ds1 > 0 )  sn2 = ss1/ds1 ;
      else            sn2 = kInfinity ;      
    }
    else   return snxt ;
  }
  else if (ss1 >= 0 && ss2 <= 0 )
  {
    // Inside Area - calculate leaving distance
    // *Don't* use exact distance to side for tolerance
    //                                          = ss1*std::cos(ang yz)
    //                                          = ss1/sqrt(1.0+tanyz*tanyz)
    sn1 = 0 ;

    if ( ds1 > 0 )
    {
      if (ss1 > 0.5*kCarTolerance) sn2 = ss1/ds1 ; // Leave +ve side extent
      else                         return snxt ;   // Leave immediately by +ve
    }
    else  sn2 = kInfinity ;
      
    if ( ds2 < 0 )
    {
      if ( ss2 < -0.5*kCarTolerance )
      {
        Dist = ss2/ds2 ; // Leave -ve side extent
        if (Dist < sn2) sn2=Dist;
      }
      else  return snxt ;
    }    
  }
  else if (ss1 < 0 && ss2 > 0 )
  {
    // Within +/- plane cross-over areas (not on boundaries ss1||ss2==0)

    if (ds1 >= 0 || ds2 <= 0 )  
    {
      return snxt ;
    }
    else  // Will intersect & stay inside
    {
      sn1 = ss1/ds1 ;
      Dist = ss2/ds2 ;
      if (Dist > sn1 ) sn1 = Dist ;
      sn2 = kInfinity ;
    }
  }
  
  // Reduce allowed range of distances as appropriate

  if ( sn1 > smin) smin = sn1 ;
  if ( sn2 < smax) smax = sn2 ;

  // Check for incompatible ranges (eg x intersects between 50 ->100 and y
  // only 10-40 -> no intersection). Set snxt if ok

  if ( smax > smin ) snxt = smin ;

  if (snxt < 0.5*kCarTolerance ) snxt = 0.0 ;

  return snxt ;
}

/////////////////////////////////////////////////////////////////////////
//
// Approximate distance to shape
// Calculate perpendicular distances to z/x/y surfaces, return largest
// which is the most fast estimation of shortest distance to Trd
//  - Safe underestimate
//  - If point within exact shape, return 0 

FQUALIFIER
G4double GPTrd_DistanceToIn( GEOMETRYLOC GPTrd *This,
			     GPThreeVector p )
{
  G4double safe=0.0;
  G4double tanxz,distx,safx;
  G4double tanyz,disty,safy;
  G4double zbase;

  safe=fabs(p.z)-This->fDz;
  if (safe<0) safe=0;      // Also used to ensure x/y distances
                           // POSITIVE 

  zbase=This->fDz+p.z;

  // Find distance along x direction to closest x plane
  //
  tanxz=(This->fDx2-This->fDx1)*0.5/This->fDz;
  //    widx=This->fDx1+tanxz*(This->fDz+p.z); // x width at p.z
  //    distx=fabs(p.x)-widx;      // distance to plane
  distx=fabs(p.x)-(This->fDx1+tanxz*zbase);
  if (distx>safe)
  {
    safx=distx/sqrt(1.0+tanxz*tanxz); // vector Dist=Dist*std::cos(ang)
    if (safx>safe) safe=safx;
  }

  // Find distance along y direction to slanted wall
  tanyz=(This->fDy2-This->fDy1)*0.5/This->fDz;
  //    widy=This->fDy1+tanyz*(This->fDz+p.z); // y width at p.z
  //    disty=fabs(p.y)-widy;      // distance to plane
  disty=fabs(p.y)-(This->fDy1+tanyz*zbase);
  if (disty>safe)    
  {
    safy=disty/sqrt(1.0+tanyz*tanyz); // distance along vector
    if (safy>safe) safe=safy;
  }
  return safe;
}

////////////////////////////////////////////////////////////////////////
//
// Calcluate distance to surface of shape from inside
// Calculate distance to x/y/z planes - smallest is exiting distance
// - z planes have std. check for tolerance
// - xz yz planes have check based on distance || to x or y axis
//   (not corrected for slope of planes)
// ?BUG? If v.z==0 are there cases when snside not set????

FQUALIFIER
G4double GPTrd_DistanceToOut2( GEOMETRYLOC  GPTrd *This,
			       GPThreeVector p,
			       GPThreeVector v,
			       const G4bool calcNorm,
			       G4bool *validNorm,
			       GPThreeVector *n )
{
  enum ESide {kUndefined,kPX,kMX,kPY,kMY,kPZ,kMZ};

  ESide side = kUndefined, snside = kUndefined;
  G4double snxt,pdist;
  G4double central,ss1,ss2,ds1,ds2,sn=0.,sn2=0.;
  G4double tanxz=0.,cosxz=0.,tanyz=0.,cosyz=0.;

  if (calcNorm) *validNorm=true; // All normals are valid

  // Calculate z plane intersection
  if (v.z>0)
  {
    pdist=This->fDz-p.z;
    if (pdist>kCarTolerance/2)
    {
      snxt=pdist/v.z;
      side=kPZ;
    }
    else
    {
      if (calcNorm)
      {
        *n=GPThreeVector_create(0,0,1);
      }
      return snxt=0;
    }
  }
  else if (v.z<0) 
  {
    pdist=This->fDz+p.z;
    if (pdist>kCarTolerance/2)
    {
      snxt=-pdist/v.z;
      side=kMZ;
    }
    else
    {
      if (calcNorm)
      {
        *n=GPThreeVector_create(0,0,-1);
      }
      return snxt=0;
    }
  }
  else
  {
    snxt=kInfinity;
  }

  //
  // Calculate x intersection
  //
  tanxz=(This->fDx2-This->fDx1)*0.5/This->fDz;
  central=0.5*(This->fDx1+This->fDx2);

  // +ve plane (1)
  //
  ss1=central+tanxz*p.z-p.x;  // distance || x axis to plane
                                  // (+ve if point inside)
  ds1=v.x-tanxz*v.z;    // component towards plane at +x
                            // (-ve if +ve -> -ve direction)
  // -ve plane (2)
  //
  ss2=-tanxz*p.z-p.x-central;  //distance || x axis to plane
                                   // (-ve if point inside)
  ds2=tanxz*v.z+v.x;    // component towards plane at -x

  if (ss1>0&&ss2<0)
  {
    // Normal case - entirely inside region
    if (ds1<=0&&ds2<0)
    {   
      if (ss2<-kCarTolerance/2)
      {
        sn=ss2/ds2;  // Leave by -ve side
        snside=kMX;
      }
      else
      {
        sn=0; // Leave immediately by -ve side
        snside=kMX;
      }
    }
    else if (ds1>0&&ds2>=0)
    {
      if (ss1>kCarTolerance/2)
      {
        sn=ss1/ds1;  // Leave by +ve side
        snside=kPX;
      }
      else
      {
        sn=0; // Leave immediately by +ve side
        snside=kPX;
      }
    }
    else if (ds1>0&&ds2<0)
    {
      if (ss1>kCarTolerance/2)
      {
        // sn=ss1/ds1;  // Leave by +ve side
        if (ss2<-kCarTolerance/2)
        {
          sn=ss1/ds1;  // Leave by +ve side
          sn2=ss2/ds2;
          if (sn2<sn)
          {
            sn=sn2;
            snside=kMX;
          }
          else
          {
            snside=kPX;
          }
        }
        else
        {
          sn=0; // Leave immediately by -ve
          snside=kMX;
        }      
      }
      else
      {
        sn=0; // Leave immediately by +ve side
        snside=kPX;
      }
    }
    else
    {
      // Must be || to both
      //
      sn=kInfinity;    // Don't leave by either side
    }
  }
  else if (ss1<=0&&ss2<0)
  {
    // Outside, in +ve Area
    
    if (ds1>0)
    {
      sn=0;       // Away from shape
                  // Left by +ve side
      snside=kPX;
    }
    else
    {
      if (ds2<0)
      {
        // Ignore +ve plane and use -ve plane intersect
        //
        sn=ss2/ds2; // Leave by -ve side
        snside=kMX;
      }
      else
      {
        // Must be || to both -> exit determined by other axes
        //
        sn=kInfinity; // Don't leave by either side
      }
    }
  }
  else if (ss1>0&&ss2>=0)
  {
    // Outside, in -ve Area

    if (ds2<0)
    {
      sn=0;       // away from shape
                  // Left by -ve side
      snside=kMX;
    }
    else
    {
      if (ds1>0)
      {
        // Ignore +ve plane and use -ve plane intersect
        //
        sn=ss1/ds1; // Leave by +ve side
        snside=kPX;
      }
      else
      {
        // Must be || to both -> exit determined by other axes
        //
        sn=kInfinity; // Don't leave by either side
      }
    }
  }

  // Update minimum exit distance

  if (sn<snxt)
  {
    snxt=sn;
    side=snside;
  }
  if (snxt>0)
  {
    // Calculate y intersection

    tanyz=(This->fDy2-This->fDy1)*0.5/This->fDz;
    central=0.5*(This->fDy1+This->fDy2);

    // +ve plane (1)
    //
    ss1=central+tanyz*p.z-p.y; // distance || y axis to plane
                                   // (+ve if point inside)
    ds1=v.y-tanyz*v.z;  // component towards +ve plane
                            // (-ve if +ve -> -ve direction)
    // -ve plane (2)
    //
    ss2=-tanyz*p.z-p.y-central; // distance || y axis to plane
                                    // (-ve if point inside)
    ds2=tanyz*v.z+v.y;  // component towards -ve plane

    if (ss1>0&&ss2<0)
    {
      // Normal case - entirely inside region

      if (ds1<=0&&ds2<0)
      {   
        if (ss2<-kCarTolerance/2)
        {
          sn=ss2/ds2;  // Leave by -ve side
          snside=kMY;
        }
        else
        {
          sn=0; // Leave immediately by -ve side
          snside=kMY;
        }
      }
      else if (ds1>0&&ds2>=0)
      {
        if (ss1>kCarTolerance/2)
        {
          sn=ss1/ds1;  // Leave by +ve side
          snside=kPY;
        }
        else
        {
          sn=0; // Leave immediately by +ve side
          snside=kPY;
        }
      }
      else if (ds1>0&&ds2<0)
      {
        if (ss1>kCarTolerance/2)
        {
          // sn=ss1/ds1;  // Leave by +ve side
          if (ss2<-kCarTolerance/2)
          {
            sn=ss1/ds1;  // Leave by +ve side
            sn2=ss2/ds2;
            if (sn2<sn)
            {
              sn=sn2;
              snside=kMY;
            }
            else
            {
              snside=kPY;
            }
          }
          else
          {
            sn=0; // Leave immediately by -ve
            snside=kMY;
          }
        }
        else
        {
          sn=0; // Leave immediately by +ve side
          snside=kPY;
        }
      }
      else
      {
        // Must be || to both
        //
        sn=kInfinity;    // Don't leave by either side
      }
    }
    else if (ss1<=0&&ss2<0)
    {
      // Outside, in +ve Area

      if (ds1>0)
      {
        sn=0;       // Away from shape
                    // Left by +ve side
        snside=kPY;
      }
      else
      {
        if (ds2<0)
        {
          // Ignore +ve plane and use -ve plane intersect
          //
          sn=ss2/ds2; // Leave by -ve side
          snside=kMY;
        }
        else
        {
          // Must be || to both -> exit determined by other axes
          //
          sn=kInfinity; // Don't leave by either side
        }
      }
    }
    else if (ss1>0&&ss2>=0)
    {
      // Outside, in -ve Area
      if (ds2<0)
      {
        sn=0;       // away from shape
                    // Left by -ve side
        snside=kMY;
      }
      else
      {
        if (ds1>0)
        {
          // Ignore +ve plane and use -ve plane intersect
          //
          sn=ss1/ds1; // Leave by +ve side
          snside=kPY;
        }
        else
        {
          // Must be || to both -> exit determined by other axes
          //
          sn=kInfinity; // Don't leave by either side
        }
      }
    }

    // Update minimum exit distance

    if (sn<snxt)
    {
      snxt=sn;
      side=snside;
    }
  }

  if (calcNorm)
  {
    switch (side)
    {
      case kPX:
        cosxz=1.0/sqrt(1.0+tanxz*tanxz);
        *n=GPThreeVector_create(cosxz,0,-tanxz*cosxz);
        break;
      case kMX:
        cosxz=-1.0/sqrt(1.0+tanxz*tanxz);
        *n=GPThreeVector_create(cosxz,0,tanxz*cosxz);
        break;
      case kPY:
        cosyz=1.0/sqrt(1.0+tanyz*tanyz);
        *n=GPThreeVector_create(0,cosyz,-tanyz*cosyz);
        break;
      case kMY:
        cosyz=-1.0/sqrt(1.0+tanyz*tanyz);
        *n=GPThreeVector_create(0,cosyz,tanyz*cosyz);
        break;
      case kPZ:
        *n=GPThreeVector_create(0,0,1);
        break;
      case kMZ:
        *n=GPThreeVector_create(0,0,-1);
        break;
      default:
	/*
        DumpInfo();
        G4Exception("GPTrd_DistanceToOut(p,v,..)",
                    "GeomSolids1002", JustWarning, 
                    "Undefined side for valid surface normal to solid.");
	*/
        break;
    }
  }
  return snxt; 
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - Returns 0 is point outside

FQUALIFIER
G4double GPTrd_DistanceToOut( GEOMETRYLOC GPTrd *This, 
			      GPThreeVector p )
{
  G4double safe=0.0;
  G4double tanxz,xdist,saf1;
  G4double tanyz,ydist,saf2;
  G4double zbase;

  safe=This->fDz-fabs(p.z);  // z perpendicular Dist

  zbase=This->fDz+p.z;

  // xdist = distance perpendicular to z axis to closest x plane from p
  //       = (x half width of shape at p.z) - fabs(p.x)
  //
  tanxz=(This->fDx2-This->fDx1)*0.5/This->fDz;
  xdist=This->fDx1+tanxz*zbase-fabs(p.x);
  saf1=xdist/sqrt(1.0+tanxz*tanxz); // x*std::cos(ang_xz) =
                                    // shortest (perpendicular)
                                    // distance to plane
  tanyz=(This->fDy2-This->fDy1)*0.5/This->fDz;
  ydist=This->fDy1+tanyz*zbase-fabs(p.y);
  saf2=ydist/sqrt(1.0+tanyz*tanyz);

  // Return minimum x/y/z distance
  //
  if (safe>saf1) safe=saf1;
  if (safe>saf2) safe=saf2;

  if (safe<0) safe=0;
  return safe;     
}

////////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Ordering [0-3] -This->fDz cross section
//          [4-7] +This->fDz cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility

FQUALIFIER
GPThreeVectorList
GPTrd_CreateRotatedVertices(GPTrd *This, 
			    GPAffineTransform pTransform )
{
  GPThreeVectorList vertices;
  //  vertices=new GPThreeVectorList();
  GPThreeVectorList_Constructor(&vertices);

  if (&vertices)
  {
    //    vertices->reserve(8);
    GPThreeVector vertex0 = GPThreeVector_create(-This->fDx1,-This->fDy1,-This->fDz);
    GPThreeVector vertex1 = GPThreeVector_create(This->fDx1,-This->fDy1,-This->fDz);
    GPThreeVector vertex2 = GPThreeVector_create(This->fDx1,This->fDy1,-This->fDz);
    GPThreeVector vertex3 = GPThreeVector_create(-This->fDx1,This->fDy1,-This->fDz);
    GPThreeVector vertex4 = GPThreeVector_create(-This->fDx2,-This->fDy2,This->fDz);
    GPThreeVector vertex5 = GPThreeVector_create(This->fDx2,-This->fDy2,This->fDz);
    GPThreeVector vertex6 = GPThreeVector_create(This->fDx2,This->fDy2,This->fDz);
    GPThreeVector vertex7 = GPThreeVector_create(-This->fDx2,This->fDy2,This->fDz);

    GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex0));
    GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex1));
    GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex2));
    GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex3));
    GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex4));
    GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex5));
    GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex6));
    GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex7));
  }
  else
  {
    ;
    /*
    DumpInfo();
    G4Exception("GPTrd_CreateRotatedVertices()",
                "GeomSolids0003", FatalException,
                "Error in allocation of vertices. Out of memory !");
    */
  }
  return vertices;
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

//G4GeometryType GPTrd_GetEntityType() const
//{
//  return G4String("GPTrd");
//}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
//GPVSolid* GPTrd_Clone() const
//{
//  return new GPTrd(*this);
//}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface
//
// Return a point (GPThreeVector) randomly and uniformly
// selected on the solid surface
/*
GPThreeVector GPTrd_GetPointOnSurface() const
{
  G4double px, py, pz, tgX, tgY, secX, secY, select, sumS, tmp;
  G4double Sxy1, Sxy2, Sxy, Sxz, Syz;

  tgX  = 0.5*(This->fDx2-This->fDx1)/This->fDz;
  secX = sqrt(1+tgX*tgX);
  tgY  = 0.5*(This->fDy2-This->fDy1)/This->fDz;
  secY = sqrt(1+tgY*tgY);

  // calculate 0.25 of side surfaces, sumS is 0.25 of total surface

  Sxy1 = This->fDx1*This->fDy1; 
  Sxy2 = This->fDx2*This->fDy2;
  Sxy  = Sxy1 + Sxy2; 
  Sxz  = (This->fDx1 + This->fDx2)*This->fDz*secY; 
  Syz  = (This->fDy1 + This->fDy2)*This->fDz*secX;
  sumS = Sxy + Sxz + Syz;

  select = sumS*G4UniformRand();
 
  if( select < Sxy )                  // Sxy1 or Sxy2
  {
    if( select < Sxy1 ) 
    {
      pz = -This->fDz;
      px = -This->fDx1 + 2*This->fDx1*G4UniformRand();
      py = -This->fDy1 + 2*This->fDy1*G4UniformRand();
    }
    else      
    {
      pz =  This->fDz;
      px = -This->fDx2 + 2*This->fDx2*G4UniformRand();
      py = -This->fDy2 + 2*This->fDy2*G4UniformRand();
    }
  }
  else if ( ( select - Sxy ) < Sxz )    // Sxz
  {
    pz  = -This->fDz  + 2*This->fDz*G4UniformRand();
    tmp =  This->fDx1 + (pz + This->fDz)*tgX;
    px  = -tmp  + 2*tmp*G4UniformRand();
    tmp =  This->fDy1 + (pz + This->fDz)*tgY;

    if(G4UniformRand() > 0.5) { py =  tmp; }
    else                      { py = -tmp; }
  }
  else                                   // Syz
  {
    pz  = -This->fDz  + 2*This->fDz*G4UniformRand();
    tmp =  This->fDy1 + (pz + This->fDz)*tgY;
    py  = -tmp  + 2*tmp*G4UniformRand();
    tmp =  This->fDx1 + (pz + This->fDz)*tgX;

    if(G4UniformRand() > 0.5) { px =  tmp; }
    else                      { px = -tmp; }
  } 
  return GPThreeVector(px,py,pz);
}
*/
///////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

 /*
void GPTrd_DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this);
}

G4Polyhedron* GPTrd_CreatePolyhedron () const
{
  return new G4PolyhedronTrd2 (This->fDx1, This->fDx2, This->fDy1, This->fDy2, This->fDz);
}

G4NURBS* GPTrd_CreateNURBS () const
{
  //  return new G4NURBSbox (fDx, fDy, This->fDz);
  return 0;
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
// $Id: GPTrd.icc,v 1.7 2006-10-19 15:33:37 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 inline definitions file
//
// GPTrd.icc
//
// Implementation of inline methods of GPTrd
// --------------------------------------------------------------------

FQUALIFIER
G4double GPTrd_GetXHalfLength1(GPTrd *This)
{
  return This->fDx1;
}

FQUALIFIER
G4double GPTrd_GetXHalfLength2(GPTrd *This)
{
  return This->fDx2;
}

FQUALIFIER
G4double GPTrd_GetYHalfLength1(GPTrd *This)
{
  return This->fDy1;
}

FQUALIFIER
G4double GPTrd_GetYHalfLength2(GPTrd *This)
{
  return This->fDy2;
}

FQUALIFIER
G4double GPTrd_GetZHalfLength(GPTrd *This)
{
  return This->fDz;
}

FQUALIFIER
void GPTrd_SetXHalfLength1(GPTrd *This,
			   G4double val)
{
  This->fDx1= val;
}

FQUALIFIER
void GPTrd_SetXHalfLength2(GPTrd *This,
			   G4double val)
{
  This->fDx2= val;
}

FQUALIFIER
void GPTrd_SetYHalfLength1(GPTrd *This,
			   G4double val)
{
  This->fDy1= val;
}

FQUALIFIER
void GPTrd_SetYHalfLength2(GPTrd *This,
			   G4double val)
{
  This->fDy2= val;
}

FQUALIFIER
void GPTrd_SetZHalfLength(GPTrd *This,
			  G4double val)
{
  This->fDz= val;
}
