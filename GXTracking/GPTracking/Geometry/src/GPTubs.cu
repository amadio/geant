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
// $Id: G4Tubs.cc,v 1.84 2010-10-19 15:42:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4Tubs
//

#include "GPTubs.h"
#include "GPGeomdefs.h"

FQUALIFIER
void GPTubs_Constructor( GPTubs *This, G4double pRMin, G4double pRMax,
			 G4double pDz, G4double pSPhi, G4double pDPhi )
{
  This->fSolid.fType = kTubs;

  This->fRMin = pRMin;
  This->fRMax = pRMax;
  This->fDz = pDz;
  This->fSPhi = pSPhi;
  This->fDPhi = pDPhi;

  GPTubs_CheckPhiAngles(This, pSPhi, pDPhi);
}

FQUALIFIER
G4bool GPTubs_CalculateExtent(GPTubs *This, 
			      const EAxis       pAxis,
			      GPVoxelLimits     pVoxelLimit,
			      GPAffineTransform pTransform,
			      G4double          *pMin, 
			      G4double          *pMax    )
{

  if ( (! GPAffineTransform_IsRotated(&pTransform) )&& (This->fDPhi == twopi) && (This->fRMin == 0) )
  {
    // Special case handling for unrotated solid tubes
    // Compute x/y/z mins and maxs fro bounding box respecting limits,
    // with early returns if outside limits. Then switch() on pAxis,
    // and compute exact x and y limit for x/y case
      
    G4double xoffset, xMin, xMax;
    G4double yoffset, yMin, yMax;
    G4double zoffset, zMin, zMax;

    G4double diff1, diff2, maxDiff, newMin, newMax;
    G4double xoff1, xoff2, yoff1, yoff2, delta;

    xoffset =  GPAffineTransform_NetTranslation(&pTransform).x ;
    xMin = xoffset - This->fRMax;
    xMax = xoffset + This->fRMax;

    if (GPVoxelLimits_IsXLimited(&pVoxelLimit))
    {
      if ( (xMin >  GPVoxelLimits_GetMaxXExtent(&pVoxelLimit))
        || (xMax <  GPVoxelLimits_GetMinXExtent(&pVoxelLimit)) )
      {
        return false;
      }
      else
      {
        if (xMin <  GPVoxelLimits_GetMinXExtent(&pVoxelLimit))
        {
          xMin =  GPVoxelLimits_GetMinXExtent(&pVoxelLimit);
        }
        if (xMax >  GPVoxelLimits_GetMaxXExtent(&pVoxelLimit))
        {
          xMax =  GPVoxelLimits_GetMaxXExtent(&pVoxelLimit);
        }
      }
    }
    yoffset =  GPAffineTransform_NetTranslation(&pTransform).y ;
    yMin    = yoffset - This->fRMax;
    yMax    = yoffset + This->fRMax;

    if ( GPVoxelLimits_IsYLimited(&pVoxelLimit) )
    {
      if ( (yMin >  GPVoxelLimits_GetMaxYExtent(&pVoxelLimit))
        || (yMax <  GPVoxelLimits_GetMinYExtent(&pVoxelLimit)) )
      {
        return false;
      }
      else
      {
        if (yMin <  GPVoxelLimits_GetMinYExtent(&pVoxelLimit))
        {
          yMin =  GPVoxelLimits_GetMinYExtent(&pVoxelLimit);
        }
        if (yMax >  GPVoxelLimits_GetMaxYExtent(&pVoxelLimit))
        {
          yMax= GPVoxelLimits_GetMaxYExtent(&pVoxelLimit);
        }
      }
    }
    zoffset =  GPAffineTransform_NetTranslation(&pTransform).z ;
    zMin    = zoffset - This->fDz;
    zMax    = zoffset + This->fDz;

    if ( GPVoxelLimits_IsZLimited(&pVoxelLimit) )
    {
      if ( (zMin >  GPVoxelLimits_GetMaxZExtent(&pVoxelLimit))
        || (zMax <  GPVoxelLimits_GetMinZExtent(&pVoxelLimit)) )
      {
        return false;
      }
      else
      {
        if (zMin <  GPVoxelLimits_GetMinZExtent(&pVoxelLimit))
        {
          zMin =  GPVoxelLimits_GetMinZExtent(&pVoxelLimit);
        }
        if (zMax >  GPVoxelLimits_GetMaxZExtent(&pVoxelLimit))
        {
          zMax =  GPVoxelLimits_GetMaxZExtent(&pVoxelLimit);
        }
      }
    }
    switch ( pAxis )  // Known to cut cylinder
    {
      case kXAxis :
      {
        yoff1 = yoffset - yMin;
        yoff2 = yMax    - yoffset;

        if ( (yoff1 >= 0) && (yoff2 >= 0) ) // Y limits cross max/min x
        {                                   // => no change
          *pMin = xMin;
          *pMax = xMax;
        }
        else
        {
          // Y limits don't cross max/min x => compute max delta x,
          // hence new mins/maxs

          delta   = This->fRMax*This->fRMax - yoff1*yoff1;
          diff1   = (delta>0.) ? sqrt(delta) : 0.;
          delta   = This->fRMax*This->fRMax - yoff2*yoff2;
          diff2   = (delta>0.) ? sqrt(delta) : 0.;
          maxDiff = (diff1 > diff2) ? diff1:diff2;
          newMin  = xoffset - maxDiff;
          newMax  = xoffset + maxDiff;
          *pMin    = (newMin < xMin) ? xMin : newMin;
          *pMax    = (newMax > xMax) ? xMax : newMax;
        }    
        break;
      }
      case kYAxis :
      {
        xoff1 = xoffset - xMin;
        xoff2 = xMax - xoffset;

        if ( (xoff1 >= 0) && (xoff2 >= 0) ) // X limits cross max/min y
        {                                   // => no change
          *pMin = yMin;
          *pMax = yMax;
        }
        else
        {
          // X limits don't cross max/min y => compute max delta y,
          // hence new mins/maxs

          delta   = This->fRMax*This->fRMax - xoff1*xoff1;
          diff1   = (delta>0.) ? sqrt(delta) : 0.;
          delta   = This->fRMax*This->fRMax - xoff2*xoff2;
          diff2   = (delta>0.) ? sqrt(delta) : 0.;
          maxDiff = (diff1 > diff2) ? diff1 : diff2;
          newMin  = yoffset - maxDiff;
          newMax  = yoffset + maxDiff;
          *pMin    = (newMin < yMin) ? yMin : newMin;
          *pMax    = (newMax > yMax) ? yMax : newMax;
        }
        break;
      }
      case kZAxis:
      {
        *pMin = zMin;
        *pMax = zMax;
        break;
      }
      default:
        break;
    }
    *pMin -= kCarTolerance;
    *pMax += kCarTolerance;
    return true;
  }
  else // Calculate rotated vertex coordinates
  {
    G4int i, noEntries, noBetweenSections4;
    G4bool existsAfterClip = false;
    GPThreeVectorList vertices = GPTubs_CreateRotatedVertices(This,pTransform);

    *pMin =  kInfinity;
    *pMax = -kInfinity;

    noEntries = GPThreeVectorList_size(&vertices);
    noBetweenSections4 = noEntries - 4;
    
    for ( i = 0 ; i < noEntries ; i += 4 )
    {
      GPVSolid_ClipCrossSection(&vertices, i, &pVoxelLimit, pAxis, pMin, pMax);
    }
    for ( i = 0 ; i < noBetweenSections4 ; i += 4 )
    {
      GPVSolid_ClipBetweenSections(&vertices, i, &pVoxelLimit, pAxis, pMin, pMax);
    }
    if ( (*pMin != kInfinity) || (*pMax != -kInfinity) )
    {
      existsAfterClip = true;
      *pMin -= kCarTolerance; // Add 2*tolerance to avoid precision troubles
      *pMax += kCarTolerance;
    }
    else
    {
      // Check for case where completely enveloping clipping volume
      // If point inside then we are confident that the solid completely
      // envelopes the clipping volume. Hence set min/max extents according
      // to clipping volume extents along the specified axis.

      GPThreeVector clipCentre = GPThreeVector_create(
             ( GPVoxelLimits_GetMinXExtent(&pVoxelLimit)+ GPVoxelLimits_GetMaxXExtent(&pVoxelLimit))*0.5,
             ( GPVoxelLimits_GetMinYExtent(&pVoxelLimit)+ GPVoxelLimits_GetMaxYExtent(&pVoxelLimit))*0.5,
             ( GPVoxelLimits_GetMinZExtent(&pVoxelLimit)+ GPVoxelLimits_GetMaxZExtent(&pVoxelLimit))*0.5 );

      GPAffineTransform invT = GPAffineTransform_Inverse( &pTransform );
      GPThreeVector aT = GPAffineTransform_TransformPoint(&invT, clipCentre);
        
      if ( GPTubs_Inside(This,aT) != kOutside )
      {
        existsAfterClip = true;
        *pMin            = GPVoxelLimits_GetMinExtent(&pVoxelLimit,pAxis);
        *pMax            = GPVoxelLimits_GetMaxExtent(&pVoxelLimit,pAxis);
      }
    }
    //    delete vertices;
    return existsAfterClip;
  }
}


///////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

FQUALIFIER
EInside GPTubs_Inside( GEOMETRYLOC GPTubs *This, 
		       GPThreeVector p )
{
  G4double r2,pPhi,tolRMin,tolRMax;
  EInside in = kOutside ;

  const G4double halfCarTolerance=kCarTolerance*0.5;
  const G4double halfRadTolerance=kRadTolerance*0.5;
  const G4double halfAngTolerance=kAngTolerance*0.5;

  if (fabs(p.z) <= This->fDz - halfCarTolerance)
  {
    r2 = p.x*p.x + p.y*p.y ;

    if (This->fRMin) { tolRMin = This->fRMin + halfRadTolerance ; }
    else       { tolRMin = 0 ; }

    tolRMax = This->fRMax - halfRadTolerance ;
      
    if ((r2 >= tolRMin*tolRMin) && (r2 <= tolRMax*tolRMax))
    {
      if ( This->fPhiFullTube )
      {
        in = kInside ;
      }
      else
      {
        // Try inner tolerant phi boundaries (=>inside)
        // if not inside, try outer tolerant phi boundaries

        if ( (tolRMin==0) && (fabs(p.x)<=halfCarTolerance)
                          && (fabs(p.y)<=halfCarTolerance) )
        {
          in=kSurface;
        }
        else
        {
          pPhi = atan2(p.y,p.x) ;

          if ( pPhi < -halfAngTolerance )  { pPhi += twopi; } // 0<=pPhi<2pi

          if ( This->fSPhi >= 0 )
          {

            if ( (fabs(pPhi) < halfAngTolerance)
              && (fabs(This->fSPhi + This->fDPhi - twopi) < halfAngTolerance) )
            { 
              pPhi += twopi ; // 0 <= pPhi < 2pi
            }
            if ( (pPhi >= This->fSPhi + halfAngTolerance)
              && (pPhi <= This->fSPhi + This->fDPhi - halfAngTolerance) )
            {
              in = kInside ;
            }
            else if ( (pPhi >= This->fSPhi - halfAngTolerance)
                   && (pPhi <= This->fSPhi + This->fDPhi + halfAngTolerance) )
            {
              in = kSurface ;
            }
          }
          else  // This->fSPhi < 0
          {
            if ( (pPhi <= This->fSPhi + twopi - halfAngTolerance)
              && (pPhi >= This->fSPhi + This->fDPhi  + halfAngTolerance) ) 
	    { ; } //kOutside
            else if ( (pPhi <= This->fSPhi + twopi + halfAngTolerance)
                   && (pPhi >= This->fSPhi + This->fDPhi  - halfAngTolerance) )
            {
              in = kSurface ;
            }
            else
            {
              in = kInside ;
            }
          }
        }                    
      }
    }
    else  // Try generous boundaries
    {
      tolRMin = This->fRMin - halfRadTolerance ;
      tolRMax = This->fRMax + halfRadTolerance ;

      if ( tolRMin < 0 )  { tolRMin = 0; }

      if ( (r2 >= tolRMin*tolRMin) && (r2 <= tolRMax*tolRMax) )
      {
        if (This->fPhiFullTube || (r2 <=halfRadTolerance*halfRadTolerance) )
        {                        // Continuous in phi or on z-axis
          in = kSurface ;
        }
        else // Try outer tolerant phi boundaries only
        {
          pPhi = atan2(p.y,p.x) ;

          if ( pPhi < -halfAngTolerance)  { pPhi += twopi; } // 0<=pPhi<2pi
          if ( This->fSPhi >= 0 )
          {
            if ( (fabs(pPhi) < halfAngTolerance)
              && (fabs(This->fSPhi + This->fDPhi - twopi) < halfAngTolerance) )
            { 
              pPhi += twopi ; // 0 <= pPhi < 2pi
            }
            if ( (pPhi >= This->fSPhi - halfAngTolerance)
              && (pPhi <= This->fSPhi + This->fDPhi + halfAngTolerance) )
            {
              in = kSurface ;
            }
          }
          else  // This->fSPhi < 0
          {
            if ( (pPhi <= This->fSPhi + twopi - halfAngTolerance)
              && (pPhi >= This->fSPhi + This->fDPhi + halfAngTolerance) ) {;} // kOutside
            else
            {
              in = kSurface ;
            }
          }
        }
      }
    }
  }
  else if (fabs(p.z) <= This->fDz + halfCarTolerance)
  {                                          // Check within tolerant r limits
    r2      = p.x*p.x + p.y*p.y ;
    tolRMin = This->fRMin - halfRadTolerance ;
    tolRMax = This->fRMax + halfRadTolerance ;

    if ( tolRMin < 0 )  { tolRMin = 0; }

    if ( (r2 >= tolRMin*tolRMin) && (r2 <= tolRMax*tolRMax) )
    {
      if (This->fPhiFullTube || (r2 <=halfRadTolerance*halfRadTolerance))
      {                        // Continuous in phi or on z-axis
        in = kSurface ;
      }
      else // Try outer tolerant phi boundaries
      {
        pPhi = atan2(p.y,p.x) ;

        if ( pPhi < -halfAngTolerance )  { pPhi += twopi; }  // 0<=pPhi<2pi
        if ( This->fSPhi >= 0 )
        {
          if ( (fabs(pPhi) < halfAngTolerance)
            && (fabs(This->fSPhi + This->fDPhi - twopi) < halfAngTolerance) )
          { 
            pPhi += twopi ; // 0 <= pPhi < 2pi
          }
          if ( (pPhi >= This->fSPhi - halfAngTolerance)
            && (pPhi <= This->fSPhi + This->fDPhi + halfAngTolerance) )
          {
            in = kSurface;
          }
        }
        else  // This->fSPhi < 0
        {
          if ( (pPhi <= This->fSPhi + twopi - halfAngTolerance)
            && (pPhi >= This->fSPhi + This->fDPhi  + halfAngTolerance) ) {;}
          else
          {
            in = kSurface ;
          }
        }      
      }
    }
  }
  return in;
}

///////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

FQUALIFIER
GPThreeVector GPTubs_SurfaceNormal( GEOMETRYLOC GPTubs *This,
				    GPThreeVector p)

{
  G4int noSurfaces = 0;
  G4double rho, pPhi;
  G4double distZ, distRMin, distRMax;
  G4double distSPhi = kInfinity, distEPhi = kInfinity;

  const G4double halfCarTolerance = 0.5*kCarTolerance;
  const G4double halfAngTolerance = 0.5*kAngTolerance;

  GPThreeVector norm = GPThreeVector_create(0.,0.,0.);
  GPThreeVector sumnorm = GPThreeVector_create(0.,0.,0.);
  GPThreeVector nZ  = GPThreeVector_create(0, 0, 1.0);
  GPThreeVector nR  = GPThreeVector_create(0.,0.,0.); 
  GPThreeVector nPs = GPThreeVector_create(0.,0.,0.);
  GPThreeVector nPe = GPThreeVector_create(0.,0.,0.);

  rho = sqrt(p.x*p.x + p.y*p.y);

  distRMin = fabs(rho - This->fRMin);
  distRMax = fabs(rho - This->fRMax);
  distZ    = fabs(fabs(p.z) - This->fDz);

  if (!This->fPhiFullTube )    // Protected against (0,0,z) 
  {
    if ( rho > halfCarTolerance )
    {
      pPhi = atan2(p.y,p.x);
    
      if(pPhi  < This->fSPhi- halfCarTolerance)           { pPhi += twopi; }
      else if(pPhi > This->fSPhi+This->fDPhi+ halfCarTolerance) { pPhi -= twopi; }

      distSPhi = fabs(pPhi - This->fSPhi);       
      distEPhi = fabs(pPhi - This->fSPhi - This->fDPhi); 
    }
    else if( !This->fRMin )
    {
      distSPhi = 0.; 
      distEPhi = 0.; 
    }
    nPs = GPThreeVector_create(sin(This->fSPhi),-cos(This->fSPhi),0);
    nPe = GPThreeVector_create(-sin(This->fSPhi+This->fDPhi),cos(This->fSPhi+This->fDPhi),0);
  }
  if ( rho > halfCarTolerance ) { nR = GPThreeVector_create(p.x/rho,p.y/rho,0); }

  if( distRMax <= halfCarTolerance )
  {
    noSurfaces ++;
    sumnorm = GPThreeVector_add(sumnorm,nR);
  }
  if( This->fRMin && (distRMin <= halfCarTolerance) )
  {
    noSurfaces ++;
    sumnorm = GPThreeVector_sub(sumnorm,nR);
  }
  if( This->fDPhi < twopi )   
  {
    if (distSPhi <= halfAngTolerance)  
    {
      noSurfaces ++;
      sumnorm = GPThreeVector_add(sumnorm,nPs);
    }
    if (distEPhi <= halfAngTolerance)  
    {
      noSurfaces ++;
      sumnorm = GPThreeVector_add(sumnorm,nPe);
    }
  }
  if (distZ <= halfCarTolerance)  
  {
    noSurfaces ++;
    if ( p.z >= 0.)  { sumnorm = GPThreeVector_add(sumnorm,nZ); }
    else             { sumnorm = GPThreeVector_sub(sumnorm,nZ); }
  }
  if ( noSurfaces == 0 )
  {
    norm = GPTubs_ApproxSurfaceNormal(This,p);
  }
  else if ( noSurfaces == 1 )  { norm = sumnorm; }
  else                         { norm = GPThreeVector_unit(sumnorm); }

  return norm;
}

/////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

FQUALIFIER
GPThreeVector GPTubs_ApproxSurfaceNormal( GEOMETRYLOC GPTubs *This, 
					  GPThreeVector p )
{
  //  ENorm side ;
  enum {kNRMin,kNRMax,kNSPhi,kNEPhi,kNZ} side;
  GPThreeVector norm = GPThreeVector_create(0,0,0);
  G4double rho, phi ;
  G4double distZ, distRMin, distRMax, distSPhi, distEPhi, distMin ;

  rho = sqrt(p.x*p.x + p.y*p.y) ;

  distRMin = fabs(rho - This->fRMin) ;
  distRMax = fabs(rho - This->fRMax) ;
  distZ    = fabs(fabs(p.z) - This->fDz) ;

  if (distRMin < distRMax) // First minimum
  {
    if ( distZ < distRMin )
    {
       distMin = distZ ;
       side    = kNZ ;
    }
    else
    {
      distMin = distRMin ;
      side    = kNRMin   ;
    }
  }
  else
  {
    if ( distZ < distRMax )
    {
      distMin = distZ ;
      side    = kNZ   ;
    }
    else
    {
      distMin = distRMax ;
      side    = kNRMax   ;
    }
  }   
  if (!This->fPhiFullTube  &&  rho ) // Protected against (0,0,z) 
  {
    phi = atan2(p.y,p.x) ;

    if ( phi < 0 )  { phi += twopi; }

    if ( This->fSPhi < 0 )
    {
      distSPhi = fabs(phi - (This->fSPhi + twopi))*rho ;
    }
    else
    {
      distSPhi = fabs(phi - This->fSPhi)*rho ;
    }
    distEPhi = fabs(phi - This->fSPhi - This->fDPhi)*rho ;
                                      
    if (distSPhi < distEPhi) // Find new minimum
    {
      if ( distSPhi < distMin )
      {
        side = kNSPhi ;
      }
    }
    else
    {
      if ( distEPhi < distMin )
      {
        side = kNEPhi ;
      }
    }
  }    
  switch ( side )
  {
    case kNRMin : // Inner radius
    {                      
      norm = GPThreeVector_create(-p.x/rho, -p.y/rho, 0) ;
      break ;
    }
    case kNRMax : // Outer radius
    {                  
      norm = GPThreeVector_create(p.x/rho, p.y/rho, 0) ;
      break ;
    }
    case kNZ :    // + or - dz
    {                              
      if ( p.z > 0 )  { norm = GPThreeVector_create(0,0,1) ; }
      else              { norm = GPThreeVector_create(0,0,-1); }
      break ;
    }
    case kNSPhi:
    {
      norm = GPThreeVector_create(sin(This->fSPhi), -cos(This->fSPhi), 0) ;
      break ;
    }
    case kNEPhi:
    {
      norm = GPThreeVector_create(-sin(This->fSPhi+This->fDPhi), cos(This->fSPhi+This->fDPhi), 0) ;
      break;
    }
    default:      // Should never reach this case ...
    {
      ;
      //      DumpInfo();
      //      G4Exception("GPTubs_ApproxSurfaceNormal()",
      //                  "GeomSolids1002", JustWarning,
      //                  "Undefined side for valid surface normal to solid.");
      break ;
    }    
  }                
  return norm;
}

////////////////////////////////////////////////////////////////////
//
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes 
//        - if at valid r, phi, return
//
// -> If point is outer outer radius, compute intersection with rmax
//        - if at valid phi,z return
//
// -> Compute intersection with inner radius, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - 'if valid' implies tolerant checking of intersection points

FQUALIFIER
G4double GPTubs_DistanceToIn2( GEOMETRYLOC GPTubs *This,
			       GPThreeVector p,
			       GPThreeVector v  )
{
  G4double snxt = kInfinity ;      // snxt = default return value
  G4double tolORMin2, tolIRMax2 ;  // 'generous' radii squared
  G4double tolORMax2, tolIRMin2, tolODz, tolIDz ;
  //  const G4double dRmax = 100.*This->fRMax;

  const G4double halfCarTolerance = 0.5*kCarTolerance;
  const G4double halfRadTolerance = 0.5*kRadTolerance;

  // Intersection point variables
  //
  G4double Dist, s, xi, yi, zi, rho2, inum, iden, cosPsi, Comp ;
  G4double t1, t2, t3, b, c, d ;     // Quadratic solver variables 
  
  // Calculate tolerant rmin and rmax

  if (This->fRMin > kRadTolerance)
  {
    tolORMin2 = (This->fRMin - halfRadTolerance)*(This->fRMin - halfRadTolerance) ;
    tolIRMin2 = (This->fRMin + halfRadTolerance)*(This->fRMin + halfRadTolerance) ;
  }
  else
  {
    tolORMin2 = 0.0 ;
    tolIRMin2 = 0.0 ;
  }
  tolORMax2 = (This->fRMax + halfRadTolerance)*(This->fRMax + halfRadTolerance) ;
  tolIRMax2 = (This->fRMax - halfRadTolerance)*(This->fRMax - halfRadTolerance) ;

  // Intersection with Z surfaces

  tolIDz = This->fDz - halfCarTolerance ;
  tolODz = This->fDz + halfCarTolerance ;

  if (fabs(p.z) >= tolIDz)
  {
    if ( p.z*v.z < 0 )    // at +Z going in -Z or visa versa
    {
      s = (fabs(p.z) - This->fDz)/fabs(v.z) ;   // Z intersect distance

      if(s < 0.0)  { s = 0.0; }

      xi   = p.x + s*v.x ;                // Intersection coords
      yi   = p.y + s*v.y ;
      rho2 = xi*xi + yi*yi ;

      // Check validity of intersection

      if ((tolIRMin2 <= rho2) && (rho2 <= tolIRMax2))
      {
	if (!This->fPhiFullTube && rho2)
        {
          // Psi = angle made with central (average) phi of shape
          //
          inum   = xi*This->cosCPhi + yi*This->sinCPhi ;
          iden   = sqrt(rho2) ;
          cosPsi = inum/iden ;
          if (cosPsi >= This->cosHDPhiIT)  { return s ; }
        }
        else
        {
          return s ;
        }
      }
    }
    else
    {
      if ( snxt<halfCarTolerance )  { snxt=0; }
      return snxt ;  // On/outside extent, and heading away
                     // -> cannot intersect
    }
  }

  // -> Can not intersect z surfaces
  //
  // Intersection with rmax (possible return) and rmin (must also check phi)
  //
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=R^2
  //
  // Hence (v.x^2+v.y^2)t^2+ 2t(p.x*v.x+p.y*v.y)+p.x^2+p.y^2-R^2=0
  //            t1                t2                t3

  t1 = 1.0 - v.z*v.z ;
  t2 = p.x*v.x + p.y*v.y ;
  t3 = p.x*p.x + p.y*p.y ;

  if ( t1 > 0 )        // Check not || to z axis
  {
    b = t2/t1 ;
    c = t3 - This->fRMax*This->fRMax ;
    if ((t3 >= tolORMax2) && (t2<0))   // This also handles the tangent case
    {
      // Try outer cylinder intersection
      //          c=(t3-This->fRMax*This->fRMax)/t1;

      c /= t1 ;
      d = b*b - c ;

      if (d >= 0)  // If real root
      {
        s = c/(-b+sqrt(d));
        if (s >= 0)  // If 'forwards'
        {
	  /*
          if ( s>dRmax ) // Avoid rounding errors due to precision issues on
          {              // 64 bits systems. Split long distances and recompute
            G4double fTerm = s-fmod(s,dRmax);
            s = fTerm + DistanceToIn(p+fTerm*v,v);
          } 
	  */
          // Check z intersection
          //
          zi = p.z + s*v.z ;
          if (fabs(zi)<=tolODz)
          {
            // Z ok. Check phi intersection if reqd
            //
            if (This->fPhiFullTube)
            {
              return s ;
            }
            else
            {
              xi     = p.x + s*v.x ;
              yi     = p.y + s*v.y ;
              cosPsi = (xi*This->cosCPhi + yi*This->sinCPhi)/This->fRMax ;
              if (cosPsi >= This->cosHDPhiIT)  { return s ; }
            }
          }  //  end if fabs(zi)
        }    //  end if (s>=0)
      }      //  end if (d>=0)
    }        //  end if (r>=This->fRMax)
    else 
    {
      // Inside outer radius :
      // check not inside, and heading through tubs (-> 0 to in)

      if ((t3 > tolIRMin2) && (t2 < 0) && (fabs(p.z) <= tolIDz))
      {
        // Inside both radii, delta r -ve, inside z extent

        if (!This->fPhiFullTube)
        {
          inum   = p.x*This->cosCPhi + p.y*This->sinCPhi ;
          iden   = sqrt(t3) ;
          cosPsi = inum/iden ;
          if (cosPsi >= This->cosHDPhiIT)
          {
            // In the old version, the small negative tangent for the point
            // on surface was not taken in account, and returning 0.0 ...
            // New version: check the tangent for the point on surface and 
            // if no intersection, return kInfinity, if intersection instead
            // return s.
            //
            c = t3-This->fRMax*This->fRMax; 
            if ( c<=0.0 )
            {
              return 0.0;
            }
            else
            {
              c = c/t1 ;
              d = b*b-c;
              if ( d>=0.0 )
              {
                snxt = c/(-b+sqrt(d)); // using safe solution
                                            // for quadratic equation 
                if ( snxt < halfCarTolerance ) { snxt=0; }
                return snxt ;
              }      
              else
              {
                return kInfinity;
              }
            }
          } 
        }
        else
        {   
          // In the old version, the small negative tangent for the point
          // on surface was not taken in account, and returning 0.0 ...
          // New version: check the tangent for the point on surface and 
          // if no intersection, return kInfinity, if intersection instead
          // return s.
          //
          c = t3 - This->fRMax*This->fRMax; 
          if ( c<=0.0 )
          {
            return 0.0;
          }
          else
          {
            c = c/t1 ;
            d = b*b-c;
            if ( d>=0.0 )
            {
              snxt= c/(-b+sqrt(d)); // using safe solution
                                         // for quadratic equation 
              if ( snxt < halfCarTolerance ) { snxt=0; }
              return snxt ;
            }      
            else
            {
              return kInfinity;
            }
          }
        } // end if   (!This->fPhiFullTube)
      }   // end if   (t3>tolIRMin2)
    }     // end if   (Inside Outer Radius) 
    if ( This->fRMin )    // Try inner cylinder intersection
    {
      c = (t3 - This->fRMin*This->fRMin)/t1 ;
      d = b*b - c ;
      if ( d >= 0.0 )  // If real root
      {
        // Always want 2nd root - we are outside and know rmax Hit was bad
        // - If on surface of rmin also need farthest root

        s =( b > 0. )? c/(-b - sqrt(d)) : (-b + sqrt(d));
        if (s >= -halfCarTolerance)  // check forwards
        {
          // Check z intersection
          //
          if(s < 0.0)  { s = 0.0; }
	  /*
          if ( s>dRmax ) // Avoid rounding errors due to precision issues seen
          {              // 64 bits systems. Split long distances and recompute
            G4double fTerm = s-fmod(s,dRmax);
            s = fTerm + DistanceToIn(p+fTerm*v,v);
          } 
	  */
          zi = p.z + s*v.z ;
          if (fabs(zi) <= tolODz)
          {
            // Z ok. Check phi
            //
            if ( This->fPhiFullTube )
            {
              return s ; 
            }
            else
            {
              xi     = p.x + s*v.x ;
              yi     = p.y + s*v.y ;
              cosPsi = (xi*This->cosCPhi + yi*This->sinCPhi)/This->fRMin ;
              if (cosPsi >= This->cosHDPhiIT)
              {
                // Good inner radius isect
                // - but earlier phi isect still possible

                snxt = s ;
              }
            }
          }        //    end if fabs(zi)
        }          //    end if (s>=0)
      }            //    end if (d>=0)
    }              //    end if (This->fRMin)
  }

  // Phi segment intersection
  //
  // o Tolerant of points inside phi planes by up to kCarTolerance*0.5
  //
  // o NOTE: Large duplication of code between sphi & ephi checks
  //         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
  //            intersection check <=0 -> >=0
  //         -> use some form of loop Construct ?
  //
  if ( !This->fPhiFullTube )
  {
    // First phi surface (Starting phi)
    //
    Comp    = v.x*This->sinSPhi - v.y*This->cosSPhi ;
                    
    if ( Comp < 0 )  // Component in outwards normal dirn
    {
      Dist = (p.y*This->cosSPhi - p.x*This->sinSPhi) ;

      if ( Dist < halfCarTolerance )
      {
        s = Dist/Comp ;

        if (s < snxt)
        {
          if ( s < 0 )  { s = 0.0; }
          zi = p.z + s*v.z ;
          if ( fabs(zi) <= tolODz )
          {
            xi   = p.x + s*v.x ;
            yi   = p.y + s*v.y ;
            rho2 = xi*xi + yi*yi ;

            if ( ( (rho2 >= tolIRMin2) && (rho2 <= tolIRMax2) )
              || ( (rho2 >  tolORMin2) && (rho2 <  tolIRMin2)
                && ( v.y*This->cosSPhi - v.x*This->sinSPhi >  0 )
                && ( v.x*This->cosSPhi + v.y*This->sinSPhi >= 0 )     )
              || ( (rho2 > tolIRMax2) && (rho2 < tolORMax2)
                && (v.y*This->cosSPhi - v.x*This->sinSPhi > 0)
                && (v.x*This->cosSPhi + v.y*This->sinSPhi < 0) )    )
            {
              // z and r intersections good
              // - check intersecting with correct half-plane
              //
              if ((yi*This->cosCPhi-xi*This->sinCPhi) <= halfCarTolerance) { snxt = s; }
            }
          }
        }
      }    
    }
      
    // Second phi surface (Ending phi)

    Comp    = -(v.x*This->sinEPhi - v.y*This->cosEPhi) ;
        
    if (Comp < 0 )  // Component in outwards normal dirn
    {
      Dist = -(p.y*This->cosEPhi - p.x*This->sinEPhi) ;

      if ( Dist < halfCarTolerance )
      {
        s = Dist/Comp ;

        if (s < snxt)
        {
          if ( s < 0 )  { s = 0; }
          zi = p.z + s*v.z ;
          if ( fabs(zi) <= tolODz )
          {
            xi   = p.x + s*v.x ;
            yi   = p.y + s*v.y ;
            rho2 = xi*xi + yi*yi ;
            if ( ( (rho2 >= tolIRMin2) && (rho2 <= tolIRMax2) )
                || ( (rho2 > tolORMin2)  && (rho2 < tolIRMin2)
                  && (v.x*This->sinEPhi - v.y*This->cosEPhi >  0)
                  && (v.x*This->cosEPhi + v.y*This->sinEPhi >= 0) )
                || ( (rho2 > tolIRMax2) && (rho2 < tolORMax2)
                  && (v.x*This->sinEPhi - v.y*This->cosEPhi > 0)
                  && (v.x*This->cosEPhi + v.y*This->sinEPhi < 0) ) )
            {
              // z and r intersections good
              // - check intersecting with correct half-plane
              //
              if ( (yi*This->cosCPhi-xi*This->sinCPhi) >= 0 ) { snxt = s; }
            }                         //?? >=-halfCarTolerance
          }
        }
      }
    }         //  Comp < 0
  }           //  !This->fPhiFullTube 
  if ( snxt<halfCarTolerance )  { snxt=0; }
  return snxt ;
}
 
//////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes 
//        - if at valid r, phi, return
//
// -> If point is outer outer radius, compute intersection with rmax
//        - if at valid phi,z return
//
// -> Compute intersection with inner radius, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - Precalculations for phi trigonometry are Done `just in time'
// - `if valid' implies tolerant checking of intersection points
//   Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

FQUALIFIER
G4double GPTubs_DistanceToIn( GEOMETRYLOC GPTubs *This,
			      GPThreeVector p )
{
  G4double safe=0.0, rho, safe1, safe2, safe3 ;
  G4double safePhi, cosPsi ;

  rho   = sqrt(p.x*p.x + p.y*p.y) ;
  safe1 = This->fRMin - rho ;
  safe2 = rho - This->fRMax ;
  safe3 = fabs(p.z) - This->fDz ;

  if ( safe1 > safe2 ) { safe = safe1; }
  else                 { safe = safe2; }
  if ( safe3 > safe )  { safe = safe3; }

  if ( (!This->fPhiFullTube) && rho )
  {
    // Psi=angle from central phi to point
    //
    cosPsi = (p.x*This->cosCPhi + p.y*This->sinCPhi)/rho ;
    
    if ( cosPsi < cos(This->fDPhi*0.5) )
    {
      // Point lies outside phi range

      if ( (p.y*This->cosCPhi - p.x*This->sinCPhi) <= 0 )
      {
        safePhi = fabs(p.x*This->sinSPhi - p.y*This->cosSPhi) ;
      }
      else
      {
        safePhi = fabs(p.x*This->sinEPhi - p.y*This->cosEPhi) ;
      }
      if ( safePhi > safe )  { safe = safePhi; }
    }
  }
  if ( safe < 0 )  { safe = 0; }
  return safe ;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

FQUALIFIER
G4double GPTubs_DistanceToOut2( GEOMETRYLOC GPTubs *This,
				GPThreeVector p,
				GPThreeVector v,
				const G4bool calcNorm,
				G4bool *validNorm,
				GPThreeVector *n    )
{  
  enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ};

  ESide side=kNull , sider=kNull, sidephi=kNull ;
  G4double snxt, sr=kInfinity, sphi=kInfinity, pdist ;
  G4double deltaR, t1, t2, t3, b, c, d2, roMin2 ;

  const G4double halfCarTolerance = kCarTolerance*0.5;
  const G4double halfAngTolerance = kAngTolerance*0.5;
 
  // Vars for phi intersection:

  G4double pDistS, compS, pDistE, compE, sphi2, xi, yi, vphi, roi2 ;
 
  // Z plane intersection

  if (v.z > 0 )
  {
    pdist = This->fDz - p.z ;
    if ( pdist > halfCarTolerance )
    {
      snxt = pdist/v.z ;
      side = kPZ ;
    }
    else
    {
      if (calcNorm)
      {
        *n         = GPThreeVector_create(0,0,1) ;
        *validNorm = true ;
      }
      return snxt = 0 ;
    }
  }
  else if ( v.z < 0 )
  {
    pdist = This->fDz + p.z ;

    if ( pdist > halfCarTolerance )
    {
      snxt = -pdist/v.z ;
      side = kMZ ;
    }
    else
    {
      if (calcNorm)
      {
        *n         = GPThreeVector_create(0,0,-1) ;
        *validNorm = true ;
      }
      return snxt = 0.0 ;
    }
  }
  else
  {
    snxt = kInfinity ;    // Travel perpendicular to z axis
    side = kNull;
  }

  // Radial Intersections
  //
  // Find intersection with cylinders at rmax/rmin
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=R^2
  //
  // Hence (v.x^2+v.y^2)t^2+ 2t(p.x*v.x+p.y*v.y)+p.x^2+p.y^2-R^2=0
  //
  //            t1                t2                    t3

  t1   = 1.0 - v.z*v.z ;      // since v normalised
  t2   = p.x*v.x + p.y*v.y ;
  t3   = p.x*p.x + p.y*p.y ;

  if ( snxt > 10*(This->fDz+This->fRMax) )  { roi2 = 2*This->fRMax*This->fRMax; }
  else  { roi2 = snxt*snxt*t1 + 2*snxt*t2 + t3; }        // radius^2 on +-This->fDz

  if ( t1 > 0 ) // Check not parallel
  {
    // Calculate sr, r exit distance
     
    if ( (t2 >= 0.0) && (roi2 > This->fRMax*(This->fRMax + kRadTolerance)) )
    {
      // Delta r not negative => leaving via rmax

      deltaR = t3 - This->fRMax*This->fRMax ;

      // NOTE: Should use rho-This->fRMax<-kRadTolerance*0.5
      // - avoid sqrt for efficiency

      if ( deltaR < -kRadTolerance*This->fRMax )
      {
        b     = t2/t1 ;
        c     = deltaR/t1 ;
        d2    = b*b-c;
        if( d2 >= 0 ) { sr = c/( -b - sqrt(d2)); }
        else          { sr = 0.; }
        sider = kRMax ;
      }
      else
      {
        // On tolerant boundary & heading outwards (or perpendicular to)
        // outer radial surface -> leaving immediately

        if ( calcNorm ) 
        {
          *n         = GPThreeVector_create(p.x/This->fRMax,p.y/This->fRMax,0) ;
          *validNorm = true ;
        }
        return snxt = 0 ; // Leaving by rmax immediately
      }
    }             
    else if ( t2 < 0. ) // i.e.  t2 < 0; Possible rmin intersection
    {
      roMin2 = t3 - t2*t2/t1 ; // min ro2 of the plane of movement 

      if ( This->fRMin && (roMin2 < This->fRMin*(This->fRMin - kRadTolerance)) )
      {
        deltaR = t3 - This->fRMin*This->fRMin ;
        b      = t2/t1 ;
        c      = deltaR/t1 ;
        d2     = b*b - c ;

        if ( d2 >= 0 )   // Leaving via rmin
        {
          // NOTE: SHould use rho-rmin>kRadTolerance*0.5
          // - avoid sqrt for efficiency

          if (deltaR > kRadTolerance*This->fRMin)
          {
            sr = c/(-b+sqrt(d2)); 
            sider = kRMin ;
          }
          else
          {
            if ( calcNorm ) { *validNorm = false; }  // Concave side
            return snxt = 0.0;
          }
        }
        else    // No rmin intersect -> must be rmax intersect
        {
          deltaR = t3 - This->fRMax*This->fRMax ;
          c     = deltaR/t1 ;
          d2    = b*b-c;
          if( d2 >=0. )
          {
            sr     = -b + sqrt(d2) ;
            sider  = kRMax ;
          }
          else // Case: On the border+t2<kRadTolerance
               //       (v is perpendicular to the surface)
          {
            if (calcNorm)
            {
              *n = GPThreeVector_create(p.x/This->fRMax,p.y/This->fRMax,0) ;
              *validNorm = true ;
            }
            return snxt = 0.0;
          }
        }
      }
      else if ( roi2 > This->fRMax*(This->fRMax + kRadTolerance) )
           // No rmin intersect -> must be rmax intersect
      {
        deltaR = t3 - This->fRMax*This->fRMax ;
        b      = t2/t1 ;
        c      = deltaR/t1;
        d2     = b*b-c;
        if( d2 >= 0 )
        {
          sr     = -b + sqrt(d2) ;
          sider  = kRMax ;
        }
        else // Case: On the border+t2<kRadTolerance
             //       (v is perpendicular to the surface)
        {
          if (calcNorm)
          {
            *n = GPThreeVector_create(p.x/This->fRMax,p.y/This->fRMax,0) ;
            *validNorm = true ;
          }
          return snxt = 0.0;
        }
      }
    }
    
    // Phi Intersection

    if ( !This->fPhiFullTube )
    {
      // add angle calculation with correction 
      // of the difference in domain of atan2 and Sphi
      //
      vphi = atan2(v.y,v.x) ;
     
      if ( vphi < This->fSPhi - halfAngTolerance  )             { vphi += twopi; }
      else if ( vphi > This->fSPhi + This->fDPhi + halfAngTolerance ) { vphi -= twopi; }


      if ( p.x || p.y )  // Check if on z axis (rho not needed later)
      {
        // pDist -ve when inside

        pDistS = p.x*This->sinSPhi - p.y*This->cosSPhi ;
        pDistE = -p.x*This->sinEPhi + p.y*This->cosEPhi ;

        // Comp -ve when in direction of outwards normal

        compS   = -This->sinSPhi*v.x + This->cosSPhi*v.y ;
        compE   =  This->sinEPhi*v.x - This->cosEPhi*v.y ;
       
        sidephi = kNull;
        
        if( ( (This->fDPhi <= pi) && ( (pDistS <= halfCarTolerance)
                              && (pDistE <= halfCarTolerance) ) )
         || ( (This->fDPhi >  pi) && !((pDistS >  halfCarTolerance)
                              && (pDistE >  halfCarTolerance) ) )  )
        {
          // Inside both phi *full* planes
          
          if ( compS < 0 )
          {
            sphi = pDistS/compS ;
            
            if (sphi >= -halfCarTolerance)
            {
              xi = p.x + sphi*v.x ;
              yi = p.y + sphi*v.y ;
              
              // Check intersecting with correct half-plane
              // (if not -> no intersect)
              //
              if( (fabs(xi)<=kCarTolerance)&&(fabs(yi)<=kCarTolerance) )
              {
                sidephi = kSPhi;
                if (((This->fSPhi-halfAngTolerance)<=vphi)
                   &&((This->fSPhi+This->fDPhi+halfAngTolerance)>=vphi))
                {
                  sphi = kInfinity;
                }
              }
              else if ( yi*This->cosCPhi-xi*This->sinCPhi >=0 )
              {
                sphi = kInfinity ;
              }
              else
              {
                sidephi = kSPhi ;
                if ( pDistS > -halfCarTolerance )
                {
                  sphi = 0.0 ; // Leave by sphi immediately
                }    
              }       
            }
            else
            {
              sphi = kInfinity ;
            }
          }
          else
          {
            sphi = kInfinity ;
          }

          if ( compE < 0 )
          {
            sphi2 = pDistE/compE ;
            
            // Only check further if < starting phi intersection
            //
            if ( (sphi2 > -halfCarTolerance) && (sphi2 < sphi) )
            {
              xi = p.x + sphi2*v.x ;
              yi = p.y + sphi2*v.y ;
              
              if ((fabs(xi)<=kCarTolerance)&&(fabs(yi)<=kCarTolerance))
              {
                // Leaving via ending phi
                //
                if( !((This->fSPhi-halfAngTolerance <= vphi)
                     &&(This->fSPhi+This->fDPhi+halfAngTolerance >= vphi)) )
                {
                  sidephi = kEPhi ;
                  if ( pDistE <= -halfCarTolerance )  { sphi = sphi2 ; }
                  else                                { sphi = 0.0 ;   }
                }
              } 
              else    // Check intersecting with correct half-plane 

              if ( (yi*This->cosCPhi-xi*This->sinCPhi) >= 0)
              {
                // Leaving via ending phi
                //
                sidephi = kEPhi ;
                if ( pDistE <= -halfCarTolerance ) { sphi = sphi2 ; }
                else                               { sphi = 0.0 ;   }
              }
            }
          }
        }
        else
        {
          sphi = kInfinity ;
        }
      }
      else
      {
        // On z axis + travel not || to z axis -> if phi of vector direction
        // within phi of shape, Step limited by rmax, else Step =0
               
        if ( (This->fSPhi - halfAngTolerance <= vphi)
           && (vphi <= This->fSPhi + This->fDPhi + halfAngTolerance ) )
        {
          sphi = kInfinity ;
        }
        else
        {
          sidephi = kSPhi ; // arbitrary 
          sphi    = 0.0 ;
        }
      }
      if (sphi < snxt)  // Order intersecttions
      {
        snxt = sphi ;
        side = sidephi ;
      }
    }
    if (sr < snxt)  // Order intersections
    {
      snxt = sr ;
      side = sider ;
    }
  }
  if (calcNorm)
  {
    switch(side)
    {
      case kRMax:
        // Note: returned vector not normalised
        // (divide by This->fRMax for unit vector)
        //
        xi = p.x + snxt*v.x ;
        yi = p.y + snxt*v.y ;
        *n = GPThreeVector_create(xi/This->fRMax,yi/This->fRMax,0) ;
        *validNorm = true ;
        break ;

      case kRMin:
        *validNorm = false ;  // Rmin is inconvex
        break ;

      case kSPhi:
        if ( This->fDPhi <= pi )
        {
          *n         = GPThreeVector_create(This->sinSPhi,-This->cosSPhi,0) ;
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;

      case kEPhi:
        if (This->fDPhi <= pi)
        {
          *n = GPThreeVector_create(-This->sinEPhi,This->cosEPhi,0) ;
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;

      case kPZ:
        *n         = GPThreeVector_create(0,0,1) ;
        *validNorm = true ;
        break ;

      case kMZ:
        *n         = GPThreeVector_create(0,0,-1) ;
        *validNorm = true ;
        break ;

      default:
	;
	/*
        G4cout << G4endl ;
        DumpInfo();
        ostringstream message;
        G4int oldprc = message.precision(16);
        message << "Undefined side for valid surface normal to solid."
                << G4endl
                << "Position:"  << G4endl << G4endl
                << "p.x = "   << p.x/mm << " mm" << G4endl
                << "p.y = "   << p.y/mm << " mm" << G4endl
                << "p.z = "   << p.z/mm << " mm" << G4endl << G4endl
                << "Direction:" << G4endl << G4endl
                << "v.x = "   << v.x << G4endl
                << "v.y = "   << v.y << G4endl
                << "v.z = "   << v.z << G4endl << G4endl
                << "Proposed distance :" << G4endl << G4endl
                << "snxt = "    << snxt/mm << " mm" << G4endl ;
        message.precision(oldprc) ;
        G4Exception("GPTubs_DistanceToOut(p,v,..)", "GeomSolids1002",
                    JustWarning, message);
	*/
        break ;
    }
  }
  if ( snxt<halfCarTolerance )  { snxt=0 ; }

  return snxt ;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

FQUALIFIER
G4double GPTubs_DistanceToOut( GEOMETRYLOC GPTubs *This,
			       GPThreeVector p )
{
  G4double safe=0.0, rho, safeR1, safeR2, safeZ, safePhi ;
  rho = sqrt(p.x*p.x + p.y*p.y) ;

  if ( This->fRMin )
  {
    safeR1 = rho   - This->fRMin ;
    safeR2 = This->fRMax - rho ;
 
    if ( safeR1 < safeR2 ) { safe = safeR1 ; }
    else                   { safe = safeR2 ; }
  }
  else
  {
    safe = This->fRMax - rho ;
  }
  safeZ = This->fDz - fabs(p.z) ;

  if ( safeZ < safe )  { safe = safeZ ; }

  // Check if phi divided, Calc distances closest phi plane
  //
  if ( !This->fPhiFullTube )
  {
    if ( p.y*This->cosCPhi-p.x*This->sinCPhi <= 0 )
    {
      safePhi = -(p.x*This->sinSPhi - p.y*This->cosSPhi) ;
    }
    else
    {
      safePhi = (p.x*This->sinEPhi - p.y*This->cosEPhi) ;
    }
    if (safePhi < safe)  { safe = safePhi ; }
  }
  if ( safe < 0 )  { safe = 0 ; }

  return safe ;  
}

/////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Ordering [0-3] -This->fDz cross section
//          [4-7] +This->fDz cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility
//  Potential improvement: For last slice, use actual ending angle
//                         to avoid rounding error problems.

FQUALIFIER
GPThreeVectorList
GPTubs_CreateRotatedVertices(GPTubs *This, 
			     GPAffineTransform pTransform ) 
{
  GPThreeVector vertex0, vertex1, vertex2, vertex3 ;

  G4double meshAngle, meshRMax, crossAngle,
           cosCrossAngle, sinCrossAngle, sAngle;
  G4double rMaxX, rMaxY, rMinX, rMinY, meshRMin ;
  G4int crossSection, noCrossSections;

  // Compute no of cross-sections necessary to mesh tube
  //
  noCrossSections = G4int(This->fDPhi/kMeshAngleDefault) + 1 ;

  if ( noCrossSections < kMinMeshSections )
  {
    noCrossSections = kMinMeshSections ;
  }
  else if (noCrossSections>kMaxMeshSections)
  {
    noCrossSections = kMaxMeshSections ;
  }
  // noCrossSections = 4 ;

  meshAngle = This->fDPhi/(noCrossSections - 1) ;
  // meshAngle = This->fDPhi/(noCrossSections) ;

  meshRMax  = (This->fRMax+100*kCarTolerance)/cos(meshAngle*0.5) ;
  meshRMin = This->fRMin - 100*kCarTolerance ; 
 
  // If complete in phi, set start angle such that mesh will be at This->fRMax
  // on the x axis. Will give better extent calculations when not rotated.

  if (This->fPhiFullTube && (This->fSPhi == 0) )  { sAngle = -meshAngle*0.5 ; }
  else                                { sAngle =  This->fSPhi ; }
    
  //  GPThreeVectorList* vertices ;
  //  vertices = new GPThreeVectorList();
  GPThreeVectorList vertices;
  GPThreeVectorList_Constructor(&vertices);
    
  if ( &vertices )
  {
    //    vertices->reserve(noCrossSections*4);
    for (crossSection = 0 ; crossSection < noCrossSections ; crossSection++ )
    {
      // Compute coordinates of cross section at section crossSection

      crossAngle    = sAngle + crossSection*meshAngle ;
      cosCrossAngle = cos(crossAngle) ;
      sinCrossAngle = sin(crossAngle) ;

      rMaxX = meshRMax*cosCrossAngle ;
      rMaxY = meshRMax*sinCrossAngle ;

      if(meshRMin <= 0.0)
      {
        rMinX = 0.0 ;
        rMinY = 0.0 ;
      }
      else
      {
        rMinX = meshRMin*cosCrossAngle ;
        rMinY = meshRMin*sinCrossAngle ;
      }

      vertex0 = GPThreeVector_create(rMinX,rMinY,-This->fDz) ;
      vertex1 = GPThreeVector_create(rMaxX,rMaxY,-This->fDz) ;
      vertex2 = GPThreeVector_create(rMaxX,rMaxY,+This->fDz) ;
      vertex3 = GPThreeVector_create(rMinX,rMinY,+This->fDz) ;

      GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex0)) ;
      GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex1)) ;
      GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex2)) ;
      GPThreeVectorList_pushback(&vertices,GPAffineTransform_TransformPoint(&pTransform,vertex3)) ;
    }
  }
  else
  {
    ;
    /*
    DumpInfo();
    G4Exception("GPTubs_CreateRotatedVertices()",
                "GeomSolids0003", FatalException,
                "Error in allocation of vertices. Out of memory !");
    */
  }
  return vertices ;
}
// $Id: G4Tubs.icc,v 1.18 2010-09-22 08:31:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 inline definitions file
//
// G4Tubs.icc
//
// Implementation of inline methods of G4Tubs
// --------------------------------------------------------------------

FQUALIFIER
G4double GPTubs_GetInnerRadius (GPTubs *This)
{
  return This->fRMin;
}

FQUALIFIER
G4double GPTubs_GetOuterRadius (GPTubs *This)
{
  return This->fRMax;
}

FQUALIFIER
G4double GPTubs_GetZHalfLength (GPTubs *This)
{
  return This->fDz;
}

FQUALIFIER
G4double GPTubs_GetStartPhiAngle (GPTubs *This)
{
  return This->fSPhi;
}

FQUALIFIER
G4double GPTubs_GetDeltaPhiAngle (GPTubs *This)
{
  return This->fDPhi;
}

FQUALIFIER 
void GPTubs_Initialize(GPTubs *This)
{
  ;
  //  fCubicVolume = 0.;
  //  fSurfaceArea = 0.;
  //  fpPolyhedron = 0;
}

FQUALIFIER 
void GPTubs_InitializeTrigonometry(GPTubs *This)
{
  G4double hDPhi = 0.5*This->fDPhi;                       // half delta phi
  G4double cPhi  = This->fSPhi + hDPhi; 
  G4double ePhi  = This->fSPhi + This->fDPhi;

  This->sinCPhi    = sin(cPhi);
  This->cosCPhi    = cos(cPhi);
  This->cosHDPhiIT = cos(hDPhi - 0.5*kAngTolerance); // inner/outer tol half dphi
  This->cosHDPhiOT = cos(hDPhi + 0.5*kAngTolerance);
  This->sinSPhi = sin(This->fSPhi);
  This->cosSPhi = cos(This->fSPhi);
  This->sinEPhi = sin(ePhi);
  This->cosEPhi = cos(ePhi);
}

FQUALIFIER 
void GPTubs_CheckSPhiAngle(GPTubs *This,G4double sPhi)
{
  // Ensure This->fSphi in 0-2PI or -2PI-0 range if shape crosses 0

  if ( sPhi < 0 )
  {
    This->fSPhi = twopi - fmod(fabs(sPhi),twopi);
  }
  else
  {
    This->fSPhi = fmod(sPhi,twopi) ;
  }
  if ( This->fSPhi+This->fDPhi > twopi )
  {
    This->fSPhi -= twopi ;
  }
}

FQUALIFIER 
void GPTubs_CheckDPhiAngle(GPTubs *This, G4double dPhi)
{
  This->fPhiFullTube = true;
  if ( dPhi >= twopi-kAngTolerance*0.5 )
  {
    This->fDPhi=twopi;
    This->fSPhi=0;
  }
  else
  {
    This->fPhiFullTube = false;
    if ( dPhi > 0 )
    {
      This->fDPhi = dPhi;
    }
    else
    {
      ;
    }
  }
}

FQUALIFIER 
void GPTubs_CheckPhiAngles(GPTubs *This, 
			   G4double sPhi, G4double dPhi)
{
  GPTubs_CheckDPhiAngle(This,dPhi);
  if ( (This->fDPhi<twopi) && (sPhi) ) { GPTubs_CheckSPhiAngle(This,sPhi); }
  GPTubs_InitializeTrigonometry(This);
}

FQUALIFIER
void GPTubs_SetInnerRadius (GPTubs *This,
			    G4double newRMin)
{
  if ( newRMin < 0 ) // Check radii
  {
    ;
  }
  This->fRMin= newRMin;
  GPTubs_Initialize(This);
}

FQUALIFIER
void GPTubs_SetOuterRadius (GPTubs *This,
			    G4double newRMax)
{
  if ( newRMax <= 0 ) // Check radii
  {
  }
  This->fRMax= newRMax;
  GPTubs_Initialize(This);
}

FQUALIFIER
void GPTubs_SetZHalfLength (GPTubs *This,
			    G4double newDz)
{
  if (newDz<=0) // Check z-len
  {
    ;
  }
  This->fDz= newDz;
   GPTubs_Initialize(This);
}

FQUALIFIER
void GPTubs_SetStartPhiAngle (GPTubs *This,
			      G4double newSPhi, G4bool compute)
{
  // Flag 'compute' can be used to explicitely avoid recomputation of
  // trigonometry in case SetDeltaPhiAngle() is invoked afterwards

  GPTubs_CheckSPhiAngle(This,newSPhi);
  //  This->fPhiFullTube = false;
  if (compute)  { GPTubs_InitializeTrigonometry(This); }
  GPTubs_Initialize(This);
}

FQUALIFIER
void GPTubs_SetDeltaPhiAngle (GPTubs *This,
			      G4double newDPhi)
{
  GPTubs_CheckPhiAngles(This,This->fSPhi, newDPhi);
  GPTubs_Initialize(This);
}

//  Older names for access functions

FQUALIFIER
G4double GPTubs_GetRMin (GPTubs *This)
{
  return GPTubs_GetInnerRadius(This);
}

FQUALIFIER
G4double GPTubs_GetRMax (GPTubs *This)
{
  return GPTubs_GetOuterRadius(This);
}

FQUALIFIER
G4double GPTubs_GetDz (GPTubs *This)
{
  return GPTubs_GetZHalfLength(This)  ;
}

FQUALIFIER
G4double GPTubs_GetSPhi (GPTubs *This)
{
  return GPTubs_GetStartPhiAngle(This);
}

FQUALIFIER
G4double GPTubs_GetDPhi (GPTubs *This)
{
  return GPTubs_GetDeltaPhiAngle(This);
}

