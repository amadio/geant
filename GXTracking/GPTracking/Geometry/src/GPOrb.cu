#include "GPOrb.h"
#include "GPUtils.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

FQUALIFIER
void GPOrb_Constructor( GPOrb *This, G4double pRmax )
{
  This->fSolid.fType = kOrb;

  This->fRmax = pRmax;
  This->fRmaxTolerance = GPfmax( kRadTolerance, 2.e-11*This->fRmax);
}

FQUALIFIER
G4bool GPOrb_CalculateExtent(const GPOrb *This,
			     const EAxis pAxis,
			     GPVoxelLimits pVoxelLimit,
			     GPAffineTransform pTransform,
			     G4double* pMin, G4double* pMax)
{
    // Compute x/y/z mins and maxs for bounding box respecting limits,
    // with early returns if outside limits. Then switch() on pAxis,
    // and compute exact x and y limit for x/y case
      
    G4double xoffset,xMin,xMax;
    G4double yoffset,yMin,yMax;
    G4double zoffset,zMin,zMax;

    G4double diff1,diff2,maxDiff,newMin,newMax;
    G4double xoff1,xoff2,yoff1,yoff2;
    
    //    const G4double kCarTolerance = kCarTolerance;

    xoffset=GPAffineTransform_NetTranslation(&pTransform).x;
    xMin=xoffset-This->fRmax;
    xMax=xoffset+This->fRmax;

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
    yoffset=GPAffineTransform_NetTranslation(&pTransform).y;
    yMin=yoffset-This->fRmax;
    yMax=yoffset+This->fRmax;

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
    zoffset=GPAffineTransform_NetTranslation(&pTransform).z;
    zMin=zoffset-This->fRmax;
    zMax=zoffset+This->fRmax;

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

    // Known to cut sphere

    switch (pAxis)
    {
      case kXAxis:
        yoff1=yoffset-yMin;
        yoff2=yMax-yoffset;

        if ( yoff1 >= 0 && yoff2 >= 0 )
        {
          // Y limits cross max/min x => no change
          //
          *pMin=xMin;
          *pMax=xMax;
        }
        else
        {
          // Y limits don't cross max/min x => compute max delta x,
          // hence new mins/maxs
          //
          diff1=sqrt(This->fRmax*This->fRmax-yoff1*yoff1);
          diff2=sqrt(This->fRmax*This->fRmax-yoff2*yoff2);
          maxDiff=(diff1>diff2) ? diff1:diff2;
          newMin=xoffset-maxDiff;
          newMax=xoffset+maxDiff;
          *pMin=(newMin<xMin) ? xMin : newMin;
          *pMax=(newMax>xMax) ? xMax : newMax;
        }
        break;
      case kYAxis:
        xoff1=xoffset-xMin;
        xoff2=xMax-xoffset;
        if (xoff1>=0&&xoff2>=0)
        {
          // X limits cross max/min y => no change
          //
          *pMin=yMin;
          *pMax=yMax;
        }
        else
        {
          // X limits don't cross max/min y => compute max delta y,
          // hence new mins/maxs
          //
          diff1=sqrt(This->fRmax*This->fRmax-xoff1*xoff1);
          diff2=sqrt(This->fRmax*This->fRmax-xoff2*xoff2);
          maxDiff=(diff1>diff2) ? diff1:diff2;
          newMin=yoffset-maxDiff;
          newMax=yoffset+maxDiff;
          *pMin=(newMin<yMin) ? yMin : newMin;
          *pMax=(newMax>yMax) ? yMax : newMax;
        }
        break;
      case kZAxis:
        *pMin=zMin;
        *pMax=zMax;
        break;
      default:
        break;
    }
    *pMin -= This->fRmaxTolerance;
    *pMax += This->fRmaxTolerance;

    return true;  
  
}

// --------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

FQUALIFIER
EInside GPOrb_Inside( GEOMETRYLOC const GPOrb *This, 
		      GPThreeVector p)
{
  G4double rad2,tolRMax;
  EInside in;


  rad2 = GPThreeVector_mag2(p); //p.x*p.x+p.y*p.y+p.z*p.z ;

  G4double rad = sqrt(rad2);

  // G4double rad = std::sqrt(rad2);
  // Check radial surface
  // sets `in'
  
  tolRMax = This->fRmax - This->fRmaxTolerance*0.5 ;
    
  if ( rad <= tolRMax )  { in = kInside ; }
  else
  {
    tolRMax = This->fRmax + This->fRmaxTolerance*0.5 ;       
    if ( rad <= tolRMax )  { in = kSurface ; }
    else                   { in = kOutside ; }
  }
  return in;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If two sides are equidistant, normal of first side (x/y/z) 
// encountered returned

FQUALIFIER
GPThreeVector GPOrb_SurfaceNormal( GEOMETRYLOC const GPOrb *This, 
				   GPThreeVector p)
{
  (void)This;
  /*ENorm side = kNRMax;
  GPThreeVector norm;
  G4double rad = sqrt(p.x*p.x+p.y*p.y+p.z*p.z);

  switch (side)
  {
    case kNRMax: 
      norm = GPThreeVector_create(p.x/rad,p.y/rad,p.z/rad);
      break;
   default:
      break;    
  } 

  return norm;*/
  return GPThreeVector_unit(p);
}



///////////////////////////////////////////////////////////////////////////
//
// Calculate distance to box from an outside point
// - return kInfinity if no intersection.
//
// ALGORITHM:
//
// Check that if point lies outside x/y/z extent of box, travel is towards
// the box (ie. there is a possibility of an intersection)
//
// Calculate pairs of minimum and maximum distances for x/y/z travel for
// intersection with the box's x/y/z extent.
// If there is a valid intersection, it is given by the maximum min distance
// (ie. distance to satisfy x/y/z intersections) *if* <= minimum max distance
// (ie. distance after which 1+ of x/y/z intersections not satisfied)
//
// NOTE:
//
// `Inside' safe - meaningful answers given if point is inside the exact
// shape.

FQUALIFIER
G4double GPOrb_DistanceToIn2( GEOMETRYLOC const GPOrb *This, 
			      GPThreeVector p,
			      GPThreeVector v)
{
  G4double snxt = kInfinity ;      // snxt = default return value

  G4double rad2, pDotV3d;
  G4double c, d2, s = kInfinity ;

  //const G4double dRmax = 100.*This->fRmax;

  // General Precalcs

  rad2    = GPThreeVector_mag2(p);// p.x*p.x + p.y*p.y + p.z*p.z ;
  pDotV3d = GPThreeVector_dot(p,v); //p.x*v.x + p.y*v.y + p.z*v.z ;

  // Radial Precalcs

  //tolORMax2 = (This->fRmax+This->fRmaxTolerance*0.5)*(This->fRmax+This->fRmaxTolerance*0.5) ;
  //tolIRMax2 = (This->fRmax-This->fRmaxTolerance*0.5)*(This->fRmax-This->fRmaxTolerance*0.5) ;

  // Outer spherical shell intersection
  // - Only if outside tolerant fRmax
  // - Check for if inside and outer GPOrb heading through solid (-> 0)
  // - No intersect -> no intersection with GPOrb
  //
  // Shell eqn: x^2+y^2+z^2 = RSPH^2
  //
  // => (px+svx)^2+(py+svy)^2+(pz+svz)^2=R^2
  //
  // => (px^2+py^2+pz^2) +2s(pxvx+pyvy+pzvz)+s^2(vx^2+vy^2+vz^2)=R^2
  // =>      rad2        +2s(pDotV3d)       +s^2                =R^2
  //
  // => s=-pDotV3d+-std::sqrt(pDotV3d^2-(rad2-R^2))


  G4double rad = sqrt(rad2);
  c = (rad - This->fRmax)*(rad + This->fRmax);

  if ( c > This->fRmaxTolerance*This->fRmax )
  {
    // If outside tolerant boundary of outer GPOrb
    // [ should be std::sqrt(rad2) - fRmax > fRmaxTolerance*0.5 ]

    d2 = pDotV3d*pDotV3d - c ;

    if ( d2 >= 0 )
    {
      s = -pDotV3d - sqrt(d2) ;
      if ( s >= 0 )
      {
		// TODO: check and fix This
        /*if ( s>dRmax ) // Avoid rounding errors due to precision issues seen on
        {              // 64 bits systems. Split long distances and recompute
          G4double fTerm = s-fmod(s,dRmax);
          s = fTerm + DistanceToIn(p+fTerm*v,v);
        } */
        return snxt = s;
      }
    }
    else    // No intersection with GPOrb
    {
      return snxt = kInfinity;
    }
  }
  else
  {
    if ( c > -This->fRmaxTolerance*This->fRmax )  // on surface  
    {
      d2 = pDotV3d*pDotV3d - c ;             
      if ( (d2 < This->fRmaxTolerance*This->fRmax) || (pDotV3d >= 0) )
      {
        return snxt = kInfinity;
      }
      else
      {
        return snxt = 0.;
      }
    }
  }
  return snxt;
}

//////////////////////////////////////////////////////////////////////////
// 
// Appoximate distance to box.
// Returns largest perpendicular distance to the closest x/y/z sides of
// the box, which is the most fast estimation of the shortest distance to box
// - If inside return 0

FQUALIFIER
G4double GPOrb_DistanceToIn( GEOMETRYLOC const GPOrb *This, 
			     GPThreeVector p)
{
  G4double safe = 0.0,
           rad  = GPThreeVector_mag(p); //sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
  safe = rad - This->fRmax;
  if( safe < 0 ) { safe = 0.; }
  return safe;
}

/////////////////////////////////////////////////////////////////////////
//
// Calcluate distance to surface of box from inside
// by calculating distances to box's x/y/z planes.
// Smallest distance is exact distance to exiting.
// - Eliminate one side of each pair by considering direction of v
// - when leaving a surface & v.close, return 0

FQUALIFIER
G4double GPOrb_DistanceToOut2( GEOMETRYLOC const GPOrb *This, 
			       GPThreeVector p,GPThreeVector v,
                               const G4bool calcNorm,
			       G4bool *validNorm,GPThreeVector *n)
{
 G4double snxt = kInfinity;     // ??? snxt is default return value
  enum {kNull,kRMax} side = kNull;
  
  G4double rad2,pDotV3d; 
  GPThreeVector ipoint;  //G4double xi,yi,zi; // Intersection point
 
  G4double c,d2;
                 
  rad2    = GPThreeVector_mag2(p); //p.x*p.x + p.y*p.y + p.z*p.z;
  pDotV3d = GPThreeVector_dot(p,v); //p.x*v.x + p.y*v.y + p.z*v.z;
    
  // Radial Intersection from GPOrb::DistanceToIn
  //
  // Outer spherical shell intersection
  // - Only if outside tolerant fRmax
  // - Check for if inside and outer GPOrb heading through solid (-> 0)
  // - No intersect -> no intersection with GPOrb
  //
  // Shell eqn: x^2+y^2+z^2=RSPH^2
  //
  // => (px+svx)^2+(py+svy)^2+(pz+svz)^2=R^2
  //
  // => (px^2+py^2+pz^2) +2s(pxvx+pyvy+pzvz)+s^2(vx^2+vy^2+vz^2)=R^2
  // =>      rad2        +2s(pDotV3d)       +s^2                =R^2
  //
  // => s=-pDotV3d+-std::sqrt(pDotV3d^2-(rad2-R^2))
  
  const G4double  Rmax_plus = This->fRmax + This->fRmaxTolerance*0.5;
  G4double rad = sqrt(rad2);

  if ( rad <= Rmax_plus )
  {
    c = (rad - This->fRmax)*(rad + This->fRmax);

    if ( c < This->fRmaxTolerance*This->fRmax ) 
    {
      // Within tolerant Outer radius 
      // 
      // The test is
      //     rad  - fRmax < 0.5*fRmaxTolerance
      // =>  rad  < fRmax + 0.5*kRadTol
      // =>  rad2 < (fRmax + 0.5*kRadTol)^2
      // =>  rad2 < fRmax^2 + 2.*0.5*fRmax*kRadTol + 0.25*kRadTol*kRadTol
      // =>  rad2 - fRmax^2    <~    fRmax*kRadTol 

      d2 = pDotV3d*pDotV3d - c;

      if( ( c > -This->fRmaxTolerance*This->fRmax) &&         // on tolerant surface
          ( ( pDotV3d >= 0 )   || ( d2 < 0 )) )   // leaving outside from Rmax 
                                                  // not re-entering
      {
        if(calcNorm)
        {
          *validNorm = true ;
          *n         = GPThreeVector_create(p.x/This->fRmax,p.y/This->fRmax,p.z/This->fRmax) ;
        }
        return snxt = 0;
      }
      else 
      {
        snxt = -pDotV3d + sqrt(d2);    // second root since inside Rmax
        side = kRMax ; 
      }
    }
  }
  else // p is outside ???
  {
  }
  if (calcNorm)    // Output switch operator
  {
    switch( side )
    {
      case kRMax:
		ipoint = GPThreeVector_saxpy(snxt,v,p);
		*n=GPThreeVector_mult(ipoint,1.0/This->fRmax);
        //*n=GPThreeVector_create(xi/This->fRmax,yi/This->fRmax,zi/This->fRmax);
        *validNorm=true;
        break;
      default:
        break;
    }
  }
  return snxt;
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - If outside return 0

FQUALIFIER
G4double GPOrb_DistanceToOut( GEOMETRYLOC const GPOrb *This, 
			      GPThreeVector p )
{
   G4double safe=0.0,rad = GPThreeVector_mag(p); //sqrt(p.x*p.x+p.y*p.y+p.z*p.z);

  safe = This->fRmax - rad;
  if ( safe < 0. ) safe = 0.;
  return safe;
}

