#include "GPBox.h"
#include "GPUtils.h"
#include "GPGeomdefs.h"

FQUALIFIER
void GPBox_Constructor(GPBox *This, G4double pX, G4double pY, G4double pZ)
{
  This->fSolid.fType = kBox;
  
  This->fDx = pX;
  This->fDy = pY;
  This->fDz = pZ;
}

FQUALIFIER
void GPBox_SetXHalfLength(GPBox *This,
                          G4double dx)
{
  if(dx > 2*kCarTolerance)  // limit to thickness of surfaces
  {
    This->fDx = dx;
  }
  else
  {
    ;
  }
} 

FQUALIFIER
void GPBox_SetYHalfLength(GPBox *This,
                          G4double dy) 
{
  if(dy > 2*kCarTolerance)  // limit to thickness of surfaces
  {
    This->fDy = dy;
  }
  else
  {
    ;
  }
} 

FQUALIFIER
void GPBox_SetZHalfLength(GPBox *This,
                          G4double dz) 
{
  if(dz > 2*kCarTolerance)  // limit to thickness of surfaces
  {
    This->fDz = dz;
  }
  else
  {
    ;
  }
} 


FQUALIFIER
G4bool GPBox_CalculateExtent(GPBox *This,
                             const EAxis pAxis,
                             GPVoxelLimits pVoxelLimit,
                             GPAffineTransform pTransform,
                             G4double *pMin, G4double *pMax)
{
  if (! GPAffineTransform_IsRotated(&pTransform)) {
    // Special case handling for unrotated boxes
    // Compute x/y/z mins and maxs respecting limits, with early returns
    // if outside limits. Then switch() on pAxis
    
    G4double xoffset,xMin,xMax;
    G4double yoffset,yMin,yMax;
    G4double zoffset,zMin,zMax;
    
    xoffset = GPAffineTransform_NetTranslation(&pTransform).x ;
    xMin    = xoffset - This->fDx ;
    xMax    = xoffset + This->fDx ;
    
    if (GPVoxelLimits_IsXLimited(&pVoxelLimit)) {
      if ((xMin > GPVoxelLimits_GetMaxXExtent(&pVoxelLimit)+kCarTolerance) || 
          (xMax < GPVoxelLimits_GetMinXExtent(&pVoxelLimit)-kCarTolerance)) { return false ; }
      else {
        xMin = GPfmax(xMin, GPVoxelLimits_GetMinXExtent(&pVoxelLimit));
        xMax = GPfmin(xMax, GPVoxelLimits_GetMaxXExtent(&pVoxelLimit));
      }
    }
    yoffset = GPAffineTransform_NetTranslation(&pTransform).y ;
    yMin    = yoffset - This->fDy ;
    yMax    = yoffset + This->fDy ;
    
    if (GPVoxelLimits_IsYLimited(&pVoxelLimit)) {
      if ((yMin > GPVoxelLimits_GetMaxYExtent(&pVoxelLimit)+kCarTolerance) ||
          (yMax < GPVoxelLimits_GetMinYExtent(&pVoxelLimit)-kCarTolerance)) { return false ; }
      else {
        yMin = GPfmax(yMin, GPVoxelLimits_GetMinYExtent(&pVoxelLimit));
        yMax = GPfmin(yMax, GPVoxelLimits_GetMaxYExtent(&pVoxelLimit));
      }
    }
    zoffset = GPAffineTransform_NetTranslation(&pTransform).z ;
    zMin    = zoffset - This->fDz ;
    zMax    = zoffset + This->fDz ;

    if (GPVoxelLimits_IsZLimited(&pVoxelLimit)) {
      if ((zMin > GPVoxelLimits_GetMaxZExtent(&pVoxelLimit)+kCarTolerance) ||
          (zMax < GPVoxelLimits_GetMinZExtent(&pVoxelLimit)-kCarTolerance)) { return false ; }
      else {
        zMin = GPfmax(zMin, GPVoxelLimits_GetMinZExtent(&pVoxelLimit));
        zMax = GPfmin(zMax, GPVoxelLimits_GetMaxZExtent(&pVoxelLimit));
      }
    }
    switch (pAxis) {
    case kXAxis:
      *pMin = xMin ;
      *pMax = xMax ;
      break ;
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
    *pMin -= kCarTolerance ;
    *pMax += kCarTolerance ;
    
    return true;
  }
  else { // General rotated case - create and clip mesh to boundaries
    G4bool existsAfterClip = false ;
    GPThreeVectorList vertices = GPBox_CreateRotatedVertices(This,pTransform);
    
    *pMin = +kInfinity ;
    *pMax = -kInfinity ;
    
    // Calculate rotated vertex coordinates
    
    GPVSolid_ClipCrossSection(&vertices,0,&pVoxelLimit,pAxis,pMin,pMax) ;
    GPVSolid_ClipCrossSection(&vertices,4,&pVoxelLimit,pAxis,pMin,pMax) ;
    GPVSolid_ClipBetweenSections(&vertices,0,&pVoxelLimit,pAxis,pMin,pMax) ;
    
    if ( GPVoxelLimits_IsLimited2(&pVoxelLimit,pAxis) == false) {  
      if ( (*pMin != kInfinity) || (*pMax != -kInfinity) ) {
        existsAfterClip = true ;
        
        // Add 2*tolerance to avoid precision troubles
        
        *pMin -= kCarTolerance;
        *pMax += kCarTolerance;
      }
    }      
    else {
      GPThreeVector clipCentre = GPThreeVector_create(
      ( GPVoxelLimits_GetMinXExtent(&pVoxelLimit)+GPVoxelLimits_GetMaxXExtent(&pVoxelLimit))*0.5,
      ( GPVoxelLimits_GetMinYExtent(&pVoxelLimit)+GPVoxelLimits_GetMaxYExtent(&pVoxelLimit))*0.5,
      ( GPVoxelLimits_GetMinZExtent(&pVoxelLimit)+GPVoxelLimits_GetMaxZExtent(&pVoxelLimit))*0.5);
      
      GPAffineTransform invT = GPAffineTransform_Inverse( &pTransform );
      GPThreeVector aT = GPAffineTransform_TransformPoint(&invT, clipCentre);
      
      if ( (*pMin != kInfinity) || (*pMax != -kInfinity)) {
        existsAfterClip = true ;
        
        
        // Check to see if endpoints are in the solid        
        //        clipCentre(pAxis) = GPVoxelLimits_GetMinExtent(&pVoxelLimit,pAxis);
        
        G4double minExt = GPVoxelLimits_GetMinExtent(&pVoxelLimit,pAxis);
        switch( pAxis )
          {
          case kXAxis: clipCentre.x = minExt; break;
          case kYAxis: clipCentre.y = minExt; break;
          case kZAxis: clipCentre.z = minExt; break;
          default:
            break;
          }
        
        aT = GPAffineTransform_TransformPoint(&invT, clipCentre);
        
        if ( GPBox_Inside(This,aT) != kOutside) {
          *pMin = GPVoxelLimits_GetMinExtent(&pVoxelLimit,pAxis);
        }
        else {
          *pMin -= kCarTolerance;
        }
        //        clipCentre(pAxis) = pVoxelLimit.GetMaxExtent(pAxis);
        G4double maxExt = GPVoxelLimits_GetMaxExtent(&pVoxelLimit,pAxis);
        switch( pAxis )
          {
          case kXAxis: clipCentre.x = maxExt; break;
          case kYAxis: clipCentre.y = maxExt; break;
          case kZAxis: clipCentre.z = maxExt; break;
          default:
            break;
          }
        
        aT = GPAffineTransform_TransformPoint(&invT, clipCentre);
        //        if (Inside(pTransform.Inverse().TransformPoint(clipCentre)) != kOutside)
        if ( GPBox_Inside(This,aT) != kOutside) {
          *pMax = GPVoxelLimits_GetMaxExtent(&pVoxelLimit,pAxis);
        }
        else {
          *pMax += kCarTolerance;
        }
      }
      
      // Check for case where completely enveloping clipping volume
      // If point inside then we are confident that the solid completely
      // envelopes the clipping volume. Hence set min/max extents according
      // to clipping volume extents along the specified axis.
      
      //      else if (Inside(pTransform.Inverse().TransformPoint(clipCentre)) != kOutside) {
      else if ( GPBox_Inside(This,aT) != kOutside) {
        existsAfterClip = true ;
        *pMin            = GPVoxelLimits_GetMinExtent(&pVoxelLimit,pAxis) ;
        *pMax            = GPVoxelLimits_GetMaxExtent(&pVoxelLimit,pAxis) ;
      }
    } 
    //    delete vertices;
    return existsAfterClip;
  } 
} 

FQUALIFIER
EInside GPBox_Inside( GEOMETRYLOC GPBox *This,
		      GPThreeVector p)
{
  //  static const 
  G4double delta=0.5*kCarTolerance;
  EInside in = kOutside ;
  GPThreeVector q = GPThreeVector_create(fabs(p.x), fabs(p.y), fabs(p.z));

  if ( q.x <= (This->fDx - delta) )
  {
    if (q.y <= (This->fDy - delta) )
    {
      if      ( q.z <= (This->fDz - delta) ) { in = kInside ;  }
      else if ( q.z <= (This->fDz + delta) ) { in = kSurface ; }
    }
    else if ( q.y <= (This->fDy + delta) )
    {
      if ( q.z <= (This->fDz + delta) ) { in = kSurface ; }
    }
  }
  else if ( q.x <= (This->fDx + delta) )
  {
    if ( q.y <= (This->fDy + delta) )
    {
      if ( q.z <= (This->fDz + delta) ) { in = kSurface ; }
    }
  }
  return in ;
}

FQUALIFIER
GPThreeVector GPBox_ApproxSurfaceNormal( GEOMETRYLOC GPBox *This,
					 GPThreeVector p )
{
  G4double distx, disty, distz ;
  GPThreeVector norm = GPThreeVector_create(0.,0.,0.);
  
  // Calculate distances as if in 1st octant
  
  distx = fabs(fabs(p.x) - This->fDx) ;
  disty = fabs(fabs(p.y) - This->fDy) ;
  distz = fabs(fabs(p.z) - This->fDz) ;

  if ( distx <= disty )
  {
    if ( distx <= distz )     // Closest to X
      {
	if ( p.x < 0 ) { norm = GPThreeVector_create(-1.0,0,0) ; }
	else           { norm = GPThreeVector_create( 1.0,0,0) ; }
      }
    else                      // Closest to Z
      {
	if ( p.z < 0 ) { norm = GPThreeVector_create(0,0,-1.0) ; }
	else           { norm = GPThreeVector_create(0,0, 1.0) ; }
      }
  }
  else
  {
    if ( disty <= distz )      // Closest to Y
      {
	if ( p.y < 0 ) { norm = GPThreeVector_create(0,-1.0,0) ; }
	else           { norm = GPThreeVector_create(0, 1.0,0) ; }
      }
    else                       // Closest to Z
      {
	if ( p.z < 0 ) { norm = GPThreeVector_create(0,0,-1.0) ; }
	else           { norm = GPThreeVector_create(0,0, 1.0) ; }
      }
  }
  return norm;
}

FQUALIFIER
GPThreeVector GPBox_SurfaceNormal( GEOMETRYLOC GPBox *This,
				   GPThreeVector p)
{
  G4double distx, disty, distz ;
  GPThreeVector norm = GPThreeVector_create(0.,0.,0.);

  // Calculate distances as if in 1st octant

  distx = fabs(fabs(p.x) - This->fDx) ;
  disty = fabs(fabs(p.y) - This->fDy) ;
  distz = fabs(fabs(p.z) - This->fDz) ;

  // New code for particle on surface including edges and corners with specific
  // normals

  //  static const 
  G4double delta    = 0.5*kCarTolerance;
  const GPThreeVector nX  = GPThreeVector_create( 1.0, 0,0  );
  const GPThreeVector nmX = GPThreeVector_create(-1.0, 0,0  );
  const GPThreeVector nY  = GPThreeVector_create( 0, 1.0,0  );
  const GPThreeVector nmY = GPThreeVector_create( 0,-1.0,0  );
  const GPThreeVector nZ  = GPThreeVector_create( 0, 0,  1.0);
  const GPThreeVector nmZ = GPThreeVector_create( 0, 0,- 1.0);

  GPThreeVector normX = GPThreeVector_create(0.,0.,0.);
  GPThreeVector normY = GPThreeVector_create(0.,0.,0.);
  GPThreeVector normZ = GPThreeVector_create(0.,0.,0.);
  GPThreeVector sumnorm = GPThreeVector_create(0.,0.,0.);

  G4int noSurfaces=0; 

  if (distx <= delta)         // on X/mX surface and around
  {
    noSurfaces ++; 
    if ( p.x >= 0. )  { normX= nX ; }       // on +X surface : (1,0,0)
    else                { normX= nmX; }       //                 (-1,0,0)
    sumnorm= normX; 
  }

  if (disty <= delta)    // on one of the +Y or -Y surfaces
  {
    noSurfaces ++; 
    if ( p.y >= 0. )  { normY= nY;  }       // on +Y surface
    else                { normY= nmY; }
    //    sumnorm += normY; 
    sumnorm = GPThreeVector_add(sumnorm,normY); 
  }

  if (distz <= delta)    // on one of the +Z or -Z surfaces
  {
    noSurfaces ++; 
    if ( p.z >= 0. )  { normZ= nZ;  }       // on +Z surface
    else                { normZ= nmZ; }
    //    sumnorm += normZ;
    sumnorm = GPThreeVector_add(sumnorm,normZ);
  }

  //  static const 
  G4double invSqrt2 = 1.0 / sqrt(2.0); 
  //  static const 
  G4double invSqrt3 = 1.0 / sqrt(3.0); 

  if( noSurfaces > 0 )
  { 
    if( noSurfaces == 1 )
    { 
      norm= sumnorm; 
    }
    else
    {
      // norm = sumnorm . unit(); 
      if( noSurfaces == 2 )
      { 
        // 2 surfaces -> on edge 
        //        norm = invSqrt2 * sumnorm; 
        norm = GPThreeVector_mult(sumnorm,invSqrt2); 
      }
      else
      { 
        // 3 surfaces (on corner)
        //        norm = invSqrt3 * sumnorm; 
        norm = GPThreeVector_mult(sumnorm,invSqrt3); 
      }
    }
  }
  else
  {
    norm = GPBox_ApproxSurfaceNormal(This,p);
  }
  
  return norm;
}

FQUALIFIER
G4double GPBox_DistanceToIn2( GEOMETRYLOC GPBox *This,
			      GPThreeVector p,
			      GPThreeVector v)
{
  G4double safx, safy, safz ;
  G4double smin=0.0, sminy, sminz ; // , sminx ;
  G4double smax=kInfinity, smaxy, smaxz ; // , smaxx ;  // they always > 0
  G4double stmp ;
  G4double sOut=kInfinity, sOuty=kInfinity, sOutz=kInfinity ;

  //  static const 
  G4double delta = 0.5*kCarTolerance;

  safx = fabs(p.x) - This->fDx ;     // minimum distance to x surface of shape
  safy = fabs(p.y) - This->fDy ;
  safz = fabs(p.z) - This->fDz ;

  // Will we intersect?
  // If safx/y/z is >-tol/2 the point is outside/on the box's x/y/z extent.
  // If both p.x/y/z and v.x/y/z repectively are both positive/negative,
  // travel is in a direction away from the shape.

  if (    ((p.x*v.x >= 0.0) && (safx > -delta)) 
       || ((p.y*v.y >= 0.0) && (safy > -delta))
       || ((p.z*v.z >= 0.0) && (safz > -delta))   ) 
  {
    return kInfinity ;  // travel away or parallel within tolerance
  }

  // Compute min / max distances for x/y/z travel:
  // X Planes

  if ( v.x )  // != 0
  {
    stmp = 1.0/fabs(v.x) ;

    if (safx >= 0.0)
    {
      smin = safx*stmp ;
      smax = (This->fDx+fabs(p.x))*stmp ;
    }
    else
    {
      if (v.x < 0)  { sOut = (This->fDx + p.x)*stmp ; }
      else            { sOut = (This->fDx - p.x)*stmp ; }
    }
  }

  // Y Planes

  if ( v.y )  // != 0
  {
    stmp = 1.0/fabs(v.y) ;

    if (safy >= 0.0)
    {
      sminy = safy*stmp ;
      smaxy = (This->fDy+fabs(p.y))*stmp ;

      if (sminy > smin) { smin=sminy ; }
      if (smaxy < smax) { smax=smaxy ; }

      if (smin >= (smax-delta))
      {
        return kInfinity ;  // touch XY corner
      }
    }
    else
    {
      if (v.y < 0)  { sOuty = (This->fDy + p.y)*stmp ; }
      else            { sOuty = (This->fDy - p.y)*stmp ; }
      if( sOuty < sOut ) { sOut = sOuty ; }
    }     
  }

  // Z planes

  if ( v.z )  // != 0
  {
    stmp = 1.0/fabs(v.z) ;

    if ( safz >= 0.0 )
    {
      sminz = safz*stmp ;
      smaxz = (This->fDz+fabs(p.z))*stmp ;

      if (sminz > smin) { smin = sminz ; }
      if (smaxz < smax) { smax = smaxz ; }

      if (smin >= (smax-delta))
      { 
        return kInfinity ;    // touch ZX or ZY corners
      }
    }
    else
    {
      if (v.z < 0)  { sOutz = (This->fDz + p.z)*stmp ; }
      else            { sOutz = (This->fDz - p.z)*stmp ; }
      if( sOutz < sOut ) { sOut = sOutz ; }
    }
  }

  if (sOut <= (smin + delta)) // travel over edge
  {
    return kInfinity ;
  }
  if (smin < delta)  { smin = 0.0 ; }

  return smin ;
}

FQUALIFIER
G4double GPBox_DistanceToIn( GEOMETRYLOC GPBox *This, 
			     GPThreeVector p)
{
  G4double safex, safey, safez, safe = 0.0 ;

  safex = fabs(p.x) - This->fDx ;
  safey = fabs(p.y) - This->fDy ;
  safez = fabs(p.z) - This->fDz ;

  if (safex > safe) { safe = safex ; }
  if (safey > safe) { safe = safey ; }
  if (safez > safe) { safe = safez ; }

  return safe ;
}

FQUALIFIER
G4double GPBox_DistanceToOut2( GEOMETRYLOC GPBox *This,
			       GPThreeVector p,
			       GPThreeVector v,
			       const G4bool calcNorm,
			       G4bool *validNorm, GPThreeVector *n)
{
  enum ESide {kUndefined,kPX,kMX,kPY,kMY,kPZ,kMZ};

  ESide side = kUndefined ;
  G4double pdist,stmp,snxt=kInfinity;

  //  static const 
  G4double delta = 0.5*kCarTolerance;

  if (calcNorm) { *validNorm = true ; }  // All normals are valid

  if (v.x > 0)   // X planes
  {
    pdist = This->fDx - p.x ;

    if (pdist > delta)
    {
      snxt = pdist/v.x ;
      side = kPX ;
    }
    else
    {
      if (calcNorm) { *n   = GPThreeVector_create(1,0,0) ; }
      return        snxt = 0 ;
    }
  }
  else if (v.x < 0)
  {
    pdist = This->fDx + p.x ;

    if (pdist > delta)
    {
      snxt = -pdist/v.x ;
      side = kMX ;
    }
    else
    {
      if (calcNorm) { *n   = GPThreeVector_create(-1,0,0) ; }
      return        snxt = 0 ;
    }
  }

  if (v.y > 0)   // Y planes
  {
    pdist = This->fDy-p.y;

    if (pdist > delta)
    {
      stmp = pdist/v.y;

      if (stmp < snxt)
      {
        snxt = stmp;
        side = kPY;
      }
    }
    else
    {
      if (calcNorm) { *n   = GPThreeVector_create(0,1,0) ; }
      return        snxt = 0 ;
    }
  }
  else if (v.y < 0)
  {
    pdist = This->fDy + p.y ;

    if (pdist > delta)
    {
      stmp = -pdist/v.y;

      if ( stmp < snxt )
      {
        snxt = stmp;
        side = kMY;
      }
    }
    else
    {
      if (calcNorm) { *n   = GPThreeVector_create(0,-1,0) ; }
      return        snxt = 0 ;
    }
  }
  if (v.z > 0)        // Z planes
  {
    pdist = This->fDz-p.z;

    if ( pdist > delta )
    {
      stmp = pdist/v.z;

      if ( stmp < snxt )
      {
        snxt = stmp;
        side = kPZ;
      }
    }
    else
    {
      if (calcNorm) { *n   = GPThreeVector_create(0,0,1) ; } 
      return        snxt = 0 ;
    }
  }
  else if (v.z < 0)
  {
    pdist = This->fDz + p.z;

    if ( pdist > delta )
    {
      stmp = -pdist/v.z;

      if ( stmp < snxt )
      {
        snxt = stmp;
        side = kMZ;
      }
    }
    else
    {
      if (calcNorm) { *n   = GPThreeVector_create(0,0,-1) ; }
      return        snxt = 0 ;
    }
  }
  if (calcNorm)
  {      
    switch (side)
    {
      case kPX:
        *n=GPThreeVector_create(1,0,0);
        break;
      case kMX:
        *n=GPThreeVector_create(-1,0,0);
        break;
      case kPY:
        *n=GPThreeVector_create(0,1,0);
        break;
      case kMY:
        *n=GPThreeVector_create(0,-1,0);
        break;
      case kPZ:
        *n=GPThreeVector_create(0,0,1);
        break;
      case kMZ:
        *n=GPThreeVector_create(0,0,-1);
        break;
      default:
        ;
        break;
    }
  }
  return snxt;
}

FQUALIFIER
G4double GPBox_DistanceToOut( GEOMETRYLOC GPBox *This,
			      GPThreeVector p)
{
  G4double safx1,safx2,safy1,safy2,safz1,safz2,safe=0.0;

  safx1 = This->fDx - p.x ;
  safx2 = This->fDx + p.x ;
  safy1 = This->fDy - p.y ;
  safy2 = This->fDy + p.y ;
  safz1 = This->fDz - p.z ;
  safz2 = This->fDz + p.z ;  
  
  // shortest Dist to any boundary now MIN(safx1,safx2,safy1..)

  if (safx2 < safx1) { safe = safx2; }
  else               { safe = safx1; }
  if (safy1 < safe)  { safe = safy1; }
  if (safy2 < safe)  { safe = safy2; }
  if (safz1 < safe)  { safe = safz1; }
  if (safz2 < safe)  { safe = safz2; }

  if (safe < 0) { safe = 0 ; }
  return safe ;  
}

FQUALIFIER
GPThreeVectorList
GPBox_CreateRotatedVertices(GPBox *This,
                            GPAffineTransform pTransform)
{
  GPThreeVectorList vertices;
  GPThreeVectorList_Constructor(&vertices);

  if (&vertices)
  {
    //    vertices->reserve(8);
    GPThreeVector vertex0 = GPThreeVector_create(-This->fDx,-This->fDy,-This->fDz) ;
    GPThreeVector vertex1 = GPThreeVector_create(This->fDx,-This->fDy,-This->fDz) ;
    GPThreeVector vertex2 = GPThreeVector_create(This->fDx,This->fDy,-This->fDz) ;
    GPThreeVector vertex3 = GPThreeVector_create(-This->fDx,This->fDy,-This->fDz) ;
    GPThreeVector vertex4 = GPThreeVector_create(-This->fDx,-This->fDy,This->fDz) ;
    GPThreeVector vertex5 = GPThreeVector_create(This->fDx,-This->fDy,This->fDz) ;
    GPThreeVector vertex6 = GPThreeVector_create(This->fDx,This->fDy,This->fDz) ;
    GPThreeVector vertex7 = GPThreeVector_create(-This->fDx,This->fDy,This->fDz) ;

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
    G4Exception("GPBox_CreateRotatedVertices()",
                "GeomSolids0003", FatalException,
                "Error in allocation of vertices. Out of memory !");
    */
  }
  return vertices;
}

FQUALIFIER
G4double GPBox_GetXHalfLength(GPBox *This)
{
  return This->fDx;
}

FQUALIFIER    
G4double GPBox_GetYHalfLength(GPBox *This)
{
  return This->fDy;
}

FQUALIFIER
G4double GPBox_GetZHalfLength(GPBox *This)
{
  return This->fDz;
}
