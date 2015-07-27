// @(#)root/geom:$Id: TGeoBBox_v.cxx 27731 2009-03-09 17:40:56Z brun $// Author: Andrei Gheata   24/10/01

// Contains() and DistFromOutside/Out() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//--------------------------------------------------------------------------
// TGeoBBox_v - box class. All shape primitives inherit from this, their
//   constructor filling automatically the parameters of the box that bounds
//   the given shape. Defined by 6 parameters :
//      fDX, fDY, fDZ - half lengths on X, Y and Z axis
//      fOrigin[3]    - position of box origin
//
//--------------------------------------------------------------------------
//
//
//--- Building boxes
//  ==================
//  Normally a box has to be build only with 3 parameters : dx, dy, dz
// representing the half lengths on X, Y and Z axis. In this case, the origin
// of the box will match the one of its reference frame. The translation of the
// origin is used only by the constructors of all other shapes in order to
// define their own bounding boxes. Users should be aware that building a
// translated box that will represent a physical shape by itself will affect any
// further positioning of other shapes inside. Therefore in order to build a
// positioned box one should follow the recipe described in class TGeoNode.
//
// Creation of boxes
// 1.   TGeoBBox_v *box = new TGeoBBox_v("BOX", 20, 30, 40);
// Begin_Html
/*
<img src="gif/t_box.gif">
*/
// End_Html
//
// 2. A volume having a box shape can be built in one step:
//      TGeoVolume *vbox = gGeoManager->MakeBox("vbox", ptrMed, 20,30,40);
//
// Divisions of boxes.
//
//   Volumes having box shape can be divided with equal-length slices on
// X, Y or Z axis. The following options are supported:
// a) Dividing the full range of one axis in N slices
//      TGeoVolume *divx = vbox->Divide("SLICEX", 1, N);
//   - here 1 stands for the division axis (1-X, 2-Y, 3-Z)
// Begin_Html
/*
<img src="gif/t_boxdivX.gif">
*/
// End_Html
//
// b) Dividing in a limited range - general case.
//      TGeoVolume *divy = vbox->Divide("SLICEY",2,N,start,step);
//   - start = starting offset within (-fDY, fDY)
//   - step  = slicing step
//
// Begin_Html
/*
<img src="gif/t_boxdivstepZ.gif">
*/
// End_Html
//
// Both cases are supported by all shapes.
//   See also class TGeoShape for utility methods provided by any particular
// shape.
//_____________________________________________________________________________

#include "Riostream.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TGeoBBox_v.h"
#include "TVirtualPad.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"
#include "TRandom.h"

ClassImp(TGeoBBox_v)

    //_____________________________________________________________________________
    TGeoBBox_v::TGeoBBox_v() {
  // Default constructor
  SetShapeBit(TGeoShape::kGeoBox);
  fDX = fDY = fDZ = 0;
  fOrigin[0] = fOrigin[1] = fOrigin[2] = 0.0;
}

//_____________________________________________________________________________
TGeoBBox_v::TGeoBBox_v(double dx, double dy, double dz, double *origin) : TGeoBBox(dx, dy, dz, origin) {}

//_____________________________________________________________________________
TGeoBBox_v::TGeoBBox_v(const char *name, double dx, double dy, double dz, double *origin)
    : TGeoBBox(name, dx, dy, dz, origin) {}

//_____________________________________________________________________________
TGeoBBox_v::TGeoBBox_v(double *param) : TGeoBBox(param) {}

//_____________________________________________________________________________
TGeoBBox_v::~TGeoBBox_v() {
  // Destructor
}

//_____________________________________________________________________________
Bool_t TGeoBBox_v::AreOverlapping(const TGeoBBox_v *box1, const TGeoMatrix *mat1, const TGeoBBox_v *box2,
                                  const TGeoMatrix *mat2) {
  // Check if 2 positioned boxes overlap.
  double master[3];
  double local[3];
  double ldir1[3], ldir2[3];
  const double *o1 = box1->GetOrigin();
  const double *o2 = box2->GetOrigin();
  // Convert center of first box to the local frame of second
  mat1->LocalToMaster(o1, master);
  mat2->MasterToLocal(master, local);
  if (TGeoBBox_v::Contains(local, box2->GetDX(), box2->GetDY(), box2->GetDZ(), o2))
    return kTRUE;
  double distsq = (local[0] - o2[0]) * (local[0] - o2[0]) + (local[1] - o2[1]) * (local[1] - o2[1]) +
                    (local[2] - o2[2]) * (local[2] - o2[2]);
  // Compute distance between box centers and compare with max value
  double rmaxsq = (box1->GetDX() + box2->GetDX()) * (box1->GetDX() + box2->GetDX()) +
                    (box1->GetDY() + box2->GetDY()) * (box1->GetDY() + box2->GetDY()) +
                    (box1->GetDZ() + box2->GetDZ()) * (box1->GetDZ() + box2->GetDZ());
  if (distsq > rmaxsq + TGeoShape::Tolerance())
    return kFALSE;
  // We are still not sure: shoot a ray from the center of "1" towards the
  // center of 2.
  double dir[3];
  mat1->LocalToMaster(o1, ldir1);
  mat2->LocalToMaster(o2, ldir2);
  distsq = sqrt(1 / distsq);
  dir[0] = (ldir2[0] - ldir1[0]) * distsq;
  dir[1] = (ldir2[1] - ldir1[1]) * distsq;
  dir[2] = (ldir2[2] - ldir1[2]) * distsq;
  mat1->MasterToLocalVect(dir, ldir1);
  mat2->MasterToLocalVect(dir, ldir2);
  // Distance to exit from o1
  double dist1 = TGeoBBox_v::DistFromInside(o1, ldir1, box1->GetDX(), box1->GetDY(), box1->GetDZ(), o1);
  // Distance to enter from o2
  double dist2 = TGeoBBox_v::DistFromOutside(local, ldir2, box2->GetDX(), box2->GetDY(), box2->GetDZ(), o2);
  if (dist1 > dist2)
    return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________
void TGeoBBox_v::ComputeNormal(double *point, double *dir, double *norm) {
  // Computes normal to closest surface from POINT.
  memset(norm, 0, 3 * sizeof(double));
  double saf[3];
  int i;
  saf[0] = fabs(fabs(point[0] - fOrigin[0]) - fDX);
  saf[1] = fabs(fabs(point[1] - fOrigin[1]) - fDY);
  saf[2] = fabs(fabs(point[2] - fOrigin[2]) - fDZ);
  i = TMath::LocMin(3, saf);
  norm[i] = (dir[i] > 0) ? 1 : (-1);
}

//_____________________________________________________________________________
Bool_t TGeoBBox_v::CouldBeCrossed(double *point, double *dir) const {
  // Decides fast if the bounding box could be crossed by a vector.
  double mind = fDX;
  if (fDY < mind)
    mind = fDY;
  if (fDZ < mind)
    mind = fDZ;
  double dx = fOrigin[0] - point[0];
  double dy = fOrigin[1] - point[1];
  double dz = fOrigin[2] - point[2];
  double do2 = dx * dx + dy * dy + dz * dz;
  if (do2 <= (mind * mind))
    return kTRUE;
  double rmax2 = fDX * fDX + fDY * fDY + fDZ * fDZ;
  if (do2 <= rmax2)
    return kTRUE;
  // inside bounding sphere
  double doct = dx * dir[0] + dy * dir[1] + dz * dir[2];
  // leaving ray
  if (doct <= 0)
    return kFALSE;
  double dirnorm = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];
  if ((doct * doct) >= (do2 - rmax2) * dirnorm)
    return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________
int TGeoBBox_v::DistancetoPrimitive(int px, int py) {
  // Compute closest distance from point px,py to each corner.
  const int numPoints = 8;
  return ShapeDistancetoPrimitive(numPoints, px, py);
}

//_____________________________________________________________________________
Bool_t TGeoBBox_v::Contains(double *point) const {
  // Test if point is inside this shape.
  if (fabs(point[2] - fOrigin[2]) > fDZ)
    return kFALSE;
  if (fabs(point[0] - fOrigin[0]) > fDX)
    return kFALSE;
  if (fabs(point[1] - fOrigin[1]) > fDY)
    return kFALSE;
  return kTRUE;
}

#define vector(elcount, type) __attribute__((vector_size((elcount) * sizeof(type)))) type
//_____________________________________________________________________________
void TGeoBBox_v::Contains_v(const double *point, Bool_t *isin, const int np) const {
  // Test if point is inside this shape.
  vector(32, double)*vx, *vy;
  double xx, yy, zz;
  for (int i = 0; i < np; ++i) {
    xx = point[3 * i] - fOrigin[0];
    yy = point[3 * i + 1] - fOrigin[1];
    zz = point[3 * i + 2] - fOrigin[2];
    isin[i] = (fabs(xx) < fDX) && (fabs(yy) < fDY) && (fabs(zz) < fDZ);
  }
}

//_____________________________________________________________________________
Bool_t TGeoBBox_v::Contains(const double *point, double dx, double dy, double dz, const double *origin) {
  // Test if point is inside this shape.
  if (fabs(point[2] - origin[2]) > dz)
    return kFALSE;
  if (fabs(point[0] - origin[0]) > dx)
    return kFALSE;
  if (fabs(point[1] - origin[1]) > dy)
    return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
double TGeoBBox_v::DistFromInside(double *point, double *dir, int iact, double step, double *safe) const {
  // Compute distance from inside point to surface of the box.
  // Boundary safe algorithm.
  double s, smin, saf[6];
  double newpt[3];
  int i;
  for (i = 0; i < 3; i++)
    newpt[i] = point[i] - fOrigin[i];
  saf[0] = fDX + newpt[0];
  saf[1] = fDX - newpt[0];
  saf[2] = fDY + newpt[1];
  saf[3] = fDY - newpt[1];
  saf[4] = fDZ + newpt[2];
  saf[5] = fDZ - newpt[2];
  if (iact < 3 && safe) {
    smin = saf[0];
    // compute safe distance
    for (i = 1; i < 6; i++)
      if (saf[i] < smin)
        smin = saf[i];
    *safe = smin;
    if (smin < 0)
      *safe = 0.0;
    if (iact == 0)
      return TGeoShape::Big();
    if (iact == 1 && step < *safe)
      return TGeoShape::Big();
  }
  // compute distance to surface
  smin = TGeoShape::Big();
  for (i = 0; i < 3; i++) {
    if (dir[i] != 0) {
      s = (dir[i] > 0) ? (saf[(i << 1) + 1] / dir[i]) : (-saf[i << 1] / dir[i]);
      if (s < 0)
        return 0.0;
      if (s < smin)
        smin = s;
    }
  }
  return smin;
}

//_____________________________________________________________________________
double TGeoBBox_v::DistFromInside(const double *point, const double *dir, double dx, double dy, double dz,
                                    const double *origin, double /*stepmax*/) {
  // Compute distance from inside point to surface of the box.
  // Boundary safe algorithm.
  double s, smin, saf[6];
  double newpt[3];
  int i;
  for (i = 0; i < 3; i++)
    newpt[i] = point[i] - origin[i];
  saf[0] = dx + newpt[0];
  saf[1] = dx - newpt[0];
  saf[2] = dy + newpt[1];
  saf[3] = dy - newpt[1];
  saf[4] = dz + newpt[2];
  saf[5] = dz - newpt[2];
  // compute distance to surface
  smin = TGeoShape::Big();
  for (i = 0; i < 3; i++) {
    if (dir[i] != 0) {
      s = (dir[i] > 0) ? (saf[(i << 1) + 1] / dir[i]) : (-saf[i << 1] / dir[i]);
      if (s < 0)
        return 0.0;
      if (s < smin)
        smin = s;
    }
  }
  return smin;
}

//_____________________________________________________________________________
double TGeoBBox_v::DistFromOutside(double *point, double *dir, int iact, double step, double *safe) const {
  // Compute distance from outside point to surface of the box.
  // Boundary safe algorithm.
  Bool_t in = kTRUE;
  double saf[3];
  double par[3];
  double newpt[3];
  int i, j;
  for (i = 0; i < 3; i++)
    newpt[i] = point[i] - fOrigin[i];
  par[0] = fDX;
  par[1] = fDY;
  par[2] = fDZ;
  for (i = 0; i < 3; i++) {
    saf[i] = fabs(newpt[i]) - par[i];
    if (saf[i] >= step)
      return TGeoShape::Big();
    if (in && saf[i] > 0)
      in = kFALSE;
  }
  if (iact < 3 && safe) {
    // compute safe distance
    if (in) {
      *safe = 0.0;
    } else {
      *safe = saf[0];
      if (saf[1] > *safe)
        *safe = saf[1];
      if (saf[2] > *safe)
        *safe = saf[2];
    }
    if (iact == 0)
      return TGeoShape::Big();
    if (iact == 1 && step < *safe)
      return TGeoShape::Big();
  }
  // compute distance from point to box
  double coord, snxt = TGeoShape::Big();
  int ibreak = 0;
  // protection in case point is actually inside box
  if (in) {
    j = 0;
    double ss = saf[0];
    if (saf[1] > ss) {
      ss = saf[1];
      j = 1;
    }
    if (saf[2] > ss)
      j = 2;
    if (newpt[j] * dir[j] > 0)
      return TGeoShape::Big(); // in fact exiting
    return 0.0;
  }
  for (i = 0; i < 3; i++) {
    if (saf[i] < 0)
      continue;
    if (newpt[i] * dir[i] >= 0)
      continue;
    snxt = saf[i] / fabs(dir[i]);
    ibreak = 0;
    for (j = 0; j < 3; j++) {
      if (j == i)
        continue;
      coord = newpt[j] + snxt * dir[j];
      if (fabs(coord) > par[j]) {
        ibreak = 1;
        break;
      }
    }
    if (!ibreak)
      return snxt;
  }
  return TGeoShape::Big();
}

//_____________________________________________________________________________
double TGeoBBox_v::DistFromOutside(const double *point, const double *dir, double dx, double dy, double dz,
                                     const double *origin, double stepmax) {
  // Compute distance from outside point to surface of the box.
  // Boundary safe algorithm.
  Bool_t in = kTRUE;
  double saf[3];
  double par[3];
  double newpt[3];
  int i, j;
  for (i = 0; i < 3; i++)
    newpt[i] = point[i] - origin[i];
  par[0] = dx;
  par[1] = dy;
  par[2] = dz;
  for (i = 0; i < 3; i++) {
    saf[i] = fabs(newpt[i]) - par[i];
    if (saf[i] >= stepmax)
      return TGeoShape::Big();
    if (in && saf[i] > 0)
      in = kFALSE;
  }
  // In case point is inside return ZERO
  if (in)
    return 0.0;
  double coord, snxt = TGeoShape::Big();
  int ibreak = 0;
  for (i = 0; i < 3; i++) {
    if (saf[i] < 0)
      continue;
    if (newpt[i] * dir[i] >= 0)
      continue;
    snxt = saf[i] / fabs(dir[i]);
    ibreak = 0;
    for (j = 0; j < 3; j++) {
      if (j == i)
        continue;
      coord = newpt[j] + snxt * dir[j];
      if (fabs(coord) > par[j]) {
        ibreak = 1;
        break;
      }
    }
    if (!ibreak)
      return snxt;
  }
  return TGeoShape::Big();
}

//_____________________________________________________________________________
Bool_t TGeoBBox_v::GetPointsOnFacet(int index, int npoints, double *array) const {
  // Fills array with n random points located on the surface of indexed facet.
  // The output array must be provided with a length of minimum 3*npoints. Returns
  // true if operation succeeded.
  // Possible index values:
  //    0 - all facets togeather
  //    1 to 6 - facet index from bottom to top Z
  if (index < 0 || index > 6)
    return kFALSE;
  double surf[6];
  double area = 0.;
  if (index == 0) {
    for (int isurf = 0; isurf < 6; isurf++) {
      surf[isurf] = TGeoBBox_v::GetFacetArea(isurf + 1);
      if (isurf > 0)
        surf[isurf] += surf[isurf - 1];
    }
    area = surf[5];
  }

  for (int i = 0; i < npoints; i++) {
    // Generate randomly a surface index if needed.
    double *point = &array[3 * i];
    int surfindex = index;
    if (surfindex == 0) {
      double val = area * gRandom->Rndm();
      surfindex = 2 + TMath::BinarySearch(6, surf, val);
      if (surfindex > 6)
        surfindex = 6;
    }
    switch (surfindex) {
    case 1:
      point[0] = -fDX + 2 * fDX * gRandom->Rndm();
      point[1] = -fDY + 2 * fDY * gRandom->Rndm();
      point[2] = -fDZ;
      break;
    case 2:
      point[0] = -fDX + 2 * fDX * gRandom->Rndm();
      point[1] = -fDY;
      point[2] = -fDZ + 2 * fDZ * gRandom->Rndm();
      break;
    case 3:
      point[0] = -fDX;
      point[1] = -fDY + 2 * fDY * gRandom->Rndm();
      point[2] = -fDZ + 2 * fDZ * gRandom->Rndm();
      break;
    case 4:
      point[0] = -fDX + 2 * fDX * gRandom->Rndm();
      point[1] = fDY;
      point[2] = -fDZ + 2 * fDZ * gRandom->Rndm();
      break;
    case 5:
      point[0] = fDX;
      point[1] = -fDY + 2 * fDY * gRandom->Rndm();
      point[2] = -fDZ + 2 * fDZ * gRandom->Rndm();
      break;
    case 6:
      point[0] = -fDX + 2 * fDX * gRandom->Rndm();
      point[1] = -fDY + 2 * fDY * gRandom->Rndm();
      point[2] = fDZ;
      break;
    }
  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t TGeoBBox_v::GetPointsOnSegments(int npoints, double *array) const {
  // Fills array with n random points located on the line segments of the shape mesh.
  // The output array must be provided with a length of minimum 3*npoints. Returns
  // true if operation is implemented.
  if (npoints < GetNmeshVertices()) {
    Error("GetPointsOnSegments", "You should require at least %d points", GetNmeshVertices());
    return kFALSE;
  }
  TBuffer3D &buff = (TBuffer3D &)GetBuffer3D(TBuffer3D::kRawSizes | TBuffer3D::kRaw, kTRUE);
  int npnts = buff.NbPnts();
  int nsegs = buff.NbSegs();
  // Copy buffered points  in the array
  memcpy(array, buff.fPnts, 3 * npnts * sizeof(double));
  int ipoints = npoints - npnts;
  int icrt = 3 * npnts;
  int nperseg = (int)(double(ipoints) / nsegs);
  double *p0, *p1;
  double x, y, z, dx, dy, dz;
  for (int i = 0; i < nsegs; i++) {
    p0 = &array[3 * buff.fSegs[3 * i + 1]];
    p1 = &array[3 * buff.fSegs[3 * i + 2]];
    if (i == (nsegs - 1))
      nperseg = ipoints;
    dx = (p1[0] - p0[0]) / (nperseg + 1);
    dy = (p1[1] - p0[1]) / (nperseg + 1);
    dz = (p1[2] - p0[2]) / (nperseg + 1);
    for (int j = 0; j < nperseg; j++) {
      x = p0[0] + (j + 1) * dx;
      y = p0[1] + (j + 1) * dy;
      z = p0[2] + (j + 1) * dz;
      array[icrt++] = x;
      array[icrt++] = y;
      array[icrt++] = z;
      ipoints--;
    }
  }
  return kTRUE;
}

//_____________________________________________________________________________
double TGeoBBox_v::Safety(double *point, Bool_t in) const {
  // Computes the closest distance from given point to this shape.

  double safe, safy, safz;
  if (in) {
    safe = fDX - fabs(point[0] - fOrigin[0]);
    safy = fDY - fabs(point[1] - fOrigin[1]);
    safz = fDZ - fabs(point[2] - fOrigin[2]);
    if (safy < safe)
      safe = safy;
    if (safz < safe)
      safe = safz;
  } else {
    safe = -fDX + fabs(point[0] - fOrigin[0]);
    safy = -fDY + fabs(point[1] - fOrigin[1]);
    safz = -fDZ + fabs(point[2] - fOrigin[2]);
    if (safy > safe)
      safe = safy;
    if (safz > safe)
      safe = safz;
  }
  return safe;
}
