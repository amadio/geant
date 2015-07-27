#include <iostream>
#include "TGeoCone_v.h"
#include "base/Global.h"
using vecgeom::kRadToDeg;
#ifdef VEC_EXTENSIONS
#include "Vc/vector.h"
#include <Vc/double_v>
typedef Vc::double_v vd;
#endif

#ifndef VEC_EXTENSIONS
//_____________________________________________________________________________
void TGeoCone_v::Contains_v(const StructOfCoord &pointi, Bool_t *isin, Int_t np) const {
  for (unsigned int i = 0; i < np; i++) {
    Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    isin[i] = TGeoCone::Contains(point);
  }
}
#else
// PUT VC CODE OR THE LIKE HERE
//_____________________________________________________________________________
void TGeoCone_v::Contains_v(const StructOfCoord &pointi, Bool_t *isin, Int_t np) const {
  static vd vfDz(fDz);
  static vd vfRmin1(fRmin1);
  static vd vfRmin2(fRmin2);
  static vd vfRmax1(fRmax1);
  static vd vfRmax2(fRmax2);

  Vc::double_m particleinside;
  std::cout << "fabio: " << Vc::double_v::Size << std::endl;

  int tailsize = np % Vc::double_v::Size;
  for (unsigned int i = 0; i < np - tailsize; i += Vc::double_v::Size) {
    vd x(&pointi.x[i]); // will copy a certain number of x's into vc vector x;
    vd y(&pointi.y[i]);
    vd z(&pointi.z[i]);

    vd me = Vc::abs(z);
    Vc::double_m c1 = (me > vfDz);
    //  if( c1 ) continue;

    vd r2 = x * x + y * y;
    vd rl = 0.5 * (vfRmin2 * (z + vfDz) + vfRmin1 * (vfDz - z)) / vfDz;
    vd rh = 0.5 * (vfRmax2 * (z + vfDz) + vfRmax1 * (vfDz - z)) / vfDz;
    Vc::double_m c2 = (r2 < rl * rl);
    // if( c2 ) continue;

    Vc::double_m c3 = (r2 > rh * rh);
    // if( c3 ) continue;

    // Vc::double_m particleinside;
    particleinside = !c1 && !c2 && !c3;
    for (unsigned int j = 0; j < Vc::double_v::Size; ++j)
      isin[i + j] = particleinside[j];
  }

  // do the tail part for the moment, we just call the old static version
  for (unsigned int i = 0; i < tailsize; ++i) {
    Double_t xx, yy, zz;
    xx = pointi.x[np - tailsize + i];
    yy = pointi.y[np - tailsize + i];
    zz = pointi.z[np - tailsize + i];

    Double_t rr2 = xx * xx + yy * yy;
    Double_t rrl = 0.5 * (fRmin2 * (zz + fDz) + fRmin1 * (fDz - zz)) / fDz;
    Double_t rrh = 0.5 * (fRmax2 * (zz + fDz) + fRmax1 * (fDz - zz)) / fDz;

    isin[i] = ((fabs(zz) <= fDz) & (rr2 >= rrl * rrl) & (rr2 <= rrh * rrh));
  }
}
#endif

//_____________________________________________________________________________

void TGeoCone_v::Safety_v(const StructOfCoord &pointi, Bool_t in, Double_t *safety, Int_t np) const {
  for (unsigned int i = 0; i < np; i++) {
    Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    safety[i] = TGeoCone::Safety(point, in);
  }
}

//_____________________________________________________________________________
void TGeoCone_v::DistFromInside_v(const StructOfCoord &pointi, const StructOfCoord &diri, Int_t /*iact*/,
                                  const Double_t *step, Double_t * /*safe*/, Double_t *distance, Int_t np) const {
  for (unsigned int i = 0; i < np; i++) {
    Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};

    distance[i] = TGeoCone::DistFromInside(point, dir, 3, step[i], 0);
  }
}

//_____________________________________________________________________________

void TGeoCone_v::DistFromOutside_v(const StructOfCoord &pointi, const StructOfCoord &diri, Int_t /*iact*/,
                                   const Double_t *step, Double_t * /*safe*/, Double_t *distance, Int_t np) const {
  for (unsigned int i = 0; i < np; i++) {
    Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};

    distance[i] = TGeoCone::DistFromOutside(point, dir, 3, step[i], 0);
  }
}

#ifndef VEC_EXTENSIONS
//_____________________________________________________________________________
void TGeoConeSeg_v::Contains_v(const StructOfCoord &pointi, Bool_t *isin, Int_t np) const {
  std::cout << "ma almeno qui? " << std::endl;

  for (unsigned int i = 0; i < np; i++) {

    Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    isin[i] = TGeoConeSeg::Contains(point);
  }
}
#else

//_____________________________________________________________________________
void TGeoConeSeg_v::Contains_v4(Vc::double_v const &x, Vc::double_v const &y, Vc::double_v const &z, Vc::double_m &c1) {
  vd vfDz(fDz);
  vd vfRmin1(fRmin1);
  vd vfRmin2(fRmin2);
  vd vfRmax1(fRmax1);
  vd vfRmax2(fRmax2);

  vd absZ = Vc::abs(z);
  vd r2(x * x + y * y);
  vd rl(0.5 * (vfRmin2 * (z + vfDz) + vfRmin1 * (vfDz - z)) / vfDz);
  vd rh(0.5 * (vfRmax2 * (z + vfDz) + vfRmax1 * (vfDz - z)) / vfDz);
  c1 = ((absZ <= vfDz) && (r2 >= rl * rl) && (r2 <= rh * rh));
}

//_____________________________________________________________________________
void TGeoConeSeg_v::Contains_v(const StructOfCoord &pointi, Bool_t *isin, Int_t np) const {
  static vd vfPhi1(fPhi1);
  static vd vfPhi2(fPhi2);
  int vectorsize = Vc::double_v::Size;
  int tailsize = np % vectorsize;
  Vc::double_m particleinside;

  for (unsigned int i = 0; i < np - tailsize; i += vectorsize) {
    vd x(&pointi.x[i]); // will copy a certain number of x's into vc vector x
    vd y(&pointi.y[i]);
    vd z(&pointi.z[i]);
    vd dphi(vfPhi2 - vfPhi1);
    vd phi = Vc::atan2(y, x) * 57.295780181884765625f; // 180/phi
    // vd phi = Vc::atan2(y, x) * kRadToDeg;

    phi(phi < 0.f) += 360.f;
    // if (phi < 0 ) phi+=360.;
    vd ddp = phi - vfPhi1;

    ddp(ddp < 0.f) += 360.f;
    // if (ddp < 0) ddp+=360.;
    Vc::double_m c1;
    TGeoConeSeg_v::Contains_v4(x, y, z, c1);

    Vc::double_m c2 = (dphi >= 360.);
    Vc::double_m c3 = (ddp > dphi);
    particleinside = c1 && (c2 || (!c2 && !c3));
    // particleinside = !c1;
    for (unsigned int j = 0; j < vectorsize; ++j)
      isin[i + j] = (bool)particleinside[j];
  }
  // do the tail part for the moment, we just call the old static version
  for (unsigned int i = 0; i < tailsize; ++i) {
    Double_t point[3];

    point[0] = pointi.x[np - tailsize + i];
    point[1] = pointi.y[np - tailsize + i];
    point[2] = pointi.z[np - tailsize + i];

    Double_t t_dphi = fPhi2 - fPhi1;
    Double_t t_phi = atan2(point[1], point[0]) * kRadToDeg;
    if (t_phi < 0)
      t_phi += 360.;
    Double_t t_ddp = t_phi - fPhi1;
    if (t_ddp < 0)
      t_ddp += 360.;
    isin[i] = ((TGeoCone::Contains(point) & (t_ddp <= t_dphi)) | t_dphi >= 360.);
  }
}
#endif

//_____________________________________________________________________________

void TGeoConeSeg_v::Safety_v(const StructOfCoord &pointi, Bool_t in, Double_t *safety, Int_t np) const {
  for (unsigned int i = 0; i < np; i++) {
    Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    safety[i] = TGeoConeSeg::Safety(point, in);
  }
}
//_____________________________________________________________________________
void TGeoConeSeg_v::DistFromInside_v(const StructOfCoord &pointi, const StructOfCoord &diri, Int_t /*iact*/,
                                     const Double_t *step, Double_t * /*safe*/, Double_t *distance, Int_t np) const {
  for (unsigned int i = 0; i < np; i++) {
    Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};

    distance[i] = TGeoConeSeg::DistFromInside(point, dir, 3, step[i], 0);
  }
}

//_____________________________________________________________________________

void TGeoConeSeg_v::DistFromOutside_v(const StructOfCoord &pointi, const StructOfCoord &diri, Int_t /*iact*/,
                                      const Double_t *step, Double_t * /*safe*/, Double_t *distance, Int_t np) const {
  for (unsigned int i = 0; i < np; i++) {
    Double_t point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    Double_t dir[3] = {diri.x[i], diri.y[i], diri.z[i]};

    distance[i] = TGeoConeSeg::DistFromOutside(point, dir, 3, step[i], 0);
  }
}
