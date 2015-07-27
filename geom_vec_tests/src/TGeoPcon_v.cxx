// Include Vc library and definitions (here or in the header file), math library, TGeoShape library
#include <iostream>

#include "TGeoPcon_v.h"

#include "Util.h"

#include "base/Global.h"
using vecgeom::kPi;

#ifdef VEC_EXTENSIONS
#include "Vc/vector.h"

typedef Vc::double_v vd;  // short for vector double
typedef Vc::double_m vdm; // short for double mask
typedef Vc::int_v vi;     // short for vector integer

static vdm true_m(true);

struct Foo {
  static Vc::double_m IsSameWithinTolerance(Vc::double_v const &a, Vc::double_v const &b) {
    Vc::double_m c = Vc::abs(a - b) < 1.E-10;
    return c;
  }
};
#endif

//_____________________________________________________________________________
#ifdef VEC_EXTENSIONS
void TGeoPcon_v::Contains_v(const StructOfCoord &pointi, Bool_t *isin, Int_t np) const {
  // declare variables
  static vd fZ1_v(fZ[0]); // vector with the first z-plane position in all its components
  static vd fZN_v(fZ[fNz - 1]);
  static vd a360_v(360.);
  static vd radToDeg_v(180. / kPi);
  static vd fPhi1_v(fPhi1);
  static vd fDphi_v(fDphi);
  static vd e10_v(1E-10);

  // copmute tailsize
  int tailsize = np % Vc::double_v::Size;

  // Vc part
  for (unsigned int i = 0; i < np - tailsize; i += Vc::double_v::Size) {
    vd particleInside_v(Vc::Zero);

    vd x_v(&pointi.x[i]); // will copy a certain number of x's into vc vector x;
    vd y_v(&pointi.y[i]);
    vd z_v(&pointi.z[i]);

    // check if the particles are inside in the z-direction
    vdm zInside_m = (z_v > fZ1_v) & (z_v < fZN_v);

    if (zInside_m.isEmpty()) // if all the particles are outside
    {
      //  isin[i] = zInside_m;
      for (unsigned int j = 0; j < Vc::double_v::Size; ++j) {
        isin[i + j] = (bool)particleInside_v[j];
      }
      continue;
    }

    // bound between two z-planes
    vd fZl_v(Vc::Zero);    // z-position of the LBP
    vd fZh_v(Vc::Zero);    // z-position of the HBP
    vd fRminl_v(Vc::Zero); // inner radius of the LBP
    vd fRminh_v(Vc::Zero); // inner radius of the HBP
    vd fRmaxl_v(Vc::Zero); // outer radius of the LBP
    vd fRmaxh_v(Vc::Zero); // outer radius of the HBP

    for (unsigned int j = 0; j < Vc::double_v::Size; ++j) {
      Int_t izl = 0;       // lower boundary plane (LBP)
      Int_t izh = fNz - 1; // higher boundary plane (HBP)
      Int_t izt = (fNz - 1) / 2;

      while ((izh - izl) > 1) {
        if (z_v[j] > fZ[izt])
          izl = izt;
        else
          izh = izt;
        izt = (izh + izl) >> 1;
      }
      fZl_v[j] = fZ[izl];
      fZh_v[j] = fZ[izh];
      fRminl_v[j] = fRmin[izl];
      fRminh_v[j] = fRmin[izh];
      fRmaxl_v[j] = fRmax[izl];
      fRmaxh_v[j] = fRmax[izh];
    }

    // compute r squared
    vd r2_v = x_v * x_v + y_v * y_v;

    // compute rmin, rmax
    vd rmin_v(Vc::Zero);
    vd rmax_v(Vc::Zero);
    vdm same_m = Foo::IsSameWithinTolerance(fZl_v, fZh_v) && Foo::IsSameWithinTolerance(z_v, fZl_v);
    vd dz_v = fZh_v - fZl_v;
    vd dz1_v = z_v - fZl_v;

    rmin_v(same_m) = Vc::min(fRminl_v, fRminh_v);
    rmax_v(same_m) = Vc::max(fRmaxl_v, fRmaxh_v);
    rmin_v(!same_m) = (fRminl_v * (dz_v - dz1_v) + fRminh_v * dz1_v) / dz_v;
    rmax_v(!same_m) = (fRmaxl_v * (dz_v - dz1_v) + fRmaxh_v * dz1_v) / dz_v;

    // check if the particles are radially in the volume
    vdm rInside_m = (r2_v > rmin_v * rmin_v) && (r2_v < rmax_v * rmax_v);
    if (rInside_m.isEmpty()) // if all the particles are outside
    {
      // isin[i] = rInside_m;
      for (unsigned int j = 0; j < Vc::double_v::Size; ++j) {
        isin[i + j] = (bool)particleInside_v[j];
      }
      continue;
    }
    // check if the particles are inside the phi-range of the volume
    if (TGeoShape::IsSameWithinTolerance(fDphi, 360)) {
      //  isin[i] = zInside_m && rInside_m;
      particleInside_v(zInside_m && rInside_m) = 1.;
      for (unsigned int j = 0; j < Vc::double_v::Size; ++j) {
        isin[i + j] = (bool)particleInside_v[j];
      }
      continue;
    }
    // THIS SECTION IS NOT TESTED BY THE BENCHMARK (THE VOLUME HAS 360 DEGREES)

    vd phi_v = Vc::atan2(y_v, x_v) * radToDeg_v;
    // a360_v.setZero(phi_v < Vc::Zero);
    // phi_v += a360_v;
    phi_v(phi_v < Vc::Zero) += a360_v;

    vd ddp_v = phi_v - fPhi1_v;
    ddp_v(ddp_v < Vc::Zero) += a360_v;

    vdm phiInside_m = (ddp_v <= fDphi_v) || (r2_v < e10_v);
    // isin[i] = zInside_m && rInside_m && phiInside_m;

    particleInside_v(zInside_m && rInside_m && phiInside_m) = 1.;
    for (unsigned int j = 0; j < Vc::double_v::Size; ++j) {
      isin[i + j] = (bool)particleInside_v[j];
    }
  }

  // tail part
  for (unsigned int i = 0; i < tailsize; ++i) {
    double point[3] = {pointi.x[np - tailsize + i], pointi.y[np - tailsize + i], pointi.z[np - tailsize + i]};
    isin[np - tailsize + i] = TGeoPcon::Contains(point);
  }
}
#else
void TGeoPcon_v::Contains_v(const StructOfCoord &pointi, Bool_t *isin, Int_t np) const {
  for (unsigned int i = 0; i < np; i++) {
    double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    isin[i] = TGeoPcon::Contains(point);
  }
}
#endif

//_____________________________________________________________________________

void TGeoPcon_v::Safety_v(const StructOfCoord &pointi, Bool_t in, double *safety, Int_t np) const {
  for (unsigned int i = 0; i < np; i++) {
    double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    safety[i] = TGeoPcon::Safety(point, in);
  }
}
//_____________________________________________________________________________
void TGeoPcon_v::DistFromInside_v(const StructOfCoord &pointi, const StructOfCoord &diri, Int_t /*iact*/,
                                  const double *step, double * /*safe*/, double *distance, Int_t np) const {
  for (unsigned int i = 0; i < np; i++) {
    double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    double dir[3] = {diri.x[i], diri.y[i], diri.z[i]};

    distance[i] = TGeoPcon::DistFromInside(point, dir, 3, step[i], 0);
  }
}

//_____________________________________________________________________________

void TGeoPcon_v::DistFromOutside_v(const StructOfCoord &pointi, const StructOfCoord &diri, Int_t /*iact*/,
                                   const double *step, double * /*safe*/, double *distance, Int_t np) const {
  for (unsigned int i = 0; i < np; i++) {
    double point[3] = {pointi.x[i], pointi.y[i], pointi.z[i]};
    double dir[3] = {diri.x[i], diri.y[i], diri.z[i]};

    distance[i] = TGeoPcon::DistFromOutside(point, dir, 3, step[i], 0);
  }
}
