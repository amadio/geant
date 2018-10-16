#ifndef GEANT_MATH
#define GEANT_MATH

#include "Geant/Config.h"

// include VecCore's math ...
#include "VecCore/VecCore"
#include <VecMath/Math.h>

namespace Math {
template <typename T>
VECCORE_ATT_HOST_DEVICE
inline T Min(T const &val1, T const &val2)
{
  return vecCore::math::Min(val1, val2);
}
template <typename T>
VECCORE_ATT_HOST_DEVICE
inline T Max(T const &val1, T const &val2)
{
  return vecCore::math::Max(val1, val2);
}
template <typename T>
VECCORE_ATT_HOST_DEVICE
inline T Sqrt(T const &val)
{
  return vecCore::math::Sqrt(val);
}
template <typename T>
VECCORE_ATT_HOST_DEVICE
inline T Abs(T const &val)
{
  return vecCore::math::Abs(val);
}
template <typename T>
VECCORE_ATT_HOST_DEVICE
inline bool AreEqualAbs(T const &val1, T const &val2, T const &epsilon)
{
  return (vecCore::math::Abs(val1 - val2) < epsilon);
}
template <typename T>
VECCORE_ATT_HOST_DEVICE
inline bool AreEqualRel(T const &val1, T const &val2, T const &relPrec)
{
  return (vecCore::math::Abs(val1 - val2) < 0.5 * relPrec * (vecCore::math::Abs(val1) + vecCore::math::Abs(val2)));
}

//  template <typename T> VECCORE_ATT_HOST_DEVICE inline T Normalize(T const &val[3]) { return
//  vecCore::math::Normalize(val); }
//  VECCORE_ATT_HOST_DEVICE vecCore::math::Precision inline TwoPi() { return vecCore::math::TwoPi(); }

// From TMath.cxx ....
VECCORE_ATT_HOST_DEVICE
inline float Normalize(float v[3])
{
  // Normalize a vector v in place.
  // Returns the norm of the original vector.

  float d = Sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  if (d != 0) {
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
  }
  return d;
}
VECCORE_ATT_HOST_DEVICE
inline double Normalize(double v[3])
{
  // Normalize a vector v in place.
  // Returns the norm of the original vector.
  // This implementation (thanks Kevin Lynch <krlynch@bu.edu>) is protected
  // against possible overflows.

  // Find the largest element, and divide that one out.

  double av0 = Abs(v[0]), av1 = Abs(v[1]), av2 = Abs(v[2]);

  double amax, foo, bar;
  // 0 >= {1, 2}
  if (av0 >= av1 && av0 >= av2) {
    amax = av0;
    foo  = av1;
    bar  = av2;
  }
  // 1 >= {0, 2}
  else if (av1 >= av0 && av1 >= av2) {
    amax = av1;
    foo  = av0;
    bar  = av2;
  }
  // 2 >= {0, 1}
  else {
    amax = av2;
    foo  = av0;
    bar  = av1;
  }

  if (amax == 0.0) return 0.;

  double foofrac = foo / amax, barfrac = bar / amax;
  double d = amax * Sqrt(1. + foofrac * foofrac + barfrac * barfrac);

  v[0] /= d;
  v[1] /= d;
  v[2] /= d;
  return d;
}
constexpr VECCORE_ATT_HOST_DEVICE inline double TwoPi()
{
  return 2 * 3.14159265358979323846;
}
constexpr VECCORE_ATT_HOST_DEVICE inline double Pi()
{
  return 3.14159265358979323846;
}

template <typename R>
R Exp(R x)
{
  return vecMath::FastExp(x);
}

template <typename R>
R Log(R x)
{
  return vecMath::FastLog(x);
}

template <typename R>
R Log10(R x)
{
  return vecMath::FastLog(x) * 0.43429448190325182;
}

template <typename R>
R Sin(R x)
{
  return vecMath::FastSin(x);
}

template <typename R>
R Cos(R x)
{
  return vecMath::FastCos(x);
}

template <typename R>
void SinCos(R x, R &s, R &c)
{
  vecMath::FastSinCos(x, s, c);
}

template <typename R>
R Pow(R x, R n)
{
  return vecMath::FastPow(x, n);
}

using vecMath::IntPow;

/**
 * @brief Rotate vector u,v,w to labframe defined by vector u1, u2, u3
 * @tparam double or Real_v type from veccore
 * @param[in,out] u, v, w - vector being rotated to lab frame
 * @param[in] u1, u2, u3 - direction equal to z-direction of this scattering frame
 */

template <typename R>
inline void RotateToLabFrame(R &u, R &v, R &w, R u1, R u2, R u3)
{
  R up                = u1 * u1 + u2 * u2;
  vecCore::Mask<R> m1 = up > 0.0;
  up                  = Math::Sqrt(up);
  R px                = u;
  R py                = v;
  R pz                = w;
  vecCore__MaskedAssignFunc(u, m1, (u3 * px * u1 - u2 * py) * (1.0 / up) + u1 * pz);
  vecCore__MaskedAssignFunc(v, m1, (u3 * px * u2 + u1 * py) * (1.0 / up) + u2 * pz);
  vecCore__MaskedAssignFunc(w, m1, -up * px + u3 * pz);
  vecCore::Mask<R> m2 = !m1 && u3 < 0.;
  if (!vecCore::MaskEmpty(m2)) { // Up zero AND u3 negative
    vecCore::MaskedAssign(u, m2, -u);
    vecCore::MaskedAssign(w, m2, -w);
  }
}
} // namespace Math

#endif // GEANT_MATH
