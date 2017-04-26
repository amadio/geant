#ifndef THREEVECTOR_H
#define THREEVECTOR_H

/* 
SIMD/SIMT version of HEP3VECTOR of CLHEP (adopted from vecGeom::Vector3D)
*/

#include "base/VecPhys.h"

#include <cstdlib>
#include <ostream>
#include <string>

namespace vecphys {

VECPHYS_HOST_FORWARD_DECLARE(template <typename T> class ThreeVector;);

inline namespace VECPHYS_IMPL_NAMESPACE {

template <typename T>
class ThreeVector 
{
private:
  T dx;
  T dy;
  T dz;

public:
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  ThreeVector(const T x, const T y, const T z)
  {
    dx = x; dy = y; dz = z;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  ThreeVector()
  {
    dx = 0; dy = 0; dz = 0;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  ThreeVector(const T c)
  {
    dx = c; dy = c; dz = c;
  }

  template <typename T1>
  VECCORE_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  ThreeVector(ThreeVector<T1> const &rhs)
  {
    dx = rhs.dx; dy = rhs.dy; dz = rhs.dz;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  ThreeVector &operator=(ThreeVector const &rhs)
  {
    dx = rhs.dx; dy = rhs.dy; dz = rhs.dz;
    return *this;
  }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T &x() { return dx; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const &x() const { return dx; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T &y() { return dy; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const &y() const { return dy; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T &z() { return dz; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T const &z() const { return dz; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetX(T const &x) { dx = x; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetY(T const &y) { dx = y; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void SetZ(T const &z) { dx = z; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void Set(T const &x, T const &y, T const &z) { dx = x; dy = y; dz = z; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  void Set(const T c) { Set(c, c, c); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Perp2() const { return dx * dx + dy * dy; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Perp() const { return math::Sqrt(Perp2()); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Mag2() const { return dx*dx + dy*dy + dz*dz; }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Mag() const { return math::Sqrt(Mag2()); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Phi() const { return math::ATan2(dy, dx); }

  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  T Theta() const { return math::ACos(dz / Mag()); }

  //Dot product
  template <typename T1>
  VECCORE_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  T Dot(ThreeVector<T1> const &p) const { return dx*p.dx + dy*p.dy + dz*p.dz; }

  //Cross product
  template <class T1>
  VECCORE_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  ThreeVector<T> Cross(ThreeVector<T1> const &p) const
  {
    return ThreeVector(dy*p.z()-p.y()*dz, dz*p.x()-p.z()*dx, dx*p.y()-p.x()*dy);
  }

  //Unit : Vector parallel to this, but of length 1.
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
  ThreeVector<T> Unit() const
  {
    ThreeVector<T> output(*this);
    const T mag2 = Mag2();
    output /= math::Sqrt(mag2);
    return output;
  }

  // Rotates reference frame from Uz to newUz (unit vector) (Geant4)
  VECCORE_ATT_HOST_DEVICE 
  VECCORE_FORCE_INLINE
  ThreeVector<T>& RotateUz(ThreeVector<T>& newUz);

// Inplace binary operators
// adopted from vecGeom::Vector3D

#define HEP3VECTOR_BINARY_OP(OPERATOR)                          \
  VECCORE_ATT_HOST_DEVICE                                       \
  VECCORE_FORCE_INLINE                                          \
  ThreeVector<T> &operator OPERATOR(const ThreeVector<T> &p)    \
  {                                                             \
    dx OPERATOR p.dx;                                           \
    dy OPERATOR p.dy;                                           \
    dz OPERATOR p.dz;                                           \
    return *this;                                               \
  }                                                             \
  template <typename T1>                                        \
  VECCORE_ATT_HOST_DEVICE                                       \
  VECCORE_FORCE_INLINE                                          \
  ThreeVector<T> &operator OPERATOR(const ThreeVector<T1> &p)   \
  {                                                             \
    dx OPERATOR p.dx;                                           \
    dy OPERATOR p.dy;                                           \
    dz OPERATOR p.dz;                                           \
    return *this;                                               \
  }                                                             \
  VECCORE_ATT_HOST_DEVICE                                       \
  VECCORE_FORCE_INLINE                                          \
  ThreeVector<T> &operator OPERATOR(const T &c)                 \
  {                                                             \
    dx OPERATOR c;                                              \
    dy OPERATOR c;                                              \
    dz OPERATOR c;                                              \
    return *this;                                               \
  }
  HEP3VECTOR_BINARY_OP(+=)
  HEP3VECTOR_BINARY_OP(-=)
  HEP3VECTOR_BINARY_OP(*=)
  HEP3VECTOR_BINARY_OP(/=)
#undef HEP3VECTOR_BINARY_OP
};

template <typename T>
VECCORE_ATT_HOST_DEVICE 
VECCORE_FORCE_INLINE
ThreeVector<T>& ThreeVector<T>::RotateUz(ThreeVector<T>& newUz)
{
  // newUzVector must be normalized !

  T u1 = newUz.x();
  T u2 = newUz.y();
  T u3 = newUz.z();
  T up = math::Sqrt(u1*u1 + u2*u2);

  Mask_v<T> positive = (up > 0.);
  Mask_v<T> negativeZ = (u3 < 0.);

  T invup = Blend(positive,1.0/up, static_cast<T>(0.0));

  T px = dx;
  T py = dy;
  T pz = dz;

  dx = (u1*u3*px - u2*py)*invup + u1*pz;
  dy = (u2*u3*px + u1*py)*invup + u2*pz;
  dz =    -up*px +                u3*pz;

  dx = Blend(positive, dx, Blend(negativeZ, -px, px));
  dy = Blend(positive, dy, Blend(negativeZ,  py, py));
  dz = Blend(positive, dz, Blend(negativeZ, -pz, pz));

  return *this;
}

template <>
VECCORE_ATT_HOST_DEVICE 
VECCORE_FORCE_INLINE
ThreeVector<double>& ThreeVector<double>::RotateUz(ThreeVector<double>& newUz)
{
  //from CLHEP Hep3Vector::rotateUz
  double u1 = newUz.x();
  double u2 = newUz.y();
  double u3 = newUz.z();
  double up = u1*u1 + u2*u2;

  if (up>0) {
    up = std::sqrt(up);
    double px = dx,  py = dy,  pz = dz;
    dx = (u1*u3*px - u2*py)/up + u1*pz;
    dy = (u2*u3*px + u1*py)/up + u2*pz;
    dz =    -up*px +             u3*pz;
  }
  else if (u3 < 0.) { dx = -dx; dz = -dz; }      // phi=0  teta=pi
  else {};
  return *this;
}

#define HEP3VECTOR_BINARY_OP(OPERATOR, INPLACE)                                          \
template <typename T, typename T1>                                                       \
VECCORE_FORCE_INLINE                                                                     \
VECCORE_ATT_HOST_DEVICE                                                                  \
ThreeVector<T> operator OPERATOR(const ThreeVector<T> &lhs, const ThreeVector<T1> &rhs)  \
{                                                                                        \
  ThreeVector<T> result(lhs);                                                            \
  result INPLACE rhs;                                                                    \
  return result;                                                                         \
}                                                                                        \
template <typename T, typename ScalarT>                                                  \
VECCORE_FORCE_INLINE                                                                     \
VECCORE_ATT_HOST_DEVICE                                                                  \
ThreeVector<T> operator OPERATOR(ThreeVector<T> const &lhs, const ScalarT rhs)           \
{                                                                                        \
  ThreeVector<T> result(lhs);                                                            \
  result INPLACE rhs;                                                                    \
  return result;                                                                         \
}                                                                                        \
template <typename T, typename ScalarT>                                                  \
VECCORE_FORCE_INLINE                                                                     \
VECCORE_ATT_HOST_DEVICE                                                                  \
ThreeVector<T> operator OPERATOR(const ScalarT lhs, ThreeVector<T> const &rhs)	         \
{                                                                                        \
  ThreeVector<T> result(lhs);                                                            \
  result INPLACE rhs;                                                                    \
  return result;                                                                         \
}
HEP3VECTOR_BINARY_OP(+, +=)
HEP3VECTOR_BINARY_OP(-, -=)
HEP3VECTOR_BINARY_OP(*, *=)
HEP3VECTOR_BINARY_OP(/, /=)
#undef HEP3VECTOR_BINARY_OP

VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
bool operator==(ThreeVector<Real_t> const &lhs, ThreeVector<Real_t> const &rhs)
{
  return math::Abs(lhs.x() - rhs.x()) < 0. && 
         math::Abs(lhs.y() - rhs.y()) < 0. && 
         math::Abs(lhs.z() - rhs.z()) < 0. ;
}

VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
ThreeVector<bool> operator!=(ThreeVector<Real_t> const &lhs, ThreeVector<Real_t> const &rhs)
{
  return !(lhs == rhs);
}

template <typename T>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
ThreeVector<T> operator-(ThreeVector<T> const &v)
{
  return ThreeVector<T>(-v.x(), -v.y(), -v.z());
}

VECCORE_ATT_HOST_DEVICE
VECCORE_FORCE_INLINE
ThreeVector<bool> operator!(ThreeVector<bool> const &v)
{
  return ThreeVector<bool>(!v.x(), !v.y(), !v.z());
}

template <typename T>
std::ostream &operator<<(std::ostream &os, ThreeVector<T> const &v)
{
  os << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
  return os;
}

} // end namespace impl
} // end namespace vecphys

#endif
