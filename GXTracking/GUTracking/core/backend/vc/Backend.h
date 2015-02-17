/// \file vc/backend.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)
//
// temporary clone for GUTracking based on VecGeom - syjun
//
#ifndef VECPHYS_BACKEND_VCBACKEND_H
#define VECPHYS_BACKEND_VCBACKEND_H

#include "base/Global.h"
#include "backend/scalar/Backend.h"

#include <Vc/Vc>

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

struct kVc {

  typedef Vc::int_v                   Int_t;
  typedef Vc::Vector<Precision>       Double_t;
  typedef Vc::Vector<Precision>::Mask Bool_t;
  typedef Vc::Vector<Precision>       Index_t;

  const static Bool_t kTrue;
  const static Bool_t kFalse;
  constexpr static bool early_returns = false;

  constexpr static int kSize = kVc::Double_t::Size;
  const static Double_t kOne;
  const static Double_t kZero;
};

typedef kVc::Int_t       VcInt;
typedef kVc::Double_t    VcPrecision;
typedef kVc::Bool_t      VcBool;
typedef kVc::Int_t       VcInside;

template <typename Type>
VECPHYS_INLINE
void CondAssign(typename Vc::Vector<Type>::Mask const &cond,
                Vc::Vector<Type> const &thenval,
                Vc::Vector<Type> const &elseval,
                Vc::Vector<Type> *const output) {
  (*output)(cond) = thenval;
  (*output)(!cond) = elseval;
}

template <typename Type>
VECPHYS_INLINE
void CondAssign(typename Vc::Vector<Type>::Mask const &cond,
                Type const &thenval,
                Type const &elseval,
                Vc::Vector<Type> *const output) {
  (*output)(cond) = thenval;
  (*output)(!cond) = elseval;
}

template <typename Type>
VECPHYS_INLINE
void MaskedAssign(typename Vc::Vector<Type>::Mask const &cond,
                  Vc::Vector<Type> const &thenval,
                  Vc::Vector<Type> *const output) {
  (*output)(cond) = thenval;
}

template <typename Type>
VECPHYS_INLINE
void MaskedAssign(typename Vc::Vector<Type>::Mask const &cond,
                  Type const &thenval,
                  Vc::Vector<Type> *const output) {
  (*output)(cond) = thenval;
}

VECPHYS_INLINE
void MaskedAssign(VcBool const &cond,
                  const Inside_t thenval,
                  VcInside *const output) {
  (*output)(VcInside::Mask(cond)) = thenval;
}


VECPHYS_INLINE
bool IsFull(VcBool const &cond) {
  return cond.isFull();
}

VECPHYS_INLINE
bool Any(VcBool const &cond) {
  return !cond.isEmpty();
}

VECPHYS_INLINE
bool IsEmpty(VcBool const &cond) {
  return cond.isEmpty();
}

VECPHYS_INLINE
VcPrecision Abs(VcPrecision const &val) {
  return Vc::abs(val);
}

VECPHYS_INLINE
VcPrecision Sqrt(VcPrecision const &val) {
  return Vc::sqrt(val);
}

VECPHYS_INLINE
VcPrecision ATan2(VcPrecision const &y, VcPrecision const &x) {
  return Vc::atan2(y, x);
}


VECPHYS_INLINE
VcPrecision sin(VcPrecision const &x) {
  return Vc::sin(x);
}

VECPHYS_INLINE
VcPrecision cos(VcPrecision const &x) {
  return Vc::cos(x);
}

VECPHYS_INLINE
VcPrecision tan(VcPrecision const &radians) {
  // apparently Vc does not have a tan function
  //  return Vc::tan(radians);
  // emulating it for the moment
  VcPrecision s,c;
  Vc::sincos(radians,&s,&c);
  return s/c;
}

VECPHYS_INLINE
Precision Pow(Precision const &x, Precision arg) {
   return std::pow(x,arg);
}

VECPHYS_INLINE
VcPrecision Min(VcPrecision const &val1, VcPrecision const &val2) {
  return Vc::min(val1, val2);
}

VECPHYS_INLINE
VcPrecision Max(VcPrecision const &val1, VcPrecision const &val2) {
  return Vc::max(val1, val2);
}

VECPHYS_INLINE
VcInt Min(VcInt const &val1, VcInt const &val2) {
  return Vc::min(val1, val2);
}

VECPHYS_INLINE
VcInt Max(VcInt const &val1, VcInt const &val2) {
  return Vc::max(val1, val2);
}


VECPHYS_INLINE
VcPrecision Floor( VcPrecision const &val ){
  return Vc::floor( val );
}

VECPHYS_INLINE
  VcPrecision UniformRandom(Random_t* state,  VcInt val){
  return kVc::Double_t::Random( );
}

} // End inline namespace

} // End global namespace


#endif // VECPHYS_BACKEND_VCBACKEND_H
