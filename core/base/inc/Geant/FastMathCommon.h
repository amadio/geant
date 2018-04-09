#ifndef GEANTV_FASTMATHCOMMON_H
#define GEANTV_FASTMATHCOMMON_H

#include <cstdint>

namespace geant {
namespace details {

union ieee754 {
  inline ieee754(){};
  inline ieee754(double thed) { d = thed; };
  inline ieee754(uint64_t thell) { ll = thell; };
  inline ieee754(float thef) { f[0] = thef; };
  inline ieee754(uint32_t thei) { i[0] = thei; };
  double d;
  float f[2];
  uint32_t i[2];
  uint64_t ll;
  uint16_t s[4];
};

inline double uint642dp(uint64_t ll)
{
  ieee754 tmp;
  tmp.ll = ll;
  return tmp.d;
}

inline uint32_t sp2uint32(float x)
{
  ieee754 tmp;
  tmp.f[0] = x;
  return tmp.i[0];
}

inline double fpfloor(const double x)
{
  // no problem since exp is defined between -708 and 708. Int is enough for it!
  int32_t ret = int32_t(x);
  ret -= (sp2uint32(x) >> 31);
  return ret;
}

inline uint64_t dp2uint64(double x)
{
  ieee754 tmp;
  tmp.d = x;
  return tmp.ll;
}

inline double getMantExponent(const double x, double &fe)
{

  uint64_t n = dp2uint64(x);

  // Shift to the right up to the beginning of the exponent.
  // Then with a mask, cut off the sign bit
  uint64_t le = (n >> 52);

  // chop the head of the number: an int contains more than 11 bits (32)
  int32_t e = le; // This is important since sums on uint64_t do not vectorise
  fe        = e - 1023;

  // This puts to 11 zeroes the exponent
  n &= 0x800FFFFFFFFFFFFFULL;
  // build a mask which is 0.5, i.e. an exponent equal to 1022
  // which means *2, see the above +1.
  const uint64_t p05 = 0x3FE0000000000000ULL; // dp2uint64(0.5);
  n |= p05;

  return uint642dp(n);
}
}
}

#endif // GEANTV_FASTMATHCOMMON_H
