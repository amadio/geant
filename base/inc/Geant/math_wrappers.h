#ifndef GEANT_MATH
#define GEANT_MATH

#include "Geant/Config.h"

// include VecCore's math ...
#include "VecCore/VecCore"

namespace Math {
   template <typename T> VECCORE_ATT_HOST_DEVICE inline T Min(T const &val1, T const &val2) { return vecCore::math::Min(val1,val2); }
   template <typename T> VECCORE_ATT_HOST_DEVICE inline T Max(T const &val1, T const &val2) { return vecCore::math::Max(val1,val2); }
   template <typename T> VECCORE_ATT_HOST_DEVICE inline T Sqrt(T const &val) { return vecCore::math::Sqrt(val); }
   template <typename T> VECCORE_ATT_HOST_DEVICE inline T Abs(T const &val) { return vecCore::math::Abs(val); }
   template <typename T> VECCORE_ATT_HOST_DEVICE inline T Cos(T const &val) { return vecCore::math::Cos(val); }
   template <typename T> VECCORE_ATT_HOST_DEVICE inline T Sin(T const &val) { return vecCore::math::Sin(val); }
   template <typename T> VECCORE_ATT_HOST_DEVICE inline bool AreEqualAbs(T const &val1, T const &val2, T const &epsilon) 
     { return ( vecCore::math::Abs(val1-val2) < epsilon ); }
   template <typename T> VECCORE_ATT_HOST_DEVICE inline bool AreEqualRel(T const &val1, T const &val2, T const &relPrec) 
     { return ( vecCore::math::Abs(val1-val2) < 0.5*relPrec*(vecCore::math::Abs(val1)+vecCore::math::Abs(val2)) ); }

//  template <typename T> VECCORE_ATT_HOST_DEVICE inline T Normalize(T const &val[3]) { return vecCore::math::Normalize(val); }
//  VECCORE_ATT_HOST_DEVICE vecCore::math::Precision inline TwoPi() { return vecCore::math::TwoPi(); }

// From TMath.cxx ....
   VECCORE_ATT_HOST_DEVICE
   inline float Normalize(float v[3])
   {
      // Normalize a vector v in place.
      // Returns the norm of the original vector.

      float d = Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
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
      if( av0 >= av1 && av0 >= av2 ) {
         amax = av0;
         foo = av1;
         bar = av2;
      }
      // 1 >= {0, 2}
      else if (av1 >= av0 && av1 >= av2) {
         amax = av1;
         foo = av0;
         bar = av2;
      }
      // 2 >= {0, 1}
      else {
         amax = av2;
         foo = av0;
         bar = av1;
      }

      if (amax == 0.0)
         return 0.;

      double foofrac = foo/amax, barfrac = bar/amax;
      double d = amax * Sqrt(1.+foofrac*foofrac+barfrac*barfrac);

      v[0] /= d;
      v[1] /= d;
      v[2] /= d;
      return d;
   }
   constexpr VECCORE_ATT_HOST_DEVICE inline double TwoPi() { return 2*3.14159265358979323846; }
   constexpr VECCORE_ATT_HOST_DEVICE inline double Pi() { return 3.14159265358979323846; }
   

}

#endif // GEANT_MATH
