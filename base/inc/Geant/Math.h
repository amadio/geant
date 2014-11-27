#ifndef GEANT_MATH
#define GEANT_MATH

#include "Geant/Config.h"

// include VecGeom's math ...
#include "backend/Backend.h"

namespace Math {
   template <typename T> GEANT_CUDA_BOTH_CODE inline T Min(T const &val1, T const &val2) { return VECGEOM_NAMESPACE::Min(val1,val2); }
   template <typename T> GEANT_CUDA_BOTH_CODE inline T Max(T const &val1, T const &val2) { return VECGEOM_NAMESPACE::Max(val1,val2); }
   template <typename T> GEANT_CUDA_BOTH_CODE inline T Sqrt(T const &val) { return VECGEOM_NAMESPACE::Sqrt(val); }
   template <typename T> GEANT_CUDA_BOTH_CODE inline T Abs(T const &val) { return VECGEOM_NAMESPACE::Abs(val); }
   template <typename T> GEANT_CUDA_BOTH_CODE inline T Cos(T const &val) { return VECGEOM_NAMESPACE::cos(val); }
   template <typename T> GEANT_CUDA_BOTH_CODE inline T Sin(T const &val) { return VECGEOM_NAMESPACE::sin(val); }

//  template <typename T> GEANT_CUDA_BOTH_CODE inline T Normalize(T const &val[3]) { return VECGEOM_NAMESPACE::Normalize(val); }
//  GEANT_CUDA_BOTH_CODE VECGEOM_NAMESPACE::Precision inline TwoPi() { return VECGEOM_NAMESPACE::TwoPi(); }

// From TMath.cxx ....
   GEANT_CUDA_BOTH_CODE
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
   GEANT_CUDA_BOTH_CODE
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
   constexpr GEANT_CUDA_BOTH_CODE inline VECGEOM_NAMESPACE::Precision TwoPi() { return 2*3.14159265358979323846; }
   constexpr GEANT_CUDA_BOTH_CODE inline VECGEOM_NAMESPACE::Precision Pi() { return 3.14159265358979323846; }
   

}

#endif // GEANT_MATH
