#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
#ifdef GEANT_NVCC
#include "base/Vector.h"
template <typename T>
using vector = vecgeom::Vector<T>;
#else
using std::vector;
#endif
namespace geant {
   inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void CreateParticle0027() {

   // Creating f0(600)
   new Particle("f0(600)", 9000221, 1, "Unknown", 100, 0, 0.8, 0.8, -100, 0, -100, -1, -1);

   // Creating f0(980)
   new Particle("f0(980)", 9010221, 1, "Unknown", 100, 0, 0.98, 0.07, -100, 0, -100, -1, -1);

   // Creating eta(1405)
   new Particle("eta(1405)", 9020221, 1, "Unknown", 100, 0, 1.4098, 0.0511, -100, 0, -100, -1, -1);

   // Creating f0(1500)
   new Particle("f0(1500)", 9030221, 1, "Unknown", 100, 0, 1.505, 0.109, -100, 0, -100, -1, -1);

   // Creating f2(1810)
   new Particle("f2(1810)", 9030225, 1, "Unknown", 100, 0, 1.815, 0.197, -100, 0, -100, -1, -1);

   // Creating f2(2010)
   new Particle("f2(2010)", 9060225, 1, "Unknown", 100, 0, 2.01, 0.2, -100, 0, -100, -1, -1);

   // Creating Cherenkov
   new Particle("Cherenkov", 50000050, 1, "Unknown", 100, 0, 0, 0, -100, 0, -100, -1, -1);

   // Creating ChargedRootino
   new Particle("ChargedRootino", 50000052, 1, "Unknown", 100, -0.333333, 0, 0, -100, 0, -100, -1, -1);

   // Creating GenericIon
   new Particle("GenericIon", 50000060, 1, "Unknown", 100, 0.333333, 0.938272, 0, -100, 0, -100, -1, -1);

   // Creating N(2220)0
   new Particle("N(2220)0", 100002110, 1, "Unknown", 100, 0, 2.25, 0.4, -100, 0, -100, -1, -1);

   // Creating N(2220)+
   new Particle("N(2220)+", 100002210, 1, "Unknown", 100, 1, 2.25, 0.4, -100, 0, -100, -1, -1);

   // Creating N(2250)0
   new Particle("N(2250)0", 100012110, 1, "Unknown", 100, 0, 2.275, 0.5, -100, 0, -100, -1, -1);

   // Creating N(2250)+
   new Particle("N(2250)+", 100012210, 1, "Unknown", 100, 1, 2.275, 0.5, -100, 0, -100, -1, -1);

   // Creating Deuteron
   new Particle("Deuteron", 1000010020, 1, "ion", 100, 1, 1.87106, 0, -100, 0, -100, -1, -1);

   // Creating Triton
   new Particle("Triton", 1000010030, 1, "ion", 100, 1, 2.80941, 1.6916e-33, -100, 0, -100, -1, -1);

   // Creating HE3
   new Particle("HE3", 1000020030, 1, "ion", 100, 2, 2.80941, 0, -100, 0, -100, -1, -1);

   // Creating Alpha
   new Particle("Alpha", 1000020040, 1, "ion", 100, 2, 3.7284, 1.6916e-33, -100, 0, -100, -1, -1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
