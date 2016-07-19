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
void CreateParticle0022() {

   // Creating delta(1900)-
   new Particle("delta(1900)-", 11112, 1, "Unknown", 100, -1, 1.9, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1700)-
   new Particle("delta(1700)-", 11114, 1, "Unknown", 100, -1, 1.7, 0.3, -100, 0, -100, -1, -1);

   // Creating delta(1930)-
   new Particle("delta(1930)-", 11116, 1, "Unknown", 100, -1, 1.96, 0.36, -100, 0, -100, -1, -1);

   // Creating delta(1900)0
   new Particle("delta(1900)0", 11212, 1, "Unknown", 100, 0, 1.9, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1930)0
   new Particle("delta(1930)0", 11216, 1, "Unknown", 100, 0, 1.96, 0.36, -100, 0, -100, -1, -1);

   // Creating N(1440)0
   new Particle("N(1440)0", 12112, 1, "Unknown", 100, 0, 1.44, 0.3, -100, 0, -100, -1, -1);

   // Creating delta(1700)0
   new Particle("delta(1700)0", 12114, 1, "Unknown", 100, 0, 1.7, 0.3, -100, 0, -100, -1, -1);

   // Creating N(1680)0
   new Particle("N(1680)0", 12116, 1, "Unknown", 100, 0, 1.685, 0.13, -100, 0, -100, -1, -1);

   // Creating N(1990)0
   new Particle("N(1990)0", 12118, 1, "Unknown", 100, 0, 1.95, 0.555, -100, 0, -100, -1, -1);

   // Creating delta(1900)+
   new Particle("delta(1900)+", 12122, 1, "Unknown", 100, 1, 1.9, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1930)+
   new Particle("delta(1930)+", 12126, 1, "Unknown", 100, 1, 1.96, 0.36, -100, 0, -100, -1, -1);

   // Creating N(1440)+
   new Particle("N(1440)+", 12212, 1, "Unknown", 100, 1, 1.44, 0.3, -100, 0, -100, -1, -1);

   // Creating delta(1700)+
   new Particle("delta(1700)+", 12214, 1, "Unknown", 100, 1, 1.7, 0.3, -100, 0, -100, -1, -1);

   // Creating N(1680)+
   new Particle("N(1680)+", 12216, 1, "Unknown", 100, 1, 1.685, 0.13, -100, 0, -100, -1, -1);

   // Creating N(1990)+
   new Particle("N(1990)+", 12218, 1, "Unknown", 100, 1, 1.95, 0.555, -100, 0, -100, -1, -1);

   // Creating delta(1900)++
   new Particle("delta(1900)++", 12222, 1, "Unknown", 100, 2, 1.9, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1700)++
   new Particle("delta(1700)++", 12224, 1, "Unknown", 100, 2, 1.7, 0.3, -100, 0, -100, -1, -1);

   // Creating delta(1930)++
   new Particle("delta(1930)++", 12226, 1, "Unknown", 100, 2, 1.96, 0.36, -100, 0, -100, -1, -1);

   // Creating sigma(1660)-
   new Particle("sigma(1660)-", 13112, 1, "Unknown", 100, -1, 1.66, 0.1, -100, 0, -100, -1, -1);

   // Creating sigma(1670)-
   new Particle("sigma(1670)-", 13114, 1, "Unknown", 100, -1, 1.67, 0.06, -100, 0, -100, -1, -1);

   // Creating sigma(1915)-
   new Particle("sigma(1915)-", 13116, 1, "Unknown", 100, -1, 1.915, 0.12, -100, 0, -100, -1, -1);

   // Creating lambda(1405)
   new Particle("lambda(1405)", 13122, 1, "Unknown", 100, 0, 1.4051, 0.05, -100, 0, -100, -1, -1);

   // Creating lambda(1690)
   new Particle("lambda(1690)", 13124, 1, "Unknown", 100, 0, 1.69, 0.06, -100, 0, -100, -1, -1);

   // Creating lambda(1830)
   new Particle("lambda(1830)", 13126, 1, "Unknown", 100, 0, 1.83, 0.095, -100, 0, -100, -1, -1);

   // Creating sigma(1660)0
   new Particle("sigma(1660)0", 13212, 1, "Unknown", 100, 0, 1.66, 0.1, -100, 0, -100, -1, -1);

   // Creating sigma(1670)0
   new Particle("sigma(1670)0", 13214, 1, "Unknown", 100, 0, 1.67, 0.06, -100, 0, -100, -1, -1);

   // Creating sigma(1915)0
   new Particle("sigma(1915)0", 13216, 1, "Unknown", 100, 0, 1.915, 0.12, -100, 0, -100, -1, -1);

   // Creating sigma(1660)+
   new Particle("sigma(1660)+", 13222, 1, "Unknown", 100, 1, 1.66, 0.1, -100, 0, -100, -1, -1);

   // Creating sigma(1670)+
   new Particle("sigma(1670)+", 13224, 1, "Unknown", 100, 1, 1.67, 0.06, -100, 0, -100, -1, -1);

   // Creating sigma(1915)+
   new Particle("sigma(1915)+", 13226, 1, "Unknown", 100, 1, 1.915, 0.12, -100, 0, -100, -1, -1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
