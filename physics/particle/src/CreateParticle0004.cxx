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
void CreateParticle0004() {

   // Creating xi(1820)0_bar
   new Particle("xi(1820)0_bar", -13324, 0, "Unknown", 100, 0, 1.823, 0.024, 100, 100, 0, 100, 1);

   // Creating xi(2030)-_bar
   new Particle("xi(2030)-_bar", -13316, 0, "Unknown", 100, 1, 2.025, 0.02, 100, 100, 0, 100, 1);

   // Creating xi(1820)-_bar
   new Particle("xi(1820)-_bar", -13314, 0, "Unknown", 100, 1, 1.823, 0.024, 100, 100, 0, 100, 1);

   // Creating sigma(1915)+_bar
   new Particle("sigma(1915)+_bar", -13226, 0, "Unknown", 100, -1, 1.915, 0.12, 100, 100, 0, 100, 1);

   // Creating sigma(1670)+_bar
   new Particle("sigma(1670)+_bar", -13224, 0, "Unknown", 100, -1, 1.67, 0.06, 100, 100, 0, 100, 1);

   // Creating sigma(1660)+_bar
   new Particle("sigma(1660)+_bar", -13222, 0, "Unknown", 100, -1, 1.66, 0.1, 100, 100, 0, 100, 1);

   // Creating sigma(1915)0_bar
   new Particle("sigma(1915)0_bar", -13216, 0, "Unknown", 100, 0, 1.915, 0.12, 100, 100, 0, 100, 1);

   // Creating sigma(1670)0_bar
   new Particle("sigma(1670)0_bar", -13214, 0, "Unknown", 100, 0, 1.67, 0.06, 100, 100, 0, 100, 1);

   // Creating sigma(1660)0_bar
   new Particle("sigma(1660)0_bar", -13212, 0, "Unknown", 100, 0, 1.66, 0.1, 100, 100, 0, 100, 1);

   // Creating lambda(1830)_bar
   new Particle("lambda(1830)_bar", -13126, 0, "Unknown", 100, 0, 1.83, 0.095, 100, 100, 0, 100, 1);

   // Creating lambda(1690)_bar
   new Particle("lambda(1690)_bar", -13124, 0, "Unknown", 100, 0, 1.69, 0.06, 100, 100, 0, 100, 1);

   // Creating lambda(1405)_bar
   new Particle("lambda(1405)_bar", -13122, 0, "Unknown", 100, 0, 1.4051, 0.05, 100, 100, 0, 100, 1);

   // Creating sigma(1915)-_bar
   new Particle("sigma(1915)-_bar", -13116, 0, "Unknown", 100, 1, 1.915, 0.12, 100, 100, 0, 100, 1);

   // Creating sigma(1670)-_bar
   new Particle("sigma(1670)-_bar", -13114, 0, "Unknown", 100, 1, 1.67, 0.06, 100, 100, 0, 100, 1);

   // Creating sigma(1660)-_bar
   new Particle("sigma(1660)-_bar", -13112, 0, "Unknown", 100, 1, 1.66, 0.1, 100, 100, 0, 100, 1);

   // Creating delta(1930)++_bar
   new Particle("delta(1930)++_bar", -12226, 0, "Unknown", 100, -2, 1.96, 0.36, 100, 100, 0, 100, 1);

   // Creating delta(1700)++_bar
   new Particle("delta(1700)++_bar", -12224, 0, "Unknown", 100, -2, 1.7, 0.3, 100, 100, 0, 100, 1);

   // Creating delta(1900)++_bar
   new Particle("delta(1900)++_bar", -12222, 0, "Unknown", 100, -2, 1.9, 0.2, 100, 100, 0, 100, 1);

   // Creating N(1990)+_bar
   new Particle("N(1990)+_bar", -12218, 0, "Unknown", 100, -1, 1.95, 0.555, 100, 100, 0, 100, 1);

   // Creating N(1680)+_bar
   new Particle("N(1680)+_bar", -12216, 0, "Unknown", 100, -1, 1.685, 0.13, 100, 100, 0, 100, 1);

   // Creating delta(1700)+_bar
   new Particle("delta(1700)+_bar", -12214, 0, "Unknown", 100, -1, 1.7, 0.3, 100, 100, 0, 100, 1);

   // Creating N(1440)+_bar
   new Particle("N(1440)+_bar", -12212, 0, "Unknown", 100, -1, 1.44, 0.3, 100, 100, 0, 100, 1);

   // Creating delta(1930)+_bar
   new Particle("delta(1930)+_bar", -12126, 0, "Unknown", 100, -1, 1.96, 0.36, 100, 100, 0, 100, 1);

   // Creating delta(1900)+_bar
   new Particle("delta(1900)+_bar", -12122, 0, "Unknown", 100, -1, 1.9, 0.2, 100, 100, 0, 100, 1);

   // Creating N(1990)0_bar
   new Particle("N(1990)0_bar", -12118, 0, "Unknown", 100, 0, 1.95, 0.555, 100, 100, 0, 100, 1);

   // Creating N(1680)0_bar
   new Particle("N(1680)0_bar", -12116, 0, "Unknown", 100, 0, 1.685, 0.13, 100, 100, 0, 100, 1);

   // Creating delta(1700)0_bar
   new Particle("delta(1700)0_bar", -12114, 0, "Unknown", 100, 0, 1.7, 0.3, 100, 100, 0, 100, 1);

   // Creating N(1440)0_bar
   new Particle("N(1440)0_bar", -12112, 0, "Unknown", 100, 0, 1.44, 0.3, 100, 100, 0, 100, 1);

   // Creating delta(1930)0_bar
   new Particle("delta(1930)0_bar", -11216, 0, "Unknown", 100, 0, 1.96, 0.36, 100, 100, 0, 100, 1);

   // Creating delta(1900)0_bar
   new Particle("delta(1900)0_bar", -11212, 0, "Unknown", 100, 0, 1.9, 0.2, 100, 100, 0, 100, 1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
