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
void CreateParticle0002() {

   // Creating rho(1450)0_bar
   new Particle("rho(1450)0_bar", -100113, 0, "Unknown", 100, 0, 1.465, 0.4, 100, 100, 0, 100, 1);

   // Creating pi(1300)0_bar
   new Particle("pi(1300)0_bar", -100111, 0, "Unknown", 100, 0, 1.3, 0.4, 100, 100, 0, 100, 1);

   // Creating lambda(1810)_bar
   new Particle("lambda(1810)_bar", -53122, 0, "Unknown", 100, 0, 1.81, 0.15, 100, 100, 0, 100, 1);

   // Creating N(2090)+_bar
   new Particle("N(2090)+_bar", -52214, 0, "Unknown", 100, -1, 2.08, 0.35, 100, 100, 0, 100, 1);

   // Creating N(2090)0_bar
   new Particle("N(2090)0_bar", -52114, 0, "Unknown", 100, 0, 2.08, 0.35, 100, 100, 0, 100, 1);

   // Creating lambda(1800)_bar
   new Particle("lambda(1800)_bar", -43122, 0, "Unknown", 100, 0, 1.8, 0.3, 100, 100, 0, 100, 1);

   // Creating N(1710)+_bar
   new Particle("N(1710)+_bar", -42212, 0, "Unknown", 100, -1, 1.71, 0.1, 100, 100, 0, 100, 1);

   // Creating N(1900)+_bar
   new Particle("N(1900)+_bar", -42124, 0, "Unknown", 100, -1, 1.9, 0.5, 100, 100, 0, 100, 1);

   // Creating N(1710)0_bar
   new Particle("N(1710)0_bar", -42112, 0, "Unknown", 100, 0, 1.71, 0.1, 100, 100, 0, 100, 1);

   // Creating N(1900)0_bar
   new Particle("N(1900)0_bar", -41214, 0, "Unknown", 100, 0, 1.9, 0.5, 100, 100, 0, 100, 1);

   // Creating xi(1950)0_bar
   new Particle("xi(1950)0_bar", -33324, 0, "Unknown", 100, 0, 1.95, 0.06, 100, 100, 0, 100, 1);

   // Creating xi(1950)-_bar
   new Particle("xi(1950)-_bar", -33314, 0, "Unknown", 100, 1, 1.95, 0.06, 100, 100, 0, 100, 1);

   // Creating lambda(1670)_bar
   new Particle("lambda(1670)_bar", -33122, 0, "Unknown", 100, 0, 1.67, 0.035, 100, 100, 0, 100, 1);

   // Creating delta(1600)++_bar
   new Particle("delta(1600)++_bar", -32224, 0, "Unknown", 100, -2, 1.6, 0.35, 100, 100, 0, 100, 1);

   // Creating delta(1600)+_bar
   new Particle("delta(1600)+_bar", -32214, 0, "Unknown", 100, -1, 1.6, 0.35, 100, 100, 0, 100, 1);

   // Creating N(1650)+_bar
   new Particle("N(1650)+_bar", -32212, 0, "Unknown", 100, -1, 1.655, 0.165, 100, 100, 0, 100, 1);

   // Creating N(1720)+_bar
   new Particle("N(1720)+_bar", -32124, 0, "Unknown", 100, -1, 1.72, 0.2, 100, 100, 0, 100, 1);

   // Creating delta(1600)0_bar
   new Particle("delta(1600)0_bar", -32114, 0, "Unknown", 100, 0, 1.6, 0.35, 100, 100, 0, 100, 1);

   // Creating N(1650)0_bar
   new Particle("N(1650)0_bar", -32112, 0, "Unknown", 100, 0, 1.655, 0.165, 100, 100, 0, 100, 1);

   // Creating N(1720)0_bar
   new Particle("N(1720)0_bar", -31214, 0, "Unknown", 100, 0, 1.72, 0.2, 100, 100, 0, 100, 1);

   // Creating delta(1600)-_bar
   new Particle("delta(1600)-_bar", -31114, 0, "Unknown", 100, 1, 1.6, 0.35, 100, 100, 0, 100, 1);

   // Creating k_star(1680)-_bar
   new Particle("k_star(1680)-_bar", -30323, 0, "Unknown", 100, -1, 1.717, 0.32, 100, 100, 0, 100, 1);

   // Creating k_star(1680)0_bar
   new Particle("k_star(1680)0_bar", -30313, 0, "Unknown", 100, 0, 1.717, 0.32, 100, 100, 0, 100, 1);

   // Creating omega(1650)_bar
   new Particle("omega(1650)_bar", -30223, 0, "Unknown", 100, 0, 1.67, 0.315, 100, 100, 0, 100, 1);

   // Creating rho(1700)-_bar
   new Particle("rho(1700)-_bar", -30213, 0, "Unknown", 100, -1, 1.72, 0.25, 100, 100, 0, 100, 1);

   // Creating rho(1700)0_bar
   new Particle("rho(1700)0_bar", -30113, 0, "Unknown", 100, 0, 1.72, 0.25, 100, 100, 0, 100, 1);

   // Creating xi(1690)0_bar
   new Particle("xi(1690)0_bar", -23324, 0, "Unknown", 100, 0, 1.69, 0.05, 100, 100, 0, 100, 1);

   // Creating xi(1690)-_bar
   new Particle("xi(1690)-_bar", -23314, 0, "Unknown", 100, 1, 1.69, 0.05, 100, 100, 0, 100, 1);

   // Creating sigma(1940)+_bar
   new Particle("sigma(1940)+_bar", -23224, 0, "Unknown", 100, -1, 1.94, 0.22, 100, 100, 0, 100, 1);

   // Creating sigma(1750)+_bar
   new Particle("sigma(1750)+_bar", -23222, 0, "Unknown", 100, -1, 1.75, 0.09, 100, 100, 0, 100, 1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
