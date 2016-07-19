#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
using std::vector;
namespace geant {
   inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
void CreateParticle0018() {

   // Creating delta(1920)++
   new Particle("delta(1920)++", 22224, 1, "Unknown", 100, 2, 1.92, 0.2, -100, 0, -100, -1, -1);

   // Creating sigma(1750)-
   new Particle("sigma(1750)-", 23112, 1, "Unknown", 100, -1, 1.75, 0.09, -100, 0, -100, -1, -1);

   // Creating sigma(1940)-
   new Particle("sigma(1940)-", 23114, 1, "Unknown", 100, -1, 1.94, 0.22, -100, 0, -100, -1, -1);

   // Creating lambda(1600)
   new Particle("lambda(1600)", 23122, 1, "Unknown", 100, 0, 1.6, 0.15, -100, 0, -100, -1, -1);

   // Creating lambda(1890)
   new Particle("lambda(1890)", 23124, 1, "Unknown", 100, 0, 1.89, 0.1, -100, 0, -100, -1, -1);

   // Creating lambda(2110)
   new Particle("lambda(2110)", 23126, 1, "Unknown", 100, 0, 2.11, 0.2, -100, 0, -100, -1, -1);

   // Creating sigma(1750)0
   new Particle("sigma(1750)0", 23212, 1, "Unknown", 100, 0, 1.75, 0.09, -100, 0, -100, -1, -1);

   // Creating sigma(1940)0
   new Particle("sigma(1940)0", 23214, 1, "Unknown", 100, 0, 1.94, 0.22, -100, 0, -100, -1, -1);

   // Creating sigma(1750)+
   new Particle("sigma(1750)+", 23222, 1, "Unknown", 100, 1, 1.75, 0.09, -100, 0, -100, -1, -1);

   // Creating sigma(1940)+
   new Particle("sigma(1940)+", 23224, 1, "Unknown", 100, 1, 1.94, 0.22, -100, 0, -100, -1, -1);

   // Creating xi(1690)-
   new Particle("xi(1690)-", 23314, 1, "Unknown", 100, -1, 1.69, 0.05, -100, 0, -100, -1, -1);

   // Creating xi(1690)0
   new Particle("xi(1690)0", 23324, 1, "Unknown", 100, 0, 1.69, 0.05, -100, 0, -100, -1, -1);

   // Creating rho(1700)0
   new Particle("rho(1700)0", 30113, 1, "Unknown", 100, 0, 1.72, 0.25, -100, 0, -100, -1, -1);

   // Creating rho(1700)+
   new Particle("rho(1700)+", 30213, 1, "Unknown", 100, 1, 1.72, 0.25, -100, 0, -100, -1, -1);

   // Creating omega(1650)
   new Particle("omega(1650)", 30223, 1, "Unknown", 100, 0, 1.67, 0.315, -100, 0, -100, -1, -1);

   // Creating k_star(1680)0
   new Particle("k_star(1680)0", 30313, 1, "Unknown", 100, 0, 1.717, 0.32, -100, 0, -100, -1, -1);

   // Creating k_star(1680)+
   new Particle("k_star(1680)+", 30323, 1, "Unknown", 100, 1, 1.717, 0.32, -100, 0, -100, -1, -1);

   // Creating delta(1600)-
   new Particle("delta(1600)-", 31114, 1, "Unknown", 100, -1, 1.6, 0.35, -100, 0, -100, -1, -1);

   // Creating N(1720)0
   new Particle("N(1720)0", 31214, 1, "Unknown", 100, 0, 1.72, 0.2, -100, 0, -100, -1, -1);

   // Creating N(1650)0
   new Particle("N(1650)0", 32112, 1, "Unknown", 100, 0, 1.655, 0.165, -100, 0, -100, -1, -1);

   // Creating delta(1600)0
   new Particle("delta(1600)0", 32114, 1, "Unknown", 100, 0, 1.6, 0.35, -100, 0, -100, -1, -1);

   // Creating N(1720)+
   new Particle("N(1720)+", 32124, 1, "Unknown", 100, 1, 1.72, 0.2, -100, 0, -100, -1, -1);

   // Creating N(1650)+
   new Particle("N(1650)+", 32212, 1, "Unknown", 100, 1, 1.655, 0.165, -100, 0, -100, -1, -1);

   // Creating delta(1600)+
   new Particle("delta(1600)+", 32214, 1, "Unknown", 100, 1, 1.6, 0.35, -100, 0, -100, -1, -1);

   // Creating delta(1600)++
   new Particle("delta(1600)++", 32224, 1, "Unknown", 100, 2, 1.6, 0.35, -100, 0, -100, -1, -1);

   // Creating lambda(1670)
   new Particle("lambda(1670)", 33122, 1, "Unknown", 100, 0, 1.67, 0.035, -100, 0, -100, -1, -1);

   // Creating xi(1950)-
   new Particle("xi(1950)-", 33314, 1, "Unknown", 100, -1, 1.95, 0.06, -100, 0, -100, -1, -1);

   // Creating xi(1950)0
   new Particle("xi(1950)0", 33324, 1, "Unknown", 100, 0, 1.95, 0.06, -100, 0, -100, -1, -1);

   // Creating N(1900)0
   new Particle("N(1900)0", 41214, 1, "Unknown", 100, 0, 1.9, 0.5, -100, 0, -100, -1, -1);

   // Creating N(1710)0
   new Particle("N(1710)0", 42112, 1, "Unknown", 100, 0, 1.71, 0.1, -100, 0, -100, -1, -1);

   // Creating N(1900)+
   new Particle("N(1900)+", 42124, 1, "Unknown", 100, 1, 1.9, 0.5, -100, 0, -100, -1, -1);

   // Creating N(1710)+
   new Particle("N(1710)+", 42212, 1, "Unknown", 100, 1, 1.71, 0.1, -100, 0, -100, -1, -1);

   // Creating lambda(1800)
   new Particle("lambda(1800)", 43122, 1, "Unknown", 100, 0, 1.8, 0.3, -100, 0, -100, -1, -1);

   // Creating N(2090)0
   new Particle("N(2090)0", 52114, 1, "Unknown", 100, 0, 2.08, 0.35, -100, 0, -100, -1, -1);

   // Creating N(2090)+
   new Particle("N(2090)+", 52214, 1, "Unknown", 100, 1, 2.08, 0.35, -100, 0, -100, -1, -1);

   // Creating lambda(1810)
   new Particle("lambda(1810)", 53122, 1, "Unknown", 100, 0, 1.81, 0.15, -100, 0, -100, -1, -1);

   // Creating pi(1300)0
   new Particle("pi(1300)0", 100111, 1, "Unknown", 100, 0, 1.3, 0.4, -100, 0, -100, -1, -1);

   // Creating rho(1450)0
   new Particle("rho(1450)0", 100113, 1, "Unknown", 100, 0, 1.465, 0.4, -100, 0, -100, -1, -1);

   // Creating pi(1300)+
   new Particle("pi(1300)+", 100211, 1, "Unknown", 100, 1, 1.3, 0.4, -100, 0, -100, -1, -1);

   // Creating rho(1450)+
   new Particle("rho(1450)+", 100213, 1, "Unknown", 100, 1, 1.465, 0.4, -100, 0, -100, -1, -1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
