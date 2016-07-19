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
void CreateParticle0009() {

   // Creating lambda(1520)_bar
   new Particle("lambda(1520)_bar", -3124, 0, "Unknown", 100, 0, 1.5195, 0.0156, 100, 100, 0, 100, 1);

   // Creating Lambda0_bar
   new Particle("Lambda0_bar", -3122, 0, "Baryon", 100, 0, 1.11568, 0, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-3122));
   part->AddDecay(Particle::Decay(0, 0.639,  vector<int>{-2212,211}));
   part->AddDecay(Particle::Decay(0, 0.358,  vector<int>{-2112,-111}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-2112,-22}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{12,-11,-2212}));

   // Creating sigma(2030)-_bar
   new Particle("sigma(2030)-_bar", -3118, 0, "Unknown", 100, 1, 2.03, 0.18, 100, 100, 0, 100, 1);

   // Creating sigma(1775)-_bar
   new Particle("sigma(1775)-_bar", -3116, 0, "Unknown", 100, 1, 1.775, 0.12, 100, 100, 0, 100, 1);

   // Creating Sigma*-_bar
   new Particle("Sigma*-_bar", -3114, 0, "Baryon", 100, 1, 1.3872, 0.0394, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3114));
   part->AddDecay(Particle::Decay(0, 0.88,  vector<int>{-3122,211}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{-3212,211}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{-3112,-111}));

   // Creating Sigma-_bar
   new Particle("Sigma-_bar", -3112, 0, "Baryon", 100, 1, 1.19744, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3112));
   part->AddDecay(Particle::Decay(0, 0.999,  vector<int>{-2112,211}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{12,-11,-2112}));

   // Creating sd_1_bar
   new Particle("sd_1_bar", -3103, 0, "Unknown", 100, 0.666667, 0.1088, 0, 100, 100, 1, 100, 1);

   // Creating sd_0_bar
   new Particle("sd_0_bar", -3101, 0, "Unknown", 100, 0.666667, 0.108, 0, 100, 100, 1, 100, 1);

   // Creating delta(1950)++_bar
   new Particle("delta(1950)++_bar", -2228, 0, "Unknown", 100, -2, 1.93, 0.28, 100, 100, 0, 100, 1);

   // Creating delta(1905)++_bar
   new Particle("delta(1905)++_bar", -2226, 0, "Unknown", 100, -2, 1.89, 0.33, 100, 100, 0, 100, 1);

   // Creating Delta--
   new Particle("Delta--", -2224, 0, "Baryon", 100, -2, 1.232, 0.12, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2224));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-2212,-211}));

   // Creating delta(1620)++_bar
   new Particle("delta(1620)++_bar", -2222, 0, "Unknown", 100, -2, 1.63, 0.145, 100, 100, 0, 100, 1);

   // Creating delta(1950)+_bar
   new Particle("delta(1950)+_bar", -2218, 0, "Unknown", 100, -1, 1.93, 0.28, 100, 100, 0, 100, 1);

   // Creating N(1675)+_bar
   new Particle("N(1675)+_bar", -2216, 0, "Unknown", 100, -1, 1.675, 0.15, 100, 100, 0, 100, 1);

   // Creating Delta+_bar
   new Particle("Delta+_bar", -2214, 0, "Baryon", 100, -1, 1.232, 0.12, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2214));
   part->AddDecay(Particle::Decay(0, 0.663,  vector<int>{-2212,-111}));
   part->AddDecay(Particle::Decay(0, 0.331,  vector<int>{-2112,-211}));
   part->AddDecay(Particle::Decay(0, 0.006,  vector<int>{-2212,-22}));

   // Creating antiproton
   new Particle("antiproton", -2212, 0, "Baryon", 100, -1, 0.938272, 0, 100, 100, 1, 100, 1);

   // Creating p_diffr+_bar
   new Particle("p_diffr+_bar", -2210, 0, "Unknown", 100, -1, 0, 0, 100, 100, 1, 100, 1);

   // Creating uu_1_bar
   new Particle("uu_1_bar", -2203, 0, "Unknown", 100, -1.33333, 0.0048, 0, 100, 100, 1, 100, 1);

   // Creating N(2190)+_bar
   new Particle("N(2190)+_bar", -2128, 0, "Unknown", 100, -1, 2.19, 0.5, 100, 100, 0, 100, 1);

   // Creating delta(1905)+_bar
   new Particle("delta(1905)+_bar", -2126, 0, "Unknown", 100, -1, 1.89, 0.33, 100, 100, 0, 100, 1);

   // Creating N(1520)+_bar
   new Particle("N(1520)+_bar", -2124, 0, "Unknown", 100, -1, 1.52, 0.115, 100, 100, 0, 100, 1);

   // Creating delta(1620)+_bar
   new Particle("delta(1620)+_bar", -2122, 0, "Unknown", 100, -1, 1.63, 0.145, 100, 100, 0, 100, 1);

   // Creating delta(1950)0_bar
   new Particle("delta(1950)0_bar", -2118, 0, "Unknown", 100, 0, 1.93, 0.28, 100, 100, 0, 100, 1);

   // Creating N(1675)0_bar
   new Particle("N(1675)0_bar", -2116, 0, "Unknown", 100, 0, 1.675, 0.15, 100, 100, 0, 100, 1);

   // Creating Delta0_bar
   new Particle("Delta0_bar", -2114, 0, "Baryon", 100, 0, 1.232, 0.12, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2114));
   part->AddDecay(Particle::Decay(0, 0.663,  vector<int>{-2112,-111}));
   part->AddDecay(Particle::Decay(0, 0.331,  vector<int>{-2212,211}));
   part->AddDecay(Particle::Decay(0, 0.006,  vector<int>{-2112,-22}));

   // Creating antineutron
   new Particle("antineutron", -2112, 0, "Baryon", 100, 0, 0.939565, 0, 100, 100, 1, 100, 1);

   // Creating n_diffr0_bar
   new Particle("n_diffr0_bar", -2110, 0, "Unknown", 100, 0, 0, 0, 100, 100, 1, 100, 1);

   // Creating ud_1_bar
   new Particle("ud_1_bar", -2103, 0, "Unknown", 100, -0.333333, 0.0072, 0, 100, 100, 1, 100, 1);

   // Creating ud_0_bar
   new Particle("ud_0_bar", -2101, 0, "Unknown", 100, -0.333333, 0.0073, 0, 100, 100, 1, 100, 1);

   // Creating N(2190)0_bar
   new Particle("N(2190)0_bar", -1218, 0, "Unknown", 100, 0, 2.19, 0.5, 100, 100, 0, 100, 1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
