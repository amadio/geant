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
void CreateParticle0017() {

   // Creating delta(1950)0
   new Particle("delta(1950)0", 2118, 1, "Unknown", 100, 0, 1.93, 0.28, -100, 0, -100, -1, -1);

   // Creating delta(1620)+
   new Particle("delta(1620)+", 2122, 1, "Unknown", 100, 1, 1.63, 0.145, -100, 0, -100, -1, -1);

   // Creating N(1520)+
   new Particle("N(1520)+", 2124, 1, "Unknown", 100, 1, 1.52, 0.115, -100, 0, -100, -1, -1);

   // Creating delta(1905)+
   new Particle("delta(1905)+", 2126, 1, "Unknown", 100, 1, 1.89, 0.33, -100, 0, -100, -1, -1);

   // Creating N(2190)+
   new Particle("N(2190)+", 2128, 1, "Unknown", 100, 1, 2.19, 0.5, -100, 0, -100, -1, -1);

   // Creating uu_1
   new Particle("uu_1", 2203, 1, "Unknown", 100, 1.33333, 0.0048, 0, -100, -1, -100, -1, -1);

   // Creating p_diffr+
   new Particle("p_diffr+", 2210, 1, "Unknown", 100, 1, 0, 0, -100, -1, -100, -1, -1);

   // Creating proton
   new Particle("proton", 2212, 1, "Baryon", 100, 1, 0.938272, 0, -100, -1, -100, -1, -1);

   // Creating Delta+
   new Particle("Delta+", 2214, 1, "Baryon", 100, 1, 1.232, 0.12, -100, -1, -100, -1, -1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(2214));
   part->AddDecay(Particle::Decay(0, 0.663,  vector<int>{2212,111}));
   part->AddDecay(Particle::Decay(0, 0.331,  vector<int>{2112,211}));
   part->AddDecay(Particle::Decay(0, 0.006,  vector<int>{2212,22}));

   // Creating N(1675)+
   new Particle("N(1675)+", 2216, 1, "Unknown", 100, 1, 1.675, 0.15, -100, 0, -100, -1, -1);

   // Creating delta(1950)+
   new Particle("delta(1950)+", 2218, 1, "Unknown", 100, 1, 1.93, 0.28, -100, 0, -100, -1, -1);

   // Creating delta(1620)++
   new Particle("delta(1620)++", 2222, 1, "Unknown", 100, 2, 1.63, 0.145, -100, 0, -100, -1, -1);

   // Creating Delta++
   new Particle("Delta++", 2224, 1, "Baryon", 100, 2, 1.232, 0.12, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(2224));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{2212,211}));

   // Creating delta(1905)++
   new Particle("delta(1905)++", 2226, 1, "Unknown", 100, 2, 1.89, 0.33, -100, 0, -100, -1, -1);

   // Creating delta(1950)++
   new Particle("delta(1950)++", 2228, 1, "Unknown", 100, 2, 1.93, 0.28, -100, 0, -100, -1, -1);

   // Creating sd_0
   new Particle("sd_0", 3101, 1, "Unknown", 100, -0.666667, 0.108, 0, -100, -1, -100, -1, -1);

   // Creating sd_1
   new Particle("sd_1", 3103, 1, "Unknown", 100, -0.666667, 0.1088, 0, -100, -1, -100, -1, -1);

   // Creating Sigma-
   new Particle("Sigma-", 3112, 1, "Baryon", 100, -1, 1.19744, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3112));
   part->AddDecay(Particle::Decay(0, 0.999,  vector<int>{2112,-211}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-12,11,2112}));

   // Creating Sigma*-
   new Particle("Sigma*-", 3114, 1, "Baryon", 100, -1, 1.3872, 0.0394, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3114));
   part->AddDecay(Particle::Decay(0, 0.88,  vector<int>{3122,-211}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{3212,-211}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{3112,111}));

   // Creating sigma(1775)-
   new Particle("sigma(1775)-", 3116, 1, "Unknown", 100, -1, 1.775, 0.12, -100, 0, -100, -1, -1);

   // Creating sigma(2030)-
   new Particle("sigma(2030)-", 3118, 1, "Unknown", 100, -1, 2.03, 0.18, -100, 0, -100, -1, -1);

   // Creating Lambda0
   new Particle("Lambda0", 3122, 1, "Baryon", 100, 0, 1.11568, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3122));
   part->AddDecay(Particle::Decay(0, 0.639,  vector<int>{2212,-211}));
   part->AddDecay(Particle::Decay(0, 0.358,  vector<int>{2112,111}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{2112,22}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-12,11,2212}));

   // Creating lambda(1520)
   new Particle("lambda(1520)", 3124, 1, "Unknown", 100, 0, 1.5195, 0.0156, -100, 0, -100, -1, -1);

   // Creating lambda(1820)
   new Particle("lambda(1820)", 3126, 1, "Unknown", 100, 0, 1.82, 0.08, -100, 0, -100, -1, -1);

   // Creating lambda(2100)
   new Particle("lambda(2100)", 3128, 1, "Unknown", 100, 0, 2.1, 0.2, -100, 0, -100, -1, -1);

   // Creating su_0
   new Particle("su_0", 3201, 1, "Unknown", 100, 0.333333, 0.1064, 0, -100, -1, -100, -1, -1);

   // Creating su_1
   new Particle("su_1", 3203, 1, "Unknown", 100, 0.333333, 0.1064, 0, -100, -1, -100, -1, -1);

   // Creating Sigma0
   new Particle("Sigma0", 3212, 1, "Baryon", 100, 0, 1.19264, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3212));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{3122,22}));

   // Creating Sigma*0
   new Particle("Sigma*0", 3214, 1, "Baryon", 100, 0, 1.3837, 0.036, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3214));
   part->AddDecay(Particle::Decay(0, 0.88,  vector<int>{3122,111}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{3222,-211}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{3112,211}));

   // Creating sigma(1775)0
   new Particle("sigma(1775)0", 3216, 1, "Unknown", 100, 0, 1.775, 0.12, -100, 0, -100, -1, -1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
