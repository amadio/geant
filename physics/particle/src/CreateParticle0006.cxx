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
void CreateParticle0006() {

   // Creating Omega*_bbc0_bar
   new Particle("Omega*_bbc0_bar", -5544, 0, "B-Baryon", 100, 0, 11.7115, 0, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-5544));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Omega_bbc0_bar
   new Particle("Omega_bbc0_bar", -5542, 0, "B-Baryon", 100, 0, 11.7077, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5542));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Omega*_bb+
   new Particle("Omega*_bb+", -5534, 0, "B-Baryon", 100, 1, 10.6143, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5534));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Omega_bb+
   new Particle("Omega_bb+", -5532, 0, "B-Baryon", 100, 1, 10.6021, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5532));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Xi*_bb0_bar
   new Particle("Xi*_bb0_bar", -5524, 0, "B-Baryon", 100, 0, 10.4414, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5524));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Xi_bb0_bar
   new Particle("Xi_bb0_bar", -5522, 0, "B-Baryon", 100, 0, 10.4227, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5522));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Xi*_bb+
   new Particle("Xi*_bb+", -5514, 0, "B-Baryon", 100, 1, 10.4414, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5514));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Xi_bb+
   new Particle("Xi_bb+", -5512, 0, "Unknown", 100, 1, 10.4227, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5512));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating bb_1_bar
   new Particle("bb_1_bar", -5503, 0, "Unknown", 100, 0.666667, 10.0735, 0, 100, 100, 1, 100, 1);

   // Creating Omega*_bcc-
   new Particle("Omega*_bcc-", -5444, 0, "B-Baryon", 100, -1, 8.31325, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5444));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Omega_bcc+_bar
   new Particle("Omega_bcc+_bar", -5442, 0, "B-Baryon", 100, -1, 8.30945, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5442));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Omega*_bc0_bar
   new Particle("Omega*_bc0_bar", -5434, 0, "B-Baryon", 100, 0, 7.219, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5434));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Omega'_bc0_bar
   new Particle("Omega'_bc0_bar", -5432, 0, "B-Baryon", 100, 0, 7.21101, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5432));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Xi*_bc-
   new Particle("Xi*_bc-", -5424, 0, "B-Baryon", 100, -1, 7.0485, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5424));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Xi'_bc-
   new Particle("Xi'_bc-", -5422, 0, "B-Baryon", 100, -1, 7.03724, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5422));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Xi*_bc0_bar
   new Particle("Xi*_bc0_bar", -5414, 0, "B-Baryon", 100, 0, 7.0485, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5414));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Xi'_bc0_bar
   new Particle("Xi'_bc0_bar", -5412, 0, "B-Baryon", 100, 0, 7.03724, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5412));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating bc_1_bar
   new Particle("bc_1_bar", -5403, 0, "Unknown", 100, -0.333333, 6.67397, 0, 100, 100, 1, 100, 1);

   // Creating bc_0_bar
   new Particle("bc_0_bar", -5401, 0, "Unknown", 100, -0.333333, 6.67143, 0, 100, 100, 1, 100, 1);

   // Creating Omega_bc0_bar
   new Particle("Omega_bc0_bar", -5342, 0, "B-Baryon", 100, 0, 7.19099, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5342));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Omega*_b+
   new Particle("Omega*_b+", -5334, 0, "B-Baryon", 100, 1, 6.13, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5334));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-5332,-22}));

   // Creating Omega_b+
   new Particle("Omega_b+", -5332, 0, "B-Baryon", 100, 1, 6.12, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5332));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Xi*_b0_bar
   new Particle("Xi*_b0_bar", -5324, 0, "B-Baryon", 100, 0, 5.97, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5324));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-5232,-22}));

   // Creating Xi'_b0_bar
   new Particle("Xi'_b0_bar", -5322, 0, "B-Baryon", 100, 0, 5.96, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5322));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-5232,-22}));

   // Creating Xi*_b+
   new Particle("Xi*_b+", -5314, 0, "B-Baryon", 100, 1, 5.97, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5314));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-5132,-22}));

   // Creating Xi'_b+
   new Particle("Xi'_b+", -5312, 0, "B-Baryon", 100, 1, 5.96, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5312));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-5132,-22}));

   // Creating bs_1_bar
   new Particle("bs_1_bar", -5303, 0, "Unknown", 100, 0.666667, 5.57536, 0, 100, 100, 1, 100, 1);

   // Creating bs_0_bar
   new Particle("bs_0_bar", -5301, 0, "Unknown", 100, 0.666667, 5.56725, 0, 100, 100, 1, 100, 1);

   // Creating Xi_bc-
   new Particle("Xi_bc-", -5242, 0, "B-Baryon", 100, -1, 7.00575, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5242));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Xi_b0_bar
   new Particle("Xi_b0_bar", -5232, 0, "B-Baryon", 100, 0, 5.7924, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5232));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
