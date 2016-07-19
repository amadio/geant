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
void CreateParticle0020() {

   // Creating Xi*_b0
   new Particle("Xi*_b0", 5324, 1, "B-Baryon", 100, 0, 5.97, 0, -100, -1, -100, -1, -1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(5324));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{5232,22}));

   // Creating Omega_b-
   new Particle("Omega_b-", 5332, 1, "B-Baryon", 100, -1, 6.12, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5332));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Omega*_b-
   new Particle("Omega*_b-", 5334, 1, "B-Baryon", 100, -1, 6.13, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5334));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{5332,22}));

   // Creating Omega_bc0
   new Particle("Omega_bc0", 5342, 1, "B-Baryon", 100, 0, 7.19099, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5342));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating bc_0
   new Particle("bc_0", 5401, 1, "Unknown", 100, 0.333333, 6.67143, 0, -100, -1, -100, -1, -1);

   // Creating bc_1
   new Particle("bc_1", 5403, 1, "Unknown", 100, 0.333333, 6.67397, 0, -100, -1, -100, -1, -1);

   // Creating Xi'_bc0
   new Particle("Xi'_bc0", 5412, 1, "B-Baryon", 100, 0, 7.03724, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5412));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Xi*_bc0
   new Particle("Xi*_bc0", 5414, 1, "B-Baryon", 100, 0, 7.0485, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5414));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Xi'_bc+
   new Particle("Xi'_bc+", 5422, 1, "B-Baryon", 100, 1, 7.03724, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5422));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Xi*_bc+
   new Particle("Xi*_bc+", 5424, 1, "B-Baryon", 100, 1, 7.0485, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5424));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Omega'_bc0
   new Particle("Omega'_bc0", 5432, 1, "B-Baryon", 100, 0, 7.21101, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5432));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Omega*_bc0
   new Particle("Omega*_bc0", 5434, 1, "B-Baryon", 100, 0, 7.219, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5434));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Omega_bcc+
   new Particle("Omega_bcc+", 5442, 1, "B-Baryon", 100, 1, 8.30945, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5442));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Omega*_bcc+
   new Particle("Omega*_bcc+", 5444, 1, "B-Baryon", 100, 1, 8.31325, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5444));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating bb_1
   new Particle("bb_1", 5503, 1, "Unknown", 100, -0.666667, 10.0735, 0, -100, -1, -100, -1, -1);

   // Creating Xi_bb-
   new Particle("Xi_bb-", 5512, 1, "Unknown", 100, -1, 10.4227, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5512));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Xi*_bb-
   new Particle("Xi*_bb-", 5514, 1, "B-Baryon", 100, -1, 10.4414, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5514));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Xi_bb0
   new Particle("Xi_bb0", 5522, 1, "B-Baryon", 100, 0, 10.4227, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5522));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Xi*_bb0
   new Particle("Xi*_bb0", 5524, 1, "B-Baryon", 100, 0, 10.4414, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5524));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Omega_bb-
   new Particle("Omega_bb-", 5532, 1, "B-Baryon", 100, -1, 10.6021, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5532));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Omega*_bb-
   new Particle("Omega*_bb-", 5534, 1, "B-Baryon", 100, -1, 10.6143, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5534));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Omega_bbc0
   new Particle("Omega_bbc0", 5542, 1, "B-Baryon", 100, 0, 11.7077, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5542));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Omega*_bbc0
   new Particle("Omega*_bbc0", 5544, 1, "B-Baryon", 100, 0, 11.7115, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5544));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Omega*_bbb-
   new Particle("Omega*_bbb-", 5554, 1, "B-Baryon", 100, -1, 15.1106, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5554));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating a_00
   new Particle("a_00", 10111, 0, "Unknown", 100, 0, 0.9835, 0.06, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10111));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{221,111}));

   // Creating b_10
   new Particle("b_10", 10113, 0, "Unknown", 100, 0, 1.2295, 0.142, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10113));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{223,111}));

   // Creating pi2(1670)0
   new Particle("pi2(1670)0", 10115, 1, "Unknown", 100, 0, 1.6722, 0.26, -100, 0, -100, -1, -1);

   // Creating a_0+
   new Particle("a_0+", 10211, 1, "Unknown", 100, 1, 0.9835, 0.06, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10211));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{221,211}));

   // Creating b_1+
   new Particle("b_1+", 10213, 1, "Unknown", 100, 1, 1.2295, 0.142, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10213));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{223,211}));

   // Creating pi2(1670)+
   new Particle("pi2(1670)+", 10215, 1, "Unknown", 100, 1, 1.6722, 0.26, -100, 0, -100, -1, -1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
