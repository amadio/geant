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
void CreateParticle0014() {

   // Creating Sigma*_c+
   new Particle("Sigma*_c+", 4214, 1, "CharmedBaryon", 100, 1, 2.5175, 0, -100, -1, -100, -1, -1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(4214));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{4122,111}));

   // Creating Sigma_c++
   new Particle("Sigma_c++", 4222, 1, "CharmedBaryon", 100, 2, 2.45402, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4222));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{4122,211}));

   // Creating Sigma*_c++
   new Particle("Sigma*_c++", 4224, 1, "CharmedBaryon", 100, 2, 2.5184, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4224));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{4122,211}));

   // Creating Xi_c+
   new Particle("Xi_c+", 4232, 1, "CharmedBaryon", 100, 1, 2.4679, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4232));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{2,-1,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-13,14,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-11,12,3,81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{2,-3,3,81}));

   // Creating cs_0
   new Particle("cs_0", 4301, 1, "CharmedBaryon", 100, 0.333333, 2.15432, 0, -100, -1, -100, -1, -1);

   // Creating cs_1
   new Particle("cs_1", 4303, 1, "Unknown", 100, 0.333333, 2.17967, 0, -100, -1, -100, -1, -1);

   // Creating Xi'_c0
   new Particle("Xi'_c0", 4312, 1, "CharmedBaryon", 100, 0, 2.578, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4312));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{4132,22}));

   // Creating Xi*_c0
   new Particle("Xi*_c0", 4314, 1, "CharmedBaryon", 100, 0, 2.6461, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4314));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{4132,111}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{4132,22}));

   // Creating Xi'_c+
   new Particle("Xi'_c+", 4322, 1, "CharmedBaryon", 100, 1, 2.5757, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4322));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{4232,22}));

   // Creating Xi*_c+
   new Particle("Xi*_c+", 4324, 1, "CharmedBaryon", 100, 1, 2.6466, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4324));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{4232,111}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{4232,22}));

   // Creating Omega_c0
   new Particle("Omega_c0", 4332, 1, "CharmedBaryon", 100, 0, 2.6975, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4332));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{2,-1,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-13,14,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-11,12,3,81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{2,-3,3,81}));

   // Creating Omega*_c0
   new Particle("Omega*_c0", 4334, 1, "CharmedBaryon", 100, 0, 2.7683, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4334));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{4332,22}));

   // Creating cc_1
   new Particle("cc_1", 4403, 1, "Unknown", 100, 1.33333, 3.27531, 0, -100, -1, -100, -1, -1);

   // Creating Xi_cc+
   new Particle("Xi_cc+", 4412, 1, "CharmedBaryon", 100, 1, 3.59798, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4412));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{2,-1,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-13,14,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-11,12,3,81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{2,-3,3,81}));

   // Creating Xi*_cc+
   new Particle("Xi*_cc+", 4414, 1, "CharmedBaryon", 100, 1, 3.65648, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4414));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{2,-1,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-13,14,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-11,12,3,81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{2,-3,3,81}));

   // Creating Xi_cc++
   new Particle("Xi_cc++", 4422, 1, "CharmedBaryon", 100, 2, 3.59798, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4422));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{2,-1,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-13,14,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-11,12,3,81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{2,-3,3,81}));

   // Creating Xi*_cc++
   new Particle("Xi*_cc++", 4424, 1, "CharmedBaryon", 100, 2, 3.65648, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4424));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{2,-1,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-13,14,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-11,12,3,81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{2,-3,3,81}));

   // Creating Omega_cc+
   new Particle("Omega_cc+", 4432, 1, "CharmedBaryon", 100, 1, 3.78663, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4432));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{2,-1,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-13,14,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-11,12,3,81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{2,-3,3,81}));

   // Creating Omega*_cc+
   new Particle("Omega*_cc+", 4434, 1, "CharmedBaryon", 100, 1, 3.82466, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4434));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{2,-1,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-13,14,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-11,12,3,81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{2,-3,3,81}));

   // Creating Omega*_ccc++
   new Particle("Omega*_ccc++", 4444, 1, "CharmedBaryon", 100, 2, 4.91594, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4444));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{2,-1,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-13,14,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-11,12,3,81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{2,-3,3,81}));

   // Creating bd_0
   new Particle("bd_0", 5101, 1, "Unknown", 100, -0.666667, 5.38897, 0, -100, -1, -100, -1, -1);

   // Creating bd_1
   new Particle("bd_1", 5103, 1, "Unknown", 100, -0.666667, 5.40145, 0, -100, -1, -100, -1, -1);

   // Creating Sigma_b-
   new Particle("Sigma_b-", 5112, 1, "B-Baryon", 100, -1, 5.8152, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5112));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{5122,-211}));

   // Creating Sigma*_b-
   new Particle("Sigma*_b-", 5114, 1, "B-Baryon", 100, -1, 5.8364, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5114));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{5122,-211}));

   // Creating Lambda_b0
   new Particle("Lambda_b0", 5122, 1, "B-Baryon", 100, 0, 5.6202, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5122));
   part->AddDecay(Particle::Decay(48, 0.4291,  vector<int>{-2,1,4,2101}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4122}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4122}));
   part->AddDecay(Particle::Decay(13, 0.08,  vector<int>{-2,4,1,2101}));
   part->AddDecay(Particle::Decay(13, 0.07,  vector<int>{-4,3,4,2101}));
   part->AddDecay(Particle::Decay(0, 0.0435,  vector<int>{4122,-433}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4122}));
   part->AddDecay(Particle::Decay(0, 0.0285,  vector<int>{4122,-431}));
   part->AddDecay(Particle::Decay(0, 0.0235,  vector<int>{4122,-20213}));
   part->AddDecay(Particle::Decay(0, 0.02,  vector<int>{4122,-213}));
   part->AddDecay(Particle::Decay(13, 0.02,  vector<int>{-4,4,3,2101}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,2101}));
   part->AddDecay(Particle::Decay(0, 0.0077,  vector<int>{4122,-211}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,2101}));
   part->AddDecay(Particle::Decay(0, 0.0044,  vector<int>{20443,3122}));
   part->AddDecay(Particle::Decay(0, 0.0022,  vector<int>{443,3122}));
   part->AddDecay(Particle::Decay(0, 0.0011,  vector<int>{441,3122}));

   // Creating Xi_b-
   new Particle("Xi_b-", 5132, 1, "B-Baryon", 100, -1, 5.7924, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5132));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Xi_bc0
   new Particle("Xi_bc0", 5142, 1, "B-Baryon", 100, 0, 7.00575, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5142));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating bu_0
   new Particle("bu_0", 5201, 1, "Unknown", 100, 0.333333, 5.38897, 0, -100, -1, -100, -1, -1);

   // Creating bu_1
   new Particle("bu_1", 5203, 1, "Unknown", 100, 0.333333, 5.40145, 0, -100, -1, -100, -1, -1);

   // Creating Sigma_b0
   new Particle("Sigma_b0", 5212, 1, "B-Baryon", 100, 0, 5.8078, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5212));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{5122,111}));

   // Creating Sigma*_b0
   new Particle("Sigma*_b0", 5214, 1, "B-Baryon", 100, 0, 5.829, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5214));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{5122,111}));

   // Creating Sigma_b+
   new Particle("Sigma_b+", 5222, 1, "B-Baryon", 100, 1, 5.8078, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5222));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{5122,211}));

   // Creating Sigma*_b+
   new Particle("Sigma*_b+", 5224, 1, "B-Baryon", 100, 1, 5.829, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5224));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{5122,211}));

   // Creating Xi_b0
   new Particle("Xi_b0", 5232, 1, "B-Baryon", 100, 0, 5.7924, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5232));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating Xi_bc+
   new Particle("Xi_bc+", 5242, 1, "B-Baryon", 100, 1, 7.00575, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5242));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating bs_0
   new Particle("bs_0", 5301, 1, "Unknown", 100, -0.666667, 5.56725, 0, -100, -1, -100, -1, -1);

   // Creating bs_1
   new Particle("bs_1", 5303, 1, "Unknown", 100, -0.666667, 5.57536, 0, -100, -1, -100, -1, -1);

   // Creating Xi'_b-
   new Particle("Xi'_b-", 5312, 1, "B-Baryon", 100, -1, 5.96, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5312));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{5132,22}));

   // Creating Xi*_b-
   new Particle("Xi*_b-", 5314, 1, "B-Baryon", 100, -1, 5.97, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5314));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{5132,22}));

   // Creating Xi'_b0
   new Particle("Xi'_b0", 5322, 1, "B-Baryon", 100, 0, 5.96, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5322));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{5232,22}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
