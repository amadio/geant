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
void CreateParticle0007() {

   // Creating Sigma*_b+_bar
   new Particle("Sigma*_b+_bar", -5224, 0, "B-Baryon", 100, -1, 5.829, 0, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-5224));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-5122,-211}));

   // Creating Sigma_b+_bar
   new Particle("Sigma_b+_bar", -5222, 0, "B-Baryon", 100, -1, 5.8078, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5222));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-5122,-211}));

   // Creating Sigma*_b0_bar
   new Particle("Sigma*_b0_bar", -5214, 0, "B-Baryon", 100, 0, 5.829, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5214));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-5122,-111}));

   // Creating Sigma_b0_bar
   new Particle("Sigma_b0_bar", -5212, 0, "B-Baryon", 100, 0, 5.8078, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5212));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-5122,-111}));

   // Creating bu_1_bar
   new Particle("bu_1_bar", -5203, 0, "Unknown", 100, -0.333333, 5.40145, 0, 100, 100, 1, 100, 1);

   // Creating bu_0_bar
   new Particle("bu_0_bar", -5201, 0, "Unknown", 100, -0.333333, 5.38897, 0, 100, 100, 1, 100, 1);

   // Creating Xi_bc0_bar
   new Particle("Xi_bc0_bar", -5142, 0, "B-Baryon", 100, 0, 7.00575, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5142));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Xi_b+
   new Particle("Xi_b+", -5132, 0, "B-Baryon", 100, 1, 5.7924, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5132));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Lambda_b0_bar
   new Particle("Lambda_b0_bar", -5122, 0, "B-Baryon", 100, 0, 5.6202, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5122));
   part->AddDecay(Particle::Decay(48, 0.4291,  vector<int>{2,-1,-4,-2101}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4122}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4122}));
   part->AddDecay(Particle::Decay(13, 0.08,  vector<int>{2,-4,-1,-2101}));
   part->AddDecay(Particle::Decay(13, 0.07,  vector<int>{4,-3,-4,-2101}));
   part->AddDecay(Particle::Decay(0, 0.0435,  vector<int>{-4122,433}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4122}));
   part->AddDecay(Particle::Decay(0, 0.0285,  vector<int>{-4122,431}));
   part->AddDecay(Particle::Decay(0, 0.0235,  vector<int>{-4122,20213}));
   part->AddDecay(Particle::Decay(0, 0.02,  vector<int>{-4122,213}));
   part->AddDecay(Particle::Decay(13, 0.02,  vector<int>{4,-4,-3,-2101}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-2101}));
   part->AddDecay(Particle::Decay(0, 0.0077,  vector<int>{-4122,211}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-2101}));
   part->AddDecay(Particle::Decay(0, 0.0044,  vector<int>{-20443,-3122}));
   part->AddDecay(Particle::Decay(0, 0.0022,  vector<int>{-443,-3122}));
   part->AddDecay(Particle::Decay(0, 0.0011,  vector<int>{-441,-3122}));

   // Creating Sigma*_b-_bar
   new Particle("Sigma*_b-_bar", -5114, 0, "B-Baryon", 100, 1, 5.8364, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5114));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-5122,211}));

   // Creating Sigma_b-_bar
   new Particle("Sigma_b-_bar", -5112, 0, "B-Baryon", 100, 1, 5.8152, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5112));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-5122,211}));

   // Creating bd_1_bar
   new Particle("bd_1_bar", -5103, 0, "Unknown", 100, 0.666667, 5.40145, 0, 100, 100, 1, 100, 1);

   // Creating bd_0_bar
   new Particle("bd_0_bar", -5101, 0, "Unknown", 100, 0.666667, 5.38897, 0, 100, 100, 1, 100, 1);

   // Creating Omega*_ccc--
   new Particle("Omega*_ccc--", -4444, 0, "CharmedBaryon", 100, -2, 4.91594, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4444));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{-2,1,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{13,-14,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{11,-12,-3,-81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{-2,3,-3,-81}));

   // Creating Omega*_cc-
   new Particle("Omega*_cc-", -4434, 0, "CharmedBaryon", 100, -1, 3.82466, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4434));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{-2,1,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{13,-14,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{11,-12,-3,-81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{-2,3,-3,-81}));

   // Creating Omega_cc-
   new Particle("Omega_cc-", -4432, 0, "CharmedBaryon", 100, -1, 3.78663, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4432));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{-2,1,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{13,-14,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{11,-12,-3,-81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{-2,3,-3,-81}));

   // Creating Xi*_cc--
   new Particle("Xi*_cc--", -4424, 0, "CharmedBaryon", 100, -2, 3.65648, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4424));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{-2,1,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{13,-14,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{11,-12,-3,-81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{-2,3,-3,-81}));

   // Creating Xi_cc--
   new Particle("Xi_cc--", -4422, 0, "CharmedBaryon", 100, -2, 3.59798, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4422));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{-2,1,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{13,-14,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{11,-12,-3,-81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{-2,3,-3,-81}));

   // Creating Xi*_cc-
   new Particle("Xi*_cc-", -4414, 0, "CharmedBaryon", 100, -1, 3.65648, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4414));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{-2,1,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{13,-14,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{11,-12,-3,-81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{-2,3,-3,-81}));

   // Creating Xi_cc-
   new Particle("Xi_cc-", -4412, 0, "CharmedBaryon", 100, -1, 3.59798, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4412));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{-2,1,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{13,-14,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{11,-12,-3,-81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{-2,3,-3,-81}));

   // Creating cc_1_bar
   new Particle("cc_1_bar", -4403, 0, "Unknown", 100, -1.33333, 3.27531, 0, 100, 100, 1, 100, 1);

   // Creating Omega*_c0_bar
   new Particle("Omega*_c0_bar", -4334, 0, "CharmedBaryon", 100, 0, 2.7683, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4334));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-4332,-22}));

   // Creating Omega_c0_bar
   new Particle("Omega_c0_bar", -4332, 0, "CharmedBaryon", 100, 0, 2.6975, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4332));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{-2,1,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{13,-14,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{11,-12,-3,-81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{-2,3,-3,-81}));

   // Creating Xi*_c-
   new Particle("Xi*_c-", -4324, 0, "CharmedBaryon", 100, -1, 2.6466, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4324));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-4232,-111}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-4232,-22}));

   // Creating Xi'_c-
   new Particle("Xi'_c-", -4322, 0, "CharmedBaryon", 100, -1, 2.5757, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4322));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-4232,-22}));

   // Creating Xi*_c0_bar
   new Particle("Xi*_c0_bar", -4314, 0, "CharmedBaryon", 100, 0, 2.6461, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4314));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-4132,-111}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-4132,-22}));

   // Creating Xi'_c0_bar
   new Particle("Xi'_c0_bar", -4312, 0, "CharmedBaryon", 100, 0, 2.578, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4312));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-4132,-22}));

   // Creating cs_1_bar
   new Particle("cs_1_bar", -4303, 0, "Unknown", 100, -0.333333, 2.17967, 0, 100, 100, 1, 100, 1);

   // Creating cs_0_bar
   new Particle("cs_0_bar", -4301, 0, "CharmedBaryon", 100, -0.333333, 2.15432, 0, 100, 100, 1, 100, 1);

   // Creating Xi_c-
   new Particle("Xi_c-", -4232, 0, "CharmedBaryon", 100, -1, 2.4679, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4232));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{-2,1,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{13,-14,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{11,-12,-3,-81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{-2,3,-3,-81}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
