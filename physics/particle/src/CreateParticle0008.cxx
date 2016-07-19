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
void CreateParticle0008() {

   // Creating Sigma*_c--
   new Particle("Sigma*_c--", -4224, 0, "CharmedBaryon", 100, -2, 2.5184, 0, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-4224));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-4122,-211}));

   // Creating Sigma_c--
   new Particle("Sigma_c--", -4222, 0, "CharmedBaryon", 100, -2, 2.45402, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4222));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-4122,-211}));

   // Creating Sigma*_c-
   new Particle("Sigma*_c-", -4214, 0, "CharmedBaryon", 100, -1, 2.5175, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4214));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-4122,-111}));

   // Creating Sigma_c-
   new Particle("Sigma_c-", -4212, 0, "CharmedBaryon", 100, -1, 2.4529, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4212));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-4122,-111}));

   // Creating cu_1_bar
   new Particle("cu_1_bar", -4203, 0, "Unknown", 100, -1.33333, 2.00808, 0, 100, 100, 1, 100, 1);

   // Creating cu_0_bar
   new Particle("cu_0_bar", -4201, 0, "Unknown", 100, -1.33333, 1.96908, 0, 100, 100, 1, 100, 1);

   // Creating Xi_c0_bar
   new Particle("Xi_c0_bar", -4132, 0, "CharmedBaryon", 100, 0, 2.471, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4132));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{-2,1,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{13,-14,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{11,-12,-3,-81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{-2,3,-3,-81}));

   // Creating Lambda_c-
   new Particle("Lambda_c-", -4122, 0, "CharmedBaryon", 100, -1, 2.28646, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4122));
   part->AddDecay(Particle::Decay(13, 0.2432,  vector<int>{-2,1,-3,-2101}));
   part->AddDecay(Particle::Decay(13, 0.15,  vector<int>{-3,-2203}));
   part->AddDecay(Particle::Decay(13, 0.075,  vector<int>{-2,-3201}));
   part->AddDecay(Particle::Decay(13, 0.075,  vector<int>{-2,-3203}));
   part->AddDecay(Particle::Decay(13, 0.057,  vector<int>{-2,1,-3,-2103}));
   part->AddDecay(Particle::Decay(13, 0.035,  vector<int>{-2,1,-1,-2101}));
   part->AddDecay(Particle::Decay(13, 0.035,  vector<int>{-2,3,-3,-2101}));
   part->AddDecay(Particle::Decay(13, 0.03,  vector<int>{-1,-2203}));
   part->AddDecay(Particle::Decay(0, 0.025,  vector<int>{-2224,323}));
   part->AddDecay(Particle::Decay(42, 0.018,  vector<int>{13,-14,-3122}));
   part->AddDecay(Particle::Decay(42, 0.018,  vector<int>{11,-12,-3122}));
   part->AddDecay(Particle::Decay(0, 0.016,  vector<int>{-2212,311}));
   part->AddDecay(Particle::Decay(13, 0.015,  vector<int>{-2,-2101}));
   part->AddDecay(Particle::Decay(13, 0.015,  vector<int>{-2,-2103}));
   part->AddDecay(Particle::Decay(0, 0.0088,  vector<int>{-2212,313}));
   part->AddDecay(Particle::Decay(0, 0.0066,  vector<int>{-2224,321}));
   part->AddDecay(Particle::Decay(42, 0.006,  vector<int>{13,-14,-2212,211}));
   part->AddDecay(Particle::Decay(42, 0.006,  vector<int>{13,-14,-2112,-111}));
   part->AddDecay(Particle::Decay(42, 0.006,  vector<int>{11,-12,-2112,-111}));
   part->AddDecay(Particle::Decay(42, 0.006,  vector<int>{11,-12,-2212,211}));
   part->AddDecay(Particle::Decay(0, 0.0058,  vector<int>{-3122,-211}));
   part->AddDecay(Particle::Decay(0, 0.0055,  vector<int>{-3212,-211}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{11,-12,-3212}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-2214,311}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-2214,313}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{13,-14,-3212}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-3122,-213}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{13,-14,-3214}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-3122,-321}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-3122,-323}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{11,-12,-3214}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{-3212,-213}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{-3214,-211}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{-3214,-213}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{-3222,-111}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{-3222,-113}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{-3222,-223}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{-3224,-113}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{-3224,-223}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{-2112,-211}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{-2112,-213}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{-2114,-211}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{-2114,-213}));
   part->AddDecay(Particle::Decay(42, 0.003,  vector<int>{13,-14,-2112}));
   part->AddDecay(Particle::Decay(42, 0.003,  vector<int>{11,-12,-2112}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{-3224,-111}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-3322,-321}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-3212,-321}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-3212,-323}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-3222,-311}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-3222,-313}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-3322,-323}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-3324,-321}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-2212,-111}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-2212,-113}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-2212,-223}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{11,-12,-2114}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-3222,-221}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-3224,-221}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-3222,-331}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{13,-14,-2114}));
   part->AddDecay(Particle::Decay(0, 0.0018,  vector<int>{-2212,-10221}));
   part->AddDecay(Particle::Decay(0, 0.0013,  vector<int>{-2212,-333}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-2224,213}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-3224,-311}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-3224,-313}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-2224,211}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-2212,-221}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-2212,-331}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-2214,-111}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-2214,-221}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-2214,-331}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-2214,-113}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-3214,-321}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-3214,-323}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-2214,-223}));

   // Creating Sigma*_c0_bar
   new Particle("Sigma*_c0_bar", -4114, 0, "CharmedBaryon", 100, 0, 2.518, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4114));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-4122,211}));

   // Creating Sigma_c0_bar
   new Particle("Sigma_c0_bar", -4112, 0, "CharmedBaryon", 100, 0, 2.45376, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4112));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-4122,211}));

   // Creating cd_1_bar
   new Particle("cd_1_bar", -4103, 0, "Unknown", 100, -0.333333, 2.00808, 0, 100, 100, 1, 100, 1);

   // Creating cd_0_bar
   new Particle("cd_0_bar", -4101, 0, "Unknown", 100, -0.333333, 1.96908, 0, 100, 100, 1, 100, 1);

   // Creating Omega+
   new Particle("Omega+", -3334, 0, "Baryon", 100, 1, 1.67245, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3334));
   part->AddDecay(Particle::Decay(0, 0.676,  vector<int>{-3122,321}));
   part->AddDecay(Particle::Decay(0, 0.234,  vector<int>{-3322,211}));
   part->AddDecay(Particle::Decay(0, 0.085,  vector<int>{-3312,-111}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{12,-11,-3322}));

   // Creating Xi*0_bar
   new Particle("Xi*0_bar", -3324, 0, "Baryon", 100, 0, 1.5318, 0.0091, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3324));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-3312,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-3322,-111}));

   // Creating Xi0_bar
   new Particle("Xi0_bar", -3322, 0, "Baryon", 100, 0, 1.31486, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3322));
   part->AddDecay(Particle::Decay(0, 0.9954,  vector<int>{-3122,-111}));
   part->AddDecay(Particle::Decay(0, 0.0035,  vector<int>{-3212,-22}));
   part->AddDecay(Particle::Decay(0, 0.0011,  vector<int>{-3122,-22}));

   // Creating Xi*+
   new Particle("Xi*+", -3314, 0, "Baryon", 100, 1, 1.535, 0.0099, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3314));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-3322,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-3312,-111}));

   // Creating Xi-_bar
   new Particle("Xi-_bar", -3312, 0, "Baryon", 100, 1, 1.32171, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3312));
   part->AddDecay(Particle::Decay(0, 0.9988,  vector<int>{-3122,211}));
   part->AddDecay(Particle::Decay(0, 0.0006,  vector<int>{12,-11,-3122}));
   part->AddDecay(Particle::Decay(0, 0.0004,  vector<int>{14,-13,-3122}));
   part->AddDecay(Particle::Decay(0, 0.0001,  vector<int>{-3112,-22}));
   part->AddDecay(Particle::Decay(0, 0.0001,  vector<int>{12,-11,-3212}));

   // Creating ss_1_bar
   new Particle("ss_1_bar", -3303, 0, "Unknown", 100, 0.666667, 2.08, 0, 100, 100, 1, 100, 1);

   // Creating sigma(2030)+_bar
   new Particle("sigma(2030)+_bar", -3228, 0, "Unknown", 100, -1, 2.03, 0.18, 100, 100, 0, 100, 1);

   // Creating sigma(1775)+_bar
   new Particle("sigma(1775)+_bar", -3226, 0, "Unknown", 100, -1, 1.775, 0.12, 100, 100, 0, 100, 1);

   // Creating Sigma*+_bar
   new Particle("Sigma*+_bar", -3224, 0, "Baryon", 100, -1, 1.3828, 0.0358, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3224));
   part->AddDecay(Particle::Decay(0, 0.88,  vector<int>{-3122,-211}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{-3222,-111}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{-3212,-211}));

   // Creating Sigma+_bar
   new Particle("Sigma+_bar", -3222, 0, "Baryon", 100, -1, 1.18937, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3222));
   part->AddDecay(Particle::Decay(0, 0.516,  vector<int>{-2212,-111}));
   part->AddDecay(Particle::Decay(0, 0.483,  vector<int>{-2112,-211}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-2212,-22}));

   // Creating sigma(2030)0_bar
   new Particle("sigma(2030)0_bar", -3218, 0, "Unknown", 100, 0, 2.03, 0.18, 100, 100, 0, 100, 1);

   // Creating sigma(1775)0_bar
   new Particle("sigma(1775)0_bar", -3216, 0, "Unknown", 100, 0, 1.775, 0.12, 100, 100, 0, 100, 1);

   // Creating Sigma*0_bar
   new Particle("Sigma*0_bar", -3214, 0, "Baryon", 100, 0, 1.3837, 0.036, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3214));
   part->AddDecay(Particle::Decay(0, 0.88,  vector<int>{-3122,-111}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{-3222,211}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{-3112,-211}));

   // Creating Sigma0_bar
   new Particle("Sigma0_bar", -3212, 0, "Baryon", 100, 0, 1.19264, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3212));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-3122,-22}));

   // Creating su_1_bar
   new Particle("su_1_bar", -3203, 0, "Unknown", 100, -0.333333, 0.1064, 0, 100, 100, 1, 100, 1);

   // Creating su_0_bar
   new Particle("su_0_bar", -3201, 0, "Unknown", 100, -0.333333, 0.1064, 0, 100, 100, 1, 100, 1);

   // Creating lambda(2100)_bar
   new Particle("lambda(2100)_bar", -3128, 0, "Unknown", 100, 0, 2.1, 0.2, 100, 100, 0, 100, 1);

   // Creating lambda(1820)_bar
   new Particle("lambda(1820)_bar", -3126, 0, "Unknown", 100, 0, 1.82, 0.08, 100, 100, 0, 100, 1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
