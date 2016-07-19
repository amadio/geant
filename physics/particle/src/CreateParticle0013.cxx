#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
using std::vector;
namespace geant {
   inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
void CreateParticle0013() {

   // Creating delta(1950)+
   new Particle("delta(1950)+", 2218, 1, "Unknown", 100, 1, 1.93, 0.28, -100, 0, -100, -1, -1);

   // Creating delta(1620)++
   new Particle("delta(1620)++", 2222, 1, "Unknown", 100, 2, 1.63, 0.145, -100, 0, -100, -1, -1);

   // Creating Delta++
   new Particle("Delta++", 2224, 1, "Baryon", 100, 2, 1.232, 0.12, -100, -1, -100, -1, -1);
   Particle *part = 0;
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

   // Creating sigma(2030)0
   new Particle("sigma(2030)0", 3218, 1, "Unknown", 100, 0, 2.03, 0.18, -100, 0, -100, -1, -1);

   // Creating Sigma+
   new Particle("Sigma+", 3222, 1, "Baryon", 100, 1, 1.18937, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3222));
   part->AddDecay(Particle::Decay(0, 0.516,  vector<int>{2212,111}));
   part->AddDecay(Particle::Decay(0, 0.483,  vector<int>{2112,211}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{2212,22}));

   // Creating Sigma*+
   new Particle("Sigma*+", 3224, 1, "Baryon", 100, 1, 1.3828, 0.0358, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3224));
   part->AddDecay(Particle::Decay(0, 0.88,  vector<int>{3122,211}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{3222,111}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{3212,211}));

   // Creating sigma(1775)+
   new Particle("sigma(1775)+", 3226, 1, "Unknown", 100, 1, 1.775, 0.12, -100, 0, -100, -1, -1);

   // Creating sigma(2030)+
   new Particle("sigma(2030)+", 3228, 1, "Unknown", 100, 1, 2.03, 0.18, -100, 0, -100, -1, -1);

   // Creating ss_1
   new Particle("ss_1", 3303, 1, "Unknown", 100, -0.666667, 2.08, 0, -100, -1, -100, -1, -1);

   // Creating Xi-
   new Particle("Xi-", 3312, 1, "Baryon", 100, -1, 1.32171, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3312));
   part->AddDecay(Particle::Decay(0, 0.9988,  vector<int>{3122,-211}));
   part->AddDecay(Particle::Decay(0, 0.0006,  vector<int>{-12,11,3122}));
   part->AddDecay(Particle::Decay(0, 0.0004,  vector<int>{-14,13,3122}));
   part->AddDecay(Particle::Decay(0, 0.0001,  vector<int>{3112,22}));
   part->AddDecay(Particle::Decay(0, 0.0001,  vector<int>{-12,11,3212}));

   // Creating Xi*-
   new Particle("Xi*-", 3314, 1, "Baryon", 100, -1, 1.535, 0.0099, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3314));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{3322,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{3312,111}));

   // Creating Xi0
   new Particle("Xi0", 3322, 1, "Baryon", 100, 0, 1.31486, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3322));
   part->AddDecay(Particle::Decay(0, 0.9954,  vector<int>{3122,111}));
   part->AddDecay(Particle::Decay(0, 0.0035,  vector<int>{3212,22}));
   part->AddDecay(Particle::Decay(0, 0.0011,  vector<int>{3122,22}));

   // Creating Xi*0
   new Particle("Xi*0", 3324, 1, "Baryon", 100, 0, 1.5318, 0.0091, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3324));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{3312,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{3322,111}));

   // Creating Omega-
   new Particle("Omega-", 3334, 1, "Baryon", 100, -1, 1.67245, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3334));
   part->AddDecay(Particle::Decay(0, 0.676,  vector<int>{3122,-321}));
   part->AddDecay(Particle::Decay(0, 0.234,  vector<int>{3322,-211}));
   part->AddDecay(Particle::Decay(0, 0.085,  vector<int>{3312,111}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-12,11,3322}));

   // Creating cd_0
   new Particle("cd_0", 4101, 1, "Unknown", 100, 0.333333, 1.96908, 0, -100, -1, -100, -1, -1);

   // Creating cd_1
   new Particle("cd_1", 4103, 1, "Unknown", 100, 0.333333, 2.00808, 0, -100, -1, -100, -1, -1);

   // Creating Sigma_c0
   new Particle("Sigma_c0", 4112, 1, "CharmedBaryon", 100, 0, 2.45376, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4112));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{4122,-211}));

   // Creating Sigma*_c0
   new Particle("Sigma*_c0", 4114, 1, "CharmedBaryon", 100, 0, 2.518, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4114));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{4122,-211}));

   // Creating Lambda_c+
   new Particle("Lambda_c+", 4122, 1, "CharmedBaryon", 100, 1, 2.28646, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4122));
   part->AddDecay(Particle::Decay(13, 0.2432,  vector<int>{2,-1,3,2101}));
   part->AddDecay(Particle::Decay(13, 0.15,  vector<int>{3,2203}));
   part->AddDecay(Particle::Decay(13, 0.075,  vector<int>{2,3201}));
   part->AddDecay(Particle::Decay(13, 0.075,  vector<int>{2,3203}));
   part->AddDecay(Particle::Decay(13, 0.057,  vector<int>{2,-1,3,2103}));
   part->AddDecay(Particle::Decay(13, 0.035,  vector<int>{2,-1,1,2101}));
   part->AddDecay(Particle::Decay(13, 0.035,  vector<int>{2,-3,3,2101}));
   part->AddDecay(Particle::Decay(13, 0.03,  vector<int>{1,2203}));
   part->AddDecay(Particle::Decay(0, 0.025,  vector<int>{2224,-323}));
   part->AddDecay(Particle::Decay(42, 0.018,  vector<int>{-13,14,3122}));
   part->AddDecay(Particle::Decay(42, 0.018,  vector<int>{-11,12,3122}));
   part->AddDecay(Particle::Decay(0, 0.016,  vector<int>{2212,-311}));
   part->AddDecay(Particle::Decay(13, 0.015,  vector<int>{2,2101}));
   part->AddDecay(Particle::Decay(13, 0.015,  vector<int>{2,2103}));
   part->AddDecay(Particle::Decay(0, 0.0088,  vector<int>{2212,-313}));
   part->AddDecay(Particle::Decay(0, 0.0066,  vector<int>{2224,-321}));
   part->AddDecay(Particle::Decay(42, 0.006,  vector<int>{-13,14,2212,-211}));
   part->AddDecay(Particle::Decay(42, 0.006,  vector<int>{-13,14,2112,111}));
   part->AddDecay(Particle::Decay(42, 0.006,  vector<int>{-11,12,2112,111}));
   part->AddDecay(Particle::Decay(42, 0.006,  vector<int>{-11,12,2212,-211}));
   part->AddDecay(Particle::Decay(0, 0.0058,  vector<int>{3122,211}));
   part->AddDecay(Particle::Decay(0, 0.0055,  vector<int>{3212,211}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-11,12,3212}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{2214,-311}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{2214,-313}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-13,14,3212}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{3122,213}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-13,14,3214}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{3122,321}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{3122,323}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-11,12,3214}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{3212,213}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{3214,211}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{3214,213}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{3222,111}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{3222,113}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{3222,223}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{3224,113}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{3224,223}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{2112,211}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{2112,213}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{2114,211}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{2114,213}));
   part->AddDecay(Particle::Decay(42, 0.003,  vector<int>{-13,14,2112}));
   part->AddDecay(Particle::Decay(42, 0.003,  vector<int>{-11,12,2112}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{3224,111}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{3322,321}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{3212,321}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{3212,323}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{3222,311}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{3222,313}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{3322,323}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{3324,321}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{2212,111}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{2212,113}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{2212,223}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{-11,12,2114}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{3222,221}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{3224,221}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{3222,331}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{-13,14,2114}));
   part->AddDecay(Particle::Decay(0, 0.0018,  vector<int>{2212,10221}));
   part->AddDecay(Particle::Decay(0, 0.0013,  vector<int>{2212,333}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{2224,-213}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{3224,311}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{3224,313}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{2224,-211}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{2212,221}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{2212,331}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{2214,111}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{2214,221}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{2214,331}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{2214,113}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{3214,321}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{3214,323}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{2214,223}));

   // Creating Xi_c0
   new Particle("Xi_c0", 4132, 1, "CharmedBaryon", 100, 0, 2.471, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4132));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{2,-1,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-13,14,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-11,12,3,81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{2,-3,3,81}));

   // Creating cu_0
   new Particle("cu_0", 4201, 1, "Unknown", 100, 1.33333, 1.96908, 0, -100, -1, -100, -1, -1);

   // Creating cu_1
   new Particle("cu_1", 4203, 1, "Unknown", 100, 1.33333, 2.00808, 0, -100, -1, -100, -1, -1);

   // Creating Sigma_c+
   new Particle("Sigma_c+", 4212, 1, "CharmedBaryon", 100, 1, 2.4529, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4212));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{4122,111}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
