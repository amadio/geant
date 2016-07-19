#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
using std::vector;
namespace geant {
   inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
void CreateParticle0004() {

   // Creating B*_00_bar
   new Particle("B*_00_bar", -10511, 0, "Unknown", 100, 0, 5.68, 0.05, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-10511));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-521,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-511,-111}));

   // Creating D_1s-
   new Particle("D_1s-", -10433, 0, "Unknown", 100, -1, 2.5353, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10433));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-423,-321}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-413,-311}));

   // Creating D*_0s-
   new Particle("D*_0s-", -10431, 0, "Unknown", 100, -1, 2.3178, 0.0046, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10431));
   part->AddDecay(Particle::Decay(0, 0.8,  vector<int>{-431,-111}));
   part->AddDecay(Particle::Decay(0, 0.2,  vector<int>{-431,-22}));

   // Creating D_10_bar
   new Particle("D_10_bar", -10423, 0, "Unknown", 100, 0, 2.4223, 0.02, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10423));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-413,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-423,-111}));

   // Creating D*_00_bar
   new Particle("D*_00_bar", -10421, 0, "Unknown", 100, 0, 2.272, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10421));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-411,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-421,-111}));

   // Creating D_1-
   new Particle("D_1-", -10413, 0, "Unknown", 100, -1, 2.424, 0.02, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10413));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-423,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-413,-111}));

   // Creating D*_0-
   new Particle("D*_0-", -10411, 0, "Unknown", 100, -1, 2.272, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10411));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-421,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-411,-111}));

   // Creating eta2(1870)_bar
   new Particle("eta2(1870)_bar", -10335, 0, "Unknown", 100, 0, 1.842, 0.225, 100, 100, 0, 100, 1);

   // Creating k2(1770)-_bar
   new Particle("k2(1770)-_bar", -10325, 0, "Unknown", 100, -1, 1.773, 0.186, 100, 100, 0, 100, 1);

   // Creating K_1-
   new Particle("K_1-", -10323, 0, "Unknown", 100, -1, 1.272, 0.09, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10323));
   part->AddDecay(Particle::Decay(0, 0.313,  vector<int>{-313,-211}));
   part->AddDecay(Particle::Decay(0, 0.28,  vector<int>{-311,-213}));
   part->AddDecay(Particle::Decay(0, 0.157,  vector<int>{-323,-111}));
   part->AddDecay(Particle::Decay(0, 0.14,  vector<int>{-321,-113}));
   part->AddDecay(Particle::Decay(0, 0.11,  vector<int>{-321,-223}));

   // Creating K*_0-
   new Particle("K*_0-", -10321, 0, "Unknown", 100, -1, 1.42, 0.287, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10321));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-311,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-321,-111}));

   // Creating k2(1770)0_bar
   new Particle("k2(1770)0_bar", -10315, 0, "Unknown", 100, 0, 1.773, 0.186, 100, 100, 0, 100, 1);

   // Creating K_10_bar
   new Particle("K_10_bar", -10313, 0, "Unknown", 100, 0, 1.272, 0.09, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10313));
   part->AddDecay(Particle::Decay(0, 0.313,  vector<int>{-323,211}));
   part->AddDecay(Particle::Decay(0, 0.28,  vector<int>{-321,213}));
   part->AddDecay(Particle::Decay(0, 0.157,  vector<int>{-313,-111}));
   part->AddDecay(Particle::Decay(0, 0.14,  vector<int>{-311,-113}));
   part->AddDecay(Particle::Decay(0, 0.11,  vector<int>{-311,-223}));

   // Creating K*_00_bar
   new Particle("K*_00_bar", -10311, 0, "Unknown", 100, 0, 1.42, 0.287, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10311));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-321,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-311,-111}));

   // Creating eta2(1645)_bar
   new Particle("eta2(1645)_bar", -10225, 0, "Unknown", 100, 0, 1.617, 0.181, 100, 100, 0, 100, 1);

   // Creating pi2(1670)-_bar
   new Particle("pi2(1670)-_bar", -10215, 0, "Unknown", 100, -1, 1.6722, 0.26, 100, 100, 0, 100, 1);

   // Creating b_1-
   new Particle("b_1-", -10213, 0, "Unknown", 100, -1, 1.2295, 0.142, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10213));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-223,-211}));

   // Creating a_0-
   new Particle("a_0-", -10211, 0, "Unknown", 100, -1, 0.9835, 0.06, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10211));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-221,-211}));

   // Creating pi2(1670)0_bar
   new Particle("pi2(1670)0_bar", -10115, 0, "Unknown", 100, 0, 1.6722, 0.26, 100, 100, 0, 100, 1);

   // Creating Omega*_bbb+
   new Particle("Omega*_bbb+", -5554, 0, "B-Baryon", 100, 1, 15.1106, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5554));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating Omega*_bbc0_bar
   new Particle("Omega*_bbc0_bar", -5544, 0, "B-Baryon", 100, 0, 11.7115, 0, 100, 100, 1, 100, 1);
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
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
