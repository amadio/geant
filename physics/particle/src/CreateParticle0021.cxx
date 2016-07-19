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
void CreateParticle0021() {

   // Creating f_0
   new Particle("f_0", 10221, 0, "Unknown", 100, 0, 1.3, 0, -100, -1, -100, -1, -1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(10221));
   part->AddDecay(Particle::Decay(0, 0.52,  vector<int>{211,-211}));
   part->AddDecay(Particle::Decay(0, 0.26,  vector<int>{111,111}));
   part->AddDecay(Particle::Decay(0, 0.11,  vector<int>{321,-321}));
   part->AddDecay(Particle::Decay(0, 0.055,  vector<int>{130,130}));
   part->AddDecay(Particle::Decay(0, 0.055,  vector<int>{310,310}));

   // Creating h_1
   new Particle("h_1", 10223, 0, "Unknown", 100, 0, 1.17, 0.36, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10223));
   part->AddDecay(Particle::Decay(0, 0.334,  vector<int>{113,111}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{213,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-213,211}));

   // Creating eta2(1645)
   new Particle("eta2(1645)", 10225, 1, "Unknown", 100, 0, 1.617, 0.181, -100, 0, -100, -1, -1);

   // Creating K*_00
   new Particle("K*_00", 10311, 1, "Unknown", 100, 0, 1.42, 0.287, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10311));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{321,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{311,111}));

   // Creating K_10
   new Particle("K_10", 10313, 1, "Unknown", 100, 0, 1.272, 0.09, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10313));
   part->AddDecay(Particle::Decay(0, 0.313,  vector<int>{323,-211}));
   part->AddDecay(Particle::Decay(0, 0.28,  vector<int>{321,-213}));
   part->AddDecay(Particle::Decay(0, 0.157,  vector<int>{313,111}));
   part->AddDecay(Particle::Decay(0, 0.14,  vector<int>{311,113}));
   part->AddDecay(Particle::Decay(0, 0.11,  vector<int>{311,223}));

   // Creating k2(1770)0
   new Particle("k2(1770)0", 10315, 1, "Unknown", 100, 0, 1.773, 0.186, -100, 0, -100, -1, -1);

   // Creating K*_0+
   new Particle("K*_0+", 10321, 1, "Unknown", 100, 1, 1.42, 0.287, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10321));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{311,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{321,111}));

   // Creating K_1+
   new Particle("K_1+", 10323, 1, "Unknown", 100, 1, 1.272, 0.09, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10323));
   part->AddDecay(Particle::Decay(0, 0.313,  vector<int>{313,211}));
   part->AddDecay(Particle::Decay(0, 0.28,  vector<int>{311,213}));
   part->AddDecay(Particle::Decay(0, 0.157,  vector<int>{323,111}));
   part->AddDecay(Particle::Decay(0, 0.14,  vector<int>{321,113}));
   part->AddDecay(Particle::Decay(0, 0.11,  vector<int>{321,223}));

   // Creating k2(1770)+
   new Particle("k2(1770)+", 10325, 1, "Unknown", 100, 1, 1.773, 0.186, -100, 0, -100, -1, -1);

   // Creating f'_0
   new Particle("f'_0", 10331, 0, "Unknown", 100, 0, 1.724, 0.25, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10331));
   part->AddDecay(Particle::Decay(0, 0.36,  vector<int>{211,-211}));
   part->AddDecay(Particle::Decay(0, 0.2,  vector<int>{211,-211,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.2,  vector<int>{211,-211,111,111}));
   part->AddDecay(Particle::Decay(0, 0.18,  vector<int>{111,111}));
   part->AddDecay(Particle::Decay(0, 0.03,  vector<int>{321,-321}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{130,130}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{310,310}));

   // Creating h'_1
   new Particle("h'_1", 10333, 0, "Unknown", 100, 0, 1.4, 0.08, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10333));
   part->AddDecay(Particle::Decay(0, 0.25,  vector<int>{313,-311}));
   part->AddDecay(Particle::Decay(0, 0.25,  vector<int>{-313,311}));
   part->AddDecay(Particle::Decay(0, 0.25,  vector<int>{323,-321}));
   part->AddDecay(Particle::Decay(0, 0.25,  vector<int>{-323,321}));

   // Creating eta2(1870)
   new Particle("eta2(1870)", 10335, 1, "Unknown", 100, 0, 1.842, 0.225, -100, 0, -100, -1, -1);

   // Creating D*_0+
   new Particle("D*_0+", 10411, 1, "Unknown", 100, 1, 2.272, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10411));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{421,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{411,111}));

   // Creating D_1+
   new Particle("D_1+", 10413, 1, "Unknown", 100, 1, 2.424, 0.02, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10413));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{423,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{413,111}));

   // Creating D*_00
   new Particle("D*_00", 10421, 1, "Unknown", 100, 0, 2.272, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10421));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{411,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{421,111}));

   // Creating D_10
   new Particle("D_10", 10423, 1, "Unknown", 100, 0, 2.4223, 0.02, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10423));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{413,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{423,111}));

   // Creating D*_0s+
   new Particle("D*_0s+", 10431, 1, "Unknown", 100, 1, 2.3178, 0.0046, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10431));
   part->AddDecay(Particle::Decay(0, 0.8,  vector<int>{431,111}));
   part->AddDecay(Particle::Decay(0, 0.2,  vector<int>{431,22}));

   // Creating D_1s+
   new Particle("D_1s+", 10433, 1, "Unknown", 100, 1, 2.5353, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10433));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{423,321}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{413,311}));

   // Creating chi_0c
   new Particle("chi_0c", 10441, 0, "Unknown", 100, 0, 3.41475, 0.014, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10441));
   part->AddDecay(Particle::Decay(12, 0.993,  vector<int>{82,-82}));
   part->AddDecay(Particle::Decay(0, 0.007,  vector<int>{443,22}));

   // Creating h_1c
   new Particle("h_1c", 10443, 0, "Unknown", 100, 0, 3.52593, 0.01, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10443));
   part->AddDecay(Particle::Decay(12, 1,  vector<int>{82,-82}));

   // Creating B*_00
   new Particle("B*_00", 10511, 1, "Unknown", 100, 0, 5.68, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10511));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{521,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{511,111}));

   // Creating B_10
   new Particle("B_10", 10513, 1, "Unknown", 100, 0, 5.73, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10513));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{523,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{513,111}));

   // Creating B*_0+
   new Particle("B*_0+", 10521, 1, "Unknown", 100, 1, 5.68, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10521));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{511,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{521,111}));

   // Creating B_1+
   new Particle("B_1+", 10523, 1, "Unknown", 100, 1, 5.73, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10523));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{513,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{523,111}));

   // Creating B*_0s0
   new Particle("B*_0s0", 10531, 1, "Unknown", 100, 0, 5.92, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10531));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{521,-321}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{511,-311}));

   // Creating B_1s0
   new Particle("B_1s0", 10533, 1, "Unknown", 100, 0, 5.97, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10533));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{523,-321}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{513,-311}));

   // Creating B*_0c+
   new Particle("B*_0c+", 10541, 1, "Unknown", 100, 1, 7.25, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10541));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{511,411}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{521,421}));

   // Creating B_1c+
   new Particle("B_1c+", 10543, 1, "Unknown", 100, 1, 7.3, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10543));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{513,411}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{523,421}));

   // Creating chi_0b
   new Particle("chi_0b", 10551, 0, "Unknown", 100, 0, 9.8594, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10551));
   part->AddDecay(Particle::Decay(32, 0.98,  vector<int>{21,21}));
   part->AddDecay(Particle::Decay(0, 0.02,  vector<int>{553,22}));

   // Creating h_1b
   new Particle("h_1b", 10553, 0, "Unknown", 100, 0, 9.875, 0.01, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10553));
   part->AddDecay(Particle::Decay(32, 1,  vector<int>{21,21}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
