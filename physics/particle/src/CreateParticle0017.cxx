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

   // Creating sigma(1915)-
   new Particle("sigma(1915)-", 13116, 1, "Unknown", 100, -1, 1.915, 0.12, -100, 0, -100, -1, -1);

   // Creating lambda(1405)
   new Particle("lambda(1405)", 13122, 1, "Unknown", 100, 0, 1.4051, 0.05, -100, 0, -100, -1, -1);

   // Creating lambda(1690)
   new Particle("lambda(1690)", 13124, 1, "Unknown", 100, 0, 1.69, 0.06, -100, 0, -100, -1, -1);

   // Creating lambda(1830)
   new Particle("lambda(1830)", 13126, 1, "Unknown", 100, 0, 1.83, 0.095, -100, 0, -100, -1, -1);

   // Creating sigma(1660)0
   new Particle("sigma(1660)0", 13212, 1, "Unknown", 100, 0, 1.66, 0.1, -100, 0, -100, -1, -1);

   // Creating sigma(1670)0
   new Particle("sigma(1670)0", 13214, 1, "Unknown", 100, 0, 1.67, 0.06, -100, 0, -100, -1, -1);

   // Creating sigma(1915)0
   new Particle("sigma(1915)0", 13216, 1, "Unknown", 100, 0, 1.915, 0.12, -100, 0, -100, -1, -1);

   // Creating sigma(1660)+
   new Particle("sigma(1660)+", 13222, 1, "Unknown", 100, 1, 1.66, 0.1, -100, 0, -100, -1, -1);

   // Creating sigma(1670)+
   new Particle("sigma(1670)+", 13224, 1, "Unknown", 100, 1, 1.67, 0.06, -100, 0, -100, -1, -1);

   // Creating sigma(1915)+
   new Particle("sigma(1915)+", 13226, 1, "Unknown", 100, 1, 1.915, 0.12, -100, 0, -100, -1, -1);

   // Creating xi(1820)-
   new Particle("xi(1820)-", 13314, 1, "Unknown", 100, -1, 1.823, 0.024, -100, 0, -100, -1, -1);

   // Creating xi(2030)-
   new Particle("xi(2030)-", 13316, 1, "Unknown", 100, -1, 2.025, 0.02, -100, 0, -100, -1, -1);

   // Creating xi(1820)0
   new Particle("xi(1820)0", 13324, 1, "Unknown", 100, 0, 1.823, 0.024, -100, 0, -100, -1, -1);

   // Creating xi(2030)0
   new Particle("xi(2030)0", 13326, 1, "Unknown", 100, 0, 2.025, 0.02, -100, 0, -100, -1, -1);

   // Creating a_10
   new Particle("a_10", 20113, 0, "Unknown", 100, 0, 1.23, 0.4, -100, -1, -100, -1, -1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(20113));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{213,-211}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-213,211}));

   // Creating a_1+
   new Particle("a_1+", 20213, 1, "Unknown", 100, 1, 1.23, 0.4, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20213));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{113,211}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{213,111}));

   // Creating f_1
   new Particle("f_1", 20223, 0, "Unknown", 100, 0, 1.2818, 0.025, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20223));
   part->AddDecay(Particle::Decay(0, 0.15,  vector<int>{113,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.146,  vector<int>{10111,111}));
   part->AddDecay(Particle::Decay(0, 0.146,  vector<int>{-10211,211}));
   part->AddDecay(Particle::Decay(0, 0.146,  vector<int>{10211,-211}));
   part->AddDecay(Particle::Decay(0, 0.066,  vector<int>{113,22}));
   part->AddDecay(Particle::Decay(0, 0.05,  vector<int>{213,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.05,  vector<int>{221,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.05,  vector<int>{113,111,111}));
   part->AddDecay(Particle::Decay(0, 0.05,  vector<int>{-213,211,111}));
   part->AddDecay(Particle::Decay(0, 0.05,  vector<int>{221,111,111}));
   part->AddDecay(Particle::Decay(0, 0.024,  vector<int>{321,-311,-211}));
   part->AddDecay(Particle::Decay(0, 0.024,  vector<int>{311,-311,111}));
   part->AddDecay(Particle::Decay(0, 0.024,  vector<int>{311,-321,211}));
   part->AddDecay(Particle::Decay(0, 0.024,  vector<int>{321,-321,111}));

   // Creating K*_10
   new Particle("K*_10", 20313, 1, "Unknown", 100, 0, 1.403, 0.174, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20313));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{323,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{313,111}));

   // Creating K*_1+
   new Particle("K*_1+", 20323, 1, "Unknown", 100, 1, 1.403, 0.174, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20323));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{313,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{323,111}));

   // Creating f'_1
   new Particle("f'_1", 20333, 0, "Unknown", 100, 0, 1.4264, 0.053, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20333));
   part->AddDecay(Particle::Decay(0, 0.25,  vector<int>{313,-311}));
   part->AddDecay(Particle::Decay(0, 0.25,  vector<int>{-313,311}));
   part->AddDecay(Particle::Decay(0, 0.25,  vector<int>{323,-321}));
   part->AddDecay(Particle::Decay(0, 0.25,  vector<int>{-323,321}));

   // Creating D*_1+
   new Particle("D*_1+", 20413, 1, "Unknown", 100, 1, 2.372, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20413));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{423,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{413,111}));

   // Creating D*_10
   new Particle("D*_10", 20423, 1, "Unknown", 100, 0, 2.372, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20423));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{413,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{423,111}));

   // Creating D*_1s+
   new Particle("D*_1s+", 20433, 1, "Unknown", 100, 1, 2.4596, 0.0055, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20433));
   part->AddDecay(Particle::Decay(0, 0.8,  vector<int>{433,111}));
   part->AddDecay(Particle::Decay(0, 0.2,  vector<int>{433,22}));

   // Creating chi_1c
   new Particle("chi_1c", 20443, 0, "Unknown", 100, 0, 3.51066, 0.0009, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20443));
   part->AddDecay(Particle::Decay(12, 0.727,  vector<int>{82,-82}));
   part->AddDecay(Particle::Decay(0, 0.273,  vector<int>{443,22}));

   // Creating B*_10
   new Particle("B*_10", 20513, 1, "Unknown", 100, 0, 5.78, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20513));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{523,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{513,111}));

   // Creating B*_1+
   new Particle("B*_1+", 20523, 1, "Unknown", 100, 1, 5.78, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20523));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{513,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{523,111}));

   // Creating B*_1s0
   new Particle("B*_1s0", 20533, 1, "Unknown", 100, 0, 6.02, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20533));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{523,-321}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{513,-311}));

   // Creating B*_1c+
   new Particle("B*_1c+", 20543, 1, "Unknown", 100, 1, 7.3, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20543));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{513,411}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{523,421}));

   // Creating chi_1b
   new Particle("chi_1b", 20553, 0, "Unknown", 100, 0, 9.8928, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20553));
   part->AddDecay(Particle::Decay(32, 0.65,  vector<int>{21,21}));
   part->AddDecay(Particle::Decay(0, 0.35,  vector<int>{553,22}));

   // Creating delta(1910)-
   new Particle("delta(1910)-", 21112, 1, "Unknown", 100, -1, 1.91, 0.25, -100, 0, -100, -1, -1);

   // Creating delta(1920)-
   new Particle("delta(1920)-", 21114, 1, "Unknown", 100, -1, 1.92, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1910)0
   new Particle("delta(1910)0", 21212, 1, "Unknown", 100, 0, 1.91, 0.25, -100, 0, -100, -1, -1);

   // Creating N(1700)0
   new Particle("N(1700)0", 21214, 1, "Unknown", 100, 0, 1.7, 0.1, -100, 0, -100, -1, -1);

   // Creating N(1535)0
   new Particle("N(1535)0", 22112, 1, "Unknown", 100, 0, 1.535, 0.15, -100, 0, -100, -1, -1);

   // Creating delta(1920)0
   new Particle("delta(1920)0", 22114, 1, "Unknown", 100, 0, 1.92, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1910)+
   new Particle("delta(1910)+", 22122, 1, "Unknown", 100, 1, 1.91, 0.25, -100, 0, -100, -1, -1);

   // Creating N(1700)+
   new Particle("N(1700)+", 22124, 1, "Unknown", 100, 1, 1.7, 0.1, -100, 0, -100, -1, -1);

   // Creating N(1535)+
   new Particle("N(1535)+", 22212, 1, "Unknown", 100, 1, 1.535, 0.15, -100, 0, -100, -1, -1);

   // Creating delta(1920)+
   new Particle("delta(1920)+", 22214, 1, "Unknown", 100, 1, 1.92, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1910)++
   new Particle("delta(1910)++", 22222, 1, "Unknown", 100, 2, 1.91, 0.25, -100, 0, -100, -1, -1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
