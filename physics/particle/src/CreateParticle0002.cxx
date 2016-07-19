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
void CreateParticle0002() {

   // Creating delta(1600)-_bar
   new Particle("delta(1600)-_bar", -31114, 0, "Unknown", 100, 1, 1.6, 0.35, 100, 100, 0, 100, 1);

   // Creating k_star(1680)-_bar
   new Particle("k_star(1680)-_bar", -30323, 0, "Unknown", 100, -1, 1.717, 0.32, 100, 100, 0, 100, 1);

   // Creating k_star(1680)0_bar
   new Particle("k_star(1680)0_bar", -30313, 0, "Unknown", 100, 0, 1.717, 0.32, 100, 100, 0, 100, 1);

   // Creating omega(1650)_bar
   new Particle("omega(1650)_bar", -30223, 0, "Unknown", 100, 0, 1.67, 0.315, 100, 100, 0, 100, 1);

   // Creating rho(1700)-_bar
   new Particle("rho(1700)-_bar", -30213, 0, "Unknown", 100, -1, 1.72, 0.25, 100, 100, 0, 100, 1);

   // Creating rho(1700)0_bar
   new Particle("rho(1700)0_bar", -30113, 0, "Unknown", 100, 0, 1.72, 0.25, 100, 100, 0, 100, 1);

   // Creating xi(1690)0_bar
   new Particle("xi(1690)0_bar", -23324, 0, "Unknown", 100, 0, 1.69, 0.05, 100, 100, 0, 100, 1);

   // Creating xi(1690)-_bar
   new Particle("xi(1690)-_bar", -23314, 0, "Unknown", 100, 1, 1.69, 0.05, 100, 100, 0, 100, 1);

   // Creating sigma(1940)+_bar
   new Particle("sigma(1940)+_bar", -23224, 0, "Unknown", 100, -1, 1.94, 0.22, 100, 100, 0, 100, 1);

   // Creating sigma(1750)+_bar
   new Particle("sigma(1750)+_bar", -23222, 0, "Unknown", 100, -1, 1.75, 0.09, 100, 100, 0, 100, 1);

   // Creating sigma(1940)0_bar
   new Particle("sigma(1940)0_bar", -23214, 0, "Unknown", 100, 0, 1.94, 0.22, 100, 100, 0, 100, 1);

   // Creating sigma(1750)0_bar
   new Particle("sigma(1750)0_bar", -23212, 0, "Unknown", 100, 0, 1.75, 0.09, 100, 100, 0, 100, 1);

   // Creating lambda(2110)_bar
   new Particle("lambda(2110)_bar", -23126, 0, "Unknown", 100, 0, 2.11, 0.2, 100, 100, 0, 100, 1);

   // Creating lambda(1890)_bar
   new Particle("lambda(1890)_bar", -23124, 0, "Unknown", 100, 0, 1.89, 0.1, 100, 100, 0, 100, 1);

   // Creating lambda(1600)_bar
   new Particle("lambda(1600)_bar", -23122, 0, "Unknown", 100, 0, 1.6, 0.15, 100, 100, 0, 100, 1);

   // Creating sigma(1940)-_bar
   new Particle("sigma(1940)-_bar", -23114, 0, "Unknown", 100, 1, 1.94, 0.22, 100, 100, 0, 100, 1);

   // Creating sigma(1750)-_bar
   new Particle("sigma(1750)-_bar", -23112, 0, "Unknown", 100, 1, 1.75, 0.09, 100, 100, 0, 100, 1);

   // Creating delta(1920)++_bar
   new Particle("delta(1920)++_bar", -22224, 0, "Unknown", 100, -2, 1.92, 0.2, 100, 100, 0, 100, 1);

   // Creating delta(1910)++_bar
   new Particle("delta(1910)++_bar", -22222, 0, "Unknown", 100, -2, 1.91, 0.25, 100, 100, 0, 100, 1);

   // Creating delta(1920)+_bar
   new Particle("delta(1920)+_bar", -22214, 0, "Unknown", 100, -1, 1.92, 0.2, 100, 100, 0, 100, 1);

   // Creating N(1535)+_bar
   new Particle("N(1535)+_bar", -22212, 0, "Unknown", 100, -1, 1.535, 0.15, 100, 100, 0, 100, 1);

   // Creating N(1700)+_bar
   new Particle("N(1700)+_bar", -22124, 0, "Unknown", 100, -1, 1.7, 0.1, 100, 100, 0, 100, 1);

   // Creating delta(1910)+_bar
   new Particle("delta(1910)+_bar", -22122, 0, "Unknown", 100, -1, 1.91, 0.25, 100, 100, 0, 100, 1);

   // Creating delta(1920)0_bar
   new Particle("delta(1920)0_bar", -22114, 0, "Unknown", 100, 0, 1.92, 0.2, 100, 100, 0, 100, 1);

   // Creating N(1535)0_bar
   new Particle("N(1535)0_bar", -22112, 0, "Unknown", 100, 0, 1.535, 0.15, 100, 100, 0, 100, 1);

   // Creating N(1700)0_bar
   new Particle("N(1700)0_bar", -21214, 0, "Unknown", 100, 0, 1.7, 0.1, 100, 100, 0, 100, 1);

   // Creating delta(1910)0_bar
   new Particle("delta(1910)0_bar", -21212, 0, "Unknown", 100, 0, 1.91, 0.25, 100, 100, 0, 100, 1);

   // Creating delta(1920)-_bar
   new Particle("delta(1920)-_bar", -21114, 0, "Unknown", 100, 1, 1.92, 0.2, 100, 100, 0, 100, 1);

   // Creating delta(1910)-_bar
   new Particle("delta(1910)-_bar", -21112, 0, "Unknown", 100, 1, 1.91, 0.25, 100, 100, 0, 100, 1);

   // Creating B*_1c-
   new Particle("B*_1c-", -20543, 0, "Unknown", 100, -1, 7.3, 0.05, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-20543));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-513,-411}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-523,-421}));

   // Creating B*_1s0_bar
   new Particle("B*_1s0_bar", -20533, 0, "Unknown", 100, 0, 6.02, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-20533));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-523,321}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-513,311}));

   // Creating B*_1-
   new Particle("B*_1-", -20523, 0, "Unknown", 100, -1, 5.78, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-20523));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-513,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-523,-111}));

   // Creating B*_10_bar
   new Particle("B*_10_bar", -20513, 0, "Unknown", 100, 0, 5.78, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-20513));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-523,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-513,-111}));

   // Creating D*_1s-
   new Particle("D*_1s-", -20433, 0, "Unknown", 100, -1, 2.4596, 0.0055, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-20433));
   part->AddDecay(Particle::Decay(0, 0.8,  vector<int>{-433,-111}));
   part->AddDecay(Particle::Decay(0, 0.2,  vector<int>{-433,-22}));

   // Creating D*_10_bar
   new Particle("D*_10_bar", -20423, 0, "Unknown", 100, 0, 2.372, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-20423));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-413,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-423,-111}));

   // Creating D*_1-
   new Particle("D*_1-", -20413, 0, "Unknown", 100, -1, 2.372, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-20413));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-423,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-413,-111}));

   // Creating K*_1-
   new Particle("K*_1-", -20323, 0, "Unknown", 100, -1, 1.403, 0.174, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-20323));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-313,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-323,-111}));

   // Creating K*_10_bar
   new Particle("K*_10_bar", -20313, 0, "Unknown", 100, 0, 1.403, 0.174, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-20313));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-323,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-313,-111}));

   // Creating a_1-
   new Particle("a_1-", -20213, 0, "Unknown", 100, -1, 1.23, 0.4, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-20213));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-113,-211}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-213,-111}));

   // Creating xi(2030)0_bar
   new Particle("xi(2030)0_bar", -13326, 0, "Unknown", 100, 0, 2.025, 0.02, 100, 100, 0, 100, 1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
