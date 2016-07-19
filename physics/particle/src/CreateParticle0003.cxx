#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
using std::vector;
namespace geant {
   inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
void CreateParticle0003() {

   // Creating xi(1820)0_bar
   new Particle("xi(1820)0_bar", -13324, 0, "Unknown", 100, 0, 1.823, 0.024, 100, 100, 0, 100, 1);

   // Creating xi(2030)-_bar
   new Particle("xi(2030)-_bar", -13316, 0, "Unknown", 100, 1, 2.025, 0.02, 100, 100, 0, 100, 1);

   // Creating xi(1820)-_bar
   new Particle("xi(1820)-_bar", -13314, 0, "Unknown", 100, 1, 1.823, 0.024, 100, 100, 0, 100, 1);

   // Creating sigma(1915)+_bar
   new Particle("sigma(1915)+_bar", -13226, 0, "Unknown", 100, -1, 1.915, 0.12, 100, 100, 0, 100, 1);

   // Creating sigma(1670)+_bar
   new Particle("sigma(1670)+_bar", -13224, 0, "Unknown", 100, -1, 1.67, 0.06, 100, 100, 0, 100, 1);

   // Creating sigma(1660)+_bar
   new Particle("sigma(1660)+_bar", -13222, 0, "Unknown", 100, -1, 1.66, 0.1, 100, 100, 0, 100, 1);

   // Creating sigma(1915)0_bar
   new Particle("sigma(1915)0_bar", -13216, 0, "Unknown", 100, 0, 1.915, 0.12, 100, 100, 0, 100, 1);

   // Creating sigma(1670)0_bar
   new Particle("sigma(1670)0_bar", -13214, 0, "Unknown", 100, 0, 1.67, 0.06, 100, 100, 0, 100, 1);

   // Creating sigma(1660)0_bar
   new Particle("sigma(1660)0_bar", -13212, 0, "Unknown", 100, 0, 1.66, 0.1, 100, 100, 0, 100, 1);

   // Creating lambda(1830)_bar
   new Particle("lambda(1830)_bar", -13126, 0, "Unknown", 100, 0, 1.83, 0.095, 100, 100, 0, 100, 1);

   // Creating lambda(1690)_bar
   new Particle("lambda(1690)_bar", -13124, 0, "Unknown", 100, 0, 1.69, 0.06, 100, 100, 0, 100, 1);

   // Creating lambda(1405)_bar
   new Particle("lambda(1405)_bar", -13122, 0, "Unknown", 100, 0, 1.4051, 0.05, 100, 100, 0, 100, 1);

   // Creating sigma(1915)-_bar
   new Particle("sigma(1915)-_bar", -13116, 0, "Unknown", 100, 1, 1.915, 0.12, 100, 100, 0, 100, 1);

   // Creating sigma(1670)-_bar
   new Particle("sigma(1670)-_bar", -13114, 0, "Unknown", 100, 1, 1.67, 0.06, 100, 100, 0, 100, 1);

   // Creating sigma(1660)-_bar
   new Particle("sigma(1660)-_bar", -13112, 0, "Unknown", 100, 1, 1.66, 0.1, 100, 100, 0, 100, 1);

   // Creating delta(1930)++_bar
   new Particle("delta(1930)++_bar", -12226, 0, "Unknown", 100, -2, 1.96, 0.36, 100, 100, 0, 100, 1);

   // Creating delta(1700)++_bar
   new Particle("delta(1700)++_bar", -12224, 0, "Unknown", 100, -2, 1.7, 0.3, 100, 100, 0, 100, 1);

   // Creating delta(1900)++_bar
   new Particle("delta(1900)++_bar", -12222, 0, "Unknown", 100, -2, 1.9, 0.2, 100, 100, 0, 100, 1);

   // Creating N(1990)+_bar
   new Particle("N(1990)+_bar", -12218, 0, "Unknown", 100, -1, 1.95, 0.555, 100, 100, 0, 100, 1);

   // Creating N(1680)+_bar
   new Particle("N(1680)+_bar", -12216, 0, "Unknown", 100, -1, 1.685, 0.13, 100, 100, 0, 100, 1);

   // Creating delta(1700)+_bar
   new Particle("delta(1700)+_bar", -12214, 0, "Unknown", 100, -1, 1.7, 0.3, 100, 100, 0, 100, 1);

   // Creating N(1440)+_bar
   new Particle("N(1440)+_bar", -12212, 0, "Unknown", 100, -1, 1.44, 0.3, 100, 100, 0, 100, 1);

   // Creating delta(1930)+_bar
   new Particle("delta(1930)+_bar", -12126, 0, "Unknown", 100, -1, 1.96, 0.36, 100, 100, 0, 100, 1);

   // Creating delta(1900)+_bar
   new Particle("delta(1900)+_bar", -12122, 0, "Unknown", 100, -1, 1.9, 0.2, 100, 100, 0, 100, 1);

   // Creating N(1990)0_bar
   new Particle("N(1990)0_bar", -12118, 0, "Unknown", 100, 0, 1.95, 0.555, 100, 100, 0, 100, 1);

   // Creating N(1680)0_bar
   new Particle("N(1680)0_bar", -12116, 0, "Unknown", 100, 0, 1.685, 0.13, 100, 100, 0, 100, 1);

   // Creating delta(1700)0_bar
   new Particle("delta(1700)0_bar", -12114, 0, "Unknown", 100, 0, 1.7, 0.3, 100, 100, 0, 100, 1);

   // Creating N(1440)0_bar
   new Particle("N(1440)0_bar", -12112, 0, "Unknown", 100, 0, 1.44, 0.3, 100, 100, 0, 100, 1);

   // Creating delta(1930)0_bar
   new Particle("delta(1930)0_bar", -11216, 0, "Unknown", 100, 0, 1.96, 0.36, 100, 100, 0, 100, 1);

   // Creating delta(1900)0_bar
   new Particle("delta(1900)0_bar", -11212, 0, "Unknown", 100, 0, 1.9, 0.2, 100, 100, 0, 100, 1);

   // Creating delta(1930)-_bar
   new Particle("delta(1930)-_bar", -11116, 0, "Unknown", 100, 1, 1.96, 0.36, 100, 100, 0, 100, 1);

   // Creating delta(1700)-_bar
   new Particle("delta(1700)-_bar", -11114, 0, "Unknown", 100, 1, 1.7, 0.3, 100, 100, 0, 100, 1);

   // Creating delta(1900)-_bar
   new Particle("delta(1900)-_bar", -11112, 0, "Unknown", 100, 1, 1.9, 0.2, 100, 100, 0, 100, 1);

   // Creating B_1c-
   new Particle("B_1c-", -10543, 0, "Unknown", 100, -1, 7.3, 0.05, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-10543));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-513,-411}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-523,-421}));

   // Creating B*_0c-
   new Particle("B*_0c-", -10541, 0, "Unknown", 100, -1, 7.25, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10541));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-511,-411}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-521,-421}));

   // Creating B_1s0_bar
   new Particle("B_1s0_bar", -10533, 0, "Unknown", 100, 0, 5.97, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10533));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-523,321}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-513,311}));

   // Creating B*_0s0_bar
   new Particle("B*_0s0_bar", -10531, 0, "Unknown", 100, 0, 5.92, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10531));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-521,321}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-511,311}));

   // Creating B_1-
   new Particle("B_1-", -10523, 0, "Unknown", 100, -1, 5.73, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10523));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-513,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-523,-111}));

   // Creating B*_0-
   new Particle("B*_0-", -10521, 0, "Unknown", 100, -1, 5.68, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10521));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-511,-211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-521,-111}));

   // Creating B_10_bar
   new Particle("B_10_bar", -10513, 0, "Unknown", 100, 0, 5.73, 0.05, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-10513));
   part->AddDecay(Particle::Decay(0, 0.667,  vector<int>{-523,211}));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-513,-111}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
