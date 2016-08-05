// This files was autogenerated by geant::Particle::ReadFile

#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
#ifdef VECCORE_CUDA
#include "base/Vector.h"
template <typename T>
using vector = vecgeom::Vector<T>;
#else
using std::vector;
#endif

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void CreateParticle0049() {
   vector<int> daughters;
   Particle *part = nullptr;

   // Creating h'_1
   new Particle("h'_1", 10333, 0, "Unknown", 100, 0, 1.4, 0.08, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10333));
   daughters.clear();
   daughters.push_back(313);
   daughters.push_back(-311);
   part->AddDecay(Particle::Decay(0, 0.25,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(311);
   part->AddDecay(Particle::Decay(0, 0.25,  daughters));
   daughters.clear();
   daughters.push_back(323);
   daughters.push_back(-321);
   part->AddDecay(Particle::Decay(0, 0.25,  daughters));
   daughters.clear();
   daughters.push_back(-323);
   daughters.push_back(321);
   part->AddDecay(Particle::Decay(0, 0.25,  daughters));

   // Creating eta2(1870)
   new Particle("eta2(1870)", 10335, 1, "Unknown", 100, 0, 1.842, 0.225, -100, 0, -100, -1, -1);

   // Creating D*_0+
   new Particle("D*_0+", 10411, 1, "Unknown", 100, 1, 2.272, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10411));
   daughters.clear();
   daughters.push_back(421);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(411);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating D_1+
   new Particle("D_1+", 10413, 1, "Unknown", 100, 1, 2.424, 0.02, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10413));
   daughters.clear();
   daughters.push_back(423);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(413);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating D*_00
   new Particle("D*_00", 10421, 1, "Unknown", 100, 0, 2.272, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10421));
   daughters.clear();
   daughters.push_back(411);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(421);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating D_10
   new Particle("D_10", 10423, 1, "Unknown", 100, 0, 2.4223, 0.02, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10423));
   daughters.clear();
   daughters.push_back(413);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(423);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating D*_0s+
   new Particle("D*_0s+", 10431, 1, "Unknown", 100, 1, 2.3178, 0.0046, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10431));
   daughters.clear();
   daughters.push_back(431);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.8,  daughters));
   daughters.clear();
   daughters.push_back(431);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.2,  daughters));

   // Creating D_1s+
   new Particle("D_1s+", 10433, 1, "Unknown", 100, 1, 2.5353, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10433));
   daughters.clear();
   daughters.push_back(423);
   daughters.push_back(321);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(413);
   daughters.push_back(311);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));

   // Creating chi_0c
   new Particle("chi_0c", 10441, 0, "Unknown", 100, 0, 3.41475, 0.014, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10441));
   daughters.clear();
   daughters.push_back(82);
   daughters.push_back(-82);
   part->AddDecay(Particle::Decay(12, 0.993,  daughters));
   daughters.clear();
   daughters.push_back(443);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.007,  daughters));

   // Creating h_1c
   new Particle("h_1c", 10443, 0, "Unknown", 100, 0, 3.52593, 0.01, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10443));
   daughters.clear();
   daughters.push_back(82);
   daughters.push_back(-82);
   part->AddDecay(Particle::Decay(12, 1,  daughters));

   // Creating B*_00
   new Particle("B*_00", 10511, 1, "Unknown", 100, 0, 5.68, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10511));
   daughters.clear();
   daughters.push_back(521);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(511);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating B_10
   new Particle("B_10", 10513, 1, "Unknown", 100, 0, 5.73, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10513));
   daughters.clear();
   daughters.push_back(523);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(513);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating B*_0+
   new Particle("B*_0+", 10521, 1, "Unknown", 100, 1, 5.68, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10521));
   daughters.clear();
   daughters.push_back(511);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(521);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating B_1+
   new Particle("B_1+", 10523, 1, "Unknown", 100, 1, 5.73, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10523));
   daughters.clear();
   daughters.push_back(513);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(523);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating B*_0s0
   new Particle("B*_0s0", 10531, 1, "Unknown", 100, 0, 5.92, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10531));
   daughters.clear();
   daughters.push_back(521);
   daughters.push_back(-321);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(511);
   daughters.push_back(-311);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));

   // Creating B_1s0
   new Particle("B_1s0", 10533, 1, "Unknown", 100, 0, 5.97, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10533));
   daughters.clear();
   daughters.push_back(523);
   daughters.push_back(-321);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(513);
   daughters.push_back(-311);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));

   // Creating B*_0c+
   new Particle("B*_0c+", 10541, 1, "Unknown", 100, 1, 7.25, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10541));
   daughters.clear();
   daughters.push_back(511);
   daughters.push_back(411);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(521);
   daughters.push_back(421);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));

   // Creating B_1c+
   new Particle("B_1c+", 10543, 1, "Unknown", 100, 1, 7.3, 0.05, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10543));
   daughters.clear();
   daughters.push_back(513);
   daughters.push_back(411);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(523);
   daughters.push_back(421);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));

   // Creating chi_0b
   new Particle("chi_0b", 10551, 0, "Unknown", 100, 0, 9.8594, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10551));
   daughters.clear();
   daughters.push_back(21);
   daughters.push_back(21);
   part->AddDecay(Particle::Decay(32, 0.98,  daughters));
   daughters.clear();
   daughters.push_back(553);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.02,  daughters));

   // Creating h_1b
   new Particle("h_1b", 10553, 0, "Unknown", 100, 0, 9.875, 0.01, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(10553));
   daughters.clear();
   daughters.push_back(21);
   daughters.push_back(21);
   part->AddDecay(Particle::Decay(32, 1,  daughters));

   // Creating delta(1900)-
   new Particle("delta(1900)-", 11112, 1, "Unknown", 100, -1, 1.9, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1700)-
   new Particle("delta(1700)-", 11114, 1, "Unknown", 100, -1, 1.7, 0.3, -100, 0, -100, -1, -1);

   // Creating delta(1930)-
   new Particle("delta(1930)-", 11116, 1, "Unknown", 100, -1, 1.96, 0.36, -100, 0, -100, -1, -1);

   // Creating delta(1900)0
   new Particle("delta(1900)0", 11212, 1, "Unknown", 100, 0, 1.9, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1930)0
   new Particle("delta(1930)0", 11216, 1, "Unknown", 100, 0, 1.96, 0.36, -100, 0, -100, -1, -1);

   // Creating N(1440)0
   new Particle("N(1440)0", 12112, 1, "Unknown", 100, 0, 1.44, 0.3, -100, 0, -100, -1, -1);

   // Creating delta(1700)0
   new Particle("delta(1700)0", 12114, 1, "Unknown", 100, 0, 1.7, 0.3, -100, 0, -100, -1, -1);

   // Creating N(1680)0
   new Particle("N(1680)0", 12116, 1, "Unknown", 100, 0, 1.685, 0.13, -100, 0, -100, -1, -1);

   // Creating N(1990)0
   new Particle("N(1990)0", 12118, 1, "Unknown", 100, 0, 1.95, 0.555, -100, 0, -100, -1, -1);

   // Creating delta(1900)+
   new Particle("delta(1900)+", 12122, 1, "Unknown", 100, 1, 1.9, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1930)+
   new Particle("delta(1930)+", 12126, 1, "Unknown", 100, 1, 1.96, 0.36, -100, 0, -100, -1, -1);

   // Creating N(1440)+
   new Particle("N(1440)+", 12212, 1, "Unknown", 100, 1, 1.44, 0.3, -100, 0, -100, -1, -1);

   // Creating delta(1700)+
   new Particle("delta(1700)+", 12214, 1, "Unknown", 100, 1, 1.7, 0.3, -100, 0, -100, -1, -1);

   // Creating N(1680)+
   new Particle("N(1680)+", 12216, 1, "Unknown", 100, 1, 1.685, 0.13, -100, 0, -100, -1, -1);

   // Creating N(1990)+
   new Particle("N(1990)+", 12218, 1, "Unknown", 100, 1, 1.95, 0.555, -100, 0, -100, -1, -1);

   // Creating delta(1900)++
   new Particle("delta(1900)++", 12222, 1, "Unknown", 100, 2, 1.9, 0.2, -100, 0, -100, -1, -1);

   // Creating delta(1700)++
   new Particle("delta(1700)++", 12224, 1, "Unknown", 100, 2, 1.7, 0.3, -100, 0, -100, -1, -1);

   // Creating delta(1930)++
   new Particle("delta(1930)++", 12226, 1, "Unknown", 100, 2, 1.96, 0.36, -100, 0, -100, -1, -1);

   // Creating sigma(1660)-
   new Particle("sigma(1660)-", 13112, 1, "Unknown", 100, -1, 1.66, 0.1, -100, 0, -100, -1, -1);

   // Creating sigma(1670)-
   new Particle("sigma(1670)-", 13114, 1, "Unknown", 100, -1, 1.67, 0.06, -100, 0, -100, -1, -1);

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
   part = const_cast<Particle*>(&Particle::Particles().at(20113));
   daughters.clear();
   daughters.push_back(213);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(-213);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));

   // Creating a_1+
   new Particle("a_1+", 20213, 1, "Unknown", 100, 1, 1.23, 0.4, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20213));
   daughters.clear();
   daughters.push_back(113);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(213);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));

   // Creating f_1
   new Particle("f_1", 20223, 0, "Unknown", 100, 0, 1.2818, 0.025, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(20223));
   daughters.clear();
   daughters.push_back(113);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.15,  daughters));
   daughters.clear();
   daughters.push_back(10111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.146,  daughters));
   daughters.clear();
   daughters.push_back(-10211);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.146,  daughters));
   daughters.clear();
   daughters.push_back(10211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.146,  daughters));
   daughters.clear();
   daughters.push_back(113);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.066,  daughters));
   daughters.clear();
   daughters.push_back(213);
   daughters.push_back(-211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.05,  daughters));
   daughters.clear();
   daughters.push_back(221);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.05,  daughters));
   daughters.clear();
   daughters.push_back(113);
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.05,  daughters));
   daughters.clear();
   daughters.push_back(-213);
   daughters.push_back(211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.05,  daughters));
   daughters.clear();
   daughters.push_back(221);
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.05,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(-311);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.024,  daughters));
   daughters.clear();
   daughters.push_back(311);
   daughters.push_back(-311);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.024,  daughters));
   daughters.clear();
   daughters.push_back(311);
   daughters.push_back(-321);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.024,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(-321);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.024,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
