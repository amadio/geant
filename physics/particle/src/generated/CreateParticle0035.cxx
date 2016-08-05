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
void CreateParticle0035() {
   vector<int> daughters;
   Particle *part = nullptr;

   // Creating k3_star(1780)0
   new Particle("k3_star(1780)0", 317, 1, "Unknown", 100, 0, 1.776, 0.159, -100, 0, -100, -1, -1);

   // Creating K+
   new Particle("K+", 321, 1, "Meson", 100, 1, 0.493677, 5.31674e-17, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(321));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   part->AddDecay(Particle::Decay(0, 0.6352,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.2116,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.0559,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(42, 0.0482,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(42, 0.0318,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.0173,  daughters));

   // Creating K*+
   new Particle("K*+", 323, 1, "Meson", 100, 1, 0.89166, 0.0498, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(323));
   daughters.clear();
   daughters.push_back(311);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(3, 0.666,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(3, 0.333,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.001,  daughters));

   // Creating K*_2+
   new Particle("K*_2+", 325, 1, "Meson", 100, 1, 1.4256, 0.098, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(325));
   daughters.clear();
   daughters.push_back(311);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.332,  daughters));
   daughters.clear();
   daughters.push_back(313);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.168,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.166,  daughters));
   daughters.clear();
   daughters.push_back(313);
   daughters.push_back(211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.086,  daughters));
   daughters.clear();
   daughters.push_back(323);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.084,  daughters));
   daughters.clear();
   daughters.push_back(311);
   daughters.push_back(213);
   part->AddDecay(Particle::Decay(0, 0.059,  daughters));
   daughters.clear();
   daughters.push_back(323);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.043,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(113);
   part->AddDecay(Particle::Decay(0, 0.029,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(223);
   part->AddDecay(Particle::Decay(0, 0.029,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(221);
   part->AddDecay(Particle::Decay(0, 0.002,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.002,  daughters));

   // Creating k3_star(1780)+
   new Particle("k3_star(1780)+", 327, 1, "Unknown", 100, 1, 1.776, 0.159, -100, 0, -100, -1, -1);

   // Creating phi_diff
   new Particle("phi_diff", 330, 0, "Meson", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating eta'
   new Particle("eta'", 331, 0, "Meson", 100, 0, 0.95766, 0.0002, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(331));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(221);
   part->AddDecay(Particle::Decay(0, 0.437,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(113);
   part->AddDecay(Particle::Decay(0, 0.302,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(111);
   daughters.push_back(221);
   part->AddDecay(Particle::Decay(0, 0.208,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(223);
   part->AddDecay(Particle::Decay(0, 0.0302,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.0212,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.0016,  daughters));

   // Creating phi
   new Particle("phi", 333, 0, "Meson", 100, 0, 1.01945, 0.00443, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(333));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(-321);
   part->AddDecay(Particle::Decay(3, 0.48947,  daughters));
   daughters.clear();
   daughters.push_back(130);
   daughters.push_back(310);
   part->AddDecay(Particle::Decay(3, 0.34,  daughters));
   daughters.clear();
   daughters.push_back(-213);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.043,  daughters));
   daughters.clear();
   daughters.push_back(113);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.043,  daughters));
   daughters.clear();
   daughters.push_back(213);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.043,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(1, 0.027,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(221);
   part->AddDecay(Particle::Decay(0, 0.0126,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.0013,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-11);
   part->AddDecay(Particle::Decay(0, 0.0003,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-13);
   part->AddDecay(Particle::Decay(0, 0.00025,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 8e-05,  daughters));

   // Creating f'_2
   new Particle("f'_2", 335, 0, "Meson", 100, 0, 1.525, 0.076, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(335));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(-321);
   part->AddDecay(Particle::Decay(0, 0.444,  daughters));
   daughters.clear();
   daughters.push_back(130);
   daughters.push_back(130);
   part->AddDecay(Particle::Decay(0, 0.222,  daughters));
   daughters.clear();
   daughters.push_back(310);
   daughters.push_back(310);
   part->AddDecay(Particle::Decay(0, 0.222,  daughters));
   daughters.clear();
   daughters.push_back(221);
   daughters.push_back(221);
   part->AddDecay(Particle::Decay(0, 0.104,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.004,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.004,  daughters));

   // Creating phi3(1850)
   new Particle("phi3(1850)", 337, 1, "Unknown", 100, 0, 1.854, 0.087, -100, 0, -100, -1, -1);

   // Creating D+
   new Particle("D+", 411, 1, "CharmedMeson", 100, 1, 1.86962, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(411));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.087,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(20213);
   part->AddDecay(Particle::Decay(0, 0.076,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   daughters.push_back(-311);
   part->AddDecay(Particle::Decay(42, 0.07,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   daughters.push_back(-311);
   part->AddDecay(Particle::Decay(42, 0.07,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(211);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.067,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(213);
   part->AddDecay(Particle::Decay(0, 0.066,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   daughters.push_back(-313);
   part->AddDecay(Particle::Decay(42, 0.065,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   daughters.push_back(-313);
   part->AddDecay(Particle::Decay(42, 0.065,  daughters));
   daughters.clear();
   daughters.push_back(-20313);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.045,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(213);
   part->AddDecay(Particle::Decay(0, 0.041,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(321);
   daughters.push_back(-311);
   part->AddDecay(Particle::Decay(0, 0.027,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(323);
   part->AddDecay(Particle::Decay(0, 0.026,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.026,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.022,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.0218,  daughters));
   daughters.clear();
   daughters.push_back(333);
   daughters.push_back(211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.019,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.019,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.012,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.012,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   daughters.push_back(-323);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(42, 0.011,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   daughters.push_back(-313);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(42, 0.011,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   daughters.push_back(-323);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(42, 0.011,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   daughters.push_back(-313);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(42, 0.011,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.009,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(213);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.008,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(321);
   part->AddDecay(Particle::Decay(0, 0.0073,  daughters));
   daughters.clear();
   daughters.push_back(221);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.0066,  daughters));
   daughters.clear();
   daughters.push_back(333);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.006,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(211);
   daughters.push_back(113);
   part->AddDecay(Particle::Decay(0, 0.0057,  daughters));
   daughters.clear();
   daughters.push_back(333);
   daughters.push_back(213);
   part->AddDecay(Particle::Decay(0, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   daughters.push_back(-311);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   daughters.push_back(-321);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   daughters.push_back(-311);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(221);
   daughters.push_back(213);
   part->AddDecay(Particle::Decay(0, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   daughters.push_back(-321);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(42, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(321);
   part->AddDecay(Particle::Decay(0, 0.0047,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(323);
   part->AddDecay(Particle::Decay(0, 0.0047,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(321);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.004,  daughters));
   daughters.clear();
   daughters.push_back(331);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.003,  daughters));
   daughters.clear();
   daughters.push_back(331);
   daughters.push_back(213);
   part->AddDecay(Particle::Decay(0, 0.003,  daughters));
   daughters.clear();
   daughters.push_back(113);
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.0028,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.0022,  daughters));
   daughters.clear();
   daughters.push_back(-313);
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.002,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(113);
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.0019,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.0015,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   daughters.push_back(221);
   part->AddDecay(Particle::Decay(42, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   daughters.push_back(331);
   part->AddDecay(Particle::Decay(42, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   daughters.push_back(113);
   part->AddDecay(Particle::Decay(42, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   daughters.push_back(223);
   part->AddDecay(Particle::Decay(42, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(223);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(223);
   daughters.push_back(213);
   part->AddDecay(Particle::Decay(0, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(42, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-321);
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(42, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   daughters.push_back(221);
   part->AddDecay(Particle::Decay(42, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-311);
   daughters.push_back(113);
   daughters.push_back(211);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   daughters.push_back(331);
   part->AddDecay(Particle::Decay(42, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   daughters.push_back(113);
   part->AddDecay(Particle::Decay(42, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   daughters.push_back(223);
   part->AddDecay(Particle::Decay(42, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(113);
   daughters.push_back(213);
   part->AddDecay(Particle::Decay(0, 0.0006,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(213);
   part->AddDecay(Particle::Decay(0, 0.0006,  daughters));
   daughters.clear();
   daughters.push_back(113);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.0006,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
