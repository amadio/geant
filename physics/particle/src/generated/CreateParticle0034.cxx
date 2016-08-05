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
void CreateParticle0034() {
   vector<int> daughters;
   Particle *part = nullptr;

   // Creating rho3(1690)0
   new Particle("rho3(1690)0", 117, 1, "Unknown", 100, 0, 1.6888, 0.161, -100, 0, -100, -1, -1);

   // Creating K_L0
   new Particle("K_L0", 130, 0, "Meson", 100, 0, 0.497614, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(130));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.2112,  daughters));
   daughters.clear();
   daughters.push_back(-12);
   daughters.push_back(11);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(42, 0.1939,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(42, 0.1939,  daughters));
   daughters.clear();
   daughters.push_back(-14);
   daughters.push_back(13);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(42, 0.1359,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(42, 0.1359,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.1256,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.002,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.001,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.0006,  daughters));

   // Creating pi_diffr+
   new Particle("pi_diffr+", 210, 1, "Meson", 100, 1, 0, 0, -100, -1, -100, -1, -1);

   // Creating pi+
   new Particle("pi+", 211, 1, "Meson", 100, 1, 0.13957, 2.52837e-17, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(211));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   part->AddDecay(Particle::Decay(0, 0.999877,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   part->AddDecay(Particle::Decay(0, 0.000123,  daughters));

   // Creating rho+
   new Particle("rho+", 213, 1, "Meson", 100, 1, 0.77549, 0.149, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(213));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(3, 0.99955,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.00045,  daughters));

   // Creating a_2+
   new Particle("a_2+", 215, 1, "Meson", 100, 1, 1.3183, 0.107, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(215));
   daughters.clear();
   daughters.push_back(213);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.34725,  daughters));
   daughters.clear();
   daughters.push_back(113);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.34725,  daughters));
   daughters.clear();
   daughters.push_back(221);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.144,  daughters));
   daughters.clear();
   daughters.push_back(223);
   daughters.push_back(211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.104,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(-311);
   part->AddDecay(Particle::Decay(0, 0.049,  daughters));
   daughters.clear();
   daughters.push_back(331);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.0057,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.0028,  daughters));

   // Creating rho3(1690)+
   new Particle("rho3(1690)+", 217, 1, "Unknown", 100, 1, 1.6888, 0.161, -100, 0, -100, -1, -1);

   // Creating omega_di
   new Particle("omega_di", 220, 0, "Meson", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating eta
   new Particle("eta", 221, 0, "Meson", 100, 0, 0.547853, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(221));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.3923,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.321,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.2317,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.0478,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(11);
   daughters.push_back(-11);
   part->AddDecay(Particle::Decay(2, 0.0049,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(11);
   daughters.push_back(-11);
   part->AddDecay(Particle::Decay(0, 0.0013,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(22);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.0007,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(13);
   daughters.push_back(-13);
   part->AddDecay(Particle::Decay(0, 0.0003,  daughters));

   // Creating omega
   new Particle("omega", 223, 0, "Meson", 100, 0, 0.78265, 0.00843, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(223));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(1, 0.89,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.08693,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(3, 0.0221,  daughters));
   daughters.clear();
   daughters.push_back(221);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.00083,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(111);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 7e-05,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-11);
   part->AddDecay(Particle::Decay(0, 7e-05,  daughters));

   // Creating f_2
   new Particle("f_2", 225, 0, "Meson", 100, 0, 1.2751, 0.185, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(225));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.564,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.282,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.072,  daughters));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.028,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(-321);
   part->AddDecay(Particle::Decay(0, 0.023,  daughters));
   daughters.clear();
   daughters.push_back(130);
   daughters.push_back(130);
   part->AddDecay(Particle::Decay(0, 0.0115,  daughters));
   daughters.clear();
   daughters.push_back(310);
   daughters.push_back(310);
   part->AddDecay(Particle::Decay(0, 0.0115,  daughters));
   daughters.clear();
   daughters.push_back(221);
   daughters.push_back(221);
   part->AddDecay(Particle::Decay(0, 0.005,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(111);
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.003,  daughters));

   // Creating omega3(1670)
   new Particle("omega3(1670)", 227, 1, "Unknown", 100, 0, 1.667, 0.168, -100, 0, -100, -1, -1);

   // Creating K_S0
   new Particle("K_S0", 310, 0, "Meson", 100, 0, 0.497614, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(310));
   daughters.clear();
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.6861,  daughters));
   daughters.clear();
   daughters.push_back(111);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.3139,  daughters));

   // Creating K0
   new Particle("K0", 311, 1, "Meson", 100, 0, 0.497614, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(311));
   daughters.clear();
   daughters.push_back(130);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));
   daughters.clear();
   daughters.push_back(310);
   part->AddDecay(Particle::Decay(0, 0.5,  daughters));

   // Creating K*0
   new Particle("K*0", 313, 1, "Meson", 100, 0, 0.896, 0.0505, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(313));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(3, 0.665,  daughters));
   daughters.clear();
   daughters.push_back(311);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(3, 0.333,  daughters));
   daughters.clear();
   daughters.push_back(311);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0.002,  daughters));

   // Creating K*_20
   new Particle("K*_20", 315, 1, "Meson", 100, 0, 1.4324, 0.109, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(315));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));
   daughters.clear();
   daughters.push_back(323);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.168,  daughters));
   daughters.clear();
   daughters.push_back(311);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.166,  daughters));
   daughters.clear();
   daughters.push_back(323);
   daughters.push_back(-211);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.087,  daughters));
   daughters.clear();
   daughters.push_back(313);
   daughters.push_back(111);
   part->AddDecay(Particle::Decay(0, 0.084,  daughters));
   daughters.clear();
   daughters.push_back(321);
   daughters.push_back(-213);
   part->AddDecay(Particle::Decay(0, 0.059,  daughters));
   daughters.clear();
   daughters.push_back(313);
   daughters.push_back(211);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.043,  daughters));
   daughters.clear();
   daughters.push_back(311);
   daughters.push_back(113);
   part->AddDecay(Particle::Decay(0, 0.029,  daughters));
   daughters.clear();
   daughters.push_back(311);
   daughters.push_back(223);
   part->AddDecay(Particle::Decay(0, 0.029,  daughters));
   daughters.clear();
   daughters.push_back(311);
   daughters.push_back(221);
   part->AddDecay(Particle::Decay(0, 0.002,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
