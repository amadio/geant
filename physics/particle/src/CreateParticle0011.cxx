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
void CreateParticle0011() {

   // Creating rho3(1690)0
   new Particle("rho3(1690)0", 117, 1, "Unknown", 100, 0, 1.6888, 0.161, -100, 0, -100, -1, -1);

   // Creating K_L0
   new Particle("K_L0", 130, 0, "Meson", 100, 0, 0.497614, 0, -100, -1, -100, -1, -1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(130));
   part->AddDecay(Particle::Decay(0, 0.2112,  vector<int>{111,111,111}));
   part->AddDecay(Particle::Decay(42, 0.1939,  vector<int>{-12,11,211}));
   part->AddDecay(Particle::Decay(42, 0.1939,  vector<int>{12,-11,-211}));
   part->AddDecay(Particle::Decay(42, 0.1359,  vector<int>{-14,13,211}));
   part->AddDecay(Particle::Decay(42, 0.1359,  vector<int>{14,-13,-211}));
   part->AddDecay(Particle::Decay(0, 0.1256,  vector<int>{211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{211,-211}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{111,111}));
   part->AddDecay(Particle::Decay(0, 0.0006,  vector<int>{22,22}));

   // Creating pi_diffr+
   new Particle("pi_diffr+", 210, 1, "Meson", 100, 1, 0, 0, -100, -1, -100, -1, -1);

   // Creating pi+
   new Particle("pi+", 211, 1, "Meson", 100, 1, 0.13957, 2.52837e-17, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(211));
   part->AddDecay(Particle::Decay(0, 0.999877,  vector<int>{-13,14}));
   part->AddDecay(Particle::Decay(0, 0.000123,  vector<int>{-11,12}));

   // Creating rho+
   new Particle("rho+", 213, 1, "Meson", 100, 1, 0.77549, 0.149, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(213));
   part->AddDecay(Particle::Decay(3, 0.99955,  vector<int>{211,111}));
   part->AddDecay(Particle::Decay(0, 0.00045,  vector<int>{211,22}));

   // Creating a_2+
   new Particle("a_2+", 215, 1, "Meson", 100, 1, 1.3183, 0.107, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(215));
   part->AddDecay(Particle::Decay(0, 0.34725,  vector<int>{213,111}));
   part->AddDecay(Particle::Decay(0, 0.34725,  vector<int>{113,211}));
   part->AddDecay(Particle::Decay(0, 0.144,  vector<int>{221,211}));
   part->AddDecay(Particle::Decay(0, 0.104,  vector<int>{223,211,111}));
   part->AddDecay(Particle::Decay(0, 0.049,  vector<int>{321,-311}));
   part->AddDecay(Particle::Decay(0, 0.0057,  vector<int>{331,211}));
   part->AddDecay(Particle::Decay(0, 0.0028,  vector<int>{211,22}));

   // Creating rho3(1690)+
   new Particle("rho3(1690)+", 217, 1, "Unknown", 100, 1, 1.6888, 0.161, -100, 0, -100, -1, -1);

   // Creating omega_di
   new Particle("omega_di", 220, 0, "Meson", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating eta
   new Particle("eta", 221, 0, "Meson", 100, 0, 0.547853, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(221));
   part->AddDecay(Particle::Decay(0, 0.3923,  vector<int>{22,22}));
   part->AddDecay(Particle::Decay(0, 0.321,  vector<int>{111,111,111}));
   part->AddDecay(Particle::Decay(0, 0.2317,  vector<int>{211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.0478,  vector<int>{22,211,-211}));
   part->AddDecay(Particle::Decay(2, 0.0049,  vector<int>{22,11,-11}));
   part->AddDecay(Particle::Decay(0, 0.0013,  vector<int>{211,-211,11,-11}));
   part->AddDecay(Particle::Decay(0, 0.0007,  vector<int>{111,22,22}));
   part->AddDecay(Particle::Decay(0, 0.0003,  vector<int>{22,13,-13}));

   // Creating omega
   new Particle("omega", 223, 0, "Meson", 100, 0, 0.78265, 0.00843, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(223));
   part->AddDecay(Particle::Decay(1, 0.89,  vector<int>{211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.08693,  vector<int>{22,111}));
   part->AddDecay(Particle::Decay(3, 0.0221,  vector<int>{211,-211}));
   part->AddDecay(Particle::Decay(0, 0.00083,  vector<int>{221,22}));
   part->AddDecay(Particle::Decay(0, 7e-05,  vector<int>{111,111,22}));
   part->AddDecay(Particle::Decay(0, 7e-05,  vector<int>{11,-11}));

   // Creating f_2
   new Particle("f_2", 225, 0, "Meson", 100, 0, 1.2751, 0.185, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(225));
   part->AddDecay(Particle::Decay(0, 0.564,  vector<int>{211,-211}));
   part->AddDecay(Particle::Decay(0, 0.282,  vector<int>{111,111}));
   part->AddDecay(Particle::Decay(0, 0.072,  vector<int>{211,-211,111,111}));
   part->AddDecay(Particle::Decay(0, 0.028,  vector<int>{211,-211,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.023,  vector<int>{321,-321}));
   part->AddDecay(Particle::Decay(0, 0.0115,  vector<int>{130,130}));
   part->AddDecay(Particle::Decay(0, 0.0115,  vector<int>{310,310}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{221,221}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{111,111,111,111}));

   // Creating omega3(1670)
   new Particle("omega3(1670)", 227, 1, "Unknown", 100, 0, 1.667, 0.168, -100, 0, -100, -1, -1);

   // Creating K_S0
   new Particle("K_S0", 310, 0, "Meson", 100, 0, 0.497614, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(310));
   part->AddDecay(Particle::Decay(0, 0.6861,  vector<int>{211,-211}));
   part->AddDecay(Particle::Decay(0, 0.3139,  vector<int>{111,111}));

   // Creating K0
   new Particle("K0", 311, 1, "Meson", 100, 0, 0.497614, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(311));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{130}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{310}));

   // Creating K*0
   new Particle("K*0", 313, 1, "Meson", 100, 0, 0.896, 0.0505, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(313));
   part->AddDecay(Particle::Decay(3, 0.665,  vector<int>{321,-211}));
   part->AddDecay(Particle::Decay(3, 0.333,  vector<int>{311,111}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{311,22}));

   // Creating K*_20
   new Particle("K*_20", 315, 1, "Meson", 100, 0, 1.4324, 0.109, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(315));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{321,-211}));
   part->AddDecay(Particle::Decay(0, 0.168,  vector<int>{323,-211}));
   part->AddDecay(Particle::Decay(0, 0.166,  vector<int>{311,111}));
   part->AddDecay(Particle::Decay(0, 0.087,  vector<int>{323,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.084,  vector<int>{313,111}));
   part->AddDecay(Particle::Decay(0, 0.059,  vector<int>{321,-213}));
   part->AddDecay(Particle::Decay(0, 0.043,  vector<int>{313,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.029,  vector<int>{311,113}));
   part->AddDecay(Particle::Decay(0, 0.029,  vector<int>{311,223}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{311,221}));

   // Creating k3_star(1780)0
   new Particle("k3_star(1780)0", 317, 1, "Unknown", 100, 0, 1.776, 0.159, -100, 0, -100, -1, -1);

   // Creating K+
   new Particle("K+", 321, 1, "Meson", 100, 1, 0.493677, 5.31674e-17, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(321));
   part->AddDecay(Particle::Decay(0, 0.6352,  vector<int>{-13,14}));
   part->AddDecay(Particle::Decay(0, 0.2116,  vector<int>{211,111}));
   part->AddDecay(Particle::Decay(0, 0.0559,  vector<int>{211,211,-211}));
   part->AddDecay(Particle::Decay(42, 0.0482,  vector<int>{12,-11,111}));
   part->AddDecay(Particle::Decay(42, 0.0318,  vector<int>{14,-13,111}));
   part->AddDecay(Particle::Decay(0, 0.0173,  vector<int>{211,111,111}));

   // Creating K*+
   new Particle("K*+", 323, 1, "Meson", 100, 1, 0.89166, 0.0498, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(323));
   part->AddDecay(Particle::Decay(3, 0.666,  vector<int>{311,211}));
   part->AddDecay(Particle::Decay(3, 0.333,  vector<int>{321,111}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{321,22}));

   // Creating K*_2+
   new Particle("K*_2+", 325, 1, "Meson", 100, 1, 1.4256, 0.098, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(325));
   part->AddDecay(Particle::Decay(0, 0.332,  vector<int>{311,211}));
   part->AddDecay(Particle::Decay(0, 0.168,  vector<int>{313,211}));
   part->AddDecay(Particle::Decay(0, 0.166,  vector<int>{321,111}));
   part->AddDecay(Particle::Decay(0, 0.086,  vector<int>{313,211,111}));
   part->AddDecay(Particle::Decay(0, 0.084,  vector<int>{323,111}));
   part->AddDecay(Particle::Decay(0, 0.059,  vector<int>{311,213}));
   part->AddDecay(Particle::Decay(0, 0.043,  vector<int>{323,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.029,  vector<int>{321,113}));
   part->AddDecay(Particle::Decay(0, 0.029,  vector<int>{321,223}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{321,221}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{321,22}));

   // Creating k3_star(1780)+
   new Particle("k3_star(1780)+", 327, 1, "Unknown", 100, 1, 1.776, 0.159, -100, 0, -100, -1, -1);

   // Creating phi_diff
   new Particle("phi_diff", 330, 0, "Meson", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating eta'
   new Particle("eta'", 331, 0, "Meson", 100, 0, 0.95766, 0.0002, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(331));
   part->AddDecay(Particle::Decay(0, 0.437,  vector<int>{211,-211,221}));
   part->AddDecay(Particle::Decay(0, 0.302,  vector<int>{22,113}));
   part->AddDecay(Particle::Decay(0, 0.208,  vector<int>{111,111,221}));
   part->AddDecay(Particle::Decay(0, 0.0302,  vector<int>{22,223}));
   part->AddDecay(Particle::Decay(0, 0.0212,  vector<int>{22,22}));
   part->AddDecay(Particle::Decay(0, 0.0016,  vector<int>{111,111,111}));

   // Creating phi
   new Particle("phi", 333, 0, "Meson", 100, 0, 1.01945, 0.00443, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(333));
   part->AddDecay(Particle::Decay(3, 0.48947,  vector<int>{321,-321}));
   part->AddDecay(Particle::Decay(3, 0.34,  vector<int>{130,310}));
   part->AddDecay(Particle::Decay(0, 0.043,  vector<int>{-213,211}));
   part->AddDecay(Particle::Decay(0, 0.043,  vector<int>{113,111}));
   part->AddDecay(Particle::Decay(0, 0.043,  vector<int>{213,-211}));
   part->AddDecay(Particle::Decay(1, 0.027,  vector<int>{211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.0126,  vector<int>{22,221}));
   part->AddDecay(Particle::Decay(0, 0.0013,  vector<int>{111,22}));
   part->AddDecay(Particle::Decay(0, 0.0003,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(0, 0.00025,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(0, 8e-05,  vector<int>{211,-211}));

   // Creating f'_2
   new Particle("f'_2", 335, 0, "Meson", 100, 0, 1.525, 0.076, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(335));
   part->AddDecay(Particle::Decay(0, 0.444,  vector<int>{321,-321}));
   part->AddDecay(Particle::Decay(0, 0.222,  vector<int>{130,130}));
   part->AddDecay(Particle::Decay(0, 0.222,  vector<int>{310,310}));
   part->AddDecay(Particle::Decay(0, 0.104,  vector<int>{221,221}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{211,-211}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{111,111}));

   // Creating phi3(1850)
   new Particle("phi3(1850)", 337, 1, "Unknown", 100, 0, 1.854, 0.087, -100, 0, -100, -1, -1);

   // Creating D+
   new Particle("D+", 411, 1, "CharmedMeson", 100, 1, 1.86962, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(411));
   part->AddDecay(Particle::Decay(0, 0.087,  vector<int>{-311,211,211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.076,  vector<int>{-311,20213}));
   part->AddDecay(Particle::Decay(42, 0.07,  vector<int>{-11,12,-311}));
   part->AddDecay(Particle::Decay(42, 0.07,  vector<int>{-13,14,-311}));
   part->AddDecay(Particle::Decay(0, 0.067,  vector<int>{-321,211,211}));
   part->AddDecay(Particle::Decay(0, 0.066,  vector<int>{-311,213}));
   part->AddDecay(Particle::Decay(42, 0.065,  vector<int>{-11,12,-313}));
   part->AddDecay(Particle::Decay(42, 0.065,  vector<int>{-13,14,-313}));
   part->AddDecay(Particle::Decay(0, 0.045,  vector<int>{-20313,211}));
   part->AddDecay(Particle::Decay(0, 0.041,  vector<int>{-313,213}));
   part->AddDecay(Particle::Decay(0, 0.027,  vector<int>{-311,321,-311}));
   part->AddDecay(Particle::Decay(0, 0.026,  vector<int>{-313,323}));
   part->AddDecay(Particle::Decay(0, 0.026,  vector<int>{-311,211}));
   part->AddDecay(Particle::Decay(0, 0.022,  vector<int>{-321,211,211,111,111}));
   part->AddDecay(Particle::Decay(0, 0.0218,  vector<int>{211,211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.019,  vector<int>{333,211,111}));
   part->AddDecay(Particle::Decay(0, 0.019,  vector<int>{-313,211}));
   part->AddDecay(Particle::Decay(0, 0.012,  vector<int>{-311,211,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.012,  vector<int>{-311,211,111}));
   part->AddDecay(Particle::Decay(42, 0.011,  vector<int>{-13,14,-323,211}));
   part->AddDecay(Particle::Decay(42, 0.011,  vector<int>{-11,12,-313,111}));
   part->AddDecay(Particle::Decay(42, 0.011,  vector<int>{-11,12,-323,211}));
   part->AddDecay(Particle::Decay(42, 0.011,  vector<int>{-13,14,-313,111}));
   part->AddDecay(Particle::Decay(0, 0.009,  vector<int>{-321,211,211,111}));
   part->AddDecay(Particle::Decay(0, 0.008,  vector<int>{-321,213,211}));
   part->AddDecay(Particle::Decay(0, 0.0073,  vector<int>{-311,321}));
   part->AddDecay(Particle::Decay(0, 0.0066,  vector<int>{221,211}));
   part->AddDecay(Particle::Decay(0, 0.006,  vector<int>{333,211}));
   part->AddDecay(Particle::Decay(0, 0.0057,  vector<int>{-313,211,113}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{333,213}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-11,12,-311,111}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-11,12,-321,211}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-13,14,-311,111}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{221,213}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-13,14,-321,211}));
   part->AddDecay(Particle::Decay(0, 0.0047,  vector<int>{-313,321}));
   part->AddDecay(Particle::Decay(0, 0.0047,  vector<int>{-311,323}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{-321,321,211}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{331,211}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{331,213}));
   part->AddDecay(Particle::Decay(0, 0.0028,  vector<int>{113,211,211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.0022,  vector<int>{211,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-313,211,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.0019,  vector<int>{-321,113,211,211,111}));
   part->AddDecay(Particle::Decay(0, 0.0015,  vector<int>{211,211,211,-211,-211}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{111,211}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{-11,12,221}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{-11,12,331}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{-11,12,113}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{-11,12,223}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{223,211}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{223,213}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{-11,12,111}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-321,211,211,211,-211}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{-13,14,111}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{-13,14,221}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-311,113,211,211,-211}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{-13,14,331}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{-13,14,113}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{-13,14,223}));
   part->AddDecay(Particle::Decay(0, 0.0006,  vector<int>{113,213}));
   part->AddDecay(Particle::Decay(0, 0.0006,  vector<int>{111,213}));
   part->AddDecay(Particle::Decay(0, 0.0006,  vector<int>{113,211}));

   // Creating D*+
   new Particle("D*+", 413, 1, "CharmedMeson", 100, 1, 2.01027, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(413));
   part->AddDecay(Particle::Decay(3, 0.683,  vector<int>{421,211}));
   part->AddDecay(Particle::Decay(3, 0.306,  vector<int>{411,111}));
   part->AddDecay(Particle::Decay(0, 0.011,  vector<int>{411,22}));

   // Creating D*_2+
   new Particle("D*_2+", 415, 1, "CharmedMeson", 100, 1, 2.4601, 0.023, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(415));
   part->AddDecay(Particle::Decay(0, 0.3,  vector<int>{421,211}));
   part->AddDecay(Particle::Decay(0, 0.16,  vector<int>{423,211}));
   part->AddDecay(Particle::Decay(0, 0.15,  vector<int>{411,111}));
   part->AddDecay(Particle::Decay(0, 0.13,  vector<int>{423,211,111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{413,111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{421,211,111}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{413,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.04,  vector<int>{411,211,-211}));

   // Creating D0
   new Particle("D0", 421, 1, "CharmedMeson", 100, 0, 1.86484, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(421));
   part->AddDecay(Particle::Decay(0, 0.0923,  vector<int>{-321,211,111,111}));
   part->AddDecay(Particle::Decay(0, 0.074,  vector<int>{-321,20213}));
   part->AddDecay(Particle::Decay(0, 0.073,  vector<int>{-321,213}));
   part->AddDecay(Particle::Decay(0, 0.067,  vector<int>{-311,211,-211,111,111}));
   part->AddDecay(Particle::Decay(0, 0.062,  vector<int>{-323,213}));
   part->AddDecay(Particle::Decay(0, 0.0511,  vector<int>{-311,113,111,111,111}));
   part->AddDecay(Particle::Decay(0, 0.045,  vector<int>{-323,211}));
   part->AddDecay(Particle::Decay(0, 0.0365,  vector<int>{-321,211}));
   part->AddDecay(Particle::Decay(42, 0.034,  vector<int>{-11,12,-321}));
   part->AddDecay(Particle::Decay(42, 0.034,  vector<int>{-13,14,-321}));
   part->AddDecay(Particle::Decay(42, 0.027,  vector<int>{-11,12,-323}));
   part->AddDecay(Particle::Decay(42, 0.027,  vector<int>{-13,14,-323}));
   part->AddDecay(Particle::Decay(0, 0.025,  vector<int>{-311,223}));
   part->AddDecay(Particle::Decay(0, 0.024,  vector<int>{-321,211,211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.022,  vector<int>{-311,211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.021,  vector<int>{-313,111}));
   part->AddDecay(Particle::Decay(0, 0.021,  vector<int>{-313,221}));
   part->AddDecay(Particle::Decay(0, 0.021,  vector<int>{-311,111}));
   part->AddDecay(Particle::Decay(0, 0.018,  vector<int>{-311,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.018,  vector<int>{-321,211,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.017,  vector<int>{211,211,-211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.016,  vector<int>{-313,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{-313,113}));
   part->AddDecay(Particle::Decay(0, 0.011,  vector<int>{-321,211,111}));
   part->AddDecay(Particle::Decay(0, 0.0109,  vector<int>{-10323,211}));
   part->AddDecay(Particle::Decay(0, 0.009,  vector<int>{-311,321,-321,111}));
   part->AddDecay(Particle::Decay(0, 0.0088,  vector<int>{-311,333}));
   part->AddDecay(Particle::Decay(0, 0.0085,  vector<int>{-311,211,211,-211,-211}));
   part->AddDecay(Particle::Decay(0, 0.0077,  vector<int>{-313,211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.0075,  vector<int>{211,211,-211,-211}));
   part->AddDecay(Particle::Decay(0, 0.0063,  vector<int>{-321,211,113}));
   part->AddDecay(Particle::Decay(0, 0.0061,  vector<int>{-311,113}));
   part->AddDecay(Particle::Decay(0, 0.0052,  vector<int>{-321,321,-311}));
   part->AddDecay(Particle::Decay(0, 0.0041,  vector<int>{-321,321}));
   part->AddDecay(Particle::Decay(42, 0.004,  vector<int>{-13,14,-323,111}));
   part->AddDecay(Particle::Decay(42, 0.004,  vector<int>{-11,12,-313,-211}));
   part->AddDecay(Particle::Decay(42, 0.004,  vector<int>{-11,12,-323,111}));
   part->AddDecay(Particle::Decay(42, 0.004,  vector<int>{-13,14,-313,-211}));
   part->AddDecay(Particle::Decay(0, 0.0036,  vector<int>{-313,321,-211}));
   part->AddDecay(Particle::Decay(0, 0.0035,  vector<int>{-321,323}));
   part->AddDecay(Particle::Decay(0, 0.0034,  vector<int>{-321,311,211}));
   part->AddDecay(Particle::Decay(0, 0.0028,  vector<int>{321,-321,211,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.0027,  vector<int>{-313,313}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{-13,14,-311,-211}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{-13,14,-321,111}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{-11,12,-211}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{-11,12,-213}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-323,321}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{-13,14,-211}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{-13,14,-213}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{-11,12,-311,-211}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{-11,12,-321,111}));
   part->AddDecay(Particle::Decay(0, 0.0018,  vector<int>{333,113}));
   part->AddDecay(Particle::Decay(0, 0.0016,  vector<int>{111,111}));
   part->AddDecay(Particle::Decay(0, 0.0016,  vector<int>{211,-211}));
   part->AddDecay(Particle::Decay(0, 0.0011,  vector<int>{-311,311}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-313,311}));
   part->AddDecay(Particle::Decay(0, 0.0009,  vector<int>{310,310,310}));
   part->AddDecay(Particle::Decay(0, 0.0006,  vector<int>{333,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.0004,  vector<int>{113,211,211,-211,-211}));

   // Creating D*0
   new Particle("D*0", 423, 1, "CharmedMeson", 100, 0, 2.00697, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(423));
   part->AddDecay(Particle::Decay(3, 0.619,  vector<int>{421,111}));
   part->AddDecay(Particle::Decay(0, 0.381,  vector<int>{421,22}));

   // Creating D*_20
   new Particle("D*_20", 425, 1, "CharmedMeson", 100, 0, 2.4611, 0.023, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(425));
   part->AddDecay(Particle::Decay(0, 0.3,  vector<int>{411,-211}));
   part->AddDecay(Particle::Decay(0, 0.16,  vector<int>{413,-211}));
   part->AddDecay(Particle::Decay(0, 0.15,  vector<int>{421,111}));
   part->AddDecay(Particle::Decay(0, 0.13,  vector<int>{413,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{423,111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{411,-211,111}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{423,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.04,  vector<int>{421,211,-211}));

   // Creating D_s+
   new Particle("D_s+", 431, 1, "CharmedMeson", 100, 1, 1.9685, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(431));
   part->AddDecay(Particle::Decay(13, 0.25,  vector<int>{2,-1,3,-3}));
   part->AddDecay(Particle::Decay(13, 0.0952,  vector<int>{2,-1}));
   part->AddDecay(Particle::Decay(0, 0.095,  vector<int>{331,213}));
   part->AddDecay(Particle::Decay(0, 0.079,  vector<int>{221,213}));
   part->AddDecay(Particle::Decay(0, 0.052,  vector<int>{333,213}));
   part->AddDecay(Particle::Decay(0, 0.05,  vector<int>{323,-313}));
   part->AddDecay(Particle::Decay(0, 0.037,  vector<int>{331,211}));
   part->AddDecay(Particle::Decay(0, 0.033,  vector<int>{323,-311}));
   part->AddDecay(Particle::Decay(42, 0.03,  vector<int>{-11,12,333}));
   part->AddDecay(Particle::Decay(42, 0.03,  vector<int>{-13,14,333}));
   part->AddDecay(Particle::Decay(0, 0.028,  vector<int>{333,211}));
   part->AddDecay(Particle::Decay(0, 0.028,  vector<int>{321,-311}));
   part->AddDecay(Particle::Decay(0, 0.026,  vector<int>{321,-313}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{-11,12,331}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{-11,12,221}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{-13,14,221}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{-13,14,331}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{221,211}));
   part->AddDecay(Particle::Decay(0, 0.01,  vector<int>{-15,16}));
   part->AddDecay(Particle::Decay(0, 0.01,  vector<int>{2212,-2112}));
   part->AddDecay(Particle::Decay(0, 0.0078,  vector<int>{10221,211}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-13,14,321,-321}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-13,14,311,-311}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{221,321}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{331,321}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{333,321}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{221,323}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-11,12,311,-311}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-11,12,321,-321}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{213,113}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{211,111}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{213,111}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{211,113}));

   // Creating D*_s+
   new Particle("D*_s+", 433, 1, "CharmedMeson", 100, 1, 2.1123, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(433));
   part->AddDecay(Particle::Decay(0, 0.94,  vector<int>{431,22}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{431,111}));

   // Creating D*_2s+
   new Particle("D*_2s+", 435, 1, "CharmedMeson", 100, 1, 2.5726, 0.015, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(435));
   part->AddDecay(Particle::Decay(0, 0.4,  vector<int>{421,321}));
   part->AddDecay(Particle::Decay(0, 0.4,  vector<int>{411,311}));
   part->AddDecay(Particle::Decay(0, 0.1,  vector<int>{423,321}));
   part->AddDecay(Particle::Decay(0, 0.1,  vector<int>{413,311}));

   // Creating J/psi_di
   new Particle("J/psi_di", 440, 0, "CharmedMeson", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating eta_c
   new Particle("eta_c", 441, 0, "CharmedMeson", 100, 0, 2.9803, 0.0013, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(441));
   part->AddDecay(Particle::Decay(12, 1,  vector<int>{82,-82}));

   // Creating J/psi
   new Particle("J/psi", 443, 0, "Meson", 100, 0, 3.09692, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(443));
   part->AddDecay(Particle::Decay(12, 0.8797,  vector<int>{82,-82}));
   part->AddDecay(Particle::Decay(0, 0.0602,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(0, 0.0601,  vector<int>{13,-13}));

   // Creating chi_2c
   new Particle("chi_2c", 445, 0, "CharmedMeson", 100, 0, 3.5562, 0.002, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(445));
   part->AddDecay(Particle::Decay(12, 0.865,  vector<int>{82,-82}));
   part->AddDecay(Particle::Decay(0, 0.135,  vector<int>{443,22}));

   // Creating B0
   new Particle("B0", 511, 1, "B-Meson", 100, 0, 5.27953, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(511));
   part->AddDecay(Particle::Decay(48, 0.4291,  vector<int>{2,-1,-4,1}));
   part->AddDecay(Particle::Decay(13, 0.08,  vector<int>{2,-4,-1,1}));
   part->AddDecay(Particle::Decay(13, 0.07,  vector<int>{4,-3,-4,1}));
   part->AddDecay(Particle::Decay(42, 0.055,  vector<int>{14,-13,-413}));
   part->AddDecay(Particle::Decay(42, 0.055,  vector<int>{12,-11,-413}));
   part->AddDecay(Particle::Decay(42, 0.03,  vector<int>{16,-15,-413}));
   part->AddDecay(Particle::Decay(0, 0.025,  vector<int>{-413,433}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{12,-11,-411}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{14,-13,-411}));
   part->AddDecay(Particle::Decay(13, 0.02,  vector<int>{4,-4,-3,1}));
   part->AddDecay(Particle::Decay(0, 0.0185,  vector<int>{-411,433}));
   part->AddDecay(Particle::Decay(0, 0.018,  vector<int>{-413,20213}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{-411,431}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,1}));
   part->AddDecay(Particle::Decay(0, 0.0135,  vector<int>{-413,431}));
   part->AddDecay(Particle::Decay(42, 0.012,  vector<int>{14,-13,-415}));
   part->AddDecay(Particle::Decay(42, 0.012,  vector<int>{12,-11,-415}));
   part->AddDecay(Particle::Decay(0, 0.011,  vector<int>{-411,213}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{16,-15,-411}));
   part->AddDecay(Particle::Decay(0, 0.009,  vector<int>{-413,213}));
   part->AddDecay(Particle::Decay(42, 0.008,  vector<int>{14,-13,-20413}));
   part->AddDecay(Particle::Decay(42, 0.008,  vector<int>{12,-11,-20413}));
   part->AddDecay(Particle::Decay(0, 0.0055,  vector<int>{-411,20213}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{12,-11,-10411}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{14,-13,-10413}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{14,-13,-10411}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{12,-11,-10413}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,1}));
   part->AddDecay(Particle::Decay(0, 0.0042,  vector<int>{-413,211}));
   part->AddDecay(Particle::Decay(0, 0.0035,  vector<int>{-411,211}));
   part->AddDecay(Particle::Decay(0, 0.0025,  vector<int>{20443,313}));
   part->AddDecay(Particle::Decay(0, 0.0019,  vector<int>{20443,311}));
   part->AddDecay(Particle::Decay(0, 0.0014,  vector<int>{443,313}));
   part->AddDecay(Particle::Decay(0, 0.0008,  vector<int>{443,311}));
   part->AddDecay(Particle::Decay(0, 0.0007,  vector<int>{441,313}));
   part->AddDecay(Particle::Decay(0, 0.0004,  vector<int>{441,311}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
