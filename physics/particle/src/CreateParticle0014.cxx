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
void CreateParticle0014() {

   // Creating nu_Rmu
   new Particle("nu_Rmu", 65, 1, "Unknown", 100, 0, 750, 0, -100, -1, -100, -1, -1);

   // Creating nu_Rtau
   new Particle("nu_Rtau", 66, 1, "Unknown", 100, 0, 750, 0, -100, -1, -100, -1, -1);

   // Creating specflav
   new Particle("specflav", 81, 0, "Generator", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating rndmflav
   new Particle("rndmflav", 82, 1, "Generator", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating phasespa
   new Particle("phasespa", 83, 0, "Generator", 100, 0, 1, 0, -100, -1, -100, -1, -1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(83));
   part->AddDecay(Particle::Decay(12, 1,  vector<int>{82,-82}));

   // Creating c-hadron
   new Particle("c-hadron", 84, 1, "Generator", 100, 0.666667, 2, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(84));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{2,-1,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-13,14,3,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-11,12,3,81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{2,-3,3,81}));

   // Creating b-hadron
   new Particle("b-hadron", 85, 1, "Generator", 100, -0.333333, 5, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(85));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{-2,1,4,81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{-4,3,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-12,11,4,81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{-14,13,4,81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{-2,4,1,81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{-16,15,4,81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-4,4,3,81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,81}));

   // Creating cluster
   new Particle("cluster", 91, 0, "Generator", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating string
   new Particle("string", 92, 0, "Generator", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating indep.
   new Particle("indep.", 93, 0, "Generator", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating CMshower
   new Particle("CMshower", 94, 0, "Generator", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating SPHEaxis
   new Particle("SPHEaxis", 95, 0, "Generator", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating THRUaxis
   new Particle("THRUaxis", 96, 0, "Generator", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating CLUSjet
   new Particle("CLUSjet", 97, 0, "Generator", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating CELLjet
   new Particle("CELLjet", 98, 0, "Generator", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating table
   new Particle("table", 99, 0, "Generator", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating rho_diff0
   new Particle("rho_diff0", 110, 0, "Unknown", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating pi0
   new Particle("pi0", 111, 0, "Meson", 100, 0, 0.134977, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(111));
   part->AddDecay(Particle::Decay(0, 0.988,  vector<int>{22,22}));
   part->AddDecay(Particle::Decay(2, 0.012,  vector<int>{22,11,-11}));

   // Creating rho0
   new Particle("rho0", 113, 0, "Meson", 100, 0, 0.77549, 0.151, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(113));
   part->AddDecay(Particle::Decay(3, 0.998739,  vector<int>{211,-211}));
   part->AddDecay(Particle::Decay(0, 0.00079,  vector<int>{111,22}));
   part->AddDecay(Particle::Decay(0, 0.00038,  vector<int>{221,22}));
   part->AddDecay(Particle::Decay(0, 4.6e-05,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(0, 4.5e-05,  vector<int>{11,-11}));

   // Creating a_20
   new Particle("a_20", 115, 0, "Meson", 100, 0, 1.3183, 0.107, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(115));
   part->AddDecay(Particle::Decay(0, 0.34725,  vector<int>{213,-211}));
   part->AddDecay(Particle::Decay(0, 0.34725,  vector<int>{-213,211}));
   part->AddDecay(Particle::Decay(0, 0.144,  vector<int>{221,111}));
   part->AddDecay(Particle::Decay(0, 0.104,  vector<int>{223,211,-211}));
   part->AddDecay(Particle::Decay(0, 0.0245,  vector<int>{321,-321}));
   part->AddDecay(Particle::Decay(0, 0.01225,  vector<int>{130,130}));
   part->AddDecay(Particle::Decay(0, 0.01225,  vector<int>{310,310}));
   part->AddDecay(Particle::Decay(0, 0.0057,  vector<int>{331,111}));
   part->AddDecay(Particle::Decay(0, 0.0028,  vector<int>{111,22}));

   // Creating rho3(1690)0
   new Particle("rho3(1690)0", 117, 1, "Unknown", 100, 0, 1.6888, 0.161, -100, 0, -100, -1, -1);

   // Creating K_L0
   new Particle("K_L0", 130, 0, "Meson", 100, 0, 0.497614, 0, -100, -1, -100, -1, -1);
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
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
