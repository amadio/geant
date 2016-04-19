// This files was autogenerated by geant::Particle::ReadFile

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
   vector<int> daughters;
   Particle *part = nullptr;

   // Creating c-hadron_bar
   new Particle("c-hadron_bar", -84, 0, "Generator", 100, -0.666667, 2, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-84));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(1);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(11, 0.76,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(42, 0.08,  daughters));
   daughters.clear();
   daughters.push_back(-2);
   daughters.push_back(3);
   daughters.push_back(-3);
   daughters.push_back(-81);
   part->AddDecay(Particle::Decay(11, 0.08,  daughters));

   // Creating rndmflav_bar
   new Particle("rndmflav_bar", -82, 0, "Generator", 100, 0, 0, 0, 100, 100, 1, 100, 1);

   // Creating nu_Rtau_bar
   new Particle("nu_Rtau_bar", -66, 0, "Unknown", 100, 0, 750, 0, 100, 100, 1, 100, 1);

   // Creating nu_Rmu_bar
   new Particle("nu_Rmu_bar", -65, 0, "Unknown", 100, 0, 750, 0, 100, 100, 1, 100, 1);

   // Creating nu_Re_bar
   new Particle("nu_Re_bar", -64, 0, "Unknown", 100, 0, 750, 0, 100, 100, 1, 100, 1);

   // Creating W_R-
   new Particle("W_R-", -63, 0, "Unknown", 100, -1, 750, 19.3391, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-63));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0.325914,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.32532,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0.314118,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0.016736,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.016735,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.000603001,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0.000554001,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 1e-05,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 9.00001e-06,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-64);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-65);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(-66);
   part->AddDecay(Particle::Decay(0, 0,  daughters));

   // Creating H_R--
   new Particle("H_R--", -62, 0, "Unknown", 100, -2, 200, 0.88001, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-62));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(15);
   part->AddDecay(Particle::Decay(0, 0.813719,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(13);
   part->AddDecay(Particle::Decay(0, 0.0904279,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(11);
   part->AddDecay(Particle::Decay(0, 0.0904279,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(13);
   part->AddDecay(Particle::Decay(0, 0.001809,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(15);
   part->AddDecay(Particle::Decay(0, 0.001808,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(15);
   part->AddDecay(Particle::Decay(0, 0.001808,  daughters));
   daughters.clear();
   daughters.push_back(-63);
   daughters.push_back(-63);
   part->AddDecay(Particle::Decay(0, 0,  daughters));

   // Creating H_L--
   new Particle("H_L--", -61, 0, "Unknown", 100, -2, 200, 0.88161, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-61));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(15);
   part->AddDecay(Particle::Decay(0, 0.812251,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(13);
   part->AddDecay(Particle::Decay(0, 0.0902641,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(11);
   part->AddDecay(Particle::Decay(0, 0.0902641,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-24);
   part->AddDecay(Particle::Decay(0, 0.001806,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(15);
   part->AddDecay(Particle::Decay(0, 0.001805,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(15);
   part->AddDecay(Particle::Decay(0, 0.001805,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(13);
   part->AddDecay(Particle::Decay(0, 0.001805,  daughters));

   // Creating rho_tech-
   new Particle("rho_tech-", -55, 0, "Unknown", 100, -1, 210, 0.64973, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-55));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-51);
   part->AddDecay(Particle::Decay(0, 0.474101,  daughters));
   daughters.clear();
   daughters.push_back(-52);
   daughters.push_back(-23);
   part->AddDecay(Particle::Decay(0, 0.176299,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-23);
   part->AddDecay(Particle::Decay(0, 0.138845,  daughters));
   daughters.clear();
   daughters.push_back(-52);
   daughters.push_back(-22);
   part->AddDecay(Particle::Decay(0, 0.109767,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0.0285839,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.0285299,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(0, 0.00966098,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(0, 0.00966098,  daughters));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(0, 0.00965998,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-53);
   part->AddDecay(Particle::Decay(0, 0.00816098,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0.00373499,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0.001468,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.001468,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 5.29999e-05,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 6.99999e-06,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 9.99998e-07,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(-52);
   daughters.push_back(-51);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(17);
   daughters.push_back(-18);
   part->AddDecay(Particle::Decay(0, 0,  daughters));

   // Creating pi_tech-
   new Particle("pi_tech-", -52, 0, "Unknown", 100, -1, 110, 0.0105, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-52));
   daughters.clear();
   daughters.push_back(-4);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(32, 0.90916,  daughters));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(0, 0.048905,  daughters));
   daughters.clear();
   daughters.push_back(-4);
   daughters.push_back(3);
   part->AddDecay(Particle::Decay(32, 0.041762,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(0, 0.000173,  daughters));
   daughters.clear();
   daughters.push_back(-24);
   daughters.push_back(-5);
   daughters.push_back(5);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif