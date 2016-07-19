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
void CreateParticle0012() {

   // Creating W-
   new Particle("W-", -24, 0, "GaugeBoson", 100, -1, 80.398, 2.07002, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-24));
   part->AddDecay(Particle::Decay(32, 0.321502,  vector<int>{1,-2}));
   part->AddDecay(Particle::Decay(32, 0.320778,  vector<int>{3,-4}));
   part->AddDecay(Particle::Decay(0, 0.108062,  vector<int>{11,-12}));
   part->AddDecay(Particle::Decay(0, 0.108062,  vector<int>{13,-14}));
   part->AddDecay(Particle::Decay(0, 0.107983,  vector<int>{15,-16}));
   part->AddDecay(Particle::Decay(32, 0.016509,  vector<int>{3,-2}));
   part->AddDecay(Particle::Decay(32, 0.016502,  vector<int>{1,-4}));
   part->AddDecay(Particle::Decay(32, 0.000591001,  vector<int>{5,-4}));
   part->AddDecay(Particle::Decay(32, 1e-05,  vector<int>{5,-2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{1,-8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{5,-6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{5,-8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-4}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{3,-6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{3,-8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{1,-6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{17,-18}));

   // Creating nu'_tau_bar
   new Particle("nu'_tau_bar", -18, 0, "Lepton", 100, 0, 0, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-18));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-18}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,-17}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-37,-17}));

   // Creating tau'+
   new Particle("tau'+", -17, 0, "Lepton", 100, 1, 400, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-17));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-17}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-17}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,-18}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-17}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{37,-18}));

   // Creating nu_tau_bar
   new Particle("nu_tau_bar", -16, 0, "Lepton", 100, 0, 0, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-16));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-16}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,-15}));

   // Creating tau+
   new Particle("tau+", -15, 0, "Lepton", 100, 1, 1.77684, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-15));
   part->AddDecay(Particle::Decay(0, 0.2494,  vector<int>{-16,213}));
   part->AddDecay(Particle::Decay(42, 0.1783,  vector<int>{12,-11,-16}));
   part->AddDecay(Particle::Decay(42, 0.1735,  vector<int>{14,-13,-16}));
   part->AddDecay(Particle::Decay(0, 0.1131,  vector<int>{-16,211}));
   part->AddDecay(Particle::Decay(41, 0.09,  vector<int>{-16,213,-111}));
   part->AddDecay(Particle::Decay(41, 0.08,  vector<int>{-16,211,-113}));
   part->AddDecay(Particle::Decay(41, 0.0191,  vector<int>{-16,211,-223}));
   part->AddDecay(Particle::Decay(41, 0.0133,  vector<int>{-16,211,-113,-111}));
   part->AddDecay(Particle::Decay(0, 0.012,  vector<int>{-16,323}));
   part->AddDecay(Particle::Decay(41, 0.011,  vector<int>{-16,211,-211,211}));
   part->AddDecay(Particle::Decay(41, 0.01,  vector<int>{-16,213,-111,-111}));
   part->AddDecay(Particle::Decay(0, 0.0071,  vector<int>{-16,321}));
   part->AddDecay(Particle::Decay(41, 0.0067,  vector<int>{-16,213,-211,211}));
   part->AddDecay(Particle::Decay(41, 0.005,  vector<int>{-16,213,-113}));
   part->AddDecay(Particle::Decay(41, 0.0035,  vector<int>{-16,213,-223}));
   part->AddDecay(Particle::Decay(41, 0.0034,  vector<int>{-16,321,-321,211}));
   part->AddDecay(Particle::Decay(41, 0.003,  vector<int>{-16,211,-111}));
   part->AddDecay(Particle::Decay(41, 0.0027,  vector<int>{-16,211,-111,-111}));
   part->AddDecay(Particle::Decay(41, 0.00205,  vector<int>{-16,211,-310,-111}));
   part->AddDecay(Particle::Decay(41, 0.00205,  vector<int>{-16,211,-130,-111}));
   part->AddDecay(Particle::Decay(41, 0.0015,  vector<int>{-16,213,-221}));
   part->AddDecay(Particle::Decay(41, 0.0014,  vector<int>{-16,211,-111,-111,-111}));
   part->AddDecay(Particle::Decay(41, 0.0012,  vector<int>{-16,213,-111,-111,-111}));
   part->AddDecay(Particle::Decay(41, 0.0011,  vector<int>{-16,213,-113,-111,-111}));
   part->AddDecay(Particle::Decay(41, 0.00078,  vector<int>{-16,321,-310}));
   part->AddDecay(Particle::Decay(41, 0.00078,  vector<int>{-16,321,-130}));
   part->AddDecay(Particle::Decay(41, 0.00075,  vector<int>{-16,211,-113,-113}));
   part->AddDecay(Particle::Decay(41, 0.00075,  vector<int>{-16,323,-111}));
   part->AddDecay(Particle::Decay(41, 0.00069,  vector<int>{-16,321,-310,-111}));
   part->AddDecay(Particle::Decay(41, 0.00069,  vector<int>{-16,321,-130,-111}));
   part->AddDecay(Particle::Decay(41, 0.0006,  vector<int>{-16,211,-223,-111}));
   part->AddDecay(Particle::Decay(41, 0.00051,  vector<int>{-16,211,-310,-130}));
   part->AddDecay(Particle::Decay(41, 0.0005,  vector<int>{-16,211,-211,211,-111}));
   part->AddDecay(Particle::Decay(41, 0.0004,  vector<int>{-16,323,-111,-111}));
   part->AddDecay(Particle::Decay(41, 0.0004,  vector<int>{-16,321,-111}));
   part->AddDecay(Particle::Decay(41, 0.00025,  vector<int>{-16,211,-130}));
   part->AddDecay(Particle::Decay(41, 0.00025,  vector<int>{-16,211,-310,-310}));
   part->AddDecay(Particle::Decay(41, 0.00025,  vector<int>{-16,211,-310}));
   part->AddDecay(Particle::Decay(41, 0.00025,  vector<int>{-16,211,-130,-130}));
   part->AddDecay(Particle::Decay(41, 0.00022,  vector<int>{-16,211,-113,-113,-111}));
   part->AddDecay(Particle::Decay(41, 0.00021,  vector<int>{-16,211,-221,-111}));
   part->AddDecay(Particle::Decay(41, 0.0002,  vector<int>{-16,211,213,-211,-111}));
   part->AddDecay(Particle::Decay(41, 0.0002,  vector<int>{-16,211,-113,-111,-111}));
   part->AddDecay(Particle::Decay(41, 0.0002,  vector<int>{-16,213,-113,-111}));
   part->AddDecay(Particle::Decay(41, 0.0002,  vector<int>{-16,211,-213,213}));
   part->AddDecay(Particle::Decay(41, 0.0002,  vector<int>{-16,211,-213,211,-111}));
   part->AddDecay(Particle::Decay(41, 0.0001,  vector<int>{-16,321,-111,-111,-111}));
   part->AddDecay(Particle::Decay(41, 0.0001,  vector<int>{-16,211,-221,-221}));
   part->AddDecay(Particle::Decay(41, 6e-05,  vector<int>{-16,323,-111,-111}));
   part->AddDecay(Particle::Decay(41, 6e-05,  vector<int>{-16,211,-221}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-15}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-15}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-16}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-15}));

   // Creating nu_mu_bar
   new Particle("nu_mu_bar", -14, 0, "Lepton", 100, 0, 0, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-14));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-14}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,-13}));

   // Creating mu+
   new Particle("mu+", -13, 0, "Lepton", 100, 1, 0.105658, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-13));
   part->AddDecay(Particle::Decay(42, 1,  vector<int>{12,-11,-14}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-13}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-13}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-14}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-13}));

   // Creating nu_e_bar
   new Particle("nu_e_bar", -12, 0, "Lepton", 100, 0, 0, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-12));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-12}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,-11}));

   // Creating e+
   new Particle("e+", -11, 0, "Lepton", 100, 1, 0.000510999, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-11));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-11}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-11}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-12}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-11}));

   // Creating t'_bar
   new Particle("t'_bar", -8, 0, "Quark", 100, -0.666667, 171.2, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-8));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-21,-8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,-1}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,-3}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,-5}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,-7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-37,-5}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-37,-7}));

   // Creating b'_bar
   new Particle("b'_bar", -7, 0, "Quark", 100, 0.333333, 468, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-7));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-21,-7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-7}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,-2}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,-4}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,-6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,-8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-7}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{37,-4}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{37,-6}));

   // Creating t_bar
   new Particle("t_bar", -6, 0, "Quark", 100, -0.666667, 171.2, 1.39883, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-6));
   part->AddDecay(Particle::Decay(0, 0.998205,  vector<int>{-24,-5}));
   part->AddDecay(Particle::Decay(0, 0.001765,  vector<int>{-24,-3}));
   part->AddDecay(Particle::Decay(0, 3e-05,  vector<int>{-24,-1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-21,-6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,-7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-37,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-1000006}));

   // Creating b_bar
   new Particle("b_bar", -5, 0, "Quark", 100, 0.333333, 4.68, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-5));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-21,-5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-5}));

   // Creating c_bar
   new Particle("c_bar", -4, 0, "Quark", 100, -0.666667, 1.27, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-21,-4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,-1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,-3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,-5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,-7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-4}));

   // Creating s_bar
   new Particle("s_bar", -3, 0, "Quark", 100, 0.333333, 0.104, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-21,-3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-3}));

   // Creating u_bar
   new Particle("u_bar", -2, 0, "Quark", 100, -0.666667, 0.0024, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-21,-2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,-1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,-3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,-5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,-7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-2}));

   // Creating d_bar
   new Particle("d_bar", -1, 0, "Quark", 100, 0.333333, 0.0048, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-21,-1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-22,-1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-23,-1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,-8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-25,-1}));

   // Creating Rootino
   new Particle("Rootino", 0, 0, "Unknown", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating d
   new Particle("d", 1, 1, "Quark", 100, -0.333333, 0.0048, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(1));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{21,1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,1}));

   // Creating u
   new Particle("u", 2, 1, "Quark", 100, 0.666667, 0.0024, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(2));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{21,2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,2}));

   // Creating s
   new Particle("s", 3, 1, "Quark", 100, -0.333333, 0.104, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(3));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{21,3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,3}));

   // Creating c
   new Particle("c", 4, 1, "Quark", 100, 0.666667, 1.27, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{21,4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,4}));

   // Creating b
   new Particle("b", 5, 1, "Quark", 100, -0.333333, 4.68, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(5));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{21,5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,5}));

   // Creating t
   new Particle("t", 6, 1, "Quark", 100, 0.666667, 171.2, 1.39883, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(6));
   part->AddDecay(Particle::Decay(0, 0.998205,  vector<int>{24,5}));
   part->AddDecay(Particle::Decay(0, 0.001765,  vector<int>{24,3}));
   part->AddDecay(Particle::Decay(0, 3e-05,  vector<int>{24,1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{21,6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{37,5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000021,1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,1000006}));

   // Creating b'
   new Particle("b'", 7, 1, "Quark", 100, -0.333333, 468, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(7));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{21,7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,7}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,2}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,4}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,7}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-37,4}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-37,6}));

   // Creating t'
   new Particle("t'", 8, 1, "Quark", 100, 0.666667, 171.2, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(8));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{21,8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,1}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,3}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,5}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{37,5}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{37,7}));

   // Creating e-
   new Particle("e-", 11, 1, "Lepton", 100, -1, 0.000510999, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(11));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,11}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,11}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,12}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,11}));

   // Creating nu_e
   new Particle("nu_e", 12, 1, "Lepton", 100, 0, 0, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(12));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,12}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,11}));

   // Creating mu-
   new Particle("mu-", 13, 1, "Lepton", 100, -1, 0.105658, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(13));
   part->AddDecay(Particle::Decay(42, 1,  vector<int>{-12,11,14}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,13}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,13}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,14}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,13}));

   // Creating nu_mu
   new Particle("nu_mu", 14, 1, "Lepton", 100, 0, 0, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(14));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,14}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,13}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
