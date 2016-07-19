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

   // Creating k3_star(1780)-_bar
   new Particle("k3_star(1780)-_bar", -327, 0, "Unknown", 100, -1, 1.776, 0.159, 100, 100, 0, 100, 1);

   // Creating K*_2-
   new Particle("K*_2-", -325, 0, "Meson", 100, -1, 1.4256, 0.098, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-325));
   part->AddDecay(Particle::Decay(0, 0.332,  vector<int>{-311,-211}));
   part->AddDecay(Particle::Decay(0, 0.168,  vector<int>{-313,-211}));
   part->AddDecay(Particle::Decay(0, 0.166,  vector<int>{-321,-111}));
   part->AddDecay(Particle::Decay(0, 0.086,  vector<int>{-313,-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.084,  vector<int>{-323,-111}));
   part->AddDecay(Particle::Decay(0, 0.059,  vector<int>{-311,-213}));
   part->AddDecay(Particle::Decay(0, 0.043,  vector<int>{-323,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.029,  vector<int>{-321,-113}));
   part->AddDecay(Particle::Decay(0, 0.029,  vector<int>{-321,-223}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-321,-221}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-321,-22}));

   // Creating K*-
   new Particle("K*-", -323, 0, "Meson", 100, -1, 0.89166, 0.0498, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-323));
   part->AddDecay(Particle::Decay(3, 0.666,  vector<int>{-311,-211}));
   part->AddDecay(Particle::Decay(3, 0.333,  vector<int>{-321,-111}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-321,-22}));

   // Creating K-
   new Particle("K-", -321, 0, "Meson", 100, -1, 0.493677, 5.31674e-17, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-321));
   part->AddDecay(Particle::Decay(0, 0.6352,  vector<int>{13,-14}));
   part->AddDecay(Particle::Decay(0, 0.2116,  vector<int>{-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.0559,  vector<int>{-211,-211,211}));
   part->AddDecay(Particle::Decay(42, 0.0482,  vector<int>{-12,11,-111}));
   part->AddDecay(Particle::Decay(42, 0.0318,  vector<int>{-14,13,-111}));
   part->AddDecay(Particle::Decay(0, 0.0173,  vector<int>{-211,-111,-111}));

   // Creating k3_star(1780)0_bar
   new Particle("k3_star(1780)0_bar", -317, 0, "Unknown", 100, 0, 1.776, 0.159, 100, 100, 0, 100, 1);

   // Creating K*_20_bar
   new Particle("K*_20_bar", -315, 0, "Meson", 100, 0, 1.4324, 0.109, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-315));
   part->AddDecay(Particle::Decay(0, 0.333,  vector<int>{-321,211}));
   part->AddDecay(Particle::Decay(0, 0.168,  vector<int>{-323,211}));
   part->AddDecay(Particle::Decay(0, 0.166,  vector<int>{-311,-111}));
   part->AddDecay(Particle::Decay(0, 0.087,  vector<int>{-323,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.084,  vector<int>{-313,-111}));
   part->AddDecay(Particle::Decay(0, 0.059,  vector<int>{-321,213}));
   part->AddDecay(Particle::Decay(0, 0.043,  vector<int>{-313,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.029,  vector<int>{-311,-113}));
   part->AddDecay(Particle::Decay(0, 0.029,  vector<int>{-311,-223}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-311,-221}));

   // Creating K*0_bar
   new Particle("K*0_bar", -313, 0, "Meson", 100, 0, 0.896, 0.0505, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-313));
   part->AddDecay(Particle::Decay(3, 0.665,  vector<int>{-321,211}));
   part->AddDecay(Particle::Decay(3, 0.333,  vector<int>{-311,-111}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{-311,-22}));

   // Creating K0_bar
   new Particle("K0_bar", -311, 0, "Meson", 100, 0, 0.497614, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-311));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-130}));
   part->AddDecay(Particle::Decay(0, 0.5,  vector<int>{-310}));

   // Creating omega3(1670)_bar
   new Particle("omega3(1670)_bar", -227, 0, "Unknown", 100, 0, 1.667, 0.168, 100, 100, 0, 100, 1);

   // Creating rho3(1690)-_bar
   new Particle("rho3(1690)-_bar", -217, 0, "Unknown", 100, -1, 1.6888, 0.161, 100, 100, 0, 100, 1);

   // Creating a_2-
   new Particle("a_2-", -215, 0, "Meson", 100, -1, 1.3183, 0.107, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-215));
   part->AddDecay(Particle::Decay(0, 0.34725,  vector<int>{-213,-111}));
   part->AddDecay(Particle::Decay(0, 0.34725,  vector<int>{-113,-211}));
   part->AddDecay(Particle::Decay(0, 0.144,  vector<int>{-221,-211}));
   part->AddDecay(Particle::Decay(0, 0.104,  vector<int>{-223,-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.049,  vector<int>{-321,311}));
   part->AddDecay(Particle::Decay(0, 0.0057,  vector<int>{-331,-211}));
   part->AddDecay(Particle::Decay(0, 0.0028,  vector<int>{-211,-22}));

   // Creating rho-
   new Particle("rho-", -213, 0, "Meson", 100, -1, 0.77549, 0.149, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-213));
   part->AddDecay(Particle::Decay(3, 0.99955,  vector<int>{-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.00045,  vector<int>{-211,-22}));

   // Creating pi-
   new Particle("pi-", -211, 0, "Meson", 100, -1, 0.13957, 2.52837e-17, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-211));
   part->AddDecay(Particle::Decay(0, 0.999877,  vector<int>{13,-14}));
   part->AddDecay(Particle::Decay(0, 0.000123,  vector<int>{11,-12}));

   // Creating pi_diffr-
   new Particle("pi_diffr-", -210, 0, "Meson", 100, -1, 0, 0, 100, 100, 1, 100, 1);

   // Creating rho3(1690)0_bar
   new Particle("rho3(1690)0_bar", -117, 0, "Unknown", 100, 0, 1.6888, 0.161, 100, 100, 0, 100, 1);

   // Creating b-hadron_bar
   new Particle("b-hadron_bar", -85, 0, "Generator", 100, 0.333333, 5, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-85));
   part->AddDecay(Particle::Decay(42, 0.5,  vector<int>{2,-1,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.14,  vector<int>{4,-3,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{12,-11,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.105,  vector<int>{14,-13,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{2,-4,-1,-81}));
   part->AddDecay(Particle::Decay(42, 0.04,  vector<int>{16,-15,-4,-81}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{2,-1,-2,-81}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{4,-4,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{4,-3,-2,-81}));

   // Creating c-hadron_bar
   new Particle("c-hadron_bar", -84, 0, "Generator", 100, -0.666667, 2, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-84));
   part->AddDecay(Particle::Decay(11, 0.76,  vector<int>{-2,1,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{13,-14,-3,-81}));
   part->AddDecay(Particle::Decay(42, 0.08,  vector<int>{11,-12,-3,-81}));
   part->AddDecay(Particle::Decay(11, 0.08,  vector<int>{-2,3,-3,-81}));

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
   part->AddDecay(Particle::Decay(32, 0.325914,  vector<int>{1,-2}));
   part->AddDecay(Particle::Decay(32, 0.32532,  vector<int>{3,-4}));
   part->AddDecay(Particle::Decay(32, 0.314118,  vector<int>{5,-6}));
   part->AddDecay(Particle::Decay(32, 0.016736,  vector<int>{3,-2}));
   part->AddDecay(Particle::Decay(32, 0.016735,  vector<int>{1,-4}));
   part->AddDecay(Particle::Decay(32, 0.000603001,  vector<int>{5,-4}));
   part->AddDecay(Particle::Decay(32, 0.000554001,  vector<int>{3,-6}));
   part->AddDecay(Particle::Decay(32, 1e-05,  vector<int>{5,-2}));
   part->AddDecay(Particle::Decay(32, 9.00001e-06,  vector<int>{1,-6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{11,-64}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{13,-65}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{15,-66}));

   // Creating H_R--
   new Particle("H_R--", -62, 0, "Unknown", 100, -2, 200, 0.88001, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-62));
   part->AddDecay(Particle::Decay(0, 0.813719,  vector<int>{15,15}));
   part->AddDecay(Particle::Decay(0, 0.0904279,  vector<int>{13,13}));
   part->AddDecay(Particle::Decay(0, 0.0904279,  vector<int>{11,11}));
   part->AddDecay(Particle::Decay(0, 0.001809,  vector<int>{11,13}));
   part->AddDecay(Particle::Decay(0, 0.001808,  vector<int>{13,15}));
   part->AddDecay(Particle::Decay(0, 0.001808,  vector<int>{11,15}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-63,-63}));

   // Creating H_L--
   new Particle("H_L--", -61, 0, "Unknown", 100, -2, 200, 0.88161, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-61));
   part->AddDecay(Particle::Decay(0, 0.812251,  vector<int>{15,15}));
   part->AddDecay(Particle::Decay(0, 0.0902641,  vector<int>{13,13}));
   part->AddDecay(Particle::Decay(0, 0.0902641,  vector<int>{11,11}));
   part->AddDecay(Particle::Decay(0, 0.001806,  vector<int>{-24,-24}));
   part->AddDecay(Particle::Decay(0, 0.001805,  vector<int>{13,15}));
   part->AddDecay(Particle::Decay(0, 0.001805,  vector<int>{11,15}));
   part->AddDecay(Particle::Decay(0, 0.001805,  vector<int>{11,13}));

   // Creating rho_tech-
   new Particle("rho_tech-", -55, 0, "Unknown", 100, -1, 210, 0.64973, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-55));
   part->AddDecay(Particle::Decay(0, 0.474101,  vector<int>{-24,-51}));
   part->AddDecay(Particle::Decay(0, 0.176299,  vector<int>{-52,-23}));
   part->AddDecay(Particle::Decay(0, 0.138845,  vector<int>{-24,-23}));
   part->AddDecay(Particle::Decay(0, 0.109767,  vector<int>{-52,-22}));
   part->AddDecay(Particle::Decay(32, 0.0285839,  vector<int>{1,-2}));
   part->AddDecay(Particle::Decay(32, 0.0285299,  vector<int>{3,-4}));
   part->AddDecay(Particle::Decay(0, 0.00966098,  vector<int>{11,-12}));
   part->AddDecay(Particle::Decay(0, 0.00966098,  vector<int>{13,-14}));
   part->AddDecay(Particle::Decay(0, 0.00965998,  vector<int>{15,-16}));
   part->AddDecay(Particle::Decay(0, 0.00816098,  vector<int>{-24,-53}));
   part->AddDecay(Particle::Decay(32, 0.00373499,  vector<int>{5,-6}));
   part->AddDecay(Particle::Decay(32, 0.001468,  vector<int>{3,-2}));
   part->AddDecay(Particle::Decay(32, 0.001468,  vector<int>{1,-4}));
   part->AddDecay(Particle::Decay(32, 5.29999e-05,  vector<int>{5,-4}));
   part->AddDecay(Particle::Decay(32, 6.99999e-06,  vector<int>{3,-6}));
   part->AddDecay(Particle::Decay(32, 9.99998e-07,  vector<int>{5,-2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{1,-8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{5,-8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-4}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{3,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-52,-51}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{1,-6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{17,-18}));

   // Creating pi_tech-
   new Particle("pi_tech-", -52, 0, "Unknown", 100, -1, 110, 0.0105, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-52));
   part->AddDecay(Particle::Decay(32, 0.90916,  vector<int>{-4,5}));
   part->AddDecay(Particle::Decay(0, 0.048905,  vector<int>{15,-16}));
   part->AddDecay(Particle::Decay(32, 0.041762,  vector<int>{-4,3}));
   part->AddDecay(Particle::Decay(0, 0.000173,  vector<int>{13,-14}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-24,-5,5}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{11,-12}));

   // Creating R0_bar
   new Particle("R0_bar", -40, 0, "Unknown", 100, 0, 5000, 417.465, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-40));
   part->AddDecay(Particle::Decay(32, 0.215134,  vector<int>{-1,3}));
   part->AddDecay(Particle::Decay(32, 0.215134,  vector<int>{-2,4}));
   part->AddDecay(Particle::Decay(32, 0.215133,  vector<int>{-3,5}));
   part->AddDecay(Particle::Decay(32, 0.214738,  vector<int>{-4,6}));
   part->AddDecay(Particle::Decay(0, 0.0699301,  vector<int>{-11,13}));
   part->AddDecay(Particle::Decay(0, 0.0699301,  vector<int>{-13,15}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-5,7}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-6,8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-15,17}));

   // Creating LQ_ue_bar
   new Particle("LQ_ue_bar", -39, 0, "Unknown", 100, 0.333333, 200, 0.39162, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-39));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-2,-11}));

   // Creating H-
   new Particle("H-", -37, 0, "GaugeBoson", 100, -1, 300, 5.75967, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-37));
   part->AddDecay(Particle::Decay(0, 0.929792,  vector<int>{-24,-25}));
   part->AddDecay(Particle::Decay(32, 0.067484,  vector<int>{5,-6}));
   part->AddDecay(Particle::Decay(0, 0.002701,  vector<int>{15,-16}));
   part->AddDecay(Particle::Decay(32, 1.3e-05,  vector<int>{3,-4}));
   part->AddDecay(Particle::Decay(0, 1e-05,  vector<int>{13,-14}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{1,-2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{17,-18}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{11,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,1000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,1000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-1000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-1000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-1000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-1000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-1000016}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000015,-1000016}));

   // Creating W'-
   new Particle("W'-", -34, 0, "GaugeBoson", 100, -1, 500, 16.6708, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-34));
   part->AddDecay(Particle::Decay(32, 0.251276,  vector<int>{1,-2}));
   part->AddDecay(Particle::Decay(32, 0.250816,  vector<int>{3,-4}));
   part->AddDecay(Particle::Decay(32, 0.215459,  vector<int>{5,-6}));
   part->AddDecay(Particle::Decay(0, 0.085262,  vector<int>{11,-12}));
   part->AddDecay(Particle::Decay(0, 0.085262,  vector<int>{13,-14}));
   part->AddDecay(Particle::Decay(0, 0.08526,  vector<int>{15,-16}));
   part->AddDecay(Particle::Decay(32, 0.012903,  vector<int>{3,-2}));
   part->AddDecay(Particle::Decay(32, 0.012903,  vector<int>{1,-4}));
   part->AddDecay(Particle::Decay(32, 0.000465,  vector<int>{5,-4}));
   part->AddDecay(Particle::Decay(32, 0.00038,  vector<int>{3,-6}));
   part->AddDecay(Particle::Decay(32, 8e-06,  vector<int>{5,-2}));
   part->AddDecay(Particle::Decay(32, 6e-06,  vector<int>{1,-6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-4}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{3,-8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{1,-8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{5,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{17,-18}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,-23}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,-22}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,-25}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
