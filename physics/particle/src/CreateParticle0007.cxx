#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
using std::vector;
namespace geant {
   inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
void CreateParticle0007() {

   // Creating Delta--
   new Particle("Delta--", -2224, 0, "Baryon", 100, -2, 1.232, 0.12, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-2224));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-2212,-211}));

   // Creating delta(1620)++_bar
   new Particle("delta(1620)++_bar", -2222, 0, "Unknown", 100, -2, 1.63, 0.145, 100, 100, 0, 100, 1);

   // Creating delta(1950)+_bar
   new Particle("delta(1950)+_bar", -2218, 0, "Unknown", 100, -1, 1.93, 0.28, 100, 100, 0, 100, 1);

   // Creating N(1675)+_bar
   new Particle("N(1675)+_bar", -2216, 0, "Unknown", 100, -1, 1.675, 0.15, 100, 100, 0, 100, 1);

   // Creating Delta+_bar
   new Particle("Delta+_bar", -2214, 0, "Baryon", 100, -1, 1.232, 0.12, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2214));
   part->AddDecay(Particle::Decay(0, 0.663,  vector<int>{-2212,-111}));
   part->AddDecay(Particle::Decay(0, 0.331,  vector<int>{-2112,-211}));
   part->AddDecay(Particle::Decay(0, 0.006,  vector<int>{-2212,-22}));

   // Creating antiproton
   new Particle("antiproton", -2212, 0, "Baryon", 100, -1, 0.938272, 0, 100, 100, 1, 100, 1);

   // Creating p_diffr+_bar
   new Particle("p_diffr+_bar", -2210, 0, "Unknown", 100, -1, 0, 0, 100, 100, 1, 100, 1);

   // Creating uu_1_bar
   new Particle("uu_1_bar", -2203, 0, "Unknown", 100, -1.33333, 0.0048, 0, 100, 100, 1, 100, 1);

   // Creating N(2190)+_bar
   new Particle("N(2190)+_bar", -2128, 0, "Unknown", 100, -1, 2.19, 0.5, 100, 100, 0, 100, 1);

   // Creating delta(1905)+_bar
   new Particle("delta(1905)+_bar", -2126, 0, "Unknown", 100, -1, 1.89, 0.33, 100, 100, 0, 100, 1);

   // Creating N(1520)+_bar
   new Particle("N(1520)+_bar", -2124, 0, "Unknown", 100, -1, 1.52, 0.115, 100, 100, 0, 100, 1);

   // Creating delta(1620)+_bar
   new Particle("delta(1620)+_bar", -2122, 0, "Unknown", 100, -1, 1.63, 0.145, 100, 100, 0, 100, 1);

   // Creating delta(1950)0_bar
   new Particle("delta(1950)0_bar", -2118, 0, "Unknown", 100, 0, 1.93, 0.28, 100, 100, 0, 100, 1);

   // Creating N(1675)0_bar
   new Particle("N(1675)0_bar", -2116, 0, "Unknown", 100, 0, 1.675, 0.15, 100, 100, 0, 100, 1);

   // Creating Delta0_bar
   new Particle("Delta0_bar", -2114, 0, "Baryon", 100, 0, 1.232, 0.12, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2114));
   part->AddDecay(Particle::Decay(0, 0.663,  vector<int>{-2112,-111}));
   part->AddDecay(Particle::Decay(0, 0.331,  vector<int>{-2212,211}));
   part->AddDecay(Particle::Decay(0, 0.006,  vector<int>{-2112,-22}));

   // Creating antineutron
   new Particle("antineutron", -2112, 0, "Baryon", 100, 0, 0.939565, 0, 100, 100, 1, 100, 1);

   // Creating n_diffr0_bar
   new Particle("n_diffr0_bar", -2110, 0, "Unknown", 100, 0, 0, 0, 100, 100, 1, 100, 1);

   // Creating ud_1_bar
   new Particle("ud_1_bar", -2103, 0, "Unknown", 100, -0.333333, 0.0072, 0, 100, 100, 1, 100, 1);

   // Creating ud_0_bar
   new Particle("ud_0_bar", -2101, 0, "Unknown", 100, -0.333333, 0.0073, 0, 100, 100, 1, 100, 1);

   // Creating N(2190)0_bar
   new Particle("N(2190)0_bar", -1218, 0, "Unknown", 100, 0, 2.19, 0.5, 100, 100, 0, 100, 1);

   // Creating delta(1905)0_bar
   new Particle("delta(1905)0_bar", -1216, 0, "Unknown", 100, 0, 1.89, 0.33, 100, 100, 0, 100, 1);

   // Creating N(1520)0_bar
   new Particle("N(1520)0_bar", -1214, 0, "Unknown", 100, 0, 1.52, 0.115, 100, 100, 0, 100, 1);

   // Creating delta(1620)0_bar
   new Particle("delta(1620)0_bar", -1212, 0, "Unknown", 100, 0, 1.63, 0.145, 100, 100, 0, 100, 1);

   // Creating delta(1950)-_bar
   new Particle("delta(1950)-_bar", -1118, 0, "Unknown", 100, 1, 1.93, 0.28, 100, 100, 0, 100, 1);

   // Creating delta(1905)-_bar
   new Particle("delta(1905)-_bar", -1116, 0, "Unknown", 100, 1, 1.89, 0.33, 100, 100, 0, 100, 1);

   // Creating Delta-_bar
   new Particle("Delta-_bar", -1114, 0, "Unknown", 100, 1, 1.232, 0.12, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1114));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-2112,211}));

   // Creating delta(1620)-_bar
   new Particle("delta(1620)-_bar", -1112, 0, "Unknown", 100, 1, 1.63, 0.145, 100, 100, 0, 100, 1);

   // Creating dd_1_bar
   new Particle("dd_1_bar", -1103, 0, "Unknown", 100, 0.666667, 0.96, 0, 100, 100, 1, 100, 1);

   // Creating B*_2c-
   new Particle("B*_2c-", -545, 0, "B-Meson", 100, -1, 7.35, 0.02, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-545));
   part->AddDecay(Particle::Decay(0, 0.3,  vector<int>{-511,-411}));
   part->AddDecay(Particle::Decay(0, 0.3,  vector<int>{-521,-421}));
   part->AddDecay(Particle::Decay(0, 0.2,  vector<int>{-513,-411}));
   part->AddDecay(Particle::Decay(0, 0.2,  vector<int>{-523,-421}));

   // Creating B*_c-
   new Particle("B*_c-", -543, 0, "B-Meson", 100, -1, 6.602, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-543));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-541,-22}));

   // Creating B_c-
   new Particle("B_c-", -541, 0, "B-Meson", 100, -1, 6.276, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-541));
   part->AddDecay(Particle::Decay(42, 0.24,  vector<int>{1,-2,-3,5}));
   part->AddDecay(Particle::Decay(42, 0.15,  vector<int>{-2,1,4,-4}));
   part->AddDecay(Particle::Decay(11, 0.122,  vector<int>{-4,3}));
   part->AddDecay(Particle::Decay(42, 0.065,  vector<int>{1,-3,-2,5}));
   part->AddDecay(Particle::Decay(42, 0.05,  vector<int>{-4,3,4,-4}));
   part->AddDecay(Particle::Decay(0, 0.047,  vector<int>{-16,15}));
   part->AddDecay(Particle::Decay(42, 0.042,  vector<int>{11,-12,-533}));
   part->AddDecay(Particle::Decay(42, 0.042,  vector<int>{13,-14,-533}));
   part->AddDecay(Particle::Decay(42, 0.037,  vector<int>{-2,4,1,-4}));
   part->AddDecay(Particle::Decay(42, 0.035,  vector<int>{-14,13,-443}));
   part->AddDecay(Particle::Decay(42, 0.035,  vector<int>{-12,11,-443}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-4,4,3,-4}));
   part->AddDecay(Particle::Decay(42, 0.014,  vector<int>{13,-14,-531}));
   part->AddDecay(Particle::Decay(42, 0.014,  vector<int>{11,-12,-531}));
   part->AddDecay(Particle::Decay(42, 0.014,  vector<int>{1,-2,-1,5}));
   part->AddDecay(Particle::Decay(42, 0.012,  vector<int>{-12,11,-441}));
   part->AddDecay(Particle::Decay(42, 0.012,  vector<int>{3,-2,-3,5}));
   part->AddDecay(Particle::Decay(42, 0.012,  vector<int>{-14,13,-441}));
   part->AddDecay(Particle::Decay(42, 0.008,  vector<int>{-2,3,4,-4}));
   part->AddDecay(Particle::Decay(42, 0.007,  vector<int>{-16,15,-443}));
   part->AddDecay(Particle::Decay(11, 0.006,  vector<int>{-4,1}));
   part->AddDecay(Particle::Decay(42, 0.003,  vector<int>{-16,15,-441}));
   part->AddDecay(Particle::Decay(42, 0.003,  vector<int>{3,-3,-2,5}));
   part->AddDecay(Particle::Decay(42, 0.003,  vector<int>{-4,1,4,-4}));
   part->AddDecay(Particle::Decay(42, 0.003,  vector<int>{1,-1,-2,5}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{13,-14,-513}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{-2,4,3,-4}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{11,-12,-513}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{11,-12,-511}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{-4,4,1,-4}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{13,-14,-511}));

   // Creating B*_2s0_bar
   new Particle("B*_2s0_bar", -535, 0, "B-Meson", 100, 0, 5.8397, 0.02, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-535));
   part->AddDecay(Particle::Decay(0, 0.3,  vector<int>{-521,321}));
   part->AddDecay(Particle::Decay(0, 0.3,  vector<int>{-511,311}));
   part->AddDecay(Particle::Decay(0, 0.2,  vector<int>{-523,321}));
   part->AddDecay(Particle::Decay(0, 0.2,  vector<int>{-513,311}));

   // Creating B*_s0_bar
   new Particle("B*_s0_bar", -533, 0, "B-Meson", 100, 0, 5.4128, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-533));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-531,-22}));

   // Creating B_s0_bar
   new Particle("B_s0_bar", -531, 0, "B-Meson", 100, 0, 5.3663, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-531));
   part->AddDecay(Particle::Decay(48, 0.4291,  vector<int>{-2,1,4,-3}));
   part->AddDecay(Particle::Decay(13, 0.08,  vector<int>{-2,4,1,-3}));
   part->AddDecay(Particle::Decay(13, 0.07,  vector<int>{-4,3,4,-3}));
   part->AddDecay(Particle::Decay(42, 0.055,  vector<int>{-14,13,433}));
   part->AddDecay(Particle::Decay(42, 0.055,  vector<int>{-12,11,433}));
   part->AddDecay(Particle::Decay(42, 0.03,  vector<int>{-16,15,433}));
   part->AddDecay(Particle::Decay(0, 0.025,  vector<int>{433,-433}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{-12,11,431}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{-14,13,431}));
   part->AddDecay(Particle::Decay(13, 0.02,  vector<int>{-4,4,3,-3}));
   part->AddDecay(Particle::Decay(0, 0.0185,  vector<int>{431,-433}));
   part->AddDecay(Particle::Decay(0, 0.018,  vector<int>{433,-20213}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{431,-431}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,-3}));
   part->AddDecay(Particle::Decay(0, 0.0135,  vector<int>{433,-431}));
   part->AddDecay(Particle::Decay(42, 0.012,  vector<int>{-14,13,435}));
   part->AddDecay(Particle::Decay(42, 0.012,  vector<int>{-12,11,435}));
   part->AddDecay(Particle::Decay(0, 0.011,  vector<int>{431,-213}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-16,15,431}));
   part->AddDecay(Particle::Decay(0, 0.009,  vector<int>{433,-213}));
   part->AddDecay(Particle::Decay(42, 0.008,  vector<int>{-14,13,20433}));
   part->AddDecay(Particle::Decay(42, 0.008,  vector<int>{-12,11,20433}));
   part->AddDecay(Particle::Decay(0, 0.0055,  vector<int>{431,-20213}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-12,11,10431}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-14,13,10433}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-14,13,10431}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-12,11,10433}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,-3}));
   part->AddDecay(Particle::Decay(0, 0.0042,  vector<int>{433,-211}));
   part->AddDecay(Particle::Decay(0, 0.0035,  vector<int>{431,-211}));
   part->AddDecay(Particle::Decay(0, 0.0025,  vector<int>{-20443,-333}));
   part->AddDecay(Particle::Decay(0, 0.0014,  vector<int>{-443,-333}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-20443,-221}));
   part->AddDecay(Particle::Decay(0, 0.0009,  vector<int>{-20443,-331}));
   part->AddDecay(Particle::Decay(0, 0.0007,  vector<int>{-441,-333}));
   part->AddDecay(Particle::Decay(0, 0.0004,  vector<int>{-443,-221}));
   part->AddDecay(Particle::Decay(0, 0.0004,  vector<int>{-443,-331}));
   part->AddDecay(Particle::Decay(0, 0.0002,  vector<int>{-441,-331}));
   part->AddDecay(Particle::Decay(0, 0.0002,  vector<int>{-441,-221}));

   // Creating B*_2-
   new Particle("B*_2-", -525, 0, "B-Meson", 100, -1, 5.7469, 0.02, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-525));
   part->AddDecay(Particle::Decay(0, 0.3,  vector<int>{-511,-211}));
   part->AddDecay(Particle::Decay(0, 0.16,  vector<int>{-513,-211}));
   part->AddDecay(Particle::Decay(0, 0.15,  vector<int>{-521,-111}));
   part->AddDecay(Particle::Decay(0, 0.13,  vector<int>{-513,-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{-523,-111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{-511,-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{-523,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.04,  vector<int>{-521,-211,211}));

   // Creating B*-
   new Particle("B*-", -523, 0, "Meson", 100, -1, 5.3251, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-523));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-521,-22}));

   // Creating B-
   new Particle("B-", -521, 0, "B-Meson", 100, -1, 5.27915, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-521));
   part->AddDecay(Particle::Decay(48, 0.4291,  vector<int>{-2,1,4,-2}));
   part->AddDecay(Particle::Decay(13, 0.08,  vector<int>{-2,4,1,-2}));
   part->AddDecay(Particle::Decay(13, 0.07,  vector<int>{-4,3,4,-2}));
   part->AddDecay(Particle::Decay(42, 0.055,  vector<int>{-14,13,423}));
   part->AddDecay(Particle::Decay(42, 0.055,  vector<int>{-12,11,423}));
   part->AddDecay(Particle::Decay(42, 0.03,  vector<int>{-16,15,423}));
   part->AddDecay(Particle::Decay(0, 0.025,  vector<int>{423,-433}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{-12,11,421}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{-14,13,421}));
   part->AddDecay(Particle::Decay(13, 0.02,  vector<int>{-4,4,3,-2}));
   part->AddDecay(Particle::Decay(0, 0.0185,  vector<int>{421,-433}));
   part->AddDecay(Particle::Decay(0, 0.018,  vector<int>{423,-20213}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{421,-431}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,-2}));
   part->AddDecay(Particle::Decay(0, 0.0135,  vector<int>{423,-431}));
   part->AddDecay(Particle::Decay(42, 0.012,  vector<int>{-14,13,425}));
   part->AddDecay(Particle::Decay(42, 0.012,  vector<int>{-12,11,425}));
   part->AddDecay(Particle::Decay(0, 0.011,  vector<int>{421,-213}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-16,15,421}));
   part->AddDecay(Particle::Decay(0, 0.009,  vector<int>{423,-213}));
   part->AddDecay(Particle::Decay(42, 0.008,  vector<int>{-14,13,20423}));
   part->AddDecay(Particle::Decay(42, 0.008,  vector<int>{-12,11,20423}));
   part->AddDecay(Particle::Decay(0, 0.0055,  vector<int>{421,-20213}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-12,11,10421}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-14,13,10423}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-14,13,10421}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-12,11,10423}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,-2}));
   part->AddDecay(Particle::Decay(0, 0.0042,  vector<int>{423,-211}));
   part->AddDecay(Particle::Decay(0, 0.0035,  vector<int>{421,-211}));
   part->AddDecay(Particle::Decay(0, 0.0025,  vector<int>{-20443,-323}));
   part->AddDecay(Particle::Decay(0, 0.0019,  vector<int>{-20443,-321}));
   part->AddDecay(Particle::Decay(0, 0.0014,  vector<int>{-443,-323}));
   part->AddDecay(Particle::Decay(0, 0.0008,  vector<int>{-443,-321}));
   part->AddDecay(Particle::Decay(0, 0.0007,  vector<int>{-441,-323}));
   part->AddDecay(Particle::Decay(0, 0.0004,  vector<int>{-441,-321}));

   // Creating B*_20_bar
   new Particle("B*_20_bar", -515, 0, "B-Meson", 100, 0, 5.7469, 0.02, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-515));
   part->AddDecay(Particle::Decay(0, 0.3,  vector<int>{-521,211}));
   part->AddDecay(Particle::Decay(0, 0.16,  vector<int>{-523,211}));
   part->AddDecay(Particle::Decay(0, 0.15,  vector<int>{-511,-111}));
   part->AddDecay(Particle::Decay(0, 0.13,  vector<int>{-523,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{-513,-111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{-521,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{-513,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.04,  vector<int>{-511,-211,211}));

   // Creating B*0_bar
   new Particle("B*0_bar", -513, 0, "B-Meson", 100, 0, 5.3251, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-513));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{-511,-22}));

   // Creating B0_bar
   new Particle("B0_bar", -511, 0, "B-Meson", 100, 0, 5.27953, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-511));
   part->AddDecay(Particle::Decay(48, 0.4291,  vector<int>{-2,1,4,-1}));
   part->AddDecay(Particle::Decay(13, 0.08,  vector<int>{-2,4,1,-1}));
   part->AddDecay(Particle::Decay(13, 0.07,  vector<int>{-4,3,4,-1}));
   part->AddDecay(Particle::Decay(42, 0.055,  vector<int>{-14,13,413}));
   part->AddDecay(Particle::Decay(42, 0.055,  vector<int>{-12,11,413}));
   part->AddDecay(Particle::Decay(42, 0.03,  vector<int>{-16,15,413}));
   part->AddDecay(Particle::Decay(0, 0.025,  vector<int>{413,-433}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{-12,11,411}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{-14,13,411}));
   part->AddDecay(Particle::Decay(13, 0.02,  vector<int>{-4,4,3,-1}));
   part->AddDecay(Particle::Decay(0, 0.0185,  vector<int>{411,-433}));
   part->AddDecay(Particle::Decay(0, 0.018,  vector<int>{413,-20213}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{411,-431}));
   part->AddDecay(Particle::Decay(42, 0.015,  vector<int>{-2,1,2,-1}));
   part->AddDecay(Particle::Decay(0, 0.0135,  vector<int>{413,-431}));
   part->AddDecay(Particle::Decay(42, 0.012,  vector<int>{-14,13,415}));
   part->AddDecay(Particle::Decay(42, 0.012,  vector<int>{-12,11,415}));
   part->AddDecay(Particle::Decay(0, 0.011,  vector<int>{411,-213}));
   part->AddDecay(Particle::Decay(42, 0.01,  vector<int>{-16,15,411}));
   part->AddDecay(Particle::Decay(0, 0.009,  vector<int>{413,-213}));
   part->AddDecay(Particle::Decay(42, 0.008,  vector<int>{-14,13,20413}));
   part->AddDecay(Particle::Decay(42, 0.008,  vector<int>{-12,11,20413}));
   part->AddDecay(Particle::Decay(0, 0.0055,  vector<int>{411,-20213}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-12,11,10411}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-14,13,10413}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-14,13,10411}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-12,11,10413}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{-4,3,2,-1}));
   part->AddDecay(Particle::Decay(0, 0.0042,  vector<int>{413,-211}));
   part->AddDecay(Particle::Decay(0, 0.0035,  vector<int>{411,-211}));
   part->AddDecay(Particle::Decay(0, 0.0025,  vector<int>{-20443,-313}));
   part->AddDecay(Particle::Decay(0, 0.0019,  vector<int>{-20443,-311}));
   part->AddDecay(Particle::Decay(0, 0.0014,  vector<int>{-443,-313}));
   part->AddDecay(Particle::Decay(0, 0.0008,  vector<int>{-443,-311}));
   part->AddDecay(Particle::Decay(0, 0.0007,  vector<int>{-441,-313}));
   part->AddDecay(Particle::Decay(0, 0.0004,  vector<int>{-441,-311}));
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
