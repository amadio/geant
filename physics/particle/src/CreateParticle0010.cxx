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
void CreateParticle0010() {

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
   Particle *part = 0;
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

   // Creating D*_2s-
   new Particle("D*_2s-", -435, 0, "CharmedMeson", 100, -1, 2.5726, 0.015, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-435));
   part->AddDecay(Particle::Decay(0, 0.4,  vector<int>{-421,-321}));
   part->AddDecay(Particle::Decay(0, 0.4,  vector<int>{-411,-311}));
   part->AddDecay(Particle::Decay(0, 0.1,  vector<int>{-423,-321}));
   part->AddDecay(Particle::Decay(0, 0.1,  vector<int>{-413,-311}));

   // Creating D*_s-
   new Particle("D*_s-", -433, 0, "CharmedMeson", 100, -1, 2.1123, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-433));
   part->AddDecay(Particle::Decay(0, 0.94,  vector<int>{-431,-22}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{-431,-111}));

   // Creating D_s-
   new Particle("D_s-", -431, 0, "CharmedMeson", 100, -1, 1.9685, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-431));
   part->AddDecay(Particle::Decay(13, 0.25,  vector<int>{-2,1,-3,3}));
   part->AddDecay(Particle::Decay(13, 0.0952,  vector<int>{-2,1}));
   part->AddDecay(Particle::Decay(0, 0.095,  vector<int>{-331,-213}));
   part->AddDecay(Particle::Decay(0, 0.079,  vector<int>{-221,-213}));
   part->AddDecay(Particle::Decay(0, 0.052,  vector<int>{-333,-213}));
   part->AddDecay(Particle::Decay(0, 0.05,  vector<int>{-323,313}));
   part->AddDecay(Particle::Decay(0, 0.037,  vector<int>{-331,-211}));
   part->AddDecay(Particle::Decay(0, 0.033,  vector<int>{-323,311}));
   part->AddDecay(Particle::Decay(42, 0.03,  vector<int>{11,-12,-333}));
   part->AddDecay(Particle::Decay(42, 0.03,  vector<int>{13,-14,-333}));
   part->AddDecay(Particle::Decay(0, 0.028,  vector<int>{-333,-211}));
   part->AddDecay(Particle::Decay(0, 0.028,  vector<int>{-321,311}));
   part->AddDecay(Particle::Decay(0, 0.026,  vector<int>{-321,313}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{11,-12,-331}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{11,-12,-221}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{13,-14,-221}));
   part->AddDecay(Particle::Decay(42, 0.02,  vector<int>{13,-14,-331}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{-221,-211}));
   part->AddDecay(Particle::Decay(0, 0.01,  vector<int>{15,-16}));
   part->AddDecay(Particle::Decay(0, 0.01,  vector<int>{-2212,2112}));
   part->AddDecay(Particle::Decay(0, 0.0078,  vector<int>{-10221,-211}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{13,-14,-321,321}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{13,-14,-311,311}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-221,-321}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-331,-321}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-333,-321}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-221,-323}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{11,-12,-311,311}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{11,-12,-321,321}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-213,-113}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-213,-111}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-211,-113}));

   // Creating D*_20_bar
   new Particle("D*_20_bar", -425, 0, "CharmedMeson", 100, 0, 2.4611, 0.023, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-425));
   part->AddDecay(Particle::Decay(0, 0.3,  vector<int>{-411,211}));
   part->AddDecay(Particle::Decay(0, 0.16,  vector<int>{-413,211}));
   part->AddDecay(Particle::Decay(0, 0.15,  vector<int>{-421,-111}));
   part->AddDecay(Particle::Decay(0, 0.13,  vector<int>{-413,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{-423,-111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{-411,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{-423,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.04,  vector<int>{-421,-211,211}));

   // Creating D*0_bar
   new Particle("D*0_bar", -423, 0, "CharmedMeson", 100, 0, 2.00697, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-423));
   part->AddDecay(Particle::Decay(3, 0.619,  vector<int>{-421,-111}));
   part->AddDecay(Particle::Decay(0, 0.381,  vector<int>{-421,-22}));

   // Creating D0_bar
   new Particle("D0_bar", -421, 0, "CharmedMeson", 100, 0, 1.86484, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-421));
   part->AddDecay(Particle::Decay(0, 0.0923,  vector<int>{321,-211,-111,-111}));
   part->AddDecay(Particle::Decay(0, 0.074,  vector<int>{321,-20213}));
   part->AddDecay(Particle::Decay(0, 0.073,  vector<int>{321,-213}));
   part->AddDecay(Particle::Decay(0, 0.067,  vector<int>{311,-211,211,-111,-111}));
   part->AddDecay(Particle::Decay(0, 0.062,  vector<int>{323,-213}));
   part->AddDecay(Particle::Decay(0, 0.0511,  vector<int>{311,-113,-111,-111,-111}));
   part->AddDecay(Particle::Decay(0, 0.045,  vector<int>{323,-211}));
   part->AddDecay(Particle::Decay(0, 0.0365,  vector<int>{321,-211}));
   part->AddDecay(Particle::Decay(42, 0.034,  vector<int>{11,-12,321}));
   part->AddDecay(Particle::Decay(42, 0.034,  vector<int>{13,-14,321}));
   part->AddDecay(Particle::Decay(42, 0.027,  vector<int>{11,-12,323}));
   part->AddDecay(Particle::Decay(42, 0.027,  vector<int>{13,-14,323}));
   part->AddDecay(Particle::Decay(0, 0.025,  vector<int>{311,-223}));
   part->AddDecay(Particle::Decay(0, 0.024,  vector<int>{321,-211,-211,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.022,  vector<int>{311,-211,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.021,  vector<int>{313,-111}));
   part->AddDecay(Particle::Decay(0, 0.021,  vector<int>{313,-221}));
   part->AddDecay(Particle::Decay(0, 0.021,  vector<int>{311,-111}));
   part->AddDecay(Particle::Decay(0, 0.018,  vector<int>{311,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.018,  vector<int>{321,-211,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.017,  vector<int>{-211,-211,211,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.016,  vector<int>{313,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{-211,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.015,  vector<int>{313,-113}));
   part->AddDecay(Particle::Decay(0, 0.011,  vector<int>{321,-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.0109,  vector<int>{10323,-211}));
   part->AddDecay(Particle::Decay(0, 0.009,  vector<int>{311,-321,321,-111}));
   part->AddDecay(Particle::Decay(0, 0.0088,  vector<int>{311,-333}));
   part->AddDecay(Particle::Decay(0, 0.0085,  vector<int>{311,-211,-211,211,211}));
   part->AddDecay(Particle::Decay(0, 0.0077,  vector<int>{313,-211,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.0075,  vector<int>{-211,-211,211,211}));
   part->AddDecay(Particle::Decay(0, 0.0063,  vector<int>{321,-211,-113}));
   part->AddDecay(Particle::Decay(0, 0.0061,  vector<int>{311,-113}));
   part->AddDecay(Particle::Decay(0, 0.0052,  vector<int>{321,-321,311}));
   part->AddDecay(Particle::Decay(0, 0.0041,  vector<int>{321,-321}));
   part->AddDecay(Particle::Decay(42, 0.004,  vector<int>{13,-14,323,-111}));
   part->AddDecay(Particle::Decay(42, 0.004,  vector<int>{11,-12,313,211}));
   part->AddDecay(Particle::Decay(42, 0.004,  vector<int>{11,-12,323,-111}));
   part->AddDecay(Particle::Decay(42, 0.004,  vector<int>{13,-14,313,211}));
   part->AddDecay(Particle::Decay(0, 0.0036,  vector<int>{313,-321,211}));
   part->AddDecay(Particle::Decay(0, 0.0035,  vector<int>{321,-323}));
   part->AddDecay(Particle::Decay(0, 0.0034,  vector<int>{321,-311,-211}));
   part->AddDecay(Particle::Decay(0, 0.0028,  vector<int>{-321,321,-211,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.0027,  vector<int>{313,-313}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{13,-14,311,211}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{13,-14,321,-111}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{11,-12,211}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{11,-12,213}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{323,-321}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{13,-14,211}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{13,-14,213}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{11,-12,311,211}));
   part->AddDecay(Particle::Decay(42, 0.002,  vector<int>{11,-12,321,-111}));
   part->AddDecay(Particle::Decay(0, 0.0018,  vector<int>{-333,-113}));
   part->AddDecay(Particle::Decay(0, 0.0016,  vector<int>{-111,-111}));
   part->AddDecay(Particle::Decay(0, 0.0016,  vector<int>{-211,211}));
   part->AddDecay(Particle::Decay(0, 0.0011,  vector<int>{311,-311}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{313,-311}));
   part->AddDecay(Particle::Decay(0, 0.0009,  vector<int>{-310,-310,-310}));
   part->AddDecay(Particle::Decay(0, 0.0006,  vector<int>{-333,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.0004,  vector<int>{-113,-211,-211,211,211}));

   // Creating D*_2-
   new Particle("D*_2-", -415, 0, "CharmedMeson", 100, -1, 2.4601, 0.023, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-415));
   part->AddDecay(Particle::Decay(0, 0.3,  vector<int>{-421,-211}));
   part->AddDecay(Particle::Decay(0, 0.16,  vector<int>{-423,-211}));
   part->AddDecay(Particle::Decay(0, 0.15,  vector<int>{-411,-111}));
   part->AddDecay(Particle::Decay(0, 0.13,  vector<int>{-423,-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{-413,-111}));
   part->AddDecay(Particle::Decay(0, 0.08,  vector<int>{-421,-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.06,  vector<int>{-413,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.04,  vector<int>{-411,-211,211}));

   // Creating D*-
   new Particle("D*-", -413, 0, "CharmedMeson", 100, -1, 2.01027, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-413));
   part->AddDecay(Particle::Decay(3, 0.683,  vector<int>{-421,-211}));
   part->AddDecay(Particle::Decay(3, 0.306,  vector<int>{-411,-111}));
   part->AddDecay(Particle::Decay(0, 0.011,  vector<int>{-411,-22}));

   // Creating D-
   new Particle("D-", -411, 0, "CharmedMeson", 100, -1, 1.86962, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-411));
   part->AddDecay(Particle::Decay(0, 0.087,  vector<int>{311,-211,-211,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.076,  vector<int>{311,-20213}));
   part->AddDecay(Particle::Decay(42, 0.07,  vector<int>{11,-12,311}));
   part->AddDecay(Particle::Decay(42, 0.07,  vector<int>{13,-14,311}));
   part->AddDecay(Particle::Decay(0, 0.067,  vector<int>{321,-211,-211}));
   part->AddDecay(Particle::Decay(0, 0.066,  vector<int>{311,-213}));
   part->AddDecay(Particle::Decay(42, 0.065,  vector<int>{11,-12,313}));
   part->AddDecay(Particle::Decay(42, 0.065,  vector<int>{13,-14,313}));
   part->AddDecay(Particle::Decay(0, 0.045,  vector<int>{20313,-211}));
   part->AddDecay(Particle::Decay(0, 0.041,  vector<int>{313,-213}));
   part->AddDecay(Particle::Decay(0, 0.027,  vector<int>{311,-321,311}));
   part->AddDecay(Particle::Decay(0, 0.026,  vector<int>{313,-323}));
   part->AddDecay(Particle::Decay(0, 0.026,  vector<int>{311,-211}));
   part->AddDecay(Particle::Decay(0, 0.022,  vector<int>{321,-211,-211,-111,-111}));
   part->AddDecay(Particle::Decay(0, 0.0218,  vector<int>{-211,-211,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.019,  vector<int>{-333,-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.019,  vector<int>{313,-211}));
   part->AddDecay(Particle::Decay(0, 0.012,  vector<int>{311,-211,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.012,  vector<int>{311,-211,-111}));
   part->AddDecay(Particle::Decay(42, 0.011,  vector<int>{13,-14,323,-211}));
   part->AddDecay(Particle::Decay(42, 0.011,  vector<int>{11,-12,313,-111}));
   part->AddDecay(Particle::Decay(42, 0.011,  vector<int>{11,-12,323,-211}));
   part->AddDecay(Particle::Decay(42, 0.011,  vector<int>{13,-14,313,-111}));
   part->AddDecay(Particle::Decay(0, 0.009,  vector<int>{321,-211,-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.008,  vector<int>{321,-213,-211}));
   part->AddDecay(Particle::Decay(0, 0.0073,  vector<int>{311,-321}));
   part->AddDecay(Particle::Decay(0, 0.0066,  vector<int>{-221,-211}));
   part->AddDecay(Particle::Decay(0, 0.006,  vector<int>{-333,-211}));
   part->AddDecay(Particle::Decay(0, 0.0057,  vector<int>{313,-211,-113}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-333,-213}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{11,-12,311,-111}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{11,-12,321,-211}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{13,-14,311,-111}));
   part->AddDecay(Particle::Decay(0, 0.005,  vector<int>{-221,-213}));
   part->AddDecay(Particle::Decay(42, 0.005,  vector<int>{13,-14,321,-211}));
   part->AddDecay(Particle::Decay(0, 0.0047,  vector<int>{313,-321}));
   part->AddDecay(Particle::Decay(0, 0.0047,  vector<int>{311,-323}));
   part->AddDecay(Particle::Decay(0, 0.004,  vector<int>{321,-321,-211}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{-331,-211}));
   part->AddDecay(Particle::Decay(0, 0.003,  vector<int>{-331,-213}));
   part->AddDecay(Particle::Decay(0, 0.0028,  vector<int>{-113,-211,-211,211,-111}));
   part->AddDecay(Particle::Decay(0, 0.0022,  vector<int>{-211,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.002,  vector<int>{313,-211,-211,211}));
   part->AddDecay(Particle::Decay(0, 0.0019,  vector<int>{321,-113,-211,-211,-111}));
   part->AddDecay(Particle::Decay(0, 0.0015,  vector<int>{-211,-211,-211,211,211}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-111,-211}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{11,-12,-221}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{11,-12,-331}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{11,-12,-113}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{11,-12,-223}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-223,-211}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{-223,-213}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{11,-12,-111}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{321,-211,-211,-211,211}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{13,-14,-111}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{13,-14,-221}));
   part->AddDecay(Particle::Decay(0, 0.001,  vector<int>{311,-113,-211,-211,211}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{13,-14,-331}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{13,-14,-113}));
   part->AddDecay(Particle::Decay(42, 0.001,  vector<int>{13,-14,-223}));
   part->AddDecay(Particle::Decay(0, 0.0006,  vector<int>{-113,-213}));
   part->AddDecay(Particle::Decay(0, 0.0006,  vector<int>{-111,-213}));
   part->AddDecay(Particle::Decay(0, 0.0006,  vector<int>{-113,-211}));

   // Creating phi3(1850)_bar
   new Particle("phi3(1850)_bar", -337, 0, "Unknown", 100, 0, 1.854, 0.087, 100, 100, 0, 100, 1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
