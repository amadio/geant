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
void CreateParticle0013() {

   // Creating tau-
   new Particle("tau-", 15, 1, "Lepton", 100, -1, 1.77684, 0, -100, -1, -100, -1, -1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(15));
   part->AddDecay(Particle::Decay(0, 0.2494,  vector<int>{16,-213}));
   part->AddDecay(Particle::Decay(42, 0.1783,  vector<int>{-12,11,16}));
   part->AddDecay(Particle::Decay(42, 0.1735,  vector<int>{-14,13,16}));
   part->AddDecay(Particle::Decay(0, 0.1131,  vector<int>{16,-211}));
   part->AddDecay(Particle::Decay(41, 0.09,  vector<int>{16,-213,111}));
   part->AddDecay(Particle::Decay(41, 0.08,  vector<int>{16,-211,113}));
   part->AddDecay(Particle::Decay(41, 0.0191,  vector<int>{16,-211,223}));
   part->AddDecay(Particle::Decay(41, 0.0133,  vector<int>{16,-211,113,111}));
   part->AddDecay(Particle::Decay(0, 0.012,  vector<int>{16,-323}));
   part->AddDecay(Particle::Decay(41, 0.011,  vector<int>{16,-211,211,-211}));
   part->AddDecay(Particle::Decay(41, 0.01,  vector<int>{16,-213,111,111}));
   part->AddDecay(Particle::Decay(0, 0.0071,  vector<int>{16,-321}));
   part->AddDecay(Particle::Decay(41, 0.0067,  vector<int>{16,-213,211,-211}));
   part->AddDecay(Particle::Decay(41, 0.005,  vector<int>{16,-213,113}));
   part->AddDecay(Particle::Decay(41, 0.0035,  vector<int>{16,-213,223}));
   part->AddDecay(Particle::Decay(41, 0.0034,  vector<int>{16,-321,321,-211}));
   part->AddDecay(Particle::Decay(41, 0.003,  vector<int>{16,-211,111}));
   part->AddDecay(Particle::Decay(41, 0.0027,  vector<int>{16,-211,111,111}));
   part->AddDecay(Particle::Decay(41, 0.00205,  vector<int>{16,-211,310,111}));
   part->AddDecay(Particle::Decay(41, 0.00205,  vector<int>{16,-211,130,111}));
   part->AddDecay(Particle::Decay(41, 0.0015,  vector<int>{16,-213,221}));
   part->AddDecay(Particle::Decay(41, 0.0014,  vector<int>{16,-211,111,111,111}));
   part->AddDecay(Particle::Decay(41, 0.0012,  vector<int>{16,-213,111,111,111}));
   part->AddDecay(Particle::Decay(41, 0.0011,  vector<int>{16,-213,113,111,111}));
   part->AddDecay(Particle::Decay(41, 0.00078,  vector<int>{16,-321,310}));
   part->AddDecay(Particle::Decay(41, 0.00078,  vector<int>{16,-321,130}));
   part->AddDecay(Particle::Decay(41, 0.00075,  vector<int>{16,-211,113,113}));
   part->AddDecay(Particle::Decay(41, 0.00075,  vector<int>{16,-323,111}));
   part->AddDecay(Particle::Decay(41, 0.00069,  vector<int>{16,-321,310,111}));
   part->AddDecay(Particle::Decay(41, 0.00069,  vector<int>{16,-321,130,111}));
   part->AddDecay(Particle::Decay(41, 0.0006,  vector<int>{16,-211,223,111}));
   part->AddDecay(Particle::Decay(41, 0.00051,  vector<int>{16,-211,310,130}));
   part->AddDecay(Particle::Decay(41, 0.0005,  vector<int>{16,-211,211,-211,111}));
   part->AddDecay(Particle::Decay(41, 0.0004,  vector<int>{16,-323,111,111}));
   part->AddDecay(Particle::Decay(41, 0.0004,  vector<int>{16,-321,111}));
   part->AddDecay(Particle::Decay(41, 0.00025,  vector<int>{16,-211,130}));
   part->AddDecay(Particle::Decay(41, 0.00025,  vector<int>{16,-211,310,310}));
   part->AddDecay(Particle::Decay(41, 0.00025,  vector<int>{16,-211,310}));
   part->AddDecay(Particle::Decay(41, 0.00025,  vector<int>{16,-211,130,130}));
   part->AddDecay(Particle::Decay(41, 0.00022,  vector<int>{16,-211,113,113,111}));
   part->AddDecay(Particle::Decay(41, 0.00021,  vector<int>{16,-211,221,111}));
   part->AddDecay(Particle::Decay(41, 0.0002,  vector<int>{16,-211,-213,211,111}));
   part->AddDecay(Particle::Decay(41, 0.0002,  vector<int>{16,-211,113,111,111}));
   part->AddDecay(Particle::Decay(41, 0.0002,  vector<int>{16,-213,113,111}));
   part->AddDecay(Particle::Decay(41, 0.0002,  vector<int>{16,-211,213,-213}));
   part->AddDecay(Particle::Decay(41, 0.0002,  vector<int>{16,-211,213,-211,111}));
   part->AddDecay(Particle::Decay(41, 0.0001,  vector<int>{16,-321,111,111,111}));
   part->AddDecay(Particle::Decay(41, 0.0001,  vector<int>{16,-211,221,221}));
   part->AddDecay(Particle::Decay(41, 6e-05,  vector<int>{16,-323,111,111}));
   part->AddDecay(Particle::Decay(41, 6e-05,  vector<int>{16,-211,221}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,15}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,15}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{-24,16}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,15}));

   // Creating nu_tau
   new Particle("nu_tau", 16, 1, "Lepton", 100, 0, 0, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(16));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,16}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{24,15}));

   // Creating tau'-
   new Particle("tau'-", 17, 1, "Lepton", 100, -1, 400, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(17));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{22,17}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,17}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-24,18}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{25,17}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-37,18}));

   // Creating nu'_tau
   new Particle("nu'_tau", 18, 1, "Lepton", 100, 0, 0, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(18));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{23,18}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,17}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{37,17}));

   // Creating g
   new Particle("g", 21, 0, "GaugeBoson", 100, 0, 0, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(21));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{1,-1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{2,-2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{3,-3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{6,-6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{7,-7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{8,-8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{21,21}));

   // Creating gamma
   new Particle("gamma", 22, 0, "GaugeBoson", 100, 0, 0, 0, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(22));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{1,-1}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{2,-2}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{3,-3}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{6,-6}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{7,-7}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{8,-8}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{15,-15}));
   part->AddDecay(Particle::Decay(102, 0,  vector<int>{17,-17}));

   // Creating Z0
   new Particle("Z0", 23, 0, "GaugeBoson", 100, 0, 91.187, 2.48009, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(23));
   part->AddDecay(Particle::Decay(32, 0.154075,  vector<int>{1,-1}));
   part->AddDecay(Particle::Decay(32, 0.154072,  vector<int>{3,-3}));
   part->AddDecay(Particle::Decay(32, 0.152196,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(32, 0.119483,  vector<int>{2,-2}));
   part->AddDecay(Particle::Decay(32, 0.119346,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(0, 0.0667521,  vector<int>{12,-12}));
   part->AddDecay(Particle::Decay(0, 0.0667521,  vector<int>{14,-14}));
   part->AddDecay(Particle::Decay(0, 0.0667521,  vector<int>{16,-16}));
   part->AddDecay(Particle::Decay(0, 0.033549,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(0, 0.033549,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(0, 0.033473,  vector<int>{15,-15}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{6,-6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-7}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{8,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{17,-17}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{18,-18}));

   // Creating W+
   new Particle("W+", 24, 1, "GaugeBoson", 100, 1, 80.398, 2.07002, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(24));
   part->AddDecay(Particle::Decay(32, 0.321502,  vector<int>{-1,2}));
   part->AddDecay(Particle::Decay(32, 0.320778,  vector<int>{-3,4}));
   part->AddDecay(Particle::Decay(0, 0.108062,  vector<int>{-11,12}));
   part->AddDecay(Particle::Decay(0, 0.108062,  vector<int>{-13,14}));
   part->AddDecay(Particle::Decay(0, 0.107983,  vector<int>{-15,16}));
   part->AddDecay(Particle::Decay(32, 0.016509,  vector<int>{-3,2}));
   part->AddDecay(Particle::Decay(32, 0.016502,  vector<int>{-1,4}));
   part->AddDecay(Particle::Decay(32, 0.000591001,  vector<int>{-5,4}));
   part->AddDecay(Particle::Decay(32, 1e-05,  vector<int>{-5,2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-1,8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-5,6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-5,8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,4}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-3,6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-3,8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-1,6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-17,18}));

   // Creating h0
   new Particle("h0", 25, 0, "GaugeBoson", 100, 0, 80, 0.00237, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(25));
   part->AddDecay(Particle::Decay(32, 0.852249,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(0, 0.06883,  vector<int>{15,-15}));
   part->AddDecay(Particle::Decay(32, 0.053489,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(0, 0.023981,  vector<int>{21,21}));
   part->AddDecay(Particle::Decay(0, 0.000879,  vector<int>{22,22}));
   part->AddDecay(Particle::Decay(32, 0.000327,  vector<int>{3,-3}));
   part->AddDecay(Particle::Decay(0, 0.000244,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(32, 1e-06,  vector<int>{1,-1}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{2,-2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{6,-6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{17,-17}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-7}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{8,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{22,23}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{23,23}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1000023}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000023}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000025}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000023}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000025}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000035}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-1000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,-2000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-2000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,2000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-1000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-2000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-2000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,2000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-1000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,-2000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-2000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,2000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-1000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-2000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-2000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,2000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-1000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,-2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-2000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-2000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,2000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-1000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000011,-2000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-2000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,2000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-1000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-2000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-2000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,2000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-1000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000013,-2000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-2000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,2000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-1000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-2000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-2000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,2000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-1000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000015,-2000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-2000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,2000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-1000016}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-2000016}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-2000016}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,2000016}));

   // Creating reggeon
   new Particle("reggeon", 28, 0, "GaugeBoson", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating pomeron
   new Particle("pomeron", 29, 0, "GaugeBoson", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating Z'0
   new Particle("Z'0", 32, 0, "GaugeBoson", 100, 0, 500, 14.5485, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(32));
   part->AddDecay(Particle::Decay(32, 0.145869,  vector<int>{1,-1}));
   part->AddDecay(Particle::Decay(32, 0.145869,  vector<int>{3,-3}));
   part->AddDecay(Particle::Decay(32, 0.14581,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(32, 0.113303,  vector<int>{2,-2}));
   part->AddDecay(Particle::Decay(32, 0.113298,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(0, 0.0636061,  vector<int>{12,-12}));
   part->AddDecay(Particle::Decay(0, 0.0636061,  vector<int>{14,-14}));
   part->AddDecay(Particle::Decay(0, 0.0636061,  vector<int>{16,-16}));
   part->AddDecay(Particle::Decay(32, 0.0490131,  vector<int>{6,-6}));
   part->AddDecay(Particle::Decay(0, 0.0320071,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(0, 0.0320071,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(0, 0.0320041,  vector<int>{15,-15}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-7}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{8,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{17,-17}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{18,-18}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,-24}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{37,-37}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{23,22}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{23,25}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{25,36}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{35,36}));

   // Creating Z"0
   new Particle("Z\"0", 33, 0, "GaugeBoson", 100, 0, 900, 0, -100, -1, -100, -1, -1);

   // Creating W'+
   new Particle("W'+", 34, 1, "GaugeBoson", 100, 1, 500, 16.6708, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(34));
   part->AddDecay(Particle::Decay(32, 0.251276,  vector<int>{-1,2}));
   part->AddDecay(Particle::Decay(32, 0.250816,  vector<int>{-3,4}));
   part->AddDecay(Particle::Decay(32, 0.215459,  vector<int>{-5,6}));
   part->AddDecay(Particle::Decay(0, 0.085262,  vector<int>{-11,12}));
   part->AddDecay(Particle::Decay(0, 0.085262,  vector<int>{-13,14}));
   part->AddDecay(Particle::Decay(0, 0.08526,  vector<int>{-15,16}));
   part->AddDecay(Particle::Decay(32, 0.012903,  vector<int>{-3,2}));
   part->AddDecay(Particle::Decay(32, 0.012903,  vector<int>{-1,4}));
   part->AddDecay(Particle::Decay(32, 0.000465,  vector<int>{-5,4}));
   part->AddDecay(Particle::Decay(32, 0.00038,  vector<int>{-3,6}));
   part->AddDecay(Particle::Decay(32, 8e-06,  vector<int>{-5,2}));
   part->AddDecay(Particle::Decay(32, 6e-06,  vector<int>{-1,6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,4}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-3,8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-1,8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-5,8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-17,18}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,23}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,22}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,25}));

   // Creating H0
   new Particle("H0", 35, 0, "GaugeBoson", 100, 0, 300, 8.42842, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(35));
   part->AddDecay(Particle::Decay(0, 0.688641,  vector<int>{24,-24}));
   part->AddDecay(Particle::Decay(0, 0.306171,  vector<int>{23,23}));
   part->AddDecay(Particle::Decay(0, 0.003799,  vector<int>{25,25}));
   part->AddDecay(Particle::Decay(32, 0.000754001,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(0, 0.000439,  vector<int>{21,21}));
   part->AddDecay(Particle::Decay(0, 7.40001e-05,  vector<int>{15,-15}));
   part->AddDecay(Particle::Decay(0, 6.10001e-05,  vector<int>{22,23}));
   part->AddDecay(Particle::Decay(32, 4.6e-05,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(0, 1.5e-05,  vector<int>{22,22}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{3,-3}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{17,-17}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{1,-1}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{2,-2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{6,-6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-7}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{8,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{23,25}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{36,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1000023}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000023}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000025}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000023}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000025}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000035}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-1000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,-2000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-2000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,2000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-1000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-2000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-2000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,2000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-1000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,-2000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-2000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,2000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-1000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-2000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-2000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,2000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-1000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,-2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-2000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-2000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,2000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-1000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000011,-2000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-2000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,2000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-1000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-2000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-2000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,2000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-1000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000013,-2000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-2000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,2000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-1000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-2000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-2000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,2000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-1000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000015,-2000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-2000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,2000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-1000016}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-2000016}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-2000016}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,2000016}));

   // Creating A0
   new Particle("A0", 36, 0, "GaugeBoson", 100, 0, 300, 4.92026, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(36));
   part->AddDecay(Particle::Decay(0, 0.996235,  vector<int>{23,25}));
   part->AddDecay(Particle::Decay(0, 0.002256,  vector<int>{21,21}));
   part->AddDecay(Particle::Decay(32, 0.001292,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(0, 0.000126,  vector<int>{15,-15}));
   part->AddDecay(Particle::Decay(32, 7.90002e-05,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(0, 1e-05,  vector<int>{22,22}));
   part->AddDecay(Particle::Decay(0, 2e-06,  vector<int>{22,23}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{8,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{3,-3}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{17,-17}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{1,-1}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{2,-2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{6,-6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{23,23}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{24,-24}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-7}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1000023}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000023}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000025}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000022}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000023}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000025}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000035}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-1000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000001,-2000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000001,-2000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,2000001}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-1000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000002,-2000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000002,-2000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,2000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-1000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000003,-2000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000003,-2000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,2000003}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-1000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000004,-2000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000004,-2000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,2000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-1000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000005,-2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000005,-2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-1000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-2000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-2000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,2000006}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-1000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000011,-2000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000011,-2000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,2000011}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-1000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000012,-2000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000012,-2000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,2000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-1000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000013,-2000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,-2000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,2000013}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-1000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-2000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-2000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000014,2000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-1000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000015,-2000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-2000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,2000015}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-1000016}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-2000016}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-2000016}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,2000016}));

   // Creating H+
   new Particle("H+", 37, 1, "GaugeBoson", 100, 1, 300, 5.75967, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(37));
   part->AddDecay(Particle::Decay(0, 0.929792,  vector<int>{24,25}));
   part->AddDecay(Particle::Decay(32, 0.067484,  vector<int>{-5,6}));
   part->AddDecay(Particle::Decay(0, 0.002701,  vector<int>{-15,16}));
   part->AddDecay(Particle::Decay(32, 1.3e-05,  vector<int>{-3,4}));
   part->AddDecay(Particle::Decay(0, 1e-05,  vector<int>{-13,14}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-1,2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-17,18}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-11,12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000024}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,1000037}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-1000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-1000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000006,-2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000006,-2000005}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,1000002}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,1000004}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,1000012}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000013,1000014}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000015,1000016}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000015,1000016}));

   // Creating eta_tech0
   new Particle("eta_tech0", 38, 0, "Unknown", 100, 0, 350, 0.10158, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(38));
   part->AddDecay(Particle::Decay(32, 0.547101,  vector<int>{21,21}));
   part->AddDecay(Particle::Decay(32, 0.452899,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{6,-6}));

   // Creating LQ_ue
   new Particle("LQ_ue", 39, 1, "Unknown", 100, -0.333333, 200, 0.39162, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(39));
   part->AddDecay(Particle::Decay(0, 1,  vector<int>{2,11}));

   // Creating R0
   new Particle("R0", 40, 1, "Unknown", 100, 0, 5000, 417.465, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(40));
   part->AddDecay(Particle::Decay(32, 0.215134,  vector<int>{1,-3}));
   part->AddDecay(Particle::Decay(32, 0.215134,  vector<int>{2,-4}));
   part->AddDecay(Particle::Decay(32, 0.215133,  vector<int>{3,-5}));
   part->AddDecay(Particle::Decay(32, 0.214738,  vector<int>{4,-6}));
   part->AddDecay(Particle::Decay(0, 0.0699301,  vector<int>{11,-13}));
   part->AddDecay(Particle::Decay(0, 0.0699301,  vector<int>{13,-15}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{5,-7}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{6,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{15,-17}));

   // Creating pi_tech0
   new Particle("pi_tech0", 51, 0, "Unknown", 100, 0, 110, 0.04104, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(51));
   part->AddDecay(Particle::Decay(32, 0.596654,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(32, 0.316112,  vector<int>{21,21}));
   part->AddDecay(Particle::Decay(0, 0.050055,  vector<int>{15,-15}));
   part->AddDecay(Particle::Decay(32, 0.036777,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(32, 0.000225,  vector<int>{3,-3}));
   part->AddDecay(Particle::Decay(0, 0.000177,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{6,-6}));

   // Creating pi_tech+
   new Particle("pi_tech+", 52, 1, "Unknown", 100, 1, 110, 0.0105, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(52));
   part->AddDecay(Particle::Decay(32, 0.90916,  vector<int>{4,-5}));
   part->AddDecay(Particle::Decay(0, 0.048905,  vector<int>{-15,16}));
   part->AddDecay(Particle::Decay(32, 0.041762,  vector<int>{4,-3}));
   part->AddDecay(Particle::Decay(0, 0.000173,  vector<int>{-13,14}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{24,5,-5}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-11,12}));

   // Creating pi'_tech0
   new Particle("pi'_tech0", 53, 0, "Unknown", 100, 0, 110, 0.02807, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(53));
   part->AddDecay(Particle::Decay(32, 0.872445,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(0, 0.0731921,  vector<int>{15,-15}));
   part->AddDecay(Particle::Decay(32, 0.0537761,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(32, 0.000328,  vector<int>{3,-3}));
   part->AddDecay(Particle::Decay(0, 0.000259,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{6,-6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{21,21}));

   // Creating rho_tech0
   new Particle("rho_tech0", 54, 0, "Unknown", 100, 0, 210, 0.82101, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(54));
   part->AddDecay(Particle::Decay(0, 0.342802,  vector<int>{24,-52}));
   part->AddDecay(Particle::Decay(0, 0.342802,  vector<int>{52,-24}));
   part->AddDecay(Particle::Decay(0, 0.153373,  vector<int>{24,-24}));
   part->AddDecay(Particle::Decay(0, 0.0868672,  vector<int>{22,51}));
   part->AddDecay(Particle::Decay(0, 0.0312801,  vector<int>{22,53}));
   part->AddDecay(Particle::Decay(32, 0.00691101,  vector<int>{2,-2}));
   part->AddDecay(Particle::Decay(32, 0.00691101,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(32, 0.00478901,  vector<int>{3,-3}));
   part->AddDecay(Particle::Decay(32, 0.00478901,  vector<int>{1,-1}));
   part->AddDecay(Particle::Decay(32, 0.00478901,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(0, 0.00307701,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(0, 0.00307701,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(0, 0.00307701,  vector<int>{15,-15}));
   part->AddDecay(Particle::Decay(0, 0.001598,  vector<int>{23,51}));
   part->AddDecay(Particle::Decay(0, 0.00103,  vector<int>{14,-14}));
   part->AddDecay(Particle::Decay(0, 0.00103,  vector<int>{12,-12}));
   part->AddDecay(Particle::Decay(0, 0.00103,  vector<int>{16,-16}));
   part->AddDecay(Particle::Decay(0, 0.000768002,  vector<int>{23,53}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-7}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{8,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{52,-52}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{6,-6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{17,-17}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{18,-18}));

   // Creating rho_tech+
   new Particle("rho_tech+", 55, 1, "Unknown", 100, 1, 210, 0.64973, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(55));
   part->AddDecay(Particle::Decay(0, 0.474101,  vector<int>{24,51}));
   part->AddDecay(Particle::Decay(0, 0.176299,  vector<int>{52,23}));
   part->AddDecay(Particle::Decay(0, 0.138845,  vector<int>{24,23}));
   part->AddDecay(Particle::Decay(0, 0.109767,  vector<int>{52,22}));
   part->AddDecay(Particle::Decay(32, 0.0285839,  vector<int>{-1,2}));
   part->AddDecay(Particle::Decay(32, 0.0285299,  vector<int>{-3,4}));
   part->AddDecay(Particle::Decay(0, 0.00966098,  vector<int>{-11,12}));
   part->AddDecay(Particle::Decay(0, 0.00966098,  vector<int>{-13,14}));
   part->AddDecay(Particle::Decay(0, 0.00965998,  vector<int>{-15,16}));
   part->AddDecay(Particle::Decay(0, 0.00816098,  vector<int>{24,53}));
   part->AddDecay(Particle::Decay(32, 0.00373499,  vector<int>{-5,6}));
   part->AddDecay(Particle::Decay(32, 0.001468,  vector<int>{-3,2}));
   part->AddDecay(Particle::Decay(32, 0.001468,  vector<int>{-1,4}));
   part->AddDecay(Particle::Decay(32, 5.29999e-05,  vector<int>{-5,4}));
   part->AddDecay(Particle::Decay(32, 6.99999e-06,  vector<int>{-3,6}));
   part->AddDecay(Particle::Decay(32, 9.99998e-07,  vector<int>{-5,2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-1,8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-5,8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,2}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,4}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,6}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-7,8}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-3,8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{52,51}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{-1,6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-17,18}));

   // Creating omega_tech
   new Particle("omega_tech", 56, 0, "Unknown", 100, 0, 210, 0.1575, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(56));
   part->AddDecay(Particle::Decay(0, 0.45294,  vector<int>{22,53}));
   part->AddDecay(Particle::Decay(0, 0.163019,  vector<int>{22,51}));
   part->AddDecay(Particle::Decay(32, 0.045908,  vector<int>{2,-2}));
   part->AddDecay(Particle::Decay(32, 0.045908,  vector<int>{4,-4}));
   part->AddDecay(Particle::Decay(0, 0.038354,  vector<int>{11,-11}));
   part->AddDecay(Particle::Decay(0, 0.038354,  vector<int>{13,-13}));
   part->AddDecay(Particle::Decay(0, 0.038354,  vector<int>{15,-15}));
   part->AddDecay(Particle::Decay(0, 0.038042,  vector<int>{52,-24}));
   part->AddDecay(Particle::Decay(0, 0.038042,  vector<int>{24,-52}));
   part->AddDecay(Particle::Decay(32, 0.017733,  vector<int>{3,-3}));
   part->AddDecay(Particle::Decay(32, 0.017733,  vector<int>{1,-1}));
   part->AddDecay(Particle::Decay(32, 0.017733,  vector<int>{5,-5}));
   part->AddDecay(Particle::Decay(0, 0.011181,  vector<int>{14,-14}));
   part->AddDecay(Particle::Decay(0, 0.011181,  vector<int>{12,-12}));
   part->AddDecay(Particle::Decay(0, 0.011181,  vector<int>{16,-16}));
   part->AddDecay(Particle::Decay(0, 0.00833401,  vector<int>{23,53}));
   part->AddDecay(Particle::Decay(0, 0.004003,  vector<int>{23,51}));
   part->AddDecay(Particle::Decay(0, 0.001999,  vector<int>{24,-24}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{7,-7}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{8,-8}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{52,-52}));
   part->AddDecay(Particle::Decay(32, 0,  vector<int>{6,-6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{17,-17}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{18,-18}));

   // Creating H_L++
   new Particle("H_L++", 61, 1, "Unknown", 100, 2, 200, 0.88161, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(61));
   part->AddDecay(Particle::Decay(0, 0.812251,  vector<int>{-15,-15}));
   part->AddDecay(Particle::Decay(0, 0.0902641,  vector<int>{-13,-13}));
   part->AddDecay(Particle::Decay(0, 0.0902641,  vector<int>{-11,-11}));
   part->AddDecay(Particle::Decay(0, 0.001806,  vector<int>{24,24}));
   part->AddDecay(Particle::Decay(0, 0.001805,  vector<int>{-13,-15}));
   part->AddDecay(Particle::Decay(0, 0.001805,  vector<int>{-11,-15}));
   part->AddDecay(Particle::Decay(0, 0.001805,  vector<int>{-11,-13}));

   // Creating H_R++
   new Particle("H_R++", 62, 1, "Unknown", 100, 2, 200, 0.88001, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(62));
   part->AddDecay(Particle::Decay(0, 0.813719,  vector<int>{-15,-15}));
   part->AddDecay(Particle::Decay(0, 0.0904279,  vector<int>{-13,-13}));
   part->AddDecay(Particle::Decay(0, 0.0904279,  vector<int>{-11,-11}));
   part->AddDecay(Particle::Decay(0, 0.001809,  vector<int>{-11,-13}));
   part->AddDecay(Particle::Decay(0, 0.001808,  vector<int>{-13,-15}));
   part->AddDecay(Particle::Decay(0, 0.001808,  vector<int>{-11,-15}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{63,63}));

   // Creating W_R+
   new Particle("W_R+", 63, 1, "Unknown", 100, 1, 750, 19.3391, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(63));
   part->AddDecay(Particle::Decay(32, 0.325914,  vector<int>{-1,2}));
   part->AddDecay(Particle::Decay(32, 0.32532,  vector<int>{-3,4}));
   part->AddDecay(Particle::Decay(32, 0.314118,  vector<int>{-5,6}));
   part->AddDecay(Particle::Decay(32, 0.016736,  vector<int>{-3,2}));
   part->AddDecay(Particle::Decay(32, 0.016735,  vector<int>{-1,4}));
   part->AddDecay(Particle::Decay(32, 0.000603001,  vector<int>{-5,4}));
   part->AddDecay(Particle::Decay(32, 0.000554001,  vector<int>{-3,6}));
   part->AddDecay(Particle::Decay(32, 1e-05,  vector<int>{-5,2}));
   part->AddDecay(Particle::Decay(32, 9.00001e-06,  vector<int>{-1,6}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-11,64}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-13,65}));
   part->AddDecay(Particle::Decay(0, 0,  vector<int>{-15,66}));

   // Creating nu_Re
   new Particle("nu_Re", 64, 1, "Unknown", 100, 0, 750, 0, -100, -1, -100, -1, -1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
