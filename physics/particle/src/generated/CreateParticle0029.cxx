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
void CreateParticle0029() {
   vector<int> daughters;
   Particle *part = nullptr;

   // Creating reggeon
   new Particle("reggeon", 28, 0, "GaugeBoson", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating pomeron
   new Particle("pomeron", 29, 0, "GaugeBoson", 100, 0, 0, 0, -100, -1, -100, -1, -1);

   // Creating Z'0
   new Particle("Z'0", 32, 0, "GaugeBoson", 100, 0, 500, 14.5485, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(32));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(32, 0.145869,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-3);
   part->AddDecay(Particle::Decay(32, 0.145869,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-5);
   part->AddDecay(Particle::Decay(32, 0.14581,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0.113303,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 0.113298,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-12);
   part->AddDecay(Particle::Decay(0, 0.0636061,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-14);
   part->AddDecay(Particle::Decay(0, 0.0636061,  daughters));
   daughters.clear();
   daughters.push_back(16);
   daughters.push_back(-16);
   part->AddDecay(Particle::Decay(0, 0.0636061,  daughters));
   daughters.clear();
   daughters.push_back(6);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0.0490131,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-11);
   part->AddDecay(Particle::Decay(0, 0.0320071,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-13);
   part->AddDecay(Particle::Decay(0, 0.0320071,  daughters));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(-15);
   part->AddDecay(Particle::Decay(0, 0.0320041,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-7);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(8);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(17);
   daughters.push_back(-17);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(18);
   daughters.push_back(-18);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(-24);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(37);
   daughters.push_back(-37);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(23);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(23);
   daughters.push_back(25);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(25);
   daughters.push_back(36);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(35);
   daughters.push_back(36);
   part->AddDecay(Particle::Decay(0, 0,  daughters));

   // Creating Z"0
   new Particle("Z\"0", 33, 0, "GaugeBoson", 100, 0, 900, 0, -100, -1, -100, -1, -1);

   // Creating W'+
   new Particle("W'+", 34, 1, "GaugeBoson", 100, 1, 500, 16.6708, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(34));
   daughters.clear();
   daughters.push_back(-1);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(32, 0.251276,  daughters));
   daughters.clear();
   daughters.push_back(-3);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(32, 0.250816,  daughters));
   daughters.clear();
   daughters.push_back(-5);
   daughters.push_back(6);
   part->AddDecay(Particle::Decay(32, 0.215459,  daughters));
   daughters.clear();
   daughters.push_back(-11);
   daughters.push_back(12);
   part->AddDecay(Particle::Decay(0, 0.085262,  daughters));
   daughters.clear();
   daughters.push_back(-13);
   daughters.push_back(14);
   part->AddDecay(Particle::Decay(0, 0.085262,  daughters));
   daughters.clear();
   daughters.push_back(-15);
   daughters.push_back(16);
   part->AddDecay(Particle::Decay(0, 0.08526,  daughters));
   daughters.clear();
   daughters.push_back(-3);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(32, 0.012903,  daughters));
   daughters.clear();
   daughters.push_back(-1);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(32, 0.012903,  daughters));
   daughters.clear();
   daughters.push_back(-5);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(32, 0.000465,  daughters));
   daughters.clear();
   daughters.push_back(-3);
   daughters.push_back(6);
   part->AddDecay(Particle::Decay(32, 0.00038,  daughters));
   daughters.clear();
   daughters.push_back(-5);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(32, 8e-06,  daughters));
   daughters.clear();
   daughters.push_back(-1);
   daughters.push_back(6);
   part->AddDecay(Particle::Decay(32, 6e-06,  daughters));
   daughters.clear();
   daughters.push_back(-7);
   daughters.push_back(2);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(-7);
   daughters.push_back(4);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(-7);
   daughters.push_back(6);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(-7);
   daughters.push_back(8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(-3);
   daughters.push_back(8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1);
   daughters.push_back(8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(-5);
   daughters.push_back(8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(-17);
   daughters.push_back(18);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(23);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(25);
   part->AddDecay(Particle::Decay(0, 0,  daughters));

   // Creating H0
   new Particle("H0", 35, 0, "GaugeBoson", 100, 0, 300, 8.42842, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(35));
   daughters.clear();
   daughters.push_back(24);
   daughters.push_back(-24);
   part->AddDecay(Particle::Decay(0, 0.688641,  daughters));
   daughters.clear();
   daughters.push_back(23);
   daughters.push_back(23);
   part->AddDecay(Particle::Decay(0, 0.306171,  daughters));
   daughters.clear();
   daughters.push_back(25);
   daughters.push_back(25);
   part->AddDecay(Particle::Decay(0, 0.003799,  daughters));
   daughters.clear();
   daughters.push_back(5);
   daughters.push_back(-5);
   part->AddDecay(Particle::Decay(32, 0.000754001,  daughters));
   daughters.clear();
   daughters.push_back(21);
   daughters.push_back(21);
   part->AddDecay(Particle::Decay(0, 0.000439,  daughters));
   daughters.clear();
   daughters.push_back(15);
   daughters.push_back(-15);
   part->AddDecay(Particle::Decay(0, 7.40001e-05,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(23);
   part->AddDecay(Particle::Decay(0, 6.10001e-05,  daughters));
   daughters.clear();
   daughters.push_back(4);
   daughters.push_back(-4);
   part->AddDecay(Particle::Decay(32, 4.6e-05,  daughters));
   daughters.clear();
   daughters.push_back(22);
   daughters.push_back(22);
   part->AddDecay(Particle::Decay(0, 1.5e-05,  daughters));
   daughters.clear();
   daughters.push_back(13);
   daughters.push_back(-13);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(3);
   daughters.push_back(-3);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(17);
   daughters.push_back(-17);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(1);
   daughters.push_back(-1);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(2);
   daughters.push_back(-2);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(6);
   daughters.push_back(-6);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(7);
   daughters.push_back(-7);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(8);
   daughters.push_back(-8);
   part->AddDecay(Particle::Decay(32, 0,  daughters));
   daughters.clear();
   daughters.push_back(23);
   daughters.push_back(25);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(11);
   daughters.push_back(-11);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(36);
   daughters.push_back(36);
   part->AddDecay(Particle::Decay(0, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000022);
   daughters.push_back(1000022);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000023);
   daughters.push_back(1000022);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000023);
   daughters.push_back(1000023);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000025);
   daughters.push_back(1000022);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000025);
   daughters.push_back(1000023);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000025);
   daughters.push_back(1000025);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000035);
   daughters.push_back(1000022);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000035);
   daughters.push_back(1000023);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000035);
   daughters.push_back(1000025);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000035);
   daughters.push_back(1000035);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000024);
   daughters.push_back(-1000024);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000024);
   daughters.push_back(-1000037);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000037);
   daughters.push_back(-1000024);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000037);
   daughters.push_back(-1000037);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000001);
   daughters.push_back(-1000001);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000001);
   daughters.push_back(-2000001);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000001);
   daughters.push_back(-2000001);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000001);
   daughters.push_back(2000001);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000002);
   daughters.push_back(-1000002);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000002);
   daughters.push_back(-2000002);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000002);
   daughters.push_back(-2000002);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000002);
   daughters.push_back(2000002);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000003);
   daughters.push_back(-1000003);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000003);
   daughters.push_back(-2000003);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000003);
   daughters.push_back(-2000003);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000003);
   daughters.push_back(2000003);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000004);
   daughters.push_back(-1000004);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000004);
   daughters.push_back(-2000004);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000004);
   daughters.push_back(-2000004);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000004);
   daughters.push_back(2000004);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000005);
   daughters.push_back(-1000005);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000005);
   daughters.push_back(-2000005);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000005);
   daughters.push_back(-2000005);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000005);
   daughters.push_back(2000005);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000006);
   daughters.push_back(-1000006);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000006);
   daughters.push_back(-2000006);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000006);
   daughters.push_back(-2000006);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000006);
   daughters.push_back(2000006);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000011);
   daughters.push_back(-1000011);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000011);
   daughters.push_back(-2000011);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000011);
   daughters.push_back(-2000011);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000011);
   daughters.push_back(2000011);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000012);
   daughters.push_back(-1000012);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000012);
   daughters.push_back(-2000012);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000012);
   daughters.push_back(-2000012);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000012);
   daughters.push_back(2000012);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000013);
   daughters.push_back(-1000013);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000013);
   daughters.push_back(-2000013);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000013);
   daughters.push_back(-2000013);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000013);
   daughters.push_back(2000013);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000014);
   daughters.push_back(-1000014);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000014);
   daughters.push_back(-2000014);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000014);
   daughters.push_back(-2000014);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000014);
   daughters.push_back(2000014);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000015);
   daughters.push_back(-1000015);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000015);
   daughters.push_back(-2000015);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000015);
   daughters.push_back(-2000015);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000015);
   daughters.push_back(2000015);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000016);
   daughters.push_back(-1000016);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(2000016);
   daughters.push_back(-2000016);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(1000016);
   daughters.push_back(-2000016);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
   daughters.clear();
   daughters.push_back(-1000016);
   daughters.push_back(2000016);
   part->AddDecay(Particle::Decay(53, 0,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif