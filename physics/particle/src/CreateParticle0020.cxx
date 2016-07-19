#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
using std::vector;
namespace geant {
   inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
void CreateParticle0020() {

   // Creating ~mu_R-
   new Particle("~mu_R-", 2000013, 1, "Sparticle", 100, -1, 500, 1, -100, -1, -100, -1, -1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(2000013));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,14}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,13}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000013,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000014,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000014,-37}));

   // Creating ~nu_muR
   new Particle("~nu_muR", 2000014, 1, "Sparticle", 100, 0, 500, 0, -100, -1, -100, -1, -1);

   // Creating ~tau_2-
   new Particle("~tau_2-", 2000015, 1, "Sparticle", 100, -1, 500, 1, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(2000015));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000039,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,16}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000022,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000023,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000025,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000035,15}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,23}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,25}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,35}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,36}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000016,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{2000016,-37}));

   // Creating ~nu_tauR
   new Particle("~nu_tauR", 2000016, 1, "Sparticle", 100, 0, 500, 0, -100, -1, -100, -1, -1);

   // Creating d*
   new Particle("d*", 4000001, 1, "Excited", 100, -0.333333, 400, 2.65171, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4000001));
   part->AddDecay(Particle::Decay(53, 0.85422,  vector<int>{21,1}));
   part->AddDecay(Particle::Decay(0, 0.096449,  vector<int>{-24,2}));
   part->AddDecay(Particle::Decay(0, 0.044039,  vector<int>{23,1}));
   part->AddDecay(Particle::Decay(0, 0.005292,  vector<int>{22,1}));

   // Creating u*
   new Particle("u*", 4000002, 1, "Excited", 100, 0.666667, 400, 2.65499, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4000002));
   part->AddDecay(Particle::Decay(0, 0.853166,  vector<int>{21,2}));
   part->AddDecay(Particle::Decay(0, 0.0963291,  vector<int>{24,1}));
   part->AddDecay(Particle::Decay(0, 0.029361,  vector<int>{23,2}));
   part->AddDecay(Particle::Decay(0, 0.021144,  vector<int>{22,2}));

   // Creating e*-
   new Particle("e*-", 4000011, 1, "Excited", 100, -1, 400, 0.42901, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4000011));
   part->AddDecay(Particle::Decay(0, 0.596149,  vector<int>{-24,12}));
   part->AddDecay(Particle::Decay(0, 0.294414,  vector<int>{22,11}));
   part->AddDecay(Particle::Decay(0, 0.109437,  vector<int>{23,11}));

   // Creating nu*_e0
   new Particle("nu*_e0", 4000012, 1, "Excited", 100, 0, 400, 0.41917, -100, -1, -100, -1, -1);
   part = const_cast<Particle*>(&Particle::Particles().at(4000012));
   part->AddDecay(Particle::Decay(0, 0.610139,  vector<int>{24,11}));
   part->AddDecay(Particle::Decay(0, 0.389861,  vector<int>{23,12}));

   // Creating a0(980)0
   new Particle("a0(980)0", 9000111, 1, "Unknown", 100, 0, 0.98, 0.075, -100, 0, -100, -1, -1);

   // Creating a0(980)+
   new Particle("a0(980)+", 9000211, 1, "Unknown", 100, 1, 0.98, 0.06, -100, 0, -100, -1, -1);

   // Creating f0(600)
   new Particle("f0(600)", 9000221, 1, "Unknown", 100, 0, 0.8, 0.8, -100, 0, -100, -1, -1);

   // Creating f0(980)
   new Particle("f0(980)", 9010221, 1, "Unknown", 100, 0, 0.98, 0.07, -100, 0, -100, -1, -1);

   // Creating eta(1405)
   new Particle("eta(1405)", 9020221, 1, "Unknown", 100, 0, 1.4098, 0.0511, -100, 0, -100, -1, -1);

   // Creating f0(1500)
   new Particle("f0(1500)", 9030221, 1, "Unknown", 100, 0, 1.505, 0.109, -100, 0, -100, -1, -1);

   // Creating f2(1810)
   new Particle("f2(1810)", 9030225, 1, "Unknown", 100, 0, 1.815, 0.197, -100, 0, -100, -1, -1);

   // Creating f2(2010)
   new Particle("f2(2010)", 9060225, 1, "Unknown", 100, 0, 2.01, 0.2, -100, 0, -100, -1, -1);

   // Creating Cherenkov
   new Particle("Cherenkov", 50000050, 1, "Unknown", 100, 0, 0, 0, -100, 0, -100, -1, -1);

   // Creating ChargedRootino
   new Particle("ChargedRootino", 50000052, 1, "Unknown", 100, -0.333333, 0, 0, -100, 0, -100, -1, -1);

   // Creating GenericIon
   new Particle("GenericIon", 50000060, 1, "Unknown", 100, 0.333333, 0.938272, 0, -100, 0, -100, -1, -1);

   // Creating N(2220)0
   new Particle("N(2220)0", 100002110, 1, "Unknown", 100, 0, 2.25, 0.4, -100, 0, -100, -1, -1);

   // Creating N(2220)+
   new Particle("N(2220)+", 100002210, 1, "Unknown", 100, 1, 2.25, 0.4, -100, 0, -100, -1, -1);

   // Creating N(2250)0
   new Particle("N(2250)0", 100012110, 1, "Unknown", 100, 0, 2.275, 0.5, -100, 0, -100, -1, -1);

   // Creating N(2250)+
   new Particle("N(2250)+", 100012210, 1, "Unknown", 100, 1, 2.275, 0.5, -100, 0, -100, -1, -1);

   // Creating Deuteron
   new Particle("Deuteron", 1000010020, 1, "ion", 100, 1, 1.87106, 0, -100, 0, -100, -1, -1);

   // Creating Triton
   new Particle("Triton", 1000010030, 1, "ion", 100, 1, 2.80941, 1.6916e-33, -100, 0, -100, -1, -1);

   // Creating HE3
   new Particle("HE3", 1000020030, 1, "ion", 100, 2, 2.80941, 0, -100, 0, -100, -1, -1);

   // Creating Alpha
   new Particle("Alpha", 1000020040, 1, "ion", 100, 2, 3.7284, 1.6916e-33, -100, 0, -100, -1, -1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
