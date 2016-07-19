#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
using std::vector;
namespace geant {
   inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
void CreateParticle0001() {

   // Creating ~nu_eL_bar
   new Particle("~nu_eL_bar", -1000012, 0, "Sparticle", 100, 0, 500, 1, 100, 100, 1, 100, 1);
   Particle *part = 0;
   part = const_cast<Particle*>(&Particle::Particles().at(-1000012));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000011,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000011,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000011,-37}));

   // Creating ~e_L+
   new Particle("~e_L+", -1000011, 0, "Sparticle", 100, 1, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000011));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-12}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-11}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000012,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000012,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000012,37}));

   // Creating ~t_1_bar
   new Particle("~t_1_bar", -1000006, 0, "Sparticle", 100, -0.666667, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000006));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000005,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000005,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000016,15,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000015,-16,-5}));

   // Creating ~b_1_bar
   new Particle("~b_1_bar", -1000005, 0, "Sparticle", 100, 0.333333, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000005));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-6}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-5}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000006,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000006,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-5}));

   // Creating ~c_L_bar
   new Particle("~c_L_bar", -1000004, 0, "Sparticle", 100, -0.666667, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000004));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000003,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000003,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-4}));

   // Creating ~s_L_bar
   new Particle("~s_L_bar", -1000003, 0, "Sparticle", 100, 0.333333, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000003));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-4}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-3}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000004,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000004,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-3}));

   // Creating ~u_L_bar
   new Particle("~u_L_bar", -1000002, 0, "Sparticle", 100, -0.666667, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000002));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000024,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000037,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,-24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000001,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000001,-37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-2}));

   // Creating ~d_L_bar
   new Particle("~d_L_bar", -1000001, 0, "Sparticle", 100, 0.333333, 500, 1, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-1000001));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000039,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000024,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{1000037,-2}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000022,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000023,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000025,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000035,-1}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,24}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000002,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-2000002,37}));
   part->AddDecay(Particle::Decay(53, 0,  vector<int>{-1000021,-1}));

   // Creating phi(1680)_bar
   new Particle("phi(1680)_bar", -100333, 0, "Unknown", 100, 0, 1.68, 0.15, 100, 100, 0, 100, 1);

   // Creating eta(1475)_bar
   new Particle("eta(1475)_bar", -100331, 0, "Unknown", 100, 0, 1.476, 0.085, 100, 100, 0, 100, 1);

   // Creating k2_star(1980)-_bar
   new Particle("k2_star(1980)-_bar", -100325, 0, "Unknown", 100, -1, 1.973, 0.373, 100, 100, 0, 100, 1);

   // Creating k_star(1410)-_bar
   new Particle("k_star(1410)-_bar", -100323, 0, "Unknown", 100, -1, 1.414, 0.232, 100, 100, 0, 100, 1);

   // Creating k(1460)-_bar
   new Particle("k(1460)-_bar", -100321, 0, "Unknown", 100, -1, 1.46, 0.26, 100, 100, 0, 100, 1);

   // Creating k2_star(1980)0_bar
   new Particle("k2_star(1980)0_bar", -100315, 0, "Unknown", 100, 0, 1.973, 0.373, 100, 100, 0, 100, 1);

   // Creating k_star(1410)0_bar
   new Particle("k_star(1410)0_bar", -100313, 0, "Unknown", 100, 0, 1.414, 0.232, 100, 100, 0, 100, 1);

   // Creating k(1460)0_bar
   new Particle("k(1460)0_bar", -100311, 0, "Unknown", 100, 0, 1.46, 0.26, 100, 100, 0, 100, 1);

   // Creating omega(1420)_bar
   new Particle("omega(1420)_bar", -100223, 0, "Unknown", 100, 0, 1.425, 0.215, 100, 100, 0, 100, 1);

   // Creating eta(1295)_bar
   new Particle("eta(1295)_bar", -100221, 0, "Unknown", 100, 0, 1.294, 0.055, 100, 100, 0, 100, 1);

   // Creating rho(1450)-_bar
   new Particle("rho(1450)-_bar", -100213, 0, "Unknown", 100, -1, 1.465, 0.4, 100, 100, 0, 100, 1);

   // Creating pi(1300)-_bar
   new Particle("pi(1300)-_bar", -100211, 0, "Unknown", 100, -1, 1.3, 0.4, 100, 100, 0, 100, 1);

   // Creating rho(1450)0_bar
   new Particle("rho(1450)0_bar", -100113, 0, "Unknown", 100, 0, 1.465, 0.4, 100, 100, 0, 100, 1);

   // Creating pi(1300)0_bar
   new Particle("pi(1300)0_bar", -100111, 0, "Unknown", 100, 0, 1.3, 0.4, 100, 100, 0, 100, 1);

   // Creating lambda(1810)_bar
   new Particle("lambda(1810)_bar", -53122, 0, "Unknown", 100, 0, 1.81, 0.15, 100, 100, 0, 100, 1);

   // Creating N(2090)+_bar
   new Particle("N(2090)+_bar", -52214, 0, "Unknown", 100, -1, 2.08, 0.35, 100, 100, 0, 100, 1);

   // Creating N(2090)0_bar
   new Particle("N(2090)0_bar", -52114, 0, "Unknown", 100, 0, 2.08, 0.35, 100, 100, 0, 100, 1);

   // Creating lambda(1800)_bar
   new Particle("lambda(1800)_bar", -43122, 0, "Unknown", 100, 0, 1.8, 0.3, 100, 100, 0, 100, 1);

   // Creating N(1710)+_bar
   new Particle("N(1710)+_bar", -42212, 0, "Unknown", 100, -1, 1.71, 0.1, 100, 100, 0, 100, 1);

   // Creating N(1900)+_bar
   new Particle("N(1900)+_bar", -42124, 0, "Unknown", 100, -1, 1.9, 0.5, 100, 100, 0, 100, 1);

   // Creating N(1710)0_bar
   new Particle("N(1710)0_bar", -42112, 0, "Unknown", 100, 0, 1.71, 0.1, 100, 100, 0, 100, 1);

   // Creating N(1900)0_bar
   new Particle("N(1900)0_bar", -41214, 0, "Unknown", 100, 0, 1.9, 0.5, 100, 100, 0, 100, 1);

   // Creating xi(1950)0_bar
   new Particle("xi(1950)0_bar", -33324, 0, "Unknown", 100, 0, 1.95, 0.06, 100, 100, 0, 100, 1);

   // Creating xi(1950)-_bar
   new Particle("xi(1950)-_bar", -33314, 0, "Unknown", 100, 1, 1.95, 0.06, 100, 100, 0, 100, 1);

   // Creating lambda(1670)_bar
   new Particle("lambda(1670)_bar", -33122, 0, "Unknown", 100, 0, 1.67, 0.035, 100, 100, 0, 100, 1);

   // Creating delta(1600)++_bar
   new Particle("delta(1600)++_bar", -32224, 0, "Unknown", 100, -2, 1.6, 0.35, 100, 100, 0, 100, 1);

   // Creating delta(1600)+_bar
   new Particle("delta(1600)+_bar", -32214, 0, "Unknown", 100, -1, 1.6, 0.35, 100, 100, 0, 100, 1);

   // Creating N(1650)+_bar
   new Particle("N(1650)+_bar", -32212, 0, "Unknown", 100, -1, 1.655, 0.165, 100, 100, 0, 100, 1);

   // Creating N(1720)+_bar
   new Particle("N(1720)+_bar", -32124, 0, "Unknown", 100, -1, 1.72, 0.2, 100, 100, 0, 100, 1);

   // Creating delta(1600)0_bar
   new Particle("delta(1600)0_bar", -32114, 0, "Unknown", 100, 0, 1.6, 0.35, 100, 100, 0, 100, 1);

   // Creating N(1650)0_bar
   new Particle("N(1650)0_bar", -32112, 0, "Unknown", 100, 0, 1.655, 0.165, 100, 100, 0, 100, 1);

   // Creating N(1720)0_bar
   new Particle("N(1720)0_bar", -31214, 0, "Unknown", 100, 0, 1.72, 0.2, 100, 100, 0, 100, 1);
}

 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
