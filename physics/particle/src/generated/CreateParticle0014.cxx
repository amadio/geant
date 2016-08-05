// This files was autogenerated by geant::Particle::ReadFile

#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize off
#endif
#include "Particle.h"
#ifdef VECCORE_CUDA
#include "base/Vector.h"
template <typename T>
using vector = vecgeom::Vector<T>;
#else
using std::vector;
#endif

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {


//________________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void CreateParticle0014() {
   vector<int> daughters;
   Particle *part = nullptr;

   // Creating Sigma*_c0_bar
   new Particle("Sigma*_c0_bar", -4114, 0, "CharmedBaryon", 100, 0, 2.518, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4114));
   daughters.clear();
   daughters.push_back(-4122);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 1,  daughters));

   // Creating Sigma_c0_bar
   new Particle("Sigma_c0_bar", -4112, 0, "CharmedBaryon", 100, 0, 2.45376, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-4112));
   daughters.clear();
   daughters.push_back(-4122);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 1,  daughters));

   // Creating cd_1_bar
   new Particle("cd_1_bar", -4103, 0, "Unknown", 100, -0.333333, 2.00808, 0, 100, 100, 1, 100, 1);

   // Creating cd_0_bar
   new Particle("cd_0_bar", -4101, 0, "Unknown", 100, -0.333333, 1.96908, 0, 100, 100, 1, 100, 1);

   // Creating Omega+
   new Particle("Omega+", -3334, 0, "Baryon", 100, 1, 1.67245, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3334));
   daughters.clear();
   daughters.push_back(-3122);
   daughters.push_back(321);
   part->AddDecay(Particle::Decay(0, 0.676,  daughters));
   daughters.clear();
   daughters.push_back(-3322);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.234,  daughters));
   daughters.clear();
   daughters.push_back(-3312);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.085,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-3322);
   part->AddDecay(Particle::Decay(0, 0.005,  daughters));

   // Creating Xi*0_bar
   new Particle("Xi*0_bar", -3324, 0, "Baryon", 100, 0, 1.5318, 0.0091, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3324));
   daughters.clear();
   daughters.push_back(-3312);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(-3322);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating Xi0_bar
   new Particle("Xi0_bar", -3322, 0, "Baryon", 100, 0, 1.31486, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3322));
   daughters.clear();
   daughters.push_back(-3122);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.9954,  daughters));
   daughters.clear();
   daughters.push_back(-3212);
   daughters.push_back(-22);
   part->AddDecay(Particle::Decay(0, 0.0035,  daughters));
   daughters.clear();
   daughters.push_back(-3122);
   daughters.push_back(-22);
   part->AddDecay(Particle::Decay(0, 0.0011,  daughters));

   // Creating Xi*+
   new Particle("Xi*+", -3314, 0, "Baryon", 100, 1, 1.535, 0.0099, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3314));
   daughters.clear();
   daughters.push_back(-3322);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.667,  daughters));
   daughters.clear();
   daughters.push_back(-3312);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.333,  daughters));

   // Creating Xi-_bar
   new Particle("Xi-_bar", -3312, 0, "Baryon", 100, 1, 1.32171, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3312));
   daughters.clear();
   daughters.push_back(-3122);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.9988,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-3122);
   part->AddDecay(Particle::Decay(0, 0.0006,  daughters));
   daughters.clear();
   daughters.push_back(14);
   daughters.push_back(-13);
   daughters.push_back(-3122);
   part->AddDecay(Particle::Decay(0, 0.0004,  daughters));
   daughters.clear();
   daughters.push_back(-3112);
   daughters.push_back(-22);
   part->AddDecay(Particle::Decay(0, 0.0001,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-3212);
   part->AddDecay(Particle::Decay(0, 0.0001,  daughters));

   // Creating ss_1_bar
   new Particle("ss_1_bar", -3303, 0, "Unknown", 100, 0.666667, 2.08, 0, 100, 100, 1, 100, 1);

   // Creating sigma(2030)+_bar
   new Particle("sigma(2030)+_bar", -3228, 0, "Unknown", 100, -1, 2.03, 0.18, 100, 100, 0, 100, 1);

   // Creating sigma(1775)+_bar
   new Particle("sigma(1775)+_bar", -3226, 0, "Unknown", 100, -1, 1.775, 0.12, 100, 100, 0, 100, 1);

   // Creating Sigma*+_bar
   new Particle("Sigma*+_bar", -3224, 0, "Baryon", 100, -1, 1.3828, 0.0358, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3224));
   daughters.clear();
   daughters.push_back(-3122);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.88,  daughters));
   daughters.clear();
   daughters.push_back(-3222);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.06,  daughters));
   daughters.clear();
   daughters.push_back(-3212);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.06,  daughters));

   // Creating Sigma+_bar
   new Particle("Sigma+_bar", -3222, 0, "Baryon", 100, -1, 1.18937, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3222));
   daughters.clear();
   daughters.push_back(-2212);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.516,  daughters));
   daughters.clear();
   daughters.push_back(-2112);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.483,  daughters));
   daughters.clear();
   daughters.push_back(-2212);
   daughters.push_back(-22);
   part->AddDecay(Particle::Decay(0, 0.001,  daughters));

   // Creating sigma(2030)0_bar
   new Particle("sigma(2030)0_bar", -3218, 0, "Unknown", 100, 0, 2.03, 0.18, 100, 100, 0, 100, 1);

   // Creating sigma(1775)0_bar
   new Particle("sigma(1775)0_bar", -3216, 0, "Unknown", 100, 0, 1.775, 0.12, 100, 100, 0, 100, 1);

   // Creating Sigma*0_bar
   new Particle("Sigma*0_bar", -3214, 0, "Baryon", 100, 0, 1.3837, 0.036, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3214));
   daughters.clear();
   daughters.push_back(-3122);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.88,  daughters));
   daughters.clear();
   daughters.push_back(-3222);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.06,  daughters));
   daughters.clear();
   daughters.push_back(-3112);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.06,  daughters));

   // Creating Sigma0_bar
   new Particle("Sigma0_bar", -3212, 0, "Baryon", 100, 0, 1.19264, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3212));
   daughters.clear();
   daughters.push_back(-3122);
   daughters.push_back(-22);
   part->AddDecay(Particle::Decay(0, 1,  daughters));

   // Creating su_1_bar
   new Particle("su_1_bar", -3203, 0, "Unknown", 100, -0.333333, 0.1064, 0, 100, 100, 1, 100, 1);

   // Creating su_0_bar
   new Particle("su_0_bar", -3201, 0, "Unknown", 100, -0.333333, 0.1064, 0, 100, 100, 1, 100, 1);

   // Creating lambda(2100)_bar
   new Particle("lambda(2100)_bar", -3128, 0, "Unknown", 100, 0, 2.1, 0.2, 100, 100, 0, 100, 1);

   // Creating lambda(1820)_bar
   new Particle("lambda(1820)_bar", -3126, 0, "Unknown", 100, 0, 1.82, 0.08, 100, 100, 0, 100, 1);

   // Creating lambda(1520)_bar
   new Particle("lambda(1520)_bar", -3124, 0, "Unknown", 100, 0, 1.5195, 0.0156, 100, 100, 0, 100, 1);

   // Creating Lambda0_bar
   new Particle("Lambda0_bar", -3122, 0, "Baryon", 100, 0, 1.11568, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3122));
   daughters.clear();
   daughters.push_back(-2212);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.639,  daughters));
   daughters.clear();
   daughters.push_back(-2112);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.358,  daughters));
   daughters.clear();
   daughters.push_back(-2112);
   daughters.push_back(-22);
   part->AddDecay(Particle::Decay(0, 0.002,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-2212);
   part->AddDecay(Particle::Decay(0, 0.001,  daughters));

   // Creating sigma(2030)-_bar
   new Particle("sigma(2030)-_bar", -3118, 0, "Unknown", 100, 1, 2.03, 0.18, 100, 100, 0, 100, 1);

   // Creating sigma(1775)-_bar
   new Particle("sigma(1775)-_bar", -3116, 0, "Unknown", 100, 1, 1.775, 0.12, 100, 100, 0, 100, 1);

   // Creating Sigma*-_bar
   new Particle("Sigma*-_bar", -3114, 0, "Baryon", 100, 1, 1.3872, 0.0394, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3114));
   daughters.clear();
   daughters.push_back(-3122);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.88,  daughters));
   daughters.clear();
   daughters.push_back(-3212);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.06,  daughters));
   daughters.clear();
   daughters.push_back(-3112);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.06,  daughters));

   // Creating Sigma-_bar
   new Particle("Sigma-_bar", -3112, 0, "Baryon", 100, 1, 1.19744, 0, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-3112));
   daughters.clear();
   daughters.push_back(-2112);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.999,  daughters));
   daughters.clear();
   daughters.push_back(12);
   daughters.push_back(-11);
   daughters.push_back(-2112);
   part->AddDecay(Particle::Decay(0, 0.001,  daughters));

   // Creating sd_1_bar
   new Particle("sd_1_bar", -3103, 0, "Unknown", 100, 0.666667, 0.1088, 0, 100, 100, 1, 100, 1);

   // Creating sd_0_bar
   new Particle("sd_0_bar", -3101, 0, "Unknown", 100, 0.666667, 0.108, 0, 100, 100, 1, 100, 1);

   // Creating delta(1950)++_bar
   new Particle("delta(1950)++_bar", -2228, 0, "Unknown", 100, -2, 1.93, 0.28, 100, 100, 0, 100, 1);

   // Creating delta(1905)++_bar
   new Particle("delta(1905)++_bar", -2226, 0, "Unknown", 100, -2, 1.89, 0.33, 100, 100, 0, 100, 1);

   // Creating Delta--
   new Particle("Delta--", -2224, 0, "Baryon", 100, -2, 1.232, 0.12, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2224));
   daughters.clear();
   daughters.push_back(-2212);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 1,  daughters));

   // Creating delta(1620)++_bar
   new Particle("delta(1620)++_bar", -2222, 0, "Unknown", 100, -2, 1.63, 0.145, 100, 100, 0, 100, 1);

   // Creating delta(1950)+_bar
   new Particle("delta(1950)+_bar", -2218, 0, "Unknown", 100, -1, 1.93, 0.28, 100, 100, 0, 100, 1);

   // Creating N(1675)+_bar
   new Particle("N(1675)+_bar", -2216, 0, "Unknown", 100, -1, 1.675, 0.15, 100, 100, 0, 100, 1);

   // Creating Delta+_bar
   new Particle("Delta+_bar", -2214, 0, "Baryon", 100, -1, 1.232, 0.12, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-2214));
   daughters.clear();
   daughters.push_back(-2212);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.663,  daughters));
   daughters.clear();
   daughters.push_back(-2112);
   daughters.push_back(-211);
   part->AddDecay(Particle::Decay(0, 0.331,  daughters));
   daughters.clear();
   daughters.push_back(-2212);
   daughters.push_back(-22);
   part->AddDecay(Particle::Decay(0, 0.006,  daughters));

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
   daughters.clear();
   daughters.push_back(-2112);
   daughters.push_back(-111);
   part->AddDecay(Particle::Decay(0, 0.663,  daughters));
   daughters.clear();
   daughters.push_back(-2212);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 0.331,  daughters));
   daughters.clear();
   daughters.push_back(-2112);
   daughters.push_back(-22);
   part->AddDecay(Particle::Decay(0, 0.006,  daughters));

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
   daughters.clear();
   daughters.push_back(-2112);
   daughters.push_back(211);
   part->AddDecay(Particle::Decay(0, 1,  daughters));

   // Creating delta(1620)-_bar
   new Particle("delta(1620)-_bar", -1112, 0, "Unknown", 100, 1, 1.63, 0.145, 100, 100, 0, 100, 1);

   // Creating dd_1_bar
   new Particle("dd_1_bar", -1103, 0, "Unknown", 100, 0.666667, 0.96, 0, 100, 100, 1, 100, 1);

   // Creating B*_2c-
   new Particle("B*_2c-", -545, 0, "B-Meson", 100, -1, 7.35, 0.02, 100, 100, 1, 100, 1);
   part = const_cast<Particle*>(&Particle::Particles().at(-545));
   daughters.clear();
   daughters.push_back(-511);
   daughters.push_back(-411);
   part->AddDecay(Particle::Decay(0, 0.3,  daughters));
   daughters.clear();
   daughters.push_back(-521);
   daughters.push_back(-421);
   part->AddDecay(Particle::Decay(0, 0.3,  daughters));
   daughters.clear();
   daughters.push_back(-513);
   daughters.push_back(-411);
   part->AddDecay(Particle::Decay(0, 0.2,  daughters));
   daughters.clear();
   daughters.push_back(-523);
   daughters.push_back(-421);
   part->AddDecay(Particle::Decay(0, 0.2,  daughters));
}

} // End of inline namespace
} // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
