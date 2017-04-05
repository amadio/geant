#ifndef GUPhysicsModelName_H
#define GUPhysicsModelName_H

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

enum GUPhysicsModelIndex {
  kNullModel = -1,
  kKleinNishina,       // Compton
  kHybridKleinNishina, // Compton with Alias + Composition and Rejection
  kBetheHeitler,       // Conversion
  kSauterGavrila,      // Photo-Electric Effect
  kMollerBhabha,       // Ionization
  kSeltzerBerger,      // Bremsstrahlung
  kNumberPhysicsModel
};

static const char *GUPhysicsModelName[kNumberPhysicsModel] = {"KleinNishina",  "HybridKleinNishina", "BetheHeitler",
                                                              "SauterGavrila", "MollerBhabha",       "SeltzerBerger"};
}
} // end namespace vecphys

#endif
