#ifndef GUPhysicsModelName_H
#define GUPhysicsModelName_H

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

enum GUPhysicsModelIndex {
  kNullModel = -1,
  kKleinNishina,      // Compton
  kHybridKleinNishina,// Compton with Alias + Composition and Rejection
  kBetheHeitler,      // Conversion
  kMollerBhabha,      // Ionization
  kSauterGavrila,     // Photo-Electric Effect
  kSeltzerBerger,     // Bremsstrahlung
  kNumberPhysicsModel
};

static const char* GUPhysicsModelName[kNumberPhysicsModel] = {
  "KleinNishina",
  "HybridKleinNishina",
  "BetheHeitler",
  "MollerBhabha",
  "SauterGavrila",
  "SeltzerBerger"
};

}
} // end namespace vecphys

#endif
