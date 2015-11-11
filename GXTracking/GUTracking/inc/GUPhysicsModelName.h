#ifndef GUPhysicsModelName_H
#define GUPhysicsModelName_H 

namespace vecphys {

enum GUPhysicsModelIndex {
  kNullModel = -1, 
  kKleinNishina,      // Compton  
  kBetheHeitler,      // Conversion
  kSauterGavrila,     // Photo-Electric Effect
  kMollerBhabha,      // Ionization
  kSeltzerBerger,     // Bremsstrahlung 
  kNumberPhysicsModel
};

static 
const char* GUPhysicsModelName[kNumberPhysicsModel] = {
  "KleinNishina",
  "BetheHeitler",
  "SauterGavrila",
  "MollerBhabha",
  "SeltzerBerger"
};

} // end namespace vecphys

#endif
