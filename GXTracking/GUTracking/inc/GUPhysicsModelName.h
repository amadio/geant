#ifndef GUPhysicsModelName_H
#define GUPhysicsModelName_H 

namespace vecphys {

enum GUPhysicsModelIndex {
  kNullModel = -1, 
  kKleinNishina,      // Compton  
  kVKleinNishina,     // Compton inherited from GUEmModelBase  
  kBetheHeitler,      // Conversion
  kSauterGavrila,     // Photo-Electric Effect
  kMollerBhabha,      // Ionization
  kSeltzerBerger,     // Bremsstrahlung 
  kNumberPhysicsModel
};

static 
const char* GUPhysicsModelName[kNumberPhysicsModel] = {
  "KleinNishina ",
  "VKleinNishina",
  "BetheHeitler ",
  "SauterGavrila",
  "MollerBhabha ",
  "SeltzerBerger"
};

} // end namespace vecphys

#endif
