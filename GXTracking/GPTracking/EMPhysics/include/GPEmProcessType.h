#ifndef GPEmProcessType_h
#define GPEmProcessType_h 1

enum GPEmProcessType {
  kNullType, 
  kCompton,                // ComptonScattering
  kPhotoElectric,          // PhotoElectricEffect
  kConversion,             // GammaConversion
  kBremsstrahlung,         // eBremsstrahlung
  kIonisation,             // eIonisation
  kMultipleScattering,     // eMultipleScattering
  kNumberEmProcess
};

#endif
