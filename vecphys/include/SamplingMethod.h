#ifndef SamplingMethod_H
#define SamplingMethod_H

namespace vecphys {

enum SamplingMethod {
  kNullMethod = -1,
  kAlias,                //Alias method
  kRejection,            //Geant4 CompositionAndRejection
  kUnpack,               //Shuffling
  kNumberSamplingMethod
};

} // end namespace vecphys

#endif
