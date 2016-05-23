#ifndef MODEL_KERNEL_H
#define MODEL_KERNEL_H 1

#include "GUTrack.h"
#include "base/VPGlobal.h"
#include "SamplingMethod.h"

namespace vecphys {

typedef Real_t (*KernelFunc_t)(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                               SamplingMethod sampleType);

// Scalar

Real_t ScalarKleinNishina(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                          SamplingMethod sampleType);

Real_t ScalarHybridCompton(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                           SamplingMethod sampleType);

Real_t ScalarBetheHeitler(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                          SamplingMethod sampleType);

Real_t ScalarSauterGavrila(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                           SamplingMethod sampleType);

Real_t ScalarMollerBhabha(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                          SamplingMethod sampleType);

Real_t ScalarSeltzerBerger(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos,
                           SamplingMethod sampleType);

KernelFunc_t ScalarKernelFunc[] = {ScalarKleinNishina,  ScalarHybridCompton, ScalarBetheHeitler,
                                   ScalarSauterGavrila, ScalarMollerBhabha,  ScalarSeltzerBerger};

// Geant4

typedef Real_t (*G4KernelFunc_t)(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos);

Real_t G4KleinNishina(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos);

Real_t G4HybridCompton(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos);

Real_t G4BetheHeitler(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos);

Real_t G4SauterGavrila(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos);

Real_t G4MollerBhabha(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos);

Real_t G4SeltzerBerger(int ntrack, GUTrack *itrack_aos, int *targetElements, GUTrack *otrack_aos);

G4KernelFunc_t G4KernelFunc[] = {G4KleinNishina,  G4HybridCompton, G4BetheHeitler,
                                 G4SauterGavrila, G4MollerBhabha,  G4SeltzerBerger};

// Vector

typedef Real_t (*VectorKernelFunc_t)(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa,
                                     SamplingMethod sampleType);

Real_t VectorKleinNishina(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa, SamplingMethod sampleType);

Real_t VectorHybridCompton(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa,
                           SamplingMethod sampleType);

Real_t VectorBetheHeitler(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa, SamplingMethod sampleType);

Real_t VectorSauterGavrila(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa,
                           SamplingMethod sampleType);

Real_t VectorMollerBhabha(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa, SamplingMethod sampleType);

Real_t VectorSeltzerBerger(GUTrack_v &itrack_soa, int *targetElements, GUTrack_v &otrack_soa,
                           SamplingMethod sampleType);

VectorKernelFunc_t VectorKernelFunc[] = {VectorKleinNishina,  VectorHybridCompton, VectorBetheHeitler,
                                         VectorSauterGavrila, VectorMollerBhabha,  VectorSeltzerBerger};

} // end namespace vecphys

#endif
