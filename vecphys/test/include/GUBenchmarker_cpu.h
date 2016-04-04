#ifndef MODEL_KERNEL_H
#define MODEL_KERNEL_H 1

#include "GUTrack.h"
#include "base/Global.h"
#include "SamplingMethod.h"

namespace vecphys {

typedef Precision (*KernelFunc_t)(int ntrack,
			          GUTrack* itrack_aos,
			          int *targetElements,
			          GUTrack* otrack_aos,
                                  SamplingMethod sampleType);

// Scalar

Precision ScalarKleinNishina(int ntrack,
			     GUTrack* itrack_aos,
			     int *targetElements,
			     GUTrack* otrack_aos,
                             SamplingMethod sampleType);

Precision ScalarHybridCompton(int ntrack,
			      GUTrack* itrack_aos,
			      int *targetElements,
			      GUTrack* otrack_aos,
                              SamplingMethod sampleType);

Precision ScalarBetheHeitler(int ntrack,
			     GUTrack* itrack_aos,
			     int *targetElements,
			     GUTrack* otrack_aos,
                             SamplingMethod sampleType);

Precision ScalarSauterGavrila(int ntrack,
			      GUTrack* itrack_aos,
			      int *targetElements,
			      GUTrack* otrack_aos,
                              SamplingMethod sampleType);

Precision ScalarMollerBhabha(int ntrack,
			     GUTrack* itrack_aos,
			     int *targetElements,
			     GUTrack* otrack_aos,
                             SamplingMethod sampleType);

Precision ScalarSeltzerBerger(int ntrack,
			      GUTrack* itrack_aos,
			      int *targetElements,
			      GUTrack* otrack_aos,
                              SamplingMethod sampleType);

KernelFunc_t ScalarKernelFunc[] = {ScalarKleinNishina,
                                   ScalarHybridCompton,
                                   ScalarBetheHeitler,
                                   ScalarSauterGavrila,
                                   ScalarMollerBhabha,
                                   ScalarSeltzerBerger};

//Geant4

typedef Precision (*G4KernelFunc_t)(int ntrack,
		  	            GUTrack* itrack_aos,
			            int *targetElements,
			            GUTrack* otrack_aos);

Precision G4KleinNishina(int ntrack,
			 GUTrack* itrack_aos,
			 int *targetElements,
			 GUTrack* otrack_aos);

Precision G4HybridCompton(int ntrack,
			  GUTrack* itrack_aos,
			  int *targetElements,
			  GUTrack* otrack_aos);

Precision G4BetheHeitler(int ntrack,
			 GUTrack* itrack_aos,
			 int *targetElements,
			 GUTrack* otrack_aos);

Precision G4SauterGavrila(int ntrack,
			  GUTrack* itrack_aos,
			  int *targetElements,
			  GUTrack* otrack_aos);

Precision G4MollerBhabha(int ntrack,
			 GUTrack* itrack_aos,
			 int *targetElements,
			 GUTrack* otrack_aos);

Precision G4SeltzerBerger(int ntrack,
			  GUTrack* itrack_aos,
			  int *targetElements,
			  GUTrack* otrack_aos);

G4KernelFunc_t G4KernelFunc[] = {G4KleinNishina,
                                 G4HybridCompton,
                                 G4BetheHeitler,
                                 G4SauterGavrila,
                                 G4MollerBhabha,
                                 G4SeltzerBerger};

// Vector

typedef Precision (*VectorKernelFunc_t)(GUTrack_v& itrack_soa,
     			                int *targetElements,
			                GUTrack_v& otrack_soa,
                                        SamplingMethod sampleType);

Precision VectorKleinNishina(GUTrack_v& itrack_soa,
     			     int *targetElements,
			     GUTrack_v& otrack_soa,
                             SamplingMethod sampleType);

Precision VectorHybridCompton(GUTrack_v& itrack_soa,
     			      int *targetElements,
			      GUTrack_v& otrack_soa,
                              SamplingMethod sampleType);

Precision VectorBetheHeitler(GUTrack_v& itrack_soa,
     			     int *targetElements,
			     GUTrack_v& otrack_soa,
                             SamplingMethod sampleType);

Precision VectorSauterGavrila(GUTrack_v& itrack_soa,
     			      int *targetElements,
			      GUTrack_v& otrack_soa,
                              SamplingMethod sampleType);

Precision VectorMollerBhabha(GUTrack_v& itrack_soa,
     			     int *targetElements,
			     GUTrack_v& otrack_soa,
                             SamplingMethod sampleType);

Precision VectorSeltzerBerger(GUTrack_v& itrack_soa,
			      int *targetElements,
			      GUTrack_v& otrack_soa,
                              SamplingMethod sampleType);

VectorKernelFunc_t VectorKernelFunc[] = {VectorKleinNishina,
                                         VectorHybridCompton,
                                         VectorBetheHeitler,
                                         VectorSauterGavrila,
                                         VectorMollerBhabha,
                                         VectorSeltzerBerger};

} // end namespace vecphys

#endif
