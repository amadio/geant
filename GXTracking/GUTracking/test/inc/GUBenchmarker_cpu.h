#ifndef MODEL_KERNEL_H
#define MODEL_KERNEL_H 1

#include "GUTrack.h"
#include "base/Global.h"

namespace vecphys {

typedef Precision (*KernelFunc_t)(int ntrack, 
			          GUTrack* itrack_aos,
			          int *targetElements,
			          GUTrack* otrack_aos);

// Scalar

Precision ScalarKleinNishina(int ntrack, 
			     GUTrack* itrack_aos,
			     int *targetElements,
			     GUTrack* otrack_aos);

Precision ScalarBetheHeitler(int ntrack, 
			     GUTrack* itrack_aos,
			     int *targetElements,
			     GUTrack* otrack_aos);

Precision ScalarSauterGavrila(int ntrack, 
			      GUTrack* itrack_aos,
			      int *targetElements,
			      GUTrack* otrack_aos);

Precision ScalarMollerBhabha(int ntrack, 
			     GUTrack* itrack_aos,
			     int *targetElements,
			     GUTrack* otrack_aos);

Precision ScalarSeltzerBerger(int ntrack, 
			      GUTrack* itrack_aos,
			      int *targetElements,
			      GUTrack* otrack_aos);

KernelFunc_t ScalarKernelFunc[] = {ScalarKleinNishina, 
                                   ScalarBetheHeitler,
                                   ScalarSauterGavrila,
                                   ScalarMollerBhabha,
                                   ScalarSeltzerBerger}; 

//Geant4

Precision G4KleinNishina(int ntrack, 
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

KernelFunc_t G4KernelFunc[] = {G4KleinNishina, 
                               G4BetheHeitler,
                               G4SauterGavrila,
                               G4MollerBhabha,
                               G4SeltzerBerger}; 

// Vector

typedef Precision (*VectorKernelFunc_t)(GUTrack_v& itrack_soa,
     			                int *targetElements,
			                GUTrack_v& otrack_soa);

Precision VectorKleinNishina(GUTrack_v& itrack_soa,
     			     int *targetElements,
			     GUTrack_v& otrack_soa);

Precision VectorBetheHeitler(GUTrack_v& itrack_soa,
     			     int *targetElements,
			     GUTrack_v& otrack_soa);

Precision VectorSauterGavrila(GUTrack_v& itrack_soa,
     			      int *targetElements,
			      GUTrack_v& otrack_soa);

Precision VectorMollerBhabha(GUTrack_v& itrack_soa,
     			     int *targetElements,
			     GUTrack_v& otrack_soa);

Precision VectorSeltzerBerger(GUTrack_v& itrack_soa,
			      int *targetElements,
			      GUTrack_v& otrack_soa);

VectorKernelFunc_t VectorKernelFunc[] = {VectorKleinNishina,
                                         VectorBetheHeitler,
                                         VectorSauterGavrila,
                                         VectorMollerBhabha, 
                                         VectorSeltzerBerger}; 

} // end namespace vecphys

#endif
