#ifndef MODEL_KERNEL_H
#define MODEL_KERNEL_H 1

#include "GUTrack.h"
#include "base/VecPhys.h"

namespace vecphys {

typedef Real_t (*KernelFunc_t)(int ntrack, 
			       GUTrack* itrack_aos,
			       int *targetElements);

// Scalar

Real_t ScalarPhotonProcess(int ntrack, 
			   GUTrack* itrack_aos,
			   int *targetElements);

Real_t ScalarElectronProcess(int ntrack, 
		  	     GUTrack* itrack_aos,
			     int *targetElements);


KernelFunc_t ScalarKernelFunc[] = {ScalarPhotonProcess, 
                                   ScalarElectronProcess};

// Vector

typedef Real_t (*VectorKernelFunc_t)(GUTrack_v& itrack_soa,
     			             int *targetElements);

Real_t VectorPhotonProcess(GUTrack_v& itrack_soa,
     			   int *targetElements);

Real_t VectorElectronProcess(GUTrack_v& itrack_soa,
     		  	     int *targetElements);

VectorKernelFunc_t VectorKernelFunc[] = {VectorPhotonProcess,
                                         VectorElectronProcess}; 

// Geant3

Real_t Geant3PhotonProcess(int ntrack, 
			   GUTrack* itrack_aos,
			   int *targetElements);

Real_t Geant3ElectronProcess(int ntrack, 
		  	     GUTrack* itrack_aos,
			     int *targetElements);


KernelFunc_t Geant3KernelFunc[] = {Geant3PhotonProcess, 
                                   Geant3ElectronProcess};


// GeantV - hybrid

typedef Real_t (*GeantVKernelFunc_t)(GUTrack_v& itrack_soa,
     			             int *targetElements);

Real_t GeantVPhotonProcess(GUTrack_v& itrack_soa,
     			   int *targetElements);

Real_t GeantVElectronProcess(GUTrack_v& itrack_soa,
     		  	     int *targetElements);

GeantVKernelFunc_t GeantVKernelFunc[] = {GeantVPhotonProcess,
                                         GeantVElectronProcess}; 

} // end namespace vecphys

#endif
