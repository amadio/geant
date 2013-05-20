#include <iostream>
#include <memory>

//Common

#include "GPThreeVector.h"
#include "GPThreeVectorList.h"
#include "GPRotationMatrix.h"
#include "GPUtils.h"

//Material
#include "GPMaterial.h"

//Geometry
#include "GPGeomManager.h"
#include "GPVPhysicalVolume.h"
#include "GPTouchableHistory.h"
#include "GPNavigator.h"
#include "GXMultiLevelLocator.h"

// for debugging
#include "GPLogicalVolume.h"
#include "GPTubs.h"

//Transportation
#include "GPFieldMap.h"
#include "GXTrack.h"

#include "GXFieldTrack.h"
#include "GXEquationOfMotion.h"
#include "GXClassicalRK4.h"
#include "GXMagIntegratorDriver.h"
#include "GXChordFinder.h"
#include "GXFieldManager.h"
#include "GXPropagatorInField.h"
#include "GXTransportation.h"

#include "GXBrem.h"
#include "GXIoni.h"
#include "GXMsc.h"

#include "GPSTLVector.h"

#include "stdio.h"

#include "util_kernel.h"

//-----------------------------------------------------------------------------
//  Transportation engine
//-----------------------------------------------------------------------------
GLOBALFUNC
void tracking_kernel(curandState* devStates,
                     GPGeomManager *geomManager,
                     GPFieldMap *magMap,
                     GXTrack *track,
                     int *logVolumeIndices,
                     int *physVolumeIndices,
                     GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
                     size_t nTrackSize)
{
   unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
   
   GPVPhysicalVolume *world = geomManager->fWorldPhysical;
   
   GXMagneticField magField;
   GXMagneticField_Constructor(&magField,magMap);
   
   GXEquationOfMotion equaOfMotion;
   GXEquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);
   
   GXClassicalRK4 rk4;
   GXClassicalRK4_Constructor(&rk4,&equaOfMotion);
   
   GXMagInt_Driver magDriver;
   GXMagInt_Driver_Constructor(&magDriver,1.0,&rk4);
   
   GXChordFinder chordFinder;
   GXChordFinder_Constructor(&chordFinder,&magDriver);
   
   GXFieldManager aFieldManager;
   GXFieldManager_Constructor(&aFieldManager,&magField,&chordFinder, true);
   
   //Navigator
   GPNavigator aNavigator;
   GPNavigator_Constructor(&aNavigator);
   GPNavigator_SetWorldVolume(&aNavigator,world);
   
   //GPMultiLevelLocator
   GXMultiLevelLocator mLocator;
   
   //Propagator
   GXPropagatorInField propagatorInField;
   
   //Transporation
   GXTransportation transp;
   
   //EMPhysics
   GXBrem brem;
   GXIoni ioni;
   GXMsc msc;
   GPForceCondition condition;
   GPGPILSelection  selection;
   
   EMPhysics_Init(tid,devStates,&brem,&ioni,&msc,
                  eBrem_table,eIoni_table,msc_table);
   
   G4double x,y,z,step;
   G4double proposedStep = 0.0;
   
   GXTrack atrack;
   
   while(tid < nTrackSize) {
      if (track[tid].q != -1) {
         tid += blockDim.x * gridDim.x;
         continue;
      }
      
      atrack = track[tid];
      
      x    = atrack.x  ;
      y    = atrack.y  ;
      z    = atrack.z  ;
      step = atrack.s   ;
      
      //    printf("Trying to navigate %d/%ld with step size %f with %f %f\n",tid,nTrackSize, step,atrack.px, atrack.q);
      // printf("Trying %d: %f %f %f\n",tid,x,y,z);
      
      //Geometry - Initialize Navigator and construct related structures

      // We know where we are coming from, so we need to upgrade the
      // following code to use the location information we already have.
      GPVPhysicalVolume *where = GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(x,y,z),NULL,false,true);
      if (where) {
         if (where->fIndex != physVolumeIndices[tid])
            printf("ERROR: %d Found the track to be in phys volume node #%d expected %d at %lf %lf %lf\n",
                   tid,where->fIndex,physVolumeIndices[tid],
                   x,y,z);
         if (where->flogical->fIndex != logVolumeIndices[tid])
            printf("ERROR: %d Found the track to be in log volume node #%d expected %d at %lf %lf %lf\n",
                   tid,where->flogical->fIndex,logVolumeIndices[tid],
                   x,y,z);
      } else printf("ERROR: Could not find where the track is ...\n");
      //     if (where && where->flogical) {
      //        printf("Logical volume's solid: %p\n",where->flogical->fSolid);
      //        printf("Logical volume's solid: SPhi=%f\n",((GPTubs*)where->flogical->fSolid)->fSPhi);
      //     }
      // GPVPhysicalVolume *where = geomManager->fVolumesIndex[volumeIndices[tid]];
      
      GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
      
      GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
                                      &aFieldManager,&mLocator);
      
      GXTransportation_Constructor(&transp,&propagatorInField,0);
      
      //EM Physics
      proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,&atrack,&condition,&selection);
      
      //Transportation
      proposedStep += step;
      GPGPILSelection  fGPILSelection;
      G4double safety = 0.0;
      
      //add energyLoss here as faking output
      track[tid].s = GXTransportation_AlongStepGetPhysicalInteractionLength(&transp,
                                                                            &atrack,
                                                                            1.0,
                                                                            proposedStep,
                                                                            safety,
                                                                            &fGPILSelection);
      double norm = sqrt(atrack.px*atrack.px+atrack.py*atrack.py+atrack.pz*atrack.pz);
      atrack.x += atrack.px * track[tid].s / norm;
      atrack.y += atrack.py * track[tid].s / norm;
      atrack.z += atrack.pz * track[tid].s / norm;
      track[tid] = atrack;
      //    printf("Got %d: %f %f %f exit: %d vol: %d %ld\n",tid,
      //           atrack.x,atrack.y,atrack.z,
      //           aNavigator.fExiting,
      //           aNavigator.fBlockedPhysicalVolume ? aNavigator.fBlockedPhysicalVolume->fIndex : -1,
      //           aNavigator.fBlockedPhysicalVolume);
      // if (aNavigator.fBlockedPhysicalVolume) {
      //    physVolumeIndices[tid] = aNavigator.fBlockedPhysicalVolume->fIndex;
      //    logVolumeIndices[tid] = aNavigator.fBlockedPhysicalVolume->flogical->fIndex;
      //    // printf("Moved from %d to in %d\n",logVolumeIndices[tid],aNavigator.fBlockedPhysicalVolume ? aNavigator.fBlockedPhysicalVolume->flogical->fIndex : -1);
      // } else {
         // humm so we have not moved ... let's check..
         where = GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(atrack.x,atrack.y,atrack.z),NULL,false,true);
         if (!where) {
            // printf("Moved out of the geometry\n");
            physVolumeIndices[tid] = -1;
            logVolumeIndices[tid] = -1;
         } else {
            // // if (where->flogical->fIndex != logVolumeIndices[tid] || where->fIndex != physVolumeIndices[tid]) {
            //    printf("Moved from %d:%d to %d:%d\n", logVolumeIndices[tid], physVolumeIndices[tid],
            //           where->flogical->fIndex, where->fIndex);
            // }
            physVolumeIndices[tid] = where->fIndex;
            logVolumeIndices[tid] = where->flogical->fIndex;
         }            
      // }
      tid += blockDim.x * gridDim.x;
      
   }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void tracking_gpu(curandState* devStates,
                  GPGeomManager *geomManager,
                  GPFieldMap *magMap,
                  GXTrack *track,
                  int *logVolumeIndices,
                  int *physVolumeIndices,
                  GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
                  size_t nTrackSize,
                  int NBLOCKS,
                  int NTHREADS,
                  cudaStream_t stream)
{
   int threadsPerBlock = NTHREADS;
   int blocksPerGrid   = NBLOCKS;
   //  int blocksPerGrid   = (int)((nTrackSize+threadsPerBlock-1)/threadsPerBlock);
   
   tracking_kernel<<< blocksPerGrid, threadsPerBlock, 0 , stream >>>(devStates,
                                                                     geomManager,magMap,
                                                                     track,logVolumeIndices,physVolumeIndices,
                                                                     eBrem_table,eIoni_table,msc_table,
                                                                     nTrackSize);
}

//-----------------------------------------------------------------------------
//  Transportation engine - Electron/Positron
//-----------------------------------------------------------------------------

GLOBALFUNC
void tracking_electron_kernel(curandState* devStates,
                              GPVPhysicalVolume *world,
                              GPFieldMap *magMap,
                              GXTrack *track,
                              GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
                              size_t nTrackSize)
{
   unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
   
   GXMagneticField magField;
   GXMagneticField_Constructor(&magField,magMap);
   
   GXEquationOfMotion equaOfMotion;
   GXEquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);
   
   GXClassicalRK4 rk4;
   GXClassicalRK4_Constructor(&rk4,&equaOfMotion);
   
   GXMagInt_Driver magDriver;
   GXMagInt_Driver_Constructor(&magDriver,1.0,&rk4);
   
   GXChordFinder chordFinder;
   GXChordFinder_Constructor(&chordFinder,&magDriver);
   
   GXFieldManager aFieldManager;
   GXFieldManager_Constructor(&aFieldManager,&magField,&chordFinder, true);
   
   //Navigator
   GPNavigator aNavigator;
   GPNavigator_Constructor(&aNavigator);
   GPNavigator_SetWorldVolume(&aNavigator,world);
   
   //GPMultiLevelLocator
   GXMultiLevelLocator mLocator;
   
   //Propagator
   GXPropagatorInField propagatorInField;
   
   //Transporation
   GXTransportation transp;
   
   //EMPhysics
   GXBrem brem;
   GXIoni ioni;
   GXMsc msc;
   GPForceCondition condition;
   GPGPILSelection  selection;
   
   EMPhysics_Init(tid,devStates,&brem,&ioni,&msc,
                  eBrem_table,eIoni_table,msc_table);
   
   G4double x,y,z,step;
   G4double proposedStep = 0.0;
   
   GXTrack atrack;
   
   while(tid < nTrackSize) {
      
      atrack = track[tid];
      
      x    = atrack.x  ;
      y    = atrack.y  ;
      z    = atrack.z  ;
      step = atrack.s   ;
      
      //Geometry - Initialize Navigator and construct related structures
      GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(x,y,z),
                                            NULL,false,true);
      
      GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
      
      GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
                                      &aFieldManager,&mLocator);
      
      GXTransportation_Constructor(&transp,&propagatorInField,0);
      
      //EM Physics
      proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,&atrack,
                                                &condition,&selection);
      
      //Transportation
      proposedStep += step;
      GPGPILSelection  fGPILSelection;
      G4double safety = 0.0;
      
      //add energyLoss here as faking output
      track[tid].s = GXTransportation_AlongStepGPIL_Electron(&transp,
                                                             &atrack,
                                                             1.0,
                                                             proposedStep,
                                                             safety,
                                                             &fGPILSelection);
      tid += blockDim.x * gridDim.x;
      
   }
   
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void tracking_electron_gpu(curandState* devStates,
                           GPVPhysicalVolume *world,
                           GPFieldMap *magMap,
                           GXTrack *track,
                           GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
                           size_t nTrackSize,
                           int NBLOCKS,
                           int NTHREADS,
                           cudaStream_t stream)
{
   int threadsPerBlock = NTHREADS;
   int blocksPerGrid   = NBLOCKS;
   
   tracking_electron_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
   (devStates,world,magMap,track,eBrem_table,eIoni_table,msc_table,nTrackSize);
   
}


//-----------------------------------------------------------------------------
//  Transportation engine - Photon
//-----------------------------------------------------------------------------

GLOBALFUNC
void tracking_photon_kernel(curandState* devStates,
                            GPVPhysicalVolume *world,
                            GPFieldMap *magMap,
                            GXTrack *track,
                            GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
                            size_t nTrackSize)
{
   unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
   
   GXMagneticField magField;
   GXMagneticField_Constructor(&magField,magMap);
   
   GXEquationOfMotion equaOfMotion;
   GXEquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);
   
   GXClassicalRK4 rk4;
   GXClassicalRK4_Constructor(&rk4,&equaOfMotion);
   
   GXMagInt_Driver magDriver;
   GXMagInt_Driver_Constructor(&magDriver,1.0,&rk4);
   
   GXChordFinder chordFinder;
   GXChordFinder_Constructor(&chordFinder,&magDriver);
   
   GXFieldManager aFieldManager;
   GXFieldManager_Constructor(&aFieldManager,&magField,&chordFinder, true);
   
   //Navigator
   GPNavigator aNavigator;
   GPNavigator_Constructor(&aNavigator);
   GPNavigator_SetWorldVolume(&aNavigator,world);
   
   //GPMultiLevelLocator
   GXMultiLevelLocator mLocator;
   
   //Propagator
   GXPropagatorInField propagatorInField;
   
   //Transporation
   GXTransportation transp;
   
   //EMPhysics
   GXBrem brem;
   GXIoni ioni;
   GXMsc msc;
   GPForceCondition condition;
   GPGPILSelection  selection;
   
   EMPhysics_Init(tid,devStates,&brem,&ioni,&msc,
                  eBrem_table,eIoni_table,msc_table);
   
   G4double x,y,z,step;
   G4double proposedStep = 0.0;
   
   GXTrack atrack;
   
   while(tid < nTrackSize) {
      
      atrack = track[tid];
      
      x    = atrack.x  ;
      y    = atrack.y  ;
      z    = atrack.z  ;
      step = atrack.s   ;
      
      //Geometry - Initialize Navigator and construct related structures
      GPNavigator_LocateGlobalPointAndSetup(&aNavigator,GPThreeVector_create(x,y,z),
                                            NULL,false,true);
      
      GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);
      
      GXPropagatorInField_Constructor(&propagatorInField,&aNavigator,
                                      &aFieldManager,&mLocator);
      
      GXTransportation_Constructor(&transp,&propagatorInField,0);
      
      //EM Physics
      proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,&atrack,
                                                &condition,&selection);
      
      //Transportation
      proposedStep += step;
      GPGPILSelection  fGPILSelection;
      G4double safety = 0.0;
      
      //add energyLoss here as faking output
      track[tid].s = GXTransportation_AlongStepGPIL_Photon(&transp,
                                                           &atrack,
                                                           1.0,
                                                           proposedStep,
                                                           safety,
                                                           &fGPILSelection);
      tid += blockDim.x * gridDim.x;
      
   }
   
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void tracking_photon_gpu(curandState* devStates,
                         GPVPhysicalVolume *world,
                         GPFieldMap *magMap,
                         GXTrack *track,
                         GPPhysicsTable* eBrem_table, GPPhysicsTable* eIoni_table, GPPhysicsTable* msc_table,
                         size_t nTrackSize,
                         int NBLOCKS,
                         int NTHREADS,
                         cudaStream_t stream)
{
   int threadsPerBlock = NTHREADS;
   int blocksPerGrid   = NBLOCKS; 
   
   tracking_photon_kernel<<< blocksPerGrid, threadsPerBlock, 0, stream >>>
   (devStates,world,magMap,track,eBrem_table,eIoni_table,msc_table,nTrackSize);
}

#if 0
//-----------------------------------------------------------------------------
// Other Kernels
//-----------------------------------------------------------------------------
#include "random_kernel.cu"
#include "util_kernel.cu"

//-----------------------------------------------------------------------------
// GPCommon
//-----------------------------------------------------------------------------

#include "GPThreeVector.cu"
#include "GPThreeVectorList.cu"
#include "GPRotationMatrix.cu"
#include "GPUtils.cu"

//-----------------------------------------------------------------------------
// GPMaterial
//-----------------------------------------------------------------------------

#include "GPMaterial.cu"

//-----------------------------------------------------------------------------
// GPGeometry
//-----------------------------------------------------------------------------
#include "GPAffineTransform.cu"
#include "GPAuxiliaryNavServices.cu"
#include "GPCombinedNavigation.cu"
#include "GPLineSection.cu"
#include "GPLogicalVolume.cu"
#include "GPNavigationLevel.cu"
#include "GPNavigator.cu"
#include "GPNavigationHistory.cu"
#include "GPNormalNavigation.cu"
#include "GPSmartVoxelHeader.cu"
#include "GPSmartVoxelNode.cu"
#include "GPSmartVoxelProxy.cu"
#include "GPTouchableHistory.cu"
#include "GPVoxelLimits.cu"
#include "GPVoxelNavigation.cu"
#include "GPVPhysicalVolume.cu"

#include "GXMultiLevelLocator.cu"

#include "GPBox.cu"
#include "GPCons.cu"
#include "GPOrb.cu"
#include "GPTrd.cu"
#include "GPTubs.cu"
#include "GPVSolid.cu"

#include "GPChargeState.cu"

#include "GPUserGeometry.cu"
#include "GPVoxelHeader.cu"

//chordFiner
#include "GXFieldTrack.cu"
#include "GXClassicalRK4.cu"
#include "GXEquationOfMotion.cu"
#include "GXMagneticField.cu"
#include "GXMagIntegratorDriver.cu"
#include "GXChordFinder.cu"
#include "GXFieldManager.cu"
#include "GXPropagatorInField.cu"
#include "GXTransportation.cu"

//EMPhysics
#include "GPPhysicsTable.cu"
#include "GXBrem.cu"
#include "GXIoni.cu"
#include "GXMsc.cu"

#endif
