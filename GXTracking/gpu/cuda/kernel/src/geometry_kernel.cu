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

//Transportation
#include "GPFieldTrack.h"
#include "GXFieldMap.h"
#include "GXTrack.h"

#include "GXFieldTrack.h"
#include "GXEquationOfMotion.h"
#include "GXClassicalRK4.h"
#include "GXMagIntegratorDriver.h"
#include "GXChordFinder.h"

GLOBALFUNC
void navigator_kernel(GPGeomManager *geomManager,
		      GXFieldMap *magMap,
		      GXTrack *track, 
		      size_t nTrackSize) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  GPVPhysicalVolume *world = geomManager->fWorldPhysical;

  //Navigator
  GPNavigator aNavigator;
  GPNavigator_Constructor(&aNavigator);
  GPNavigator_SetWorldVolume(&aNavigator,world);

  GXTrack atrack;
  G4double x,y,z,px,py,pz,step;

  while(tid < nTrackSize) {

    //input
    atrack = track[tid];

    x    = atrack.x  ;
    y    = atrack.y  ;
    z    = atrack.z  ;
    px   = atrack.px ;
    py   = atrack.py ;
    pz   = atrack.pz ;
    step = atrack.s  ;

    GPThreeVector pos = GPThreeVector_create(x,y,z);
    GPThreeVector mom = GPThreeVector_create(px,py,pz);
    GPThreeVector dir = GPThreeVector_unit(mom);

    double safety = 0.0;
    GPThreeVector end = GPThreeVector_create(x + dir.x*step,
					     y + dir.y*step,
					     z + dir.z*step);

    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,pos,NULL,false,true);
    step = GPNavigator_ComputeStep(&aNavigator, pos, dir, step, &safety); 
    //output
    track[tid].s = step;

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for navigator kernel
//-----------------------------------------------------------------------------
void navigator_gpu(GPGeomManager *geomManager,
		   GXFieldMap *magMap,
		   GXTrack *track,  
		   size_t nTrackSize,
		   int NBLOCKS,
		   int NTHREADS) 
{
  int kstatus = 0;

  int threadsPerBlock = NTHREADS;
  int blocksPerGrid   = NBLOCKS; 

  navigator_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (geomManager,magMap,track,nTrackSize);

  kstatus = cudaThreadSynchronize();
  if(kstatus) std::cout << "navigator_gpu status = " << kstatus << "\n";
}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------
void navigator_cpu(GPGeomManager *geomManager,
		   GXFieldMap *magMap, 
		   GXTrack *track, 
		   size_t nTrackSize)
{
  GPVPhysicalVolume *world = geomManager->fWorldPhysical;

  //Navigator
  GPNavigator aNavigator;
  GPNavigator_Constructor(&aNavigator);
  GPNavigator_SetWorldVolume(&aNavigator,world);

  GXTrack atrack;
  G4double x,y,z,px,py,pz,step;

  for (int tid = 0 ; tid < nTrackSize ; tid++) {

    //input
    atrack = track[tid];

    x    = atrack.x  ;
    y    = atrack.y  ;
    z    = atrack.z  ;
    px   = atrack.px ;
    py   = atrack.py ;
    pz   = atrack.pz ;
    step = atrack.s  ;

    GPThreeVector pos = GPThreeVector_create(x,y,z);
    GPThreeVector mom = GPThreeVector_create(px,py,pz);
    GPThreeVector dir = GPThreeVector_unit(mom);

    double safety = 0.0;
    GPThreeVector end = GPThreeVector_create(x + dir.x*step,
					     y + dir.y*step,
					     z + dir.z*step);

    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,pos,NULL,false,true);
    step = GPNavigator_ComputeStep(&aNavigator, pos, dir,step, &safety); 

    //output
    track[tid].s = step;
  }
}

//--------------------------------------------
// MultilevelLocator
//--------------------------------------------
GLOBALFUNC
void mllocator_kernel(GPGeomManager *geomManager,
		      GXFieldMap *magMap,
		      GXTrack *track, 
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

  //Navigator
  GPNavigator aNavigator;
  GPNavigator_Constructor(&aNavigator);
  GPNavigator_SetWorldVolume(&aNavigator,world);

  //GPMultiLevelLocator
  GXMultiLevelLocator mLocator;

  const G4double e_mass = 0.510998910;

  GXFieldTrack aTrack;
  GXFieldTrack bTrack;
  GXFieldTrack IntersectPointVelct_G;

  GXTrack atrack;
  G4double x,y,z,px,py,pz,step,q,E;

  while(tid < nTrackSize) {

    //input
    atrack = track[tid];

    x    = atrack.x  ;
    y    = atrack.y  ;
    z    = atrack.z  ;
    px   = atrack.px ;
    py   = atrack.py ;
    pz   = atrack.pz ;
    step = atrack.s  ;
    q    = atrack.q;

    GPThreeVector pos = GPThreeVector_create(x,y,z);
    GPThreeVector mom = GPThreeVector_create(px,py,pz);
    GPThreeVector dir = GPThreeVector_unit(mom);
    G4double pinv = 1.0/GPThreeVector_mag(mom);

    //initialize the navigator and instantiate the multilevel locator
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,pos,NULL,false,true);

    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);

    GPVIntersectionLocator_SetChordFinderFor(&mLocator,
					     &chordFinder);

    //arguments
    GPThreeVector xend = GPThreeVector_create(x + dir.x*step,
					      y + dir.y*step,
					      z + dir.z*step);
    GPThreeVector pend = GPThreeVector_create(px + px*pinv,
					      py + py*pinv,
					      pz + py*pinv);

    E = sqrt(px*px+py*py+pz*pz+e_mass*e_mass);//calcuate 
    GXFieldTrack_Constructor(&aTrack,pos,mom,E,e_mass,q,100.0);
    GXFieldTrack_Constructor(&bTrack,xend,pend,E,e_mass,q,100.0+step);

    GPThreeVector intPos = GPThreeVector_create(0.5*(xend.x-pos.x),
						0.5*(xend.y-pos.y),
						0.5*(xend.z-pos.z));

    bool recalculatedEndPt = false;
    G4double previousSafety = 0.0;
    GPThreeVector previousSftOrigin = GPThreeVector_create(0,0,0);

    G4bool found_intersection =
      GXMultiLevelLocator_EstimateIntersectionPoint(&mLocator,
						    aTrack, bTrack, intPos,
						    IntersectPointVelct_G,
						    recalculatedEndPt,
						    previousSafety,
						    previousSftOrigin);
    //output
    pos = GXFieldTrack_GetPosition(&IntersectPointVelct_G);
    mom = GXFieldTrack_GetMomentum(&IntersectPointVelct_G);
    atrack.x = pos.x;
    atrack.y = pos.y;
    atrack.z = pos.z;
    atrack.px = mom.x;
    atrack.py = mom.y;
    atrack.pz = mom.z;

    track[tid] = atrack;

    tid += blockDim.x * gridDim.x;
  }
}

//-----------------------------------------------------------------------------
//  cuda wrapper for kernel
//-----------------------------------------------------------------------------
void mllocator_gpu( GPGeomManager *geomManager,
		  GXFieldMap *magMap,
		  GXTrack *track,  
		  size_t nTrackSize,
		  int NBLOCKS,
		  int NTHREADS) 
{
  int kstatus = 0;

  int threadsPerBlock = NTHREADS;
  int blocksPerGrid   = NBLOCKS; 

  mllocator_kernel<<< blocksPerGrid, threadsPerBlock >>>
    (geomManager,magMap,track,nTrackSize);

  kstatus = cudaThreadSynchronize();
  if(kstatus) std::cout << "geometry_gpu status = " << kstatus << "\n";

}

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void mllocator_cpu( GPGeomManager *geomManager,
		    GXFieldMap *magMap, 
		    GXTrack *track, 
		    size_t nTrackSize)
{
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

  //Navigator
  GPNavigator aNavigator;
  GPNavigator_Constructor(&aNavigator);
  GPNavigator_SetWorldVolume(&aNavigator,world);

  //GPMultiLevelLocator
  GXMultiLevelLocator mLocator;

  const G4double e_mass = 0.510998910;

  GXFieldTrack aTrack;
  GXFieldTrack bTrack;
  GXFieldTrack IntersectPointVelct_G;

  GXTrack atrack;
  G4double x,y,z,px,py,pz,step,q,E;

  for (int tid = 0 ; tid < nTrackSize ; tid++) {

    //input
    atrack = track[tid];

    x    = atrack.x  ;
    y    = atrack.y  ;
    z    = atrack.z  ;
    px   = atrack.px ;
    py   = atrack.py ;
    pz   = atrack.pz ;
    step = atrack.s  ;
    q    = atrack.q;

    GPThreeVector pos = GPThreeVector_create(x,y,z);
    GPThreeVector mom = GPThreeVector_create(px,py,pz);
    GPThreeVector dir = GPThreeVector_unit(mom);
    G4double pinv = 1.0/GPThreeVector_mag(mom);

    //initialize the navigator and instantiate the multilevel locator
    GPNavigator_LocateGlobalPointAndSetup(&aNavigator,pos,NULL,false,true);

    GXMultiLevelLocator_Constructor(&mLocator, &aNavigator);

    GPVIntersectionLocator_SetChordFinderFor(&mLocator,
					     &chordFinder);

    //arguments
    GPThreeVector xend = GPThreeVector_create(x + dir.x*step,
					      y + dir.y*step,
					      z + dir.z*step);
    GPThreeVector pend = GPThreeVector_create(px + px*pinv,
					      py + py*pinv,
					      pz + py*pinv);

    E = sqrt(px*px+py*py+pz*pz+e_mass*e_mass);//calcuate 
    GXFieldTrack_Constructor(&aTrack,pos,mom,E,e_mass,q,100.0);
    GXFieldTrack_Constructor(&bTrack,xend,pend,E,e_mass,q,100.0+step);

    GPThreeVector intPos = GPThreeVector_create(0.5*(xend.x-pos.x),
						0.5*(xend.y-pos.y),
						0.5*(xend.z-pos.z));

    bool recalculatedEndPt = false;
    G4double previousSafety = 0.0;
    GPThreeVector previousSftOrigin = GPThreeVector_create(0,0,0);

    G4bool found_intersection =
      GXMultiLevelLocator_EstimateIntersectionPoint(&mLocator,
						    aTrack, bTrack, intPos,
						    IntersectPointVelct_G,
						    recalculatedEndPt,
						    previousSafety,
						    previousSftOrigin);
    //output
    pos = GXFieldTrack_GetPosition(&IntersectPointVelct_G);
    mom = GXFieldTrack_GetMomentum(&IntersectPointVelct_G);
    atrack.x = pos.x;
    atrack.y = pos.y;
    atrack.z = pos.z;
    atrack.px = mom.x;
    atrack.py = mom.y;
    atrack.pz = mom.z;

    track[tid] = atrack;
  }
}

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
#include "GPElement.cu"

//-----------------------------------------------------------------------------
// GPGeometry
//-----------------------------------------------------------------------------
#include "GXMultiLevelLocator.cu"

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

#include "GPBox.cu"
#include "GPCons.cu"
#include "GPOrb.cu"
#include "GPTrd.cu"
#include "GPTubs.cu"
#include "GPVSolid.cu"

#include "GPFieldTrack.cu"
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
