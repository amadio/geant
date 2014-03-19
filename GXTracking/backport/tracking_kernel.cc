#include <iostream>
#include <memory>

//Common
#include "G4ThreeVector.hh"
// #include "GXThreeVectorList.h"
#include "G4RotationMatrix.hh"
// #include "GPUtils.hh"

//Material
#include "G4Material.hh"

//Geometry
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Navigator.hh"
#include "G4MultiLevelLocator.hh"

//Transportation
#include "GPFieldMap.h"
#include "G4Track.hh"

#include "G4FieldTrack.hh"
#include "G4EquationOfMotion.hh"
#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"
#include "G4PropagatorInField.hh"
#include "G4Transportation.hh"

// #include "GXBrem.h"
// #include "GXIoni.h"
// #include "GXMsc.h"

//#include "GPSTLVector.hh"
#include <vector>

#include <stdio.h>

// #include "util_kernel.h"


#include "G4MagneticField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4Electron.hh"
#include "G4ProductionCutsTable.hh"
#include "G4UrbanMscModel95.hh"

int EMPhysics_Init(unsigned int /* tid */,
		    void /* curandState*/ * /* devStates */,
		    G4VEnergyLossProcess /* GXBrem */ *brem,
		    G4VEnergyLossProcess /* GXIoni */ *ioni,
                    G4VMultipleScattering /* GXMsc */  *msc,
		    G4PhysicsTable* eBrem_table, 
		    G4PhysicsTable* eIoni_table, 
		    G4PhysicsTable*  msc_table) 
{
   // G4ProductionCutsTable::GetProductionCutsTable()->RetrieveCutsTable("./data",true); 

  //Initialize EMPhysics Processes
  // brem->SetCurandState(&devStates[tid]);
  brem->SetIntegral(true);
  brem->PreparePhysicsTable(*G4Electron::Definition());
  brem->SetLambdaTable(eBrem_table);
  //brem->UseLambdaTable(true);

  // ioni->SetCurandState(&devStates[tid]);
  ioni->SetIntegral(true);
  ioni->PreparePhysicsTable(*G4Electron::Definition());
  ioni->SetLambdaTable(eIoni_table);
  //ioni->UseLambdaTable(true);

  // msc->SetCurandState(&devStates[tid]);
  //msc->SetIntegral(true);
  msc->AddEmModel(1, new G4UrbanMscModel95());
  msc->PreparePhysicsTable(*G4Electron::Definition());
   msc->RetrievePhysicsTable(G4Electron::Definition(),"data",true);
  //msc->SetLambdaTable(msc_table);
  //msc->UseLambdaTable(true);
  //msc->SetBuildLambdaTable(true);
   return 0;
}

G4double EMPhysics_DefineStepLength(G4eBremsstrahlung     *brem,
                                    G4eIonisation         *ioni,
                                    G4eMultipleScattering *msc,
                                    G4Track *atrack,
                                    G4ForceCondition *condition,
                                    G4GPILSelection  *selection,
                                    G4double &safety) 
{
   //G4double px,py,pz;
   G4double step;
   // G4double E;
   
   G4double energyLoss = 0.0;
   G4double step_brem = 0.0;
   G4double step_ioni = 0.0;
   G4double step_msc  = 0.0;
   G4double proposedStep = 0.0;
   //EM Physics
   G4ThreeVector momentum(atrack->GetMomentum());
   // px   = momentum.x();
   // py   = momentum.y(); //atrack->py;
   // pz   = momentum.z(); // atrack->pz;
   step = atrack->GetStepLength();

   // E = sqrt(px*px+py*py+pz*pz+0.510998910*0.510998910);//calcuate 

   step_brem = brem->PostStepGetPhysicalInteractionLength(*atrack, step, condition);
   step_ioni = ioni->PostStepGetPhysicalInteractionLength(*atrack, step, condition);
   step_msc = msc->PostStepGetPhysicalInteractionLength(*atrack, step, condition);

   //physics model defining the current step
   unsigned int model =0;
   if(step_brem > step_ioni) model = 1;   
   if(step_ioni > step_msc) model = 2;  
   
   G4Step stepinfo;
   G4VParticleChange *pchange = 0;
   switch (model) {
   case 0 : 
      proposedStep = step_brem;
      pchange = brem->PostStepDoIt(*atrack,stepinfo);
      energyLoss += pchange->GetLocalEnergyDeposit();
      break; 
   case 1 :
      proposedStep = step_ioni;
      pchange = ioni->PostStepDoIt(*atrack,stepinfo);
      energyLoss += pchange->GetLocalEnergyDeposit();
      break; 
   case 2 :
      proposedStep = step_msc;
      pchange = msc->PostStepDoIt(*atrack,stepinfo);
      energyLoss += pchange->GetLocalEnergyDeposit();
      break; 
   }
   pchange->Clear();

   //alongstep
   // its does not seems right to pass step twice ... the cude kernal pass only currentMinimalStep
   step_msc = msc->AlongStepGetPhysicalInteractionLength(*atrack, step, step, safety, selection);
   pchange = msc->AlongStepDoIt(*atrack,stepinfo);
   energyLoss += pchange->GetLocalEnergyDeposit();
   if(step_msc < proposedStep) proposedStep = step_msc;
   
   return proposedStep;
}

//-----------------------------------------------------------------------------
//  Transportation engine
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void tracking_cpu(G4VPhysicalVolume *world, G4Navigator &aNavigator,
         G4MagneticField *magMap, 
         G4Track **track, 
         G4PhysicsTable* eBrem_table, G4PhysicsTable* eIoni_table, G4PhysicsTable* msc_table,
         size_t nTrackSize)
{
   // cmsExpG4MagneticField magField(magMap);

   G4Mag_UsualEqRhs equaOfMotion(magMap);
   // G4EquationOfMotion equaOfMotion(magMap);
   // G4EquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);
   
   G4ClassicalRK4 rk4( &equaOfMotion );
   //G4ClassicalRK4_Constructor(&rk4,&equaOfMotion);

   G4MagInt_Driver *magDriver = new G4MagInt_Driver(1.0, &rk4);
   //G4MagInt_Driver_Constructor(&magDriver,1.0,&rk4);

   G4ChordFinder chordFinder(magDriver);
   //G4ChordFinder_Constructor(&chordFinder,&magDriver);
   
   G4FieldManager aFieldManager(magMap,&chordFinder, true);
   //G4FieldManager_Constructor(&aFieldManager,&magField,&chordFinder, true);
   
   // //Navigator
   // G4Navigator aNavigator;
   // aNavigator.SetWorldVolume(world);
   // //G4Navigator_Constructor(&aNavigator);
   // //G4Navigator_SetWorldVolume(&aNavigator,world);
   
   //GPMultiLevelLocator
   G4MultiLevelLocator mLocator(&aNavigator);
   
   //Propagator
   G4PropagatorInField propagatorInField(&aNavigator,&aFieldManager,&mLocator);
   
   //Transporation
   G4Transportation transp(0);
   transp.SetPropagatorInField(&propagatorInField);
   
   //EMPhysics
   static G4eBremsstrahlung brem; // GXBrem brem;
   static G4eIonisation ioni; // GXIoni ioni;
   static G4eMultipleScattering msc; // GXMsc msc;
   
   G4ProductionCutsTable::GetProductionCutsTable()->UpdateCoupleTable(world);
   static int s = EMPhysics_Init(0,NULL,&brem,&ioni,&msc,
                                 eBrem_table,eIoni_table,msc_table);

   
   G4ForceCondition condition;
   G4GPILSelection  selection;
   

   // G4double x,y,z;
   G4double step = 0.0;
   G4double proposedStep = 0.0;
   
   for (size_t tid = 0 ; tid < nTrackSize ; tid++) {
      
      //Geometry - Initialize Navigator and construct related structures
      aNavigator.LocateGlobalPointAndSetup(track[tid]->GetPosition(),
                                           NULL,false,true);
      
      // G4MultiLevelLocator_Constructor(&mLocator, &aNavigator);
      
      // G4PropagatorInField_Constructor(&propagatorInField,&aNavigator,
      //                                 &aFieldManager,&mLocator); 
      
      // G4Transportation_Constructor(&transp,&propagatorInField,0);
      
      //EM Physics
      G4double safety = 0.0;
      proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,track[tid],
                                                &condition,&selection,safety);
      
      //Transportation
      proposedStep += step;
      G4GPILSelection  fGPILSelection;
      
      //add energyLoss here as faking output
      track[tid]->SetStepLength(transp.AlongStepGetPhysicalInteractionLength(*track[tid],
                                                                             1.0,
                                                                             proposedStep,
                                                                             safety,
                                                                             &fGPILSelection));
      // track[tid].s = G4Transportation_AlongStepGetPhysicalInteractionLength(&transp,
      //    								  &atrack,
      //    								  1.0,
      //    								  proposedStep,
      //    								  safety,
      //    								  &fGPILSelection);
   }
  
}

//-----------------------------------------------------------------------------
//  Transportation engine - Electron/Positron
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void tracking_electron_cpu( G4VPhysicalVolume *world, G4Navigator &aNavigator,
                            G4MagneticField *magMap, 
			    G4Track **track, 
			    G4PhysicsTable* eBrem_table, G4PhysicsTable* eIoni_table, G4PhysicsTable* msc_table,
			    size_t nTrackSize)
{
  // G4MagneticField magField;
  // G4MagneticField_Constructor(&magField,magMap);

   G4Mag_UsualEqRhs equaOfMotion(magMap);
   //  G4EquationOfMotion equaOfMotion;
   //G4EquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);

   G4ClassicalRK4 rk4(&equaOfMotion);
   //G4ClassicalRK4_Constructor(&rk4,&equaOfMotion);

   G4MagInt_Driver *magDriver = new G4MagInt_Driver(1.0, &rk4);
   //G4MagInt_Driver_Constructor(&magDriver,1.0,&rk4);
   
   G4ChordFinder chordFinder(magDriver);
   //G4ChordFinder_Constructor(&chordFinder,&magDriver);
   
   G4FieldManager aFieldManager(magMap,&chordFinder, true);
   //G4FieldManager_Constructor(&aFieldManager,&magField,&chordFinder, true);

   // //Navigator
   // G4Navigator aNavigator;
   // aNavigator.SetWorldVolume(world);
   // //G4Navigator_Constructor(&aNavigator);
   // //G4Navigator_SetWorldVolume(&aNavigator,world);
   
   //G4MultiLevelLocator
   G4MultiLevelLocator mLocator(&aNavigator);
   
   //Propagator
   G4PropagatorInField propagatorInField(&aNavigator,&aFieldManager,&mLocator);

  //Transporation
  G4Transportation transp(0);
  transp.SetPropagatorInField(&propagatorInField);

  //EMPhysics
  static G4eBremsstrahlung brem; // GXBrem brem;
  static G4eIonisation ioni; // GXIoni ioni;
  static G4eMultipleScattering msc; // GXMsc msc;
  G4ForceCondition condition;
  G4GPILSelection  selection;

   G4ProductionCutsTable::GetProductionCutsTable()->UpdateCoupleTable(world);
  static int init = EMPhysics_Init(0,NULL,&brem,&ioni,&msc,
                                   eBrem_table,eIoni_table,msc_table);

  // G4double x,y,z;
  G4double step = 0.0;
  G4double proposedStep = 0.0;

  G4Track atrack;

  for (size_t tid = 0 ; tid < nTrackSize ; tid++) {

    //Geometry - Initialize Navigator and construct related structures
     aNavigator.LocateGlobalPointAndSetup(track[tid]->GetPosition(),
                                          NULL,false,true);

     // G4MultiLevelLocator_Constructor(&mLocator, &aNavigator);
     
     // G4PropagatorInField_Constructor(&propagatorInField,&aNavigator,
     //    			    &aFieldManager,&mLocator); 
    
     // G4Transportation_Constructor(&transp,&propagatorInField,0);
    
     //EM Physics
    G4double safety = 0.0;

    proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,track[tid],&condition,&selection,safety);
    
    //Transportation
    proposedStep += step;
    G4GPILSelection  fGPILSelection;

    //add energyLoss here as faking output
    track[tid]->SetStepLength(transp.AlongStepGetPhysicalInteractionLength(*track[tid],
                                                                           1.0,
                                                                           proposedStep,
                                                                           safety,
                                                                           &fGPILSelection));
    // track[tid].s = G4Transportation_AlongStepGPIL_Electron(&transp,
    //     						   &atrack,
    //     						   1.0,
    //     						   proposedStep,
    //     						   safety,
    //     						   &fGPILSelection);
  }

}

//-----------------------------------------------------------------------------
//  Transportation engine - Photon
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//  cuda wrapper for CPU
//-----------------------------------------------------------------------------

void tracking_photon_cpu(G4VPhysicalVolume *world, G4Navigator &aNavigator,
                         G4MagneticField *magMap, 
			 G4Track **track, 
			 G4PhysicsTable* eBrem_table, G4PhysicsTable* eIoni_table, G4PhysicsTable* msc_table,
			 size_t nTrackSize)
{
   // // G4MagneticField magField;
   // G4MagneticField_Constructor(&magField,magMap);
   
   G4Mag_UsualEqRhs equaOfMotion(magMap);
   // G4EquationOfMotion equaOfMotion(magMap);
   // G4EquationOfMotion_Constructor(&equaOfMotion,&magField,-1.0);
   
   G4ClassicalRK4 rk4( &equaOfMotion );
   // G4ClassicalRK4_Constructor(&rk4,&equaOfMotion);
   
   G4MagInt_Driver *magDriver = new G4MagInt_Driver(1.0, &rk4);
   //G4MagInt_Driver_Constructor(&magDriver,1.0,&rk4);
   
   G4ChordFinder chordFinder(magDriver);
   //G4ChordFinder_Constructor(&chordFinder,&magDriver);
   
   G4FieldManager aFieldManager(magMap,&chordFinder, true);
   //G4FieldManager_Constructor(&aFieldManager,&magField,&chordFinder, true);
   
   // //Navigator
   // G4Navigator aNavigator;
   // aNavigator.SetWorldVolume(world);
   // //G4Navigator_Constructor(&aNavigator);
   // //G4Navigator_SetWorldVolume(&aNavigator,world);
   
   //G4MultiLevelLocator
   G4MultiLevelLocator mLocator(&aNavigator);
   
   //Propagator
   G4PropagatorInField propagatorInField(&aNavigator,&aFieldManager,&mLocator);
   
   //Transporation
   G4Transportation transp(0);
   transp.SetPropagatorInField(&propagatorInField);
   
   //EMPhysics
   static G4eBremsstrahlung brem; // GXBrem brem;
   static G4eIonisation ioni; // GXIoni ioni;
   static G4eMultipleScattering msc; // GXMsc msc;
   G4ForceCondition condition;
   G4GPILSelection  selection;
   
   G4ProductionCutsTable::GetProductionCutsTable()->UpdateCoupleTable(world);
   static int init = EMPhysics_Init(0,NULL,&brem,&ioni,&msc,
                                    eBrem_table,eIoni_table,msc_table);
   
   // G4double x,y,z;
   G4double step = 0.0;
   G4double proposedStep = 0.0;
   
   // G4Track atrack;
   
   //  #pragma omp parallel for
   for (size_t tid = 0 ; tid < nTrackSize ; tid++) {
      
      //Geometry - Initialize Navigator and construct related structures
      aNavigator.LocateGlobalPointAndSetup(track[tid]->GetPosition(),
                                           NULL,false,true);
      // G4Navigator_LocateGlobalPointAndSetup(&aNavigator,G4ThreeVector_create(x,y,z),
      //      				  NULL,false,true);
      
      // G4MultiLevelLocator_Constructor(&mLocator, &aNavigator);
      
      // G4PropagatorInField_Constructor(&propagatorInField,&aNavigator,
      //     			    &aFieldManager,&mLocator); 
      
      // G4Transportation_Constructor(&transp,&propagatorInField,0);
      
      //EM Physics
      G4double safety = 0.0;
      proposedStep = EMPhysics_DefineStepLength(&brem,&ioni,&msc,track[tid],
                                                &condition,&selection,safety);
      
      //Transportation
      proposedStep += step;
      G4GPILSelection  fGPILSelection;
      
      //add energyLoss here as faking output
      track[tid]->SetStepLength(transp.AlongStepGetPhysicalInteractionLength(*track[tid],
                                                                             1.0,
                                                                             proposedStep,
                                                                             safety,
                                                                             &fGPILSelection));
      // track[tid].s = G4Transportation_AlongStepGPIL_Photon(&transp,
      //     						 &atrack,
      //     						 1.0,
      //     						 proposedStep,
      //     						 safety,
      //     						 &fGPILSelection);
   }
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
