#include "GXBrem.h"
#include "GXIoni.h"
#include "GXMsc.h"

#include "GXTrack.h"
#include "util_kernel.h"

//-----------------------------------------------------------------------------
//  EMPhysics 
//-----------------------------------------------------------------------------
FQUALIFIER
void EMPhysics_Init(unsigned int tid,
		    curandState* devStates,
		    GXBrem *brem,
		    GXIoni *ioni,
		    GXMsc  *msc,
		    GPPhysicsTable* eBrem_table, 
		    GPPhysicsTable* eIoni_table, 
		    GPPhysicsTable* msc_table) 
{
  //Initialize EMPhysics Processes
  brem->SetCurandState(&devStates[tid]);
  brem->UseIntegral(true);
  brem->SetLambdaTable(eBrem_table);
  brem->UseLambdaTable(true);

  ioni->SetCurandState(&devStates[tid]);
  ioni->UseIntegral(true);
  ioni->SetLambdaTable(eIoni_table);
  ioni->UseLambdaTable(true);

  msc->SetCurandState(&devStates[tid]);
  msc->UseIntegral(true);
  msc->SetLambdaTable(msc_table);
  msc->UseLambdaTable(true);

}

FQUALIFIER
G4double EMPhysics_DefineStepLength(GXBrem *brem,
				GXIoni *ioni,
				GXMsc  *msc,
				GXTrack *atrack,
                                GPForceCondition *condition,
				GPGPILSelection  *selection) 
{
  G4double px,py,pz,step,E;

  G4double energyLoss = 0.0;
  G4double step_brem = 0.0;
  G4double step_ioni = 0.0;
  G4double step_msc  = 0.0;
  G4double proposedStep = 0.0;
  //EM Physics
  px   = atrack->px;
  py   = atrack->py;
  pz   = atrack->pz;
  step = atrack->s;

  E = sqrt(px*px+py*py+pz*pz+0.510998910*0.510998910);//calcuate 

  step_brem = brem->PostStepGetPhysicalInteractionLength(E, step, condition);
  step_ioni = ioni->PostStepGetPhysicalInteractionLength(E, step, condition);
  step_msc = msc->PostStepGetPhysicalInteractionLength(E, step, condition);

  //physics model defining the current step
  unsigned int model =0;
  if(step_brem > step_ioni) model = 1;   
  if(step_ioni > step_msc) model = 2;  

  switch (model) {
  case 0 : 
    proposedStep = step_brem;
    energyLoss += brem->PostStepDoIt(atrack);
    break; 
  case 1 :
    proposedStep = step_ioni;
    energyLoss += ioni->PostStepDoIt(atrack);
    break; 
  case 2 :
    proposedStep = step_msc;                         
    energyLoss += msc->PostStepDoIt(atrack);
    break;
  default:
    printf("oups dont know the model, we have %d %f %f %f\n",step_brem,step_ioni,step_msc);
  }

  //alongstep
  step_msc = msc->AlongStepGetPhysicalInteractionLength(atrack, step, selection);
  energyLoss += msc->AlongStepDoIt(atrack);
  if(step_msc < proposedStep) proposedStep = step_msc;

  return proposedStep;
}
