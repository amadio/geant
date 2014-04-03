#ifndef GPStepStatus_h
#define GPStepStatus_h 1

enum GPStepStatus
{
  fWorldBoundary,           
    // Step reached the world boundary
  fGeomBoundary,            
    // Step defined by a geometry boundary
  fAtRestDoItProc,          
    // Step defined by a PreStepDoItVector
  fAlongStepDoItProc,       
    // Step defined by a AlongStepDoItVector
  fPostStepDoItProc,        
    // Step defined by a PostStepDoItVector
  fUserDefinedLimit,
    // Step defined by the user Step limit in the logical volume
  fExclusivelyForcedProc,   
    // Step defined by an exclusively forced PostStepDoIt process 
  fUndefined                
    // Step not defined yet
};

#endif


