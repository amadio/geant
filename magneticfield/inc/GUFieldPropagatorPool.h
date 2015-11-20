//  Create and maintain a 'pool' of Field Propagator instances
//  Each thread will use the same instance, indexed by its thread-id (0<tid<=maxTID)
//

//  For Multi-threaded version only -- ie not for CUDA

//  An implementation on GPU will require a different approach, potentially
//   - a dynamically created set of classes GUFieldPropagator/Driver/Stepper ...
//   - a static method for each class (revised GUFieldPropagator/Driver/..)
//        without the need for an object
//   - a mixture of the two approaches.
//  Current status: to be investigated, once classes are more stable.

#include <vector>
class GUFieldPropagator;
// #include "GUFieldPropagator.h"

class GUFieldPropagatorPool
{
  public:
    // Access methods
    // static GUFieldPropagator* CreateOrFind(int numThreads);
      // It can be called from many threads -- same value must be returned
      //  numThreads must be constant between calls

    static GUFieldPropagatorPool* Instance();

    bool   RegisterPrototype( GUFieldPropagator* prototype );
     // prototype for the propagators for each thread
    
    bool Initialize(unsigned int numThreads); 
     // Create new propagators for each thread !
    
    void CheckIndex(int num){
       assert(num>=0);
       assert(num< fFieldPropagatorVec.size());
    }
    
    GUFieldPropagator* GetPropagator(int num) {
       CheckIndex(num);       
       return fFieldPropagatorVec[num];
    }

  private:

    GUFieldPropagatorPool( GUFieldPropagator* prototype = 0); // , void** banks=0 );  // Ensure one per thread
    ~GUFieldPropagatorPool() {} 

    void Extend(unsigned int Num);
     // Create additional propagators, so that total is 'Num'
  private:
    unsigned int fNumberPropagators;
    const  GUFieldPropagator* fPrototype; //  
    
    static std::vector<GUFieldPropagator*> fFieldPropagatorVec;
};
