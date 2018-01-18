#include <cassert>

#include "GUFieldPropagatorPool.h"

// For implementation
#include "GUFieldPropagator.h"
// #include "ScalarIntegrationDriver.h"
#include "VScalarEquationOfMotion.h"
#include "VVectorField.h"

#include <iostream>

// static
std::vector<GUFieldPropagator*> // GUFieldPropagation:: // namespace ...
   GUFieldPropagatorPool::fFieldPropagatorVec;

/// --------------  GUFieldPropagatorPool ------------------------------------
// #include "GUFieldPropagatorPool.h"   // For now, not a separate file

// static
GUFieldPropagatorPool* 
GUFieldPropagatorPool::Instance()
{
   // A lock is REQUIRED for the next line - TODO
   static GUFieldPropagatorPool sInstance;

   return &sInstance;
}

GUFieldPropagatorPool::GUFieldPropagatorPool( GUFieldPropagator* prototype )
   : fInitialisedRKIntegration(false),
     fNumberPropagators(0),
     fPrototype(prototype),
     fFieldPrototype(nullptr)
{
   // prototype can be null initially
}

GUFieldPropagatorPool::~GUFieldPropagatorPool()
{
   delete fPrototype;
}

bool
GUFieldPropagatorPool::RegisterPrototype( GUFieldPropagator* prototype )
{
   bool ok = ((fNumberPropagators > 0)
              && (fPrototype) && (prototype!=fPrototype) ); 
   if( !ok){
      std::cerr << "WARNING from GUFieldPropagatorPool:  "
                << "Changing prototype propagator after having created "
                << fNumberPropagators << " instances. " << std::endl;
      std::cerr << "     prototype =   "  << prototype
                << " old-prototype= " << fPrototype << std::endl;
   }
   assert( prototype );
   fPrototype = prototype;

   fFieldPrototype= prototype->GetField();
  
   fInitialisedRKIntegration = true;
   return ok;
}

bool
GUFieldPropagatorPool::Initialize( unsigned int numThreads )
{
   if( ! fPrototype ){
       std::cerr << "ERROR> from GUFieldPropagatorPool::Initialize:  "
                 << "Must register prototype propagator before calling Initialize. "
                 << std::endl
                 << "    # propagators= " << fNumberPropagators << " instances. "
                 << std::endl;
       exit(1);
   }
   bool goodExpansion= true;
   if( numThreads > fNumberPropagators )
   {
      // std::cout << "GUFieldPropagatorPool::Initialize  calling Extend for "
      //           << numThreads - fNumberPropagators << " new propagators. " << std::endl;
      Extend( numThreads - fNumberPropagators );
   }

   size_t revSize= fFieldPropagatorVec.size();
   std::cout << " Pool:  revised size= " << revSize
             << " request= " << numThreads << std::endl;

   goodExpansion = ( fFieldPropagatorVec.size() >= numThreads );
   assert (goodExpansion);

   fNumberPropagators= revSize;
   
   return (fPrototype != 0) && goodExpansion;
}

void
GUFieldPropagatorPool::Extend(size_t noNeeded)
{
    size_t num= fFieldPropagatorVec.size();
    // size_t originalNum= num;
    // const char* method= "GUFieldPropagatorPool::Extend";
    assert( fPrototype );
    assert( num < noNeeded );

    // printf("%s method called.  Num needed = %ld,  existing= %ld\n", method,
    //        noNeeded, num );

    while ( num < noNeeded )
    {
      auto prop= fPrototype->Clone();
      fFieldPropagatorVec.push_back( prop );

      num++;

      // printf("            Created propagator %p for slot %ld\n", prop, num );

      // fFieldVec.push_back( fFieldPrototype->CloneOrSafeSelf() );

      // Extension idea - also keep fields in 
      // auto field= prop->GetField();
      // fFieldVec.push_back( field );
    }
    // printf("%s method ended.  Created %ld propagators.  New total= %ld\n", method,
    //     num - originalNum, num );
}

#if 0

//// ---------------------  Postpone handling of multiple 
GUFieldPropagator* 
GUFieldPropagatorPool::CreateOrFind( int noNeeded ) // , void** banks )
{
  static int numberCreated= -1;
  static GUFieldPropagatorPool* pInstance= Instance();

  // A lock is REQUIRED for this section - TODO
  if( numberCreated < noNeeded)
  {
    Extend(noNeeded);
    assert( fFieldPropagatorVec.size() == noNeeded );
    // fNum = fFieldPropagatorVec.size();
    numberCreated= noNeeded;
  }
}
#endif
