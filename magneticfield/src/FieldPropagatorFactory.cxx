#include "FieldPropagatorFactory.h"

// #include "Geant/Error.h"

// Source file is required to aid compiler/linker in placing inline methods into a library.

//______________________________________________________________________________
void
FieldPropagatorFactory::RegisterPropagator(GUFieldPropagator* fieldPropagator)
{
  GUFieldPropagatorPool* fpPool= GUFieldPropagatorPool::Instance();
  assert( fpPool );  // Cannot be zero
  if( fpPool ) {
     fpPool->RegisterPrototype( fieldPropagator );
     // printf( "FieldPropagatorFactory: Registered Prototype field-prop %p\n", fieldPropagator );
  } else {
     // Geant::Error("PrepareRkIntegration","Cannot find GUFieldPropagatorPool Instance.");
     std::cerr << "ERROR in PrepareRkIntegration: "
               << "Cannot find GUFieldPropagatorPool Instance." << std::endl;
  }
}
