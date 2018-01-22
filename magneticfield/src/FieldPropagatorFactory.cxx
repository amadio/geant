#include "FieldPropagatorFactory.h"

// #include "Geant/Error.h"
bool  FieldPropagatorFactory::fVerboseConstruct= false;

// Source file is required to aid compiler/linker in placing inline methods into a library.

//______________________________________________________________________________
void
FieldPropagatorFactory::RegisterPropagator(GUFieldPropagator* fieldPropagator)
{
  GUFieldPropagatorPool* fpPool= GUFieldPropagatorPool::Instance();
  assert( fpPool );  // Cannot be zero
  if( fpPool ) {
     fpPool->RegisterPrototype( fieldPropagator );
     // Not complete until   fpPool->Initialize( numThreads ); is called
     // Geant::Printf( "FieldPropagatorFactory: Registered Prototype field-prop %p\n", fieldPropagator );
     std::cout << "FieldPropagatorFactory: Registered Prototype field-prop " << fieldPropagator
               << std::endl;
  } else {
     // Geant::Error("PrepareRkIntegration","Cannot find GUFieldPropagatorPool Instance.");
     std::cerr << "ERROR in PrepareRkIntegration: "
               << "Cannot find GUFieldPropagatorPool Instance." << std::endl;
  }
}
