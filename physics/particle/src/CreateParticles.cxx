#include "Particle.h"
namespace geant {
   inline namespace GEANT_IMPL_NAMESPACE {

void CreateParticle0000();
void CreateParticle0001();
void CreateParticle0002();
void CreateParticle0003();
void CreateParticle0004();
void CreateParticle0005();
void CreateParticle0006();
void CreateParticle0007();
void CreateParticle0008();
void CreateParticle0009();
void CreateParticle0010();
void CreateParticle0011();
void CreateParticle0012();
void CreateParticle0013();
void CreateParticle0014();
void CreateParticle0015();
void CreateParticle0016();
void CreateParticle0017();
void CreateParticle0018();
void CreateParticle0019();
void CreateParticle0020();

//________________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void Particle::CreateParticles() {
   static bool initDone=false;
   if(initDone) return;
   initDone = true;
    CreateParticle0000();
    CreateParticle0001();
    CreateParticle0002();
    CreateParticle0003();
    CreateParticle0004();
    CreateParticle0005();
    CreateParticle0006();
    CreateParticle0007();
    CreateParticle0008();
    CreateParticle0009();
    CreateParticle0010();
    CreateParticle0011();
    CreateParticle0012();
    CreateParticle0013();
    CreateParticle0014();
    CreateParticle0015();
    CreateParticle0016();
    CreateParticle0017();
    CreateParticle0018();
    CreateParticle0019();
    CreateParticle0020();
}
 } // End of inline namespace
 } // End of geant namespace
#if defined(__clang__) && !defined(__APPLE__)
#pragma clang optimize on
#endif
