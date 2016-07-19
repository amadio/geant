#include "Particle.h"
namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

GEANT_CUDA_BOTH_CODE
void CreateParticle0000();
GEANT_CUDA_BOTH_CODE
void CreateParticle0001();
GEANT_CUDA_BOTH_CODE
void CreateParticle0002();
GEANT_CUDA_BOTH_CODE
void CreateParticle0003();
GEANT_CUDA_BOTH_CODE
void CreateParticle0004();
GEANT_CUDA_BOTH_CODE
void CreateParticle0005();
GEANT_CUDA_BOTH_CODE
void CreateParticle0006();
GEANT_CUDA_BOTH_CODE
void CreateParticle0007();
GEANT_CUDA_BOTH_CODE
void CreateParticle0008();
GEANT_CUDA_BOTH_CODE
void CreateParticle0009();
GEANT_CUDA_BOTH_CODE
void CreateParticle0010();
GEANT_CUDA_BOTH_CODE
void CreateParticle0011();
GEANT_CUDA_BOTH_CODE
void CreateParticle0012();
GEANT_CUDA_BOTH_CODE
void CreateParticle0013();
GEANT_CUDA_BOTH_CODE
void CreateParticle0014();
GEANT_CUDA_BOTH_CODE
void CreateParticle0015();
GEANT_CUDA_BOTH_CODE
void CreateParticle0016();
GEANT_CUDA_BOTH_CODE
void CreateParticle0017();
GEANT_CUDA_BOTH_CODE
void CreateParticle0018();
GEANT_CUDA_BOTH_CODE
void CreateParticle0019();
GEANT_CUDA_BOTH_CODE
void CreateParticle0020();

#ifdef GEANT_NVCC
GEANT_CUDA_DEVICE_CODE bool fgCreateParticlesInitDoneDev = false;
#endif

//________________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void Particle::CreateParticles() {
#ifndef GEANT_CUDA_DEVICE_BUILD
   static bool fgCreateParticlesInitDone = false;
#else
   bool &fgCreateParticlesInitDone(fgCreateParticlesInitDoneDev);
#endif
   if(fgCreateParticlesInitDone) return;
   fgCreateParticlesInitDone = true;
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
