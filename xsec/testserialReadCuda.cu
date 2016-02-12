#include "backend/cuda/Interface.h"
//#include "TPartIndex.h"
#include "materials/Particle.h"
using vecgeom::Particle; 
/*
#include "TEXsec.h"
#include "TEFstate.h"
#include "TPDecay.h"

#include "TTabPhysMgr.h"
*/
__global__
void expandPhysics(char *buf) {
   printf("Rebuilding TPartIndex store\n");
   Particle::CreateParticles();
/*
   TPartIndex::I()->RebuildClass(buf);

   int sizet = TPartIndex::I()->SizeOf();
   printf("Number of bytes for TPartIndex %d\n",sizet);
   buf += sizet;
   printf("Rebuilding x-sec store\n");
   TEXsec::RebuildStore(buf);
   int sizex = TEXsec::SizeOfStore();
   printf("Number of bytes for x-sec %d\n",sizex);
   buf += sizex;
   printf("Rebuilding decay store\n");
   TPDecay *dec = (TPDecay *) buf;
   dec->RebuildClass();
   TEFstate::SetDecayTable(dec);
   int sized = dec->SizeOf();
   printf("Number of bytes for decay %d\n",sized);
   buf += sized;
   printf("Rebuilding final state store");
   TEFstate::RebuildStore(buf);
   int sizef = TEFstate::SizeOfStore();
   printf("Number of bytes for final state %d\n",sizef);
*/
}

namespace vecgeom {
namespace cxx {

template size_t DevicePtr<char>::SizeOf();
template void DevicePtr<char>::Construct() const;

} // End cxx namespace
}

void launchExpandPhysicsOnDevice(vecgeom::cxx::DevicePtr<char> devBuf, int nBlocks, int nThreads) {
 int threadsPerBlock = nThreads;
 int blocksPerGrid   = nBlocks;
 expandPhysics<<< blocksPerGrid, threadsPerBlock >>>(devBuf);

}
