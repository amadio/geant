#include "backend/cuda/Interface.h"

#include "TPartIndex.h"
#include "materials/Particle.h"
using vecgeom::Particle; 
#include "TEXsec.h"
#include "TPDecay.h"

#include "TEFstate.h"
#include "base/RNG.h"
using vecgeom::RNG;
#define UNIFORM() RNG::Instance().uniform()
#define FNPROC 18
/*
#include "TTabPhysMgr.h"
*/

constexpr int kIEL = 1;
__global__
void expandPhysics(char *buf, double *idSampled, int* idPart, float* idEnergy, int nrep) {
   printf("Rebuild TPartIndex class\n");
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
   printf("Rebuilding final state store\n");
   TEFstate::RebuildStore(buf);
   int sizef = TEFstate::SizeOfStore();
   printf("Number of bytes for final state %d\n",sizef);
//   for(auto iel=0; iel<TEXsec::NLdElems(); ++iel) 
   for(int irep=0; irep<nrep; irep++) { 
       idPart[irep] = 1;
       idEnergy[irep] = 1.0;
       printf("%f\n",idSampled[irep]);
   }
   for(int irep=0; irep<nrep; irep++) {

     idPart[irep] = (int) (idSampled[irep] * TPartIndex::I()->NPartReac());
     int ireac = idSampled[irep] * FNPROC;
        // int ireac = idSampled[irep] * 18;
     float en =  idSampled[irep] * (TPartIndex::I()->Emax() - TPartIndex::I()->Emin())
	    + TPartIndex::I()->Emin();
          //cout<<"using RNG "<<ipart<<endl;
     float xs = TEXsec::Element(kIEL)->XS(idPart[irep], ireac, en);
     if(xs < 0) continue;
     int npart=0;
     float weight=0;
     float kerma=0;
     float enr = 0.;
     const int *pid=0; 
     const float *mom=0;
     //int ebinindx=0;
     TEFstate::Element(kIEL)->GetReac(idPart[irep], ireac, en, 1,npart, weight, kerma, enr, pid, mom);
     //TEFstate::Element(iel)->SampleReac(idPart[irep], ireac, en, npart, weight, kerma, enr, pid, mom, ebinindx);
     idEnergy[irep]=enr;
     if (npart>0)
        printf("idPart %d, PDG %d, xs %f,ireac %d,  energy %f, npart %d, enr %f, pid %d, mom %f %f %f \n", idPart[irep], TPartIndex::I()->PDG(idPart[irep]), xs,ireac, en, npart, enr, pid[0],mom[0],mom[1],mom[2]);
       //  printf("dev energy for part %s is %f\n", TPartIndex::I()->PartName(idPart[irep]), idEnergy[irep]);
	// if(npart <= 0) continue;
  }
 // }
  return;
}
namespace vecgeom {
namespace cxx {

template size_t DevicePtr<char>::SizeOf();
template void DevicePtr<char>::Construct() const;

} // End cxx namespace
}

void launchExpandPhysicsOnDevice(vecgeom::cxx::DevicePtr<char> &devBuf, int nBlocks, int nThreads, double *idSampled, int* idPart, float* idEnergy, int nrep) {
//void launchExpandPhysicsOnDevice(vecgeom::cxx::DevicePtr<char> &devBuf, int nBlocks, int nThreads, int nrep, vecgeom::cxx::DevicePtr<double> idPart, vecgeom::cxx::DevicePtr<double> idEnergy) {
   int threadsPerBlock = nThreads;
   int blocksPerGrid   = nBlocks;
   printf("Launching expandPhysics threads: %d, blocks: %d\n",nThreads,nBlocks);

   expandPhysics<<< blocksPerGrid, threadsPerBlock >>>(devBuf, idSampled, idPart, idEnergy,nrep);
   cudaDeviceSynchronize();
}


