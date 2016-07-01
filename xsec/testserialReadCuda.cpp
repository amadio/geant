
/*
#include "TTabPhysMgr.h"
#include "GeantPropagator.h"
*/

#ifdef USE_VECGEOM_NAVIGATOR
#include "base/RNG.h"
using vecgeom::RNG;
#define UNIFORM() RNG::Instance().uniform()
#elif USE_ROOT
#include <TRandom.h>
#define UNIFORM() gRandom->Uniform()
#else
#define UNIFORM() ((double)rand())/RAND_MAX
#endif

#include <iostream>
#include <fstream>
#include "backend/cuda/Interface.h"

constexpr unsigned int kNREP =100;
void launchExpandPhysicsOnDevice(vecgeom::DevicePtr<char>&, int nBlocks, int nThreads, double* iSampled,int* devIPart, float* devIEnergy, int kNREP);
/*
void expandPhysicsLocal(char *buf) {
   std::cout << "Rebuilding TPartIndex store" << std::endl;
   TPartIndex::I()->RebuildClass(buf);
   int sizet = TPartIndex::I()->SizeOf();
   std::cout << "Number of bytes for TPartIndex " << sizet << std::endl;
   buf += sizet;
   std::cout << "Rebuilding x-sec store" << std::endl;
   TEXsec::RebuildStore(buf);
   int sizex = TEXsec::SizeOfStore();
   std::cout << "Number of bytes for x-sec " << sizex << std::endl;
   buf += sizex;
   std::cout << "Rebuilding decay store" << std::endl;
   TPDecay *dec = (TPDecay *) buf;
   dec->RebuildClass();
   TEFstate::SetDecayTable(dec);
   int sized = dec->SizeOf();
   std::cout << "Number of bytes for decay " << sized << std::endl;
   buf += sized;
   std::cout << "Rebuilding final state store" << std::endl;
   TEFstate::RebuildStore(buf);
   int sizef = TEFstate::SizeOfStore();
   std::cout << "Number of bytes for final state -- HOST --" << sizef << std::endl;
}
*/
int main()
{
   char *hostBuf=nullptr;
   int totsize;
   // read from file
   std::ifstream fin("xfphys.bin", std::ios::binary);
   fin.read(reinterpret_cast<char*>(&totsize), sizeof(totsize));
   //buf = new char[totsize];
   hostBuf = (char*)_mm_malloc(totsize,sizeof(double));
   fin.read(reinterpret_cast<char*>(hostBuf), totsize);
   fin.close();
   vecgeom::DevicePtr<char> devBuf;
   devBuf.Allocate(totsize);
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR ALLOC buffer\n");
      return 0;
   } 
   devBuf.ToDevice(hostBuf,totsize);
   if (cudaSuccess!=cudaGetLastError()) {
      printf("ERROR MEMCPY buffer\n");
      return 0;
   }

   printf("Total size of store %d\n", totsize);

   double *iSampled;
   float *iEnergy; 
   int *iPart; 
   cudaMallocHost(&iSampled, kNREP*sizeof(double));
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating iSample on Host: %s \n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }
   cudaMallocHost(&iEnergy, kNREP*sizeof(float));
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating iEnergy on Host: %s \n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }
   cudaMallocHost(&iPart, kNREP*sizeof(int));
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating iPart on Host: %s \n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }

   for(unsigned int irep=0; irep<kNREP; irep++) {
      iSampled[irep] = (double) UNIFORM();
   }

   double  *devISampled;
   int  *devIPart;
   float  *devIEnergy;

   cudaMalloc((void**) &devISampled, kNREP*sizeof(double)); 
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating devIPart to Device: %s \n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }
   cudaMalloc((void**) &devIPart, kNREP*sizeof(int)); 
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating devIPart to Device: %s \n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }
   cudaMalloc((void**) &devIEnergy, kNREP*sizeof(float)); 
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR allocating devIEnergy to Device: %s\n",cudaGetErrorString(cudaGetLastError()));
      return 0;
   }
 
   cudaMemcpy(devISampled,iSampled, kNREP*sizeof(double),cudaMemcpyHostToDevice);
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR copying iSampled to devIPart on Device:%s \n", cudaGetErrorString(cudaGetLastError()));
      return 0;
   }

   launchExpandPhysicsOnDevice(devBuf, 1, 1,devISampled, devIPart,devIEnergy,kNREP);

   cudaThreadSynchronize();
   cudaMemcpy(iSampled,devISampled, kNREP*sizeof(double),cudaMemcpyDeviceToHost);
   cudaError_t error=cudaGetLastError(); 

   if (error != cudaSuccess) {
      printf(" ERROR copy iSampled from Device: %s\n", cudaGetErrorString(error));
      return 0;
   }

   cudaMemcpy(iPart,devIPart, kNREP*sizeof(int),cudaMemcpyDeviceToHost);
   error=cudaGetLastError(); 
   if (error != cudaSuccess) {
      printf(" ERROR copy iPart from Device: %s\n", cudaGetErrorString(error));
      return 0;
   }
   cudaMemcpy(iEnergy,devIEnergy, kNREP*sizeof(float),cudaMemcpyDeviceToHost);
   if (cudaGetLastError() != cudaSuccess) {
      printf(" ERROR copy iEnergy from Device\n");
      return 0;
   }
    

   
   cudaMemcpy(iEnergy,devIEnergy, kNREP*sizeof(float),cudaMemcpyDeviceToHost);
   cudaMemcpy(iPart,devIPart, kNREP*sizeof(float),cudaMemcpyDeviceToHost);

   std::ofstream fftest("xphysR.txt");

   for(unsigned int irep=0; irep<kNREP; ++irep)
      if (iEnergy[irep] >0.)
        fftest << "idPart "<<  iPart[irep] << ", energy " << iEnergy[irep] <<std::endl;
   fftest.close();

   return 0;
}
