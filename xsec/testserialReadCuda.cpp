/*
#include "TEXsec.h"
#include "TEFstate.h"
#include "TPDecay.h"

#include "TTabPhysMgr.h"
#include "GeantPropagator.h"
*/
#ifdef USE_ROOT
#include "TGeoManager.h"
#endif
#include <iostream>
#include <fstream>

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

#include "backend/cuda/Interface.h"

void launchExpandPhysicsOnDevice(vecgeom::DevicePtr<char>, int nBlocks, int nThreads);

using std::cout;
using std::endl;

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

   std::cout << "Total size of store " << totsize << std::endl;
   launchExpandPhysicsOnDevice(devBuf, 1, 1);

  /*  
   expandPhysicsKernel<<<>>>(devBuf);

   const char *fxsec = "/dev/null";
   const char *ffins = "/dev/null";
   #ifdef USE_ROOT
   GeantPropagator::Instance(1,1,1);
   TGeoManager *geom = TGeoManager::Import("http://root.cern.ch/files/cms.root");

   #endif
   TTabPhysMgr::Instance(fxsec, ffins );

   constexpr int nrep = 1000;

   #ifndef USE_VECGEOM_NAVIGATOR
   #ifdef USE_ROOT
   gRandom->SetSeed(12345);
   #else
   srand(12345);
   #endif
   #endif

   std::ofstream fftest("xphysR.txt");
   for(auto iel=0; iel<TEXsec::NLdElems(); ++iel) {
      for(auto irep=0; irep<nrep; ++irep) {

	 int ipart = UNIFORM() * TPartIndex::I()->NPartReac();
         int ireac = UNIFORM() * FNPROC;
	 float en =  UNIFORM() * (TPartIndex::I()->Emax() - TPartIndex::I()->Emin())
	    + TPartIndex::I()->Emin();
          //cout<<"using RNG "<<ipart<<endl;
         float xs = TEXsec::Element(iel)->XS(ipart, ireac, en);
 	 if(xs < 0) continue;
	 int npart=0;
	 float weight=0;
	 float kerma=0;
	 float enr=0;
	 const int *pid=0;
	 const float *mom=0;
	 int ebinindx=0;
	 TEFstate::Element(iel)->SampleReac(ipart, ireac, en, npart, weight, kerma, enr, pid, mom, ebinindx);
	 if(npart <= 0) continue;
	 fftest <<  iel << ":" << TPartIndex::I()->PartName(ipart) << ":" << ireac << ":" << en
		<< ":" << xs << ":" << npart << ":" << weight << ":" << kerma << ":" << enr << ":";
	 for(auto i=0; i<npart; ++i)
	    fftest << pid[i] << ":" << mom[i*3] << ":" << mom[i*3+1] << ":" << mom[i*3+2];
	 fftest <<":" << ebinindx << std::endl;
      }
   }
   fftest.close();
   #ifdef USE_ROOT
   delete geom;
   #endif
*/
   return 0;
}
/*
void expandPhysicsKernel<<<>>>(char *buf) {
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
   std::cout << "Number of bytes for final state " << sizef << std::endl;
}
*/
