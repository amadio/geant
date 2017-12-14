#include "TSystem.h"
#include "TGeoManager.h"
#include "TEXsec.h"
#include "TEFstate.h"
#include "TPDecay.h"
#include "TRandom.h"

#include "TTabPhysMgr.h"
#include "GeantPropagator.h"

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

using std::cout;
using std::endl;

int serializePhysics(char *&buf);

void testserialWrite()
{
   gSystem->Load("libPhysics.so");
   gSystem->Load("libHist.so");
   gSystem->Load("libThread.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("../../lib/libXsec");
   gSystem->Load("../../lib/libGeant_v.so");
   gSystem->Load("../../lib/libUser.so");
   TGeoManager *geom = nullptr;
   //	geom = TGeoManager::Import("Al_H2O_H.root");
   const char *fxsec = "../../data/xsec_FTFP_BERT_G496p02_1mev.root";
   const char *ffins = "../../data/fstate_FTFP_BERT_G496p02_1mev.root";
   GeantPropagator::Instance(1,1,1);
   geom = TGeoManager::Import("http://root.cern.ch/files/cms.root");
   TTabPhysMgr::Instance(fxsec, ffins );

   constexpr int nrep = 1000;
   #ifndef USE_VECGEOM_NAVIGATOR
   #ifdef USE_ROOT
   gRandom->SetSeed(12345);
   #else
   srand(12345);
   #endif
   #endif
   std::ofstream fftest("xphysW.txt");
   for(auto iel=0; iel<TEXsec::NLdElems(); ++iel) {
      for(auto irep=0; irep<nrep; ++irep) {
	 // Get a random particle & reaction & energy
	    int ipart = UNIFORM() * TPartIndex::I()->NPartReac();
      int ireac = UNIFORM() * FNPROC;
	    float en =  UNIFORM() * (TPartIndex::I()->Emax() - TPartIndex::I()->Emin())
	    + TPartIndex::I()->Emin();
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

   char *buf = nullptr;
   int totsize = serializePhysics(buf);

   // write to file
   std::ofstream fxsecs("xfphys.bin", std::ios::binary);
   fxsecs.write(reinterpret_cast<char*>(&totsize), sizeof(totsize));
   fxsecs.write(reinterpret_cast<char*>(buf), totsize);
   fxsecs.close();

   delete [] buf;

   delete geom;
   delete TTabPhysMgr::Instance();
}

int serializePhysics(char *&buf) {
   TPDecay *decay = TTabPhysMgr::Instance()->GetDecayTable();
   int totsize = TPartIndex::I()->SizeOf() + TEXsec::SizeOfStore()
               + TEFstate::SizeOfStore() + decay->SizeOf();
   cout << "Total size of store " << totsize << endl;
   delete [] buf;
   buf= new char[totsize];
   char *cur = buf;
   int sizet = TPartIndex::I()->MakeCompactBuffer(buf);
   cout << "Size of the TPartIndex buffer = " << sizet << " bytes " << endl;
   cur += sizet;
   int sizex = TEXsec::MakeCompactBuffer(cur);
   cout << "Size of the X-sec buffer = " << sizex << " bytes " << endl;
   cur += sizex;
   int sized = decay->MakeCompactBuffer(cur);
   cout << "Size of the decay buffer = " << sized << " bytes " << endl;
   cur += sized;
   int sizef = TEFstate::MakeCompactBuffer(cur);
   cout << "Size of the fin state buffer = " << sizef << " bytes " << endl;
   cur += sizef;
   return totsize;
}
