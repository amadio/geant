#include "TSystem.h"
#include "TGeoManager.h"
#include "TEXsec.h"
#include "TEFstate.h"
#include "TPDecay.h"
#include "TRandom.h"
#include "TFile.h"

#include "TTabPhysMgr.h"
#include "GeantPropagator.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

void testserialRead()
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
   GeantPropagator::Instance(1,1,1);
   geom = TGeoManager::Import("http://root.cern.ch/files/cms.root");

   /*
     const char *tpf = "../../data/xsec_FTFP_BERT_G496p02_1mev.root";
     TFile *f = new TFile(tpf);
     TPartIndex *p = (TPartIndex *) f->Get("PartIndex");
     f->Close();delete f;
   */

   char *b=nullptr;
   char *t=nullptr;
   { // read from file
      std::ifstream fin("xsec.bin", std::ios::binary);
      int nb;
      fin.read(reinterpret_cast<char*>(&nb), sizeof(nb));
      std::cout << "Number of bytes for TPartIndex " << nb << std::endl;
      t = (char *) malloc(nb);
      fin.read(reinterpret_cast<char*>(t), nb);
      fin.close();
      std::cout << "Rebuilding TPartIndex store" << std::endl;
      TPartIndex::I()->RebuildClass(t);
      t += TPartIndex::I()->SizeOf();
      std::cout << "Number of bytes for x-sec " << nb - TPartIndex::I()->SizeOf() << std::endl;
      std::cout << "Rebuilding x-sec store" << std::endl;
      TEXsec::RebuildStore(t);
   }

   char *d=nullptr;
   char *de=nullptr;
   { // read from file
      std::ifstream fin("fins.bin", std::ios::binary);
      int nb;
      fin.read(reinterpret_cast<char*>(&nb), sizeof(nb));
      std::cout << "Number of bytes for fins " << nb << std::endl;
      d = (char *) malloc(nb);
      fin.read(reinterpret_cast<char*>(d), nb);
      fin.close();
      std::cout << "Rebuilding fins store" << std::endl;
      TEFstate::RebuildStore(d);
      d += TEFstate::SizeOfStore();

      std::cout << "Number of bytes for decay " << nb - TEFstate::SizeOfStore() << std::endl;
      std::cout << "Rebuilding decay store" << std::endl;
      TPDecay *dec = (TPDecay *) d;
      dec->RebuildClass();
      TEFstate::SetDecayTable(dec);
   }

   const char *fxsec = "/dev/null";
   const char *ffins = "/dev/null";
   TTabPhysMgr::Instance(fxsec, ffins );

   std::cout << TPartIndex::I()->NPartReac() << std::endl;
   constexpr int nrep = 1000;
   srand(12345);
   std::ofstream fout("xsecsR.txt");
   for(auto iel=0; iel<TEXsec::NLdElems(); ++iel) {
      for(auto irep=0; irep<nrep; ++irep) {
	 // Get a random particle & reaction & energy
	 int ipart = (((double) rand())/RAND_MAX) * TPartIndex::I()->NPartReac();
	 int ireac = (((double) rand())/RAND_MAX) * FNPROC;
	 float en = (((double) rand())/RAND_MAX) * (TPartIndex::I()->Emax() - TPartIndex::I()->Emin())
	    + TPartIndex::I()->Emin();
	 float xs = TEXsec::Element(iel)->XS(ipart, ireac, en);
	 if(xs < 0) continue;
	 fout <<  iel << " " << TPartIndex::I()->PartName(ipart) << " " << ireac << " " << en << " " << xs << std::endl;
      }
   }
   fout.close();

   std::ofstream fftest("xfinsR.txt");
   for(auto iel=0; iel<TEXsec::NLdElems(); ++iel) {
      for(auto irep=0; irep<nrep; ++irep) {
	 // Get a random particle & reaction & energy
	 int ipart = (((double) rand())/RAND_MAX) * TPartIndex::I()->NPartReac();
	 int ireac = (((double) rand())/RAND_MAX) * FNPROC;
	 float en = (((double) rand())/RAND_MAX) * (TPartIndex::I()->Emax() - TPartIndex::I()->Emin())
	    + TPartIndex::I()->Emin();
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
		<< ":" << npart << ":" << weight << ":" << kerma << ":" << enr << ":";
	 for(auto i=0; i<npart; ++i)
	    fftest << pid[i] << ":" << mom[i*3] << ":" << mom[i*3+1] << ":" << mom[i*3+2];
	 fftest <<":" << ebinindx << std::endl;
      }
   }
   fftest.close();

   delete geom;
//   delete TTabPhysMgr::Instance();
}
