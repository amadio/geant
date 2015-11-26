#include "TSystem.h"
#include "TGeoManager.h"
#include "TEXsec.h"
#include "TEFstate.h"
#include "TPDecay.h"
#include "TRandom.h"

#include "TTabPhysMgr.h"
#include "GeantPropagator.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

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
   srand(12345);
   std::ofstream fxtest("xsecsW.txt");
   for(auto iel=0; iel<TEXsec::NLdElems(); ++iel) {
      for(auto irep=0; irep<nrep; ++irep) {
	 // Get a random particle & reaction & energy
	 int ipart = (((double) rand())/RAND_MAX) * TPartIndex::I()->NPartReac();
	 int ireac = (((double) rand())/RAND_MAX) * FNPROC;
	 float en = (((double) rand())/RAND_MAX) * (TPartIndex::I()->Emax() - TPartIndex::I()->Emin())
	    + TPartIndex::I()->Emin();
	 float xs = TEXsec::Element(iel)->XS(ipart, ireac, en);
	 if(xs < 0) continue;
	 fxtest <<  iel << " " << TPartIndex::I()->PartName(ipart) << " " << ireac << " " << en << " " << xs << std::endl;
      }
   }
   fxtest.close();

   std::ofstream fftest("xfinsW.txt");
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

   char *b=nullptr;

   // write to file
   std::ofstream fxsecs("xsec.bin", std::ios::binary);
   size_t sizex = TPartIndex::I()->MakeCompactBuffer(b);
   cout << "Size of the TPartIndex buffer = " << sizex << " bytes " << endl;
   fxsecs.write(reinterpret_cast<char*>(&sizex), sizeof(sizex));
   fxsecs.write(reinterpret_cast<char*>(b), sizex);

   delete [] b;

   sizex = TEXsec::MakeCompactBuffer(b);
   cout << "Size of the X-sec buffer = " << sizex << " bytes " << endl;
   fxsecs.write(reinterpret_cast<char*>(&sizex), sizeof(sizex));
   int nelem = TEXsec::NLdElems();
   fxsecs.write(reinterpret_cast<char*>(&nelem), sizeof(nelem));
   fxsecs.write(reinterpret_cast<char*>(b), sizex);

   fxsecs.close();
   
   delete [] b;

   char *d=nullptr;
   size_t sizef = TEFstate::MakeCompactBuffer(d);
   cout << "Size of the fin state buffer = " << sizef << " bytes " << endl;

   // write to file
   std::ofstream ffinst("fins.bin", std::ios::binary);
   ffinst.write(reinterpret_cast<char*>(&sizef), sizeof(sizef));
   ffinst.write(reinterpret_cast<char*>(&nelem), sizeof(nelem));
   ffinst.write(reinterpret_cast<char*>(d), sizef);
   
   delete [] d;

   // now the decay table... dirty... 
   TPDecay *decay = TTabPhysMgr::Instance()->GetDecayTable();
   sizef = decay->MakeCompactBuffer(d);
   cout << "Size of the decay buffer = " << sizef << " bytes " << endl;
   // write to file
   ffinst.write(reinterpret_cast<char*>(&sizef), sizeof(sizef));
   ffinst.write(reinterpret_cast<char*>(d), sizef);
   ffinst.close();

   delete [] d;

   delete geom;
   delete TTabPhysMgr::Instance();
}
