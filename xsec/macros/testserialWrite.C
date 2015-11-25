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
   TRandom *r = new TRandom(12345);
   std::ofstream ftest("xsecsW.txt");
   for(auto iel=0; iel<TEXsec::NLdElems(); ++iel) {
      for(auto irep=0; irep<nrep; ++irep) {
	 // Get a random particle & reaction & energy
	 int ipart = r->Uniform() * TPartIndex::I()->NPartReac();
	 int ireac = r->Uniform() * FNPROC;
	 float en = r->Uniform() * (TPartIndex::I()->Emax() - TPartIndex::I()->Emin())
	    + TPartIndex::I()->Emin();
	 float xs = TEXsec::Element(iel)->XS(ipart, ireac, en);
	 if(xs < 0) continue;
	 ftest <<  iel << " " << TPartIndex::I()->PartName(ipart) << " " << ireac << " " << en << " " << xs << std::endl;
      }
   }
   ftest.close();


   char *b=nullptr;
   size_t sizex = TEXsec::MakeCompactBuffer(b);
   cout << "Size of the X-sec buffer = " << sizex << " bytes " << endl;

   // write to file
   std::ofstream fxsecs("xsec.bin", std::ios::binary);
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
