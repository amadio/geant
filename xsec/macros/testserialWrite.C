#include "TSystem.h"
#include "TGeoManager.h"
#include "TEXsec.h"
#include "TEFstate.h"

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

   char *b=nullptr;
   size_t sizex = TEXsec::MakeCompactBuffer(b);
   cout << "Size of the X-sec buffer = " << sizex << " bytes " << endl;

   { // write to file
      std::ofstream fout("xsec.bin", std::ios::binary);
      fout.write(reinterpret_cast<char*>(&sizex), sizeof(sizex));
      int nelem = TEXsec::NLdElems();
      fout.write(reinterpret_cast<char*>(&nelem), sizeof(nelem));
      fout.write(reinterpret_cast<char*>(b), sizex);
      fout.close();
   }
   
   delete [] b;

   char *d=nullptr;
   size_t sizef = TEFstate::MakeCompactBuffer(d);
   cout << "Size of the fin state buffer = " << sizef << " bytes " << endl;

   { // write to file
      std::ofstream fout("fins.bin", std::ios::binary);
      fout.write(reinterpret_cast<char*>(&sizef), sizeof(sizef));
      int nelem = TEXsec::NLdElems();
      fout.write(reinterpret_cast<char*>(&nelem), sizeof(nelem));
      fout.write(reinterpret_cast<char*>(d), sizef);
      fout.close();
   }
   
   delete [] d;
   
   delete geom;
   delete TTabPhysMgr::Instance();
}
