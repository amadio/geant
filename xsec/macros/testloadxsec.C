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

void testloadxsec()
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
   //   geom = TGeoManager::Import("http://root.cern.ch/files/cms.root");
   geom = TGeoManager::Import("cms.root");
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

   
   { // read from file
      std::ifstream fin("xsec.bin", std::ios::binary);
      size_t nb;
      fin.read(reinterpret_cast<char*>(&nb), sizeof(nb));
      std::cout << "Number of bytes for x-sec " << nb << std::endl;
      int nel;
      fin.read(reinterpret_cast<char*>(&nel), sizeof(nel));
      std::cout << "Number of elements in x-sec " << nel << std::endl;
      b = (char *) malloc(nb);
      fin.read(reinterpret_cast<char*>(b), nb);
      std::cout << "Rebuilding x-sec store" << std::endl;
      TEXsec::RebuildStore(nb,nel,b);
      fin.close();
   }

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

   
   { // read from file
      std::ifstream fin("fins.bin", std::ios::binary);
      size_t nb;
      fin.read(reinterpret_cast<char*>(&nb), sizeof(nb));
      std::cout << "Number of bytes for fins " << nb << std::endl;
      int nel;
      fin.read(reinterpret_cast<char*>(&nel), sizeof(nel));
      std::cout << "Number of elements in fins " << nel << std::endl;
      d = (char *) malloc(nb);
      fin.read(reinterpret_cast<char*>(d), nb);
      std::cout << "Rebuilding fins store" << std::endl;
      TEFstate::RebuildStore(sizef,TEXsec::NLdElems(),d);
      fin.close();
   }



   delete geom;
   delete TTabPhysMgr::Instance();
}
