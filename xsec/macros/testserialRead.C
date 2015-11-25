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

   char *b=nullptr;
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
      TEFstate::RebuildStore(nb,TEXsec::NLdElems(),d);
      fin.close();
   }

   const char *fxsec = "/dev/null";
   const char *ffins = "/dev/null";
   TTabPhysMgr::Instance(fxsec, ffins );

   delete geom;
   delete TTabPhysMgr::Instance();
}
