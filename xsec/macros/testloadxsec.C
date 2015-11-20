#include "TSystem.h"
#include "TGeoManager.h"
#include "TEXsec.h"

#include "TTabPhysMgr.h"
#include "GeantPropagator.h"

#include <iostream>

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
   geom = TGeoManager::Import("http://root.cern.ch/files/cms.root");
   TTabPhysMgr::Instance(fxsec, ffins );
   for(auto i=0; i<TEXsec::NLdElems(); ++i)
      TEXsec::Element(i)->GetPartSize();
   char *b=nullptr;
   size_t size = TEXsec::MakeCompactBuffer(b);
   cout << "Size of the buffer = " << size << " bytes " << endl;
   TEXsec::RebuildStore(size,TEXsec::NLdElems(),b);
   delete geom;
   delete TTabPhysMgr::Instance();
}
