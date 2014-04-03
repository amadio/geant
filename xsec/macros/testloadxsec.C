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
	TGeoManager *geom;
//	geom = TGeoManager::Import("Al_H2O_H.root");
	geom = TGeoManager::Import("http://root.cern.ch/files/cms.root");
	TTabPhysMgr::Instance(geom, "xsec_FTFP_BERT.root", "fstate_FTFP_BERT.root" );
   delete geom;
   delete TTabPhysMgr::Instance();
}
