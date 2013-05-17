/*
 int nphi = 4;
int nz   = 3;
double density = 8.28;
*/

TGeoVolume *VP_SimpleECal(int nphi = 4, int nz = 3, double density = 8.28)
{
   const double world_x =  9000.;
   const double world_y =  9000.;
   const double world_z = 16000;

   TGeoManager *geom = new TGeoManager("SimpleECal","Simplification of the CMS ECal");

   TGeoMaterial *world_mat = new TGeoMaterial("world",1,2,0);
   TGeoMedium *world_med=new TGeoMedium("world_medium",0,world_mat);
   TGeoVolume *world = geom->MakeBox("top",world_med,world_x, world_y, world_z);
   
   geom->SetTopVolume(world);
   geom->SetTopVisible(1);
   
   int crystal_nphi  = nphi;
   int crystal_nz    = nz;
   double ecal_density = density;

   const int crystal_n     = crystal_nphi*crystal_nz;
   const double ecal_zmin  = -3000.;
   const double ecal_zmax  =  3000.;
   
   const double ecal_rmin  =  10.;
   const double ecal_rmax  =  5000.;
   const double ecal_dz    =  0.5*(ecal_zmax-ecal_zmin)/crystal_nz;
   const double ecal_sphi  =     0.;
   //    const G4double ecal_dphi  =  2.0*M_PI/crystal_nphi;
   // G4 seems to be in radian while TGeo seems to be in degree.
   // const double ecal_dphi  =  2.0*TMath::Pi()/crystal_nphi;
   const double ecal_dphi  =  2.0*180/crystal_nphi;

   int iptr = 0;
   
   TGeoElementTable *table = gGeoManager->GetElementTable();
   //TGeoElement* elPb = new TGeoElement( "Lead", "Pb", 82., 207.19*g/mole );
   //TGeoElement* elW = new TGeoElement( "Tungstenm", "W",74., 183.85*g/mole);
   //TGeoElement* elO = new TGeoElement( "Oxygen", "O2", 8., 16.*g/mole );
   TGeoElement* elPb = table->GetElement(82);
   TGeoElement* elW = table->GetElement(74);
   TGeoElement* elO = table->GetElement(8);
   //TGeoMaterial *ecal_mat = new TGeoMaterial("ecal",90,120,density);
   TGeoMixture *ecal_mat = new TGeoMixture("ecal_mat", 3, density);
   ecal_mat->AddElement(elPb,1);
   ecal_mat->AddElement(elW,1);
   ecal_mat->AddElement(elO,4);
   TGeoMedium *ecal_med = new TGeoMedium("ecal_med",0,ecal_mat);
   
   for ( int j = 0; j < crystal_nz ; ++j ) {
      for ( int i = 0; i < crystal_nphi ; ++i ) {
         
         iptr = i+j*crystal_nphi;
         
         TGeoVolume *ecal = geom->MakeTubs(TString::Format("ecal-%d-%d",j,i),ecal_med,
                                           ecal_rmin, ecal_rmax,
                                           ecal_dz,
                                           ecal_sphi+i*ecal_dphi,
                                           ecal_sphi+(i+1)*ecal_dphi);
         ecal->SetLineColor(iptr);
         // top->AddNode(ecal,1,new TGeoCombiTrans(0,0,0,new TGeoRotation("ecal",0,0,0)));

         // GPThreeVector ecal_trans = GPThreeVector_create(0,0,ecal_zmin+(2.0*j+1.0)*ecal_dz);
         double dx = 0.0;
         double dy = 0.0;
         double dz = ecal_zmin+(2.0*j+1.0)*ecal_dz;
         TGeoTranslation *ecal_trans = new TGeoTranslation("",dx,dy,dz);

//         GPLogicalVolume_Constructor(ecal_log+iptr, (GPVSolid*)ecal, ecal_mat);
         // GPVPhysicalVolume_Constructor(ecal_phy+iptr, idRot, ecal_trans, ecal_log+iptr);
         world->AddNode(ecal,iptr,ecal_trans);
         
         //Set mother
//         GPVPhysicalVolume_SetMotherLogical(ecal_phy+iptr, world_log);
         
//         addLogicalVolumePointers( ecal_log+iptr);
//         addPhysicalVolumePointers( ecal_phy+iptr);
         
      }
   }
   
   //add daughter volume
//   for ( int j=0; j < crystal_nz ; ++j ) {
//      for ( int i=0; i < crystal_nphi ; ++i ) {
//         iptr = i+j*crystal_nphi;
//         addLogicalVolumeDaughter( world_log, ecal_phy+iptr);
//      }
//   }
   
   // Register world volume pointers for relocation
//   addLogicalVolumePointers( world_log );
//   addPhysicalVolumePointers( world_phy );

   geom->CloseGeometry();

   return world;
}





void runGPU(Int_t nthreads=10, Bool_t graphics=kFALSE, const char *geomfile="http://root.cern.ch/files/cms.root")
{
   gSystem->Load("libPhysics.so");
   gSystem->Load("libHist.so");
   gSystem->Load("libThread.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libGeantCuda.so");
   gSystem->Load("libGeant.so");
   
   GeantPropagator *prop = GeantPropagator::Instance();
   WorkloadManager *wmgr = WorkloadManager::Instance(nthreads);
   wmgr->SetNminThreshold(5*nthreads);
   CoprocessorBroker *gpuBroker = new CoprocessorBroker();
   gpuBroker->CudaSetup(32,128,1);
   wmgr->SetCoprocessorBroker(gpuBroker);
   
   prop->fNtotal   = 150; // 0;  // Number of events to be transported
   prop->fNevents  = 100;   // Number of buffered events
   prop->fNaverage = 50; // 0;   // Average number of tracks per event
   prop->fNperBasket = 100; // added zero
//   prop->fUseDebug = kTRUE;
//   prop->fDebugTrk = 1;
   
   // This sets gGeomManager and hence superseeds the filename.
   VP_SimpleECal();
   prop->PropagatorGeom("", nthreads, graphics);
   
   delete gGeoManager;
   
}   
