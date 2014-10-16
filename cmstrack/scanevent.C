void scanevent(Int_t event=0, Bool_t verbose=kFALSE)
{

   TGeoManager::Import("cms.gdml");
   gGeoManager->Export("cms.root");

   double fX;            // x position
   double fY;            // y position
   double fZ;            // z position
   double fPx;           // x momentum
   double fPy;           // y momentum
   double fPz;           // z momentum
   Short_t fPID;         // PDG particle id
   UShort_t  fLVid;      // logical volume id
   double fSafety;       // safety
   double fSnext;        // snext
   double fStep;         // step
   UChar_t fSurfid;      // surface id
   UChar_t fProcess;     // Process
   UChar_t fBegEnd;      // Beginning or end of track

   TFile *f = new TFile("cmstrack.root","read");
   char evname[10];
   snprintf(evname,10,"Event%4.4d",event);
   printf("%s\n",evname);
   TTree *ev = (TTree*) f->Get(evname);
   // ev->Print();

   ev->SetBranchAddress("x",&fX);
   ev->SetBranchAddress("y",&fY);
   ev->SetBranchAddress("z",&fZ);
   ev->SetBranchAddress("px",&fPx);
   ev->SetBranchAddress("py",&fPy);
   ev->SetBranchAddress("pz",&fPz);
   ev->SetBranchAddress("pid",&fPID);
   ev->SetBranchAddress("lvid",&fLVid);
   ev->SetBranchAddress("safety",&fSafety);
   ev->SetBranchAddress("snext",&fSnext);
   ev->SetBranchAddress("step",&fStep);
   ev->SetBranchAddress("surfid",&fSurfid);
   ev->SetBranchAddress("process",&fProcess);
   ev->SetBranchAddress("begend",&fBegEnd);

   TObjArray lv;
   lv.Read("LogicalVolumes");
   //   lv.Print();
   TObjArray lp;
   lp.Read("ProcessDictionary");
   // lp.Print();
   TDatabasePDG * pdg = TDatabasePDG::Instance();
   // pdg->Print();

   int nmiss = 0;
   TH1F *hstep = new TH1F("hstep","Difference of distance to boundary",100,0,0);
   for(int ipoint=0; ipoint<ev->GetEntries(); ++ipoint) {
      ev->GetEntry(ipoint);
      TParticlePDG *particle = pdg->GetParticle(fPID);
      const char* gvol = ((TObjString*)lv->At(fLVid))->GetString()->Data();

      double point[3];
      double dir[3];
      double norm = sqrt(1./(fPx*fPx+fPy*fPy+fPz*fPz));
      point[0]=fX*0.1;point[1]=fY*0.1;point[2]=fZ*0.1;
      dir[0]=fPx*norm;dir[1]=fPy*norm;dir[2]=fPz*norm;
      const char* rvol = gGeoManager->InitTrack(point,dir)->GetVolume()->GetName();
      if(strcmp(rvol,gvol)) {
	 printf("ROOT vol %s != Geant4 vol %s sf %9.3g sn %9.3g st %9.3g\n",rvol,gvol,fSafety, fSnext, fStep);
	 nmiss++;
      } else if(strcmp(rvol,"SIC1")) {
	 gGeoManager->FindNextBoundaryAndStep();
	 hstep->Fill(gGeoManager->GetStep()/10-fStep);
      }

      if(verbose) {
	 if(fBegEnd==1 || fBegEnd==3) {
	    printf("========================================================= start ");
	    if(particle) printf("%-8.8s ",particle->GetName());
	    else printf("Unknown   ");
	    printf(" =========================================================\n");
	 }
	 printf("x(%9.3g,%9.3g,%9.3g) ", fX, fY, fZ);
	 printf("p(%9.3g,%9.3g,%9.3g) ", fPx, fPy, fPz);
	 printf("%-9.9s",((TObjString*)lp->At(fProcess))->GetString()->Data());
	 printf(" sf %9.3g sn %9.3g st %9.3g",fSafety, fSnext, fStep);
	 printf(" lvl %s (%s)",gvol,
		&(gGeoManager->GetVolume(gvol)->GetShape()->ClassName())[4]);
	 printf("\n");
	 if(fBegEnd==2 || fBegEnd==3) {
	    printf("========================================================= end   ");
	    if(particle) printf("%-8.8s ",particle->GetName());
	    else printf("Unknown   ");
	    printf(" =========================================================\n");
	 }
      }
   }
   printf("Mismatch %10.3g\%\n",(100.*nmiss)/ev->GetEntries());
   hstep->Draw();
}
