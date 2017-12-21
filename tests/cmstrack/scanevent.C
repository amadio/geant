#include <stdio.h>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TH1F.h>
#include <THashList.h>
#include <TTree.h>
#include <TRandom.h>
#include <TMath.h>

void scanevent(Int_t event=0, Bool_t verbose=false, double zlim=1e10)
{

   double fX;            // x position
   double fY;            // y position
   double fZ;            // z position
   double fPx;           // x momentum
   double fPy;           // y momentum
   double fPz;           // z momentum
   Short_t fPID;         // PDG particle id
   UShort_t fLVid;       // logical volume id
   UShort_t fShapeid;    // shape id 
   double fSafety;       // safety
   double fSnext;        // snext
   double fStep;         // step
   UChar_t fSurfid;      // surface id
   UChar_t fProcess;     // Process
   UChar_t fBegEnd;      // Beginning or end of track
   UInt_t  fTrid;        // Track ID
   UInt_t  fTrPid;       // Track Parend ID
   Double_t fCPUtime;    // CPU time used since start of track
   Double_t fCPUstep;    // CPU time used for current step

   TGeoManager::Import("cmstrack.root");

   TFile *f = new TFile("cmstrack.root","read");

   char evname[10];
   snprintf(evname,10,"Event%4.4d",event);
   printf("%s\n",evname);
   TTree *ev = (TTree*) f->Get(evname);
   ev->Print();

   ev->SetBranchAddress("x",&fX);
   ev->SetBranchAddress("y",&fY);
   ev->SetBranchAddress("z",&fZ);
   ev->SetBranchAddress("px",&fPx);
   ev->SetBranchAddress("py",&fPy);
   ev->SetBranchAddress("pz",&fPz);
   ev->SetBranchAddress("pid",&fPID);
   ev->SetBranchAddress("lvid",&fLVid);
   ev->SetBranchAddress("shapeid",&fShapeid);
   ev->SetBranchAddress("safety",&fSafety);
   ev->SetBranchAddress("snext",&fSnext);
   ev->SetBranchAddress("step",&fStep);
   ev->SetBranchAddress("surfid",&fSurfid);
   ev->SetBranchAddress("process",&fProcess);
   ev->SetBranchAddress("begend",&fBegEnd);
   ev->SetBranchAddress("trid",&fTrid);
   ev->SetBranchAddress("trpid",&fTrPid);
   ev->SetBranchAddress("cputime",&fCPUtime);
   ev->SetBranchAddress("cpustep",&fCPUstep);

   THashList *lv = (THashList *) f->Get("LogicalVolumes");
   // lv->Print();
   THashList *lp = (THashList *) f->Get("ProcessDictionary");
   // lp->Print();
   THashList *ls = (THashList *) f->Get("ShapeDictionary");
   // ls->Print();
   TDatabasePDG * pdg = TDatabasePDG::Instance();
   // pdg->Print();

   int nmiss = 0;
   TH1F *hstep = new TH1F("hstep","Difference of distance to boundary",100,0,0);
   int npoint = ev->GetEntries();
   double dperc = 2;
   int dpoint = 0.01*2*npoint+0.5;
   for(int ipoint=0; ipoint<npoint; ++ipoint) {
	
      if(!(ipoint%dpoint)) printf("Now at %3d%%\n",(int)(100.*ipoint/npoint+0.5));
	 
      ev->GetEntry(ipoint);
      TParticlePDG *particle = pdg->GetParticle(fPID);
      const char *gvol = lv->At(fLVid)->GetName();

      double point[3];
      double dir[3];
      double norm = sqrt(fPx*fPx+fPy*fPy+fPz*fPz);
      if(fSafety > 5e-10 && TMath::Abs(fZ)<zlim) {
	 if(norm>1e-10) {
	    norm=1./norm;
	    dir[0]=fPx*norm;dir[1]=fPy*norm;dir[2]=fPz*norm;
	 } else {
	    const double phi = gRandom->Rndm()*TMath::TwoPi();
	    double costh = 1-2*gRandom->Rndm();
	    double sinth = TMath::Sqrt((1+costh)*(1-costh));
	    dir[0] = sinth*TMath::Cos(phi);
	    dir[1] = sinth*TMath::Sin(phi);
	    dir[2] = costh;
	 }
	 point[0]=fX*0.1;point[1]=fY*0.1;point[2]=fZ*0.1;
	 const char* rvol = gGeoManager->InitTrack(point,dir)->GetVolume()->GetName();

	 char rvols [1024];
	 int lrvol = strlen(rvol);
	 //
	 // problem: sometimes root adds an hex number at the end of the name
	 // if there is a hex number at the end of the file, we create a new
	 // name without the last four chars
	 // if the root name does not match either G4 live or G4 cold, we 
	 // check the root name without the four last hex chars
	 //
	 strcpy(rvols,rvol);
	 bool numbers=lrvol>4;
	 if(numbers)
	    for(int i=lrvol-4; i<lrvol; ++i) 
	       numbers &= (('0'<=rvol[i] && rvol[i]<='9') || ('a'<=rvol[i] && rvol[i]<='f'));
	 if(numbers) rvols[strlen(rvol)-4]='\0';


	 //	 printf("Now in %s\n",rvol);
	 if(strcmp(rvol,gvol)&&strcmp(rvols,gvol)) {
	    ++nmiss;
	    printf("ROOT vol %s != Geant4 vol %s sf %9.3g sn %9.3g st %9.3g mmiss ppm %7.1f\n",rvol,gvol,fSafety, fSnext, fStep, 1000000.*nmiss/npoint);
	    const char *rpath = gGeoManager->GetPath();
	    printf("ROOT    path %s\n",rpath);
	    char rpoint[2048];
	    rpoint[0]='\0';
	    strcat(rpoint,Form("/%8.03g,%8.03g,%8.03g",point[0],point[1],point[2]));
	    for ( Int_t i=0; i<gGeoManager->GetLevel(); ++i) {
	       Double_t plocal[3];
	       gGeoManager->GetMotherMatrix(gGeoManager->GetLevel()-1-i)->MasterToLocal(point,plocal);
	       strcat(rpoint,Form("/%8.03g,%8.03g,%8.03g",plocal[0],plocal[1],plocal[2]));
	    }
	    printf("%s\n",rpoint);

	 } else if(norm>1e-10) {
	    //	    printf("%f %f %f %f\n",fPx,fPy,fPz,norm);
	    gGeoManager->FindNextBoundaryAndStep();
	    hstep->Fill(gGeoManager->GetStep()*10-fStep);
	 }
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
	 printf("%-9.9s",lp->At(fProcess)->GetName());
	 printf(" sf %9.3g sn %9.3g st %9.3g",fSafety, fSnext, fStep);
	 printf(" lvl %s(%d) (%s)",gvol,fLVid,
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
   printf("Mismatch %10.3g%%\n",(100.*nmiss)/ev->GetEntries());
   hstep->Draw();
}
