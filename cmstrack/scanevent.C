#include <stdio.h>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TH1F.h>
#include <THashList.h>
#include <TTree.h>

void scanevent(Int_t event=0, Bool_t verbose=kFALSE)
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
   for(int ipoint=0; ipoint<ev->GetEntries(); ++ipoint) {
      ev->GetEntry(ipoint);
      TParticlePDG *particle = pdg->GetParticle(fPID);
      const char *gvol = lv->At(fLVid)->GetName();

      double point[3];
      double dir[3];
      double norm = sqrt(fPx*fPx+fPy*fPy+fPz*fPz);
      if(norm>1e-10 && fSafety > 5e-10) {
	 norm=1./norm;
	 point[0]=fX*0.1;point[1]=fY*0.1;point[2]=fZ*0.1;
	 dir[0]=fPx*norm;dir[1]=fPy*norm;dir[2]=fPz*norm;
	 const char* rvol = gGeoManager->InitTrack(point,dir)->GetVolume()->GetName();
	 //	 printf("Now in %s\n",rvol);
	 if(strcmp(rvol,gvol)) {
	    printf("ROOT vol %s != Geant4 vol %s sf %9.3g sn %9.3g st %9.3g\n",rvol,gvol,fSafety, fSnext, fStep);
	    nmiss++;
	 } else {
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
