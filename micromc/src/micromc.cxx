
#include <TFile.h>
#include <TGeoExtension.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TList.h>
#include <TRandom.h>
#include <TGeoBBox.h>
#include <GeantTrack.h>

#include <TPartIndex.h>
#include <TMXsec.h>
#include <TPXsec.h>

void GenerateEvent(Double_t avemult, Double_t energy, Double_t fVertex[3]);
Double_t SampleMaxwell(Double_t emean);
void IncreaseStack();
void VertexIn(TGeoBBox *bbox, Double_t ori[3]);

static Int_t stacksize=100;
static Int_t hwmark=0;
static GeantTrack *particleStack=new GeantTrack[stacksize];
static TGeoManager *geom;

Int_t main (int argc, char *argv[]) {

   for(Int_t i=0; i<argc; ++i) {
      printf("argv[%d] = %s\n",i,argv[i]);
   }
/*
   Int_t nevent=1;
   if(argc>1) sscanf(argv[1],"%d",&nevent);

   Double_t avemult = 10.;
   if(argc>2) sscanf(argv[2],"%lf",&avemult);
   
   Double_t energy = 10.;
   if(argc>3) sscanf(argv[3],"%lf",&energy);
		 
   printf("Generating %d events with ave multiplicity %f and energy %f\n",nevent,avemult, energy);

   const Char_t *geofile="http://root.cern.ch/files/cms.root";
   geom = TGeoManager::Import(geofile);

   // loop materials
   TFile *f = new TFile("xsec.root");
   f->Get("PartIndex");

   TPXsec::SetVerbose(1);
   TPartIndex::I()->SetEnergyGrid(1e-3,1e3,100);
   TList *matlist = (TList*) geom->GetListOfMaterials();
   TIter next(matlist);
   TGeoMaterial *mat=0;

   TList *matXS = new TList();
   matXS->Add(TPartIndex::I());
   TMXsec *mxs=0;
   printf("Total of %d materials\n",matlist->GetSize());
   while((mat = (TGeoMaterial*) next())) {
      if(!mat->IsUsed()) continue;
      Int_t nelem = mat->GetNelements();
      Int_t *z = new Int_t[nelem];
      Int_t *a = new Int_t[nelem];
      Float_t *w = new Float_t[nelem];
      for(Int_t iel=0; iel<nelem; ++iel) {
	 Double_t ad;
	 Double_t zd;
	 Double_t wd;
	 mat->GetElementProp(ad,zd,wd,iel);
	 a[iel]=ad;
	 z[iel]=zd;
	 w[iel]=wd;
	 //	 printf("Mixture %s element %s z %d a %d\n",
	 //	mat->GetName(), mat->GetElement(iel)->GetName(),
	 //	z[iel],a[iel]);
      }
      mxs = new TMXsec(mat->GetName(),mat->GetTitle(),
		       z,a,w,nelem,mat->GetDensity(),kTRUE);
      matXS->Add(mxs);
      mat->SetFWExtension(new TGeoRCExtension(mxs));
      //      myObject = mat->GetExtension()->GetUserObject();
      delete [] a;
      delete [] z;
      delete [] w;
   }

   TMXsec::Prune();

   TFile *fmxs = new TFile("mxs.root","recreate");
   fmxs->SetCompressionLevel(0);
   matXS->Write();
   fmxs->Close();
   
   TGeoVolume *top = geom->GetTopVolume();
   TGeoShape *shape = top->GetShape();
   TGeoBBox *bbox = (TGeoBBox*) shape;
   Double_t dx = bbox->GetDX();
   Double_t dy = bbox->GetDY();
   Double_t dz = bbox->GetDZ();
   const Double_t *origin = bbox->GetOrigin();
   printf("Top volume is %s shape %s\n",top->GetName(),shape->GetName());
   printf("BBox dx %f dy %f dz %f origin %f %f %f\n",dx,dy,dz,origin[0],origin[1],origin[2]);

   for(Int_t iev=0; iev<nevent; ++iev) {
      // should define a vertex, origin for the moment
      Double_t vertex[3]={0,0,0};
      VertexIn(bbox,vertex);
      GenerateEvent(avemult, energy, vertex);
      while(hwmark) {
	 GeantTrack *track = &particleStack[--hwmark];
	 printf("Transporting particle #%d %s energy %g\n",hwmark+1,
		TDatabasePDG::Instance()->GetParticle(track->fPDG)->GetName(),
		track->fE-track->fMass);
	 Double_t pos[3] = {track->fXpos,track->fYpos,track->fZpos};
	 Double_t dir[3] = {track->fXdir,track->fYdir,track->fZdir};
	 TGeoNode *current = geom->InitTrack(pos,dir);
	 Double_t pintl = -TMath::Log(gRandom->Rndm());
	 printf("Initial point %f %f %f in %s\n",pos[0],pos[1],pos[2],current->GetName());
	 while(!geom->IsOutside()) {
	    mat = current->GetVolume()->GetMaterial();
	    const Double_t *cpos = geom->GetCurrentPoint();
	    printf("Point now %f %f %f in %s made of %s\n",cpos[0],cpos[1],cpos[2],
		   current->GetName(),current->GetVolume()->GetMaterial()->GetName());
	    Double_t ken = track->fE-track->fMass;
	    TMXsec *mx = ((TMXsec *)
			  ((TGeoRCExtension*) 
			   mat->GetFWExtension())->GetUserObject());
	    Double_t xlen = mx->Xlength(track->fG5code,ken,TMath::Sqrt( (track->fE+track->fMass)*(track->fE-track->fMass) )  );
	    Double_t pnext = pintl*xlen;
	    current = geom->FindNextBoundaryAndStep(pnext);
	    Double_t snext = geom->GetStep();
	    printf("pnext = %f snext = %f\n",pnext,snext);
	    if(pnext<=snext) {
	       //phys wins
	       Int_t reac;
	       mat->Print();
	       TEXsec *el = mx->SampleInt(track->fG5code,ken,reac);
	       printf("particle does a %s on %s\n",TPartIndex::I()->ProcName(reac),el->GetName());
	       break;
	    } else {
	       // geom wins
	       pintl-=snext/xlen;
	       //	       printf("Geom wins\n");
	    }
	    //	    dirnew = something;
	    //	       geom->SetCurrentDirection(dirnew);
	 }
	 if(geom->IsOutside()) printf("Particle exited setup\n");
      }
   }
*/   
   /*
   Double_t dir[3];
   Double_t pos[3];
   TGeoNode *current=0;
   TGeoNode *nexnode=0;
   TIter next(particleStack);
   for(Int_t iev=0; iev<nevent; ++iev) {
      GenerateEvent();
      next.Reset();
      GeantTrack *tr=0;
      while((tr=(GeantTrack*)next())) {
	 Int_t G5index = TPartIndex::I()->PartIndex(tr->pdg);
	 tr->Direction(dir);
	 x[0]=tr->xpos;
	 x[1]=tr->ypos;
	 x[2]=tr->zpos;
	 // where am I
	    
      }
   */
   return 0;
}

#define NPART 11

void GenerateEvent(Double_t avemult, Double_t energy, Double_t fVertex[3]) {
   static Bool_t first=kTRUE;
   static const Int_t kMaxPart=NPART;
   static const Char_t* G5name[NPART] = {"pi+","pi-","proton","antiproton","neutron","antineutron","e-","e+",
					 "gamma", "mu+","mu-"};
   static const Species_t G5species[NPART] = {kHadron, kHadron, kHadron, kHadron, kHadron, kHadron, 
					      kLepton, kLepton, kLepton, kLepton, kLepton};
   static Int_t G5part[NPART];
   static Float_t G5prob[NPART] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

   const Double_t etamin = -3, etamax = 3;

   // Initialise simple generator
   if(first) {
      Double_t sumprob=0;
      for(Int_t ip=0; ip<kMaxPart; ++ip) {
	 G5part[ip] = TPartIndex::I()->PartIndex(G5name[ip]);
	 printf("part %s code %d\n",G5name[ip],G5part[ip]);
	 sumprob += G5prob[ip];
      }
      for(Int_t ip=0; ip<kMaxPart; ++ip) {
	 G5prob[ip]/=sumprob;
	 if(ip) G5prob[ip]+=G5prob[ip-1];
      }
      first=kFALSE;
   }
   
   Int_t ntracks = gRandom->Poisson(avemult)+0.5;
   
   hwmark=0;
   for (Int_t i=0; i<ntracks; i++) {
      if(hwmark==stacksize) IncreaseStack();
      GeantTrack *track=&particleStack[hwmark++];
      Double_t prob = gRandom->Uniform();
      for(Int_t j=0; j<kMaxPart; ++j) {
	 if(prob <= G5prob[j]) {
	    track->fG5code = G5part[j];
	    track->fPDG = TPartIndex::I()->PDG(G5part[j]);
	    track->fSpecies = G5species[j];
	    printf("Generating a %s\n",TDatabasePDG::Instance()->GetParticle(track->fPDG)->GetName());
	    //	    pdgCount[j]++;
	    break;
	 }
      }   
      if(!track->fPDG) Fatal("ImportTracks","No particle generated!");
      TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(track->fPDG);
      track->fCharge = part->Charge()/3.;
      track->fMass   = part->Mass();
      track->fXpos = fVertex[0];
      track->fYpos = fVertex[1];
      track->fZpos = fVertex[2];
      Double_t ekin = SampleMaxwell(energy/avemult);
      track->fE = ekin+track->fMass;
      Double_t p = TMath::Sqrt(ekin*(2*ekin+track->fMass));
      Double_t eta = gRandom->Uniform(etamin,etamax);  //multiplicity is flat in rapidity
      Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
      //Double_t theta = TMath::ACos((1.-2.*gRandom->Rndm()));
      Double_t phi = TMath::TwoPi()*gRandom->Rndm();
      track->fP = p;
      track->fXdir = TMath::Sin(theta)*TMath::Cos(phi);
      track->fYdir = TMath::Sin(theta)*TMath::Sin(phi);
      track->fZdir = TMath::Cos(theta);
      track->fFrombdr = kFALSE;
   }
//      Printf("Event #%d: Generated species for %6d particles:", event, ntracks);
}

Double_t SampleMaxwell(Double_t emean) 
{
   Double_t th = gRandom->Uniform()*TMath::TwoPi();
   Double_t rho = TMath::Sqrt(-TMath::Log(gRandom->Uniform()));
   Double_t mx = rho*TMath::Sin(th);
   return 2*emean*(-TMath::Log(gRandom->Uniform())+mx*mx)/3.;
}

void IncreaseStack() {
   Int_t newstacksize = stacksize*1.5;
   GeantTrack *tmp = new GeantTrack[newstacksize];
   memcpy((void*)tmp,(void*)particleStack,stacksize*sizeof(GeantTrack));
   delete [] particleStack;
   particleStack=tmp;
   stacksize = newstacksize;
}

void VertexIn(TGeoBBox *bbox, Double_t ori[3])
{
   Double_t eta=0;
   do {
      eta = gRandom->Rndm();
      ori[0] = bbox->GetDX()*eta*eta*(1?gRandom->Rndm()>0.5:-1);
      eta = gRandom->Rndm();
      ori[1] = bbox->GetDY()*eta*eta*(1?gRandom->Rndm()>0.5:-1);
      eta = gRandom->Rndm();
      ori[2] = bbox->GetDZ()*eta*eta*(1?gRandom->Rndm()>0.5:-1);
      geom->SetCurrentPoint(ori);
   } while(geom->IsOutside());
}
