#include <TGeoExtension.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TList.h>

#include <TMXsec.h>

void GenerateEvent();

static TList particleStack;

Int_t main (int argc, char *argv[]) {

   for(Int_t i=0; i<argc; ++i) {
      printf("argv[%d] = %s\n",i,argv[i]);
   }

   const Char_t *geofile="http://root.cern.ch/files/cms.root";
   TGeoManager *geom = TGeoManager::Import(geofile);
   
   // loop materials

   TList *matlist = (TList*) geom->GetListOfMaterials();
   TIter next(matlist);
   TGeoMaterial *mat=0;
   TGeoMixture *mix=0;
   Int_t nmater = matlist->GetEntries();
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
	 printf("Mixture %s element %s z %d a %d\n",
		mat->GetName(), mat->GetElement(iel)->GetName(),
		z[iel],a[iel]);
      }
      mat->SetFWExtension(
	new TGeoRCExtension(
	   new TMXsec(mat->GetName(),mat->GetTitle(),
		      z,a,w,nelem,mat->GetDensity(),kTRUE)));
      //      myObject = mat->GetExtension()->GetUserObject();
      delete [] a;
      delete [] z;
      delete [] w;
   }

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
	 current = geom->InitTrack(pos,dir);
	 Double_t ken = 
	 while(!geom->IsOutside()) {
	    mat = current->GetVolume()->GetMaterial();
	    Double_t xlen = mat->GetUserField()->Xlength(G5index,;
	    nexnode = geom->FindNextBoundaryAndStep(xlen);
	    Double_t snext = geom->GetStep();
	    if(snext>xlen) {
	       //phys wins
	    } else {
	       // geom wins
	       dirnew = something;
	       geom->SetCurrentDirection(dirnew);
		  
	    }
	 }
	 
	    
      }
   */
   return 0;
}
/*
#define NPART 11

void GenerateEvent(Double_t avemult) {
   static Bool_t first=kTRUE;
   const npart=NPART;
   static const Char_t* partnam[NPART] = {"pi+","pi-","proton","antiproton","neutron","antineutron","e-","e+","gamma"
				    "mu+","mu-"};
   static Int_t G5part[NPART];
   static Float_t G5prob[NPART] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

   const Double_t etamin = -3, etamax = 3;

   // Initialise simple generator
   if(first) {
      Double_t sumprob=0;
      for(Int_t ip=0; ip<npart; ++ip) {
	 G5part[ip] = TPartIndex::I()->PartIndex(partnam[ip]);
	 sumprob += G5prob[ip];
      }
      for(Int_t ip=0; ip<npart; ++ip) {
	 G5prob[ip]/=sumprob;
	 if(ip) G5Prob[i]+=G5prob[ip-1];
      }
   }
   
   Int_t ntracks = ntracks = td->fRndm->Poisson(average);
   for (Int_t i=0; i<ntracks; i++) {
      Double_t prob = gRandom->Uniform();
      for(Int_t j=0; j<kMaxPart; ++j) {
	 if(prob <= pdgProb[j]) {
	    track->pdg = pdgGen[j];
	    track->species = pdgSpec[j];
	    //            Printf("Generating a %s",TDatabasePDG::Instance()->GetParticle(track->pdg)->GetName());
	    pdgCount[j]++;
	    break;
	 }
      }   
         if(!track->pdg) Fatal("ImportTracks","No particle generated!");
         TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(track->pdg);
         track->charge = part->Charge()/3.;
         track->mass   = part->Mass();
         track->xpos = fVertex[0];
         track->ypos = fVertex[1];
         track->zpos = fVertex[2];
         track->e = fKineTF1->GetRandom()+track->mass;
         Double_t p = TMath::Sqrt((track->e-track->mass)*(track->e+track->mass));
         Double_t eta = td->fRndm->Uniform(etamin,etamax);  //multiplicity is flat in rapidity
         Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
         //Double_t theta = TMath::ACos((1.-2.*gRandom->Rndm()));
         Double_t phi = TMath::TwoPi()*td->fRndm->Rndm();
         track->px = p*TMath::Sin(theta)*TMath::Cos(phi);
         track->py = p*TMath::Sin(theta)*TMath::Sin(phi);
         track->pz = p*TMath::Cos(theta);
         track->frombdr = kFALSE;
         Int_t itrack = track->particle;
	 
         fCollections[tid]->AddTrack(itrack, basket);
//         gPropagator->fCollections[tid]->AddTrack(itrack, basket);
    //     basket->AddTrack(fNstart);
         fNstart++;
      }
//      Printf("Event #%d: Generated species for %6d particles:", event, ntracks);
      event++;
      for (Int_t i=0; i<kMaxPart; i++) {
//         Printf("%15s : %6d particles", TDatabasePDG::Instance()->GetParticle(pdgGen[i])->GetName(), pdgCount[i]);
         pdgCount[i] = 0;
      }   
   }
//   Printf("Injecting %d events...", nevents);
   InjectCollection(tid);      
   return basket;
}
   */
