#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoShape.h"
#include "TMath.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoPgon.h"
#include "TGeoCompositeShape.h"
#include "TGeoSphere.h"
#include "TGeoXtru.h"

const Int_t NG = 33;
const Int_t NS = 15;
const char *exps[NG] = {"aleph", "barres", "felix", "phenix", "chambers", "p326",
                        "bes", "dubna", "ganil", "e907", "phobos2", "hermes", "na35",
                        "na47", "na49", "wa91", "sdc", "integral", "ams", "brahms",
                        "gem", "tesla", "btev", "cdf", "hades2", "lhcbfull", "star", 
                        "sld", "cms", "alice2", "babar2", "belle", "atlas" };
const char *shp[NS] = {"box", "tube", "tubeseg", "cone", "coneseg", "trd1", "trd2",
                       "para", "pcon", "pgon", "trap", "gtra", "sphere", "composite", "xtru"};
void perfshapes();
void countshapes();

Double_t performance(TGeoVolume* vol) {
   TGeoShape *shape = vol->GetShape();
   TStopwatch timer;
   Int_t ntr = 10000000;
   Double_t theta, phi;
   TGeoBBox *box = (TGeoBBox *)shape;
   Double_t dx = box->GetDX();
   Double_t dy = box->GetDY();
   Double_t dz = box->GetDZ();
   Double_t ox = (box->GetOrigin())[0];
   Double_t oy = (box->GetOrigin())[1];
   Double_t oz = (box->GetOrigin())[2];
   Double_t *vec = new Double_t[6*ntr];
   for (Int_t i=0; i<ntr; i++) {
      vec[6*i] = ox-dx+2*dx*gRandom->Rndm();
      vec[6*i+1] = oy-dy+2*dy*gRandom->Rndm();
      vec[6*i+2] = oz-dz+2*dz*gRandom->Rndm();
      phi = 2*TMath::Pi()*gRandom->Rndm();
      theta= TMath::ACos(1.-2.*gRandom->Rndm());
      vec[6*i+3]=TMath::Sin(theta)*TMath::Cos(phi);
      vec[6*i+4]=TMath::Sin(theta)*TMath::Sin(phi);
      vec[6*i+5]=TMath::Cos(theta);
   }
   timer.Start(kTRUE);
   for (Int_t i=0; i<ntr; i++) shape->DistFromOutside(&vec[6*i], &vec[6*i+3],3);
   timer.Stop();
//   timer.Print();
   Double_t time = timer.CpuTime()/ntr;
   printf("%s: Time per DistFromOutside: %f [ms]\n", shape->ClassName(), time*1.E6);
   delete [] vec;
   return time;
}

void countshapes()
{
   Double_t counters[NS] = {0};
   TString fname;
   Int_t nnodes;
   TGeoNode *node;
   TGeoShape *shape;
   for (Int_t i=0; i<NG; i++) {
      nnodes = 0;
      fname = TString::Format("%s.root", exps[i]);
      TGeoManager *geom = TGeoManager::Import(Form("http://root.cern.ch/files/%s",fname.Data()));
      TGeoIterator *next = new TGeoIterator(gGeoManager->GetTopVolume());
      while ((node=next->Next())) {
         nnodes++;
         shape = node->GetVolume()->GetShape();
         if (!strcmp(shape->ClassName(), "TGeoBBox")) counters[0] += 1;
         if (!strcmp(shape->ClassName(), "TGeoTube")) counters[1]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoTubeSeg")) counters[2]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoCone")) counters[3]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoConeSeg")) counters[4]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoTrd1")) counters[5]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoTrd2")) counters[6]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoPara")) counters[7]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoPcon")) counters[8]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoPgon")) counters[9]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoTrap")) counters[10]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoGtra")) counters[11]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoSphere")) counters[12]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoCompositeShape")) counters[13]+= 1;
         if (!strcmp(shape->ClassName(), "TGeoXtru")) counters[14]+= 1;
      }
      printf("geometry %s with %d nodes\n", exps[i], nnodes+1);
      delete next;
      delete geom;
   }
   
      
   Double_t sum = 0;
   for (Int_t i=0; i<NS; ++i) sum += counters[i];
   printf("=========  SHAPES  USAGE ==========\n");
   for (Int_t i=0; i<NS; ++i)  {counters[i] /= sum; printf("%s: %g %%\n", shp[i], 100.*counters[i]);}
   perfshapes();
}
   
void perfshapes()
{
   TGeoVolume *vol;
   TGeoManager *geom = new TGeoManager("test", "test");
   TGeoMaterial *mat = new TGeoMaterial("Al", 26.98,13,2.7);
   TGeoMedium *med = new TGeoMedium("MED",1,mat);
   Double_t t[NS] = {0.};
   vol = gGeoManager->MakeBox("BOX",med, 20,30,40);
   t[0] = performance(vol);
   vol = gGeoManager->MakeTube("TUBE",med, 20,30,40);
   t[1] = performance(vol);
   vol = gGeoManager->MakeTubs("TUBESEG",med, 20,30,40,-30,270);
   t[2] = performance(vol);
   vol = gGeoManager->MakeCone("CONE",med, 40,10,20,35,45);
   t[3] = performance(vol);
   vol = gGeoManager->MakeCons("CONESEG",med, 40,30,40,10,20,-30,250);
   t[4] = performance(vol);
   vol = gGeoManager->MakeTrd1("Trd1",med, 10,20,30,40);
   t[5] = performance(vol);
   vol = gGeoManager->MakeTrd2("Trd2",med, 10,20,30,10,40);
   t[6] = performance(vol);
   vol = gGeoManager->MakePara("PARA",med, 20,30,40,30,15,30);
   t[7] = performance(vol);
   vol = gGeoManager->MakePcon("PCON",med, -30.0,300,4);
   TGeoPcon *pcon = (TGeoPcon*)(vol->GetShape());
   pcon->DefineSection(0,0,15,20);
   pcon->DefineSection(1,20,15,20);
   pcon->DefineSection(2,20,15,25);
   pcon->DefineSection(3,50,15,20);
   t[8] = performance(vol);
   vol = gGeoManager->MakePgon("PGON",med, -45.0,270.0,4,4);
   TGeoPgon *pgon = (TGeoPgon*)(vol->GetShape());
   pgon->DefineSection(0,-70,45,50);
   pgon->DefineSection(1,0,35,40);
   pgon->DefineSection(2,0,30,35);
   pgon->DefineSection(3,70,90,100);
   t[9] = performance(vol);
   vol = gGeoManager->MakeTrap("Trap",med, 30,15,30,20,10,15,0,20,10,15,0);
   t[10] = performance(vol);
   vol = gGeoManager->MakeGtra("Gtra",med, 30,15,30,30,20,10,15,0,20,10,15,0);
   t[11] = performance(vol);
   vol = gGeoManager->MakeSphere("SPHERE",med, 30,40,60,120,30,240);
   t[12] = performance(vol);
   pgon = new TGeoPgon("pg",0.,360.,6,2); 
   pgon->DefineSection(0,0,0,20);
   pgon->DefineSection(1, 30,0,20);

   new TGeoSphere("sph", 40., 45.);
   // define named geometrical transformations with names
   TGeoTranslation *tr = new TGeoTranslation(0., 0., 45.);
   tr->SetName("tr");
   // register all used transformations
   tr->RegisterYourself();
   // create the composite shape based on a Boolean expression
   TGeoCompositeShape *cs = new TGeoCompositeShape("mir", "sph:tr*pg");

   vol = new TGeoVolume("COMP",cs);
   t[13] = performance(vol);
   vol = gGeoManager->MakeXtru("XTRU",med,4);
   TGeoXtru *xtru = (TGeoXtru*)vol->GetShape();
   Double_t x[8] = {-30,-30,30,30,15,15,-15,-15};
   Double_t y[8] = {-30,30,30,-30,-30,15,15,-30};
   xtru->DefinePolygon(8,x,y);
   xtru->DefineSection(0,-40, -20., 10., 1.5);
   xtru->DefineSection(1, 10, 0., 0., 0.5);
   xtru->DefineSection(2, 10, 0., 0., 0.7);
   xtru->DefineSection(3, 40, 10., 20., 0.9);
   t[14] = performance(vol);
   Double_t tot = 0.;
   printf("=========  SHAPES  BUDGET (DistFromOutside) ==========\n");
   for (Int_t i=0; i<NS; ++i) tot += t[i];
   for (Int_t i=0; i<NS; ++i) {t[i] = t[i]/tot; printf("%s relative budget=%g %%\n", shp[i], t[i]*100);}
   delete geom;
}
