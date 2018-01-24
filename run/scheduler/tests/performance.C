#include "TMath.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoShape.h"
#include "TGeoBBox.h"

void performance() {
   TStopwatch timer;
   Int_t ntr = 10000000;
   Double_t theta, phi;
   TGeoVolume *vol = gGeoManager->GetTopVolume()->GetNode(0)->GetVolume();
   const TGeoShape *shape = vol->GetShape();
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
   timer.Start(true);
   for (Int_t i=0; i<ntr; i++) shape->DistFromOutside(&vec[6*i], &vec[6*i+3],3);
   timer.Stop();
   timer.Print();
   Double_t time = timer.CpuTime()/ntr;
   printf("Time per DistFromOutside: %f [ms]\n", time*1.E6);
   delete [] vec;
}
