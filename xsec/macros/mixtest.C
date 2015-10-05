#include "TSystem.h"
#include "TFile.h"

#include "TPartIndex.h"
#include "TMXsec.h"
#include "TPDecay.h"

#include "iostream"

void mixtest() {
   const char *fxsec = "../../data/xsec_FTFP_BERT_G496p02_1mev.root";
   const char *ffins = "../../data/fstate_FTFP_BERT_G496p02_1mev.root";
   gSystem->Load("libXsec");
   TFile *ff = new TFile(ffins,"r");
   TPDecay *dt = (TPDecay*) ff->Get("DecayTable");
   TFile *fx = new TFile(fxsec,"r");

   int z[2]={1,16};
   float w[2]={2,1};
   int a[2];
   TMXsec *tm = nullptr;
   tm = new TMXsec("h2o","water",z,a,w,2,1,kFALSE,dt);
   tm->Print();
   w[0] = 2*TPartIndex::I()->WEle(1)/(2*TPartIndex::I()->WEle(1)+TPartIndex::I()->WEle(16));
   w[1] = TPartIndex::I()->WEle(16)/(2*TPartIndex::I()->WEle(1)+TPartIndex::I()->WEle(16));
   tm = new TMXsec("h2o","water",z,a,w,2,1,kTRUE,dt);
   tm->Print();
}
