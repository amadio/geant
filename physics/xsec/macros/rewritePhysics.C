#include <iostream>
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TPartIndex.h"
#include "TEXsec.h"
#include "TPDecay.h"
#include "TEFstate.h"

void rewritePhysics()
{
/*
  const char *fxsec = "../../data/xsec_FTFP_BERT_G496p02_1mev.root";
  const char *ffins = "../../data/fstate_FTFP_BERT_G496p02_1mev.root";
  const char *fxsecOut="xsec_FTFP_BERT_G496p02_1mev_phase1";
  const char *ffinsOut="fstate_FTFP_BERT_G496p02_1mev_phase1";
*/
  const char *fxsec="xsec_FTFP_BERT_G496p02_1mev_phase1";
  const char *ffins="fstate_FTFP_BERT_G496p02_1mev_phase1";
  const char *fxsecOut="xsec_FTFP_BERT_G496p02_1mev_phase2";
  const char *ffinsOut="fstate_FTFP_BERT_G496p02_1mev_phase2";

  TFile *xsec = new TFile(fxsec,"read");
  TFile *xsecOut = new TFile(fxsecOut,"recreate");
// Deal with xsec first
TIter next1(xsec->GetListOfKeys());
TKey *fkey=nullptr;
while (( fkey = (TKey*) next1())) {
   std::cout << fkey->GetName() << " class " << fkey->GetClassName() << std::endl;
   void *obj=nullptr;
   obj = xsec->GetObjectChecked(fkey->GetName(),fkey->GetClassName());
   if(!strcmp(fkey->GetClassName(),"TPartIndex"))
    xsecOut->WriteObject((TPartIndex*) obj,fkey->GetName());
    if(!strcmp(fkey->GetClassName(),"TEXsec"))
     xsecOut->WriteObject((TEXsec*) obj,fkey->GetName());
}

xsec->Close();
xsecOut->Close();

TFile *fins = new TFile(ffins,"read");
TFile *finsOut = new TFile(ffinsOut,"recreate");
TIter next2(fins->GetListOfKeys());

while (( fkey = (TKey*) next2())) {
   std::cout << fkey->GetName() << " class " << fkey->GetClassName() << std::endl;
   void *obj=nullptr;
   obj = fins->GetObjectChecked(fkey->GetName(),fkey->GetClassName());
   if(!strcmp(fkey->GetClassName(),"TPDecay"))
    finsOut->WriteObject((TPDecay*) obj,fkey->GetName());
  if(!strcmp(fkey->GetClassName(),"TEFstate"))
    finsOut->WriteObject((TEFstate*) obj,fkey->GetName());
}

 fins->Close();
 finsOut->Close();
}
