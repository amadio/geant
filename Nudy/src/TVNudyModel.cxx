#include "TNudyCore.h"
#include "TVNudyModel.h"
//Header files required to display data
#include "TMath.h"
#include "TFile.h"
#include "TRandom3.h"
#include "Math/SpecFuncMathMore.h"
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TGraph2D.h>
#include <TStyle.h>

ClassImp(TVNudyModel)


//______________________________________________________________________________
TVNudyModel::TVNudyModel(){
  fMAT = 0;
  fTemp = 0;
  fReaction = (Reaction_t)0;
  fMaterial = NULL;
  fProjectile = NULL;
  fE_file3 = NULL;
  fXSect_file3 = NULL;
  fEXSect_length = 0;
  f4nens = 0;  
  fAPAlias = NULL;
  fEPtable = NULL;
  f5Tein = -1;
  f5Tel = 0;
  fPerc = NULL;
  fEPAlias = NULL;
  maxpop = 50;
  nens = 0;
  nperc = 25;
}

//_______________________________________________________________________________
TVNudyModel::TVNudyModel(TGeoElementRN *mat, Reaction_t reac, ULong_t temp,TParticlePDG* projectile,TNudyEndfMat *material){
  SetName(TNudyCore::Instance()->GetKey(mat,reac,temp));
  fEndf = mat->ENDFCode();
  fPdg = projectile->PdgCode();
  fMaterial = mat;
  fReaction = (Reaction_t)reac;
  fTemp = temp;    
  fProjectile = projectile;
  // Make TGeoElementRN -> MAT map?
  fMAT = 0;
  fE_file3 = NULL;
  fXSect_file3 = NULL;
  fEXSect_length = 0;
  fEPtable = NULL;
  fPerc = NULL;
  fEPAlias = NULL;
  maxpop = 50;
  f5Tein = -1;
  f5Tel = 0;
  nens = 0;
  nperc = 25;
  if(material){
    fMAT = material->GetMAT();
    ReadFile(material);
  }
}

//______________________________________________________________________________
TVNudyModel::~TVNudyModel(){
  if(fEPAlias){
    delete [] fEPAlias;
    fEPAlias = 0;
  }
  if(fAPAlias){
    delete [] fAPAlias;
    fAPAlias = 0;
  }
  if(fE_file3){
    delete [] fE_file3;
    fE_file3 = 0;
  }
  if(fXSect_file3){
    delete [] fXSect_file3;
    fXSect_file3 = 0;
  }
  if(fPerc){
    delete [] fPerc;
    fPerc = 0;
  }
}

//______________________________________________________________________________
Double_t TVNudyModel::GetEo(Double_t ein){
  Int_t el,eh;
  if(f5Tein != ein){
    Info("EIN","%e --- %e",f5Tein, ein);
   f5Tein = ein; 
   if(ein <= xengr[0]){
     Warning("GetEo","Incident energy %e below minimum threshold %e",ein,xengr[0]);
      el = 0;
    }else if(ein >= xengr[xengr.GetSize()-1]){
      Warning("GetEo","Incident energy %e above maximum threshold %e",ein,xengr[xengr.GetSize()-1]);
      el = xengr.GetSize()-1;
    }else{
      el = TNudyCore::Instance()->BinarySearch(xengr.GetArray(),nens,ein);
      Info("Searching","%d",el);
    }
   Info("","%d - > %d",f5Tel,el);
    f5Tel = el;
  }
  else{
    el = f5Tel;
  }
/*  eh = el+1;
  Double_t pIndex = fRnd.Rndm();
  Int_t pl = (pIndex * 25.0);
  Int_t ph = pl + 1;
  return (TNudyCore::Instance()->BilinearInterploation(xengr[el],pl/25.0,xengr[eh],ph/25.0,
						       fPerc[el].GetAt(pl),fPerc[el].GetAt(ph),fPerc[eh].GetAt(pl),fPerc[eh].GetAt(ph),
						       ein,pIndex));
*/
  return fEPAlias[el].ImprovedInterpolation(f5Tein);
}
//_____________________________________________________________________________
Double_t TVNudyModel::GetAo(Double_t ein){
   Int_t el,eh;
  if(f4Tein != ein){
    Info("EIN","%e --- %e",f5Tein, ein);
   f4Tein = ein; 
   if(ein <= fAPAlias[0].GetAlpha()){
     Warning("GetAo","Incident energy %e below minimum threshold %e",ein,fAPAlias[0].GetAlpha());
      el = 0;
    }else if(ein >= fAPAlias[f4nens-1].GetAlpha()){
      Warning("GetAo","Incident energy %e above maximum threshold %e",ein,fAPAlias[f4nens-1].GetAlpha());
      el = f4nens-2;
    }else{
      el = TNudyCore::Instance()->BinarySearch(f4eins.GetArray(),f4nens,ein);
      Info("Searching","%d",el);
    }
    f4Tel = el;
  }
  else{
    el = f4Tel;
  }
  return fAPAlias[el].ImprovedInterpolation(f4Tein);
}
//______________________________________________________________________________
TGeoElementRN* TVNudyModel::GetMaterial(){
  if(!fMaterial) fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fMaterial;
}

//______________________________________________________________________________
const char* TVNudyModel::GetMaterialName(){
  if(!fMaterial) fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fMaterial->GetName();
}

//______________________________________________________________________________
Int_t TVNudyModel::GetZ(){
  if(!fMaterial) fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fMaterial->Z();
}

//______________________________________________________________________________
Int_t TVNudyModel::GetA(){
  if(!fMaterial) fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return (Int_t)fMaterial->A();
}

//______________________________________________________________________________
Int_t TVNudyModel::GetISO(){
  if(!fMaterial) fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fMaterial->IsoNo();
}

//______________________________________________________________________________
Reaction_t TVNudyModel::GetReaction(){
  if(!fMaterial) fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fReaction;
}

//______________________________________________________________________________
ULong_t TVNudyModel::GetTemp(){
  if(!fMaterial) fMaterial = TNudyCore::Instance()->GetMaterial(fEndf);
  return fTemp;
}

//______________________________________________________________________________
void TVNudyModel::ReadFile(TNudyEndfMat *material){
  //Function to Read and Process File in ENDF Tape
  if(material->GetMAT() != fMAT){
    printf("Wrong Material");
    return;
  }
  TIter iter(material->GetFiles());
  TNudyEndfFile *file;
  while((file = (TNudyEndfFile*)iter.Next())){
    //Read File data into class structures
    SetTitle(Form("%s:%d",GetTitle(),file->GetMF()));
    switch(file->GetMF()){
    case 3:
      ReadFile3(file);
      break;
    case 4:
      ReadFile4(file);
      break;
    case 5:
      SetTitle(Form("%s(",GetTitle()));
      ReadFile5(file);
      SetTitle(Form("%s)",GetTitle()));
      break;
    }
  }
}

//______________________________________________________________________________
void TVNudyModel::ReadFile5(TNudyEndfFile *file){
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  
  //nens = 0;
  //maxpop = 50;
  while((sec = (TNudyEndfSec*)secIter.Next())){
    if(sec->GetMT() == (Int_t)fReaction){
      File5_Pass1(sec);
      File5_Pass2(sec);
      File5_Pass3();
    }
  }
} 
//______________________________________________________________________________
void TVNudyModel::ReadFile4(TNudyEndfFile *file){
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while((sec = (TNudyEndfSec*)secIter.Next())){
    if(sec->GetMT() == (Int_t)fReaction){
      TIter recIter(sec->GetRecords());
      TNudyEndfCont* header = (TNudyEndfCont*)recIter.Next();
      Int_t LTT = sec->GetL2();
      Int_t LI = header->GetL1();
      printf("LTT = %d LI = %d\n",LTT,LI);
      if(LTT == 1 && LI == 0){
        TNudyEndfTab2* subheader = (TNudyEndfTab2*)recIter.Next();
        TArrayD ein(subheader->GetN2());
        TArrayD *lCoef = new TArrayD[subheader->GetN2()]; 
        for(Int_t i = 0; i < subheader->GetN2();i++){
          TNudyEndfList* tab = (TNudyEndfList*)recIter.Next();
          lCoef[i].Set(tab->GetNPL());
          ein[i] = tab->GetC2();
          for(Int_t j = 0; j < tab->GetNPL(); j++)
            lCoef[i].SetAt(tab->GetLIST(j),j);
        }
        for(Int_t i = 0; i < ein.GetSize(); i++){
          printf("Ein = %e\n",ein[i]);
          for(Int_t j = 0; j < lCoef[i].GetSize();j++){
            printf("a%d = %e\n",j,lCoef[i].At(j));
          }
        }
      }
      else if(LTT == 3 && LI == 0){
        TNudyEndfTab2* lowE = (TNudyEndfTab2*)recIter.Next();
        TArrayD ein(lowE->GetN2());
        TArrayD *lCoef = new TArrayD[lowE->GetN2()]; 
        for(Int_t i = 0; i < lowE->GetN2();i++){
          TNudyEndfList* tab = (TNudyEndfList*)recIter.Next();
          lCoef[i].Set(tab->GetNPL());
          ein[i] = tab->GetC2();
          for(Int_t j = 0; j < tab->GetNPL(); j++)
            lCoef[i].SetAt(tab->GetLIST(j),j);
        }
        TNudyEndfTab2* highE = (TNudyEndfTab2*)recIter.Next();
        f4nens = highE->GetN2();
        f4eins.Set(f4nens);
        fAPAlias = new TNudyAliasCont[f4nens];
        TArrayD *tVal = new TArrayD[highE->GetN2()];
        ein.Set(ein.GetSize()+highE->GetN2());
        for(Int_t i = 0; i < highE->GetN2(); i++){
          TNudyEndfTab1* tab = (TNudyEndfTab1*)recIter.Next();
          f4eins.SetAt(tab->GetC2(),i);
          ein[lowE->GetN2() + i] = tab->GetC2();
          tVal[i].Set(tab->GetNP()*2);
          fAPAlias[i].Initialize(tab->Y(),tab->X(),tab->GetNP(),tab->GetC2(),i);
          for(Int_t j = 0; j < tab->GetNP(); j++){
            tVal[i].SetAt(tab->GetX(j),j*2);
            tVal[i].SetAt(tab->GetY(j),j*2+1);
          }
        }

        for(Int_t i = 0; i < lowE->GetN2(); i++){
          printf("Ein = %e\n",ein[i]);
          for(Int_t j = 0; j < lCoef[i].GetSize();j++){
            printf("a%d = %e\n",j,lCoef[i].At(j));
          }
        } 
        Double_t h = 0;
        for(Int_t i = 0; i < highE->GetN2(); i++){
          printf("Ein = %e\n",ein[i+lowE->GetN2()]);
          h = 0;
          for(Int_t j = 0; j < tVal[i].GetSize()/2;j++){
            if(j>1)
              h+= 0.5 * (tVal[i].At(j*2)-tVal[i].At(j*2-2)) * (tVal[i].At(j*2+1) + tVal[i].At(j*2-1));
            printf("X = %e Y = %e\n",tVal[i].At(j*2), tVal[i].At(j*2+1));
          }
          printf("Integral = %e",h);
        } 
        TNudyAliasCont::BuildIntermediate(fAPAlias,f4nens);
      }
    }
  }
}
//______________________________________________________________________________
void TVNudyModel::File5_Pass3(){
  nperc = 25;
  fPerc = new TArrayD[nens];
  fEPAlias = new TNudyAliasCont[nens];
  for(Int_t jen = 1; jen <= nens; jen++){
    fPerc[jen-1].Set(nperc+1);
    Int_t ibeg = 0,iend = 0;
    Int_t ipoint = 0;
    Double_t hintno = 0;
    Double_t aX[nEout[jen-1]];
    Double_t aP[nEout[jen-1]];
    for(Int_t jp = 1; jp <= nEout[jen-1]-1; jp++){
      hintno += 0.5*(fEPtable[jen-1].GetAt(2*jp-1)+fEPtable[jen-1].GetAt(2*jp+1))*(fEPtable[jen-1].GetAt(2*jp)-fEPtable[jen-1].GetAt(2*jp-2));
    }
    if(TMath::Abs(hintno-1) > 2e-2){
      Warning("File5_Pass3","Integral is not normalized at %e Integral = %e",xengr[jen-1],hintno);
    }
    if(hintno <= 0){
      for(Int_t jp = 1; jp <= nperc+1; jp++){
        fPerc[jen-1].SetAt(0,jp-1);
      }
    }
    else{
      Int_t np;
      for( np = 1; np <= nEout[jen-1]-1; np++){
        if(fEPtable[jen-1].GetAt(2*np+3) > 0){
          //    Info("Pass3","First Value = %e E = %e",fEPtable[jen-1].GetAt(2*np+2),xengr[jen-1]);
          fPerc[jen-1].SetAt(fEPtable[jen-1].GetAt(2*np),0);
          ibeg = np;
          break;
        }
      }
      if(!(np <= nEout[jen-1]-1)) Error("File5_Pass3", "No first point different from 0");
      for( np = nEout[jen-1]; np >= 1; np--){
        if(fEPtable[jen-1].GetAt(2*np-1) > 0){
          //  Info("Pass3","Last Value = %e E = %e",fEPtable[jen-1].GetAt(2*np-2),xengr[jen-1]);
          fPerc[jen-1].SetAt(fEPtable[jen-1].GetAt(2*np-2),nperc);
          iend = np;
          break;
        }
      }
      if(!(np >= 0)) Error("File5_Pass3","No last point different from 0");
      Double_t hinteg = 0;
      Double_t hintol = 0;
      ipoint = ibeg - 1;
      for(Int_t jperc = 1; jperc <= nperc-1 ; jperc++){
        Double_t percen = 1.0*jperc/nperc;
        while(hinteg < percen && ipoint < iend){
          hintol = hinteg;
          //      Info("F(x)","
          //  Info("int","ibeg = %d ipoint = %d, iend = %d , Int = %e",ibeg,ipoint,iend, hintol);
          ipoint++;
          hinteg+=0.5*(fEPtable[jen-1].GetAt(2*ipoint-1)+fEPtable[jen-1].GetAt(2*ipoint+1)) * (fEPtable[jen-1].GetAt(2*ipoint) - fEPtable[jen-1].GetAt(2*ipoint-2))/hintno;
        }
        Double_t f1 = fEPtable[jen-1].GetAt(2*ipoint-1)/hintno;
        Double_t f2 = fEPtable[jen-1].GetAt(2*ipoint+1)/hintno;
        Double_t e1 = fEPtable[jen-1].GetAt(2*ipoint-2);
        Double_t e2 = fEPtable[jen-1].GetAt(2*ipoint);
        Double_t delf = f2 - f1;
        Double_t hintrp;
        if(delf != 0){
          Double_t dele = e2 -e1;
          Double_t temp = 2*(percen-hintol)*delf/dele;
          hintrp = e1 + dele*(TMath::Sqrt(TMath::Power(f1,2)+temp)-f1)/delf;
        }
        else{
          hintrp = e1 + (percen-hintol)/f1;
        }
        fPerc[jen-1].SetAt(hintrp,jperc);
//        printf("Eout - %e, P - %f\n",hintrp,percen);

      }
    }
//    if(iend >= nEout[jen-1])
//      iend = nEout[jen-1] - 1;
//    while(iend > ibeg && fEPtable[jen-1].GetAt(2*iend+1) == 0) iend--;
    while(fEPtable[jen-1].GetAt(2*iend+1) == 0) iend--;

    printf("ENERGY = %e\n",xengr[jen-1]);
    for(Int_t ja = 0 ; ja <= iend;ja++){
      aX[ja] = fEPtable[jen-1].GetAt(2*ja);
      aP[ja] = fEPtable[jen-1].GetAt(2*ja + 1);
      printf("FX = %e FP = %e\n",aX[ja],aP[ja]);
    }
    fEPAlias[jen-1].Initialize(aP,aX,iend+1,xengr[jen-1]);
  }
  TNudyAliasCont::BuildIntermediate(fEPAlias,nens);
}

//______________________________________________________________________________
void TVNudyModel::File5_Pass2(TNudyEndfSec *sec){
  TIter recIter(sec->GetRecords());
  for(Int_t k = 0; k < sec->GetN1();k++){
    TNudyEndfTab1* header = (TNudyEndfTab1*)recIter.Next();
    CheckLinear(header);
    Double_t u = header->GetC1();
    SetTitle(Form("%s,%d",GetTitle(),header->GetL2()));
    if(header->GetL2() == 1) {
      TNudyEndfTab2* record = (TNudyEndfTab2*)recIter.Next();
      Int_t nep = record->GetN2();
      Double_t* engr = new Double_t[nep];
      TList temp;
      for(Int_t j = 1; j <= nep; j++){
        TNudyEndfTab1* row = (TNudyEndfTab1*)recIter.Next();
        temp.Add(row);
        engr[j-1]=row->GetC2();
      }	
      for(Int_t jg = 1; jg <= nens; jg++){
        Double_t energy = xengr[jg-1];
        Double_t ppe = TNudyCore::Instance()->Interpolate(header->NBT(),header->INT(),header->GetN1(),header->X(),header->Y(),header->GetN2(),energy);
        if(ppe <= 0) continue;
        Int_t ichan = TNudyCore::Instance()->BinarySearch(engr,nep,energy);
        if(engr[ichan] > energy || engr[ichan+1] < energy){
          Error("File5_Pass2","Error in Interpolation %e does not lie between %e and %e",energy,engr[ichan],engr[ichan+1]);
        }
        Double_t xrat = (engr[ichan+1]-energy)/(engr[ichan+1]-engr[ichan]);
        Double_t xra1 = 1-xrat;
        if(TMath::Min(xrat,xra1) < 0 || TMath::Max(xrat,xra1) > 1){
          Error("File5_Pass2","Error in XRAT, XRA1 = %e, %e\nE = %e Ei = %e Eo = %e",xrat,xra1,energy,engr[ichan],engr[ichan+1]);

        }
        Double_t hint = 0;
        Double_t eold = 0;
        Double_t pold = 0;
        TNudyEndfTab1* lfunc = (TNudyEndfTab1*)temp.At(ichan);
        TNudyEndfTab1* ufunc = (TNudyEndfTab1*)temp.At(ichan+1);
        for(Int_t jef = 1; jef <= nEout[jg-1]; jef++){
          Double_t efen = fEPtable[jg-1].GetAt(2*jef-2);
          Double_t pb1 = TNudyCore::Instance()->Interpolate(lfunc->NBT(),lfunc->INT(),lfunc->GetN1(),lfunc->X(),lfunc->Y(),lfunc->GetN2(),efen);
          Double_t pb2 = TNudyCore::Instance()->Interpolate(ufunc->NBT(),ufunc->INT(),ufunc->GetN1(),ufunc->X(),ufunc->Y(),ufunc->GetN2(),efen);
          Double_t prob = xrat * pb1 + xra1 * pb2;
          //    Info("File5_Pass2:Integration","e1,p1 = %e,%e , e2,p2 = %e,%e at E = %e Int = %e A=%e B=%e",eold,pold,efen,prob,energy, hint,pb1,pb2); 
          fEPtable[jg-1].GetArray()[2*jef-1] += ppe*prob;
          if(jef > 1){
            hint=hint+(0.5*(efen-eold)*(prob+pold));
          }
          eold = efen;
          pold = prob;
        }
        if(TMath::Abs(1-hint) > 2e-2){
          Error("File5_Pass2","Integral error = %e at E = %e PPE = %e",hint*ppe, energy,ppe);
        }
        else{
          Info("File5_Pass2","Succesful Integral %e at E = %e PPE = %e",hint,energy,ppe);
        }
      }
    }
    else if(header->GetL1() == 5){
      Double_t u = header->GetC1();
      TNudyEndfTab1* temptab =(TNudyEndfTab1*) recIter.Next();
      TNudyEndfTab1* probtab = (TNudyEndfTab1*)recIter.Next();
      for(Int_t jg = 1; jg <= nens; jg++){
        Double_t energy = xengr[jg-1];
        if(energy <= u ) continue;
        Double_t ppe = TNudyCore::Instance()->Interpolate(header->NBT(),header->INT(),header->GetN1(),header->X(),header->Y(),header->GetN2(),energy);
        if(ppe <= 0) continue;
        Double_t teta = TNudyCore::Instance()->Interpolate(temptab->NBT(),temptab->INT(),temptab->GetN1(),temptab->X(),temptab->Y(),temptab->GetN2(),energy);
        Double_t hint = 0;
        Double_t pold = 0;
        Double_t eold = 0;
        for(Int_t jef = 1; jef <= nEout[jg-1]; jef++){
          Double_t efen = fEPtable[jg-1].GetAt(2*jef-2);
          if(!(efen >= (energy - u))){
            Double_t prob = TNudyCore::Instance()->Interpolate(probtab->NBT(),probtab->INT(),probtab->GetNR(),probtab->X(),probtab->Y(),probtab->GetN2(),efen/teta);
            fEPtable[jg-1].GetArray()[2*jef-2]+= ppe*prob;
            if(jef > 1)
              hint += 0.5*(efen-eold)*(prob+pold);
            eold = efen;
            pold = prob;
          }else{
            break;
          }

        }
        if(TMath::Abs(1-hint) > 2e-2){
          Error("File5_Pass2","Integral error = %e at E = %e",hint*ppe, energy);
        }
        else{
          Info("File5_Pass2","Succesful Integral %e at E = %e PPE = %e",hint,energy,ppe);
        }
      }
    }
    else if(header->GetL2() == 7){
      Double_t u = header->GetC1();
      TNudyEndfTab1* temptab = (TNudyEndfTab1*)recIter.Next();
      for(Int_t jg = 1; jg <= nens; jg++){
        Double_t energy = xengr[jg-1];
        Double_t ppe = TNudyCore::Instance()->Interpolate(header->NBT(),header->INT(),header->GetN1(),header->X(),header->Y(),header->GetN2(),energy);
        if(ppe <= 0) continue;
        Double_t teta = TNudyCore::Instance()->Interpolate(temptab->NBT(),temptab->INT(),temptab->GetN1(),temptab->X(),temptab->Y(),temptab->GetN2(),energy);
        Double_t rede = (energy-u)/teta;
        Double_t hnorm = TMath::Power(teta,1.5)*(TMath::Sqrt(TMath::Pi())/2*TMath::Erf(TMath::Sqrt(rede)) - TMath::Sqrt(rede)*TMath::Exp(-rede));
        Double_t hint = 0;
        Double_t pold = 0;
        Double_t eold = 0;
        for(Int_t jef = 1; jef <= nEout[jg-1]; jef++){
          Double_t efen = fEPtable[jg-1].GetAt(2*jef-2);
          if(!(efen >= (energy-u))){
            Double_t prob = TMath::Sqrt(efen)*TMath::Exp(-efen/teta)/hnorm;
            fEPtable[jg-1].GetArray()[2*jef-1]+= ppe*prob;
            if(jef > 1)
              hint += 0.5*(efen-eold)*(prob+pold);
            eold = efen;
            pold = prob;
          }
          else{
            break;}
        }
        if(TMath::Abs(1-hint) > 2e-2){
          Error("File5_Pass2","Integral error = %e at E = %e",hint*ppe, energy);
        }
        else{
          Info("File5_Pass2","Succesful Integral %e at E = %e PPE = %e",hint,energy,ppe);
        }
      }
    }
    else if(header->GetL2() == 9){
      Double_t u = header->GetC1();
      TNudyEndfTab1* temptab = (TNudyEndfTab1*)recIter.Next();
      for(Int_t jg = 1; jg <= nens; jg++){
        Double_t energy = xengr[jg-1];
        if(energy <= u) continue;
        Double_t ppe = TNudyCore::Instance()->Interpolate(header->NBT(),header->INT(),header->GetN1(),header->X(),header->Y(),header->GetN2(),energy);
        if(ppe <= 0) continue;
        Double_t teta = TNudyCore::Instance()->Interpolate(temptab->NBT(),temptab->INT(),temptab->GetN1(),temptab->X(),temptab->Y(),temptab->GetN2(),energy);
        Double_t rede = (energy-u)/teta;
        Double_t hnorm;
        if(rede > 1e-6)
          hnorm = TMath::Power(teta,2)*(1-TMath::Exp(-rede)*(1+rede));
        else
          hnorm = TMath::Power(teta,2)*0.5*TMath::Power(rede,2)*(1-rede);
        Double_t hint = 0;
        Double_t pold = 0;
        Double_t eold = 0;
        for(Int_t jef = 1; jef <= nEout[jg-1]; jef++){
          Double_t efen = fEPtable[jg-1].GetAt(2*jef-2);
          if(!(efen >= (energy-u))){
            Double_t prob = TMath::Sqrt(efen)*TMath::Exp(-efen/teta)/hnorm;
            fEPtable[jg-1].GetArray()[2*jef-1]+= ppe*prob;
            if(jef > 1)
              hint += 0.5*(efen-eold)*(prob+pold);
            eold = efen;
            pold = prob;
          }
          else{
            break;}
        }
        if(TMath::Abs(1-hint) > 2e-2){
          Error("File5_Pass2","Integral error = %e at E = %e",hint*ppe, energy);
        }
        else{
          Info("File5_Pass2","Succesful Integral %e at E = %e PPE = %e",hint,energy,ppe);
        }		
      }
    }
    else if(header->GetL2() == 11){
      Double_t u = header->GetC1();
      TNudyEndfTab1* atab = (TNudyEndfTab1*)recIter.Next();
      TNudyEndfTab1* btab = (TNudyEndfTab1*)recIter.Next();
      for(Int_t jg = 1; jg <= nens; jg++){
        Double_t energy = xengr[jg-1];
        if(energy <= u) continue;
        Double_t ppe = TNudyCore::Instance()->Interpolate(header->NBT(),header->INT(),header->GetN1(),header->X(),header->Y(),header->GetN2(),energy);
        if(ppe <= 0) continue;
        Double_t a = TNudyCore::Instance()->Interpolate(atab->NBT(),atab->INT(),atab->GetN1(),atab->X(),atab->Y(),atab->GetN2(),energy);
        Double_t b = TNudyCore::Instance()->Interpolate(btab->NBT(),btab->INT(),btab->GetN1(),btab->X(),btab->Y(),btab->GetN2(),energy);
        Double_t ab4 = 0.25*a*b;
        Double_t sab4 = TMath::Sqrt(ab4);
        Double_t elim = energy-u;
        Double_t rede = (elim)/a;
        Double_t srede = TMath::Sqrt(rede);
        Double_t hnorm = 0.5*a*TMath::Sqrt(TMath::Pi())*sab4*TMath::Exp(ab4)*(TMath::Erf(srede-sab4)+TMath::Erf(srede+sab4))-a*TMath::Exp(-rede)*TMath::SinH(TMath::Sqrt(b*elim));
        Double_t hint = 0;
        Double_t pold = 0;
        Double_t eold = 0;
        for(Int_t jef = 1; jef <= nEout[jg-1]; jef++){
          Double_t efen = fEPtable[jg-1].GetAt(2*jef-2);
          if(!(efen >= (energy-u))){
            Double_t prob = TMath::Exp(-efen/a)*TMath::SinH(TMath::Sqrt(b*efen))/hnorm;
            fEPtable[jg-1].GetArray()[2*jef-1]+= ppe*prob;
            if(jef > 1)
              hint += 0.5*(efen-eold)*(prob+pold);
            eold = efen;
            pold = prob;
          }
          else{
            break;}
        }
        if(TMath::Abs(1-hint) > 2e-2){
          Error("File5_Pass2","Integral error = %e at E = %e",hint*ppe, energy);
        }
        else{
          Info("File5_Pass2","Succesful Integral %e at E = %e PPE = %e",hint,energy,ppe);
        }
      }
    }
    else if(header->GetL2() == 12){
      Double_t u = header->GetC1();
      TNudyEndfTab1* temptab = (TNudyEndfTab1*)recIter.Next();
      Double_t ef[2] = {temptab->GetC1(),temptab->GetC2()};
      for(Int_t jg = 1; jg <= nens; jg++){
        Double_t energy = xengr[jg-1];
        if(energy <= u) continue;
        Double_t elim = energy-u;
        Double_t ppe = TNudyCore::Instance()->Interpolate(header->NBT(),header->INT(),header->GetN1(),header->X(),header->Y(),header->GetN2(),energy);
        if(ppe <= 0) continue;
        Double_t tm = TNudyCore::Instance()->Interpolate(temptab->NBT(),temptab->INT(),temptab->GetN1(),temptab->X(),temptab->Y(),temptab->GetN2(),energy);
        Double_t hint = 0;
        Double_t *tmppb = new Double_t[nEout[jg-1]];
        for(Int_t jef = 1; jef <= nEout[jg-1]; jef++){
          Double_t efen = fEPtable[jg-1].GetAt(2*jef-2);
          if(!(efen >= (elim))){
            Double_t prob = 0;
            for(Int_t jhl = 0; jhl < 2; jhl++){
              Double_t u1 = TMath::Power(TMath::Sqrt(efen)-TMath::Sqrt(ef[jhl]),2)/tm;
              Double_t u2 = TMath::Power(TMath::Sqrt(efen)+TMath::Sqrt(ef[jhl]),2)/tm;
              prob += 0.5*((TMath::Power(u2,1.5)*ROOT::Math::expint(u2)+TMath::Gamma(1.5,u2))-(TMath::Power(u1,1.5)*ROOT::Math::expint(u1)+TMath::Gamma(1.5,u1)))/(3*TMath::Sqrt(ef[jhl]*tm));
            }
            tmppb[jef-1] = prob;

            if(jef > 1)
              hint += 0.5*(tmppb[jef-2]+tmppb[jef-1])*(fEPtable[jg-1].GetAt(2*jef-1) - fEPtable[jg-1].GetAt(2*jef-3));
          }
          else{
            break;}
        }
        for(Int_t j = 1; j <= nEout[jg-1]; j++)
          fEPtable[jg-1].GetArray()[2*j-1]+=ppe*tmppb[j-1]/hint;

      }

    }
  }
}

//______________________________________________________________________________
void TVNudyModel::File5_Pass1(TNudyEndfSec* sec){
  Int_t index = 0;
  TIter recIter(sec->GetRecords());
  fEPtable = new TArrayD[500];
  for(Int_t k = 0; k < sec->GetN1();k++){
    Info("error","%d < %d",k,sec->GetN1());
    TNudyEndfTab1* header = (TNudyEndfTab1*)recIter.Next();
    CheckLinear(header);
    Double_t u = header->GetC1();
    Info("File5_Pass1","%s - File 5 - LF %d",GetName(),header->GetL2());
    if(header->GetL2() == 1) {
      TNudyEndfTab2* range = (TNudyEndfTab2*)recIter.Next();
      CheckLinear(range);
      Int_t nz;
      TNudyEndfTab1 *row;
      for(nz = 0; nz < range->GetN2(); nz++){
	row = (TNudyEndfTab1*)recIter.Next();
	//Linearize Interplolation if it is histogram like
	Linearize(row);
	//Add upper threshold value
/*	if(row->GetX(0) == 0 && row->GetY(0) == 0){
	  TNudyEndfTab1 *temp = new TNudyEndfTab1(row->GetC1(),row->GetC2(),row->GetL1(),row->GetL2(),row->GetN1(),row->GetN2()-1);
	  //memcpy(temp->X(),row->X()+1,sizeof(Double_t)*row->GetN2()-1);
	  //		    memcpy(temp->Y(),row->Y()+1,sizeof(Double_t)*row->GetN2()-1);
	  //		    memcpy(temp->INT(),row->INT(),sizeof(Int_t)*row->GetN1());
	  for(Int_t j = 1; j < row->GetN2(); j++){
	    temp->SetX(row->GetX(j),j-1);
	    temp->SetY(row->GetY(j),j-1);
	  }
	  for(Int_t i = 0; i < row->GetN1(); i++){
	    temp->SetNBT(row->GetNBT(i)-1,i);
	    temp->SetINT(row->GetINT(i),i);
	  }
	  row->Equate(temp);
	  delete temp;
	}*/
	if(row->GetY(row->GetN2()-1) > 0){
	  TNudyEndfTab1 *temp = new TNudyEndfTab1(row,row->GetN1(),row->GetN2()+1);
	  temp->SetX(temp->GetX(row->GetN2()-1)*(1+1e-4),row->GetN2());
	  temp->SetY(1e-10,row->GetN2());
	  row->Equate(temp);
	  delete temp;
	}
	Double_t tempen = row->GetC2();
	Double_t ppe = TNudyCore::Instance()->Interpolate(header->NBT(),header->INT(),header->GetNR(),header->X(),header->Y(),header->GetN2(),tempen);
	if(ppe <= 0)
	  continue;
	Int_t exists = 0;
	for(Int_t jen = 0; jen < nens; jen++){
	  if(tempen == xengr[jen]){
	    index = jen;
	    exists = 1;
	    break;
	  }
	}
	if(!exists){
	  nens++;	    
	  index=nens-1;
	  xengr.Set(nens);
	  nEout.Set(nens);
	  xengr[index]=tempen;
	  Info("Pass1","Adding %e energy at %d",tempen,index);
	  nEout[index]=0;
	  Double_t emin = TMath::Max(row->GetX(0)*(1-1e-5),1e-5);
	  Double_t emax = row->GetX(row->GetN2()-1)*(1+1e-5);
	  Double_t fact = TMath::Exp(TMath::Log(emax/emin)/(maxpop-1));
	  Double_t ef = emin/fact;
	  for(Int_t jef = 0; jef < maxpop; jef++){
	    ef = ef*fact;
	    if(!EoExists(index,ef)){	
	      //		    Info("Pass1","Adding outgoing energy %e at %d",ef,index);
	      nEout[index]++;
	      fEPtable[index].Set(2*nEout[index]);
	      fEPtable[index].SetAt(ef,2*nEout[index]-2);
	      fEPtable[index].SetAt(0,2*nEout[index]-1);
	    }
	  }
	  for(Int_t jeps = 1; jeps <= row->GetNP(); jeps++){
	    ef = row->GetX(jeps-1);
	    //Separating equal probability secondary energies
	    if(jeps-1 < row->GetNP() && ef == row->GetX(jeps)){
	      ef = ef*(1-5e-6);
	      row->SetX(row->GetX(jeps)*(1+5e-6),jeps);
	    }
	    if(!EoExists(index,ef)){
	      nEout[index]++;
	      fEPtable[index].Set(2*nEout[index]);
	      fEPtable[index].SetAt(ef,2*nEout[index]-2);
	      fEPtable[index].SetAt(0,2*nEout[index]-1);
	    }
	  }
	}
      }
    }
    else if(header->GetL2() == 5){	 
      TNudyEndfTab1 *theta = (TNudyEndfTab1*)recIter.Next();
      FillGrid(u,theta->GetN2(),theta,header);
      TNudyEndfTab1 *dEdTheta = (TNudyEndfTab1*)recIter.Next();
      for(index = 0;  index < nens; index++){
	Double_t energy = xengr[index];
	Double_t tet = TNudyCore::Instance()->Interpolate(theta->NBT(),theta->INT(),theta->GetN1(),theta->X(),theta->Y(),theta->GetN2(),energy);
	for(Int_t jef = 1; jef <= dEdTheta->GetN2();jef++){
	  Double_t ef = tet*dEdTheta->GetX(jef-1);
	  if(!EoExists(index,ef)){
	    nEout[index]++;
	    fEPtable[index].Set(2*nEout[index]);
	    fEPtable[index].SetAt(ef,2*nEout[index]-2);
	    fEPtable[index].SetAt(0,2*nEout[index]-1);
	  }
	}
      }
    }
    else if(header->GetL2() == 7){	      
      TNudyEndfTab1 *theta =  (TNudyEndfTab1*)recIter.Next();
      FillGrid(u,theta->GetN2(),theta,header);
    }
    else if(header->GetL2() == 9){	      
      TNudyEndfTab1 *theta =  (TNudyEndfTab1*)recIter.Next();
      FillGrid(u,theta->GetN2(),theta,header);
    }
    else if(header->GetL2() == 11){	      
      TNudyEndfTab1 *a =  (TNudyEndfTab1*)recIter.Next();
      FillGrid(u,a->GetN2(),a,header);
      TNudyEndfTab1 *b =  (TNudyEndfTab1*)recIter.Next();
      FillGrid(u,b->GetN2(),b,header);
    }
    else if(header->GetL2() == 12){	      
      TNudyEndfTab1 *tm =  (TNudyEndfTab1*)recIter.Next();
      FillGrid(u,tm->GetN2(),tm,header);
    }
  }
  for(Int_t i = 1; i <= nens-1; i++){
    for(Int_t j = i+1; j <= nens; j++){
      if(xengr[i-1] > xengr[j-1]){
	Double_t temp = xengr[i-1];
	xengr[i-1] = xengr[j-1];
	xengr[j-1] = temp;
	Int_t ntemp = nEout[i-1];
	nEout[i-1] = nEout[j-1];
	nEout[j-1] = ntemp;
	char buffer[sizeof(TArrayD)];
	memcpy(buffer, &fEPtable[i-1], sizeof(TArrayD));
	memcpy(&fEPtable[i-1], &fEPtable[j-1], sizeof(TArrayD));
	memcpy(&fEPtable[j-1], buffer, sizeof(TArrayD));
      }
    }
  }
  for(Int_t jen = 0; jen < nens; jen++){
    for(Int_t j1 = 0; j1 < nEout[jen]-1; j1++){
      for(Int_t j2 = j1+1; j2 < nEout[jen]; j2++){
	if(fEPtable[jen].GetAt(2*j1) > fEPtable[jen].GetAt(2*j2)){
	  Double_t tmp = fEPtable[jen].GetAt(2*j1);
	  fEPtable[jen].GetArray()[2*j1] = fEPtable[jen].GetAt(2*j2);
	  fEPtable[jen].GetArray()[2*j2] = tmp;
	}
      }
    }
  }
}

//______________________________________________________________________________
void TVNudyModel::PopulateGrid(Int_t index){
  Int_t max = 50;
  Double_t emin = 1e-5;
  Double_t emax = 2e7;
  Double_t fact,ef;
  if(emax > 0){
    fact = TMath::Exp(TMath::Log(emax/emin)/(max-1));
    ef = emin/fact;
  } else {

  }
  if(index > nEout.GetSize()-1)
    nEout.Set(index+1);
  nEout[index]=max;
  fEPtable[index].Set(2*max);
  for(Int_t jeps = 1; jeps <= max;jeps++){
    ef = ef*fact;
    fEPtable[index].GetArray()[2*jeps-2] = ef;
    fEPtable[index].GetArray()[2*jeps-1] = 0;
  }
}

//______________________________________________________________________________
void TVNudyModel::FillGrid(Double_t u, Int_t nep,TNudyEndfTab1 *tab ,TNudyEndfTab1 *pe){
  Int_t nene = 0;
  TArrayD tmpene;
  Int_t j,keps,index;
  Double_t ratmax = 2;
  Double_t ratio,diff,ratiol,nadd,eadd,fact,ef,ethre1,ethre2;
  Int_t maxene = 200;
  tmpene.Set(maxene);
  for( j = 1; j <= nep; j++){
    Double_t et = tab->GetX(j-1);
    Double_t ppe = TNudyCore::Instance()->Interpolate(pe->NBT(),pe->INT(),pe->GetNR(),pe->X(),pe->Y(),pe->GetN2(),et);
    if(ppe <= 0) continue;
    if(et - u < 1e-5){
      et = TMath::Max(u*1.001,u+1e-4);
    }
    nene++;
    tmpene[nene-1] = et;
  }
  Int_t ntota = 0;
  for( j = 1; j <= nene-1; j++){
    ratio = tmpene[j]/tmpene[j-1];
    diff = 2*TMath::Abs(tab->GetY(j-1)-tab->GetY(j))/TMath::Max(1e-10,tab->GetY(j-1)+tab->GetY(j));
    if(ratio > ratmax && diff > 1e-6){
      ratiol = TMath::Log(ratio);
      nadd = ratiol/TMath::Log(ratmax)+1;
      fact = TMath::Exp(ratiol/nadd);
      eadd = tmpene[j-1];
      for(Int_t jad = 1; jad <= nadd-1; jad++){
	eadd = eadd*fact;
	tmpene[nene+ntota+jad-2] = eadd;
      }
      ntota = ntota+nadd-1;
    }
  }
  nene = nene+ntota;
  for(j = 1 ; j <= nene; j++){
    Int_t check = 0;
    for(Int_t jen = 1;jen <= nens; jen++){
      if(TMath::Abs(2*(tmpene[j-1]-xengr[jen-1])/(tmpene[j-1]+xengr[jen-1])) < 1e-7){
	index = jen-1;
	check = 1;
	break;
      }
    } 
    if(!check){
      nens++;
      index = nens-1;
      xengr.Set(nens);
      xengr[index]=tmpene[j-1];
      PopulateGrid(index);
    }
    Double_t emin = 1e-5;
    Double_t emax = TMath::Max(xengr[index]-u,emin+1e-5);
    if(emax > 0){
      fact = TMath::Exp(TMath::Log(emax/emin)/(maxene-1));
      ef = emin/fact;
    } else {
      fact = 1;
      ef = 0;
    }
    for(Int_t jeps = 1; jeps <= maxene;jeps++){
      ef = ef*fact;
      if(!EoExists(index,ef)){
	nEout[index]++;
	fEPtable[index].Set(2*nEout[index]);
	fEPtable[index].SetAt(ef,2*nEout[index]-2);
	fEPtable[index].SetAt(0,2*nEout[index]-1);
      }
    }
    check = 1;
    ethre1 = emax * (1-1e-4);
    
    for( keps = 1; keps <= nEout[index]; keps++){
      if(fEPtable[index].GetAt(2*keps-2) == ethre1){
	check=0;
	break;
      }
    }
    if(check){
      nEout[index]++;
      fEPtable[index].Set(2*nEout[index]);
      fEPtable[index].SetAt(ethre1,2*(nEout[index])-2);
      fEPtable[index].SetAt(0,2*(nEout[index])-1);
    }
    ethre2 = emax * (1+1e-4);
    check = 1;
    for( keps = 1; keps <= nEout[index]; keps++){
      if(fEPtable[index].GetAt(2*keps-2) == ethre2){
	check=0;
	break;
      }
    }
    if(check){
      nEout[index]++;
      fEPtable[index].Set(2*nEout[index]);
      fEPtable[index].SetAt(ethre2,2*nEout[index]-2);
      fEPtable[index].SetAt(0,2*nEout[index]-1);
    }
  }
}

//______________________________________________________________________________
Int_t TVNudyModel::EoExists(Int_t index,Double_t ef){
  Int_t check = 0;
  for(Int_t keps = 1; keps <= nEout[index]; keps++){
    if(fEPtable[index].GetAt(2*keps-2)==ef){
      check = 1;
      break;
    }
  }
  return check;
}

//______________________________________________________________________________
void TVNudyModel::Linearize(TNudyEndfTab1 *tab){
  if(!CheckLinear(tab)){
    Int_t islin = 0;
    Int_t start = 1;

    Int_t n2 = tab->GetN2();
    Double_t *epval = new Double_t[tab->GetN2()*4];
    for(Int_t i = 0; i < tab->GetN2(); i++){
      epval[2*i] = tab->GetX(i);
      epval[2*i+1] = tab->GetY(i);
    }
    for(Int_t jr = 1; jr <= tab->GetN1(); jr++){
      if(tab->GetINT(jr-1) == 1){
	islin = 1;
	for(Int_t jp = start; jp <= tab->GetNBT(jr-1)-1; jp++){
	  n2 = n2 +1;
	  epval[2*n2-2] = TMath::Max((1-1e-6)*epval[2*jp],epval[jp*2-2]);
	  epval[2*n2-1] = epval[jp*2-1];
	}
      }
      else{
	start = tab->GetNBT(jr-1);
	if(tab->GetINT(jr)!=2){
	  Warning("Linearize","Unexpected interpolation law");
	}
      }
    }
    if(!islin) return;
    for(Int_t j1 = 1; j1 <= n2; j1++){
      for(Int_t j2 = j1+1; j2 <= n2; j2++){
	if(epval[2*j1-2] > epval[2*j2-2]){
	  Double_t tmp1 = epval[2*j1-2];
	  Double_t tmp2 = epval[2*j1-1];
	  epval[2*j1-2] = epval[2*j2-2];
	  epval[2*j1-1] = epval[2*j2-1];
	  epval[2*j2-2] = tmp1;
	  epval[2*j2-1] = tmp2;
	}
      }
    }
    TNudyEndfTab1 *linear = new TNudyEndfTab1(tab->GetC1(),tab->GetC2(),tab->GetL1(),tab->GetL2(),1,n2);
    linear->SetINT(2,0);
    linear->SetNBT(n2,0);
    for(Int_t i = 0; i < n2; i++){
      linear->SetX(epval[2*i],i);
      linear->SetY(epval[2*i+1],i);
    }
    tab->Equate(linear);
    delete linear;
  }
}

//_______________________________________________________________________________
Int_t TVNudyModel::CheckLinear(TNudyEndfTab1* tab){
  for(Int_t nr = 0; nr < tab->GetNR(); nr++){
    if(tab->GetINT(nr) != 2){
      Error("CheckLinear","Tab1 record data is not interpolated linearly");
      return 0;
    }
  }
  return 1;
}
Int_t TVNudyModel::CheckLinear(TNudyEndfTab2* tab){
  for(Int_t nr = 0; nr < tab->GetNR(); nr++){
    if(tab->GetINT(nr) != 2){
      Error("CheckLinear","Tab2 record data is not interpolated linearly");
      return 0;
    }
  }
  return 1;
}

//_______________________________________________________________________________
void TVNudyModel::DisplayData(FileData_t file){
  //This function displays data for a particular file stored in this model
  //Points drawn in graph are those from data.
  TCanvas *c1 = new TCanvas("c1","Display Data",200,10,700,500);
  c1->SetFillColor(42);
  c1->SetGrid();
  TCanvas *canvas;
  TGraph2D *gr2;
  TGraph *gr = NULL;
  //Draw data
  switch(file) {
  case kReac_XSect: 
    if(!fE_file3 || !fXSect_file3){
      Error("DisplayData()","This model does not contain any data about file 3");
      return;
    }
    c1->SetLogy();
    c1->SetLogx();
    gr = new TGraph(fEXSect_length,fE_file3,fXSect_file3);
    gr->SetTitle("File 3");
    gr->GetXaxis()->SetTitle("Energy");
    gr->GetYaxis()->SetTitle("Cross Section");
    break;
  case kEnergy_Dist:
    if(!fPerc){
      Error("DisplayData()","This model does not contain any data about file 5");
      return;
    }
    gr2 = new TGraph2D();
    for(Int_t i = 0; i < nens; i++){
      for(Int_t j = 1; j <= 24; j++){
        gr2->SetPoint(i*24+j,xengr[i],j/24.0,fPerc[i].GetAt(j));
      }
    }
    gStyle->SetPalette(1);
    gr2->SetName(GetMaterialName());
    gr2->SetTitle("Energy Probability Distribution");
    gr2->GetXaxis()->SetTitle("Primary Energy (Ein)");
    gr2->GetYaxis()->SetTitle("Probability (P)");
    gr2->GetZaxis()->SetTitle("Secondary Energy (Eo)");
    gr2->Draw("surf1");
    break;
  default:
    break;
  }
  if(gr){
    gr->SetLineColor(2);
    gr->SetLineWidth(1);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->Draw("al");
  }
  if(gr2){
    gr2->Draw("surf1");
  }
  c1->Update();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(5);
  c1->Modified();
  c1->SaveAs(Form("%d_%d_%d_%d_%d.eps",fMAT,fTemp,fPdg,fReaction,file));
  c1->Close();
  if(!gr){
    delete gr;
  }
  if(!gr2){
    delete gr2;
  }
  delete c1;
}

//_______________________________________________________________________________
void TVNudyModel::DumpData(FileData_t file){

  switch (file){
  case kEnergy_Dist:
    printf("Energy Probability Distribution for %s\n\n",fMaterial->GetName());
    break;
  default:return;
  }
}

//______________________________________________________________________________
TArrayD* TVNudyModel::GetFile5Data(Double_t ein) {
  register Int_t lo,hi,mid;
  Int_t found = -1;
  lo = 0;
  hi = nens - 1;
  while(hi-lo>1) {
    mid = (lo+hi)/2;
    if(ein>xengr.GetAt(mid)) lo = mid;
    else if(ein<xengr.GetAt(mid)) hi = mid;
    else {
      found = mid;
      break;
    }
  }
  if(found>=0)
    return &fEPtable[found];
  else {
    printf("hi=%d,lo=%d\n",hi,lo);
    Error("getFile5Data","No such Ein in data");
    return NULL;
  }
}

//______________________________________________________________________________
TArrayD* TVNudyModel::GetFile5ProcessedData(Double_t ein) {
  register Int_t lo,hi,mid;
  Int_t found = -1;
  lo = 0;
  hi = nens - 1;
  while(hi-lo>1) {
    mid = (lo+hi)/2;
    if(ein>xengr.GetAt(mid)) lo = mid;
    else if(ein<xengr.GetAt(mid)) hi = mid;
    else {
      found = mid;
      break;
    }
  }
  if(found>=0)
    return &fPerc[found];
  else {
    printf("hi=%d,lo=%d\n",hi,lo);
    Error("getFile5Data","No such Ein in data");
    return NULL;
  }
}

