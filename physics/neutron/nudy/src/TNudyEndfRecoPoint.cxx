// This class is reconstructing ENDF data and creating probability tables for the angle
// and energy distributions of the secondatries
// Author: Dr. Harphool Kumawat
// Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// date of creation: March 22, 2016

#include "TList.h"
#include "TFile.h"
#include "TKey.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfTape.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfMat.h"
#include "Geant/TNudyEndfNuPh.h"
#include "Geant/TNudyEndfAng.h"
#include "Geant/TNudyEndfEnergy.h"
#include "Geant/TNudyEndfEnergyAng.h"
#include "Geant/TNudyEndfFissionYield.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfPhYield.h"
#include "Geant/TNudyEndfPhProd.h"
#include "Geant/TNudyEndfPhAng.h"
#include "Geant/TNudyEndfPhEnergy.h"
#include "Geant/TNudyEndfRecoPoint.h"
#include "TTree.h"
#include "TH1D.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfRecoPoint)
#endif

#ifdef USE_ROOT
#include "TRandom3.h"
#endif

    TNudyEndfRecoPoint::TNudyEndfRecoPoint()
    : fElemId(0), rENDF()
{
}
TNudyEndfRecoPoint::TNudyEndfRecoPoint(int ielemId, const char *irENDF) : fElemId(ielemId), rENDF(irENDF)
{
  GetData(ielemId, irENDF);
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ReadFile3(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  int mt1multi = -1;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    int MT = sec->GetMT();
    if (MT != 1 && MT != 3 && MT != 4 && MT != 27 && MT != 19 && MT != 20 && MT != 21 && MT != 38 && MT != 101 &&
        (MT < 120 || MT >= 600)) {
      fMtNumbers.push_back(MT);
      TIter recIter(sec->GetRecords());
      TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
      double QI             = header->GetC2();
      fQValue[sec->GetMT()] = QI;
      fQvalueTemp.push_back(QI);
      fQvalueTemp.push_back(MT);
      fNR                 = header->GetN1();
      fNP                 = header->GetN2();
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab1->GetX(crs));
        fXLinearFile3.push_back(tab1->GetY(crs));
      }
      fEneTemp.insert(std::end(fEneTemp), std::begin(fELinearFile3), std::end(fELinearFile3));
      fEneTemp.insert(std::end(fEneTemp), std::begin(fXLinearFile3), std::end(fXLinearFile3));
      fELinearFile3.clear();
      fXLinearFile3.clear();
      fSigmaOfMts.push_back(fEneTemp);
      fEneTemp.clear();
    }
    if (MT == 1 && mt1multi == -1) {
      mt1multi = 0;
      TIter recIter(sec->GetRecords());
      TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
      fNR                   = header->GetN1();
      fNP                   = header->GetN2();
      TNudyEndfTab1 *tab1   = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab1->GetX(crs));
        fXLinearFile3.push_back(tab1->GetY(crs));
      }
      fEnergyMts.insert(std::end(fEnergyMts), std::begin(fELinearFile3), std::end(fELinearFile3));
      fSigmaMts.insert(std::end(fSigmaMts), std::begin(fXLinearFile3), std::end(fXLinearFile3));
      fELinearFile3.clear();
      fXLinearFile3.clear();
    }
  }
  MtValues.push_back(fMtNumbers);
  fQvalue.push_back(fQvalueTemp);
  fQvalueTemp.clear();
  fMtNumbers.clear();
}
void TNudyEndfRecoPoint::GetData(int ielemId, const char *rENDF)
{
  // TFile *f = new TFile("testetau235.root", "recreate");
  // TTree *eta =new TTree("eta","a Tree with fast fission data");
  // TH1D *h1 = new TH1D("h1", "", 150, -5, 8);
  // TH1D *h2 = new TH1D("h2", "", 150, -5, 8);
  // TH1D *h3 = new TH1D("h3", "", 150, -5, 8);
  fElemId     = ielemId;
  TFile *rEND = TFile::Open(rENDF);
  //  TFile *rEND = TFile::Open(rENDF,"UPDATE");
  if (!rEND || rEND->IsZombie()) printf("Error: TFile :: Cannot open file %s\n", rENDF);
  TKey *rkey              = (TKey *)rEND->GetListOfKeys()->First();
  TNudyEndfTape *rENDFVol = (TNudyEndfTape *)rkey->ReadObj();
  TNudyEndfMat *tMat      = 0;
  TList *mats             = (TList *)rENDFVol->GetMats();
  int nmats               = mats->GetEntries();
  int mt455               = 1000;
  for (int iMat = 0; iMat < nmats; iMat++) {
    tMat = (TNudyEndfMat *)mats->At(iMat);
    for (int inf = 0; inf < tMat->GetNXC(); inf++) {
      switch (tMat->GetMFn(inf)) {
      case 4:
        fMtNum4.push_back(tMat->GetMFn(inf));
        fMtNum4.push_back(tMat->GetMTn(inf));
        break;
      case 5:
        fMtNum5.push_back(tMat->GetMFn(inf));
        if (tMat->GetMTn(inf) == 455) {
          fMtNum5.push_back(mt455);
          mt455++;
        } else {
          fMtNum5.push_back(tMat->GetMTn(inf));
        }
        break;
      case 6:
        fMtNum6.push_back(tMat->GetMFn(inf));
        fMtNum6.push_back(tMat->GetMTn(inf));
        break;
      }
    }
    fMt4.push_back(fMtNum4);
    fMt5.push_back(fMtNum5);
    fMt6.push_back(fMtNum6);
    std::vector<int>().swap(fMtNum4);
    std::vector<int>().swap(fMtNum5);
    std::vector<int>().swap(fMtNum6);
    //    std::cout<<" MAT number "<< tMat->GetMAT() <<"  "<< nmats << std::endl;
    TIter iter(tMat->GetFiles());
    TNudyEndfFile *file;
    while ((file = (TNudyEndfFile *)iter.Next())) {
      // std::cout <<" mf "<< file->GetMF() << std::endl;
      // Read File data into class structures
      switch (file->GetMF()) {
      case 1:
        if (fFlagRead != 1) {
          recoNuPh  = new TNudyEndfNuPh(file);
          fFlagRead = 1;
          // std::cout << "file 1 OK: Should be printed only one time " << std::endl;
        }
        break;
      case 2:
        break;
      case 3: {
        ReadFile3(file);
        TIter secIter(file->GetSections());
        TNudyEndfSec *sec;
        while ((sec = (TNudyEndfSec *)secIter.Next())) {
          int MT = sec->GetMT();
          if (MT == 1) {
            fEneUni.push_back(fEnergyMts);
            fSigUniT.push_back(fSigmaMts);
            break;
          }
        }
        FixupTotal(fEnergyMts);
        //	  double sigfis, sigcap;
        //	  std::vector<double> energ, etavalue;
        //	  std::vector<double> energjunk, etajunk;
        //	  eta->Branch("energy",&energ);
        //	  eta->Branch("eta",&etavalue);
        //	  std::cout<<"hello rrebin  "<<fEnergyMts.size()<<std::endl;
        //	    for (unsigned int j1 = 1; j1 < fEnergyMts.size() ; j1++){
        //	       for(unsigned int j2 = 0; j2 < MtValues[0].size(); j2++){
        //	          if(MtValues[0][j2] ==18)sigfis = GetSigmaPartial(0, j2, fEnergyMts[j1]);
        //	          if(MtValues[0][j2] ==102)sigcap = GetSigmaPartial(0, j2, fEnergyMts[j1]);
        //	       }
        // 	       if (recoNuPh->GetNuTotal(0, fEnergyMts[j1])/(1+sigcap/sigfis)>0){
        // 		h1->Fill(log10(fEnergyMts[j1]),recoNuPh->GetNuTotal(0, fEnergyMts[j1])/(1+sigcap/sigfis));
        // 		h2->Fill(log10(fEnergyMts[j1]));
        // 	       }
        // 	       energ.push_back(fEnergyMts[j1]);
        // 	       etavalue.push_back(recoNuPh->GetNuTotal(0, fEnergyMts[j1])/(1+sigcap/sigfis));
        //	      std::cout<< fEnergyMts[j1] <<"  "<< recoNuPh->GetNuTotal(0, fEnergyMts[j1])/(1+sigcap/sigfis) <<
        // std::endl;
        //    }
        // 	    for(int j = 0; j < 150; j++){
        // 	      double num = h1->GetBinContent(j);
        // 	      int    num2 = h2->GetBinContent(j);
        // 	      if ( num >0 && num2 > 0)h3->SetBinContent(j,num/num2);
        // 	      std::cout<<j<<"  "<<num<<"  "<<num2<<"  "<<num/num2 << std::endl;
        // 	    }
        // 	  std::cout<<"hello eta "<<std::endl;
        // 	       eta->Fill();
        // 	       f->Write();
        fEnergyMts.clear();
        fSigmaMts.clear();
      } break;
      case 4:
        // std::cout << "before file 4 " << std::endl;
        recoAng = new TNudyEndfAng(file);
        // std::cout << "file 4 OK " << std::endl;
        break;
      case 5:
        // std::cout << "before file 5 " << std::endl;
        recoEnergy = new TNudyEndfEnergy(file);
        // std::cout << "file 5 OK " << std::endl;
        break;
      case 6:
        // std::cout << "before file 6 " << std::endl;
        recoEnergyAng = new TNudyEndfEnergyAng(file, fQValue);
        // std::cout << "file 6 OK " << std::endl;
        break;
      case 8:
        // std::cout << "before file 8 " << std::endl;
        recoFissY = new TNudyEndfFissionYield(file);
        // std::cout << "file 8 OK " << std::endl;
        break;
      // //       case 12:
      // //         std::cout << "before file 12 " << std::endl;
      // // 	recoPhYield = new TNudyEndfPhYield(file);
      // // 	std::cout<<"file 12 OK "<<std::endl;
      // // 	break;
      // //       case 13:
      // //         std::cout << "before file 13 " << std::endl;
      // // 	recoPhProd = new TNudyEndfPhProd(file);
      // // 	std::cout<<"file 13 OK "<<std::endl;
      // // 	break;
      case 14:
        // std::cout << "before file 14 " << std::endl;
        recoPhAng = new TNudyEndfPhAng(file);
        // std::cout<<"file 14 OK "<<std::endl;
        break;
        // //       case 15:
        // //         std::cout << "before file 15 " << std::endl;
        // // 	recoPhEnergy = new TNudyEndfPhEnergy(file);
        // // 	std::cout<<"file 15 OK "<<std::endl;
        // // 	break;
      }
    }
  }
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::FixupTotal(std::vector<double> &x1)
{
  for (unsigned long i = 0; i < fSigmaOfMts.size(); i++) {
    int size    = fSigmaOfMts[i].size() / 2;
    int flagLoc = -1;
    for (unsigned long k = 0; k < x1.size(); k++) {
      int min = 0;
      int max = size - 1;
      int mid = 0;
      if (x1[k] <= fSigmaOfMts[i][min])
        min = 0;
      else if (x1[k] >= fSigmaOfMts[i][max])
        min = max;
      else {
        while (max - min > 1) {
          mid = (min + max) / 2;
          if (x1[k] < fSigmaOfMts[i][mid])
            max = mid;
          else
            min = mid;
        }
      }
      if (x1[k] == fSigmaOfMts[i][min] && fSigmaOfMts[i][size + min] > 0) {
        fEneTemp.push_back(x1[k]);
        fSigTemp.push_back(fSigmaOfMts[i][size + min]);
      }
      if (fEneTemp.size() == 1 && flagLoc == -1) {
        fEnergyLocationMts.push_back(k);
        flagLoc = 0;
      }
    }
    fSigmaUniOfMts.push_back(fSigTemp);
    fEneTemp.clear();
    fSigTemp.clear();
  }
  fSigmaOfMts.clear();
  fSigUniOfMt.push_back(fSigmaUniOfMts);
  fEnergyLocMtId.push_back(fEnergyLocationMts);
  fSigmaUniOfMts.clear();
  fEnergyLocationMts.clear();
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetSigmaTotal(int ielemId, double energyK)
{
  int min = 0;
  int max = fEneUni[ielemId].size() - 1;
  int mid = 0;
  if (energyK <= fEneUni[ielemId][min]) {
    min = 0;
  } else if (energyK >= fEneUni[ielemId][max]) {
    min = max - 1;
  } else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEneUni[ielemId][mid]) {
        max = mid;
      } else {
        min = mid;
      }
    }
  }
  return fSigUniT[ielemId][min] +
         (fSigUniT[ielemId][min + 1] - fSigUniT[ielemId][min]) * (energyK - fEneUni[ielemId][min]) /
             (fEneUni[ielemId][min + 1] - fEneUni[ielemId][min]);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetSigmaPartial(int ielemId, int i, double energyK)
{
  // std::cout<< fEneUni[ielemId].size() <<std::endl;
  int min = 0;
  int max = fEneUni[ielemId].size() - 1;
  int mid = 0;
  if (energyK <= fEneUni[ielemId][min])
    min = 0;
  else if (energyK >= fEneUni[ielemId][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEneUni[ielemId][mid])
        max = mid;
      else
        min = mid;
    }
  }
  int minp = min - fEnergyLocMtId[ielemId][i];
  //     for (unsigned int j1 = 0; j1 < fEneUni[ielemId].size(); j1++)
  //       std::cout << fEneUni[ielemId][j1] <<"  "<< fSigUniOfMt[ielemId][i][j1] <<"  "<<
  //       fSigUniOfMt[ielemId][i].size() << std::endl;
  if (minp <= 0) return 0;
  if (minp >= (int)fSigUniOfMt[ielemId][i].size()) return 0;
  return fSigUniOfMt[ielemId][i][minp] +
         (fSigUniOfMt[ielemId][i][minp + 1] - fSigUniOfMt[ielemId][i][minp]) * (energyK - fEneUni[ielemId][min]) /
             (fEneUni[ielemId][min + 1] - fEneUni[ielemId][min]);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetCos6(int ielemId, int mt, double energyK)
{
  return recoEnergyAng->GetCos6(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetEnergy6(int ielemId, int mt, double energyK)
{
  return recoEnergyAng->GetEnergy6(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetCos64(int ielemId, int mt, double energyK)
{
  return recoEnergyAng->GetCos64(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetCos4(int ielemId, int mt, double energyK)
{
  return recoAng->GetCos4(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetEnergy5(int ielemId, int mt, double energyK)
{
  return recoEnergy->GetEnergy5(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetAd6(int ielemId, int mt)
{
  return recoEnergyAng->GetAd6(ielemId, mt);
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetZd6(int ielemId, int mt)
{
  return recoEnergyAng->GetZd6(ielemId, mt);
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetLaw6(int ielemId, int mt)
{
  return recoEnergyAng->GetLaw6(ielemId, mt);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetQValue(int ielemId, int mt)
{
  int size = fQvalue[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (fQvalue[ielemId][2 * i + 1] == mt) {
      return fQvalue[ielemId][2 * i];
    }
  }
  return 0;
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetMt4(int ielemId, int mt)
{
  if (fMt4.size() <= 0) return 99;
  int size = fMt4[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (fMt4[ielemId][2 * i + 1] == mt) {
      return fMt4[ielemId][2 * i];
    }
  }
  return 99;
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetMt5(int ielemId, int mt)
{
  if (fMt5.size() <= 0) return -1;
  int size = fMt5[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (fMt5[ielemId][2 * i + 1] == mt) {
      return fMt5[ielemId][2 * i];
    }
  }
  return -1;
}
//-------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetMt6Neutron(int ielemId, int mt)
{
  return recoEnergyAng->GetMt6Neutron(ielemId, mt);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetMt6(int ielemId, int mt)
{
  if (fMt6.size() <= 0) return 99;
  int size = fMt6[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (fMt6[ielemId][2 * i + 1] == mt) {
      return fMt6[ielemId][2 * i];
    }
  }
  return 99;
}
//-------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetCos4Lct(int ielemId, int mt)
{
  return recoAng->GetCos4Lct(ielemId, mt);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetNuTotal(int elemid, double energyK)
{
  return recoNuPh->GetNuTotal(elemid, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetNuPrompt(int elemid, double energyK)
{
  return recoNuPh->GetNuPrompt(elemid, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetNuDelayed(int elemid, double energyK)
{
  return recoNuPh->GetNuDelayed(elemid, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetFissHeat(int elemid, double energyK)
{
  return recoNuPh->GetFissHeat(elemid, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetFisYield(int elemid, double energyK)
{
  return recoFissY->GetFisYield(elemid, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetLambdaD(int ielemId, int time)
{
  return recoNuPh->GetLambdaD(ielemId, time);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetDelayedFraction(int ielemId, int mt, double energyK)
{
  return recoEnergy->GetDelayedFraction(ielemId, mt, energyK);
}
//-------------------------------------------------------------------------------------------------------
TNudyEndfRecoPoint::~TNudyEndfRecoPoint()
{
  fSigmaOfMts.shrink_to_fit();
  fSigmaUniOfMts.shrink_to_fit();
  fEnergyLocationMts.shrink_to_fit();
  fMtNumbers.shrink_to_fit();
  fMtNum4.shrink_to_fit();
  fMtNum5.shrink_to_fit();
  fMtNum6.shrink_to_fit();
  fSigmaMts.shrink_to_fit();
  fQvalueTemp.shrink_to_fit();
  fELinearFile3.shrink_to_fit();
  fXLinearFile3.shrink_to_fit();
  fEneTemp.shrink_to_fit();
  fSigTemp.shrink_to_fit();
}
