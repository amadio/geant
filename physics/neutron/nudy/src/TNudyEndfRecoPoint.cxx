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
#include "Geant/TNudyEndfRecoPoint.h"

#ifdef USE_ROOT
ClassImp(TNudyEndfRecoPoint)
#endif

#ifdef USE_ROOT
#include "TRandom3.h"
#endif

    TNudyEndfRecoPoint::TNudyEndfRecoPoint()
    : elemId(0), rENDF()
{
}
TNudyEndfRecoPoint::TNudyEndfRecoPoint(int ielemId, const char *irENDF)
    : elemId(ielemId), rENDF(irENDF)
{
  GetData(ielemId, irENDF);
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ReadFile3(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    // std::cout <<"MT  "<<sec->GetMT() <<std::endl;
    int MT = sec->GetMT();
    if (MT != 1 && MT != 3 && MT != 4 && MT != 27 && MT != 19 && MT != 20 && MT != 21 && MT != 38 && MT != 101 &&
        MT < 250) {
      MtNumbers.push_back(MT);
      TIter recIter(sec->GetRecords());
      TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
      double QI            = header->GetC2();
      QValue[sec->GetMT()] = QI;
      qvaluetemp.push_back(QI);
      qvaluetemp.push_back(MT);
      NR = header->GetN1();
      NP = header->GetN2();
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
      for (int crs = 0; crs < NP; crs++) {
        eLinearFile3.push_back(tab1->GetX(crs));
        xLinearFile3.push_back(tab1->GetY(crs));
      }
      if (MT==1) {
         energyMts.insert(std::end(energyMts),std::begin(eLinearFile3),std::end(eLinearFile3));
         sigmaMts.insert(std::end(sigmaMts),std::begin(xLinearFile3),std::end(xLinearFile3));
      } else {
         sigmaOfMts.push_back (eLinearFile3);
         sigmaOfMts.push_back (xLinearFile3);
      }
      eLinearFile3.clear();
      xLinearFile3.clear();
    }
  }
  MtValues.push_back(MtNumbers);
  qvalue.push_back(qvaluetemp);
  qvaluetemp.clear();
  MtNumbers.clear();
}
void TNudyEndfRecoPoint::GetData(int ielemId, const char *rENDF)
{
  elemId      = ielemId;
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
        MtNum4.push_back(tMat->GetMFn(inf));
        MtNum4.push_back(tMat->GetMTn(inf));
        break;
      case 5:
        MtNum5.push_back(tMat->GetMFn(inf));
        if (tMat->GetMTn(inf) == 455) {
          MtNum5.push_back(mt455);
          mt455++;
        } else {
          MtNum5.push_back(tMat->GetMTn(inf));
        }
        break;
      case 6:
        MtNum6.push_back(tMat->GetMFn(inf));
        MtNum6.push_back(tMat->GetMTn(inf));
        break;
      }
    }
    Mt4.push_back(MtNum4);
    Mt5.push_back(MtNum5);
    Mt6.push_back(MtNum6);
    std::vector<int>().swap(MtNum4);
    std::vector<int>().swap(MtNum5);
    std::vector<int>().swap(MtNum6);
    //    std::cout<<" MAT number "<< tMat->GetMAT() <<"  "<< nmats << std::endl;
    TIter iter(tMat->GetFiles());
    TNudyEndfFile *file;
    while ((file = (TNudyEndfFile *)iter.Next())) {
      // std::cout <<" mf "<< file->GetMF() << std::endl;
      // Read File data into class structures
      switch (file->GetMF()) {
      case 1:
        if (flagRead != 1) {
          recoNuPh = new TNudyEndfNuPh(file);
          flagRead = 1;
          std::cout << "file 1 OK: Should be printed only one time " << std::endl;
        }
        break;
      case 2:
        break;
      case 3:
        { 
          ReadFile3(file);
          TIter secIter(file->GetSections());
          TNudyEndfSec *sec;
          while ((sec = (TNudyEndfSec *)secIter.Next())) {
    //    std::cout <<"MT  "<<sec->GetMT() <<std::endl;
             int MT = sec->GetMT();
             if (MT != 3 && MT != 4 && MT != 27 && MT != 19 && MT != 20 && MT != 21 && MT != 38 && MT != 101 && MT < 250) {
                if (MT == 1) {
                   eneUni.push_back (energyMts);
                   sigUniT.push_back (sigmaMts);
                   continue;
                }
                fixupTotal(energyMts);
             }
          }
                 energyMts.clear();
                 sigmaMts.clear();
        }break;
      case 4:
        recoAng = new TNudyEndfAng(file);
        std::cout << "file 4 OK " << std::endl;
        break;
      case 5:
        recoEnergy = new TNudyEndfEnergy(file);
        std::cout << "file 5 OK " << std::endl;
        break;
      case 6:
        recoEnergyAng = new TNudyEndfEnergyAng(file, QValue);
        std::cout << "file 6 OK " << std::endl;
        break;
      case 8:
        recoFissY = new TNudyEndfFissionYield(file);
        std::cout << "file 8 OK " << std::endl;
        break;
      }
    }
  }
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::fixupTotal(std::vector<double> &x1)
{

  for (unsigned long i = 0; i < sigmaOfMts.size(); i++) {
    int size = sigmaOfMts[i].size()/ 2;
    for (unsigned long k = 0; k < x1.size(); k++) {
      int min = 0;
      int max = size - 1;
      int mid = 0;
      if (x1[k] <= sigmaOfMts[i][min])
        min = 0;
      else if (x1[k] >= sigmaOfMts[i][max])
        min = max - 1;
      else {
        while (max - min > 1) {
          mid = (min + max) / 2;
          if (x1[k] < sigmaOfMts[i][mid])
            max = mid;
          else
            min = mid;
        }
      }
      if (x1[k] == sigmaOfMts[i][min] && sigmaOfMts[i][size + min] > 1E-20) {
        eneTemp.push_back(x1[k]);
        sigTemp.push_back(sigmaOfMts[i][size + min]);
      }
      if (eneTemp.size() == 1) {
        energyLocationMts.push_back(k);
      }
    }
    sigmaUniOfMts.push_back(sigTemp);
    eneTemp.clear();
    sigTemp.clear();
  }
  sigmaOfMts.clear();
  sigUniOfMt.push_back(sigmaUniOfMts);
  energyLocMtId.push_back(energyLocationMts);
  sigmaUniOfMts.clear();
  energyLocationMts.clear();
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetSigmaTotal(int ielemId, double energyK)
{
  int min = 0;
  int max = eneUni[ielemId].size() - 1;
  int mid = 0;
  if (energyK <= eneUni[ielemId][min])
    min = 0;
  else if (energyK >= eneUni[ielemId][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < eneUni[ielemId][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return sigUniT[ielemId][min] +
         (sigUniT[ielemId][min + 1] - sigUniT[ielemId][min]) * (energyK - eneUni[ielemId][min]) /
             (eneUni[ielemId][min + 1] - eneUni[ielemId][min]);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetSigmaPartial(int ielemId, int i, double energyK)
{
  int min = 0;
  int max = eneUni[ielemId].size() - 1;
  int mid = 0;
  if (energyK <= eneUni[ielemId][min])
    min = 0;
  else if (energyK >= eneUni[ielemId][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < eneUni[ielemId][mid])
        max = mid;
      else
        min = mid;
    }
  }
  int minp = min - energyLocMtId[ielemId][i];
  // for (unsigned int j1 = 0; j1 < eneUni[ielemId].size(); j1++)
  // std::cout << eneUni[ielemId][j1] <<"  "<<  sigUniOfMt[ielemId][i][minp] << std::endl;
  if (minp < 0) return 0;
  if (minp >= (int)sigUniOfMt[ielemId][i].size()) return 0;
  return sigUniOfMt[ielemId][i][minp] +
         (sigUniOfMt[ielemId][i][minp + 1] - sigUniOfMt[ielemId][i][minp]) * (energyK - eneUni[ielemId][min]) /
             (eneUni[ielemId][min + 1] - eneUni[ielemId][min]);
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
  int size = qvalue[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (qvalue[ielemId][2 * i + 1] == mt) {
      return qvalue[ielemId][2 * i];
    }
  }
  return 0;
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetMt4(int ielemId, int mt)
{
  if (Mt4.size() <= 0) return 99;
  int size = Mt4[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (Mt4[ielemId][2 * i + 1] == mt) {
      return Mt4[ielemId][2 * i];
    }
  }
  return 99;
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetMt5(int ielemId, int mt)
{
  if (Mt5.size() <= 0) return -1;
  int size = Mt5[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (Mt5[ielemId][2 * i + 1] == mt) {
      return Mt5[ielemId][2 * i];
    }
  }
  return -1;
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetMt6(int ielemId, int mt)
{
  if (Mt6.size() <= 0) return 99;
  int size = Mt6[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (Mt6[ielemId][2 * i + 1] == mt) {
      return Mt6[ielemId][2 * i];
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
  sigmaOfMts.shrink_to_fit();
  sigmaUniOfMts.shrink_to_fit();
  energyLocationMts.shrink_to_fit();
  MtNumbers.shrink_to_fit();
  MtNum4.shrink_to_fit();
  MtNum5.shrink_to_fit();
  MtNum6.shrink_to_fit();
  sigmaMts.shrink_to_fit();
  qvaluetemp.shrink_to_fit();
  eLinearFile3.shrink_to_fit();
  xLinearFile3.shrink_to_fit();
  eneTemp.shrink_to_fit();
  sigTemp.shrink_to_fit();
}
