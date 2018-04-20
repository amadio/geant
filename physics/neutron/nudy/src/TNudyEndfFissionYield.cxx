// 	This class is reconstructing Fission yields
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 24, 2016

#include "TList.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfFissionYield.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfFissionYield)
#include "TRandom3.h"
#endif

    TNudyEndfFissionYield::TNudyEndfFissionYield()
{
}

//______________________________________________________________________________
TNudyEndfFissionYield::TNudyEndfFissionYield(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    // double ZA   = sec->GetC1();
    // double AWR  = sec->GetC2();
    // div_t divr;
    int MT = sec->GetMT();
    int LE = sec->GetL1();
    if (MT == 454) { // Neutron induced independent fission yield
      for (int i = 0; i < LE; i++) {
        TNudyEndfList *list1 = (TNudyEndfList *)recIter.Next();
        fEin.push_back(list1->GetC1());
        // int NN   = list1->GetN1();
        int NFP = list1->GetN2();
        // std::cout<<"energy i " <<fEin[i] << std::endl;
        double sum = 0;
        for (int j = 0; j < NFP; j++) {
          if (list1->GetLIST(4 * j + 2) > 0) {
            fZafp1.push_back(10 * list1->GetLIST(4 * j + 0) + list1->GetLIST(4 * j + 1));
            fFps1.push_back(list1->GetLIST(4 * j + 1));
            fYi1.push_back(list1->GetLIST(4 * j + 2) / 2);
            fDyi1.push_back(list1->GetLIST(4 * j + 3));
            sum += list1->GetLIST(4 * j + 2);
            fCyi1.push_back(sum / 2);
          }
          // divr = div(fZafp1[j],1000);
          // std::cout<< 10*list1->GetLIST(4 * j + 0) + list1->GetLIST(4 * j + 1) <<"  "<<  list1->GetLIST(4 * j + 2)/2
          // << std::endl;
        }
        //        TNudyCore::Instance()->cdfGenerateT(fZafp1, fYi1, fCyi1);
        // std::cout << "sum fission yield \t" << sum << std::endl;
        fZafp.push_back(fZafp1);
        fFps.push_back(fFps1);
        fYi.push_back(fYi1);
        fCyi.push_back(fCyi1);
        fDyi.push_back(fDyi1);
        fZafp1.clear();
        fFps1.clear();
        fYi1.clear();
        fCyi1.clear();
        fDyi1.clear();
      }
      fEinfId.push_back(fEin);
      fZafId.push_back(fZafp);
      fPdfYieldId.push_back(fYi);
      fCdfYieldId.push_back(fCyi);
      fYi.clear();
      fCyi.clear();
      fZafp.clear();
      fEin.clear();
    } else if (MT == 459) { // Neutron induced cummulative fission yield
      for (int i = 0; i < LE; i++) {
        TNudyEndfList *list1 = (TNudyEndfList *)recIter.Next();
        fEinc.push_back(list1->GetC1());
        // int NN   = list1->GetN1();
        int NFP = list1->GetN2();
        //          std::cout<<"energy c " <<fEinc[i] << std::endl;
        for (int j = 0; j < NFP; j++) {
          fZafpc1.push_back(10 * list1->GetLIST(4 * j + 0) + list1->GetLIST(4 * j + 1));
          fFpsc1.push_back(list1->GetLIST(4 * j + 1));
          fYc1.push_back(list1->GetLIST(4 * j + 2));
          fDyc1.push_back(list1->GetLIST(4 * j + 3));
          //	    std::cout<<"fission yield c "<< fZafpc[i].At(j) << std::endl;
        }
        fZafpc.push_back(fZafpc1);
        fFpsc.push_back(fFpsc1);
        fYc.push_back(fYc1);
        fDyc.push_back(fDyc1);
        fZafpc1.clear();
        fFpsc1.clear();
        fYc1.clear();
        fDyc1.clear();
      }
    }
  }
}

TNudyEndfFissionYield::~TNudyEndfFissionYield()
{
  fEin.shrink_to_fit();
  fEinc.shrink_to_fit();
  fZafp.shrink_to_fit();
  fFps.shrink_to_fit();
  fZafpc.shrink_to_fit();
  fFpsc.shrink_to_fit();
  fYi.shrink_to_fit();
  fCyi.shrink_to_fit();
  fDyi.shrink_to_fit();
  fYc.shrink_to_fit();
  fDyc.shrink_to_fit();
  fZafp1.shrink_to_fit();
  fFps1.shrink_to_fit();
  fZafpc1.shrink_to_fit();
  fFpsc1.shrink_to_fit();
  fYi1.shrink_to_fit();
  fCyi1.shrink_to_fit();
  fDyi1.shrink_to_fit();
  fYc1.shrink_to_fit();
  fDyc1.shrink_to_fit();
}

double TNudyEndfFissionYield::GetFisYield(int ielemId, double energyK)
{
  fRnd = new TRandom3(0);
  // std::cout<<"element "<< fEinfId.size() << std::endl;
  // std::cout<<"energies "<< fEinfId[ielemId].size() << std::endl;
  int min = 0;
  int max = fEinfId[ielemId].size() - 1;
  int mid = 0;
  if (energyK <= fEinfId[ielemId][min])
    min = 0;
  else if (energyK >= fEinfId[ielemId][max])
    min = max;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEinfId[ielemId][mid])
        max = mid;
      else
        min = mid;
    }
  }
  double fraction          = (energyK - fEinfId[ielemId][min]) / (fEinfId[ielemId][min + 1] - fEinfId[ielemId][min]);
  double rnd1              = fRnd->Uniform(1);
  double rnd2              = fRnd->Uniform(1);
  if (rnd2 < fraction) min = min + 1;
  int k                    = 0;
  int size                 = fPdfYieldId[ielemId][min].size();
  for (int j = 1; j < size; j++) {
    if (rnd1 <= fCdfYieldId[ielemId][min][j]) {
      k                    = j - 1;
      if (k >= size - 1) k = size - 1;
      break;
    }
  }
  double plk = (fPdfYieldId[ielemId][min][k + 1] - fPdfYieldId[ielemId][min][k]) /
               (fZafId[ielemId][min][k + 1] - fZafId[ielemId][min][k]);
  double plk2                      = fPdfYieldId[ielemId][min][k] * fPdfYieldId[ielemId][min][k];
  double plsq                      = plk2 + 2 * plk * (rnd1 - fCdfYieldId[ielemId][min][k]);
  double zaf                       = 0;
  if (plk == 0 && rnd1 < 0.5) zaf  = fZafId[ielemId][min][k];
  if (plk == 0 && rnd1 >= 0.5) zaf = fZafId[ielemId][min][k + 1];
  if (plk != 0 && plsq > 0)
    zaf = fZafId[ielemId][min][k] + (sqrt(std::fabs(plsq)) - fPdfYieldId[ielemId][min][k]) / plk;
  return zaf;
}
