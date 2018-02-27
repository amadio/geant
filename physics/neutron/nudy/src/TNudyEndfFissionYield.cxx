// 	This class is reconstructing Fission yields
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 24, 2016

#include "TList.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfFissionYield.h"

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
        ein.push_back(list1->GetC1());
        // int NN   = list1->GetN1();
        int NFP = list1->GetN2();
        // std::cout<<"energy i " <<ein[i] << std::endl;
        for (int j = 0; j < NFP; j++) {
          zafp1.push_back(list1->GetLIST(4 * j + 0));
          fps1.push_back(list1->GetLIST(4 * j + 1));
          yi1.push_back(list1->GetLIST(4 * j + 2));
          dyi1.push_back(list1->GetLIST(4 * j + 3));
          // divr = div(zafp1[j],1000);
          // std::cout<< divr.rem <<"  "<<  yi1[j] << std::endl;
        }
        TNudyCore::Instance()->cdfGenerateT(zafp1, yi1, cyi1);
        zafp.push_back(zafp1);
        fps.push_back(fps1);
        yi.push_back(yi1);
        cyi.push_back(cyi1);
        dyi.push_back(dyi1);
        zafp1.clear();
        fps1.clear();
        yi1.clear();
        cyi1.clear();
        dyi1.clear();
      }
      einfId.push_back(ein);
      zafId.push_back(zafp);
      pdfYieldId.push_back(yi);
      cdfYieldId.push_back(cyi);
      yi.clear();
      cyi.clear();
      zafp.clear();
      ein.clear();
    } else if (MT == 459) { // Neutron induced cummulative fission yield
      for (int i = 0; i < LE; i++) {
        TNudyEndfList *list1 = (TNudyEndfList *)recIter.Next();
        einc.push_back(list1->GetC1());
        // int NN   = list1->GetN1();
        int NFP = list1->GetN2();
        //          std::cout<<"energy c " <<einc[i] << std::endl;
        for (int j = 0; j < NFP; j++) {
          zafpc1.push_back(list1->GetLIST(4 * j + 0));
          fpsc1.push_back(list1->GetLIST(4 * j + 1));
          yc1.push_back(list1->GetLIST(4 * j + 2));
          dyc1.push_back(list1->GetLIST(4 * j + 3));
          //	    std::cout<<"fission yield c "<< zafpc[i].At(j) << std::endl;
        }
        zafpc.push_back(zafpc1);
        fpsc.push_back(fpsc1);
        yc.push_back(yc1);
        dyc.push_back(dyc1);
        zafpc1.clear();
        fpsc1.clear();
        yc1.clear();
        dyc1.clear();
      }
    }
  }
}

TNudyEndfFissionYield::~TNudyEndfFissionYield()
{
  ein.shrink_to_fit();
  einc.shrink_to_fit();
  zafp.shrink_to_fit();
  fps.shrink_to_fit();
  zafpc.shrink_to_fit();
  fpsc.shrink_to_fit();
  yi.shrink_to_fit();
  cyi.shrink_to_fit();
  dyi.shrink_to_fit();
  yc.shrink_to_fit();
  dyc.shrink_to_fit();
  zafp1.shrink_to_fit();
  fps1.shrink_to_fit();
  zafpc1.shrink_to_fit();
  fpsc1.shrink_to_fit();
  yi1.shrink_to_fit();
  cyi1.shrink_to_fit();
  dyi1.shrink_to_fit();
  yc1.shrink_to_fit();
  dyc1.shrink_to_fit();
}

double TNudyEndfFissionYield::GetFisYield(int ielemId, double energyK)
{
  fRnd = new TRandom3(0);
  // std::cout<<"element "<< einfId.size() << std::endl;
  // std::cout<<"energies "<< einfId[ielemId].size() << std::endl;
  int min = 0;
  int max = einfId[ielemId].size() - 1;
  int mid = 0;
  if (energyK <= einfId[ielemId][min])
    min = 0;
  else if (energyK >= einfId[ielemId][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < einfId[ielemId][mid])
        max = mid;
      else
        min = mid;
    }
  }
  double fraction          = (energyK - einfId[ielemId][min]) / (einfId[ielemId][min + 1] - einfId[ielemId][min]);
  double rnd1              = fRnd->Uniform(1);
  double rnd2              = fRnd->Uniform(1);
  if (rnd2 < fraction) min = min + 1;
  int k                    = 0;
  int size                 = pdfYieldId[ielemId][min].size();
  for (int j = 1; j < size; j++) {
    if (rnd1 < cdfYieldId[ielemId][min][j]) {
      k                    = j - 1;
      if (k >= size - 1) k = size - 1;
      break;
    }
  }
  double plk = (pdfYieldId[ielemId][min][k + 1] - pdfYieldId[ielemId][min][k]) /
               (zafId[ielemId][min][k + 1] - zafId[ielemId][min][k]);
  double plk2 = pdfYieldId[ielemId][min][k] * pdfYieldId[ielemId][min][k];
  double plsq = plk2 + 2 * plk * (rnd1 - cdfYieldId[ielemId][min][k]);
  // std::cout <<"plk "<< plk <<" plk2 "<< plk2 <<"  "<< k << std::endl;
  double zaf                    = 0;
  if (plk != 0 && plsq > 0) zaf = zafId[ielemId][min][k] + (sqrt(std::fabs(plsq)) - pdfYieldId[ielemId][min][k]) / plk;
  return zaf;
}
