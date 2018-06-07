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
#endif

TNudyEndfFissionYield::TNudyEndfFissionYield()
{
}

//------------------------------------------------------------------------------------------------------
TNudyEndfFissionYield::TNudyEndfFissionYield(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    int MT = sec->GetMT();
    int LE = sec->GetL1();
    if (MT == 454) { // Neutron induced independent fission yield
      for (int i = 0; i < LE; i++) {
        TNudyEndfList *list1 = (TNudyEndfList *)recIter.Next();
        fEin.push_back(list1->GetC1());
        int NFP = list1->GetN2();
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
        }
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
        int NFP = list1->GetN2();
        for (int j = 0; j < NFP; j++) {
          fZafpc1.push_back(10 * list1->GetLIST(4 * j + 0) + list1->GetLIST(4 * j + 1));
          fFpsc1.push_back(list1->GetLIST(4 * j + 1));
          fYc1.push_back(list1->GetLIST(4 * j + 2));
          fDyc1.push_back(list1->GetLIST(4 * j + 3));
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

