#include "TList.h"
#include "TFile.h"
#include "TKey.h"
#include "Geant/TNudyENDF.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfTape.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyEndfMat.h"
#include "Geant/TNudyEndfDoppler.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfSigma.h"
#include "Math/SpecFuncMathMore.h"
#include "Geant/TNudyEndfEnergy.h"
#include "Geant/TNudyEndfEnergyAng.h"
#include "Geant/TNudyEndfFissionYield.h"
#include "Geant/TNudyEndfPhYield.h"
#include "Geant/TNudyEndfPhProd.h"
#include "Geant/TNudyEndfPhAng.h"
#include "Geant/TNudyEndfPhEnergy.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfSigma)
#endif

#ifdef USE_ROOT
#include "TRandom3.h"
#endif

#include <algorithm>

    TNudyEndfSigma::TNudyEndfSigma()
    : rENDF(), fSigDiff(0)
{
}

TNudyEndfSigma::TNudyEndfSigma(const char *irENDF, double isigDiff) : rENDF(irENDF), fSigDiff(isigDiff)
{
  // GetData(irENDF, isigDiff);
}

//____________________________________________________________________________________________________________________
double TNudyEndfSigma::RecursionLinearNuPh(double x1, double x2, double sig1, double sig2, std::vector<double> x,
                                           std::vector<double> sig)
{
  double siga;
  double mid = 0.5 * (x1 + x2);
  if ((sig1 == 0.0 && sig2 == 0.0) || x1 == x2 || x1 < 1E-5 || x2 < 1E-5) return 0;
  siga           = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, x, sig, fNP, mid);
  double sigmid1 = sig1 + (sig2 - sig1) * (mid - x1) / (x2 - x1);
  if (fabs((siga - sigmid1) / sigmid1) <= fSigDiff) {
    return 0;
  }
  x.push_back(mid);
  sig.push_back(siga);
  RecursionLinearNuPh(x1, mid, sig1, siga, x, sig);
  RecursionLinearNuPh(mid, x2, siga, sig2, x, sig);
  return 0;
}
//______________________________________________________________________________
void TNudyEndfSigma::ReadFile1(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    int fMT = sec->GetMT();
    if (fMT == 452) {
      // Total fission neutron multiplicity polynomial expansion
      int LNU = sec->GetL2();
      if (LNU == 1) {
        TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
        for (int j = 0, eN1 = list->GetN1(); j != eN1; ++j) {
          fCnc.push_back(list->GetLIST(j));
        }
        double ein = 1E-5;
        do {
          double nun = 0;
          for (int i = 0, eN1 = list->GetN1(); i != eN1; ++i) {
            nun += fCnc[i] * pow(ein, i);
          }
          fEintFile1.push_back(ein);
          fNutFile1.push_back(nun);
          ein *= 2;
        } while (ein < 21E8);
        fCnc.clear();
      } else {
        // Total fission neutron multiplicity tabulated representation
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        fNR                 = tab1->GetN1();
        fNP                 = tab1->GetN2();
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fEintFile1.push_back(tab1->GetX(crs));
          fNutFile1.push_back(tab1->GetY(crs));
        }
        for (int cr = 0; cr < fNP - 1; cr++) {
          RecursionLinearNuPh(fEintFile1[cr], fEintFile1[cr + 1], fNutFile1[cr], fNutFile1[cr + 1], fEintFile1,
                              fNutFile1);
        }
        TNudyCore::Instance()->Sort(fEintFile1, fNutFile1);
      }
      fEint.push_back(fEintFile1);
      fNut.push_back(fNutFile1);
      fEintFile1.clear();
      fNutFile1.clear();
      fNbt1.clear();
      fInt1.clear();
    } else if (fMT == 455) {
      // delayed neutron multiplicity
      int LDG = sec->GetL1();
      int LNU = sec->GetL2();
      if (LNU == 1 && LDG == 0) {

      } else if (LNU == 1 && LDG == 1) {

      } else if (LNU == 2 && LDG == 0) {
        TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
        int NNF             = list->GetNPL();
        for (int i = 0; i < NNF; i++) {
          fNui.push_back(list->GetLIST(i));
        }
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        fNR                 = tab1->GetN1();
        fNP                 = tab1->GetN2();
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fEindFile1.push_back(tab1->GetX(crs));
          fNudFile1.push_back(tab1->GetY(crs));
        }
        for (int cr = 0; cr < fNP - 1; cr++) {
          RecursionLinearNuPh(fEindFile1[cr], fEindFile1[cr + 1], fNudFile1[cr], fNudFile1[cr + 1], fEindFile1,
                              fNudFile1);
        }
        TNudyCore::Instance()->Sort(fEindFile1, fNudFile1);
        fEind.push_back(fEindFile1);
        fNud.push_back(fNudFile1);
        fLambdaD.push_back(fNui);
        fEindFile1.clear();
        fNudFile1.clear();
        fNui.clear();
        fNbt1.clear();
        fInt1.clear();
      } else if (LNU == 2 && LDG == 1) {
      }
    } else if (fMT == 456) {
      // prompt neutron multiplicity
      int LNU = sec->GetL2();
      if (LNU == 1) {
      } else {
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        fNR                 = tab1->GetN1();
        fNP                 = tab1->GetN2();
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fEinFile1.push_back(tab1->GetX(crs));
          fNuFile1.push_back(tab1->GetY(crs));
        }
        for (int cr = 0; cr < fNP - 1; cr++) {
          RecursionLinearNuPh(fEinFile1[cr], fEinFile1[cr + 1], fNuFile1[cr], fNuFile1[cr + 1], fEinFile1, fNuFile1);
        }
        TNudyCore::Instance()->Sort(fEinFile1, fNuFile1);
        fNbt1.clear();
        fInt1.clear();
      }
      fEinp.push_back(fEinFile1);
      fNup.push_back(fNuFile1);
      fEinFile1.clear();
      fNuFile1.clear();
    } else if (fMT == 458) {
      TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
      int NPLY            = list->GetL2();
      if (NPLY == 0) {
        double EFR = list->GetLIST(0);
        double ENP = list->GetLIST(2);
        double END = list->GetLIST(4);
        double EGP = list->GetLIST(6);
        double EGD = list->GetLIST(8);
        double EB  = list->GetLIST(10);
        double ENU = list->GetLIST(12);
        double ein = 1E-5;
        do {
          double EFis = 0;
          EFis        = EFR + ENP + END + EGP + EGD + EB + ENU;
          EFis -= 0.100 * ein;
          // nuetrino energy dependence
          EFis -= 0.075 * ein;
          // delayed gamma energy dependence
          EFis -= 0.075 * ein;
          // delayed beta energy dependence
          int max = fEintFile1.size() - 1;
          int n0  = 0;
          max     = fEintFile1.size() - 1;
          int mid = 0;
          if (ein <= fEintFile1[n0]) {
            n0 = 0;
          } else if (ein >= fEintFile1[max]) {
            n0 = max - 1;
          } else {
            while (max - n0 > 1) {
              mid = (n0 + max) / 2;
              if (ein < fEintFile1[mid])
                max = mid;
              else
                n0 = mid;
              //      std::cout<<"n0 "<< n0 <<" max "<< max << std::endl;
            }
          }
          double nue = TNudyCore::Instance()->LinearInterpolation(fEintFile1[n0], fNutFile1[n0], fEintFile1[n0 + 1],
                                                                  fNutFile1[n0 + 1], ein);
          double nu0 = fNutFile1[0];
          EFis -= -1.307 * ein + 8.07 * 1E6 * (nue - nu0); // prompt neutron energy dependence
          fEinfFile1.push_back(ein);
          fHeatFile1.push_back(EFis);
          ein *= 2;
        } while (ein < 21E8);
      } else {
        int nply1 = 9 * (NPLY + 1);
        double c0[nply1], c1[nply1];
        for (int i = 0; i < nply1; i++) {
          c0[i] = list->GetLIST(i);
        }
        for (int i = nply1; i < 2 * nply1; i++) {
          c1[i - nply1] = list->GetLIST(i);
        }
        double ein = 1E-5;
        do {
          double EFis = 0;
          for (int i = 0; i < 9 * NPLY / 2 - 2; i++) {
            EFis += c0[i * 2] + c1[i * 2] * ein;
          }
          fEinfFile1.push_back(ein);
          fHeatFile1.push_back(EFis);
          ein *= 2;
        } while (ein < 21E8);
      }
      fEinFissHeat.push_back(fEinfFile1);
      fFissHeat.push_back(fHeatFile1);
      fEinfFile1.clear();
      fHeatFile1.clear();
    } else if (fMT == 460) {
      int LO = sec->GetL1();
      int NG = sec->GetN1();
      if (LO == 1) {
        for (int ng = 0; ng < NG; ng++) {
          TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
          fNR                 = tab1->GetN1();
          fNP                 = tab1->GetN2();
          for (int cr = 0; cr < fNR; cr++) {
            fNbt1.push_back(tab1->GetNBT(cr));
            fInt1.push_back(tab1->GetINT(cr));
          }
          for (int crs = 0; crs < fNP; crs++) {
            fEinPhFile1.push_back(tab1->GetX(crs));
            fPhFile1.push_back(tab1->GetY(crs));
          }
          // linearzation is stopped due to discrete photons which are to be matched with file 12 later
          // for(int cr=0; cr < fNP - 1 ; cr ++){
          //  recursionLinearNuPh(fEinPhFile1[cr], fEinPhFile1[cr+1], fPhFile1[cr], fPhFile1[cr+1],fEinPhFile1,
          //  fPhFile1);
          //}
          // TNudyCore::Instance()->Sort(fEinPhFile1, fPhFile1);
          fEinPhFile1.clear();
          fPhFile1.clear();
          fNbt1.clear();
          fInt1.clear();
        }
      } else {
        TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
        int NNF             = list->GetN1();
        // double lambda[NNF];
        for (int i = 0; i < NNF; i++) {
          // lambda[i] = list->GetLIST(i);
          //	     std::cout <<"lambda  "<< lambda[i] << std::endl;
        }
      }
    }
  }
}
//______________________________________________________________________________
void TNudyEndfSigma::ReadFile2(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  NLS           = 0;
  double gjdeno = 1;
  LSSF          = 0;
  int nrsl0     = 0;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    for (int k = 0, eN1 = sec->GetN1(); k != eN1; ++k) {
      TIter recIter(sec->GetRecords());
      ZAI                  = sec->GetC1();
      AWRI                 = sec->GetC2();
      NIS                  = sec->GetN1();
      TNudyEndfCont *cont1 = (TNudyEndfCont *)recIter.Next();
      LFW                  = cont1->GetL2();
      NER                  = cont1->GetN1();
      for (int j = 0; j < NER; j++) {
        TNudyEndfCont *cont2 = (TNudyEndfCont *)recIter.Next();
        fElo                 = cont2->GetC1();
        fEhi                 = cont2->GetC2();
        LRU                  = cont2->GetL1();
        LRF                  = cont2->GetL2();
        NRO                  = cont2->GetN1();
        NAPS                 = cont2->GetN2();
        if (LRU == 2 && LRF == 1 && LFW == 1) {
        } else {
          TNudyEndfCont *cont3 = (TNudyEndfCont *)recIter.Next();
          fSpi                 = cont3->GetC1();
          fAp                  = cont3->GetC2();
          if (LRU == 2) LSSF   = cont3->GetL1();
          // if(LRF==3)int LAD = cont3->GetL1();
          if (LRU == 1) NLS = cont3->GetN1();
          NLS2              = cont3->GetN1();
          fNRS.resize(NLS, 0);
          gjdeno = 4.0 * fSpi + 2.0;
        }
        switch (LRU) {
        case 0: {
        } break;
        case 1: { // resolved resonance region
	  std::cout <<"LRF "<< LRF << std::endl;
          switch (LRF) {
          case 1: // single level resonance region
          case 2: // multi level resonance region
          case 3: // RM resonance region
          case 4: // Adler - Adler resonance region
          {
            fFlagResolve      = 1;
            if (j == 0) fElo1 = fElo;
            fEhi1             = fEhi;
            for (int lseries = 0; lseries < NLS; lseries++) {
              TNudyEndfList *list                      = (TNudyEndfList *)recIter.Next();
              AWRI                                     = list->GetC1();
              fQX                                      = list->GetC2();
              if (fQX == 0.0) fApl[lseries]            = fAp;
              if (LRF == 3 && fQX > 0.0) fApl[lseries] = list->GetC2();
              l.push_back(list->GetL1());
              fLrx                = list->GetL2();
              fNRS[lseries]       = list->GetN2();
              fA                  = Mn * AWRI;
              Z                   = (int)(ZA - fA) / 1000.0;
              fRad_a              = 0.08 + 0.123 * pow(1.00866491578 * AWRI, (1. / 3.));
              if (fAp == 0.0) fAp = fRad_a;
              fFactor_k           = kconst * (AWRI / (AWRI + 1.0));
              fJMin               = (std::fabs(fSpi - l[lseries]) - 0.5);
              fJMax               = (fSpi + l[lseries] + 0.5);
              nrsl0               = (!lseries) ? 0 : nrsl0 + fNRS[lseries - 1]; // last fNRS value
              double jtemp[fNRS[lseries]];
              if (LRF == 1 || LRF == 2) {
                for (int ii = 0, eNRS = fNRS[lseries]; ii != eNRS; ++ii) {
                  fCueMat = nrsl0 + ii;
                  fEr.push_back(list->GetLIST(ii * 6 + 0));
                  J.push_back(list->GetLIST(ii * 6 + 1));                      // J value
                  fGamma_r.push_back(list->GetLIST(ii * 6 + 2));               // total width
                  fGamma_n.push_back(list->GetLIST(ii * 6 + 3));               // neutron width
                  fGamma_g.push_back(list->GetLIST(ii * 6 + 4));               // gamma width
                  fGamma_f.push_back(list->GetLIST(ii * 6 + 5));               // fission width
                  fGJ.push_back((2.0 * std::fabs(J[fCueMat]) + 1.0) / gjdeno); //(2J+1)/2(2I+1)
                  fPhiEr.push_back(CalcPene(std::fabs(list->GetLIST(ii * 6 + 0)), lseries));
                  fShiftEr.push_back(CalcShift(std::fabs(list->GetLIST(ii * 6 + 0)), lseries));
                  fGamma_n[nrsl0 + ii] = fGamma_n[nrsl0 + ii] / fPhiEr[nrsl0 + ii];
                  jtemp[ii]            = list->GetLIST(ii * 6 + 1);
                }
                double jmis = -99, temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9;
                for (int isort = 0, eNRS = fNRS[lseries]; isort != eNRS; ++isort) {
                  for (int isort1 = 1; isort1 != eNRS; ++isort1) {
                    if (jtemp[isort1] < jtemp[isort1 - 1]) {
                      temp              = jtemp[isort1];
                      jtemp[isort1]     = jtemp[isort1 - 1];
                      jtemp[isort1 - 1] = temp;

                      temp1                   = fEr[nrsl0 + isort1];
                      fEr[nrsl0 + isort1]     = fEr[nrsl0 + isort1 - 1];
                      fEr[nrsl0 + isort1 - 1] = temp1;

                      temp2                 = J[nrsl0 + isort1];
                      J[nrsl0 + isort1]     = J[nrsl0 + isort1 - 1];
                      J[nrsl0 + isort1 - 1] = temp2;

                      temp3                        = fGamma_r[nrsl0 + isort1];
                      fGamma_r[nrsl0 + isort1]     = fGamma_r[nrsl0 + isort1 - 1];
                      fGamma_r[nrsl0 + isort1 - 1] = temp3;

                      temp4                        = fGamma_n[nrsl0 + isort1];
                      fGamma_n[nrsl0 + isort1]     = fGamma_n[nrsl0 + isort1 - 1];
                      fGamma_n[nrsl0 + isort1 - 1] = temp4;

                      temp5                        = fGamma_g[nrsl0 + isort1];
                      fGamma_g[nrsl0 + isort1]     = fGamma_g[nrsl0 + isort1 - 1];
                      fGamma_g[nrsl0 + isort1 - 1] = temp5;

                      temp6                        = fGamma_f[nrsl0 + isort1];
                      fGamma_f[nrsl0 + isort1]     = fGamma_f[nrsl0 + isort1 - 1];
                      fGamma_f[nrsl0 + isort1 - 1] = temp6;

                      temp7                   = fGJ[nrsl0 + isort1];
                      fGJ[nrsl0 + isort1]     = fGJ[nrsl0 + isort1 - 1];
                      fGJ[nrsl0 + isort1 - 1] = temp7;

                      temp8                      = fPhiEr[nrsl0 + isort1];
                      fPhiEr[nrsl0 + isort1]     = fPhiEr[nrsl0 + isort1 - 1];
                      fPhiEr[nrsl0 + isort1 - 1] = temp8;

                      temp9                        = fShiftEr[nrsl0 + isort1];
                      fShiftEr[nrsl0 + isort1]     = fShiftEr[nrsl0 + isort1 - 1];
                      fShiftEr[nrsl0 + isort1 - 1] = temp9;
                    }
                  }
                }
                int ju               = 0;
                fNJValue[l[lseries]] = 0;
                fMisGj[l[lseries]]   = 0;
                for (int j1 = 0, eNRS = fNRS[lseries]; j1 != eNRS; ++j1) {
                  if (jtemp[j1] != jmis) {
                    fMissingJ[l[lseries]][ju] = jtemp[j1];
                    fNJValue[l[lseries]] += 1;
                    jmis = jtemp[j1];
                    ju += 1;
                  }
                }
              } else if (LRF == 3) {
                for (int ii = 0, eNRS = fNRS[lseries]; ii != eNRS; ++ii) {
                  fCueMat = nrsl0 + ii;
                  fEr.push_back(list->GetLIST(ii * 6 + 0));
                  J.push_back(list->GetLIST(ii * 6 + 1));                                  // J value
                  fGamma_n.push_back(list->GetLIST(ii * 6 + 2));                           // total width
                  fGamma_g.push_back(list->GetLIST(ii * 6 + 3));                           // neutron width
                  fGamma_fa.push_back(list->GetLIST(ii * 6 + 4));                          // gamma width
                  fGamma_fb.push_back(list->GetLIST(ii * 6 + 5));                          // fission width
                  fGamma_fasq.push_back(sqrt(0.5 * std::fabs(list->GetLIST(ii * 6 + 4)))); // gamma width
                  fGamma_fbsq.push_back(sqrt(0.5 * std::fabs(list->GetLIST(ii * 6 + 5)))); // fission width
                  if (fGamma_fa[fCueMat] < 0.0) fGamma_fasq[fCueMat] = -fGamma_fasq[fCueMat];
                  if (fGamma_fb[fCueMat] < 0.0) fGamma_fbsq[fCueMat] = -fGamma_fbsq[fCueMat];
                  fGJ.push_back((2.0 * std::fabs(J[fCueMat]) + 1.0) / gjdeno); //(2J+1)/2(2I+1)
                  jtemp[ii] = list->GetLIST(ii * 6 + 1);
                  fPhiEr.push_back(CalcPene(std::fabs(list->GetLIST(ii * 6 + 0)), lseries));
                  fShiftEr.push_back(CalcShift(std::fabs(list->GetLIST(ii * 6 + 0)), lseries));
                  fGamma_n[nrsl0 + ii] = fGamma_n[nrsl0 + ii] / fPhiEr[nrsl0 + ii];
                }
                double gjfound = 0.0, jmis = -99, temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9,
                       temp10, temp11;
                for (int isort = 0, eNRS = fNRS[lseries]; isort != eNRS; ++isort) {
                  for (int isort1 = 1; isort1 != eNRS; ++isort1) {
                    if (jtemp[isort1] < jtemp[isort1 - 1]) {
                      temp              = jtemp[isort1];
                      jtemp[isort1]     = jtemp[isort1 - 1];
                      jtemp[isort1 - 1] = temp;

                      temp1                   = fEr[nrsl0 + isort1];
                      fEr[nrsl0 + isort1]     = fEr[nrsl0 + isort1 - 1];
                      fEr[nrsl0 + isort1 - 1] = temp1;

                      temp2                 = J[nrsl0 + isort1];
                      J[nrsl0 + isort1]     = J[nrsl0 + isort1 - 1];
                      J[nrsl0 + isort1 - 1] = temp2;

                      temp3                         = fGamma_fa[nrsl0 + isort1];
                      fGamma_fa[nrsl0 + isort1]     = fGamma_fa[nrsl0 + isort1 - 1];
                      fGamma_fa[nrsl0 + isort1 - 1] = temp3;

                      temp4                        = fGamma_n[nrsl0 + isort1];
                      fGamma_n[nrsl0 + isort1]     = fGamma_n[nrsl0 + isort1 - 1];
                      fGamma_n[nrsl0 + isort1 - 1] = temp4;

                      temp5                        = fGamma_g[nrsl0 + isort1];
                      fGamma_g[nrsl0 + isort1]     = fGamma_g[nrsl0 + isort1 - 1];
                      fGamma_g[nrsl0 + isort1 - 1] = temp5;

                      temp6                         = fGamma_fb[nrsl0 + isort1];
                      fGamma_fb[nrsl0 + isort1]     = fGamma_fb[nrsl0 + isort1 - 1];
                      fGamma_fb[nrsl0 + isort1 - 1] = temp6;

                      temp7                   = fGJ[nrsl0 + isort1];
                      fGJ[nrsl0 + isort1]     = fGJ[nrsl0 + isort1 - 1];
                      fGJ[nrsl0 + isort1 - 1] = temp7;

                      temp8                      = fPhiEr[nrsl0 + isort1];
                      fPhiEr[nrsl0 + isort1]     = fPhiEr[nrsl0 + isort1 - 1];
                      fPhiEr[nrsl0 + isort1 - 1] = temp8;

                      temp9                        = fShiftEr[nrsl0 + isort1];
                      fShiftEr[nrsl0 + isort1]     = fShiftEr[nrsl0 + isort1 - 1];
                      fShiftEr[nrsl0 + isort1 - 1] = temp9;

                      temp10                          = fGamma_fasq[nrsl0 + isort1];
                      fGamma_fasq[nrsl0 + isort1]     = fGamma_fasq[nrsl0 + isort1 - 1];
                      fGamma_fasq[nrsl0 + isort1 - 1] = temp10;

                      temp11                          = fGamma_fbsq[nrsl0 + isort1];
                      fGamma_fbsq[nrsl0 + isort1]     = fGamma_fbsq[nrsl0 + isort1 - 1];
                      fGamma_fbsq[nrsl0 + isort1 - 1] = temp11;
                    }
                  }
                }
                int ju               = 0;
                fNJValue[l[lseries]] = 0;
                fMisGj[l[lseries]]   = 0;
                for (int j1 = 0, eNRS = fNRS[lseries]; j1 != eNRS; ++j1) {
                  if (jtemp[j1] != jmis) {
                    fMissingJ[l[lseries]][ju] = jtemp[j1];
                    fNJValue[l[lseries]] += 1;
                    jmis = jtemp[j1];
                    gjfound += (2 * std::fabs(jtemp[j1]) + 1) / gjdeno;
                    ju += 1;
                  }
                }
                fMisGj[l[lseries]] = 2 * l[lseries] + 1 - gjfound;
              } else if (LRF == 4) {                         // Adler -Adler resonance region
                for (int ii = 0, eNRS = fNRS[lseries]; ii != eNRS; ++ii) { // line 6 onwards data
                  fCueMat = nrsl0 + ii;
                  fAt1.push_back(list->GetLIST(ii * 6 + 0));
                  fAt2.push_back(list->GetLIST(ii * 6 + 1));
                  fAt3.push_back(list->GetLIST(ii * 6 + 2));
                  fAt4.push_back(list->GetLIST(ii * 6 + 3));
                  fBt1.push_back(list->GetLIST(ii * 6 + 4));
                  fBt2.push_back(list->GetLIST(ii * 6 + 5));
                }
                TNudyEndfList *list1 = (TNudyEndfList *)recIter.Next();
                for (int ii = 0, eN2 = list1->GetN2(); ii != eN2; ++ii) {
                  fDet1.push_back(list1->GetLIST(ii * 6 + 0));
                  fDwt1.push_back(list1->GetLIST(ii * 6 + 1));
                  fGrt1.push_back(list1->GetLIST(ii * 6 + 2));
                  fGit1.push_back(list1->GetLIST(ii * 6 + 3));
                  fDef1.push_back(list1->GetLIST(ii * 6 + 4));
                  fDwf1.push_back(list1->GetLIST(ii * 6 + 5));
                  fGrf1.push_back(list1->GetLIST(ii * 6 + 6 + 0));
                  fGif1.push_back(list1->GetLIST(ii * 6 + 6 + 1));
                  fDec1.push_back(list1->GetLIST(ii * 6 + 6 + 2));
                  fDwc1.push_back(list1->GetLIST(ii * 6 + 6 + 3));
                  fGrc1.push_back(list1->GetLIST(ii * 6 + 6 + 4));
                  fGic1.push_back(list1->GetLIST(ii * 6 + 6 + 5));
                }
              }
            } // loop for L values
          } break;
          case 7: { // R-Matrix
          } break;
          } // break;
          Linearize(j);
          // clearing vectors for different energy regions
          fNRS.clear();
          l.clear();
          fEr.clear();
          J.clear();
          fGamma_r.clear();
          fGamma_n.clear();
          fGamma_g.clear();
          fGamma_f.clear();
          fGJ.clear();
          fPhiEr.clear();
          fShiftEr.clear();
          fGamma_fa.clear();
          fGamma_fasq.clear();
          fGamma_fb.clear();
          fGamma_fbsq.clear();
          fAt1.clear();
          fAt2.clear();
          fAt3.clear();
          fAt4.clear();
          fBt1.clear();
          fBt2.clear();
          fDet1.clear();
          fDwt1.clear();
          fGrt1.clear();
          fGit1.clear();
          fDef1.clear();
          fDwf1.clear();
          fGrf1.clear();
          fGif1.clear();
          fDec1.clear();
          fDwc1.clear();
          fGrc1.clear();
          fGic1.clear();

        } break;
        case 2: { // unresolved resonance region
          fFlagUnResolve = 1;
          switch (LRF) {
          case 2: {
            for (int lseries = 0; lseries < NLS2; lseries++) {
              TNudyEndfList *list2 = (TNudyEndfList *)recIter.Next();
              AWRI                 = list2->GetC1();
              l.push_back(list2->GetL1());
              NJS = list2->GetN1();
              fNRJ.push_back(NJS);
              fApl[lseries] = fAp;
              for (int jseries = 0; jseries < NJS; jseries++) {
                TNudyEndfList *list3 = (TNudyEndfList *)recIter.Next();
                fJSM.push_back(list3->GetN2());
                fRad_a              = 0.08 + 0.123 * pow(1.00866491578 * AWRI, (1. / 3.));
                if (fAp == 0.0) fAp = fRad_a;
                fFactor_k           = kconst * (AWRI / (AWRI + 1.0));
                fJMin               = list3->GetC1();
                INT                 = list3->GetL1();
                fAmux.push_back(list3->GetLIST(2));           // total width
                fAmun.push_back(list3->GetLIST(3));           // neutron width
                fAmug.push_back(list3->GetLIST(4));           // gamma width
                fAmuf.push_back(list3->GetLIST(5));           // fission width
                fGJ.push_back((2.0 * fJMin + 1.0) / gjdeno);  //(2J+1)/2(2I+1)
                for (int ii = 0, eN2 = list3->GetN2(); ii != eN2; ++ii) { // line 6 onwards data
                  fEs.push_back(list3->GetLIST((ii + 1) * 6 + 0));
                  fD.push_back(list3->GetLIST((ii + 1) * 6 + 1));   // J value
                  fGX.push_back(list3->GetLIST((ii + 1) * 6 + 2));  // total width
                  fGNO.push_back(list3->GetLIST((ii + 1) * 6 + 3)); // neutron width
                  fGG.push_back(list3->GetLIST((ii + 1) * 6 + 4));  // gamma width
                  fGF.push_back(list3->GetLIST((ii + 1) * 6 + 5));  // fission width
                }
              } // complete NJS reading
            }   // complete NLS2 reading
            fElo2 = fElo;
            fEhi2 = fEhi;
            // std::cout << " fElo2 "<< fElo2 <<" fEhi2 "<< fEhi2 << std::endl;
            if (LSSF == 0) Linearize(j);
            l.clear();
            fNRS.clear();
            fNRJ.clear();
            fJSM.clear();
            fAmux.clear();
            fAmun.clear();
            fAmug.clear();
            fAmuf.clear();
            fGJ.clear();
            fEs.clear();
            fD.clear();
            fGX.clear();
            fGNO.clear();
            fGG.clear();
            fGF.clear();
          } break;
          case 1: {
            switch (LFW) {
            case 0: {
              for (int lseries = 0; lseries < NLS2; lseries++) {
                TNudyEndfList *list2 = (TNudyEndfList *)recIter.Next();
                AWRI                 = list2->GetC1();
                l.push_back(list2->GetL1());
                NJS = list2->GetN2();
                fNRJ.push_back(NJS);
                fApl[lseries] = fAp;
                for (int ii = 0; ii < NJS; ii++) { // line 6 onwards data
                  fJSM.push_back(1);
                  fRad_a              = 0.08 + 0.123 * pow(1.00866491578 * AWRI, (1. / 3.));
                  if (fAp == 0.0) fAp = fRad_a;
                  fFactor_k           = kconst * (AWRI / (AWRI + 1.0));
                  INT                 = 2;
                  fEs.push_back(fElo);
                  fD.push_back(list2->GetLIST((ii)*6 + 0));    // Average level spacing for resonances with spin J
                  fAmun.push_back(list2->GetLIST((ii)*6 + 2)); // degree of freedom neutron width
                  fGNO.push_back(list2->GetLIST((ii)*6 + 3));  // neutron width
                  fGG.push_back(list2->GetLIST((ii)*6 + 4));   // gamma width
                  fGF.push_back(0.0);                          // fission width
                  fGX.push_back(0.0);                          // fission width
                  fAmug.push_back(0.0);                        // degree of freedom gamma width
                  fAmuf.push_back(0.0);                        // degree of freedom fission width
                  fAmux.push_back(0.0);                        // degree of freedom competitive width
                  fGJ.push_back((2.0 * list2->GetLIST((ii)*6 + 1) + 1.0) / gjdeno); //(2J+1)/2(2I+1)
                }
              }
              fElo2 = fElo;
              fEhi2 = fEhi;
              if (LSSF == 0) Linearize(j);
              l.clear();
              fNRS.clear();
              fNRJ.clear();
              fJSM.clear();
              fAmux.clear();
              fAmun.clear();
              fAmug.clear();
              fAmuf.clear();
              fGJ.clear();
              fEs.clear();
              fD.clear();
              fGX.clear();
              fGNO.clear();
              fGG.clear();
              fGF.clear();
            } break;
            case 1: {
              TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
              fSpi                = list->GetC1();
              fAp                 = list->GetC2();
              LSSF                = list->GetL1();
              NLS2                = list->GetN2();
              fNE                 = list->GetN1();
              fNRS.resize(NLS, 0);
              gjdeno = 4.0 * fSpi + 2.0;
              for (int lseries = 0; lseries < NLS2; lseries++) {
                TNudyEndfCont *list2 = (TNudyEndfCont *)recIter.Next();
                AWRI                 = list2->GetC1();
                l.push_back(list2->GetL1());
                NJS = list2->GetN1();
                fNRJ.push_back(NJS);
                fApl[lseries] = fAp;
                for (int jseries = 0; jseries < NJS; jseries++) {
                  TNudyEndfList *list3 = (TNudyEndfList *)recIter.Next();
                  fJSM.push_back(fNE);
                  fRad_a              = 0.08 + 0.123 * pow(1.00866491578 * AWRI, (1. / 3.));
                  if (fAp == 0.0) fAp = fRad_a;
                  fFactor_k           = kconst * (AWRI / (AWRI + 1.0));
                  INT                 = 2;
                  fAmux.push_back(0.0);                                      // total width
                  fAmun.push_back(list3->GetLIST(2));                        // neutron width
                  fAmug.push_back(0.0);                                      // gamma width
                  fAmuf.push_back(list3->GetL2());                           // fission width
                  fGJ.push_back((2.0 * list3->GetLIST((1)) + 1.0) / gjdeno); //(2J+1)/2(2I+1)
                  for (int ii = 0; ii < fNE; ii++) {                         // line 6 onwards data
                    fEs.push_back(list->GetLIST(ii));
                    fGF.push_back(list3->GetLIST(6 + ii)); // fission width
                    fGG.push_back(list3->GetLIST(4));      // gamma width
                    fGX.push_back(0.0);                    // gamma width
                    fD.push_back(list3->GetLIST(0));       // J value
                    fGNO.push_back(list3->GetLIST(3));     // neutron width
                  }
                } // complete NJS reading
              }
              fElo2 = fElo;
              fEhi2 = fEhi;
              if (LSSF == 0) Linearize(j);
              l.clear();
              fNRS.clear();
              fNRJ.clear();
              fJSM.clear();
              fAmux.clear();
              fAmun.clear();
              fAmug.clear();
              fAmuf.clear();
              fGJ.clear();
              fEs.clear();
              fD.clear();
              fGX.clear();
              fGNO.clear();
              fGG.clear();
              fGF.clear();
            } break;
            }
          } break;
          }
        } break;
        }
      }
    }
  }
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::RecursionLinearFile3(double x1, double x2, double sig1, double sig2, std::vector<double> ene,
                                            std::vector<double> sig)
{
  double siga;
  double mid = 0.5 * (x1 + x2);
  if (sig1 == 0.0 || fMloop > 1000) return 0;
  if (sig1 == 0.0 && sig2 == 0.0) return 0;
  if (x1 == x2 || x1 < 1E-5 || x2 < 1E-5) {
    return 0;
  }
  siga           = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, ene, sig, fNP, mid);
  double sigmid1 = sig1 + (sig2 - sig1) * (mid - x1) / (x2 - x1);
  if (siga == 0.0 || x1 == x2) return 0;
  // std::cout <<"x \t"<< x1 <<"  \t"<< x2 <<"  \t"<< sig1 <<"  \t"<< sig2 <<"  \t"<< sigmid1 <<"   \t"<< siga <<
  // std::endl;
  if (std::fabs((siga - sigmid1) / siga) <= fSigDiff) {
    return 0;
  } else {
    // std::cout <<"mid \t"<< x1 <<"  \t"<< x2 <<"  \t"<< sig1 <<"  \t"<< sig2 <<"  \t"<< sigmid1 <<"   \t"<< siga <<"
    // \t"<< std::fabs((siga - sigmid1) / siga )<< std::endl;
    fELinearFile3.push_back(mid);
    fXLinearFile3.push_back(siga);
  }
  RecursionLinearFile3(x1, mid, sig1, siga, ene, sig);
  RecursionLinearFile3(mid, x2, siga, sig2, ene, sig);
  return 0;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::AddFile3Resonance(double &x1, double &x2, std::vector<double> &x3, std::vector<double> &x4)
{
  int linsize = x3.size();
  int size1   = fELinearFile3.size();
  for (int cr1 = 0; cr1 < linsize; cr1++) {
    if (x3[cr1] >= x1 && x3[cr1] <= x2) {
      double crs = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fELinearFile3, fXLinearFile3, size1, x3[cr1]);
      x4[cr1] += crs;
    }
  }
  return 0;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::InsertFile3(std::vector<double> &x1, std::vector<double> &x2)
{
  int size1 = fELinearFile3.size();
  for (int cr1 = 0; cr1 < size1; cr1++) {
    if ((fELinearFile3[cr1] > fElo2 && fELinearFile3[cr1] <= fEhi2)) {
      x1.push_back(fELinearFile3[cr1]);
      x2.push_back(fXLinearFile3[cr1]);
    }
  }
  return 0;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::InsertFile3High(std::vector<double> &x1, std::vector<double> &x2)
{
  int eELinearFile3Size = fELinearFile3.size();
  if (fFlagResolve != 0) {
    for (int cr = 0; cr != eELinearFile3Size; ++cr) {
      if (fXLinearFile3[cr] > 0.0 && fELinearFile3[cr] <= fElo1) {
        x1.push_back(fELinearFile3[cr]);
        x2.push_back(fXLinearFile3[cr]);
      }
    }
  } else if (fFlagResolve == 0 && fFlagUnResolve != 0) {
    for (int cr = 0; cr != eELinearFile3Size; ++cr) {
      if (fXLinearFile3[cr] > 0.0 && fELinearFile3[cr] <= fElo2) {
        x1.push_back(fELinearFile3[cr]);
        x2.push_back(fXLinearFile3[cr]);
      }
    }
  }
  if (fFlagUnResolve != 0) {
    for (int cr = 0; cr != eELinearFile3Size; ++cr) {
      if (fXLinearFile3[cr] > 0.0 && fELinearFile3[cr] > fEhi2) {
        x1.push_back(fELinearFile3[cr]);
        x2.push_back(fXLinearFile3[cr]);
      }
    }
  } else if (fFlagUnResolve == 0) {
    for (int cr = 0; cr != eELinearFile3Size; ++cr) {
      if (fXLinearFile3[cr] > 0.0 && fELinearFile3[cr] > fEhi1) {
        x1.push_back(fELinearFile3[cr]);
        x2.push_back(fXLinearFile3[cr]);
      }
    }
  }
  return 0;
}

//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::ReadFile3(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  std::vector<double> eneExtra;
  if (LRF == 0) fEhi = 0.0;
  double eneLow      = 1E-5;
  do {
    eneExtra.push_back(eneLow);
    eneLow *= 1.15;
  } while (eneLow < 1E3);

  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    int fMT = sec->GetMT();
    /*
    //     if (fMT != 1 && fMT != 3 && fMT != 4 && fMT != 5 && fMT != 27 && fMT != 19 && fMT != 20 && fMT != 21 && fMT
    != 38 && fMT != 101 &&
    //         fMT < 250) {
    //       MtNumbers.push_back(fMT);
    */
    TIter recIter(sec->GetRecords());
    TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
    ZA                    = sec->GetC1();
    //      fAWR	  	= sec->GetC2();
    fMAT = sec->GetMAT();
    /*
          // double QM   = header->GetC1();
    //       double QI            = header->GetC2();
    //       fQValue[sec->GetMT()] = QI;
    //       fQvalueTemp.push_back(QI);
    //       fQvalueTemp.push_back(fMT);
          // int LR = header->GetL2();
    */
    fNR = header->GetN1();
    fNP = header->GetN2();
    //       std::cout<<sec->GetMT() <<" fMT "<< fNR <<"  "<< fNP << std::endl;
    TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
    for (int cr = 0; cr < fNR; cr++) {
      fNbt1.push_back(tab1->GetNBT(cr));
      fInt1.push_back(tab1->GetINT(cr));
    }
    int npp = 0;
    for (int crs = 0; crs < fNP; crs++) {
      fEneTemp.push_back(tab1->GetX(crs));
      fSigTemp.push_back(tab1->GetY(crs));
      if (fPrepro == 0) {
        if (npp > 0 && fEneTemp[fEneTemp.size() - 1] == fEneTemp[fEneTemp.size() - 2]) {
          fEneTemp[fEneTemp.size() - 1] += 1E-9; // adding small difference for boundary points
        }
        if (npp > 0 && fEneTemp[fEneTemp.size() - 1] == fEneTemp[fEneTemp.size() - 3]) {
          fEneTemp[fEneTemp.size() - 1] += 2E-9; // adding small difference for boundary points
        }
        if (npp > 0 && fEneTemp[fEneTemp.size() - 1] == fEneTemp[fEneTemp.size() - 4]) {
          fEneTemp[fEneTemp.size() - 1] += 3E-9; // adding small difference for boundary points
        }
      }
      npp += 1;
    }
    if (fPrepro == 0) {
      TNudyCore::Instance()->Sort(fEneTemp, fSigTemp);
      TNudyCore::Instance()->ThinningDuplicate(fEneTemp, fSigTemp);
    }
    fNP = fEneTemp.size();
    for (int crs = 0; crs < fNP; crs++) {
      fELinearFile3.push_back(fEneTemp[crs]);
      fXLinearFile3.push_back(fSigTemp[crs]);
    }
    // Linearization of file 3 data
    if (fPrepro == 0) {
      for (int cr = 0; cr < fNP - 1; cr++) {
        fMloop = 0;
        RecursionLinearFile3(fELinearFile3[cr], fELinearFile3[cr + 1], fXLinearFile3[cr], fXLinearFile3[cr + 1],
                             fELinearFile3, fXLinearFile3);
      }
      double enehigh = 1E+3;
      eneLow         = 1E+3;
      int multifac   = 2;
      do {
        eneLow += multifac * enehigh;
        eneExtra.push_back(eneLow);
        multifac += 1;
      } while (eneLow < fEneTemp[fNP - 1]);
    
      for (int ie = 0, eneExtraSize = eneExtra.size(); ie != eneExtraSize; ++ie) {
        double sigExtra =
            TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fELinearFile3, fXLinearFile3, fNP, eneExtra[ie]);
        if (sigExtra > 0.0) {
          fELinearFile3.push_back(eneExtra[ie]);
          fXLinearFile3.push_back(sigExtra);
        }
      }
      TNudyCore::Instance()->Sort(fELinearFile3, fXLinearFile3);
      TNudyCore::Instance()->ThinningDuplicate(fELinearFile3, fXLinearFile3);
    }
    //	Filling of array to interpolate and add with file 2 data if it is given
    fNbt1.clear();
    fInt1.clear();

    fNbt1.push_back(fELinearFile3.size());
    fInt1.push_back(2);
    fNR = 1;
    /////////////////////////////////////////////////////////////////////
    //	resolved range is always added in file 3
    if (fPrepro == 0) {
      if (LRU != 0) {
        if (fFlagResolve != 0) {
          switch (fMT) {
          case 2: {
            AddFile3Resonance(fElo1, fEhi1, fELinElastic, fXLinElastic);
          } break;
          case 18: {
            AddFile3Resonance(fElo1, fEhi1, fELinFission, fXLinFission);
          } break;
          case 102: {
            AddFile3Resonance(fElo1, fEhi1, fELinCapture, fXLinCapture);
          } break;
          }
        }
        // unresolved resonance region is added if LSSF = 0
        if (fFlagUnResolve != 0) {
          switch (LSSF) {
          case 0: {
            switch (fMT) {
            case 2: {
              AddFile3Resonance(fElo2, fEhi2, fELinElastic, fXLinElastic);
            } break;
            case 18: {
              AddFile3Resonance(fElo2, fEhi2, fELinFission, fXLinFission);
            } break;
            case 102: {
              AddFile3Resonance(fElo2, fEhi2, fELinCapture, fXLinCapture);
            } break;
            }
          } break;
          case 1: {
            switch (fMT) {
            case 2: {
              InsertFile3(fELinElastic, fXLinElastic);
            } break;
            case 18: {
              InsertFile3(fELinFission, fXLinFission);
            } break;
            case 102: {
              InsertFile3(fELinCapture, fXLinCapture);
            } break;
            }
          } break;
          }
        }
        // This adds additional points to the elastic, fission and capture below resolved range and above uresolved
        // range
        switch (fMT) {
        case 2: {
          InsertFile3High(fELinElastic, fXLinElastic);
        } break;
        case 18: {
          InsertFile3High(fELinFission, fXLinFission);
        } break;
        case 102: {
          InsertFile3High(fELinCapture, fXLinCapture);
        } break;
        }
      } else { // no resonance parameters are given
           int ELinearFile3Size = fELinearFile3.size();
        switch (fMT) {
        case 2: {
          for (int cr = 0; cr != ELinearFile3Size; ++cr) {
            fELinElastic.push_back(fELinearFile3[cr]);
            fXLinElastic.push_back(fXLinearFile3[cr]);
          }
        } break;
        case 18: {
          for (int cr = 0; cr != ELinearFile3Size; ++cr) {
            fELinFission.push_back(fELinearFile3[cr]);
            fXLinFission.push_back(fXLinearFile3[cr]);
          }
        } break;
        case 102: {
          for (int cr = 0; cr != ELinearFile3Size; ++cr) {
            fELinCapture.push_back(fELinearFile3[cr]);
            fXLinCapture.push_back(fXLinearFile3[cr]);
          }
        } break;
        }
      }
      int fNR = 1;
      if (fMT == 2) {
        int ELinElasticSize = fELinElastic.size(); 
        header->SetCont(header->GetC1(), header->GetC2(), header->GetL1(), header->GetL2(), fNR,
                        ELinElasticSize);
        tab1->SetNBT(ELinElasticSize, 0);
        tab1->SetINT(2, 0);
        for (int crs = 0; crs != ELinElasticSize; ++crs) {
          tab1->SetX(fELinElastic[crs], crs);
          tab1->SetY(fXLinElastic[crs], crs);
        }
      } else if (fMT == 18) {
        int ELinFissionSize = fELinFission.size();
        header->SetCont(header->GetC1(), header->GetC2(), header->GetL1(), header->GetL2(), fNR,
                        ELinFissionSize);
        tab1->SetNBT(ELinFissionSize, 0);
        tab1->SetINT(2, 0);
        for (int crs = 0; crs != ELinFissionSize; ++crs) {
          tab1->SetX(fELinFission[crs], crs);
          tab1->SetY(fXLinFission[crs], crs);
        }
      } else if (fMT == 102) {
        int ELinCaptureSize = fELinCapture.size();
        header->SetCont(header->GetC1(), header->GetC2(), header->GetL1(), header->GetL2(), fNR,
                        ELinCaptureSize);
        tab1->SetNBT(ELinCaptureSize, 0);
        tab1->SetINT(2, 0);
        for (int crs = 0; crs != ELinCaptureSize; ++crs) {
          tab1->SetX(fELinCapture[crs], crs);
          tab1->SetY(fXLinCapture[crs], crs);
        }
      } else if (fMT != 2 && fMT != 18 && fMT != 102) {
	int ELinearFile3Size = fELinearFile3.size();
        header->SetCont(header->GetC1(), header->GetC2(), header->GetL1(), header->GetL2(), fNR,
                        ELinearFile3Size);
        tab1->SetNBT(ELinearFile3Size, 0);
        tab1->SetINT(2, 0);
        for (int crs = 0; crs != ELinearFile3Size; ++crs) {
          tab1->SetX(fELinearFile3[crs], crs);
          tab1->SetY(fXLinearFile3[crs], crs);
          //	    std::cout <<" fMT "<< sec->GetMT()<<" energy "<< fELinearFile3 [crs] <<" fSigma "<< fXLinearFile3
          //[crs] << std::endl;
        }
      }

      // 	std::cout <<" proton flag \t" <<  fMTChargeFlag [0] <<" proton flag \t" <<  fMTChargeFlag [5] <<" fMT \t" <<
      // fMT << std::endl;
      // if data are given only in fMT = 600-850
      if (fMTChargeFlag[0] != -1 && (fMT >= 600 && fMT < 650)) {
        //  summing for n,p reaction
        // 	  std::cout<<" first case0 != -1 fMT= 600-649 "<< std::endl;
        fEneUniP.insert(std::end(fEneUniP), std::begin(fELinearFile3), std::end(fELinearFile3));
        fEneUniPAll.insert(std::end(fEneUniPAll), std::begin(fELinearFile3), std::end(fELinearFile3));
      } else if (fMTChargeFlag[1] != -1 && fMT >= 650 && fMT < 700) {
        //  summing for n,d reaction
        // 	  std::cout<<" first case0 != -1 fMT= 650-700 "<< std::endl;
        fEneUniD.insert(std::end(fEneUniD), std::begin(fELinearFile3), std::end(fELinearFile3));
        fEneUniDAll.insert(std::end(fEneUniDAll), std::begin(fELinearFile3), std::end(fELinearFile3));
      } else if (fMTChargeFlag[2] != -1 && fMT >= 700 && fMT < 750) {
        //  summing for n,t reaction
        // 	  std::cout<<" first case0 != -1 fMT= 700-749 "<< std::endl;
        fEneUniT.insert(std::end(fEneUniT), std::begin(fELinearFile3), std::end(fELinearFile3));
        fEneUniTAll.insert(std::end(fEneUniTAll), std::begin(fELinearFile3), std::end(fELinearFile3));
      } else if (fMTChargeFlag[3] != -1 && fMT >= 750 && fMT < 800) {
        //  summing for n,He3 reaction
        // 	  std::cout<<" first case0 != -1 fMT= 749-800 "<< std::endl;
        fEneUniHe3.insert(std::end(fEneUniHe3), std::begin(fELinearFile3), std::end(fELinearFile3));
        fEneUniHe3All.insert(std::end(fEneUniHe3All), std::begin(fELinearFile3), std::end(fELinearFile3));
      } else if (fMTChargeFlag[4] != -1 && fMT >= 800 && fMT < 850) {
        //  summing for n,He4 reaction
        // 	  std::cout<<" first case0 != -1 fMT= 800-849 "<< std::endl;
        fEneUniHe4.insert(std::end(fEneUniHe4), std::begin(fELinearFile3), std::end(fELinearFile3));
        fEneUniHe4All.insert(std::end(fEneUniHe4All), std::begin(fELinearFile3), std::end(fELinearFile3));
      } else if (fMT == 28 || (fMT >= 41 && fMT < 46) || fMT == 111 || fMT == 112 || fMT == 115 || fMT == 116 ||
                 fMT == 103) {
        //  summing for n,px reaction
        fEneUniPAll.insert(std::end(fEneUniPAll), std::begin(fELinearFile3), std::end(fELinearFile3));
      } else if (fMT == 11 || fMT == 32 || fMT == 35 || fMT == 114 || fMT == 115 || fMT == 117 || fMT == 104) {
        //  summing for n,dx reaction
        fEneUniDAll.insert(std::end(fEneUniDAll), std::begin(fELinearFile3), std::end(fELinearFile3));
      } else if (fMT == 33 || fMT == 113 || fMT == 116 || fMT == 105) {
        //  summing for n,tx reaction
        fEneUniTAll.insert(std::end(fEneUniTAll), std::begin(fELinearFile3), std::end(fELinearFile3));
      } else if (fMT == 34 || fMT == 113 || fMT == 116 || fMT == 106) {
        //  summing for n,He3x reaction
        fEneUniHe3All.insert(std::end(fEneUniHe3All), std::begin(fELinearFile3), std::end(fELinearFile3));
      } else if ((fMT >= 22 && fMT < 26) || fMT == 29 || fMT == 30 || fMT == 36 || fMT == 36 || fMT == 45 ||
                 fMT == 108 || fMT == 109 || (fMT >= 112 && fMT < 115) || fMT == 117 || fMT == 107) {
        //  summing for n,He4x reaction
        fEneUniHe4All.insert(std::end(fEneUniHe4All), std::begin(fELinearFile3), std::end(fELinearFile3));
      }
      //	check the sum of cross sections with total cross-section given in the fMT=1 + resonance
      //	    int size = fELinearFile3.size();
      if (fMT == 2) {
        fSigmaMts.insert(std::end(fSigmaMts), std::begin(fELinElastic), std::end(fELinElastic));
        fSigmaMts.insert(std::end(fSigmaMts), std::begin(fXLinElastic), std::end(fXLinElastic));
        // 	  energyUni.insert(std::end(energyUni), std::begin(fELinElastic), std::end(fELinElastic));
      } else if (fMT == 18) {
        fSigmaMts.insert(std::end(fSigmaMts), std::begin(fELinFission), std::end(fELinFission));
        fSigmaMts.insert(std::end(fSigmaMts), std::begin(fXLinFission), std::end(fXLinFission));
        // 	  energyUni.insert(std::end(energyUni), std::begin(fELinFission), std::end(fELinFission));
      } else if (fMT == 102) {
        fSigmaMts.insert(std::end(fSigmaMts), std::begin(fELinCapture), std::end(fELinCapture));
        fSigmaMts.insert(std::end(fSigmaMts), std::begin(fXLinCapture), std::end(fXLinCapture));
        // 	  energyUni.insert(std::end(energyUni), std::begin(fELinCapture), std::end(fELinCapture));
      } else if (fMT != 2 && fMT != 18 && fMT != 102) {
        fSigmaMts.insert(std::end(fSigmaMts), std::begin(fELinearFile3), std::end(fELinearFile3));
        fSigmaMts.insert(std::end(fSigmaMts), std::begin(fXLinearFile3), std::end(fXLinearFile3));
        // 	  energyUni.insert(std::end(energyUni), std::begin(fELinearFile3), std::end(fELinearFile3));
      }
    }
    if (fPrepro == 1) {
      if (fMT == 2) {
        fELinElastic.insert(std::end(fELinElastic), std::begin(fELinearFile3), std::end(fELinearFile3));
        fXLinElastic.insert(std::end(fXLinElastic), std::begin(fXLinearFile3), std::end(fXLinearFile3));
      } else if (fMT == 18) {
        fELinFission.insert(std::end(fELinFission), std::begin(fELinearFile3), std::end(fELinearFile3));
        fXLinFission.insert(std::end(fXLinFission), std::begin(fXLinearFile3), std::end(fXLinearFile3));
      } else if (fMT == 102) {
        fELinCapture.insert(std::end(fELinCapture), std::begin(fELinearFile3), std::end(fELinearFile3));
        fXLinCapture.insert(std::end(fXLinCapture), std::begin(fXLinearFile3), std::end(fXLinearFile3));
      }
      fSigmaMts.insert(std::end(fSigmaMts), std::begin(fELinearFile3), std::end(fELinearFile3));
      fSigmaMts.insert(std::end(fSigmaMts), std::begin(fXLinearFile3), std::end(fXLinearFile3));
      //         energyUni.insert(std::end(energyUni), std::begin(fELinearFile3), std::end(fELinearFile3));
    }
    fELinearFile3.clear();
    fXLinearFile3.clear();
    fEneTemp.clear();
    fSigTemp.clear();
    fNbt1.clear();
    fInt1.clear();
    fSigmaOfMts.push_back(fSigmaMts);
    fSigmaMts.clear();
    //     }
  }
  // sorting and removing duplicate points
  std::sort(fEneUniP.begin(), fEneUniP.end());
  TNudyCore::Instance()->ThinningDuplicate(fEneUniP);

  std::sort(fEneUniD.begin(), fEneUniD.end());
  TNudyCore::Instance()->ThinningDuplicate(fEneUniD);

  std::sort(fEneUniT.begin(), fEneUniT.end());
  TNudyCore::Instance()->ThinningDuplicate(fEneUniT);

  std::sort(fEneUniHe3.begin(), fEneUniHe3.end());
  TNudyCore::Instance()->ThinningDuplicate(fEneUniHe3);

  std::sort(fEneUniHe4.begin(), fEneUniHe4.end());
  TNudyCore::Instance()->ThinningDuplicate(fEneUniHe4);

  std::sort(fEneUniPAll.begin(), fEneUniPAll.end());
  TNudyCore::Instance()->ThinningDuplicate(fEneUniPAll);

  std::sort(fEneUniDAll.begin(), fEneUniDAll.end());
  TNudyCore::Instance()->ThinningDuplicate(fEneUniDAll);

  std::sort(fEneUniTAll.begin(), fEneUniTAll.end());
  TNudyCore::Instance()->ThinningDuplicate(fEneUniTAll);

  std::sort(fEneUniHe3All.begin(), fEneUniHe3All.end());
  TNudyCore::Instance()->ThinningDuplicate(fEneUniHe3All);

  std::sort(fEneUniHe4All.begin(), fEneUniHe4All.end());
  TNudyCore::Instance()->ThinningDuplicate(fEneUniHe4All);

  // summing cross-section for making total from ground and excited states for n,p
  TIter secIter1(file->GetSections());
  TNudyEndfSec *sec1;
  while ((sec1 = (TNudyEndfSec *)secIter1.Next())) {
    int fMT = sec1->GetMT();
    //     std::cout<<" fMT "<< fMT << std::endl;
    TIter recIter1(sec1->GetRecords());
    TNudyEndfCont *header = (TNudyEndfCont *)recIter1.Next();
    fNR                   = header->GetN1();
    fNP                   = header->GetN2();
    TNudyEndfTab1 *tab1   = (TNudyEndfTab1 *)(sec1->GetRecords()->At(0));
    if ((fMTChargeFlag[0] != -1 && fMT >= 600 && fMT < 650) /* || (fMTChargeFlag [5] != -1 && fMT == 103)
        || (fMTChargeFlag [0] == -1 && fMTChargeFlag [5] == -1 && fMT == 103)*/) {
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab1->GetX(crs));
        fXLinearFile3.push_back(tab1->GetY(crs));
      }
      for (int k = 0, EneUniPSize = fEneUniP.size(); k != EneUniPSize; ++k) {
        int min = BinarySearch(fEneUniP[k], fELinearFile3);
        if (fEneUniP[k] == fELinearFile3[min]) {
          fEneTemp.push_back(fEneUniP[k]);
          fSigTemp.push_back(fXLinearFile3[min]);
        } else {
          double sigmaAdd = fXLinearFile3[min] +
                            (fXLinearFile3[min + 1] - fXLinearFile3[min]) * (fEneUniP[k] - fELinearFile3[min]) /
                                (fELinearFile3[min + 1] - fELinearFile3[min]); // linear interpolation
          if (sigmaAdd > 1E-20) {
            fEneTemp.push_back(fEneUniP[k]);
            fSigTemp.push_back(sigmaAdd);
          }
        }
        if (fEneTemp.size() == 1) {
          fEneLocP.push_back(k);
        }
      }
      fSigUniOfP.push_back(fSigTemp);
    } else if ((fMTChargeFlag[1] != -1 && fMT >= 650 && fMT < 700) /* || (fMTChargeFlag [6] != -1 && fMT == 104)
        || (fMTChargeFlag [1] == -1 && fMTChargeFlag [6] == -1 && fMT == 104)*/) {
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab1->GetX(crs));
        fXLinearFile3.push_back(tab1->GetY(crs));
      }
      for (int k = 0, EneUniDSize = fEneUniD.size(); k != EneUniDSize; ++k) {
        int min = BinarySearch(fEneUniD[k], fELinearFile3);
        if (fEneUniD[k] == fELinearFile3[min]) {
          fEneTemp.push_back(fEneUniD[k]);
          fSigTemp.push_back(fXLinearFile3[min]);
        } else {
          double sigmaAdd = fXLinearFile3[min] +
                            (fXLinearFile3[min + 1] - fXLinearFile3[min]) * (fEneUniD[k] - fELinearFile3[min]) /
                                (fELinearFile3[min + 1] - fELinearFile3[min]); // linear interpolation
          fEneTemp.push_back(fEneUniD[k]);
          fSigTemp.push_back(sigmaAdd);
        }
        if (fEneTemp.size() == 1) {
          fEneLocD.push_back(k);
        }
      }
      fSigUniOfD.push_back(fSigTemp);
    } else if ((fMTChargeFlag[2] != -1 && fMT >= 700 && fMT < 750) /* || (fMTChargeFlag [7] != -1 && fMT == 105)
        || (fMTChargeFlag [2] == -1 && fMTChargeFlag [7] == -1 && fMT == 105)*/) {
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab1->GetX(crs));
        fXLinearFile3.push_back(tab1->GetY(crs));
      }
      for (int k = 0, EneUniTSize = fEneUniT.size(); k != EneUniTSize; ++k) {
        int min = BinarySearch(fEneUniT[k], fELinearFile3);
        if (fEneUniT[k] == fELinearFile3[min]) {
          fEneTemp.push_back(fEneUniT[k]);
          fSigTemp.push_back(fXLinearFile3[min]);
        } else {
          double sigmaAdd = fXLinearFile3[min] +
                            (fXLinearFile3[min + 1] - fXLinearFile3[min]) * (fEneUniT[k] - fELinearFile3[min]) /
                                (fELinearFile3[min + 1] - fELinearFile3[min]); // linear interpolation
          fEneTemp.push_back(fEneUniT[k]);
          fSigTemp.push_back(sigmaAdd);
        }
        if (fEneTemp.size() == 1) {
          fEneLocT.push_back(k);
        }
      }
      fSigUniOfT.push_back(fSigTemp);
    } else if ((fMTChargeFlag[3] != -1 && fMT >= 750 && fMT < 800) /* || (fMTChargeFlag [8] != -1 && fMT == 106)
        || (fMTChargeFlag [3] == -1 && fMTChargeFlag [8] == -1 && fMT == 106)*/) {
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab1->GetX(crs));
        fXLinearFile3.push_back(tab1->GetY(crs));
      }
      for (int k = 0, EneUniHe3Size = fEneUniHe3.size(); k != EneUniHe3Size; ++k) {
        int min = BinarySearch(fEneUniHe3[k], fELinearFile3);
        if (fEneUniHe3[k] == fELinearFile3[min]) {
          fEneTemp.push_back(fEneUniHe3[k]);
          fSigTemp.push_back(fXLinearFile3[min]);
        } else {
          double sigmaAdd = fXLinearFile3[min] +
                            (fXLinearFile3[min + 1] - fXLinearFile3[min]) * (fEneUniHe3[k] - fELinearFile3[min]) /
                                (fELinearFile3[min + 1] - fELinearFile3[min]); // linear interpolation
          fEneTemp.push_back(fEneUniHe3[k]);
          fSigTemp.push_back(sigmaAdd);
        }
        if (fEneTemp.size() == 1) {
          fEneLocHe3.push_back(k);
        }
      }
      fSigUniOfHe3.push_back(fSigTemp);
    } else if ((fMTChargeFlag[4] != -1 && fMT >= 800 && fMT < 850) /* || (fMTChargeFlag [9] != -1 && fMT == 107)
        || (fMTChargeFlag [4] == -1 && fMTChargeFlag [9] == -1 && fMT == 107)*/) {
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab1->GetX(crs));
        fXLinearFile3.push_back(tab1->GetY(crs));
      }
      for (int k = 0, EneUniHe4Size = fEneUniHe4.size(); k != EneUniHe4Size; ++k) {
        int min = BinarySearch(fEneUniHe4[k], fELinearFile3);
        if (fEneUniHe4[k] == fELinearFile3[min]) {
          fEneTemp.push_back(fEneUniHe4[k]);
          fSigTemp.push_back(fXLinearFile3[min]);
        } else {
          double sigmaAdd = fXLinearFile3[min] +
                            (fXLinearFile3[min + 1] - fXLinearFile3[min]) * (fEneUniHe4[k] - fELinearFile3[min]) /
                                (fELinearFile3[min + 1] - fELinearFile3[min]); // linear interpolation
          fEneTemp.push_back(fEneUniHe4[k]);
          fSigTemp.push_back(sigmaAdd);
        }
        if (fEneTemp.size() == 1) {
          fEneLocHe4.push_back(k);
        }
      }
      fSigUniOfHe4.push_back(fSigTemp);
    }
    fEneTemp.clear();
    fSigTemp.clear();
    fELinearFile3.clear();
    fXLinearFile3.clear();
  }
  //	Adding charge particle reaction cross-sections
  fSigUniP.resize(fEneUniP.size());
  for (int i = 0, SigUniOfPSize = fSigUniOfP.size(); i != SigUniOfPSize; ++i) {
    int size = fSigUniOfP[i].size();
    for (int j = 0; j < size; j++) {
      fSigUniP[fEneLocP[i] + j] += fSigUniOfP[i][j];
    }
  }
  fSigUniD.resize(fEneUniD.size());
  for (int i = 0, SigUniOfDSize = fSigUniOfD.size(); i != SigUniOfDSize; ++i) {
    int size = fSigUniOfD[i].size();
    for (int j = 0; j < size; j++) {
      fSigUniD[fEneLocD[i] + j] += fSigUniOfD[i][j];
    }
  }
  fSigUniT.resize(fEneUniT.size());
  for (int i = 0, SigUniOfTSize = fSigUniOfT.size(); i != SigUniOfTSize; ++i) {
    int size = fSigUniOfT[i].size();
    for (int j = 0; j < size; j++) {
      fSigUniT[fEneLocT[i] + j] += fSigUniOfT[i][j];
    }
  }
  fSigUniHe3.resize(fEneUniHe3.size());
  for (int i = 0, SigUniOfHe3Size = fSigUniOfHe3.size(); i != SigUniOfHe3Size; ++i) {
    int size = fSigUniOfHe3[i].size();
    for (int j = 0; j < size; j++) {
      fSigUniHe3[fEneLocHe3[i] + j] += fSigUniOfHe3[i][j];
    }
  }
  fSigUniHe4.resize(fEneUniHe4.size());
  for (int i = 0, SigUniOfHe4Size = fSigUniOfHe4.size(); i != SigUniOfHe4Size; ++i) {
    int size = fSigUniOfHe4[i].size();
    for (int j = 0; j < size; j++) {
      fSigUniHe4[fEneLocHe4[i] + j] += fSigUniOfHe4[i][j];
    }
  }
  TIter secIter2(file->GetSections());
  TNudyEndfSec *sec2;
  double QM[5] = {0, 0, 0, 0, 0};
  double QI[5] = {0, 0, 0, 0, 0};
  while ((sec2 = (TNudyEndfSec *)secIter2.Next())) {
    TIter recIter2(sec2->GetRecords());
    TNudyEndfCont *header2 = (TNudyEndfCont *)recIter2.Next();
    int fMT                = sec2->GetMT();
    if (fMT == 600) {
      QM[0] = header2->GetC1();
      QI[0] = header2->GetC2();
    }
    if (fMT == 650) {
      QM[1] = header2->GetC1();
      QI[1] = header2->GetC2();
    }
    if (fMT == 700) {
      QM[2] = header2->GetC1();
      QI[2] = header2->GetC2();
    }
    if (fMT == 750) {
      QM[3] = header2->GetC1();
      QI[3] = header2->GetC2();
    }
    if (fMT == 800) {
      QM[4] = header2->GetC1();
      QI[4] = header2->GetC2();
    }
  }
  if (fMTChargeFlag[0] != -1 && fMTChargeFlag[5] == -1) {
    AddSecFile3(file, QM[0], QI[0], 103, fEneUniP, fSigUniP);
  }
  if (fMTChargeFlag[1] != -1 && fMTChargeFlag[6] == -1) {
    AddSecFile3(file, QM[1], QI[1], 104, fEneUniD, fSigUniD);
  }
  if (fMTChargeFlag[2] != -1 && fMTChargeFlag[7] == -1) {
    AddSecFile3(file, QM[2], QI[2], 105, fEneUniT, fSigUniT);
  }
  if (fMTChargeFlag[3] != -1 && fMTChargeFlag[8] == -1) {
    AddSecFile3(file, QM[3], QI[3], 106, fEneUniHe3, fSigUniHe3);
  }
  if (fMTChargeFlag[4] != -1 && fMTChargeFlag[9] == -1) {
    AddSecFile3(file, QM[4], QI[4], 107, fEneUniHe4, fSigUniHe4);
  }

  // Total charge production for the gas production cross-section
  TIter secIter0(file->GetSections());
  TNudyEndfSec *sec0;
  while ((sec0 = (TNudyEndfSec *)secIter0.Next())) {
    int fMT = sec0->GetMT();
    TIter recIter0(sec0->GetRecords());
    TNudyEndfCont *header0 = (TNudyEndfCont *)recIter0.Next();
    fNR                    = header0->GetN1();
    fNP                    = header0->GetN2();
    TNudyEndfTab1 *tab0    = (TNudyEndfTab1 *)(sec0->GetRecords()->At(0));
    if (fMT == 103 || fMT == 28 || (fMT >= 41 && fMT < 46) || fMT == 111 || fMT == 112 || fMT == 115 || fMT == 116) {
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab0->GetX(crs));
        fXLinearFile3.push_back(tab0->GetY(crs));
      }
      for (int k = 0, EneUniPAllSize = fEneUniPAll.size(); k != EneUniPAllSize; ++k) {
        int min = BinarySearch(fEneUniPAll[k], fELinearFile3);
        if (fEneUniPAll[k] == fELinearFile3[min]) {
          fEneTemp.push_back(fEneUniPAll[k]);
          fSigTemp.push_back(fXLinearFile3[min]);
        } else {
          double sigmaAdd = 0;
          if (fEneUniPAll[k] > fELinearFile3[min]) {
            sigmaAdd = fXLinearFile3[min] +
                       (fXLinearFile3[min + 1] - fXLinearFile3[min]) * (fEneUniPAll[k] - fELinearFile3[min]) /
                           (fELinearFile3[min + 1] - fELinearFile3[min]); // linear interpolation
          }
          fEneTemp.push_back(fEneUniPAll[k]);
          fSigTemp.push_back(sigmaAdd);
        }
        if (fEneTemp.size() == 1) {
          fEneLocPAll.push_back(k);
        }
      }
      fSigUniOfPAll.push_back(fSigTemp);
    } else if (fMT == 11 || fMT == 32 || fMT == 35 || fMT == 114 || fMT == 115 || fMT == 117 || fMT == 104) {
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab0->GetX(crs));
        fXLinearFile3.push_back(tab0->GetY(crs));
      }
      for (int k = 0, EneUniDAllSize = fEneUniDAll.size(); k != EneUniDAllSize; ++k) {
        int min = BinarySearch(fEneUniDAll[k], fELinearFile3);
        if (fEneUniDAll[k] == fELinearFile3[min]) {
          fEneTemp.push_back(fEneUniDAll[k]);
          fSigTemp.push_back(fXLinearFile3[min]);
        } else {
          double sigmaAdd = 0;
          if (fEneUniDAll[k] > fELinearFile3[min]) {
            sigmaAdd = fXLinearFile3[min] +
                       (fXLinearFile3[min + 1] - fXLinearFile3[min]) * (fEneUniDAll[k] - fELinearFile3[min]) /
                           (fELinearFile3[min + 1] - fELinearFile3[min]); // linear interpolation
          }
          fEneTemp.push_back(fEneUniDAll[k]);
          fSigTemp.push_back(sigmaAdd);
        }
        if (fEneTemp.size() == 1) {
          fEneLocDAll.push_back(k);
        }
      }
      fSigUniOfDAll.push_back(fSigTemp);
    } else if (fMT == 33 || fMT == 113 || fMT == 116 || fMT == 105) {
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab0->GetX(crs));
        fXLinearFile3.push_back(tab0->GetY(crs));
      }
      for (int k = 0, EneUniTAllSize = fEneUniTAll.size(); k != EneUniTAllSize; ++k) {
        int min = BinarySearch(fEneUniTAll[k], fELinearFile3);
        if (fEneUniTAll[k] == fELinearFile3[min]) {
          fEneTemp.push_back(fEneUniTAll[k]);
          fSigTemp.push_back(fXLinearFile3[min]);
        } else {
          double sigmaAdd = 0;
          if (fEneUniTAll[k] > fELinearFile3[min]) {
            sigmaAdd = fXLinearFile3[min] +
                       (fXLinearFile3[min + 1] - fXLinearFile3[min]) * (fEneUniTAll[k] - fELinearFile3[min]) /
                           (fELinearFile3[min + 1] - fELinearFile3[min]); // linear interpolation
          }
          fEneTemp.push_back(fEneUniTAll[k]);
          fSigTemp.push_back(sigmaAdd);
        }
        if (fEneTemp.size() == 1) {
          fEneLocTAll.push_back(k);
        }
      }
      fSigUniOfTAll.push_back(fSigTemp);
    } else if (fMT == 34 || fMT == 113 || fMT == 116 || fMT == 106) {
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab0->GetX(crs));
        fXLinearFile3.push_back(tab0->GetY(crs));
      }
      for (int k = 0, EneUniHe3AllSize = fEneUniHe3All.size(); k != EneUniHe3AllSize; ++k) {
        int min = BinarySearch(fEneUniHe3All[k], fELinearFile3);
        if (fEneUniHe3All[k] == fELinearFile3[min]) {
          fEneTemp.push_back(fEneUniHe3All[k]);
          fSigTemp.push_back(fXLinearFile3[min]);
        } else {
          double sigmaAdd = 0;
          if (fEneUniHe3All[k] > fELinearFile3[min]) {
            sigmaAdd = fXLinearFile3[min] +
                       (fXLinearFile3[min + 1] - fXLinearFile3[min]) * (fEneUniHe3All[k] - fELinearFile3[min]) /
                           (fELinearFile3[min + 1] - fELinearFile3[min]); // linear interpolation
          }
          fEneTemp.push_back(fEneUniHe3All[k]);
          fSigTemp.push_back(sigmaAdd);
        }
        if (fEneTemp.size() == 1) {
          fEneLocHe3All.push_back(k);
        }
      }
      fSigUniOfHe3All.push_back(fSigTemp);
    } else if ((fMT >= 22 && fMT < 26) || fMT == 29 || fMT == 30 || fMT == 36 || fMT == 36 || fMT == 45 || fMT == 108 ||
               fMT == 109 || (fMT >= 112 && fMT < 115) || fMT == 117 || fMT == 107) {
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab0->GetX(crs));
        fXLinearFile3.push_back(tab0->GetY(crs));
      }
      for (int k = 0, EneUniHe4AllSize = fEneUniHe4All.size(); k != EneUniHe4AllSize; ++k) {
        int min = BinarySearch(fEneUniHe4All[k], fELinearFile3);
        if (fEneUniHe4All[k] == fELinearFile3[min]) {
          fEneTemp.push_back(fEneUniHe4All[k]);
          fSigTemp.push_back(fXLinearFile3[min]);
        } else {
          double sigmaAdd = 0;
          if (fEneUniHe4All[k] > fELinearFile3[min]) {
            sigmaAdd = fXLinearFile3[min] +
                       (fXLinearFile3[min + 1] - fXLinearFile3[min]) * (fEneUniHe4All[k] - fELinearFile3[min]) /
                           (fELinearFile3[min + 1] - fELinearFile3[min]); // linear interpolation
          }
          fEneTemp.push_back(fEneUniHe4All[k]);
          fSigTemp.push_back(sigmaAdd);
        }
        if (fEneTemp.size() == 1) {
          fEneLocHe4All.push_back(k);
        }
      }
      fSigUniOfHe4All.push_back(fSigTemp);
    }
    fEneTemp.clear();
    fSigTemp.clear();
    fELinearFile3.clear();
    fXLinearFile3.clear();
  }
  //	Adding all charge particle reaction cross-sections
  fSigUniPAll.resize(fEneUniPAll.size());
  for (int i = 0, SigUniOfPAllSize = fSigUniOfPAll.size(); i != SigUniOfPAllSize; ++i) {
    int size = fSigUniOfPAll[i].size();
    for (int j = 0; j < size; j++) {
      fSigUniPAll[fEneLocPAll[i] + j] += fSigUniOfPAll[i][j];
    }
  }
  fSigUniDAll.resize(fEneUniDAll.size());
  for (int i = 0, SigUniOfDAllSize = fSigUniOfDAll.size(); i != SigUniOfDAllSize; ++i) {
    int size = fSigUniOfDAll[i].size();
    for (int j = 0; j < size; j++) {
      fSigUniDAll[fEneLocDAll[i] + j] += fSigUniOfDAll[i][j];
    }
  }
  fSigUniTAll.resize(fEneUniTAll.size());
  for (int i = 0, SigUniOfTAllSize = fSigUniOfTAll.size(); i != SigUniOfTAllSize; ++i) {
    int size = fSigUniOfTAll[i].size();
    for (int j = 0; j < size; j++) {
      fSigUniTAll[fEneLocTAll[i] + j] += fSigUniOfTAll[i][j];
    }
  }
  fSigUniHe3All.resize(fEneUniHe3All.size());
  for (int i = 0, SigUniOfHe3AllSize = fSigUniOfHe3All.size(); i != SigUniOfHe3AllSize; ++i) {
    int size = fSigUniOfHe3All[i].size();
    for (int j = 0; j < size; j++) {
      fSigUniHe3All[fEneLocHe3All[i] + j] += fSigUniOfHe3All[i][j];
    }
  }
  fSigUniHe4All.resize(fEneUniHe4All.size());
  for (int i = 0, SigUniOfHe4AllSize = fSigUniOfHe4All.size(); i != SigUniOfHe4AllSize; ++i) {
    int size = fSigUniOfHe4All[i].size();
    for (int j = 0; j < size; j++) {
      fSigUniHe4All[fEneLocHe4All[i] + j] += fSigUniOfHe4All[i][j];
    }
  }
  // Adding new fMF for charge particles (p, d, t, He3, He4)

  if ((fMTChargeFlag[0] == -1 || fMTChargeFlag[5] == -1)) {
    AddSecFile3(file, 0, 0, 203, fEneUniPAll, fSigUniPAll);
  }
  if ((fMTChargeFlag[1] == -1 || fMTChargeFlag[6] == -1)) {
    AddSecFile3(file, 0, 0, 204, fEneUniDAll, fSigUniDAll);
  }
  if ((fMTChargeFlag[2] == -1 || fMTChargeFlag[7] == -1)) {
    AddSecFile3(file, 0, 0, 205, fEneUniTAll, fSigUniTAll);
  }
  if ((fMTChargeFlag[3] == -1 || fMTChargeFlag[8] == -1)) {
    AddSecFile3(file, 0, 0, 206, fEneUniHe3All, fSigUniHe3All);
  }
  if ((fMTChargeFlag[4] == -1 || fMTChargeFlag[9] == -1)) {
    AddSecFile3(file, 0, 0, 207, fEneUniHe4All, fSigUniHe4All);
  }

  TIter secIter20(file->GetSections());
  TNudyEndfSec *sec20;
  while ((sec20 = (TNudyEndfSec *)secIter20.Next())) {
    MtNumbers.push_back(sec20->GetMT());
    //       std::cout<<"fMT push "<< sec20->GetMT() << std::endl;
    TIter recIter21(sec20->GetRecords());
    TNudyEndfCont *header21 = (TNudyEndfCont *)recIter21.Next();
    fQvalueTemp.push_back(header21->GetC2());
    fQvalueTemp.push_back(sec20->GetMT());
  }

  /*
  TIter secIterx(file->GetSections());
  TNudyEndfSec *secx;
  while ((secx = (TNudyEndfSec *)secIterx.Next())) {
    int fMT = secx->GetMT();
    TIter recItery(secx->GetRecords());
    TNudyEndfCont *headerx = (TNudyEndfCont *)recItery.Next();
    fNR = headerx->GetN1();
    fNP = headerx->GetN2();
    std::cout<<"fMT "<< fMT <<" fNP "<< fNP <<"  "<< headerx->GetC1()<<"  "<< headerx->GetC2()<<"  "<<
  headerx->GetL1()<<"  "<< headerx->GetL2()
    <<"  "<< headerx->GetN1()<<"  "<< headerx->GetN2()<< std::endl;
    TNudyEndfTab1 *tabx = (TNudyEndfTab1 *)(secx->GetRecords()->At(0));
       for (int crs = 0; crs < fNP; crs++) {
        std::cout<<"MTX "<< fMT<<" fNP "<< fNP <<"  "<< tabx->GetX(crs) <<"  "<< tabx->GetY(crs) << std::endl;
      }
  }  */
  // fill fSigmaMts
  TIter secIter22(file->GetSections());
  TNudyEndfSec *sec22;
  while ((sec22 = (TNudyEndfSec *)secIter22.Next())) {
    TIter recIter23(sec22->GetRecords());
    TNudyEndfCont *header23 = (TNudyEndfCont *)recIter23.Next();
    if ((sec22->GetMT() >= 103 && sec22->GetMT() < 108) || (sec22->GetMT() >= 203 && sec22->GetMT() < 208)) {
      if (fMTChargeFlag[0] == -1 && sec22->GetMT() == 103) continue;
      if (fMTChargeFlag[1] == -1 && sec22->GetMT() == 104) continue;
      if (fMTChargeFlag[2] == -1 && sec22->GetMT() == 105) continue;
      if (fMTChargeFlag[3] == -1 && sec22->GetMT() == 106) continue;
      if (fMTChargeFlag[4] == -1 && sec22->GetMT() == 107) continue;
      fNP                  = header23->GetN2();
      TNudyEndfTab1 *tab23 = (TNudyEndfTab1 *)(sec22->GetRecords()->At(0));
      for (int crs = 0; crs < fNP; crs++) {
        fELinearFile3.push_back(tab23->GetX(crs));
        fXLinearFile3.push_back(tab23->GetY(crs));
        //  	  std::cout <<" fMT "<< sec22->GetMT()<<" energy "<< fELinearFile3 [crs] <<" fSigma "<< fXLinearFile3
        //  [crs] << std::endl;
      }
      fSigmaMts.insert(std::end(fSigmaMts), std::begin(fELinearFile3), std::end(fELinearFile3));
      fSigmaMts.insert(std::end(fSigmaMts), std::begin(fXLinearFile3), std::end(fXLinearFile3));
      fELinearFile3.clear();
      fXLinearFile3.clear();
      fSigmaOfMts.push_back(fSigmaMts);
      fSigmaMts.clear();
    }
  }
  fQvalue.push_back(fQvalueTemp);
  fQvalueTemp.clear();
  // empty vectors
  fSigUniOfP.clear();
  fEneLocP.clear();
  fSigUniOfD.clear();
  fEneLocD.clear();
  fSigUniOfT.clear();
  fEneLocT.clear();
  fSigUniOfHe3.clear();
  fEneLocHe3.clear();
  fSigUniOfHe4.clear();
  fEneLocHe4.clear();
  fSigUniOfPAll.clear();
  fEneLocPAll.clear();
  fSigUniOfDAll.clear();
  fEneLocDAll.clear();
  fSigUniOfTAll.clear();
  fEneLocTAll.clear();
  fSigUniOfHe3All.clear();
  fEneLocHe3All.clear();
  fSigUniOfHe4All.clear();
  fEneLocHe4All.clear();
}
int TNudyEndfSigma::BinarySearch(double x1, rowd &x2)
{
  int min  = 0;
  int size = x2.size();
  min      = 0;
  int max  = size - 1;
  int mid  = 0;
  if (x1 <= x2[min])
    min = 0;
  else if (x1 >= x2[max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (x1 < x2[mid])
        max = mid;
      else
        min = mid;
    }
  }
  return min;
}
// Adding new fMF for charge particles
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::AddSecFile3(TNudyEndfFile *file1, double a, double b, int MF2Add, rowd &x1, rowd &x2)
{

  TNudyEndfSec *sec3 = new TNudyEndfSec(fMAT, 3, MF2Add, ZA, AWRI, 0, 0, 0, 0);
  //  TNudyEndfCont *secCont = new TNudyEndfCont (a,b,0,0,1,(int)x1.size());
  //  sec3->Add(secCont);
  TNudyEndfTab1 *sec3Tab1 = new TNudyEndfTab1();
  sec3Tab1->SetCont(a, b, 0, 0, 1, (int)x1.size());
  sec3Tab1->SetNBT((int)x1.size(), 0);
  sec3Tab1->SetINT(2, 0);
  for (int j = 0, x1Size = x1.size(); j != x1Size; ++j) {
    sec3Tab1->SetX(x1[j], j);
    sec3Tab1->SetY(x2[j], j);
    //     std::cout << "x2 " << x2[j] << std::endl;
  }
  sec3->Add(sec3Tab1);
  //   TNudyEndfCont *secCont1 = new TNudyEndfCont (0,0,0,0,0,0);
  //   secCont1->SetContMF(fMAT,MF2Add,0);
  //   sec3->Add(secCont1);
  file1->Add(sec3);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::K_wnum(double x)
{
  double k;
  k = fFactor_k * sqrt(std::fabs(x));
  return k;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::GetRho(double x1, int lVal)
{
  if (!NAPS && !NRO)
    return fFactor_k * sqrt(std::fabs(x1)) *
           fRad_a; // Use fRad_a in the penetrabilities Pl and shift factors Sl , and fAp
                   // in the hard-sphere phase shifts φl
  if (!NAPS && NRO)
    return fFactor_k * sqrt(std::fabs(x1)) * fRad_a; // Use fRad_a in the penetrabilities Pl and shift factors Sl,fAp(E)
                                                     // in the hard-sphere phase shifts φl
  if (NAPS && !NRO)
    return fFactor_k * sqrt(std::fabs(x1)) *
           fApl[lVal]; // use fAp in the penetrabilities and shift factor as well as in the phase shifts
  if (NAPS && NRO)
    return fFactor_k * sqrt(std::fabs(x1)) * fApl[lVal]; // read fAp(E) and use it in all three places, Pl , Sl , φl
  return 0;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::GetRhoC(double x, int isDiff, int lVal)
{
  if (!isDiff)
    return fFactor_k * sqrt(std::fabs(x)) * fApl[lVal]; // read fAp(E) and use it in all three places, Pl , Sl , φl
  if (isDiff)
    return fFactor_k * sqrt(std::fabs(x)) * fApl[lVal]; // read fAp(E) and use it in all three places, Pl , Sl , φl
  return 0;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::CalcPhi(double x1, int l)
{
  double x = GetRhoC(x1, 0, l);
  switch (l) {
  case 0:
    return x;
  case 1:
    return (x - atan(x));
  case 2:
    return (x - atan(3.0 * x / (3.0 - x2(x))));
  case 3:
    return (x - atan(x * (15.0 - x * x) / (15.0 - 6.0 * x2(x))));
  case 4:
    return (x - atan(x * (105.0 - 10.0 * x2(x)) / (105.0 - 45.0 * x2(x) + x4(x))));
  }
  return 0;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::CalcShift(double x1, int l)
{
  double x = GetRho(x1, l);
  switch (l) {
  case 0:
    return 0.0;
  case 1:
    return (-1.0 / (1.0 + x * x));
  case 2:
    return (-(18.0 + 3.0 * x2(x)) / Fac2(x));
  case 3:
    return (-(675.0 + 90.0 * x2(x) + 6.0 * x4(x)) / (x6(x) + 6.0 * x4(x) + 45.0 * x2(x) + 225.0));
  case 4:
    return (-(44100.0 + 4725.0 * x2(x) + 270.0 * x4(x) + 10.0 * x6(x)) /
            (x8(x) + 10.0 * x6(x) + 135.0 * x4(x) + 1575.0 * x2(x) + 11025.0));
  }
  return 0;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::CalcPene(double x1, int l)
{
  double x = GetRho(x1, l);
  switch (l) {
  case 0:
    return x;
  case 1:
    return (x * x * x) / (1.0 + x * x);
  case 2:
    return (x5(x) / Fac2(x));
  case 3:
    return ((x * x6(x)) / (225.0 + 45.0 * x2(x) + 6.0 * x4(x) + x6(x)));
  case 4:
    return ((x * x8(x)) / (x8(x) + 10.0 * x6(x) + 135.0 * x4(x) + 1575.0 * x2(x) + 11025.0));
  }
  return 0;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::GetERP(double x, int r, int lVal)
{
  double er = fEr[r];
  if (lVal == 0) return er;
  return er + (fShiftEr[r] - CalcShift(x, lVal)) / (2.0 * fPhiEr[r]) * Gamma_nrE(std::fabs(er), r, lVal);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::Gamma_nrE(double x, int ii, int lval)
{
  return CalcPene(x, lval) * fGamma_n[ii];
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::Gamma_xrE(int ii, int lrx)
{
  if (!lrx) {
    return 0;
  } else {
    return fGamma_r[ii] - (fGamma_n[ii] + fGamma_g[ii] + fGamma_f[ii]);
  }
}
//------------------------------------------------------------------------------------------------------

double TNudyEndfSigma::Gamma_rE(double x, int ii, int lval, int lrx)
{
  return Gamma_nrE(x, ii, lval) + Gamma_xrE(ii, lrx) + fGamma_g[ii] + fGamma_f[ii];
}
// MC2 formalism is inspired from M. Kalki
//------------------------------------------------------------------------------------------------------
int TNudyEndfSigma::WidthFluctuation(double GNR, double gx, double gg, double gf, int jval)
{
  double XX[10][5] = {{3.0013465E-03, 1.3219203E-02, 1.0004488E-03, 1.3219203E-02, 1.0E+0},
                      {7.8592886E-02, 7.2349624E-02, 2.6197629E-02, 7.2349624E-02, 0.0E+0},
                      {4.3282415E-01, 1.9089473E-01, 1.4427472E-01, 1.9089473E-01, 0.0E+0},
                      {1.3345267E+00, 3.9528842E-01, 4.4484223E-01, 3.9528842E-01, 0.0E+0},
                      {3.0481846E+00, 7.4083443E-01, 1.0160615E+00, 7.4083443E-01, 0.0E+0},
                      {5.8263198E+00, 1.3498293E+00, 1.9421066E+00, 1.3498293E+00, 0.0E+0},
                      {9.9452656E+00, 2.5297983E+00, 3.3150885E+00, 2.5297983E+00, 0.0E+0},
                      {1.5782128E+01, 5.2384894E+00, 5.2607092E+00, 5.2384894E+00, 0.0E+0},
                      {2.3996824E+01, 1.3821772E+01, 7.9989414E+00, 1.3821772E+01, 0.0E+0},
                      {3.6216208E+01, 7.5647525E+01, 1.2072069E+01, 7.5647525E+01, 0.0E+0}};

  double WW[10][5] = {{1.1120413E-01, 3.3773418E-02, 3.3376214E-04, 1.7623788E-03, 1.0E+0},
                      {2.3546798E-01, 7.9932171E-02, 1.8506108E-02, 2.1517749E-02, 0.0E+0},
                      {2.8440987E-01, 1.2835937E-01, 1.2309946E-01, 8.0979849E-02, 0.0E+0},
                      {2.2419127E-01, 1.7652616E-01, 2.9918923E-01, 1.8797998E-01, 0.0E+0},
                      {1.0967668E-01, 2.1347043E-01, 3.3431475E-01, 3.0156335E-01, 0.0E+0},
                      {3.0493789E-02, 2.1154965E-01, 1.7766657E-01, 2.9616091E-01, 0.0E+0},
                      {4.2930874E-03, 1.3365186E-01, 4.2695894E-02, 1.0775649E-01, 0.0E+0},
                      {2.5827047E-04, 2.2630659E-02, 4.0760575E-03, 2.5171914E-03, 0.0E+0},
                      {4.9031965E-06, 1.6313638E-05, 1.1766115E-04, 8.9630388E-10, 0.0E+0},
                      {1.4079206E-08, 0.0000000E+00, 5.0989546E-07, 0.0000000E+00, 0.0E+0}};

  fRN     = 0.0;
  fRG     = 0.0;
  fRF     = 0.0;
  fRX     = 0.0;
  int mux = (int)fAmux[jval];
  int mun = (int)fAmun[jval];
  int muf = (int)fAmuf[jval];
  if (GNR <= 0.0 || gg <= 0.0) {
    return 0;
  }
  if (mun < 1 || mun > 4) mun = 5;
  if (muf < 1 || muf > 4) muf = 5;
  if (mux < 1 || mux > 4) mux = 5;
  if (gf < 0.0) {
    return 0;
  }
  if (gf == 0.0 && gx == 0.0) {
    for (int j = 0; j < 10; j++) {
      double XJ  = XX[j][mun - 1];
      double fak = XJ * WW[j][mun - 1] / (GNR * XJ + gg);
      fRN += fak * XJ;
      fRG += fak;
    } // for loop
    return 0;
  } // if
  if (gf == 0.0 && gx > 0.0) {
    for (int j = 0; j < 10; j++) {
      double XJ   = XX[j][mun - 1];
      double WJXJ = WW[j][mun - 1] * XJ;
      double EFFJ = GNR * XJ + gg;
      for (int k = 0; k < 10; k++) {
        double XK = XX[k][mux - 1];
        double fk = WW[k][mux - 1] * WJXJ / (EFFJ + gx * XK);
        fRN += XJ * fk;
        fRG += fk;
        fRX += XK * fk;
      } // 2 for loop
    }   // 1 for loop
  }     // if
  if (gf > 0.0 && gx == 0.0) {
    for (int j = 0; j < 10; j++) {
      double XJ   = XX[j][mun - 1];
      double WJXJ = WW[j][mun - 1] * XJ;
      double EFFJ = GNR * XJ + gg;
      for (int k = 0; k < 10; k++) {
        double fk = WW[k][muf - 1] * WJXJ / (EFFJ + gf * XX[k][muf - 1]);
        fRN += XJ * fk;
        fRG += fk;
        fRF += XX[k][muf - 1] * fk;
      } // 2 for loop
    }   // 1 for loop
  }     // if
  if (gf > 0.0 && gx > 0.0) {
    for (int j = 0; j < 10; j++) {
      double XJ   = XX[j][mun - 1];
      double WJXJ = WW[j][mun - 1] * XJ;
      double EFFJ = GNR * XJ + gg;
      for (int k = 0; k < 10; k++) {
        double XK     = XX[k][muf - 1];
        double WKWJXJ = WW[k][muf - 1] * WJXJ;
        double EFFJK  = gf * XK + EFFJ;
        for (int l = 0; l < 10; l++) {
          double XL = XX[l][mux - 1];
          double fk = WW[l][mux - 1] * WKWJXJ / (EFFJK + gx * XL);
          fRN += XJ * fk;
          fRG += fk;
          fRX += fk * XL;
          fRF += XK * fk;
        } // 3 for loop
      }   // 2 for loop
    }     // 1 for loop
  }       // if
  return 0;
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::InverseMatrix()
{
  double A[3][3], C[3][3], C2[3][3], C3[3][3], AI[3][3], AIB[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A[i][j] = fR[i][j];
    }
    A[i][i] = A[i][i] + 1.0;
  }
  double det1 = A[1][1] * A[2][2] - A[1][2] * A[2][1];
  double det2 = A[1][2] * A[2][0] - A[1][0] * A[2][2];
  double det3 = A[1][0] * A[2][1] - A[1][1] * A[2][0];
  double det  = 1./(A[0][0] * det1 + A[0][1] * det2 + A[0][2] * det3);
  AI[0][0]    = det1 * det;
  AI[1][0]    = det2 * det;
  AI[2][0]    = det3 * det;
  AI[0][1]    = AI[1][0];
  AI[1][1]    = (A[0][0] * A[2][2] - A[2][0] * A[0][2]) * det;
  AI[2][1]    = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) * det;
  AI[0][2]    = AI[2][0];
  AI[1][2]    = AI[2][1];
  AI[2][2]    = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * det;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double sum1 = 0.0;
      for (int k = 0; k < 3; k++) {
        sum1 += AI[i][k] * fS[k][j];
      }
      AIB[i][j] = sum1;
    }
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double sum1 = 0.0;
      for (int k = 0; k < 3; k++) {
        sum1 += fS[i][k] * AIB[k][j];
      }
      C[i][j]  = A[i][j] + sum1;
      C2[i][j] = fR[i][j] + sum1;
    }
  }
  det1     = C[1][1] * C[2][2] - C[1][2] * C[2][1];
  det2     = C[1][2] * C[2][0] - C[1][0] * C[2][2];
  det3     = C[1][0] * C[2][1] - C[1][1] * C[2][0];
  det      = 1./(C[0][0] * det1 + C[0][1] * det2 + C[0][2] * det3);
  C3[0][0] = det1 * det;
  C3[1][0] = det2 * det;
  C3[2][0] = det3 * det;
  C3[0][1] = C3[1][0];
  C3[1][1] = (C[0][0] * C[2][2] - C[2][0] * C[0][2]) * det;
  C3[2][1] = (C[0][1] * C[2][0] - C[0][0] * C[2][1]) * det;
  C3[0][2] = C3[2][0];
  C3[1][2] = C3[2][1];
  C3[2][2] = (C[0][0] * C[1][1] - C[0][1] * C[1][0]) * det;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double sum1 = 0.0;
      for (int k = 0; k < 3; k++) {
        sum1 += C3[i][k] * C2[k][j];
      }
      C[i][j]   = -sum1;
      fRI[i][j] = -sum1;
    }
    C[i][i] = C[i][i] + 1.0;
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double sum1 = 0.0;
      for (int k = 0; k < 3; k++) {
        sum1 += AIB[i][k] * C[k][j];
      }
      fSI[i][j] = -sum1;
    }
  }
}
// The Reich Moore formalism is inspired from D.E. Cullen
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::GetSigmaRMP(double x, double &siga, double &sigb, double &sigc)
{
  std::vector<double> row;
  int nrsl, nrs = 0, jVal, nrstemp;
  double fac2;
  double phil, sfil, cfil, sfil2, s2fil;
  double deno = 0.0;
  double ERP, GNR;
  double sumCrsElastic = 0.0, sumCrsCapture = 0.0, sumCrsFission = 0.0;
  double GR;
  double reduced_n, reduced_g, reduced_fa, reduced_fb;
  x             = std::fabs(x);
  double pibyk2 = PI / x2(K_wnum(x));
  double rt11   = 0.0,   st11 = 0.0;
  double sigrm1 = 0.0, sigab1 = 0.0, sigf1 = 0.0, sigel = 0.0;
  double gj;
  for (int lcnt = 0; lcnt < NLS; lcnt++) {
    fac2     = 0.0;
    int lval = l[lcnt];
    nrsl     = fNRS[lval];
    phil     = CalcPhi(x, l[lcnt]);
    sfil     = sin(phil);
    cfil     = cos(phil);
    sfil2    = sfil * sfil;
    s2fil    = 2. * sfil * cfil;
    fac2     = 4.0 * fMisGj[l[lcnt]] * sfil2 * pibyk2;
    sigrm1   = 0.0;
    sigel    = 0.0, sigab1 = 0.0, sigf1 = 0.0;
    jVal      = 0;
    nrstemp   = nrs;
    int nrend = nrs + nrsl;
    do {
      rt11 = 0.0;
      st11 = 0.0;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          fR[i][j]  = 0.0;
          fS[i][j]  = 0.0;
          fRI[i][j] = 0.0;
          fSI[i][j] = 0.0;
        }
      }
      for (int r = nrs; r < nrend; r++) {
        if (fMissingJ[l[lcnt]][jVal] != (J[r])) {
          nrs = r;
          break;
        }
        // fGamma_n,r term
        reduced_n     = 0.5 * Gamma_nrE(x, r, l[lcnt]);
        reduced_g     = 0.5 * std::fabs(fGamma_g[r]); //
        ERP           = fEr[r];
        GR            = reduced_g;
        GNR           = reduced_n;
        double DE     = ERP - x;
        deno          = DE * DE + GR * GR;
        double de2    = DE / deno;
        double gamma2 = GR / deno;
	// average fission widths are not given
        if (LFW == 0) { 
          rt11       += gamma2 * GNR;
          st11       += de2 * GNR;
        } else {
	  // average fission widths are given
          reduced_fa  = fGamma_fasq[r]; 
          reduced_fb  = fGamma_fbsq[r]; 
          double GR2  = sqrt(GNR);
          double GFA  = reduced_fa * reduced_fa;
          double GFB  = reduced_fb * reduced_fb;
          double GFAN = reduced_fa * GR2;
          double GFBN = reduced_fb * GR2;
          double GFAB = reduced_fa * reduced_fb;
          fR[0][0]   += gamma2 * GNR;
          fR[0][1]   += gamma2 * GFAN;
          fR[0][2]   += gamma2 * GFBN;
          fR[1][1]   += gamma2 * GFA;
          fR[1][2]   += gamma2 * GFAB;
          fR[2][2]   += gamma2 * GFB;
          fS[0][0]   += de2 * GNR;
          fS[0][1]   += de2 * GFAN;
          fS[0][2]   += de2 * GFBN;
          fS[1][1]   += de2 * GFA;
          fS[1][2]   += de2 * GFAB;
          fS[2][2]   += de2 * GFB;
          fR[1][0]    = fR[0][1];
          fS[1][0]    = fS[0][1];
          fR[2][0]    = fR[0][2];
          fS[2][0]    = fS[0][2];
          fR[2][1]    = fR[1][2];
          fS[2][1]    = fS[1][2];
        }
      }
      if (LFW == 0) {
        gj          = (2. * std::fabs(fMissingJ[l[lcnt]][jVal]) + 1.) / (4. * fSpi + 2.0);
        double det  = 1./((rt11 + 1.0) * (rt11 + 1.0) + st11 * st11);
        double SI11 = -st11 * det;
        double RI11 = -(rt11 * (rt11 + 1.0) + st11 * st11) * det;
        sigrm1 += -4. * gj * (RI11 + (RI11 * RI11 + SI11 * SI11));
        sigel +=
            gj * ((2.0 * sfil2 + 2. * RI11) * (2.0 * sfil2 + 2. * RI11) + (s2fil + 2.0 * SI11) * (s2fil + 2.0 * SI11));
        sigf1 = 0.0;
      } else {
        InverseMatrix();
        gj          = (2. * std::fabs(fMissingJ[l[lcnt]][jVal]) + 1.) / (4. * fSpi + 2.0);
        double SI11 = fSI[0][0];
        double RI11 = fRI[0][0];
        double SI12 = fSI[0][1];
        double RI12 = fRI[0][1];
        double SI13 = fSI[0][2];
        double RI13 = fRI[0][2];
        sigab1 += -4.0 * gj * (RI11 + (RI11 * RI11 + SI11 * SI11));
        sigel  +=
            gj * ((2.0 * sfil2 + 2. * RI11) * (2.0 * sfil2 + 2. * RI11) + (s2fil + 2.0 * SI11) * (s2fil + 2.0 * SI11));
        sigf1  += 4.0 * gj * (RI12 * RI12 + RI13 * RI13 + SI12 * SI12 + SI13 * SI13);
        sigrm1  = sigab1 - sigf1;
      }
      jVal += 1;
    } while (jVal < fNJValue[l[lcnt]]);
    nrs = nrstemp;
    nrs += nrsl;
    sumCrsElastic += pibyk2 * (sigel) + fac2;
    sumCrsCapture += sigrm1 * pibyk2;
    sumCrsFission += sigf1 * pibyk2;
  }
  siga = sumCrsElastic;
  sigb = sumCrsCapture;
  sigc = sumCrsFission;
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::GetSigma(int lrfa, double x, double &siga, double &sigb, double &sigc)
{
  std::vector<double> row;
  int nrsl, nrs = 0, lVal = 0, jVal, nrstemp;
  double fac1 = 0.0;
  double phil, sfil, cfil, sfil2, s2fil;
  double deno = 0.0;
  double ERP, ERPP1, GNR;
  double fachNumElastic = 0.0, fachNumCapture = 0.0, fachNumFission = 0.0;
  double facElastic = 0.0, facCapture = 0.0, facFission = 0.0;
  double sumElastic = 0.0, sumCapture = 0.0, sumFission = 0.0;
  double sumCrsElastic = 0.0, sumCrsCapture = 0.0, sumCrsFission = 0.0;
  double GNRS, GR, GS;
  double fachNumElasticG, fachNumElasticH, fachNumElasticM = 0.0;
  double facElasticM = 0.0, sumElasticM = 0.0;
  x             = std::fabs(x);
  double pibyk2;
  switch (LRU) {
  case 1:
    switch (lrfa) {
    // SLBW--------
    case 1: {
      pibyk2 = PI / x2(K_wnum(x));
      for (int lcnt = 0; lcnt < NLS; lcnt++) {
        sumElastic  = 0.0;
        sumCapture  = 0.0;
        sumFission  = 0.0;
        sumElasticM = 0.0;
        lVal        = l[lcnt];
        nrsl        = fNRS[lVal];
        phil        = CalcPhi(x, lVal);
        sfil        = sin(phil);
        cfil        = cos(phil);
        sfil2       = sfil * sfil;
        // cfil2 = cfil*cfil;
        s2fil = 2. * sfil * cfil;
        // c2fil = 1. - 2.*sfil2;
        fac1      = 4.0 * (2.0 * lVal + 1.0) * sfil2;
        jVal      = 0;
        nrstemp   = nrs;
        int nrend = nrs + nrsl;
        do {
          facElastic = 0.0;
          facCapture = 0.0;
          facFission = 0.0;
          for (int r = nrs; r < nrend; r++) {
            if (fMissingJ[l[lcnt]][jVal] != (J[r])) {
              nrs = r;
              break;
            }
            GNR            = Gamma_nrE(x, r, lVal); // Gamma_nr(E)
            GR             = Gamma_rE(x, r, lVal, fLrx);
            ERP            = GetERP(x, r, lVal);
            double xerp    = x - ERP;
            deno           = 1./(xerp * xerp + 0.25 * GR * GR);
            fachNumElastic = GNR * ((GNR - 2.0 * GR * sfil2 + 2.0 * xerp * s2fil) * deno);
            fachNumCapture = (GNR * fGamma_g[r] * deno);
            fachNumFission = (GNR * fGamma_f[r] * deno);
            facElastic += fGJ[r] * fachNumElastic;
            facCapture += fGJ[r] * fachNumCapture;
            facFission += fGJ[r] * fachNumFission;
          }
          sumElastic += facElastic;
          sumCapture += facCapture;
          sumFission += facFission;
          jVal += 1;
        } while (jVal < fNJValue[l[lcnt]]);
        nrs = nrstemp;
        nrs += nrsl;
        sumCrsElastic += pibyk2 * (fac1 + sumElastic);
        sumCrsCapture += pibyk2 * (sumCapture);
        sumCrsFission += pibyk2 * (sumFission);
      }
      siga = sumCrsElastic;
      sigb = sumCrsCapture;
      sigc = sumCrsFission;
    } break;
    case 2: {
      pibyk2 = PI / x2(K_wnum(x));
      for (int lcnt = 0; lcnt < NLS; lcnt++) {
        sumElastic  = 0.0;
        sumCapture  = 0.0;
        sumFission  = 0.0;
        sumElasticM = 0.0;
        lVal        = l[lcnt];
        nrsl        = fNRS[lVal];
        phil        = CalcPhi(x, lVal);
        sfil        = sin(phil);
        cfil        = cos(phil);
        sfil2       = sfil * sfil;
        // cfil2 = cfil * cfil;
        s2fil = 2. * sfil * cfil;
        // c2fil = 1. - 2.*sfil2;
        fac1      = 4.0 * (2.0 * lVal + 1.0) * sfil2;
        jVal      = 0;
        nrstemp   = nrs;
        int nrend = nrs + nrsl;
        do {
          facElastic  = 0.0;
          facCapture  = 0.0;
          facFission  = 0.0;
          facElasticM = 0.0;
          for (int r = nrs; r < nrend; r++) {
            fachNumElastic  = 0.0;
            fachNumCapture  = 0.0;
            fachNumFission  = 0.0;
            fachNumElasticM = 0.0;
            if (fMissingJ[l[lcnt]][jVal] != (J[r])) {
              nrs = r;
              break;
            }
            GNR            = Gamma_nrE(x, r, lVal); // Gamma_nr(E)
            GR             = Gamma_rE(x, r, lVal, fLrx);
            ERP            = (GetERP(x, r, lVal));
            double xerp    = x - ERP;
            deno           = 1./(xerp * xerp + 0.25 * GR * GR);
            fachNumElastic = GNR * ((GNR - 2.0 * GR * sfil2 + 2.0 * xerp * s2fil) * deno);
            fachNumCapture = (GNR * fGamma_g[r] * deno);
            fachNumFission = (GNR * fGamma_f[r] * deno);
            facElastic += fGJ[r] * fachNumElastic;
            facCapture += fGJ[r] * fachNumCapture;
            facFission += fGJ[r] * fachNumFission;
            fachNumElasticG = 0.0;
            fachNumElasticH = 0.0;
            for (int rs = nrs; rs < nrend; rs++) {
              if (fMissingJ[l[lcnt]][jVal] != J[rs]) continue;
              GNRS  = Gamma_nrE(x, rs, lVal);
              GS    = Gamma_rE(x, rs, lVal, fLrx);
              ERPP1 = (GetERP(x, rs, lVal));
              if (r != rs) {
                double grgs      = GR + GS;
                double erep1     = ERP - ERPP1;
                double deno1     = 1./(erep1 * erep1 + 0.25 * grgs * grgs);
                fachNumElasticG += 0.5 * GNRS * GNR * grgs *deno1;
                fachNumElasticH += GNRS * GNR * erep1 * deno1;
              }
            }
            fachNumElasticM      = (fachNumElasticG * GR + 2. * fachNumElasticH * xerp) * deno;
            facElasticM         += fGJ[r] * fachNumElasticM;
          }
          sumElastic += facElastic;
          sumCapture += facCapture;
          sumFission += facFission;
          sumElasticM += facElasticM;
          jVal += 1;
        } while (jVal < fNJValue[l[lcnt]]);
        nrs = nrstemp;
        nrs += nrsl;
        sumCrsElastic += pibyk2 * (fac1 + sumElastic + sumElasticM);
        sumCrsCapture += pibyk2 * (sumCapture);
        sumCrsFission += pibyk2 * (sumFission);
      }
      siga = sumCrsElastic;
      sigb = sumCrsCapture;
      sigc = sumCrsFission;
    } break;
    case 3: {
      GetSigmaRMP(x, siga, sigb, sigc);
    } break;
    // Adler-Adler
    case 4: {
    } break;
    // R-Matrix
    case 7: {
    } break;
    }
    break;
  case 2: {
    int nrsj = 0;
    nrs      = 0;
    pibyk2   = PI / x2(K_wnum(x));
    for (int lcnt = 0; lcnt < NLS2; lcnt++) {
      sumElastic  = 0.0;
      sumCapture  = 0.0;
      sumFission  = 0.0;
      sumElasticM = 0.0;
      lVal        = l[lcnt];
      phil        = CalcPhi(x, lVal);
      sfil        = sin(phil);
      cfil        = cos(phil);
      sfil2       = sfil * sfil;
      s2fil       = 2. * sfil * cfil;
      fac1        = 4.0 * (2.0 * lVal + 1.0) * sfil2;
      NJS         = fNRJ[lcnt];
      for (int jVal = 0; jVal < NJS; jVal++) {
        nrsl        = fJSM[nrsj];
        facElastic  = 0.0;
        facCapture  = 0.0;
        facFission  = 0.0;
        facElasticM = 0.0;
        int min     = 0;
        int max     = nrsl - 1;
        int mid     = 0;
        if (x <= fEs[min])
          min = 0;
        else if (x >= fEs[max])
          min = max - 1;
        else {
          while (max - min > 1) {
            mid = (min + max) / 2;
            if (x < fEs[mid])
              max = mid;
            else
              min = mid;
          }
        }
        double gnr = 0, gx = 0, gg = 0, gf = 0, dr = 0;
        double gnop[2]       = {fGNO[nrs + min], fGNO[nrs + min + 1]};
        double gxp[2]        = {fGX[nrs + min], fGX[nrs + min + 1]};
        double ggp[2]        = {fGG[nrs + min], fGG[nrs + min + 1]};
        double gfp[2]        = {fGF[nrs + min], fGF[nrs + min + 1]};
        double drp[2]        = {fD[nrs + min], fD[nrs + min + 1]};
        double esp[2]        = {fEs[nrs + min], fEs[nrs + min + 1]};
        gnr                  = TNudyCore::Instance()->InterpolateScale(esp, gnop, INT, x);
        gx                   = TNudyCore::Instance()->InterpolateScale(esp, gxp, INT, x);
        gg                   = TNudyCore::Instance()->InterpolateScale(esp, ggp, INT, x);
        gf                   = TNudyCore::Instance()->InterpolateScale(esp, gfp, INT, x);
        dr                   = TNudyCore::Instance()->InterpolateScale(esp, drp, INT, x);
        if (gnr < 1E-19) gnr = 0.0;
        if (gg < 1E-19) gg   = 0.0;
        if (gf < 1E-19) gf   = 0.0;
        if (gx < 1E-19) gx   = 0.0;
        if (fGNO[nrs + min] == 0.0 && INT > 3) {
          gnr = gnop[0] + (gnop[1] - gnop[0]) * (x - esp[0]) / (esp[1] - esp[0]);
        }
        if (fGG[nrs + min] == 0.0 && INT > 3) {
          gg = ggp[0] + (ggp[1] - ggp[0]) * (x - esp[0]) / (esp[1] - esp[0]);
        }
        if (fGX[nrs + min] == 0.0 && INT > 3) {
          gx = gxp[0] + (gxp[1] - gxp[0]) * (x - esp[0]) / (esp[1] - esp[0]);
        }
        if (fGF[nrs + min] == 0.0 && INT > 3) {
          gf = gfp[0] + (gfp[1] - gfp[0]) * (x - esp[0]) / (esp[1] - esp[0]);
        }
        if (fD[nrs + min] == 0.0 && INT > 3) {
          dr = drp[0] + (drp[1] - drp[0]) * (x - esp[0]) / (esp[1] - esp[0]);
        }
        GNR = fAmun[nrsj] * gnr * CalcPene(x, lVal) * sqrt(x) / GetRho(x, lVal);
        WidthFluctuation(GNR, gx, gg, gf, nrsj);
        double temp1 = 2. * PI * GNR * fGJ[nrsj] / dr;
        facElastic   = (GNR * fRN - 2. * sfil2) * temp1;
        facCapture   = gg * temp1 * fRG;
        facFission   = gf * temp1 * fRF;

        sumElastic += facElastic;
        sumCapture += facCapture;
        sumFission += facFission;
        nrs += nrsl;
        nrsj += 1;
      }
      sumCrsElastic += pibyk2 * (fac1 + sumElastic);
      sumCrsCapture += pibyk2 * (sumCapture);
      sumCrsFission += pibyk2 * (sumFission);
    }
    siga = sumCrsElastic;
    sigb = sumCrsCapture;
    sigc = sumCrsFission;
  } break;
  }
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::BackCrsAdler(double x, int nx)
{
  double crs1 = 0.0, crs2 = 0.0, backCrs = 0.0, pibyk2;
  double phil, sfil, sfil2, s2fil, c2fil, deno;
  phil    = CalcPhi(x, 0);
  sfil    = sin(phil);
  sfil2   = x2(sfil);
  s2fil   = sin(2.0 * phil);
  c2fil   = cos(2.0 * phil);
  pibyk2  = PI / x2(K_wnum(x));
  backCrs = sqrt(x) * pibyk2 /
            (fAt1[nx] + fAt2[nx] / x + fAt3[nx] / x2(x) + fAt4[nx] / (x * x * x) + fBt1[nx] * x + fBt2[nx] * x2(x));
  crs1 = 4.0 * pibyk2 * sfil2;
  switch (nx) {
  case 0: {
    for (int nrs = 0; nrs < totalAdler; nrs++) {
      deno = (fDet1[nrs] - x) * (fDet1[nrs] - x) + fDwt1[nrs] * fDwt1[nrs];
      crs2 = sqrt(x) * (fDwt1[nrs] * (fGrt1[nrs] * c2fil + fGit1[nrs] * s2fil) +
                        (fDet1[nrs] - x) * (fGit1[nrs] * c2fil - fGrt1[nrs] * s2fil)) /
             deno;
    }
    crsAdler[nx] = crs1 + crs2 + backCrs;
  } break;
  case 1: {
    for (int nrs = 0; nrs < totalAdler; nrs++) {
      deno = (fDef1[nrs] - x) * (fDef1[nrs] - x) + fDwf1[nrs] * fDwf1[nrs];
      crs2 = sqrt(x) * (fDwf1[nrs] * fGrf1[nrs] + (fDef1[nrs] - x) * fGif1[nrs]) / deno;
    }
    crsAdler[nx] = crs2 + backCrs;
  } break;
  case 2: {
    for (int nrs = 0; nrs < totalAdler; nrs++) {
      deno = (fDec1[nrs] - x) * (fDec1[nrs] - x) + fDwc1[nrs] * fDwc1[nrs];
      crs2 = sqrt(x) * (fDwc1[nrs] * fGrc1[nrs] + (fDec1[nrs] - x) * fGic1[nrs]) / deno;
    }
    crsAdler[nx] = crs2 + backCrs;
  } break;
  }
  return 0;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::RecursionLinear(double x1, double x2, double sig1, double sig2, double sig3, double sig4,
                                       double sig5, double sig6)
{
  double siga, sigb, sigc, elRatio = 1E-5, capRatio = 1E-5, fisRatio = 1E-5;
  if (fMloop > 10000) return 0;
//   double x10 = x1 + 0.001 * (x2 - x1);
//   double siga0, sigb0, sigc0;
//   GetSigma(LRF, x10, siga0, sigb0, sigc0);
//   double slope10 = (siga0 - sig1) / (x10 - x1);
//   double slope11 = (sigb0 - sig3) / (x10 - x1);
//   double slope12 = (sigc0 - sig5) / (x10 - x1);
// 
//   double x20 = x2 - 0.001 * (x2 - x1);
//   double siga1, sigb1, sigc1;
//   GetSigma(LRF, x20, siga1, sigb1, sigc1);
//   double slope20 = (sig2 - siga1) / (x2 - x20);
//   double slope21 = (sig4 - sigb1) / (x2 - x20);
//   double slope22 = (sig6 - sigc1) / (x2 - x20);

  double mid = 0.5 * (x1 + x2);
  if ((sig1 == 0.0 && sig2 == 0.0) || x1 == x2 || (x1 < 1E-5 || x2 < 1E-5)) return 0;
  GetSigma(LRF, mid, siga, sigb, sigc);
  double sigmid1         = sig1 + (sig2 - sig1) * (mid - x1) / (x2 - x1);
  double sigmid2         = sig3 + (sig4 - sig3) * (mid - x1) / (x2 - x1);
  double sigmid3         = sig5 + (sig6 - sig5) * (mid - x1) / (x2 - x1);
  if (siga > 0) elRatio  = std::fabs((siga - sigmid1) / siga);
  if (sigb > 0) capRatio = std::fabs((sigb - sigmid2) / sigb);
  if (sigc > 0) fisRatio = std::fabs((sigc - sigmid3) / sigc);
  fMloop++;
//   int slope = -1;
// 
//   if (slope10 / slope20 < 0.0) {
//     slope = 1;
//   }
//   if (slope11 / slope21 < 0.0) {
//     slope = 1;
//   }
//   if (slope12 / slope22 < 0.0) {
//     slope = 1;
//   }
    if (elRatio <= fSigDiff && capRatio <= fSigDiff && fisRatio <= fSigDiff ) {
  // if (elRatio <= fSigDiff && capRatio <= fSigDiff && fisRatio <= fSigDiff && slope == -1) {
  // if (elRatio <= fSigDiff && capRatio <= fSigDiff && fisRatio <= fSigDiff && fabs (sig2/siga -1) < 0.01
  //     && fabs (sig4/sigb -1) < 0.01 && slope == -1) {
    return 0;
  } else {
    fELinElastic.push_back(mid);
    fXLinElastic.push_back(siga);
    fELinCapture.push_back(mid);
    fXLinCapture.push_back(sigb);
    fELinFission.push_back(mid);
    fXLinFission.push_back(sigc);
  }
  RecursionLinear(x1, mid, sig1, siga, sig3, sigb, sig5, sigc);
  RecursionLinear(mid, x2, siga, sig2, sigb, sig4, sigc, sig6);
  return 0;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::RecursionLinear(double x1, double x2, double sig1, double sig2)
{
  double siga, sigb, sigc, elRatio;
  double mid = 0.5 * (x1 + x2);
  if ((sig1 == 0.0 && sig2 == 0.0) || x1 == x2 || (x1 < 1E-5 || x2 < 1E-5)) return 0;
  GetSigma(LRF, mid, siga, sigb, sigc);
  double sigmid1 = sig1 + (sig2 - sig1) * (mid - x1) / (x2 - x1);
  elRatio        = std::fabs((siga - sigmid1) / siga);
  if (elRatio <= fSigDiff) {
    return 0;
  } else {
    fELinElastic.push_back(mid);
    fXLinElastic.push_back(siga);
    fELinCapture.push_back(mid);
    fXLinCapture.push_back(sigb);
    fELinFission.push_back(mid);
    fXLinFission.push_back(sigc);
  }
  RecursionLinear(x1, mid, sig1, siga);
  RecursionLinear(mid, x2, siga, sig2);
  return 0;
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::AdditionalSigma(int LRF, double x)
{
  double siga, sigb, sigc;
  GetSigma(LRF, x, siga, sigb, sigc);
  fELinElastic.push_back(x);
  fELinCapture.push_back(x);
  fELinFission.push_back(x);
  fXLinElastic.push_back(siga);
  fXLinCapture.push_back(sigb);
  fXLinFission.push_back(sigc);
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::RecoPlusBroad(int flagNer)
{
  int nrs = 0, nrsl;
  if (LRU == 1) {
    for (int j = 0; j < NLS; j++) {
      nrsl = fNRS[l[j]];
      for (int i = nrs; i < nrs + fNRS[l[j]]; i++) {
        if (fEr[i] < fElo || fEr[i] > fEhi1) {
          continue;
        }
        AdditionalSigma(LRF, fEr[i]);
      }
      nrs += nrsl;
    }
    AdditionalSigma(LRF, fEhi1);
    if (flagNer == 0) {
      double eneLow = 1.0E-5;
      do {
        if (eneLow >= fElo) AdditionalSigma(LRF, eneLow);
        eneLow *= 1.03;
      } while (eneLow < 1.0E3 && eneLow < fEhi1);
      eneLow = 1.0E3;
      do {
        for (int incrE = 0; incrE < 10; incrE++) {
          if (eneLow >= 1.0E3 && eneLow + (incrE * eneLow) < fEhi1) {
            AdditionalSigma(LRF, eneLow + (incrE * eneLow));
          }
        }
        eneLow *= 10;
      } while (eneLow < 1.0E3 && eneLow < fEhi1);
    }
    TNudyCore::Instance()->Sort(fELinElastic, fXLinElastic);
    TNudyCore::Instance()->Sort(fELinCapture, fXLinCapture);
    TNudyCore::Instance()->Sort(fELinFission, fXLinFission);
    TNudyCore::Instance()->ThinningDuplicate(fELinElastic, fXLinElastic);
    TNudyCore::Instance()->ThinningDuplicate(fELinCapture, fXLinCapture);
    TNudyCore::Instance()->ThinningDuplicate(fELinFission, fXLinFission);
    int nvectorend = fELinElastic.size();
    for (int ju = intLinLru1; ju < nvectorend - 1; ju++) {
      fMloop = 0;
    //  std::cout << fELinElastic[ju] <<"  "<< fELinElastic[ju+1] <<"  "<< fXLinElastic[ju]<<"  "<<
    //  fXLinElastic[ju+1] <<"  "<< fXLinCapture[ju]<<"  "<< fXLinCapture[ju+1] <<"  "<< 
    //  fXLinFission[ju]<<"  "<< fXLinFission[ju+1]<< std::endl ;
      RecursionLinear(fELinElastic[ju], fELinElastic[ju + 1], fXLinElastic[ju], fXLinElastic[ju + 1], 
                      fXLinCapture[ju], fXLinCapture[ju + 1], fXLinFission[ju], fXLinFission[ju + 1]);
    //  std::cout<<" "<< std::endl;
    //  std::cout<<"No. of points created "<<fMloop << std::endl;
    //  std::cout<<" "<< std::endl;
      
    }
    intLinLru1 = fELinElastic.size();
  } else {
    intLinLru1 = fELinElastic.size();
    for (int l = 0; l < fJSM[0]; l++) {
      double ei = fEs[l];
      if (ei == fEs[0]) {
        ei += 0.01;
        AdditionalSigma(LRF, ei);
      } else {
        double gape = (fEs[l] - fEs[l - 1]) / 100;
        ei          = fEs[l - 1] + gape;
        do {
          if (ei >= fEs[l]) {
            ei = fEs[l];
          }
          AdditionalSigma(LRF, ei);
          ei += gape;
        } while (ei <= fEs[l]);
      }
    }
    AdditionalSigma(LRF, fEhi2);
  }
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::Linearize(int flagNer)
{
  double a = 0.0;
  switch (LRF) {
  case 1:
  case 2:
  case 3: {
    RecoPlusBroad(flagNer);
  } break;
  case 4: {
    for (int i = 1; i < (int)fEhi; i++) {
      a                      = (double)i;
      if (fEhi > 100000.0) a = (double)i * 10.;
      if (fEhi < 10000.0) a  = (double)i / 10.;
      if (fEhi < 5000.0) a   = (double)i / 100.;
      a += fElo;
      crsAdler[0] = BackCrsAdler(a, 0);
      crsAdler[1] = BackCrsAdler(a, 1);
      crsAdler[2] = BackCrsAdler(a, 2);
      crsAdler[3] = crsAdler[0] - crsAdler[1] - crsAdler[2];
    }
  } break;
  }
}
//******************** Get Resonant PArameter and cross-section data from ENDF file **********
void TNudyEndfSigma::GetData(const char *rENDF, double isigDiff)
{
  fSigDiff = isigDiff;
  for (int i = 0; i < 10; i++) {
    fMTChargeFlag[i] = 0;
  }
  //  std::cout << " file " << rENDF << " to be opened." << std::endl;
  TFile *rEND = TFile::Open(rENDF, "UPDATE");
  if (!rEND || rEND->IsZombie()) {
    printf("Error: TFile :: Cannot open file %s\n", rENDF);
  }
  TKey *rkey              = (TKey *)rEND->GetListOfKeys()->First();
  TNudyEndfTape *rENDFVol = (TNudyEndfTape *)rkey->ReadObj();
  TNudyEndfMat *tMat      = 0;
  TList *mats             = (TList *)rENDFVol->GetMats();
  int nmats               = mats->GetEntries();
  for (int iMat = 0; iMat < nmats; iMat++) {
    tMat = (TNudyEndfMat *)mats->At(iMat);
    TIter iter(tMat->GetFiles());
    TNudyEndfFile *file, *file2;
    std::vector<int>().swap(MtNumAng4Photon);
    while ((file = (TNudyEndfFile *)iter.Next())) {
      /*
            // commented due to break in photon file h. kumawat 05/04/2018 photon is to be completed
      // //       if (file->GetMF() > 13 && (file->GetMF() <= 15 || file->GetMF() <= 33)){
      // // 	       //std::cout<< "add file 2 called "<< file->GetMF() << std::endl;
      // // 	       tMat->Add(file2);
      // //       }
      */
      switch (file->GetMF()) {
      case 1:
        if (fFlagRead != 1) {
          for (int nxc = 0; nxc < tMat->GetNXC(); nxc++) {
            int mfn = tMat->GetMFn(nxc);
            int mtn = tMat->GetMTn(nxc);
            if (mfn == 3) {
              if (mtn == 103) {
                fMTChargeFlag[0] = -1;
              } else if (mtn == 104) {
                fMTChargeFlag[1] = -1;
              } else if (mtn == 105) {
                fMTChargeFlag[2] = -1;
              } else if (mtn == 106) {
                fMTChargeFlag[3] = -1;
              } else if (mtn == 107) {
                fMTChargeFlag[4] = -1;
              } else if (mtn >= 600 && mtn < 650) {
                fMTChargeFlag[5] = -1;
              } else if (mtn >= 650 && mtn < 700) {
                fMTChargeFlag[6] = -1;
              } else if (mtn >= 700 && mtn < 750) {
                fMTChargeFlag[7] = -1;
              } else if (mtn >= 750 && mtn < 800) {
                fMTChargeFlag[8] = -1;
              } else if (mtn >= 800 && mtn < 850) {
                fMTChargeFlag[9] = -1;
              }
            } else if (mfn == 12 || mfn == 13) {
              MtNumSig4Photon.push_back(mtn);
              // std::cout <<" fMF = \t"<< mfn <<" fMT = \t"<< mtn << std::endl;
            } else if (mfn == 14) {
              MtNumAng4Photon.push_back(mtn);
              // std::cout <<" fMF = \t"<< mfn <<" fMT = \t"<< mtn << std::endl;
            } else if (mfn == 15) {
              MtNumEng4Photon.push_back(mtn);
              // std::cout <<" fMF = \t"<< mfn <<" fMT = \t"<< mtn << std::endl;
            }
          }
          if (MtNumAng4Photon.size() == 0) {
            file2 = new TNudyEndfFile(tMat->GetMAT(), 14);
          } else {
            file2 = file;
          }
          ReadFile1(file);
          fFlagRead = 1;
          std::cout << "file 1 OK: Should be printed only one time " << std::endl;
          // ReWriteFile1(file);
        }
        break;
      case 2:
        if (fPrepro == 0) {
          // std::cout<<"before file 2 reading "<< std::endl;
          ReadFile2(file);
          // std::cout<<"file 2 reading OK "<< std::endl;
          if (LRU != 0) {
	  //  std::cout << fELinElastic.size() <<"   "<< fELinCapture.size() <<"  "<< fELinFission.size() << std::endl;
            TNudyCore::Instance()->Sort(fELinElastic, fXLinElastic);
            TNudyCore::Instance()->Sort(fELinCapture, fXLinCapture);
            TNudyCore::Instance()->Sort(fELinFission, fXLinFission);
	  //  std::cout << fELinElastic.size() <<"   "<< fELinCapture.size() <<"  "<< fELinFission.size() << std::endl;
//             TNudyCore::Instance()->ThinningDuplicate(fELinElastic, fXLinElastic);
//             TNudyCore::Instance()->ThinningDuplicate(fELinCapture, fXLinCapture);
//             TNudyCore::Instance()->ThinningDuplicate(fELinFission, fXLinFission);
	  //  std::cout << fELinElastic.size() <<"   "<< fELinCapture.size() <<"  "<< fELinFission.size() << std::endl;
//             Thinning(fELinElastic, fXLinElastic);
//             Thinning(fELinCapture, fXLinCapture);
//             Thinning(fELinFission, fXLinFission);
	  //  std::cout << fELinElastic.size() <<"   "<< fELinCapture.size() <<"  "<< fELinFission.size() << std::endl;
          }
        }
        /*
        //         double siga, sigb, sigc;
        //	for (unsigned int i = 0; i < fELinCapture.size() ; i++) {
        // 	    GetSigma(3, fELinCapture[i], siga, sigb, sigc);
        // 	    std::cout <<std::setprecision(12)<< fELinCapture[i] <<"   "<< fXLinCapture[i] <<"   "<< sigb <<"  "<<
        fXLinCapture[i] - sigb <<  std::endl;
        // 	}
        //               std::cout << fELinElastic.size() << std::endl;
        //            for (unsigned long j = 0 ; j < fELinElastic.size() ; j++)
        // 	     std::cout << std::setprecision(12) << fELinElastic [ j ] << "  " << fXLinElastic [ j ] << std::endl;
        //               std::cout << fELinFission.size() << std::endl;
        //            for (unsigned long j = 0 ; j < fELinFission.size() ; j++)
        //              std::cout << std::setprecision(12) << fELinFission [ j ] << "  " << fXLinFission [ j ] <<
        std::endl;
        //              std::cout << fELinCapture.size() << std::endl;
        //             for (unsigned long j = 0 ; j < fELinCapture.size() ; j++)
        //               std::cout << std::setprecision(12) << fELinCapture [ j ] << "  " << fXLinCapture [ j ] <<
        std::endl;
                   std::cout<<"file 2 OK "<< std::endl;
        */
        break;
      case 3: {
        //           for (unsigned long j = 0 ; j < fELinCapture.size() ; j++)
        //             std::cout << std::setprecision(12) << fELinCapture [ j ] << "  \t" << fXLinCapture [ j ] <<
        //             std::endl;
        ReadFile3(file);
        fSigma.clear();
        /*
                     std::cout << fELinElastic.size() << std::endl;
               std::cout << fELinCapture.size() << std::endl;
               std::cout << fELinFission.size() << std::endl;
                    // std::cout << "file 3 OK \t" << fSigmaOfMts.size() << std::endl;
        */
        TNudyCore::Instance()->ThinningDuplicate(fELinElastic, fXLinElastic);
        TNudyCore::Instance()->ThinningDuplicate(fELinCapture, fXLinCapture);
        TNudyCore::Instance()->ThinningDuplicate(fELinFission, fXLinFission);
        /*
        std::cout << "before elstic Doppler begins " << fELinElastic.size() <<"  \t"<< fXLinElastic.size() <<
        std::endl;
        //           for (unsigned long j = 0 ; j < fELinElastic.size() ; j++)
        //             std::cout << std::setprecision(12) << fELinElastic [ j ] << "  \t" << fXLinElastic [ j ] <<
        //std::endl;
        BroadSigma(fELinElastic, fXLinElastic, fXBroadElastic);
        TNudyCore::Instance()->Sort(fELinElastic, fXBroadElastic);
        //           std::cout << "after elstic Doppler " << fELinElastic.size() <<"  \t"<< fXBroadElastic.size() <<
        //std::endl;
        //           for (unsigned long j = 0 ; j < fELinFission.size() ; j++)
        //             std::cout << std::setprecision(12) << fELinFission [ j ] << "  " << fXLinFission [ j ] <<
        //std::endl;
        if(fPrepro==0)Thinning(fELinElastic, fXBroadElastic);
        // std::cout << fELinElastic.size() << std::endl;
        std::cout << "before capture Doppler begins " << fELinCapture.size() <<"  "<< fXLinCapture.size() <<
        std::endl;
        //            for (unsigned long j = 0 ; j < fELinCapture.size() ; j++)
        //              std::cout << std::setprecision(12) << fELinCapture [ j ] << "  " << fXLinCapture [ j ] <<
        //std::endl;
        BroadSigma(fELinCapture, fXLinCapture, fXBroadCapture);
        TNudyCore::Instance()->Sort(fELinCapture, fXBroadCapture);
        //           std::cout << "after capture Doppler begins " << fELinCapture.size() <<"  "<< fXBroadCapture.size()
        //<< std::endl;
        if(fPrepro==0)Thinning(fELinCapture, fXBroadCapture);
        // std::cout << fELinCapture.size() << std::endl;
        std::cout << "before fission Doppler begins "<< fELinFission.size() <<"  "<< fXLinFission.size() <<
        std::endl;
        BroadSigma(fELinFission, fXLinFission, fXBroadFission);
        TNudyCore::Instance()->Sort(fELinFission, fXBroadFission);
        if(fPrepro==0)Thinning(fELinFission, fXBroadFission);
        //            std::cout << "after Fission Doppler begins " << fELinFission.size() <<"  "<< fXBroadFission.size()
        //<< std::endl;
        //             for (unsigned long j = 0 ; j < fELinFission.size() ; j++)
        //               std::cout << std::setprecision(12) << fELinFission [ j ] << "  " << fXBroadFission [ j ] <<
        //std::endl;
        // std::cout << fELinFission.size() << std::endl;
        // std::cout<<"Doppler done "<<outstring << std::endl;
        fDopplerBroad = 0;

        out << fELinElastic.size() << std::endl;
        int ELinElasticSize = fELinElastic.size();
        for (int j = 0; j != ELinElasticSize; ++j)
           out << std::setprecision(12)<< fELinElastic[j] << "  " << fXBroadElastic[j] << std::endl;
        fELinElastic.clear();
        fXLinElastic.clear();
        fXBroadElastic.clear();
        out << fELinCapture.size() << std::endl;
        int ELinCaptureSize = fELinCapture.size();
        for (int j = 0; j != ELinCaptureSize; ++j)
           out << std::setprecision(12) << fELinCapture[j] << "  " << fXBroadCapture[j] << std::endl;
        fELinCapture.clear();
        fXLinCapture.clear();
        fXBroadCapture.clear();
        out << fELinFission.size() << std::endl;
        int ELinFissionSize = fELinFission.size();
        for (int j = 0; j != ELinFissionSize; ++j)
           out << std::setprecision(12) << fELinFission[j] << "  " << fXBroadFission[j] << std::endl;
        fELinFission.clear();
        fXLinFission.clear();
        fXBroadFission.clear();

        //std::cout << "dopplerall started \t" << std::endl;
        */
        DopplerAll();
        // std::cout << "dopplerall OK \t" << std::endl;
        FixupTotal(file, energyUni, sigmaUniTotal);
        // std::cout << "fixup total OK \t" <<energyUni.size() << std::endl;

        ReWriteFile3(file);
        // std::cout << "after rewrite file3 OK \t" << std::endl;
        MtNumbers.clear();
      } break;

      // //       case 4:
      // // //        ReadFile4(file);
      // //         std::cout << "file 4 OK " << std::endl;
      // //         //ReWriteFile4(file);
      // //         MtNumbers.clear();
      // //         break;
      // // 	case 5:
      // // //	recoEnergy = new TNudyEndfEnergy(file);
      // // 	std::cout<<"file 5 OK "<<std::endl;
      // // 	break;
      case 6: {
        TIter iter1(tMat->GetFiles());
        TNudyEndfFile *file1;
        while ((file1 = (TNudyEndfFile *)iter1.Next())) {
          if (file1->GetMF() == 6) {
            // std::cout << "file 6 found " << std::endl ;
            TIter secIter1(file1->GetSections());
            TNudyEndfSec *sec1;
            while ((sec1 = (TNudyEndfSec *)secIter1.Next())) {
              int NK = sec1->GetN1();
              TIter recIter1(sec1->GetRecords());
              for (int k = 0; k < NK; k++) {
                // std::cout<<k<<" NK "<<NK<< std::endl;
                TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter1.Next();
                div_t divr;
                AWRI    = sec1->GetC2();
                int fMT = sec1->GetMT();
                //      int LCT = sec->GetL2();
                divr = div(sec1->GetC1(), 1000);
                //		double ZA  = divr.quot;
                //		double AA  = divr.rem;
                //		double NA1 = AA - ZA;
                int ZAP = tab1->GetC1();
                //		double AWP = tab1->GetC2();
                // std::cout<<" ZAP = \t"<< ZAP <<" fMT "<< fMT <<  std::endl;
                // int LIP =tab1->GetL1();
                int LAW = tab1->GetL2();
                if (LAW == 3 || LAW == 4 || LAW == 0) continue;
                divr          = div(ZAP, 1000);
                int particleA = divr.rem;
                int particleZ = divr.quot;
                // std::cout<<"LAW = "<< LAW <<" fMT "<<fMT <<" ZAP = "<< ZAP <<" AWP "<<AWP <<" parZ "<< particleZ
                // <<" parA "<< particleA << std::endl;
                if (LAW == 2 && particleZ == 0 && particleA == 0) {
                  int LANG, NL;
                  TNudyEndfTab2 *tab2   = (TNudyEndfTab2 *)recIter1.Next();
                  int ne2               = tab2->GetN2();
                  TNudyEndfList *header = (TNudyEndfList *)recIter1.Next();
                  LANG                  = header->GetL1();
                  TNudyEndfSec *sec3;
                  if (LANG == 0) {
                    sec3 = new TNudyEndfSec(sec1->GetMAT(), 14, fMT, sec1->GetC1(), AWRI, 0, 1, 1, 0);
                  } else {
                    sec3 = new TNudyEndfSec(sec1->GetMAT(), 14, fMT, sec1->GetC1(), AWRI, 0, 2, 1, 0);
                  }
                  TNudyEndfTab2 *sec3Tab2 = new TNudyEndfTab2();
                  sec3Tab2->SetCont(tab2->GetC1(), tab2->GetC2(), 0, 0, tab2->GetN1(), ne2);
                  sec3Tab2->SetNBT(ne2, 0);
                  sec3Tab2->SetINT(tab2->GetINT(0), 0);
                  sec3->Add(sec3Tab2);
                  // std::cout<<"tab2->GetINT(0) "<< sec3Tab2->GetINT(0) <<" fNR "<< sec3Tab2->GetN1() <<" fNE "<<
                  // sec3Tab2->GetN2() << std::endl;

                  TNudyEndfList *sec3List = new TNudyEndfList();
                  sec3List->SetCont(header->GetC1(), header->GetC2(), 0, 0, header->GetN2(), 0);
                  // std::cout<<"C1 "<< header->GetC1() <<" C2 "<< header->GetC2() <<" N2 "<< header->GetN2() <<
                  // std::endl;
                  for (int j = 0, eN2 = header->GetN2(); j != eN2; ++j) {
                    sec3List->SetLIST(header->GetLIST(j), j);
                  }
                  sec3->Add(sec3List);
                  // std::cout<<"sec3List C1 "<< sec3List->GetC1()<<" sec3List C2 "<< sec3List->GetC2() <<" nl "<<
                  // sec3List->GetN1() << std::endl;
                  // for (int j = 0; j < sec3List->GetN1(); ++j) {
                  // std::cout <<"list "<< sec3List->GetLIST(j) << std::endl;
                  //}
                  for (int lis = 1; lis < ne2; lis++) {
                    TNudyEndfList *header   = (TNudyEndfList *)recIter1.Next();
                    TNudyEndfList *sec3List = new TNudyEndfList();
                    sec3List->SetCont(header->GetC1(), header->GetC2(), 0, 0, header->GetN2(), 0);
                    NL = header->GetN2();
                    if (LANG == 0) {
                      for (int j = 0; j < NL; j++) {
                        sec3List->SetLIST(header->GetLIST(j), j);
                      }
                      //  std::cout<<"sec3List C1 "<< sec3List->GetC1()<<" sec3List C2 "<< sec3List->GetC2() <<" nl "<<
                      //  sec3List->GetN1() << std::endl;
                      // for (int j = 0; j < sec3List->GetN1(); ++j) {
                      // std::cout <<"list "<< sec3List->GetLIST(j) << std::endl;
                      // }
                    } else if (LANG > 0) {
                      for (int i = 0; i < NL; i++) {
                        sec3List->SetLIST(header->GetLIST(2 * i + 0), 2 * i + 0);
                        sec3List->SetLIST(header->GetLIST(2 * i + 1), 2 * i + 1);
                      }
                    }
                    sec3->Add(sec3List);
                  }
                  file2->Add(sec3);
                } else if (LAW == 1 || LAW == 5) {
                  TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter1.Next();
                  for (int lis = 0, eN2 = tab2->GetN2(); lis != eN2; ++lis) {
                    //		      TNudyEndfList *header = (TNudyEndfList *)recIter1.Next();
                  }
                } else if (LAW == 6) {
                  //		    TNudyEndfCont *header = (TNudyEndfCont *)recIter1.Next();
                } else if (LAW == 7) {
                  TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter1.Next();
                  for (int cr1 = 0, eN2 = tab2->GetN2(); cr1 != eN2; ++cr1) {
                    TNudyEndfTab2 *tab3 = (TNudyEndfTab2 *)recIter1.Next();
                    for (int cr2 = 0, e3N2 = tab3->GetN2(); cr2 != e3N2; ++cr2) {
                      //			TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter1.Next();
                    }
                  }
                }
              }
            }
          }
        }
        //	recoEnergyAng = new TNudyEndfEnergyAng(file,fQValue);
        //	std::cout<<"file 6 OK "<<std::endl;
      } break;
        /*
        // // 	case 8:
        // // //	recoFissY = new TNudyEndfFissionYield(file);
        // // 	std::cout<<"file 8 OK "<<std::endl;
        // // 	break;
        //	case 12:
        //	recoPhYield = new TNudyEndfPhYield(file);
        // 	std::cout<<"file 12 OK "<<std::endl;
        //	break;
        //	case 13:
        //	recoPhProd = new TNudyEndfPhProd(file);
        // 	std::cout<<"file 13 OK "<<std::endl;
        //	break;
        //	case 14:
        //	recoPhAng = new TNudyEndfPhAng(file);
          // std::cout<<"file 14 OK "<<std::endl;
        //	break;
        //	case 15:
        //	recoPhEnergy = new TNudyEndfPhEnergy(file);
        // 	std::cout<<"file 15 OK "<<std::endl;
        //	break;
        */
      }
    }
  }
  // for (unsigned int i = 0; i < MtNumSig4Photon.size () ; i++){
  // std::cout <<"Photon Sig fMT = \t"<< MtNumSig4Photon[i] << std::endl;
  //}
  // for (unsigned int i = 0; i < MtNumAng4Photon.size () ; i++){
  // std::cout <<"Photon Ang fMT = \t"<< MtNumAng4Photon[i] << std::endl;
  //}
  // for (unsigned int i = 0; i < MtNumEng4Photon.size () ; i++){
  // std::cout <<"Photon Ene fMT = \t"<< MtNumEng4Photon[i] << std::endl;
  //}
  rENDFVol->Write();
  rEND->Close();
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::BroadSigma(std::vector<double> &x1, std::vector<double> &x2, std::vector<double> &x3)
{
  //     std::cout<<"before sort \t"<< x1.size() << std::endl;
  TNudyCore::Instance()->Sort(x1, x2);
  if (x1.size() > 0) {
    doppler = new TNudyEndfDoppler(fSigDiff, AWRI, fDoppTemp1, fDoppTemp2, x1, x2);
    //      std::cout<<"after doppler \t"<< x1.size() << std::endl;
    fDopplerBroad = 1;
    for (int j = 0, x1Size = x1.size(); j != x1Size; ++j) {
      //      std::cout<<"size2 "<< x1[j] <<"  "<< x2[j] <<"  "<< doppler->fSigma[j]<< std::endl;
      x3.push_back(doppler->fSigma[j]);
    }
  }
  doppler->fSigma.clear();
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::FixupTotal(TNudyEndfFile *file1, std::vector<double> &x1, std::vector<double> &x2)
{
  std::sort(x1.begin(), x1.end()); // Unionization of energy grid for cross-section
  TNudyCore::Instance()->ThinningDuplicate(x1);
  TIter secIter(file1->GetSections());
  TNudyEndfSec *sec;
  for (int i = 0, SigmaOfMtsDopSize =  fSigmaOfMtsDop.size(); i != SigmaOfMtsDopSize; ++i) {
    while ((sec = (TNudyEndfSec *)secIter.Next())) {
      //std::cout << sec->GetMT() <<" fixuptotal \t"<< MtNumbers[i] << std::endl;
      int fMT = sec->GetMT();
      if (fMT == MtNumbers[i]) {
        // 	if (fMT != 1 && fMT != 3 && fMT != 4 && fMT != 27 && fMT != 19 && fMT != 20 && fMT != 21 && fMT != 38 && fMT
        // != 101 && fMT < 120) {
        int size = fSigmaOfMtsDop[i].size() / 2;
        for (int k = 0, x1Size = x1.size(); k != x1Size; ++k) {
          int min = 0;
          int max = size - 1;
          int mid = 0;
          if (x1[k] <= fSigmaOfMtsDop[i][min])
            min = 0;
          else if (x1[k] >= fSigmaOfMtsDop[i][max])
            min = max - 1;
          else {
            while (max - min > 1) {
              mid = (min + max) / 2;
              if (x1[k] < fSigmaOfMtsDop[i][mid])
                max = mid;
              else
                min = mid;
            }
          }
          if (x1[k] == fSigmaOfMtsDop[i][min] && fSigmaOfMtsDop[i][min] > 0) {
            fEneTemp.push_back(x1[k]);
            fSigTemp.push_back(fSigmaOfMtsDop[i][size + min]);
          } else {
            double sigmaAdd = 0;
            if (x1[k] > fSigmaOfMtsDop[i][min]) {
              sigmaAdd = fSigmaOfMtsDop[i][size + min] +
                         (fSigmaOfMtsDop[i][size + min + 1] - fSigmaOfMtsDop[i][size + min]) *
                             (x1[k] - fSigmaOfMtsDop[i][min]) /
                             (fSigmaOfMtsDop[i][min + 1] - fSigmaOfMtsDop[i][min]); // linear interpolation
            }
            if (sigmaAdd > 0) {
              fEneTemp.push_back(x1[k]);
              fSigTemp.push_back(sigmaAdd);
            }
          }
          if (fEneTemp.size() == 1) {
            fEnergyLocationMts.push_back(k);
          }
        }

        fSigmaUniOfMts.push_back(fSigTemp);
        fEneTemp.clear();
        fSigTemp.clear();
        // 	}
      }
      break;
    }
  }
  x2.resize(x1.size());
  fSigmaOfMtsDop.clear();
  for (int i = 0, SigmaUniOfMtsSize = fSigmaUniOfMts.size(); i != SigmaUniOfMtsSize; ++i) {
    int size = fSigmaUniOfMts[i].size();
    int fMT  = MtNumbers[i];
    if (fMT != 1 && fMT != 3 && fMT != 4 && fMT != 27 && fMT != 19 && fMT != 20 && fMT != 21 && fMT != 38 &&
        fMT != 101 && fMT < 120) {
      // std::cout<<"size "<<  fEnergyLocationMts.size  << std::endl;
      for (int j = 0; j < size; j++) {
        x2[fEnergyLocationMts[i] + j] += fSigmaUniOfMts[i][j];
      }
    }
  }
  outtotal << x1.size() << std::endl;
  for (int j = 0, x1Size = x1.size(); j != x1Size; ++j) {
    outtotal << x1[j] << "  " << x2[j] << std::endl;
  }
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::DopplerAll()
{
  for (int i = 0, SigmaOfMtsSize = fSigmaOfMts.size(); i != SigmaOfMtsSize; ++i) {
    //std::cout<<"MT number \t"<< MtNumbers[i] <<"  \t"<< fSigmaOfMts.size() << std::endl;
    // std::cout <<"  "<< std::endl;
    int size = fSigmaOfMts[i].size() / 2;
    for (int k = 0; k < size; k++) {
      if (fSigmaOfMts[i][size + k] > 1E-20) {
        fEneTemp.push_back(fSigmaOfMts[i][k]);
        fSigTemp.push_back(fSigmaOfMts[i][size + k]);
        //	std::cout<<std::setprecision(12)<<  fSigmaOfMts [ i ][ k ]<<"  \t"<< fSigmaOfMts [ i ][ size + k ]
        //<<std::endl;
      }
    }
    /* future testing
    //      if(MtNumbers[i]==1){
    //            std::cout << " " << fEneTemp.size() <<"  "<< fSigTemp.size() << std::endl;
    //             for (unsigned long j = 0 ; j < fEneTemp.size() ; j++)
    //               std::cout << std::setprecision(12) << fEneTemp [ j ] << "  " << fSigTemp [ j ] << std::endl;
    //     }
    */
    BroadSigma(fEneTemp, fSigTemp, fSigma);
    //    std::cout<<"after broad \t"<< fSigma.size() << std::endl;
    if (fPrepro == 0) {
      TNudyCore::Instance()->Sort(fEneTemp, fSigma);
      Thinning(fEneTemp, fSigma);
    }
    //   std::cout<<"after thinning \t"<< fSigma.size() << std::endl;
    fSigmaMts.insert(std::end(fSigmaMts), std::begin(fEneTemp), std::end(fEneTemp));
    fSigmaMts.insert(std::end(fSigmaMts), std::begin(fSigma), std::end(fSigma));
    fSigmaOfMtsDop.push_back(fSigmaMts);
    energyUni.insert(std::end(energyUni), std::begin(fEneTemp), std::end(fEneTemp));
    //   std::cout<<"after energyuni \t"<< energyUni.size() << std::endl;
    fEneTemp.clear();
    fSigTemp.clear();
    fSigma.clear();
    fSigmaMts.clear();
  }
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfSigma::Thinning(std::vector<double> &x1, std::vector<double> &x2)
{
  int size  = x1.size();
  int size1 = x1.size();
  if (size <= 0) return 0;
  for (int i = 0; i < size1 - 2; i++) {
    if (x1[i] < 1E3 && fDopplerBroad == 0) continue;
    double sigmid1 = x2[i] + (x2[i + 2] - x2[i]) * (x1[i + 1] - x1[i]) / (x1[i + 2] - x1[i]);
    if (std::fabs((x2[i + 1] - sigmid1) / sigmid1) <= fSigDiff) {
      x1.erase(x1.begin() + i + 1);
      x2.erase(x2.begin() + i + 1);
    }
    if (x2[i] <= 1E-20) {
      x1.erase(x1.begin() + i);
      x2.erase(x2.begin() + i);
    }
    size1 = x1.size();
  }
  size1 = x1.size();
  if (size == size1) return 0;
  Thinning(x1, x2);
  return 0;
}
//______________________________________________________________________________
void TNudyEndfSigma::ReadFile4(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
    int fMT               = sec->GetMT();
    MtNumbers.push_back(fMT);
    int LTT = sec->GetL2();
    int LI  = header->GetL1();
    int LCT = header->GetL2();
    fMtLct.push_back(LCT);
    // printf("LCT = %d LTT = %d LI = %d\n",LCT, LTT, LI);
    // Legendre polynomial coefficients
    if (LTT == 1 && LI == 0) {
      TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0, e2N2 = tab2->GetN2(); i != e2N2; ++i) {
        TNudyEndfList *tab = (TNudyEndfList *)recIter.Next();
        fEin.push_back(tab->GetC2());
        // std::cout<<"energy "<< tab->GetC2() << std::endl;
        for (int j = 0, eNPL = tab->GetNPL(); j != eNPL; ++j) {
          fLegendCoef1.push_back(tab->GetLIST(j));
        }
        fLegendCoef.push_back(fLegendCoef1);
        fLegendCoef1.clear();
      }
      for (int i = 0, EinSize = fEin.size(); i != EinSize; ++i) {
        // printf("Ein = %e\n", fEin[i]);
        int k1     = 0;
        double fme = 0.0;
        do {
          fme      = 1.0;
          double x = -1. + k1 * 0.02;
          for (int j = 0, LegendCoefSize = fLegendCoef[i].size(); j != LegendCoefSize; ++j) {
            double leg = ROOT::Math::legendre(j + 1, x);
            fme += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[i][j] * leg;
            // printf("a%d = %e leg= %e\n", j, fLegendCoef[i].At(j),leg);
          }
          if (fme > 0.0) {
            fCosFile4.push_back(x);
            fCosPdfFile4.push_back(fme);
          }
          // printf("%e %e\n", x, fme);
          k1++;
        } while (k1 < 101);
        for (int l = 0; l < 100; l++) {
          RecursionLinearLeg(i, fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
        }
        FillPdf1d();
      }
      FillPdf2d();
      fLegendCoef.clear();
      // Tabulated probability tables
    } else if (LTT == 2 && LI == 0) {
      TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0, e2N2 = tab2->GetN2(); i != e2N2; ++i) {
        TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
        fEin.push_back(tab->GetC2());
        fNR = tab->GetNR();
        fNP = tab->GetNP();
        // std::cout<<"energy "<< tab->GetC2() << std::endl;
        for (int i = 0; i < fNR; i++) {
          fNbt1.push_back(tab->GetNBT(i));
          fInt1.push_back(tab->GetINT(i));
        }
        for (int j = 0; j < fNP; j++) {
          fCosFile4.push_back(tab->GetX(j));
          fCosPdfFile4.push_back(tab->GetY(j));
        }
        for (int l = 0; l < fNP - 1; l++) {
          RecursionLinearProb(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
        }
        FillPdf1d();
        fNbt1.clear();
        fInt1.clear();
      }
      FillPdf2d();
      // Low energy given by legendre polynomial and high energy by tabulated probability tables
    } else if (LTT == 3 && LI == 0) {
      TNudyEndfTab2 *lowE = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0,elN2 = lowE->GetN2(); i != elN2; ++i) {
        TNudyEndfList *tab = (TNudyEndfList *)recIter.Next();
        fEin.push_back(tab->GetC2());
        for (int j = 0, eNPL = tab->GetNPL(); j != eNPL; ++j) {
          fLegendCoef1.push_back(tab->GetLIST(j));
        }
        fLegendCoef.push_back(fLegendCoef1);
        fLegendCoef1.clear();
      }
      for (int i = 0, EinSize = fEin.size(); i != EinSize; ++i) {
        // printf("Ein = %e\n", fEin[i]);
        int k1     = 0;
        double fme = 0.0;
        do {
          fme      = 1.0;
          double x = -1. + k1 * 0.02;
          for (int j = 0, LegendCoefSize = fLegendCoef[i].size(); j != LegendCoefSize; ++j) {
            double leg = ROOT::Math::legendre(j + 1, x);
            fme += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[i][j] * leg;
            // printf("a%d = %e leg= %e\n", j, fLegendCoef[i][j],leg);
          }
          if (fme > 0.0) {
            fCosFile4.push_back(x);
            fCosPdfFile4.push_back(fme);
          }
          // printf("%e %e\n", x, fme);
          k1++;
        } while (k1 < 101);

        for (int l = 0; l < 100; l++) {
          RecursionLinearLeg(i, fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
        }
        FillPdf1d();
      }
      fLegendCoef.clear();
      TNudyEndfTab2 *highE = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0, eEN2 = highE->GetN2(); i != eEN2; ++i) {
        TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
        fEin.push_back(tab->GetC2());
        // std::cout <<"energy "<< fEin[fEin.size() - 1] << std::endl;
        fNR = tab->GetNR();
        fNP = tab->GetNP();
        for (int i = 0; i < fNR; i++) {
          fNbt1.push_back(tab->GetNBT(i));
          fInt1.push_back(tab->GetINT(i));
        }
        for (int j = 0; j < fNP; j++) {
          fCosFile4.push_back(tab->GetX(j));
          fCosPdfFile4.push_back(tab->GetY(j));
        }
        if (fNP > 2) {
          for (int l = 0; l < fNP - 1; l++) {
            RecursionLinearProb(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
          }
        }
        FillPdf1d();
        fNbt1.clear();
        fInt1.clear();
      }
      FillPdf2d();
    } else if (LTT == 0 && LI == 1) {
      fEin.push_back(1E-14);
      fEin.push_back(1.5E8);
      for (int j = 0; j < 2; j++) {
        fCosFile4.push_back(1);
        fCosPdfFile4.push_back(0.5);
        fCosFile4.push_back(-1);
        fCosPdfFile4.push_back(0.5);
        FillPdf1d();
      }
      FillPdf2d();
    }
    // Low energy given by legendre polynomial and high energy by tabulated probability tables
  } // end while loop
}
//--------------------------------------------------------------------------------------
double TNudyEndfSigma::RecursionLinearLeg(int i, double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  for (int j = 0, LegendCoefSize = fLegendCoef[i].size(); j != LegendCoefSize; ++j) {
    double leg = ROOT::Math::legendre(j + 1, mid);
    pdf += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[i][j] * leg;
  }
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (pdf > 0 && fabs((pdf - pdfmid1) / pdf) <= 2E-4) {
    return 0;
  }
  fCosFile4.push_back(mid);
  fCosPdfFile4.push_back(pdf);
  RecursionLinearLeg(i, x1, mid, pdf1, pdf);
  RecursionLinearLeg(i, mid, x2, pdf, pdf2);
  return 0;
}
//--------------------------------------------------------------------------------------
double TNudyEndfSigma::RecursionLinearProb(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  pdf            = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fCosFile4, fCosPdfFile4, fNP, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (pdf > 0 && fabs((pdf - pdfmid1) / pdf) <= 2E-4) {
    return 0;
  }
  fCosFile4.push_back(mid);
  fCosPdfFile4.push_back(pdf);
  RecursionLinearProb(x1, mid, pdf1, pdf);
  RecursionLinearProb(mid, x2, pdf, pdf2);
  return 0;
}
//-----------------------------------------------------------------------------------------
void TNudyEndfSigma::FillPdf1d()
{
  TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
  TNudyCore::Instance()->cdfGenerateT(fCosFile4, fCosPdfFile4, fCosCdfFile4);
  for (int i = 0, CosFile4Size = fCosFile4.size(); i != CosFile4Size; ++i) {
    fCos4.push_back(fCosFile4[i]);
    fPdf.push_back(fCosPdfFile4[i]);
    fCdf.push_back(fCosCdfFile4[i]);
  }
  fCos2D.push_back(fCos4);
  fPdf2D.push_back(fPdf);
  fCdf2D.push_back(fCdf);
  fCosFile4.clear();
  fCosPdfFile4.clear();
  fCosCdfFile4.clear();
  fCos4.clear();
  fPdf.clear();
  fCdf.clear();
}
//--------------------------------------------------------------------------------------
void TNudyEndfSigma::FillPdf2d()
{
  fEin2D.push_back(fEin);
  fCos3D.push_back(fCos2D);
  fPdf3D.push_back(fPdf2D);
  fCdf3D.push_back(fCdf2D);
  fEin.clear();
  fCos2D.clear();
  fPdf2D.clear();
  fCdf2D.clear();
}
//______________________________________________________________________________
void TNudyEndfSigma::ReWriteFile1(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    int fMT = sec->GetMT();
    if (fMT == 452) { // Total fission neutron multiplicity polynomial expansion
      int LNU = sec->GetL2();

      if (LNU == 1) {
        TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
        for (int j = 0, eN1 = list->GetN1(); j != eN1; ++j) {
          fCnc.push_back(list->GetLIST(j));
        }
        double ein = 1E-5;
        do {
          double nun = 0;
          for (int i = 0, eN1 = list->GetN1(); i != eN1; ++i) {
            nun += fCnc[i] * pow(ein, i);
          }
          fEintFile1.push_back(ein);
          fNutFile1.push_back(nun);
          ein *= 2;
        } while (ein < 21E8);
        fCnc.clear();
      } else { // Total fission neutron multiplicity tabulated representation
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        fNR                 = 1;
        fNP                 = fEintFile1.size();
        tab1->SetNBT(fNP, 0);
        tab1->SetINT(2, 0);
        for (int crs = 0; crs < fNP; crs++) {
          tab1->SetX(fEintFile1[crs], crs);
          tab1->SetY(fNutFile1[crs], crs);
        }
      }
    } else if (fMT == 455) { // delayed neutron multiplicity
      int LDG = sec->GetL1();
      int LNU = sec->GetL2();
      if (LNU == 1 && LDG == 0) {

      } else if (LNU == 1 && LDG == 1) {

      } else if (LNU == 2 && LDG == 0) {
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        fNR                 = 1;
        fNP                 = fEindFile1.size();
        tab1->SetNBT(fNP, 0);
        tab1->SetINT(2, 0);
        for (int crs = 0; crs < fNP; crs++) {
          tab1->SetX(fEindFile1[crs], crs);
          tab1->SetY(fNudFile1[crs], crs);
        }
      } else if (LNU == 2 && LDG == 1) {
      }
    } else if (fMT == 456) { // prompt neutron multiplicity
      int LNU = sec->GetL2();
      if (LNU == 1) {
        //      std::cout<<"prompt nu = "<< ZA <<" LNU "<< LNU << std::endl;
        //	TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
        // double fNup = list->GetLIST(0);
      } else {
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        fNR                 = 1;
        fNP                 = fEinFile1.size();
        tab1->SetNBT(fNP, 0);
        tab1->SetINT(2, 0);
        for (int crs = 0; crs < fNP; crs++) {
          tab1->SetX(fEinFile1[crs], crs);
          tab1->SetY(fNuFile1[crs], crs);
        }
      }
    } else if (fMT == 458) {
      TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
      int NPLY            = list->GetL2();
      if (NPLY == 0) {
        double EFR = list->GetLIST(0);
        double ENP = list->GetLIST(2);
        double END = list->GetLIST(4);
        double EGP = list->GetLIST(6);
        double EGD = list->GetLIST(8);
        double EB  = list->GetLIST(10);
        double ENU = list->GetLIST(12);
        //	double ER  = list->GetLIST(14);
        // double ET  = list->GetLIST(16);
        double ein = 1E-5;
        do {
          double EFis = 0;
          EFis        = EFR + ENP + END + EGP + EGD + EB + ENU;
          EFis -= 0.100 * ein; // nuetrino energy dependence
          EFis -= 0.075 * ein; // delayed gamma energy dependence
          EFis -= 0.075 * ein; // delayed beta energy dependence
          //	      for(unsigned long i = 0; i < fEintFile1.size(); i++)
          //		std::cout<< fEintFile1[i] <<"  "<<fNutFile1[i] << std::endl;
          int n0  = 0;
          int max = fEintFile1.size() - 1;
          int mid = 0;
          if (ein <= fEintFile1[n0]) {
            n0 = 0;
          } else if (ein >= fEintFile1[max]) {
            n0 = max - 1;
          } else {
            while (max - n0 > 1) {
              mid = (n0 + max) / 2;
              if (ein < fEintFile1[mid])
                max = mid;
              else
                n0 = mid;
              //      std::cout<<"n0 "<< n0 <<" max "<< max << std::endl;
            }
          }
          double nue = TNudyCore::Instance()->LinearInterpolation(fEintFile1[n0], fNutFile1[n0], fEintFile1[n0 + 1],
                                                                  fNutFile1[n0 + 1], ein);
          double nu0 = fNutFile1[0];
          EFis -= -1.307 * ein + 8.07 * 1E6 * (nue - nu0); // prompt neutron energy dependence
          fEinfFile1.push_back(ein);
          fHeatFile1.push_back(EFis);
          ein *= 2;
        } while (ein < 21E8);
      } else {
        int nply2 = 9 * (NPLY + 1);
        double c0[nply2], c1[nply2];
        for (int i = 0; i < nply2; i++) {
          c0[i] = list->GetLIST(i);
        }
        for (int i = nply2; i < 2 * nply2; i++) {
          c1[i - nply2] = list->GetLIST(i);
        }
        double ein = 1E-5;
        do {
          double EFis = 0;
          for (int i = 0; i < 9 * NPLY / 2 - 2; i++) {
            EFis += c0[i * 2] + c1[i * 2] * ein;
          }
          fEinfFile1.push_back(ein);
          fHeatFile1.push_back(EFis);
          ein *= 2;
        } while (ein < 21E8);
      }
      fEinFissHeat.push_back(fEinfFile1);
      fFissHeat.push_back(fHeatFile1);
      fEinfFile1.clear();
      fHeatFile1.clear();
    } else if (fMT == 460) {
      int LO = sec->GetL1();
      int NG = sec->GetN1();
      if (LO == 1) {
        for (int ng = 0; ng < NG; ng++) {
          TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
          fNR                 = 1;
          fNP                 = fEinPhFile1.size();
          tab1->SetNBT(fNP, 0);
          tab1->SetINT(2, 0);
          for (int crs = 0; crs < fNP; crs++) {
            tab1->SetX(fEinPhFile1[crs], crs);
            tab1->SetY(fPhFile1[crs], crs);
          }
        }
      } else {
        TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
        int NNF             = list->GetN1();
        // double lambda[NNF];
        for (int i = 0; i < NNF; i++) {
          // lambda[i] = list->GetLIST(i);
          //	     std::cout <<"lambda  "<< lambda[i] << std::endl;
        }
      }
    }
  }
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::ReWriteFile3(TNudyEndfFile *file)
{
  TIter secIter0(file->GetSections());
  TNudyEndfSec *sec0;
  int mtcount = 0;
  while ((sec0 = (TNudyEndfSec *)secIter0.Next())) {
    int fMT = sec0->GetMT();
    if (fMT != 1 && mtcount == 0) {
      AddSecFile3(file, 0, 0, 1, energyUni, sigmaUniTotal);
      break;
    }
  }
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    //     std::cout <<"fMT  "<<sec->GetMT() <<"  \t"<< MtNumbers[mtcount] << std::endl;
    int fMT = sec->GetMT();
    //     if (fMT != 3 && fMT != 4 && fMT != 27 && fMT != 19 && fMT != 20 && fMT != 21 && fMT != 38 && fMT != 101 &&
    //     fMT < 250) {
    TIter recIter(sec->GetRecords());
    TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
    int fNR               = 1;
    int energyUniSize     = energyUni.size();
    if (fMT == 1)
      header->SetCont(header->GetC1(), header->GetC2(), header->GetL1(), header->GetL2(), fNR, energyUniSize);
    TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
    if (fMT == 1) {
      for (int cr = 0; cr < fNR; cr++) {
        tab1->SetNBT(energyUniSize, 0);
        tab1->SetINT(2, 0);
      }
      for (int crs = 0; crs != energyUniSize; ++crs) {
        tab1->SetX(energyUni[crs], crs);
        tab1->SetY(sigmaUniTotal[crs], crs);
        //  	  std::cout<<energyUni[crs] <<"  \t"<< sigmaUniTotal[crs] << std::endl;
      }
      mtcount++;
      continue;
    }
    if (fMT == MtNumbers[mtcount]) {
      int SigmaUniOfMtsSize = fSigmaUniOfMts[mtcount].size();
      header->SetCont(header->GetC1(), header->GetC2(), header->GetL1(), header->GetL2(), fNR,
                      SigmaUniOfMtsSize);
      tab1->SetNBT(SigmaUniOfMtsSize, 0);
      tab1->SetINT(2, 0);
      for (int crs = 0; crs != SigmaUniOfMtsSize; ++crs) {
        tab1->SetX(energyUni[fEnergyLocationMts[mtcount] + crs], crs);
        tab1->SetY(fSigmaUniOfMts[mtcount][crs], crs);
        //   	  std::cout << energyUni[fEnergyLocationMts[mtcount] + crs] <<"  \t"<< fSigmaUniOfMts[mtcount][crs] <<
        //   std::endl;
      }
    }
    mtcount++;
    //     }
  }
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfSigma::ReWriteFile4(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    // if (sec->GetMT() == (int)fReaction) {
    TIter recIter(sec->GetRecords());
  }
}
//-------------------------------------------------------------------------------------------------------
TNudyEndfSigma::~TNudyEndfSigma()
{
  fEintFile1.shrink_to_fit();
  fNutFile1.shrink_to_fit();
  fEinFile1.shrink_to_fit();
  fNuFile1.shrink_to_fit();
  fEindFile1.shrink_to_fit();
  fNudFile1.shrink_to_fit();
  fEinPhFile1.shrink_to_fit();
  fPhFile1.shrink_to_fit();
  fEinfFile1.shrink_to_fit();
  fHeatFile1.shrink_to_fit();
  fCnc.shrink_to_fit();
  fNui.shrink_to_fit();
  energyUni.shrink_to_fit();
  sigmaUniTotal.shrink_to_fit();
  fSigmaOfMts.shrink_to_fit();
  fSigmaOfMtsDop.shrink_to_fit();
  fSigmaUniOfMts.shrink_to_fit();
  fEnergyLocationMts.shrink_to_fit();
  MtNumbers.shrink_to_fit();
  fSigmaMts.shrink_to_fit();
  fQvalueTemp.shrink_to_fit();
  fELinElastic.shrink_to_fit();
  fELinCapture.shrink_to_fit();
  fELinFission.shrink_to_fit();
  fXLinElastic.shrink_to_fit();
  fXLinCapture.shrink_to_fit();
  fXLinFission.shrink_to_fit();
  fXBroadElastic.shrink_to_fit();
  fXBroadCapture.shrink_to_fit();
  fXBroadFission.shrink_to_fit();
  fNbt1.shrink_to_fit();
  fInt1.shrink_to_fit();
  fELinearFile3.shrink_to_fit();
  fXLinearFile3.shrink_to_fit();
  fSigma.shrink_to_fit();
  l.shrink_to_fit();
  fNRS.shrink_to_fit();
  fNRJ.shrink_to_fit();
  fJSM.shrink_to_fit();
  fEr.shrink_to_fit();
  J.shrink_to_fit();
  fGJ.shrink_to_fit();
  fGamma_r.shrink_to_fit();
  fGamma_n.shrink_to_fit();
  fGamma_g.shrink_to_fit();
  fGamma_f.shrink_to_fit();
  fGamma_x.shrink_to_fit();
  fGamma_fa.shrink_to_fit();
  fGamma_fasq.shrink_to_fit();
  fGamma_fb.shrink_to_fit();
  fGamma_fbsq.shrink_to_fit();
  fAt1.shrink_to_fit();
  fAt2.shrink_to_fit();
  fAt3.shrink_to_fit();
  fAt4.shrink_to_fit();
  fBt1.shrink_to_fit();
  fBt2.shrink_to_fit();
  fDet1.shrink_to_fit();
  fDwt1.shrink_to_fit();
  fGrt1.shrink_to_fit();
  fGit1.shrink_to_fit();
  fDef1.shrink_to_fit();
  fDwf1.shrink_to_fit();
  fGrf1.shrink_to_fit();
  fGif1.shrink_to_fit();
  fDec1.shrink_to_fit();
  fDwc1.shrink_to_fit();
  fGrc1.shrink_to_fit();
  fGic1.shrink_to_fit();
  fAmux.shrink_to_fit();
  fAmun.shrink_to_fit();
  fAmug.shrink_to_fit();
  fAmuf.shrink_to_fit();
  fEs.shrink_to_fit();
  fD.shrink_to_fit();
  fGX.shrink_to_fit();
  fGNO.shrink_to_fit();
  fGG.shrink_to_fit();
  fGF.shrink_to_fit();
  fPhiEr.shrink_to_fit();
  fShiftEr.shrink_to_fit();
  fEneTemp.shrink_to_fit();
  fSigTemp.shrink_to_fit();
}
