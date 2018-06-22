// 	This class is reconstructing probability tables for nu multiplicity, fission heat and photon energy distribution
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 24, 2016

#include "TList.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfNuPh.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfNuPh)
#endif

    TNudyEndfNuPh::TNudyEndfNuPh()
{
}

//______________________________________________________________________________
TNudyEndfNuPh::TNudyEndfNuPh(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    //	std::cout <<"file  "<<sec->GetMT() <<std::endl;
    TIter recIter(sec->GetRecords());
    //	double ZA   = sec->GetC1();
    //	double AWR  = sec->GetC2();
    int MT = sec->GetMT();
    if (MT == 452) { // Total fission neutron multiplicity polynomial expansion
      int LNU = sec->GetL2();
      if (LNU == 1) {
        //      std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
        TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
        //	      std::cout<<list->GetNPL()<<"  "<< list->GetN1() << std::endl;
        for (int j = 0; j < list->GetN1(); j++) {
          fCnc.push_back(list->GetLIST(j));
        }
        double ein = 1E-5;
        do {
          double nun = 0;
          for (int i = 0; i < list->GetN1(); i++) {
            nun += fCnc[i] * pow(ein, i);
          }
          fEintFile1.push_back(ein);
          fNutFile1.push_back(nun);
          ein *= 2;
        } while (ein < 21E8);
        fCnc.clear();
      } else { // Total fission neutron multiplicity tabulated representation
        //	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        fNR                 = tab1->GetN1();
        fNP                 = tab1->GetN2();
        //	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fEintFile1.push_back(tab1->GetX(crs));
          fNutFile1.push_back(tab1->GetY(crs));
          //    std::cout<<fE_file1[crs]<<"  "<<fnu_file1[crs]<< std::endl;
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
      // std::cout <<"fEint size "<< fEint.size() <<"  "<< fNut.size() << std::endl;
      fNbt1.clear();
      fInt1.clear();
    } else if (MT == 455) { // delayed neutron multiplicity
      int LDG = sec->GetL1();
      int LNU = sec->GetL2();
      //	 std::cout<<" LNU "<< LNU <<"  LDG "<< LDG << std::endl;
      if (LNU == 1 && LDG == 0) {

      } else if (LNU == 1 && LDG == 1) {

      } else if (LNU == 2 && LDG == 0) {
        TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
        int NNF             = list->GetNPL();
        for (int i = 0; i < NNF; i++) {
          fNui.push_back(list->GetLIST(i));
          //	 std::cout<<" lambda "<< fNui[i] << std::endl;
        }
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        fNR                 = tab1->GetN1();
        fNP                 = tab1->GetN2();
        //	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fEindFile1.push_back(tab1->GetX(crs));
          fNudFile1.push_back(tab1->GetY(crs));
          //    std::cout<<fE_file1[crs]<<"  "<<fnu_file1[crs]<< std::endl;
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
    } else if (MT == 456) { // prompt neutron multiplicity

      int LNU = sec->GetL2();
      if (LNU == 1) {
        //      std::cout<<"prompt nu = "<< ZA <<" LNU "<< LNU << std::endl;
        //	TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
        // double fNup = list->GetLIST(0);
      } else {
        //	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        fNR                 = tab1->GetN1();
        fNP                 = tab1->GetN2();
        //	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fEinFile1.push_back(tab1->GetX(crs));
          fNuFile1.push_back(tab1->GetY(crs));
          //    std::cout<<fE_file1[crs]<<"  "<<fnu_file1[crs]<< std::endl;
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
    } else if (MT == 458) {
      TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
      int NPLY            = list->GetL2();
      //	  std::cout<< "NPLY "<< NPLY << std::endl;
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
          //  std::cout<< fEintFile1[n0] <<std::endl;
          double nue = TNudyCore::Instance()->LinearInterpolation(fEintFile1[n0], fNutFile1[n0], fEintFile1[n0 + 1],
                                                                  fNutFile1[n0 + 1], ein);
          double nu0 = fNutFile1[0];
          //   std::cout<<"n0 "<< n0 <<" nu0 "<< nu0 << std::endl;
          EFis -= -1.307 * ein + 8.07 * 1E6 * (nue - nu0); // prompt neutron energy dependence
          fEinfFile1.push_back(ein);
          fHeatFile1.push_back(EFis);
          //      std::cout<<"ein "<< ein <<"  "<< EFis <<std::endl;
          ein *= 2;
        } while (ein < 21E8);
        //	std::cout<< "NPLY "<< NPLY <<" ER "<< ER <<" ET "<< ET << std::endl;
      } else {
        double c0[9 * (NPLY + 1)], c1[9 * (NPLY + 1)];
        for (int i = 0; i < 9 * (NPLY + 1); i++) {
          c0[i] = list->GetLIST(i);
          //	     std::cout <<"c0  "<< c0[i] << std::endl;
        }
        for (int i = 9 * (NPLY + 1); i < 18 * (NPLY + 1); i++) {
          c1[i - 9 * (NPLY + 1)] = list->GetLIST(i);
          //	     std::cout <<"c1  "<< c1[i-9*(NPLY+1)] << std::endl;
        }
        double ein = 1E-5;
        do {
          double EFis = 0;
          for (int i = 0; i < 9 * NPLY / 2 - 2; i++) {
            EFis += c0[i * 2] + c1[i * 2] * ein;
            //      std::cout<<"ein "<< ein <<"  "<< EFis <<std::endl;
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
    } else if (MT == 460) {
      int LO = sec->GetL1();
      int NG = sec->GetN1();
      //	std::cout<<" Lo "<< LO <<" NG "<< NG <<std::endl;
      if (LO == 1) {
        for (int ng = 0; ng < NG; ng++) {
          TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
          fNR                 = tab1->GetN1();
          fNP                 = tab1->GetN2();
          for (int cr = 0; cr < fNR; cr++) {
            fNbt1.push_back(tab1->GetNBT(cr));
            fInt1.push_back(tab1->GetINT(cr));
            //	std::cout<<"fNR = "<< fNR <<" fNP "<< fNP << " NBT1 "<< fNbt1[cr]<< " INT1 "<<fInt1[cr] <<std::endl;
          }
          for (int crs = 0; crs < fNP; crs++) {
            fEinPhFile1.push_back(tab1->GetX(crs));
            fPhFile1.push_back(tab1->GetY(crs));
            //    std::cout<<crs <<"  "<< fEinPhFile1[crs]<<"  "<<fPhFile1[crs]<< std::endl;
          }
          // linearzation is stopped due to discrete photons which are to be matched with file 12 later
          // for(int cr=0; cr < fNP - 1 ; cr ++){
          //  RecursionLinearNuPh(fEinPhFile1[cr], fEinPhFile1[cr+1], fPhFile1[cr], fPhFile1[cr+1],fEinPhFile1,
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
        // std::cout<< "Photon NNF "<< NNF <<" LO "<< LO << std::endl;
        // double lambda[NNF];
        for (int i = 0; i < NNF; i++) {
          // lambda[i] = list->GetLIST(i);
          //	     std::cout <<"lambda  "<< lambda[i] << std::endl;
        }
      }
    }
  }
  // for(unsigned int crs=0; crs < fEint.size() ; crs ++)
  // for(unsigned int cr=0; cr < fEint[crs].size() ; cr ++)
  // std::cout<< fEint[crs][cr] <<"  "<< fNut[crs][cr] << std::endl;
}

TNudyEndfNuPh::~TNudyEndfNuPh(){}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::RecursionLinearNuPh(double x1, double x2, double sig1, double sig2, std::vector<double> x,
                                          std::vector<double> sig)
{
  double siga;
  double mid = 0.5 * (x1 + x2);
  if ((sig1 == 0.0 && sig2 == 0.0) || x1 == x2 || x1 < 1E-5 || x2 < 1E-5) return 0;
  //  std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< sig1 <<"  "<< sig2 << std::endl;
  siga = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, x, sig, fNP, mid);
  //  std::cout<< mid <<" mid  "<< siga1 <<std::endl;
  //  double sigmid1 = TNudyCore::Instance()->LinearInterpolation(x1,sig1,x2,sig2,mid);
  double sigmid1 = sig1 + (sig2 - sig1) * (mid - x1) / (x2 - x1);
  //  std::cout << mid <<" linear "<< siga <<"  "<< sigmid1 << std::endl;

  if (fabs((siga - sigmid1) / sigmid1) <= fSigDiff) {
    return 0;
  }
  x.push_back(mid);
  sig.push_back(siga);
  RecursionLinearNuPh(x1, mid, sig1, siga, x, sig);
  RecursionLinearNuPh(mid, x2, siga, sig2, x, sig);
  return 0;
}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::GetNuTotal(int ielemid, double energyK)
{
  // std::cout<<"fEint size "<< fEint.size() <<"  "<< ielemid<<"  "<< energyK << std::endl;
  int min = 0;
  int max = fEint[ielemid].size() - 1;
  int mid = 0;
  if (energyK <= fEint[ielemid][min])
    min = 0;
  else if (energyK >= fEint[ielemid][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEint[ielemid][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return fNut[ielemid][min] +
         (fNut[ielemid][min + 1] - fNut[ielemid][min]) * (energyK - fEint[ielemid][min]) /
             (fEint[ielemid][min + 1] - fEint[ielemid][min]);
}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::GetNuDelayed(int ielemid, double energyK)
{
  // std::cout<<"fEind size "<< fEind.size() <<"  "<< ielemid<<"  "<< energyK << std::endl;
  int min = 0;
  int max = fEind[ielemid].size() - 1;
  int mid = 0;
  if (energyK <= fEind[ielemid][min])
    min = 0;
  else if (energyK >= fEind[ielemid][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEind[ielemid][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return fNud[ielemid][min] +
         (fNud[ielemid][min + 1] - fNud[ielemid][min]) * (energyK - fEind[ielemid][min]) /
             (fEind[ielemid][min + 1] - fEind[ielemid][min]);
}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::GetNuPrompt(int ielemid, double energyK)
{
  // std::cout<<"fEinp size "<< fEinp.size() <<"  "<< ielemid<<"  "<< energyK << std::endl;
  int min = 0;
  int max = fEinp[ielemid].size() - 1;
  int mid = 0;
  if (energyK <= fEinp[ielemid][min])
    min = 0;
  else if (energyK >= fEinp[ielemid][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEinp[ielemid][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return fNup[ielemid][min] +
         (fNup[ielemid][min + 1] - fNup[ielemid][min]) * (energyK - fEinp[ielemid][min]) /
             (fEinp[ielemid][min + 1] - fEinp[ielemid][min]);
}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::GetFissHeat(int ielemid, double energyK)
{
  // std::cout<<"fFissHeat size "<< fEinFissHeat.size() <<"  "<< ielemid<<"  "<< energyK << std::endl;
  int min = 0;
  int max = fEinFissHeat[ielemid].size() - 1;
  int mid = 0;
  if (energyK <= fEinFissHeat[ielemid][min])
    min = 0;
  else if (energyK >= fEinFissHeat[ielemid][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEinFissHeat[ielemid][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return fFissHeat[ielemid][min] +
         (fFissHeat[ielemid][min + 1] - fFissHeat[ielemid][min]) * (energyK - fEinFissHeat[ielemid][min]) /
             (fEinFissHeat[ielemid][min + 1] - fEinFissHeat[ielemid][min]);
}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::GetLambdaD(int ielemid, int time)
{
  // std::cout<<"fEint size "<< fEint.size() <<"  "<< ielemid<<"  "<< energyK << std::endl;
  if (time > (int)fLambdaD[ielemid].size() - 1) return -99.0;
  return fLambdaD[ielemid][time];
}
