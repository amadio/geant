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
          cnc.push_back(list->GetLIST(j));
        }
        double ein = 1E-5;
        do {
          double nun = 0;
          for (int i = 0; i < list->GetN1(); i++) {
            nun += cnc[i] * pow(ein, i);
          }
          eintFile1.push_back(ein);
          nutFile1.push_back(nun);
          ein *= 2;
        } while (ein < 21E8);
        cnc.clear();
      } else { // Total fission neutron multiplicity tabulated representation
        //	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        NR                  = tab1->GetN1();
        NP                  = tab1->GetN2();
        //	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
        for (int cr = 0; cr < NR; cr++) {
          nbt1.push_back(tab1->GetNBT(cr));
          int1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < NP; crs++) {
          eintFile1.push_back(tab1->GetX(crs));
          nutFile1.push_back(tab1->GetY(crs));
          //    std::cout<<fE_file1[crs]<<"  "<<fnu_file1[crs]<< std::endl;
        }
        for (int cr = 0; cr < NP - 1; cr++) {
          recursionLinearNuPh(eintFile1[cr], eintFile1[cr + 1], nutFile1[cr], nutFile1[cr + 1], eintFile1, nutFile1);
        }
        TNudyCore::Instance()->Sort(eintFile1, nutFile1);
      }
      eint.push_back(eintFile1);
      nut.push_back(nutFile1);
      eintFile1.clear();
      nutFile1.clear();
      // std::cout <<"eint size "<< eint.size() <<"  "<< nut.size() << std::endl;
      nbt1.clear();
      int1.clear();
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
          nui.push_back(list->GetLIST(i));
          //	 std::cout<<" lambda "<< nui[i] << std::endl;
        }
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        NR                  = tab1->GetN1();
        NP                  = tab1->GetN2();
        //	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
        for (int cr = 0; cr < NR; cr++) {
          nbt1.push_back(tab1->GetNBT(cr));
          int1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < NP; crs++) {
          eindFile1.push_back(tab1->GetX(crs));
          nudFile1.push_back(tab1->GetY(crs));
          //    std::cout<<fE_file1[crs]<<"  "<<fnu_file1[crs]<< std::endl;
        }
        for (int cr = 0; cr < NP - 1; cr++) {
          recursionLinearNuPh(eindFile1[cr], eindFile1[cr + 1], nudFile1[cr], nudFile1[cr + 1], eindFile1, nudFile1);
        }
        TNudyCore::Instance()->Sort(eindFile1, nudFile1);
        eind.push_back(eindFile1);
        nud.push_back(nudFile1);
        lambdaD.push_back(nui);
        eindFile1.clear();
        nudFile1.clear();
        nui.clear();
        nbt1.clear();
        int1.clear();
      } else if (LNU == 2 && LDG == 1) {
      }
    } else if (MT == 456) { // prompt neutron multiplicity

      int LNU = sec->GetL2();
      if (LNU == 1) {
        //      std::cout<<"prompt nu = "<< ZA <<" LNU "<< LNU << std::endl;
        //	TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
        // double nup = list->GetLIST(0);
      } else {
        //	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        NR                  = tab1->GetN1();
        NP                  = tab1->GetN2();
        //	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
        for (int cr = 0; cr < NR; cr++) {
          nbt1.push_back(tab1->GetNBT(cr));
          int1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < NP; crs++) {
          einFile1.push_back(tab1->GetX(crs));
          nuFile1.push_back(tab1->GetY(crs));
          //    std::cout<<fE_file1[crs]<<"  "<<fnu_file1[crs]<< std::endl;
        }
        for (int cr = 0; cr < NP - 1; cr++) {
          recursionLinearNuPh(einFile1[cr], einFile1[cr + 1], nuFile1[cr], nuFile1[cr + 1], einFile1, nuFile1);
        }
        TNudyCore::Instance()->Sort(einFile1, nuFile1);
        nbt1.clear();
        int1.clear();
      }
      einp.push_back(einFile1);
      nup.push_back(nuFile1);
      einFile1.clear();
      nuFile1.clear();
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
          //	      for(unsigned long i = 0; i < eintFile1.size(); i++)
          //		std::cout<< eintFile1[i] <<"  "<<nutFile1[i] << std::endl;
          int n0  = 0;
          int max = eintFile1.size() - 1;
          int mid = 0;
          if (ein <= eintFile1[n0]) {
            n0 = 0;
          } else if (ein >= eintFile1[max]) {
            n0 = max - 1;
          } else {
            while (max - n0 > 1) {
              mid = (n0 + max) / 2;
              if (ein < eintFile1[mid])
                max = mid;
              else
                n0 = mid;
              //      std::cout<<"n0 "<< n0 <<" max "<< max << std::endl;
            }
          }
          //  std::cout<< eintFile1[n0] <<std::endl;
          double nue = TNudyCore::Instance()->LinearInterpolation(eintFile1[n0], nutFile1[n0], eintFile1[n0 + 1],
                                                                  nutFile1[n0 + 1], ein);
          double nu0 = nutFile1[0];
          //   std::cout<<"n0 "<< n0 <<" nu0 "<< nu0 << std::endl;
          EFis -= -1.307 * ein + 8.07 * 1E6 * (nue - nu0); // prompt neutron energy dependence
          einfFile1.push_back(ein);
          heatFile1.push_back(EFis);
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
          einfFile1.push_back(ein);
          heatFile1.push_back(EFis);
          ein *= 2;
        } while (ein < 21E8);
      }
      einFissHeat.push_back(einfFile1);
      fissHeat.push_back(heatFile1);
      einfFile1.clear();
      heatFile1.clear();
    } else if (MT == 460) {
      int LO = sec->GetL1();
      int NG = sec->GetN1();
      //	std::cout<<" Lo "<< LO <<" NG "<< NG <<std::endl;
      if (LO == 1) {
        for (int ng = 0; ng < NG; ng++) {
          TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
          NR                  = tab1->GetN1();
          NP                  = tab1->GetN2();
          for (int cr = 0; cr < NR; cr++) {
            nbt1.push_back(tab1->GetNBT(cr));
            int1.push_back(tab1->GetINT(cr));
            //	std::cout<<"NR = "<< NR <<" NP "<< NP << " NBT1 "<< nbt1[cr]<< " INT1 "<<int1[cr] <<std::endl;
          }
          for (int crs = 0; crs < NP; crs++) {
            einphFile1.push_back(tab1->GetX(crs));
            phFile1.push_back(tab1->GetY(crs));
            //    std::cout<<crs <<"  "<< einphFile1[crs]<<"  "<<phFile1[crs]<< std::endl;
          }
          // linearzation is stopped due to discrete photons which are to be matched with file 12 later
          // for(int cr=0; cr < NP - 1 ; cr ++){
          //  recursionLinearNuPh(einphFile1[cr], einphFile1[cr+1], phFile1[cr], phFile1[cr+1],einphFile1, phFile1);
          //}
          // TNudyCore::Instance()->Sort(einphFile1, phFile1);
          einphFile1.clear();
          phFile1.clear();
          nbt1.clear();
          int1.clear();
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
  // for(unsigned int crs=0; crs < eint.size() ; crs ++)
  // for(unsigned int cr=0; cr < eint[crs].size() ; cr ++)
  // std::cout<< eint[crs][cr] <<"  "<< nut[crs][cr] << std::endl;
}

TNudyEndfNuPh::~TNudyEndfNuPh()
{
  eintFile1.shrink_to_fit();
  nutFile1.shrink_to_fit();
  einFile1.shrink_to_fit();
  nuFile1.shrink_to_fit();
  eindFile1.shrink_to_fit();
  nudFile1.shrink_to_fit();
  einphFile1.shrink_to_fit();
  phFile1.shrink_to_fit();
  einfFile1.shrink_to_fit();
  heatFile1.shrink_to_fit();
  cnc.shrink_to_fit();
  nui.shrink_to_fit();
  nbt1.shrink_to_fit();
  int1.shrink_to_fit();
}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::recursionLinearNuPh(double x1, double x2, double sig1, double sig2, std::vector<double> x,
                                          std::vector<double> sig)
{
  double siga;
  double mid = 0.5 * (x1 + x2);
  if ((sig1 == 0.0 && sig2 == 0.0) || x1 == x2 || x1 < 1E-5 || x2 < 1E-5) return 0;
  //  std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< sig1 <<"  "<< sig2 << std::endl;
  siga = TNudyCore::Instance()->Interpolate(nbt1, int1, NR, x, sig, NP, mid);
  //  std::cout<< mid <<" mid  "<< siga1 <<std::endl;
  //  double sigmid1 = TNudyCore::Instance()->LinearInterpolation(x1,sig1,x2,sig2,mid);
  double sigmid1 = sig1 + (sig2 - sig1) * (mid - x1) / (x2 - x1);
  //  std::cout << mid <<" linear "<< siga <<"  "<< sigmid1 << std::endl;

  if (fabs((siga - sigmid1) / sigmid1) <= sigDiff) {
    return 0;
  }
  x.push_back(mid);
  sig.push_back(siga);
  recursionLinearNuPh(x1, mid, sig1, siga, x, sig);
  recursionLinearNuPh(mid, x2, siga, sig2, x, sig);
  return 0;
}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::GetNuTotal(int ielemid, double energyK)
{
  // std::cout<<"eint size "<< eint.size() <<"  "<< ielemid<<"  "<< energyK << std::endl;
  int min = 0;
  int max = eint[ielemid].size() - 1;
  int mid = 0;
  if (energyK <= eint[ielemid][min])
    min = 0;
  else if (energyK >= eint[ielemid][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < eint[ielemid][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return nut[ielemid][min] +
         (nut[ielemid][min + 1] - nut[ielemid][min]) * (energyK - eint[ielemid][min]) /
             (eint[ielemid][min + 1] - eint[ielemid][min]);
}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::GetNuDelayed(int ielemid, double energyK)
{
  // std::cout<<"eind size "<< eind.size() <<"  "<< ielemid<<"  "<< energyK << std::endl;
  int min = 0;
  int max = eind[ielemid].size() - 1;
  int mid = 0;
  if (energyK <= eind[ielemid][min])
    min = 0;
  else if (energyK >= eind[ielemid][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < eind[ielemid][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return nud[ielemid][min] +
         (nud[ielemid][min + 1] - nud[ielemid][min]) * (energyK - eind[ielemid][min]) /
             (eind[ielemid][min + 1] - eind[ielemid][min]);
}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::GetNuPrompt(int ielemid, double energyK)
{
  // std::cout<<"einp size "<< einp.size() <<"  "<< ielemid<<"  "<< energyK << std::endl;
  int min = 0;
  int max = einp[ielemid].size() - 1;
  int mid = 0;
  if (energyK <= einp[ielemid][min])
    min = 0;
  else if (energyK >= einp[ielemid][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < einp[ielemid][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return nup[ielemid][min] +
         (nup[ielemid][min + 1] - nup[ielemid][min]) * (energyK - einp[ielemid][min]) /
             (einp[ielemid][min + 1] - einp[ielemid][min]);
}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::GetFissHeat(int ielemid, double energyK)
{
  // std::cout<<"fissHeat size "<< einFissHeat.size() <<"  "<< ielemid<<"  "<< energyK << std::endl;
  int min = 0;
  int max = einFissHeat[ielemid].size() - 1;
  int mid = 0;
  if (energyK <= einFissHeat[ielemid][min])
    min = 0;
  else if (energyK >= einFissHeat[ielemid][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < einFissHeat[ielemid][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return fissHeat[ielemid][min] +
         (fissHeat[ielemid][min + 1] - fissHeat[ielemid][min]) * (energyK - einFissHeat[ielemid][min]) /
             (einFissHeat[ielemid][min + 1] - einFissHeat[ielemid][min]);
}
//____________________________________________________________________________________________________________________
double TNudyEndfNuPh::GetLambdaD(int ielemid, int time)
{
  // std::cout<<"eint size "<< eint.size() <<"  "<< ielemid<<"  "<< energyK << std::endl;
  if (time > (int)lambdaD[ielemid].size() - 1) return -99.0;
  return lambdaD[ielemid][time];
}
