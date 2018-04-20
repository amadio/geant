// 	This class is reconstructing probability tables for Energy distribution
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 24, 2016

#include "TList.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfEnergy.h"
#include "Math/SpecFuncMathMore.h"
#include "TMath.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfEnergy)
#include "TRandom3.h"
#endif

    TNudyEndfEnergy::TNudyEndfEnergy()
{
}

//______________________________________________________________________________
TNudyEndfEnergy::TNudyEndfEnergy(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  int mt455 = 1000;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    for (int k = 0; k < sec->GetN1(); k++) {
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      int MT              = sec->GetMT();
      if (MT == 455) {
        MT = mt455;
        mt455++;
      }
      fMtNumbers.push_back(MT);
      int LF = tab1->GetL2();
      // std::cout << " LF = " << LF << " MT " << MT << "  k " << k << "  " << sec->GetN1() << std::endl;
      fNR = tab1->GetN1();
      fNP = tab1->GetN2();
      //****************************************************************************
      // arbitrary tabulated function
      if (LF == 1) {
        //         std::cout << " fNR = " << fNR << " fNP " << fNP << std::endl;
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
        }
        TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
        fNr2                = tab2->GetN1();
        fNp2                = tab2->GetN2();

        for (int cr = 0; cr < fNr2; cr++) {
          fNbt2.push_back(tab2->GetNBT(cr));
          fInt2.push_back(tab2->GetINT(cr));
        }
        for (int cr = 0; cr < fNp2; cr++) {
          TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
          fEin.push_back(tab12->GetC2());
          //  	  std::cout<<"energy "<< tab12->GetC2() << std::endl;
          fNr3 = tab12->GetNR();
          fNp3 = tab12->GetNP();
          for (int i = 0; i < fNr3; i++) {
            fNbt3.push_back(tab12->GetNBT(i));
            fInt3.push_back(tab12->GetINT(i));
          }
          for (int crs = 0; crs < fNp3; crs++) {
            fEnergyFile5.push_back(tab12->GetX(crs));
            fEnergyPdfFile5.push_back(tab12->GetY(crs));
          }
          for (int cr = 0; cr < fNp3 - 1; cr++) {
            //               std::cout << fEnergyFile5[cr] <<"  "<< fEnergyPdfFile5[cr] << std::endl;
            RecursionLinearFile5Prob(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr],
                                     fEnergyPdfFile5[cr + 1]);
          }
          //           std::cout<<"linearization complete "<< std::endl;
          FillPdf1D();
          fNbt3.clear();
          fInt3.clear();
        }
        fNbt1.clear();
        fInt1.clear();
        fNbt2.clear();
        fInt2.clear();
        fE1.clear();
        //****************************************************************************
        // general evaporation spectrum
      } else if (LF == 5) {
        double u = tab1->GetC1();
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
        }
        TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
        fNr2                 = tab11->GetN1();
        fNp2                 = tab11->GetN2();
        for (int cr = 0; cr < fNr2; cr++) {
          fNbt2.push_back(tab11->GetNBT(cr));
          fInt2.push_back(tab11->GetINT(cr));
        }
        for (int crs = 0; crs < fNp2; crs++) {
          fE2.push_back(tab11->GetX(crs));
          fP2.push_back(tab11->GetY(crs));
          // fEin.push_back(tab11->GetX(crs));
        }
        TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
        fNr3                 = tab12->GetN1();
        fNp3                 = tab12->GetN2();
        for (int cr = 0; cr < fNr3; cr++) {
          fNbt3.push_back(tab12->GetNBT(cr));
          fInt3.push_back(tab12->GetINT(cr));
        }
        for (int crs = 0; crs < fNp3; crs++) {
          fE3.push_back(tab12->GetX(crs));
          fP3.push_back(tab12->GetY(crs));
        }
        double energy, eout;
        for (int i = 0; i < fNp2; i++) {
          energy = fE2[i];
          // std::cout<<"energy "<<energy<<std::endl;
          double sumprob         = 0.0;
          eout                   = fE3[1] / 100;
          if (eout < 1E-05) eout = 1E-05;
          do {
            double gx     = 0.0;
            double pe     = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
            double thetae = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
            if (thetae > 0.0)
              gx       = TNudyCore::Instance()->Interpolate(fNbt3, fInt3, fNr3, fE3, fP3, fNp3, eout / thetae);
            double ppe = 0.0;
            ppe        = pe * gx;
            if (ppe > 0.0) {
              fEnergyFile5.push_back(eout);
              fEnergyPdfFile5.push_back(ppe);
              sumprob += ppe;
            }
            // std:: cout<<"i = "<< i <<" pe "<< pe <<" thetae "<< thetae <<"  prob "<< ppe <<"  fEin "<< energy <<"
            // eout "<< eout <<"  "<< u << std::endl;
            eout *= 2;
          } while (eout < energy - u);
          if (sumprob > 0.0) fEin.push_back(energy);
          int size = fEnergyFile5.size();
          for (int cr = 0; cr < size - 1; cr++) {
            //            std::cout << fEnergyFile5[cr] <<"  "<< fEnergyPdfFile5[cr] << std::endl;
            RecursionLinearFile5GenEva(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr],
                                       fEnergyPdfFile5[cr + 1], energy);
          }
          //           for(unsigned int cr = 0 ; cr < fEnergyFile5.size() ; cr ++){
          //             std::cout <<"fill1d "<< energy <<"  "<< fEnergyFile5[cr] <<"  "<< fEnergyPdfFile5[cr] <<
          //             std::endl;
          // 	  }
          FillPdf1D();
        }
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fNbt2.clear();
        fInt2.clear();
        fE2.clear();
        fP2.clear();
        fNbt3.clear();
        fInt3.clear();
        fE3.clear();
        fP3.clear();
        //****************************************************************************
        // simple maxwellian fission spectrum
      } else if (LF == 7) {
        double u = tab1->GetC1();
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
        }
        TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
        fNr2                 = tab11->GetN1();
        fNp2                 = tab11->GetN2();
        for (int cr = 0; cr < fNr2; cr++) {
          fNbt2.push_back(tab11->GetNBT(cr));
          fInt2.push_back(tab11->GetINT(cr));
        }
        for (int crs = 0; crs < fNp2; crs++) {
          fE2.push_back(tab11->GetX(crs));
          fP2.push_back(tab11->GetY(crs));
          double sq  = (fE2[crs] - u) / fP2[crs];
          double sqt = sqrt(sq);
          INorm.push_back(pow(fP2[crs], 1.5) * (0.5 * 1.7724538529055 * erf(sqt) - sqt * exp(-sq)));
          // fEin.push_back(tab11->GetX(crs));
        }
        double energy, eout;
        for (int i = 0; i < fNp2; i++) {
          energy                 = fE2[i];
          eout                   = fE2[0] / 100;
          if (eout < 1E-05) eout = 1E-05;
          double sumprob         = 0.0;
          do {
            double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
            double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
            double thetae    = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
            double ppe       = 0.0;
            if (I > 0.0) ppe = pe * sqrt(eout) * exp(-eout / thetae) / I;
            if (ppe > 0.0) {
              fEnergyFile5.push_back(eout);
              fEnergyPdfFile5.push_back(ppe);
              sumprob += ppe;
            }
            //    std:: cout<<"I = "<< I <<" pe "<< pe <<" thetae "<< thetae <<"  prob "<< ppe <<"  fEin "<< energy <<"
            //    eout "<< eout << std::endl;
            eout *= 2;
          } while (eout < energy - u);
          if (sumprob > 0.0) fEin.push_back(energy);
          int size = fEnergyFile5.size();
          for (int cr = 0; cr < size - 1; cr++) {
            RecursionLinearFile5Maxwell(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr],
                                        fEnergyPdfFile5[cr + 1], energy);
          }
          FillPdf1D();
        }
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fNbt2.clear();
        fInt2.clear();
        fE2.clear();
        fP2.clear();
        INorm.clear();
        ///////////////////////////////////////////////////////////////////////////////////
        // evaporation spectrum
      } else if (LF == 9) {
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fNbt2.clear();
        fInt2.clear();
        fE2.clear();
        fP2.clear();
        INorm.clear();
        double u = tab1->GetC1();
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
        }
        TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
        fNr2                 = tab11->GetN1();
        fNp2                 = tab11->GetN2();
        for (int cr = 0; cr < fNr2; cr++) {
          fNbt2.push_back(tab11->GetNBT(cr));
          fInt2.push_back(tab11->GetINT(cr));
        }
        for (int crs = 0; crs < fNp2; crs++) {
          fE2.push_back(tab11->GetX(crs));
          fP2.push_back(tab11->GetY(crs));
          INorm.push_back(fP2[crs] * fP2[crs] *
                          (1. - exp(-(fE2[crs] - u) / fP2[crs]) * (1.0 + (fE2[crs] - u) / fP2[crs])));
          // fEin.push_back(tab11->GetX(crs));
        }
        double energy, eout;
        for (int i = 0; i < fNp2; i++) {
          energy                 = fE2[i];
          eout                   = fE2[0] / 100;
          if (eout < 1E-05) eout = 1E-05;
          double sumprob         = 0.0;
          do {
            double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
            double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
            double thetae    = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
            double ppe       = 0.0;
            if (I > 0.0) ppe = pe * eout * exp(-eout / thetae) / I;
            if (ppe > 0.0) {
              fEnergyFile5.push_back(eout);
              fEnergyPdfFile5.push_back(ppe);
              sumprob += ppe;
              // std::cout << eout <<"  "<< ppe << std::endl;
            }
            //    std:: cout<<"I = "<< I <<" pe "<< pe <<" thetae "<< thetae <<"  prob "<< ppe <<"  fEin "<< energy <<"
            //    eout "<< eout << std::endl;
            eout *= 2;
          } while (eout < u);
          if (sumprob > 0.0) fEin.push_back(energy);
          if (fEin.size() > 1 && fEin[fEin.size() - 1] == fEin[fEin.size() - 2]) {
            fEin.erase(fEin.begin() + fEin.size() - 1);
            fMtNumbers.erase(fMtNumbers.begin() + fMtNumbers.size() - 1);
            fEnergyFile5.clear();
            fEnergyPdfFile5.clear();
            continue;
          }
          // std::cout<<energy <<"  "<<sumprob<<std::endl;
          int size = fEnergyFile5.size();
          for (int cr = 0; cr < size - 1; cr++) {
            RecursionLinearFile5Maxwell(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr],
                                        fEnergyPdfFile5[cr + 1], energy);
          }
          FillPdf1D();
        }
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fNbt2.clear();
        fInt2.clear();
        fE2.clear();
        fP2.clear();
        INorm.clear();
        /////////////////////////////////////////////////////////////////////////////////////////
        // energy dependent watt spectrum
      } else if (LF == 11) {
        double u = tab1->GetC1();

        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
        }
        TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
        fNr2                 = tab11->GetN1();
        fNp2                 = tab11->GetN2();
        for (int cr = 0; cr < fNr2; cr++) {
          fNbt2.push_back(tab11->GetNBT(cr));
          fInt2.push_back(tab11->GetINT(cr));
        }

        for (int crs = 0; crs < fNp2; crs++) {
          fE2.push_back(tab11->GetX(crs));
          fP2.push_back(tab11->GetY(crs));
          // fEin.push_back(tab11->GetX(crs));
        }
        TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
        fNr3                 = tab12->GetN1();
        fNp3                 = tab12->GetN2();
        for (int cr = 0; cr < fNr3; cr++) {
          fNbt3.push_back(tab12->GetNBT(cr));
          fInt3.push_back(tab12->GetINT(cr));
        }
        for (int crs = 0; crs < fNp3; crs++) {
          fE3.push_back(tab12->GetX(crs));
          fP3.push_back(tab12->GetY(crs));

          double a   = fP2[crs];
          double b   = fP3[crs];
          double eua = (fE3[crs] - u) / a;

          INorm.push_back(0.5 * sqrt(0.25 * PI * a * a * a * b) * exp(0.25 * a * b) *
                              (erf(sqrt(eua) - sqrt(0.25 * a * b)) + erf(sqrt(eua) + sqrt(0.25 * a * b))) -
                          a * exp(-eua) * sinh(sqrt(b * (fE3[crs] - u))));
        }
        double energy, eout;
        for (int i = 0; i < fNp2; i++) {
          energy                 = fE2[i];
          eout                   = fE2[0] / 100;
          if (eout < 1E-05) eout = 1E-05;
          double sumprob         = 0.0;
          do {
            double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
            double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
            double ae        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
            double be        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE3, fP3, fNp3, energy);
            double ppe       = 0.0;
            if (I > 0.0) ppe = pe * sinh(sqrt(be * eout)) * exp(-eout / ae) / I;
            if (ppe > 0.0) {
              fEnergyFile5.push_back(eout);
              fEnergyPdfFile5.push_back(ppe);
              sumprob += ppe;
            }
            //    std:: cout<<"I = "<< I <<" pe "<< pe <<" ae "<< ae <<"  prob "<< ppe <<"  fEin "<< energy <<"  eout
            //    "<<
            //    eout << std::endl;
            eout *= 2;
          } while (eout < energy - u);
          if (sumprob > 0.0) fEin.push_back(energy);
          int size = fEnergyFile5.size();
          for (int cr = 0; cr < size - 1; cr++) {
            RecursionLinearFile5Watt(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr],
                                     fEnergyPdfFile5[cr + 1], energy);
          }
          FillPdf1D();
        }
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fNbt2.clear();
        fInt2.clear();
        fE2.clear();
        fP2.clear();
        fNbt3.clear();
        fInt3.clear();
        fE3.clear();
        fP3.clear();
        INorm.clear();
        /////////////////////////////////////////////////////////////////////////////////////////
      } else if (LF == 12) {
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
        }
        TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
        double efl           = tab11->GetC1();
        double efh           = tab11->GetC2();

        fNr2 = tab11->GetN1();
        fNp2 = tab11->GetN2();
        for (int cr = 0; cr < fNr2; cr++) {
          fNbt2.push_back(tab11->GetNBT(cr));
          fInt2.push_back(tab11->GetINT(cr));
        }
        for (int crs = 0; crs < fNp2; crs++) {
          fE2.push_back(tab11->GetX(crs));
          fP2.push_back(tab11->GetY(crs));
          fEin.push_back(tab11->GetX(crs));
        }
        double energy, eout;
        for (int i = 0; i < fNp2; i++) {
          energy                 = fE2[i];
          eout                   = fE2[0] / 100;
          if (eout < 1E-05) eout = 1E-05;
          do {
            double pe  = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
            double tm  = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
            double u1l = (sqrt(eout) - sqrt(efl)) * (sqrt(eout) - sqrt(efl)) / tm;
            double u2l = (sqrt(eout) + sqrt(efl)) * (sqrt(eout) + sqrt(efl)) / tm;
            double u1h = (sqrt(eout) - sqrt(efh)) * (sqrt(eout) - sqrt(efh)) / tm;
            double u2h = (sqrt(eout) + sqrt(efh)) * (sqrt(eout) + sqrt(efh)) / tm;
            // std::cout<<" u1l "<< u1l <<" u2l "<< u2l <<" u1h "<< u1h <<" u2h "<< u2h << std::endl;
            double e1ul = ROOT::Math::expint(u1l);
            double e2ul = ROOT::Math::expint(u2l);
            double e1uh = ROOT::Math::expint(u1h);
            double e2uh = ROOT::Math::expint(u2h);
            // std::cout<<" e1ul "<< e1ul <<" e2ul "<< e2ul <<" e1uh "<< e1uh <<" e2uh "<< e2uh << std::endl;
            double a1      = 1.5;
            double gamau1l = TMath::Gamma(a1, u1l);
            double gamau2l = TMath::Gamma(a1, u2l);
            double gamau1h = TMath::Gamma(a1, u1h);
            double gamau2h = TMath::Gamma(a1, u2h);
            // std::cout<<" gamau2l "<< gamau2l <<" gamau1l "<< gamau1l <<" gamau2h "<< gamau2h <<" gamau1h "<< gamau1h
            // << std::endl;
            double gl = (1. / (3 * sqrt(efl * tm))) * (pow(u2l, 1.5) * e2ul - pow(u1l, 1.5) * e1ul + gamau2l - gamau1l);
            double gh = (1. / (3 * sqrt(efh * tm))) * (pow(u2h, 1.5) * e2uh - pow(u1h, 1.5) * e1uh + gamau2h - gamau1h);
            /*
            double alp = sqrt(tm);
            double bet = sqrt(efl);
            double a1 = (sqrt(eout) + bet) * (sqrt(eout) + bet)/tm;
            double b1 = (sqrt(30E+6) + bet) * (sqrt(3E7) + bet)/tm;
            double a2 = (sqrt(eout) - bet) * (sqrt(eout) - bet)/tm;
            double b2 = (sqrt(30E+6) - bet) * (sqrt(3E7) - bet)/tm;
            */
            double ppe = 0.5 * pe * (gl + gh);
            fEnergyFile5.push_back(eout);
            fEnergyPdfFile5.push_back(ppe);
            //    std:: cout<< eout <<"  "<< ppe <<"  "<< energy <<"  "<< eout << std::endl;
            //    std:: cout<<"I = "<< I <<" pe "<< pe <<" ae "<< ae <<"  prob "<< ppe <<"  fEin "<< energy <<"  eout
            //    "<<
            //    eout << std::endl;
            eout *= 2;
          } while (eout < fE1[fNP - 1]);
          FillPdf1D();
        }
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fNbt2.clear();
        fInt2.clear();
        fE2.clear();
        fP2.clear();
      }
      // fFrac2D.push_back(fP1);
      // fP1.clear();
      /*
            fEin2D.push_back(fEin);
            fEne3D.push_back(fEne2D);
            fPdf3D.push_back(fPdf2D);
            fCdf3D.push_back(fCdf2D);
            fEin.clear();
            fEne2D.clear();
            fPdf2D.clear();
            fCdf2D.clear();
            */
    }
    fEin2D.push_back(fEin);
    fFrac2D.push_back(fP1);
    fEne3D.push_back(fEne2D);
    fPdf3D.push_back(fPdf2D);
    fCdf3D.push_back(fCdf2D);
    fEin.clear();
    fEne2D.clear();
    fPdf2D.clear();
    fCdf2D.clear();
    fP1.clear();
  }
  fMt5Values.push_back(fMtNumbers);
  fEnergy5OfMts.push_back(fEin2D);
  fFraction5OfMts.push_back(fFrac2D);
  fEnergyOut5OfMts.push_back(fEne3D);
  fEnergyPdf5OfMts.push_back(fPdf3D);
  fEnergyCdf5OfMts.push_back(fCdf3D);
  fMtNumbers.clear();
  fEin2D.clear();
  fEne3D.clear();
  fPdf3D.clear();
  fCdf3D.clear();
  fFrac2D.clear();
  /*
  for(unsigned long i = 0; i < fEnergy5OfMts[0].size() ; i++){
      std::cout <<" mt "<<fMt5Values[0][i]<<" size "<< fEnergy5OfMts[0][i].size() << std::endl;
    for(unsigned long j = 0; j < fEnergy5OfMts[0][i].size(); j++){
      std::cout << fEnergy5OfMts[0][i][j] <<"  "<< fFraction5OfMts[0][i][j] << std::endl;
     // for(unsigned long k =0; k < fEnergyPdf5OfMts[i][j].size()/2; k++){
  //std::cout << fEnergyPdf5OfMts[i][j][2*k] <<"  "<< fEnergyPdf5OfMts[i][j][2*k + 1] <<"  "<<
  fEnergyCdf5OfMts[i][j][2*k +
  1]<< std::endl;
      //}
    }
  }
  */
}

TNudyEndfEnergy::~TNudyEndfEnergy()
{
  fMtNumbers.shrink_to_fit();
  fE1.shrink_to_fit();
  fP1.shrink_to_fit();
  fE2.shrink_to_fit();
  fP2.shrink_to_fit();
  fE3.shrink_to_fit();
  fP3.shrink_to_fit();
  INorm.shrink_to_fit();
  fNbt1.shrink_to_fit();
  fInt1.shrink_to_fit();
  fNbt2.shrink_to_fit();
  fInt2.shrink_to_fit();
  fNbt3.shrink_to_fit();
  fInt3.shrink_to_fit();
  fEnergyFile5.shrink_to_fit();
  fEnergyPdfFile5.shrink_to_fit();
  fEnergyCdfFile5.shrink_to_fit();
  fEin.shrink_to_fit();
  fEneE.shrink_to_fit();
  fCdf.shrink_to_fit();
  fPdf.shrink_to_fit();
  fEne2D.shrink_to_fit();
  fFrac2D.shrink_to_fit();
  fCdf2D.shrink_to_fit();
  fPdf2D.shrink_to_fit();
  fEin2D.shrink_to_fit();
  fEne3D.shrink_to_fit();
  fCdf3D.shrink_to_fit();
  fPdf3D.shrink_to_fit();
}

double TNudyEndfEnergy::RecursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  pdf            = TNudyCore::Instance()->Interpolate(fNbt3, fInt3, fNr3, fEnergyFile5, fEnergyPdfFile5, fNp3, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5Prob(x1, mid, pdf1, pdf);
  RecursionLinearFile5Prob(mid, x2, pdf, pdf2);
  return 0;
}

double TNudyEndfEnergy::RecursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy)
{
  //   std::cout << std::setprecision(12) << x1 <<"  "<< std::setprecision(12) << x2 <<"  "<<
  //   std::setprecision(12) << pdf1 <<"  "<< std::setprecision(12) << pdf2 << std::endl;
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  double gx            = 0.0;
  double pe            = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
  double thetae        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
  if (thetae > 0.0) gx = TNudyCore::Instance()->Interpolate(fNbt3, fInt3, fNr3, fE3, fP3, fNp3, mid / thetae);
  pdf                  = 0.0;
  pdf                  = pe * gx;
  double pdfmid1       = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5GenEva(x1, mid, pdf1, pdf, energy);
  RecursionLinearFile5GenEva(mid, x2, pdf, pdf2, energy);
  return 0;
}

double TNudyEndfEnergy::RecursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
  double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
  double thetae    = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
  pdf              = 0.0;
  if (I > 0.0) pdf = pe * sqrt(mid) * exp(-mid / thetae) / I;
  double pdfmid1   = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5Maxwell(x1, mid, pdf1, pdf, energy);
  RecursionLinearFile5Maxwell(mid, x2, pdf, pdf2, energy);
  return 0;
}

double TNudyEndfEnergy::RecursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  double pe        = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fE1, fP1, fNP, energy);
  double I         = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, INorm, fNp2, energy);
  double ae        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE2, fP2, fNp2, energy);
  double be        = TNudyCore::Instance()->Interpolate(fNbt2, fInt2, fNr2, fE3, fP3, fNp3, energy);
  pdf              = 0.0;
  if (I > 0.0) pdf = pe * sinh(sqrt(be * mid)) * exp(-mid / ae) / I;
  double pdfmid1   = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  fEnergyFile5.push_back(mid);
  fEnergyPdfFile5.push_back(pdf);
  RecursionLinearFile5Watt(x1, mid, pdf1, pdf, energy);
  RecursionLinearFile5Watt(mid, x2, pdf, pdf2, energy);
  return 0;
}
void TNudyEndfEnergy::FillPdf1D()
{
  TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
  TNudyCore::Instance()->ThinningDuplicate(fEnergyFile5, fEnergyPdfFile5);
  TNudyCore::Instance()->cdfGenerateT(fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5);
  for (unsigned long i = 0; i < fEnergyFile5.size(); i++) {
    if (fEnergyCdfFile5[i] > 0) {
      fEneE.push_back(fEnergyFile5[i]);
      fPdf.push_back(fEnergyPdfFile5[i]);
      fCdf.push_back(fEnergyCdfFile5[i]);
      //         std::cout << fEnergyFile5[i] <<"  "<< fEnergyPdfFile5[i] <<"  "<< fEnergyCdfFile5[i] << std::endl ;
    }
  }
  fEne2D.push_back(fEneE);
  fPdf2D.push_back(fPdf);
  fCdf2D.push_back(fCdf);
  fEnergyFile5.clear();
  fEnergyPdfFile5.clear();
  fEnergyCdfFile5.clear();
  fEneE.clear();
  fPdf.clear();
  fCdf.clear();
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergy::GetEnergy5(int ielemId, int mt, double energyK)
{
  fRnd  = new TRandom3(0);
  int i = -1;
  // std::cout<<"mt "<< mt <<"  "<< fMt5Values[ielemId].size() <<"  "<<energyK << std::endl;
  for (unsigned int l = 0; l < fMt5Values[ielemId].size(); l++) {
    if (fMt5Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  // std::cout<<"i "<< i <<"  "<< fMt5Values[ielemId][i] <<"  "<<fEnergy5OfMts[ielemId][i].size() << std::endl;
  if (i < 0) return 99;
  int min = 0;
  int max = fEnergy5OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= fEnergy5OfMts[ielemId][i][min])
    min = 0;
  else if (energyK >= fEnergy5OfMts[ielemId][i][max])
    min = max;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEnergy5OfMts[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  if (min < 0) min = 0;
  // for(unsigned int kk = 0; kk < fEnergy5OfMts[ielemId][i].size(); kk++)
  // std::cout << min <<"  "<< fEnergy5OfMts[ielemId][i][kk] <<"  "<< fEnergy5OfMts[ielemId][i].size() << std::endl;
  double fraction = (energyK - fEnergy5OfMts[ielemId][i][min]) /
                    (fEnergy5OfMts[ielemId][i][min + 1] - fEnergy5OfMts[ielemId][i][min]);
  double rnd1              = fRnd->Uniform(1);
  double rnd2              = fRnd->Uniform(1);
  if (rnd2 < fraction) min = min + 1;
  // std::cout<<"min "<< min <<"  "<< fEnergy5OfMts[ielemId][i][min] << std::endl;
  int k = 0;
  //   for (unsigned int j = 0; j < fEnergyPdf5OfMts[ielemId][i][min].size(); j++) {
  //     std::cout << j <<"  "<< fEnergyPdf5OfMts[ielemId][i][min].size() <<"  "<< fEnergyCdf5OfMts[ielemId][i][min][j]
  //     <<"  "<< rnd1 << std::endl;
  //   }
  int size = fEnergyCdf5OfMts[ielemId][i][min].size();
  for (unsigned int j = 0; j < fEnergyPdf5OfMts[ielemId][i][min].size(); j++) {
    //     std::cout << j <<"  "<< fEnergyPdf5OfMts[ielemId][i][min].size() <<"  "<<
    //     fEnergyCdf5OfMts[ielemId][i][min][j] <<"  "<< rnd1 << std::endl;
    if (rnd1 <= fEnergyCdf5OfMts[ielemId][i][min][j]) {
      k                    = j;
      if (k >= size - 1) k = size - 1;
      break;
    }
  }
  //    for (unsigned int j1 = 0; j1 < fEnergyOut5OfMts[ielemId][i][min].size(); j1++){
  //     std::cout<< fEnergyOut5OfMts[ielemId][i][min][j1] <<"  "<< fEnergyPdf5OfMts[ielemId][i][min][j1] <<std::endl;
  //   }
  double plk = (fEnergyPdf5OfMts[ielemId][i][min][k] - fEnergyPdf5OfMts[ielemId][i][min][k - 1]) /
               (fEnergyOut5OfMts[ielemId][i][min][k] - fEnergyOut5OfMts[ielemId][i][min][k - 1]);
  double plk2 = fEnergyPdf5OfMts[ielemId][i][min][k - 1] * fEnergyPdf5OfMts[ielemId][i][min][k - 1];

  double edes = 0;
  if (plk != 0)
    edes = fEnergyOut5OfMts[ielemId][i][min][k - 1] +
           (sqrt(plk2 + 2 * plk * (rnd1 - fEnergyCdf5OfMts[ielemId][i][min][k - 1])) -
            fEnergyPdf5OfMts[ielemId][i][min][k - 1]) /
               plk;
  double emin = fEnergyOut5OfMts[ielemId][i][min][1] +
                fraction * (fEnergyOut5OfMts[ielemId][i][min + 1][1] - fEnergyOut5OfMts[ielemId][i][min][1]);
  double emax = fEnergyOut5OfMts[ielemId][i][min][size - 1] +
                fraction * (fEnergyOut5OfMts[ielemId][i][min + 1][fEnergyCdf5OfMts[ielemId][i][min + 1].size() - 1] -
                            fEnergyOut5OfMts[ielemId][i][min][size - 1]);
  edes = fEnergyOut5OfMts[ielemId][i][min][1] +
         (edes - emin) * (fEnergyOut5OfMts[ielemId][i][min][size - 1] - fEnergyOut5OfMts[ielemId][i][min][1]) /
             (emax - emin);
  //	       std::cout<<plk <<"  "<< edes << std::endl;
  return edes;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergy::GetDelayedFraction(int ielemId, int mt, double energyK)
{
  fRnd  = new TRandom3(0);
  int i = -1;
  // std::cout<<"mt "<< mt <<"  "<< fMt5Values[ielemId].size() <<"  "<<energyK << std::endl;
  for (unsigned int l = 0; l < fMt5Values[ielemId].size(); l++) {
    if (fMt5Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  // std::cout<<"i "<< i <<"  "<< fMt5Values[ielemId][i] << std::endl;
  if (i < 0) return 99;
  int min = 0;
  int max = fEnergy5OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= fEnergy5OfMts[ielemId][i][min])
    min = 0;
  else if (energyK >= fEnergy5OfMts[ielemId][i][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEnergy5OfMts[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  return fFraction5OfMts[ielemId][i][min];
}
