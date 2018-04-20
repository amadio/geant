// 	This class is reconstructing probability tables for Energy distribution
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 22, 2016

#include "TList.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfSec.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfEnergyAng.h"
#include "Math/SpecFuncMathMore.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfEnergyAng)
#include "TRandom3.h"
#endif

    TNudyEndfEnergyAng::TNudyEndfEnergyAng()
{
}

//______________________________________________________________________________
TNudyEndfEnergyAng::TNudyEndfEnergyAng(TNudyEndfFile *file, double iQValue[])
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  fRnd = new TRandom3(0);
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    int NK = sec->GetN1();
    TIter recIter(sec->GetRecords());
    for (int k = 0; k < NK; k++) {
      // std::cout<<k<<" NK "<<NK<< std::endl;
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      div_t divr;
      AWRI   = sec->GetC2();
      int MT = sec->GetMT();
      //      int LCT = sec->GetL2();
      divr       = div(sec->GetC1(), 1000);
      double ZA  = divr.quot;
      double AA  = divr.rem;
      double NA1 = AA - ZA;
      int ZAP    = tab1->GetC1();
      double AWP = tab1->GetC2();
      // std::cout<<"ZAP = \t"<< ZAP <<" MT "<< MT <<  std::endl;
      // int LIP =tab1->GetL1();
      int LAW = tab1->GetL2();
      if (LAW == 3 || LAW == 4 || LAW == 0) continue;
      fLaw.push_back(LAW);
      fNR           = tab1->GetN1();
      fNP           = tab1->GetN2();
      divr          = div(ZAP, 1000);
      int particleA = divr.rem;
      int particleZ = divr.quot;
      fZd.push_back(particleZ);
      fAd.push_back(particleA);
      fMtNumbers6.push_back(MT);
      if (particleZ == 0 && particleA == 1) {
        fMtNumNeutron.push_back(MT);
      } else if (particleZ == 0 && particleA == 0) {
        fMtNumPhoton.push_back(MT);
      } else if (particleZ != 0 && particleA >= 1) {
        fMtNumCharge.push_back(MT);
      }
      //  std::cout<<"LAW = "<< LAW <<" MT "<<MT <<" fNP "<< fNP <<" ZAP = "<< ZAP <<" AWP "<<AWP <<" parZ "<< particleZ
      //   <<" parA "<< particleA << std::endl;
      if (LAW == 2) {
        fMtNumbers4.push_back(MT);
        int LANG, NL;
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
        fNp2                = tab2->GetN2();
        // std::cout<<"LANG "<< LANG <<" LEP "<< LEP <<" fNE2 "<< fNE2 << std::endl;
        for (int lis = 0; lis < fNp2; lis++) {
          TNudyEndfList *header = (TNudyEndfList *)recIter.Next();
          fEin.push_back(header->GetC2());
          LANG = header->GetL1();
          // NW   = header->GetN1();
          NL  = header->GetN2();
          fNP = NL;
          // std::cout<<"energy " <<fEin[lis] <<" LANG "<< LANG << std::endl;
          if (LANG == 0) {
            // std::cout<<"energy "<< fEin[lis] << std::endl;
            for (int j = 0; j < NL; j++) {
              fLegendCoef1.push_back(header->GetLIST(j));
            }
            fLegendCoef.push_back(fLegendCoef1);

            int k1     = 0;
            double fme = 0.0;
            do {
              fme      = 0.5;
              double x = -1. + k1 * 0.02;
              for (unsigned long j = 0; j < fLegendCoef[lis].size(); j++) {
                double leg = ROOT::Math::legendre(j + 1, x);
                fme += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[lis][j] * leg;
              }
              if (fme > 0.0) {
                fCosFile4.push_back(x);
                fCosPdfFile4.push_back(fme);
              }
              // std::cout << x << "   "<< fme << std::endl;
              k1++;
            } while (k1 < 101);

            // for(int l = 0; l < 100; l++){
            // RecursionLinearFile4(lis, fCosFile4[l], fCosFile4[l+1], fCosPdfFile4[l], fCosPdfFile4[l+1]);
            //}
            TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
            TNudyCore::Instance()->cdfGenerateT(fCosFile4, fCosPdfFile4, fCosCdfFile4);
            for (unsigned long i = 0; i < fCosFile4.size(); i++) {
              fCosc.push_back(fCosFile4[i]);
              fPdfc.push_back(fCosPdfFile4[i]);
              fCdfc.push_back(fCosCdfFile4[i]);
              // std::cout << fCosFile4[i] << "  "<< fCosPdfFile4[i] <<"  "<< fCosCdfFile4[i] << std::endl;
            }
            fCos2dc.push_back(fCosc);
            fPdf2dc.push_back(fPdfc);
            fCdf2dc.push_back(fCdfc);
            fCosFile4.clear();
            fCosPdfFile4.clear();
            fCosCdfFile4.clear();
            fCosc.clear();
            fPdfc.clear();
            fCdfc.clear();
            fLegendCoef1.clear();
          } else if (LANG > 0) {
            for (int i = 0; i < NL; i++) {
              fCosFile4.push_back(header->GetLIST(2 * i + 0));
              fCosPdfFile4.push_back(header->GetLIST(2 * i + 1));
            }
            int size = fCosFile4.size();
            for (int l = 0; l < size - 1; l++) {
              RecursionLinearProb(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
            }
            TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
            TNudyCore::Instance()->cdfGenerateT(fCosFile4, fCosPdfFile4, fCosCdfFile4);
            for (unsigned long i = 0; i < fCosFile4.size(); i++) {
              fCosc.push_back(fCosFile4[i]);
              fPdfc.push_back(fCosPdfFile4[i]);
              fCdfc.push_back(fCosCdfFile4[i]);
            }
            fCos2dc.push_back(fCosc);
            fPdf2dc.push_back(fPdfc);
            fCdf2dc.push_back(fCdfc);
            fCosFile4.clear();
            fCosPdfFile4.clear();
            fCosCdfFile4.clear();
            fCosc.clear();
            fPdfc.clear();
            fCdfc.clear();
          }
        }
        fEin2dc.push_back(fEin);
        fCos3dc.push_back(fCos2dc);
        fPdf3dc.push_back(fPdf2dc);
        fCdf3dc.push_back(fCdf2dc);
        fEin.clear();
        fCos2dc.clear();
        fPdf2dc.clear();
        fCdf2dc.clear();
        fLegendCoef.clear();
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fP1.clear();
        continue;
      }
      fMtNumbers.push_back(MT);
      // std::cout<<"ZAP = "<< ZAP <<" AWP "<<AWP <<" parZ "<< particleZ <<" parA "<< particleA << std::endl;
      //---------------------------------------------
      switch (LAW) {
      case 1: {
        int NA, NEP, LANG;
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        fMtLct.push_back(1);
        for (int crs = 0; crs < fNP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
        }
        TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
        LANG                = tab2->GetL1();
        // LEP  = tab2->GetL2();
        fNr2 = tab2->GetN1();
        fNp2 = tab2->GetN2();
        // std::cout<<"LANG "<< LANG <<"  "<< NA<<"  "<< MT<< std::endl;
        for (int lis = 0; lis < fNp2; lis++) {
          TNudyEndfList *header = (TNudyEndfList *)recIter.Next();
          fEin.push_back(header->GetC2());
          double sumein = 0;
          // ND   = header->GetL1();
          NA  = header->GetL2();
          NEP = header->GetN2();
          // std::cout <<"fEin  "<<header->GetC2() <<" NA "<< NA <<" lang "<< LANG << std::endl;
          if (LANG == 2 && NA < 2) {
            for (int lis1 = 0; lis1 < NEP; lis1++) {
              fEdes6.push_back(header->GetLIST(lis1 * 3 + 0));
              fF06.push_back(header->GetLIST(lis1 * 3 + 1));
              fR6.push_back(header->GetLIST(lis1 * 3 + 2));
              // std::cout<<lis1<<"  "<<fEdes6[lis1] << std::endl;
            }
            double epsa = header->GetC2() * AWRI / (1. + AWRI);
            double AC = 1. + AA, ZC = ZA, NC = AC - ZC;
            double AB = AC - particleA, ZB = ZA - particleZ, NB = AB - ZB;
            double Ia = 0.0;
            double Ib = 0.0;
            double Sa = 15.68 * (AC - AA) - 28.07 * ((NC - ZC) * (NC - ZC) / AC - (NA1 - ZA) * (NA1 - ZA) / AA) -
                        -18.56 * (pow(AC, 2. / 3.) - pow(AA, 2. / 3.)) +
                        33.22 * ((NC - ZC) * (NC - ZC) / pow(AC, 4. / 3) - (NA1 - ZA) * (NA1 - ZA) / pow(AA, 4. / 3)) -
                        0.717 * (ZC * ZC / pow(AC, 1. / 3.) - ZA * ZA / pow(AA, 1. / 3.)) +
                        1.211 * (ZC * ZC / AC - ZA * ZA / AA) - Ia;
            double Sb = 15.68 * (AC - AB) - 28.07 * ((NC - ZC) * (NC - ZC) / AC - (NB - ZB) * (NB - ZB) / AB) -
                        -18.56 * (pow(AC, 2. / 3.) - pow(AB, 2. / 3.)) +
                        33.22 * ((NC - ZC) * (NC - ZC) / pow(AC, 4. / 3) - (NB - ZB) * (NB - ZB) / pow(AB, 4. / 3)) -
                        0.717 * (ZC * ZC / pow(AC, 1. / 3.) - ZB * ZB / pow(AB, 1. / 3.)) +
                        1.211 * (ZC * ZC / AC - ZB * ZB / AB) - Ib;
            double C1 = 0.04, C2 = 1.8E-6, C3 = 6.7E-7, Et1 = 130.0, Et3 = 41.0;
            double Mn1 = 1;   //, Mp = 1, Md = 1,Mt = 1, Mhe3 = 1, Malpha = 0;
            double mn  = 0.5; //, mp = 1, md = 1, mt = 1, mhe3 = 1, malpha = 2;

            int k1 = 0;
            do {
              double x = -1. + k1 * 0.02;
              // std::cout<<" k1 "<< k1 <<" cos "<< y << std::endl;
              double sumprob = 0;
              for (unsigned long j = 0; j < fEdes6.size(); j++) {
                // std::cout<<"j "<< j <<" cosine "<< x <<" eoutcm "<< fEdes6[j] << std::endl;
                double epsb = fEdes6[j] * (AWP + AWRI) / (AWRI + 1. - AWP);
                double ea = epsa + Sa, eb = epsb + Sb;
                double R1 = std::min(ea, Et1), R3 = std::min(eb, Et3);
                double X1 = R1 * eb / ea, X3 = R3 * eb / ea;
                double a0   = C1 * X1 + C2 * X1 * X1 * X1 + C3 * Mn1 * mn * pow(X3, 4);
                double prob = (0.5 * a0 * fF06[j] / sinh(a0)) * (cosh(a0 * x) + fR6[j] * sinh(a0 * x));
                // std::cout<<"flaglab "<< flaglab <<"  "<< xlab <<"   "<< fEnergyFile5.size() << std::endl;
                // if(prob > 1E-15){
                fEnergyFile5.push_back(fEdes6[j]);
                fEnergyPdfFile5.push_back(prob);
                sumprob += prob;
                // std::cout<< eoutlab<<"  "<< prob <<"  "<<xlab<<"  "<< y <<  std::endl;
                //}
              }
              // std::cout<<" energies "<< fEnergyFile5.size() << std::endl;
              // for(unsigned long i = 0; i < fEnergyFile5.size(); i++){
              // std::cout << fEnergyFile5[i] << "  "<< fEnergyPdfFile5[i] <<" k1 "<< k1 << std::endl;
              //}
              double fsum = 0.0;
              for (unsigned int cr = 1; cr < fEnergyFile5.size(); cr++) {
                fsum +=
                    0.5 * (fEnergyPdfFile5[cr] + fEnergyPdfFile5[cr - 1]) * (fEnergyFile5[cr] - fEnergyFile5[cr - 1]);
              }
              // if(fsum > 1E-20){
              sumein += fsum;
              fCosIn.push_back(x);
              fCosInPdf.push_back(fsum);
              // std::cout<<"fsum "<< fsum <<" sumprob "<< sumprob << std::endl;
              TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
              TNudyCore::Instance()->ThinningDuplicate(fEnergyFile5, fEnergyPdfFile5);
              TNudyCore::Instance()->cdfGenerateT(fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5);
              for (unsigned long i = 0; i < fEnergyFile5.size(); i++) {
                fEoute.push_back(fEnergyFile5[i]);
                fPdfe.push_back(fEnergyPdfFile5[i]);
                fCdfe.push_back(fEnergyCdfFile5[i]);
                // std::cout <<"energy pdf "<< fEnergyFile5[i] << "  "<< fEnergyPdfFile5[i] <<"  "<< fEnergyCdfFile5[i]
                // << std::endl;
              }
              fEout2de.push_back(fEoute);
              fPdf2de.push_back(fPdfe);
              fCdf2de.push_back(fCdfe);
              //}
              fEnergyFile5.clear();
              fEnergyPdfFile5.clear();
              fEnergyCdfFile5.clear();
              fEoute.clear();
              fPdfe.clear();
              fCdfe.clear();
              k1++;
            } while (k1 < 101);
            TNudyCore::Instance()->cdfGenerateT(fCosIn, fCosInPdf, fCosInCdf);
            // std::cout<<"cos size "<<  fCosIn.size() <<std::endl;
            for (unsigned long i = 0; i < fCosIn.size(); i++) {
              // std::cout << " cospdf "<< fCosIn[i] <<"  "<< fCosInPdf[i] <<"  "<< fCosInCdf[i] << std::endl;
            }
            // if(sumein > 1E-50){
            fCos2d.push_back(fCosIn);
            fCosinpdf2d.push_back(fCosInPdf);
            fCosincdf2d.push_back(fCosInCdf);
            fEout3de.push_back(fEout2de);
            fPdf3de.push_back(fPdf2de);
            fCdf3de.push_back(fCdf2de);
            fEout2de.clear();
            fPdf2de.clear();
            fCdf2de.clear();
            fCosIn.clear();
            fCosInPdf.clear();
            fCosInCdf.clear();
            //}
            // if(sumein < 1E-50 && fEin.size()>0)fEin.erase(fEin.begin()+ fEin.size()-1);
            // std::cout << "fEin end "<< fEin.size() << std::endl;
            fEdes6.clear();
            fF06.clear();
            fR6.clear();
          } else if (LANG == 2 && NA == 2) {
            //	    std::cout<<"NA = "<< NA <<std::endl;
            for (int lis1 = 0; lis1 < NEP; lis1++) {
              fEdes6.push_back(header->GetLIST(lis1 * 4 + 0));
              fF06.push_back(header->GetLIST(lis1 * 4 + 1));
              fR6.push_back(header->GetLIST(lis1 * 4 + 2));
              fA6.push_back(header->GetLIST(lis1 * 4 + 3));
              //	  std::cout<<lis1<<"  "<<fEdes6[lis1]<<"  "<<fF06[lis1] <<"  "<<fR6[lis1]<<"  "<<fA6[lis1]<<
              // std::endl;
            }
            fEdes6.clear();
            fF06.clear();
            fR6.clear();
            fA6.clear();
          } else if (LANG == 1) {
            for (int lis1 = 0; lis1 < NEP; lis1++) {
              fEdes6.push_back(header->GetLIST(lis1 * (NA + 2) + 0));
              for (int j = 0; j < NA + 1; j++) {
                fLegendCoef1.push_back(header->GetLIST(lis1 * (NA + 2) + j + 1));
                // std::cout <<  header->GetLIST(lis1*(NA+2) + j + 1) << std::endl;
              }
              fLegendCoef.push_back(fLegendCoef1);
              fLegendCoef1.clear();
            }
            int k1     = 0;
            double fme = 0;
            do {
              double sumprob = 0;
              double x       = -1. + k1 * 0.01;
              for (unsigned long m = 0; m < fEdes6.size(); m++) {
                fme = 0.5;
                for (unsigned long j = 0; j < fLegendCoef[m].size(); j++) {
                  double leg = ROOT::Math::legendre(j + 1, x);
                  fme += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[m][j] * leg;
                  // std::cout << x <<"  "<< leg<<" leg "<<fLegendCoef[m][j] <<"  "<<fLegendCoef[m].size() <<"  "<< fme
                  // << std::endl;
                }
                if (fme > 1E-15) {
                  fEnergyFile5.push_back(fEdes6[m]);
                  fEnergyPdfFile5.push_back(fme);
                  sumprob += fme;
                  // std::cout<<"eout "<< fEdes6[m] <<"  "<< fme << std::endl;
                }
              }
              double fsum = 0.0;
              for (unsigned int cr = 1; cr < fEnergyFile5.size(); cr++) {
                fsum +=
                    0.5 * (fEnergyPdfFile5[cr] + fEnergyPdfFile5[cr - 1]) * (fEnergyFile5[cr] - fEnergyFile5[cr - 1]);
              }
              // if(fsum > 1E-20){
              sumein += fsum;
              fCosIn.push_back(x);
              fCosInPdf.push_back(fsum);
              TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
              TNudyCore::Instance()->ThinningDuplicate(fEnergyFile5, fEnergyPdfFile5);
              TNudyCore::Instance()->cdfGenerateT(fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5);
              for (unsigned long i = 0; i < fEnergyFile5.size(); i++) {
                fEoute.push_back(fEnergyFile5[i]);
                fPdfe.push_back(fEnergyPdfFile5[i]);
                fCdfe.push_back(fEnergyCdfFile5[i]);
                // std::cout <<sumprob<<" energy "<< fEnergyFile5[i] << "  "<< fEnergyPdfFile5[i] <<"  "<<
                // fEnergyCdfFile5[i] << std::endl;
              }
              fEout2de.push_back(fEoute);
              fPdf2de.push_back(fPdfe);
              fCdf2de.push_back(fCdfe);
              //}
              fEnergyFile5.clear();
              fEnergyPdfFile5.clear();
              fEnergyCdfFile5.clear();
              fEoute.clear();
              fPdfe.clear();
              fCdfe.clear();
              k1++;
            } while (k1 < 201);
            // if(sumein > 1E-50){
            TNudyCore::Instance()->cdfGenerateT(fCosIn, fCosInPdf, fCosInCdf);
            for (unsigned long i = 0; i < fCosIn.size(); i++) {
              // std::cout << fCosIn[i] << " cospdf "<< fCosInPdf[i] <<"  "<< fCosInCdf[i] << std::endl;
            }
            fCos2d.push_back(fCosIn);
            fCosinpdf2d.push_back(fCosInPdf);
            fCosincdf2d.push_back(fCosInCdf);
            fEout3de.push_back(fEout2de);
            fPdf3de.push_back(fPdf2de);
            fCdf3de.push_back(fCdf2de);
            fEout2de.clear();
            fPdf2de.clear();
            fCdf2de.clear();
            fCosIn.clear();
            fCosInPdf.clear();
            fCosInCdf.clear();
            //}
            fEdes6.clear();
            fLegendCoef.clear();
            // if(sumein < 1E-50 && fEin.size()>0)fEin.erase(fEin.begin()+ fEin.size()-1);
            // std::cout << "fEin "<< fEin.size() << std::endl;
          }
        }
        fCos3d.push_back(fCos2d);
        fCosinpdf3d.push_back(fCosinpdf2d);
        fCosincdf3d.push_back(fCosincdf2d);
        fEin2d.push_back(fEin);
        fEout4de.push_back(fEout3de);
        fPdf4de.push_back(fPdf3de);
        fCdf4de.push_back(fCdf3de);
        fEin.clear();
        fEout3de.clear();
        fPdf3de.clear();
        fCdf3de.clear();
        fCos2d.clear();
        fCosinpdf2d.clear();
        fCosincdf2d.clear();
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fP1.clear();
        //---------------------------------------------
      } break;
      case 2: {
        fMtLct.push_back(2);
        int LANG, NL;
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        TNudyEndfTab2 *tab2   = (TNudyEndfTab2 *)recIter.Next();
        fNp2                  = tab2->GetN2();
        TNudyEndfList *header = (TNudyEndfList *)recIter.Next();
        fEin.push_back(header->GetC2());
        LANG = header->GetL1();
        // std::cout<<"LANG "<< LANG << std::endl;
        for (int lis = 0; lis < fNp2; lis++) {
          // NW   = header->GetN1();
          NL  = header->GetN2();
          fNP = NL;
          // std::cout<<"energy " <<fEin[lis] <<" LANG "<< LANG << std::endl;
          if (LANG == 0) {
            // std::cout<<"energy "<< fEin[lis] << std::endl;
            for (int j = 0; j < NL; j++) {
              fLegendCoef1.push_back(header->GetLIST(j));
            }
            fLegendCoef.push_back(fLegendCoef1);

            int k1     = 0;
            double fme = 0.0;
            do {
              fme      = 0.5;
              double x = -1. + k1 * 0.02;
              for (unsigned long j = 0; j < fLegendCoef[lis].size(); j++) {
                double leg = ROOT::Math::legendre(j + 1, x);
                fme += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[lis][j] * leg;
              }
              fCosFile4.push_back(x);
              fCosPdfFile4.push_back(fme);
              // std::cout<<"cos "<< x <<"  "<< fme << std::endl;
              k1++;
            } while (k1 < 101);
            // for(int l = 0; l < 100; l++){
            // RecursionLinearFile4(lis, fCosFile4[l], fCosFile4[l+1], fCosPdfFile4[l], fCosPdfFile4[l+1]);
            //}
            TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
            TNudyCore::Instance()->cdfGenerateT(fCosFile4, fCosPdfFile4, fCosCdfFile4);
            for (unsigned long i = 0; i < fCosFile4.size(); i++) {
              fCosc.push_back(fCosFile4[i]);
              fPdfc.push_back(fCosPdfFile4[i]);
              fCdfc.push_back(fCosCdfFile4[i]);
              // std::cout << fCosFile4[i] << "  "<< fCosPdfFile4[i] <<"  "<< fCosCdfFile4[i] << std::endl;
            }
            fCos2dc.push_back(fCosc);
            fPdf2dc.push_back(fPdfc);
            fCdf2dc.push_back(fCdfc);
            fCosFile4.clear();
            fCosPdfFile4.clear();
            fCosCdfFile4.clear();
            fCosc.clear();
            fPdfc.clear();
            fCdfc.clear();
            fLegendCoef1.clear();
          } else if (LANG > 0) {
            for (int i = 0; i < NL; i++) {
              fCosFile4.push_back(header->GetLIST(2 * i + 0));
              fCosPdfFile4.push_back(header->GetLIST(2 * i + 1));
            }
            int size = fCosFile4.size();
            for (int l = 0; l < size - 1; l++) {
              RecursionLinearProb(fCosFile4[l], fCosFile4[l + 1], fCosPdfFile4[l], fCosPdfFile4[l + 1]);
            }
            TNudyCore::Instance()->Sort(fCosFile4, fCosPdfFile4);
            TNudyCore::Instance()->cdfGenerateT(fCosFile4, fCosPdfFile4, fCosCdfFile4);
            for (unsigned long i = 0; i < fCosFile4.size(); i++) {
              fCosc.push_back(fCosFile4[i]);
              fPdfc.push_back(fCosPdfFile4[i]);
              fCdfc.push_back(fCosCdfFile4[i]);
            }
            fCos2dc.push_back(fCosc);
            fPdf2dc.push_back(fPdfc);
            fCdf2dc.push_back(fCdfc);
            fCosFile4.clear();
            fCosPdfFile4.clear();
            fCosCdfFile4.clear();
            fCosc.clear();
            fPdfc.clear();
            fCdfc.clear();
          }
        }
        fEin2d.push_back(fEin);
        fCos3dc.push_back(fCos2dc);
        fPdf3dc.push_back(fPdf2dc);
        fCdf3dc.push_back(fCdf2dc);
        fEin.clear();
        fCos2dc.clear();
        fPdf2dc.clear();
        fCdf2dc.clear();
        fLegendCoef.clear();
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fP1.clear();
        //---------------------------------------------
      } break;
      case 3:
      case 4: {
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
          fEin.push_back(tab1->GetX(crs));
          // std::cout<<"energy1 "<< fE1[crs] <<" multip1 "<< fP1[crs] << std::endl;
        }
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fP1.clear();
        fEin2d.push_back(fEin);
        fEin.clear();
        //---------------------------------------------
      } break;
      case 5: {

        //---------------------------------------------
      } break;
      case 6: {
        fMtLct.push_back(2);
        int NPSX;
        double APSX;
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
        }
        TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
        APSX                  = header->GetC1();
        NPSX                  = header->GetN2();
        double energy = std::fabs(-(1. + AWRI) * iQValue[sec->GetMT()] / AWRI), eout = 1E-5;
        energy += 0.001 * (energy);
        do {
          fEin.push_back(energy);
          double Ea    = AWRI * energy / (1. + AWRI) + iQValue[sec->GetMT()];
          double EiMax = (APSX - AWP) * Ea / APSX;
          // std::cout<<"Incident energy "<< energy <<"  "<< iQValue[sec->GetMT()] <<"  "<< Ea << std::endl;
          // double EStar = energy/(AWP + AWRI);
          double C[3] = {4. / (PI * EiMax * EiMax), 105. / (32. * pow(EiMax, 3.5)), 256. / (14. * PI * pow(EiMax, 5))};
          // std::cout << energy <<"  "<< fE1[fNP-1] <<"  "<< EiMax <<"  \t"<< AWP << std::endl;
          // std::cout<<"c[0] "<< C[0] <<" c[1] "<< C[1] <<" c[2] "<< C[2] << std::endl;
          // std::cout<<"APSX "<< APSX <<" NPSX "<< NPSX  <<" qval "<< iQValue[sec->GetMT()] << std::endl;
          int k1     = 0;
          double rn1 = 0, rn2 = 0, rn3 = 0, rn4 = 0, rn5 = 0, rn6 = 0, rn7 = 0, rn8 = 0, rn9 = 0;
          double p = 0, xp = 0, yp = 0, tp = 0, rn12 = 0, rn34 = 0;
          // sampling of the outgoing energy is sampled based on R28 of LA-9721-MS.
          do {
            double cosCM = -1 + k1 * 0.02;
            fCosIn.push_back(cosCM);
            int iei = 0;
            do {
              do {
                rn1  = fRnd->Uniform(1);
                rn2  = fRnd->Uniform(1);
                rn12 = rn1 * rn1 + rn2 * rn2;
              } while (rn12 > 1);
              do {
                rn3  = fRnd->Uniform(1);
                rn4  = fRnd->Uniform(1);
                rn34 = rn3 * rn3 + rn4 * rn4;
              } while (rn34 > 1);
              rn5              = fRnd->Uniform(1);
              rn6              = fRnd->Uniform(1);
              rn7              = fRnd->Uniform(1);
              rn8              = fRnd->Uniform(1);
              if (NPSX == 3) p = rn5;
              if (NPSX == 4) p = rn5 * rn6;
              if (NPSX == 5) p = rn5 * rn6 * rn7 * rn8;
              rn9              = fRnd->Uniform(1);
              xp               = -rn1 * log(rn12) / rn12 - log(rn9);
              yp               = -rn3 * log(rn34) / rn34 - log(p);
              tp               = xp / (xp + yp);
              eout             = tp * EiMax;
              double ppeCm     = 0.0;
              ppeCm            = 2 * C[NPSX - 3] * sqrt(eout) * pow((EiMax - eout), 1.5 * NPSX - 4);
              fEnergyFile5.push_back(eout);
              fEnergyPdfFile5.push_back(ppeCm);
              iei++;
            } while (iei < 500);
            TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
            fCosInPdf.push_back(0.5);
            TNudyCore::Instance()->cdfGenerateT(fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5);
            for (unsigned long i = 0; i < fEnergyFile5.size(); i++) {
              fEoute.push_back(fEnergyFile5[i]);
              fPdfe.push_back(fEnergyPdfFile5[i]);
              fCdfe.push_back(fEnergyCdfFile5[i]);
              // std::cout<< fEnergyFile5[i] << "  "<< fEnergyPdfFile5[i] << "  "<< fEnergyCdfFile5[i] <<"  "<< cosCM <<
              // std::endl;
            }
            fEout2de.push_back(fEoute);
            fPdf2de.push_back(fPdfe);
            fCdf2de.push_back(fCdfe);
            fEnergyFile5.clear();
            fEnergyPdfFile5.clear();
            fEnergyCdfFile5.clear();
            fEoute.clear();
            fPdfe.clear();
            fCdfe.clear();
            k1++;
          } while (k1 < 101);
          TNudyCore::Instance()->cdfGenerateT(fCosIn, fCosInPdf, fCosInCdf);
          fCos2d.push_back(fCosIn);
          fCosinpdf2d.push_back(fCosInPdf);
          fCosincdf2d.push_back(fCosInCdf);
          fEout3de.push_back(fEout2de);
          fPdf3de.push_back(fPdf2de);
          fCdf3de.push_back(fCdf2de);
          fEout2de.clear();
          fPdf2de.clear();
          fCdf2de.clear();
          fCosIn.clear();
          fCosInPdf.clear();
          fCosInCdf.clear();
          energy *= 2;
        } while (energy < fE1[fNP - 1] + iQValue[sec->GetMT()]);
        fCos3d.push_back(fCos2d);
        fCosinpdf3d.push_back(fCosinpdf2d);
        fCosincdf3d.push_back(fCosincdf2d);
        fEin2d.push_back(fEin);
        fEout4de.push_back(fEout3de);
        fPdf4de.push_back(fPdf3de);
        fCdf4de.push_back(fCdf3de);
        fEin.clear();
        fEout3de.clear();
        fPdf3de.clear();
        fCdf3de.clear();
        fCos2d.clear();
        fCosinpdf2d.clear();
        fCosincdf2d.clear();
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fP1.clear();
        // angle is isotropic in the CM System for this law
        //---------------------------------------------
      } break;
      case 7: {
        fMtLct.push_back(1);
        for (int cr = 0; cr < fNR; cr++) {
          fNbt1.push_back(tab1->GetNBT(cr));
          fInt1.push_back(tab1->GetINT(cr));
        }
        for (int crs = 0; crs < fNP; crs++) {
          fE1.push_back(tab1->GetX(crs));
          fP1.push_back(tab1->GetY(crs));
        }

        TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
        fNR2                = tab2->GetN1();
        fNE2                = tab2->GetN2();
        for (int cr = 0; cr < fNR2; cr++) {
          fNbt2.push_back(tab2->GetNBT(cr));
          fInt2.push_back(tab2->GetINT(cr));
        }
        for (int cr1 = 0; cr1 < fNE2; cr1++) {
          TNudyEndfTab2 *tab3 = (TNudyEndfTab2 *)recIter.Next();
          fEin.push_back(tab3->GetC2());
          // std::cout <<"energy "<< tab3->GetC2() << std::endl;
          int NMU = tab3->GetN2();
          for (int cr2 = 0; cr2 < NMU; cr2++) {
            TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
            fCosIn.push_back(tab12->GetC2());
            // std::cout<<"cosine "<< fCosIn[cr2] <<"  "<< fEin[cr1] << std::endl;
            fNr3 = tab12->GetN1();
            fNp3 = tab12->GetN2();
            for (int cr = 0; cr < fNr3; cr++) {
              fNbt3.push_back(tab12->GetNBT(cr));
              fInt3.push_back(tab12->GetINT(cr));
            }
            for (int crs = 0; crs < fNp3; crs++) {
              fEnergyFile5.push_back(tab12->GetX(crs));
              fEnergyPdfFile5.push_back(tab12->GetY(crs));
            }
            double fsum = 0.0;
            for (int cr = 1; cr < fNp3; cr++) {
              fsum += 0.5 * (fEnergyPdfFile5[cr] + fEnergyPdfFile5[cr - 1]) * (fEnergyFile5[cr] - fEnergyFile5[cr - 1]);
            }
            fCosInPdf.push_back(fsum);
            for (int cr = 0; cr < fNp3 - 1; cr++) {
              RecursionLinearFile5Prob(fEnergyFile5[cr], fEnergyFile5[cr + 1], fEnergyPdfFile5[cr],
                                       fEnergyPdfFile5[cr + 1]);
            }
            TNudyCore::Instance()->Sort(fEnergyFile5, fEnergyPdfFile5);
            TNudyCore::Instance()->ThinningDuplicate(fEnergyFile5, fEnergyPdfFile5);
            TNudyCore::Instance()->cdfGenerateT(fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5);
            for (unsigned long i = 0; i < fEnergyFile5.size(); i++) {
              fEoute.push_back(fEnergyFile5[i]);
              fPdfe.push_back(fEnergyPdfFile5[i]);
              fCdfe.push_back(fEnergyCdfFile5[i]);
              // std::cout << fEnergyFile5[i] << "  "<< fEnergyPdfFile5[i] <<"  "<< fEnergyCdfFile5[i] << std::endl;
            }
            fEout2de.push_back(fEoute);
            fPdf2de.push_back(fPdfe);
            fCdf2de.push_back(fCdfe);
            fEnergyFile5.clear();
            fEnergyPdfFile5.clear();
            fEnergyCdfFile5.clear();
            fEoute.clear();
            fPdfe.clear();
            fCdfe.clear();
            fNbt3.clear();
            fInt3.clear();
          }
          TNudyCore::Instance()->cdfGenerateT(fCosIn, fCosInPdf, fCosInCdf);
          // for(unsigned long i = 0; i < fCosIn.size(); i++){
          // std::cout << fCosIn[i] << "  "<< fCosInPdf[i] <<"  "<< fCosInCdf[i] << std::endl;
          //}
          fCos2d.push_back(fCosIn);
          fCosinpdf2d.push_back(fCosInPdf);
          fCosincdf2d.push_back(fCosInCdf);
          fEout3de.push_back(fEout2de);
          fPdf3de.push_back(fPdf2de);
          fCdf3de.push_back(fCdf2de);
          fEout2de.clear();
          fPdf2de.clear();
          fCdf2de.clear();
          fCosIn.clear();
          fCosInPdf.clear();
          fCosInCdf.clear();
        }
        fCos3d.push_back(fCos2d);
        fCosinpdf3d.push_back(fCosinpdf2d);
        fCosincdf3d.push_back(fCosincdf2d);
        fEin2d.push_back(fEin);
        fEout4de.push_back(fEout3de);
        fPdf4de.push_back(fPdf3de);
        fCdf4de.push_back(fCdf3de);
        fEin.clear();
        fEout3de.clear();
        fPdf3de.clear();
        fCdf3de.clear();
        fCos2d.clear();
        fCosinpdf2d.clear();
        fCosincdf2d.clear();
        fNbt1.clear();
        fInt1.clear();
        fE1.clear();
        fP1.clear();
        fNbt2.clear();
        fInt2.clear();
        fNbt3.clear();
        fInt3.clear();
      } break;
      }
    }
  }
  fCos6OfMts.push_back(fCos3d);
  fCosin6Pdf.push_back(fCosinpdf3d);
  fCosin6Cdf.push_back(fCosincdf3d);
  fEnergy5OfMts.push_back(fEin2d);
  fEnergyOut6OfMts.push_back(fEout4de);
  fEnergyPdf6OfMts.push_back(fPdf4de);
  fEnergyCdf6OfMts.push_back(fCdf4de);
  fMt6Values.push_back(fMtNumbers6);
  fMt5Values.push_back(fMtNumbers);
  fMt6Neutron.push_back(fMtNumNeutron);
  fMt4Lct.push_back(fMtLct);
  fLaw6.push_back(fLaw);
  fZD6.push_back(fZd);
  fAD6.push_back(fAd);
  fLaw.clear();
  fMtNumbers.clear();
  fMtNumbers6.clear();
  fZd.clear();
  fAd.clear();
  fMtLct.clear();
  fEin2d.clear();
  fEout4de.clear();
  fPdf4de.clear();
  fCdf4de.clear();
  fCos3d.clear();
  fMt4Values.push_back(fMtNumbers4);

  fEnergy4OfMts.push_back(fEin2dc);
  fCos4OfMts.push_back(fCos3dc);
  fCosPdf4OfMts.push_back(fPdf3dc);
  fCosCdf4OfMts.push_back(fCdf3dc);
  std::vector<int>().swap(fMtNumbers4);
  fEin2dc.clear();
  fCos3dc.clear();
  fPdf3dc.clear();
  fCdf3dc.clear();

  //  std::cout << fMt4Values[0].size() << std::endl;
  //  std::cout << fMt5Values[0].size() << std::endl;
  //  std::cout << fMt6Values[0].size() << std::endl;
  //  std::cout << fLaw6[0].size() << std::endl;
  //  for(unsigned int i = 0; i< fMt6Values[0].size(); i++ ){
  //  std::cout<<"MT = "<< fMt6Values[0][i]<<"  "<<fMt6Values[0].size()<<"  "<< fLaw6[0][i] <<"  "<< fZD6[0][i] <<"  "<<
  //  fAD6[0][i] << std::endl;
  //  }
  //     for(unsigned i = 0; i < fEnergy5OfMts[0].size(); i++)
  //       std::cout<<fEnergy5OfMts[0][i].size()<< std::endl;
  //  }
  // for(unsigned int i = 0; i< fMt4Values[0].size(); i++ ){
  // std::cout<<"MT = "<< fMt4Values[0][i]<<"  "<<fMt4Values[0].size()<< std::endl;
  //   for(unsigned j = 0; j < fEnergy4OfMts[0][i].size(); j++)
  //     std::cout<<fEnergy4OfMts[0][i][j]<<"  "<<fEnergy4OfMts[0][i].size()<< std::endl;
  //}
}
//------------------------------------------------------------------------------------------------------

TNudyEndfEnergyAng::~TNudyEndfEnergyAng()
{
  fCosFile4.shrink_to_fit();
  fEnergyFile5.shrink_to_fit();
  fCosPdfFile4.shrink_to_fit();
  fEnergyPdfFile5.shrink_to_fit();
  fCosCdfFile4.shrink_to_fit();
  fEnergyCdfFile5.shrink_to_fit();
  fEdes6.shrink_to_fit();
  fF06.shrink_to_fit();
  fR6.shrink_to_fit();
  fA6.shrink_to_fit();
  fE1.shrink_to_fit();
  fP1.shrink_to_fit();
  fE2.shrink_to_fit();
  fP2.shrink_to_fit();
  fE3.shrink_to_fit();
  fP3.shrink_to_fit();
  INorm.shrink_to_fit();
  fLaw.shrink_to_fit();
  fMtNumbers.shrink_to_fit();
  fZd.shrink_to_fit();
  fAd.shrink_to_fit();
  fMtLct.shrink_to_fit();
  fNbt1.shrink_to_fit();
  fInt1.shrink_to_fit();
  fNbt2.shrink_to_fit();
  fInt2.shrink_to_fit();
  fNbt3.shrink_to_fit();
  fInt3.shrink_to_fit();
  fEin.shrink_to_fit();
  fCosc.shrink_to_fit();
  fCdfc.shrink_to_fit();
  fPdfc.shrink_to_fit();
  fLegendCoef1.shrink_to_fit();
  fCosIn.shrink_to_fit();
  fCosInPdf.shrink_to_fit();
  fCosInCdf.shrink_to_fit();
  fEoute.shrink_to_fit();
  fCdfe.shrink_to_fit();
  fPdfe.shrink_to_fit();
  fCos2d.shrink_to_fit();
  fCosinpdf2d.shrink_to_fit();
  fCosincdf2d.shrink_to_fit();
  fCos2dc.shrink_to_fit();
  fPdf2dc.shrink_to_fit();
  fCdf2dc.shrink_to_fit();
  fLegendCoef.shrink_to_fit();
  fEin2d.shrink_to_fit();
  fCos3d.shrink_to_fit();
  fCosinpdf3d.shrink_to_fit();
  fCosincdf3d.shrink_to_fit();
  fCos3dc.shrink_to_fit();
  fPdf3dc.shrink_to_fit();
  fCdf3dc.shrink_to_fit();
  fEout2de.shrink_to_fit();
  fPdf2de.shrink_to_fit();
  fCdf2de.shrink_to_fit();
  fEout3de.shrink_to_fit();
  fPdf3de.shrink_to_fit();
  fCdf3de.shrink_to_fit();
  fEout4de.shrink_to_fit();
  fPdf4de.shrink_to_fit();
  fCdf4de.shrink_to_fit();
  fMt4Values.shrink_to_fit();
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::RecursionLinearFile4(int i, double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 0.5;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  //  std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< pdf1 <<"  "<< pdf2 << std::endl;
  for (unsigned long j = 0; j < fLegendCoef[i].size(); j++) {
    double leg = ROOT::Math::legendre(j + 1, mid);
    pdf += 0.5 * (2. * (j + 1) + 1.) * fLegendCoef[i][j] * leg;
  }
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  //  std::cout << mid <<" linear "<< pdf <<"  "<< pdfmid1 << std::endl;

  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  fCosFile4.push_back(mid);
  fCosPdfFile4.push_back(pdf);
  RecursionLinearFile4(i, x1, mid, pdf1, pdf);
  RecursionLinearFile4(i, mid, x2, pdf, pdf2);
  return 0;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::RecursionLinearProb(double x1, double x2, double pdf1, double pdf2)
{
  double pdf = 1.0;
  double mid = 0.5 * (x1 + x2);
  if ((pdf1 == 0.0 && pdf2 == 0.0) || x1 == x2) return 0;
  // std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< pdf1 <<"  "<< pdf2 << std::endl;
  pdf            = TNudyCore::Instance()->Interpolate(fNbt1, fInt1, fNR, fCosFile4, fCosPdfFile4, fNP, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1) * (mid - x1) / (x2 - x1);
  //  std::cout << mid <<" linear "<< pdf <<"  "<< pdfmid1 << std::endl;

  if (fabs((pdf - pdfmid1) / pdfmid1) <= 1E-3) {
    return 0;
  }
  fCosFile4.push_back(mid);
  fCosPdfFile4.push_back(pdf);
  RecursionLinearProb(x1, mid, pdf1, pdf);
  RecursionLinearProb(mid, x2, pdf, pdf2);
  return 0;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::RecursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2)
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
//------------------------------------------------------------------------------------------------------
int TNudyEndfEnergyAng::GetMt6Neutron(int ielemId, int mt)
{
  if (fMt6Neutron.size() <= 0) return -1;
  int size = fMt6Neutron[ielemId].size();
  for (int i = 0; i < size; i++) {
    if (fMt6Neutron[ielemId][i] == mt) {
      return fMt6Neutron[ielemId][i];
    }
  }
  return -1;
}
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
int TNudyEndfEnergyAng::GetAd6(int ielemId, int mt)
{
  int size = fAD6[ielemId].size();
  for (int i = 0; i < size; i++) {
    if (fMt6Values[ielemId][i] == mt) {
      return fAD6[ielemId][i];
    }
  }
  return 0;
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfEnergyAng::GetZd6(int ielemId, int mt)
{
  int size = fZD6[ielemId].size();
  for (int i = 0; i < size; i++) {
    if (fMt6Values[ielemId][i] == mt) {
      return fZD6[ielemId][i];
    }
  }
  return 0;
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfEnergyAng::GetLaw6(int ielemId, int mt)
{
  int size = fLaw6[ielemId].size();
  for (int i = 0; i < size; i++) {
    if (fMt6Values[ielemId][i] == mt) {
      return fLaw6[ielemId][i];
    }
  }
  return -1;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::GetCos64(int ielemId, int mt, double energyK)
{
  fRnd  = new TRandom3(0);
  int i = -1;
  //   std::cout <<"mt4 values in file 6  "<< fMt4Values[ielemId].size() << std::endl;
  for (unsigned int l = 0; l < fMt4Values[ielemId].size(); l++) {
    if (fMt4Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  int min = 0;
  int max = fEnergy4OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= fEnergy4OfMts[ielemId][i][min])
    min = 0;
  else if (energyK >= fEnergy4OfMts[ielemId][i][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEnergy4OfMts[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  double fraction = (energyK - fEnergy4OfMts[ielemId][i][min]) /
                    (fEnergy4OfMts[ielemId][i][min + 1] - fEnergy4OfMts[ielemId][i][min]);
  // std::cout <<" fraction "<< fraction <<"  "<< energyK <<"  "<< fEnergy4OfMts[ielemId][i][min] << std::endl;
  double rnd1              = fRnd->Uniform(1);
  double rnd2              = fRnd->Uniform(1);
  if (rnd1 < fraction) min = min + 1;
  int k                    = 0;
  // std::cout<<" pdf size "<< fCosPdf4OfMts[ielemId][i][min].size() << std::endl;
  int size = fCosCdf4OfMts[ielemId][i][min].size();
  for (int j = 0; j < size; j++) {
    // std::cout<<"cdf "<< fCosCdf4OfMts[ielemId][i][min][2 * j ] <<"  "<< fCosCdf4OfMts[ielemId][i][min][2 * j + 1] <<
    // std::endl;
    if (rnd2 <= fCosCdf4OfMts[ielemId][i][min][j]) {
      k                    = j - 1;
      if (k >= size - 1) k = size - 1;
      break;
    }
  }
  //    for(int j = 0; j < size; j++){
  //      std::cout << fEnergy4OfMts[ielemId][i][min]<< "  "<< fMt4Values[ielemId][i]
  //      <<"  "<< fCos4OfMts[ielemId][i][min][j] <<"  "<<fCosPdf4OfMts[ielemId][i][min][j]
  //      <<"  "<<fCosCdf4OfMts[ielemId][i][min][j] << std::endl;
  //    }
  // std::cout<< k <<"  "<<fCos4OfMts[ielemId][i][min][k]<<"  "<<fCosPdf4OfMts[ielemId][i][min][k] <<"
  // "<<fCosCdf4OfMts[ielemId][i][min][ k] << std::endl;
  // std::cout<<" pdf "<<k<<"  "<< fCosPdf4OfMts[ielemId][i][min][2 * k + 3]<<"  "<< fCosPdf4OfMts[ielemId][i][min][2 *
  // k
  // + 1] << std::endl;
  // std::cout<<" cos "<< fCosPdf4OfMts[ielemId][i][min][2 * k + 2]<<"  "<< fCosPdf4OfMts[ielemId][i][min][2 * k ] <<
  // std::endl;
  double plk = (fCosPdf4OfMts[ielemId][i][min][k + 1] - fCosPdf4OfMts[ielemId][i][min][k]) /
               (fCos4OfMts[ielemId][i][min][k + 1] - fCos4OfMts[ielemId][i][min][k]);
  double plk2 = fCosPdf4OfMts[ielemId][i][min][k] * fCosPdf4OfMts[ielemId][i][min][k];
  double plsq = plk2 + 2 * plk * (rnd2 - fCosCdf4OfMts[ielemId][i][min][k]);
  double Ang  = 0;
  if (plk != 0 && plsq > 0) {
    Ang = fCos4OfMts[ielemId][i][min][k] + (sqrt(plsq) - fCosPdf4OfMts[ielemId][i][min][k]) / plk;
  } else {
    // std::cout <<"uniform angle "<< mt <<"  "<< rnd1 <<"   "<< rnd2 << std::endl;
    Ang = 2 * fRnd->Uniform(1) - 1;
  }
  return Ang;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::GetCos6(int ielemId, int mt, double energyK)
{
  fRnd  = new TRandom3(0);
  int i = -1;
  for (unsigned int l = 0; l < fMt5Values[ielemId].size(); l++) {
    // std::cout<< fMt5Values[ielemId][l] << std::endl;
    if (fMt5Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  // std::cout<<" i "<<i <<std::endl;
  int min = 0;
  int max = fEnergy5OfMts[ielemId][i].size() - 1;
  // std::cout<<"mt "<< mt <<"  "<< i <<"  "<<fLaw <<"  "<< fEnergy5OfMts[ielemId][i].size()<<"
  // "<<fEnergy5OfMts[ielemId][i][max] << std::endl;
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
  double fraction = (energyK - fEnergy5OfMts[ielemId][i][min]) /
                    (fEnergy5OfMts[ielemId][i][min + 1] - fEnergy5OfMts[ielemId][i][min]);
  // std::cout <<" fraction "<< fraction <<"  "<< energyK <<"  "<< fEnergy4OfMts[ielemId][i][min] << std::endl;
  double rnd1              = fRnd->Uniform(1);
  double rnd2              = fRnd->Uniform(1);
  if (rnd2 < fraction) min = min + 1;
  double Ang               = 0;
  // std::cout<<"energy bin "<< min <<"   "<< fEnergy5OfMts[ielemId][i][min] << std::endl;
  int k    = 0;
  int size = fCos6OfMts[ielemId][i][min].size();
  // std::cout<<"cosine size "<<  size << std::endl;
  for (unsigned int j = 0; j < fCos6OfMts[ielemId][i][min].size(); j++) {
    if (rnd1 < fCosin6Cdf[ielemId][i][min][j]) {
      k                    = j - 1;
      if (k >= size - 2) k = size - 2;
      break;
    }
  }
  // std::cout<<"cosine bin "<< k <<"  "<< fCos6OfMts[ielemId][i][min][k] << std::endl;
  double plk = (fCosin6Pdf[ielemId][i][min][k + 1] - fCosin6Pdf[ielemId][i][min][k]) /
               (fCos6OfMts[ielemId][i][min][k + 1] - fCos6OfMts[ielemId][i][min][k]);
  double plk2 = fCosin6Pdf[ielemId][i][min][k] * fCosin6Pdf[ielemId][i][min][k];
  double plsq = plk2 + 2 * plk * (rnd1 - fCosin6Cdf[ielemId][i][min][k]);
  // std::cout <<"plk "<< plk <<" plk2 "<< plk2 <<"  "<< k << std::endl;
  if (plk != 0 && plsq > 0) {
    Ang = fCos6OfMts[ielemId][i][min][k] + (sqrt(std::fabs(plsq)) - fCosin6Pdf[ielemId][i][min][k]) / plk;
    // std::cout<< Ang <<" first "<< rnd1 << std::endl;
  } else {
    Ang = 2 * rnd1 - 1;
    // std::cout<< Ang <<" sec  "<< rnd1 << std::endl;
  }
  // std::cout<< Ang << std::endl;
  return Ang;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::GetEnergy6(int ielemId, int mt, double energyK)
{
  fRnd  = new TRandom3(0);
  int i = -1;
  //    std::cout <<"mt5 values in file 6  "<< fMt5Values[ielemId].size() <<"  "<< fZD6[ielemId][0] <<"  "<<
  //    fZD6[ielemId][1] << std::endl;
  for (unsigned int l = 0; l < fMt5Values[ielemId].size(); l++) {
    if (fMt5Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  // std::cout<<"mt "<< mt <<"  "<< energyK << std::endl;
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
  double fraction = (energyK - fEnergy5OfMts[ielemId][i][min]) /
                    (fEnergy5OfMts[ielemId][i][min + 1] - fEnergy5OfMts[ielemId][i][min]);
  double rnd1              = fRnd->Uniform(1);
  double rnd2              = fRnd->Uniform(1);
  double rnd3              = fRnd->Uniform(1);
  if (rnd2 < fraction) min = min + 1;
  // std::cout<<"energy bin "<< min <<"   "<< fEnergy5OfMts[ielemId][i][min] << std::endl;
  int k    = 0;
  int size = fCos6OfMts[ielemId][i][min].size();
  for (int j = 0; j < size; j++) {
    if (rnd3 < fCosin6Cdf[ielemId][i][min][j]) {
      k                    = j - 1;
      if (k >= size - 2) k = size - 2;
      break;
    }
  }
  // std::cout<<"cosine bin "<< k <<"  "<< fCos6OfMts[ielemId][i][min][k] << std::endl;
  int m = 0;
  size  = fEnergyCdf6OfMts[ielemId][i][min][k].size();
  for (unsigned int j = 1; j < fEnergyPdf6OfMts[ielemId][i][min][k].size(); j++) {
    if (rnd1 <= fEnergyCdf6OfMts[ielemId][i][min][k][j]) {
      m                    = j - 1;
      if (m >= size - 2) m = size - 2;
      break;
    }
  }
  // std::cout<<"energy bin "<< m <<"  "<< fEnergyOut6OfMts[ielemId][i][min][k][m] << std::endl;
  double plk = (fEnergyPdf6OfMts[ielemId][i][min][k][m + 1] - fEnergyPdf6OfMts[ielemId][i][min][k][m]) /
               (fEnergyOut6OfMts[ielemId][i][min][k][m + 1] - fEnergyOut6OfMts[ielemId][i][min][k][m]);
  double plk2 = fEnergyPdf6OfMts[ielemId][i][min][k][m] * fEnergyPdf6OfMts[ielemId][i][min][k][m];
  // std::cout <<"plk "<< plk <<" plk2 "<< plk2  << std::endl;
  double edes = 0;
  if (plk != 0)
    edes = fEnergyOut6OfMts[ielemId][i][min][k][m] +
           (sqrt(plk2 + 2 * plk * (rnd1 - fEnergyCdf6OfMts[ielemId][i][min][k][m])) -
            fEnergyPdf6OfMts[ielemId][i][min][k][m]) /
               plk;
  return edes;
}
