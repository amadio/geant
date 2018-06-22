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
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyEndfMat.h"
#include "Geant/TNudyEndfNuPh.h"
#include "Geant/TNudyEndfEnergyAng.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfPhYield.h"
#include "Geant/TNudyEndfPhProd.h"
#include "Geant/TNudyEndfPhAng.h"
#include "Geant/TNudyEndfPhEnergy.h"
#include "Geant/TNudyEndfRecoPoint.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfRecoPoint)
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
  fMtValues.push_back(fMtNumbers);
  fQvalue.push_back(fQvalueTemp);
  fQvalueTemp.clear();
  fMtNumbers.clear();
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::FillPdf1D(rowd &x1, rowd &x2, rowd &x3, matrixd2 &x4, matrixd2 &x5, matrixd2 &x6)
{
  TNudyCore::Instance()->Sort(x1, x2);
//  TNudyCore::Instance()->ThinningDuplicate(x1, x2);
  TNudyCore::Instance()->cdfGenerateT(x1, x2, x3);
  for (int i = 0, x1Size = x1.size(); i != x1Size; ++i) {
    fEneE.push_back(x1[i]);
    fPdf.push_back(x2[i]);
    fCdf.push_back(x3[i]);
  }
  x4.push_back(fEneE);
  x5.push_back(fPdf);
  x6.push_back(fCdf);
  x1.clear();
  x2.clear();
  x3.clear();
  fEneE.clear();
  fPdf.clear();
  fCdf.clear();
}
//-------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::FillPdf2D()
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
//-------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ReadFile4(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
    int fMT               = sec->GetMT();
    fMtNumbers.push_back(fMT);
    int LCT = header->GetL2();
    fMtLct.push_back(LCT);
    // Tabulated probability tables
    TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
    for (int i = 0, e2N2 = tab2->GetN2(); i != e2N2; ++i) {
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      fEin.push_back(tab1->GetC2());
      ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fCosFile4, fCosPdfFile4);
      FillPdf1D(fCosFile4, fCosPdfFile4, fCosCdfFile4, fCos2D, fPdf2D, fCdf2D);
      fNbt1.clear();
      fInt1.clear();
    }
    FillPdf2D();
  } // end while loop
  fMt4Values.push_back(fMtNumbers);
  fMt4Lct.push_back(fMtLct);
  fEnergy4OfMts.push_back(fEin2D);
  fCos4OfMts.push_back(fCos3D);
  fCosPdf4OfMts.push_back(fPdf3D);
  fCosCdf4OfMts.push_back(fCdf3D);
  std::vector<int>().swap(fMtNumbers);
  std::vector<int>().swap(fMtLct);
  fEin2D.clear();
  fCos3D.clear();
  fPdf3D.clear();
  fCdf3D.clear();
  fMtNumbers.clear();
} 
//-------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ReadFile5(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    for (int k = 0; k < sec->GetN1(); k++) {
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      int MT              = sec->GetMT();
      fMtNumbers.push_back(MT);
      // arbitrary tabulated functions
      ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
      TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
      ProcessTab2(tab2, fNr2, fNp2, fNbt2, fInt2);
      for (int cr = 0; cr < fNp2; cr++) {
        TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
        fEin.push_back(tab12->GetC2());
        ProcessTab1(tab12, fNr3, fNp3, fNbt3, fInt3, fEnergyFile5, fEnergyPdfFile5);
        FillPdf1D(fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5, fEin2D, fPdf2D, fCdf2D);
        fNbt3.clear();
        fInt3.clear();
      }
      fNbt1.clear();
      fInt1.clear();
      fNbt2.clear();
      fInt2.clear();
      fE1.clear();
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
}
//-------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ReadFile6(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    int NK = sec->GetN1();
    TIter recIter(sec->GetRecords());
    for (int k = 0; k < NK; k++) {
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      div_t divr;
      AWRI        = sec->GetC2();
      int MT      = sec->GetMT();
      int ZAP     = tab1->GetC1();
      int LCT     = sec->GetL2();
      int LAW     = tab1->GetL2();
      // flag to identify if it is angular from file (4, 14) or energy distribution from file (5, 15)
      int CosOrEn = sec->GetN2(); 
      // There is no structure for these laws
      if (LAW == 3 || LAW == 4 || LAW == 0) continue;
      fLaw.push_back(LAW);
      fNR           = tab1->GetN1();
      fNP           = tab1->GetN2();
      divr          = div(ZAP, 1000);
      int particleA = divr.rem;
      int particleZ = divr.quot;
      fZd.push_back(particleZ);
      fAd.push_back(particleA);
      switch (ZAP) {
        case 0: // data for photons
        {
          fMtNumPhoton.push_back(MT); 
        }break;
        case 1: // data for neutron 
        {
          fMtNumNeutron.push_back(MT); 
        }break;
        default: // data for charge particles
        {
          fMtNumCharge.push_back(MT);
        }break;
      }
      if (ZAP != 1) continue;
      switch (LAW) {
        case 2: // this law testedwith n-005_B_010.endf file
        { 
          switch (CosOrEn) {
            case 0: // from file 6
            case 4: // from file 4 
            case 14: // from file 14 
            {
              fMtNumbers4Cos.push_back(MT);
              fMtLct4Cos.push_back(LCT);
              ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
              TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
              ProcessTab2(tab2, fNr2, fNp2, fNbt2, fInt2);
              for (int lis = 0; lis < fNp2; lis++) {
                TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
                fEin.push_back(tab11->GetC2());
                ProcessTab1(tab11, fNR, fNP, fNbt3, fInt3, fCosFile4, fCosPdfFile4);
                TNudyCore::Instance()->cdfGenerateT(fCosFile4, fCosPdfFile4, fCosCdfFile4);
                for (int i = 0, CosFile4Size = fCosFile4.size(); i != CosFile4Size; ++i) {
                  fCos.push_back(fCosFile4[i]);
                  fPdf.push_back(fCosPdfFile4[i]);
                  fCdf.push_back(fCosCdfFile4[i]);
                }
                fCos2D.push_back(fCos);
                fPdf2D.push_back(fPdf);
                fCdf2D.push_back(fCdf);
                fCosFile4.clear();
                fCosPdfFile4.clear();
                fCosCdfFile4.clear();
                fCos.clear();
                fPdf.clear();
                fCdf.clear();
                fNbt3.clear();
                fInt3.clear();
              }
              fEin2D.push_back(fEin);
              fCos3D.push_back(fCos2D);
              fPdf3D.push_back(fPdf2D);
              fCdf3D.push_back(fCdf2D);
              fEin.clear();
              fCos2D.clear();
              fPdf2D.clear();
              fCdf2D.clear();
              fNbt1.clear();
              fInt1.clear();
              fE1.clear();
              fP1.clear();
              fNbt2.clear();
              fInt2.clear();
            }break;
            case 5:
            {
              fMtNumbers5.push_back(MT);
              ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
              TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
              ProcessTab2(tab2, fNr2, fNp2, fNbt2, fInt2);
              matrixd2 Ein2D, Pdf2D, Cdf2D;
              for (int cr = 0; cr < fNp2; cr++) {
                TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
                fEin.push_back(tab12->GetC2());
                ProcessTab1(tab12, fNr3, fNp3, fNbt3, fInt3, fEnergyFile5, fEnergyPdfFile5);
                FillPdf1D(fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5, Ein2D, Pdf2D, Cdf2D);
                fNbt3.clear();
                fInt3.clear();
              }
              fEin2D5.push_back(fEin);
              fFrac2D5.push_back(fP1);
              fEne3D5.push_back(Ein2D);
              fPdf3D5.push_back(Pdf2D);
              fCdf3D5.push_back(Cdf2D);
              fEin.clear();
              Ein2D.clear();
              Pdf2D.clear();
              Cdf2D.clear();
              fP1.clear();
              fNbt1.clear();
              fInt1.clear();
              fNbt2.clear();
              fInt2.clear();
              fE1.clear();
            }break;
            case 15:
            {
              fMtNumbers15.push_back(MT);
              ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
              TNudyEndfTab2 *tab2  = (TNudyEndfTab2 *)recIter.Next();
              ProcessTab2(tab2, fNr2, fNp2, fNbt2, fInt2);
              for (int cr = 0; cr < fNp2; cr++) {
                TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
                fEin.push_back(tab12->GetC2());
                ProcessTab1(tab12, fNr3, fNp3, fNbt3, fInt3, fEnergyFile5, fEnergyPdfFile5);
                FillPdf1D(fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5, fEin2D, fPdf2D, fCdf2D);
                fNbt3.clear();
                fInt3.clear();
              }
              fEin2D15.push_back(fEin);
              fFrac2D15.push_back(fP1);
              fEin3D15.push_back(fEin2D);
              fPdf3D15.push_back(fPdf2D);
              fCdf3D15.push_back(fCdf2D);
              fEin.clear();
              fEin2D.clear();
              fPdf2D.clear();
              fCdf2D.clear();
              fNbt1.clear();
              fInt1.clear();
              fNbt2.clear();
              fInt2.clear();
              fE1.clear();
              fP1.clear();
            }break;
          }
        }break;
        case 3:
        case 4: 
        {
          for (int cr = 0; cr < fNR; cr++) {
            fNbt1.push_back(tab1->GetNBT(cr));
            fInt1.push_back(tab1->GetINT(cr));
          }
          for (int crs = 0; crs < fNP; crs++) {
            fE1.push_back(tab1->GetX(crs));
            fP1.push_back(tab1->GetY(crs));
            // std::cout<<"energy1 "<< fE1[crs] <<" multip1 "<< fP1[crs] << std::endl;
          }
          fNbt1.clear();
          fInt1.clear();
          fE1.clear();
          fP1.clear();
          //---------------------------------------------
        } break;
        case 5: // charge particle scattering law is not implimented yet
        {
        } break;
        case 7: // laboratory angle and energy law
        {
          fMtNumbers6.push_back(MT);
          fMtLct.push_back(1);
          ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
          TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
          ProcessTab2(tab2, fNr2, fNp2, fNbt2, fInt2);
          for (int cr1 = 0; cr1 < fNp2; cr1++) {
            TNudyEndfTab2 *tab3 = (TNudyEndfTab2 *)recIter.Next();
            fEin.push_back(tab3->GetC2());
            int NMU = tab3->GetN2();
            for (int cr2 = 0; cr2 < NMU; cr2++) {
              TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
              ProcessTab1(tab12, fNr3, fNp3, fNbt3, fInt3, fEnergyFile5, fEnergyPdfFile5);
              fCosIn.push_back(tab12->GetC2());
              double fsum = 0.0;
              for (int cr = 1; cr < fNp3; cr++) {
                fsum += 0.5 * (fEnergyPdfFile5[cr] + fEnergyPdfFile5[cr - 1]) 
                * (fEnergyFile5[cr] - fEnergyFile5[cr - 1]);
              }
              fCosInPdf.push_back(fsum);
              TNudyCore::Instance()->cdfGenerateT(fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5);
              for (int i = 0, EnergyFile5Size = fEnergyFile5.size(); i != EnergyFile5Size; ++i) {
                fEoute.push_back(fEnergyFile5[i]);
                fPdfe.push_back(fEnergyPdfFile5[i]);
                fCdfe.push_back(fEnergyCdfFile5[i]);
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
}
//-------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ReadFile8(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    int MT = sec->GetMT();
    int LE = sec->GetL1();
    switch (MT) {
      case 454: // Neutron induced independent fission yield
      { 
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
      }break;
      case 459: // Neutron induced cummulative fission yield
      { 
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
      }break;
    }
  }
}
//-------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ReadFile14(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    int MT      = sec->GetMT();
    fMtNumbers.push_back(MT);
    mtf[0]  = sec->GetMAT();
    mtf[1]  = sec->GetMT();
    mtf[2]  = sec->GetMF();
    TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
    for (int i = 0, eN2 = tab2->GetN2(); i != eN2; ++i) {
      TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
      ProcessTab1(tab, fNR, fNP, fNbt1, fInt1, fCosFile4, fCosPdfFile4);
      FillPdf1D(fCosFile4, fCosPdfFile4, fCosCdfFile4, fCos2D, fPdf2D, fCdf2D);
      fNbt1.clear();
      fInt1.clear();;
    }
    FillPdf2D();            
  } // end of while loop
  fMt4ValuesPhoton.push_back(fMtNumbers);
  fMt4LctPhoton.push_back(fMtLct);
  fEnergy4OfMtsPhoton.push_back(fEin2D);
  fCos4OfMtsPhoton.push_back(fCos3D);
  fCosPdf4OfMtsPhoton.push_back(fPdf3D);
  fCosCdf4OfMtsPhoton.push_back(fCdf3D);
  std::vector<int>().swap(fMtNumbers);
  std::vector<int>().swap(fMtLct);
  fEin2D.clear();
  fCos3D.clear();
  fPdf3D.clear();
  fCdf3D.clear();
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ReadFile15(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    int NC                = sec->GetN1();
    for (int k = 0; k < NC; k++) {
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      int MT              = sec->GetMT();
      fMtNumbers.push_back(MT);
      fNR       = tab1->GetN1();
      fNP       = tab1->GetN2();
      // arbitrary tabulated function (only LF=1 is defined at this moment)
      ProcessTab1(tab1, fNR, fNP, fNbt1, fInt1, fE1, fP1);
      TNudyEndfTab2 *tab2  = (TNudyEndfTab2 *)recIter.Next();
      ProcessTab2(tab2, fNr2, fNp2, fNbt2, fInt2);
      for (int cr = 0; cr < fNp2; cr++) {
        TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
        fEin.push_back(tab12->GetC2());
        ProcessTab1(tab12, fNr3, fNp3, fNbt3, fInt3, fEnergyFile5, fEnergyPdfFile5);
        FillPdf1D(fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5, fEin2D, fPdf2D, fCdf2D);
        fNbt3.clear();
        fInt3.clear();
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
    }
    fE1.clear();
    fP1.clear();
    fNbt1.clear();
    fInt1.clear();
    fNbt2.clear();
    fInt2.clear();
  }
  fMt5ValuesPhoton.push_back(fMtNumbers);
  fEnergy5OfMtsPhoton.push_back(fEin2D);
  fFraction5OfMtsPhoton.push_back(fFrac2D);
  fEnergyOut5OfMtsPhoton.push_back(fEne3D);
  fEnergyPdf5OfMtsPhoton.push_back(fPdf3D);
  fEnergyCdf5OfMtsPhoton.push_back(fCdf3D);
  fMtNumbers.clear();
  fEin2D.clear();
  fEne3D.clear();
  fPdf3D.clear();
  fCdf3D.clear();
  fFrac2D.clear();
}
// ------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::GetData(int ielemId, const char *rENDF)
{
  fElemId     = ielemId;
  if (fElementId.size() > 0) {
    for (int i = 0, ElementIdSize = fElementId.size(); i != ElementIdSize; ++i) {
      if (ielemId == fElementId[i]) {
        return;
      }
    }
  }
  fElementId.push_back(ielemId);
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
      // Read File data into class structures
      switch (file->GetMF()) {
        case 1:
          if (fFlagRead != 1) {
            fRecoNuPh  = new TNudyEndfNuPh(file);
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

          fEnergyMts.clear();
          fSigmaMts.clear();
        } break;
        case 4:
          ReadFile4(file);
          // std::cout << "file 4 OK " << std::endl;
          break;
        case 5:
          ReadFile5(file);
          // std::cout << "file 5 OK " << std::endl;
          break;
        case 6:
          // std::cout << "before file 6 " << std::endl;
          ReadFile6(file);
          // std::cout << "file 6 OK " << std::endl;
          break;
        case 8:
          ReadFile8(file);
          break;
          ///*
          // case 12:
          // //         std::cout << "before file 12 " << std::endl;
          // //     fRecoPhYield = new TNudyEndfPhYield(file);
          // //     std::cout<<"file 12 OK "<<std::endl;
          // //     break;
          // //       case 13:
          // //         std::cout << "before file 13 " << std::endl;
          // //     fRecoPhProd = new TNudyEndfPhProd(file);
          // //     std::cout<<"file 13 OK "<<std::endl;
          // //     break;
        case 14:
          ReadFile14(file);
          std::cout<<"file 14 OK "<<std::endl;
          break;
        case 15:
          //ReadFile15(file);
          // std::cout<<"file 15 OK "<<std::endl;
          break;
      }
    }
    // cosine distributions
    fMt4Values.push_back(fMtNumbers4Cos);
    fMt4Lct.push_back(fMtLct4Cos);
    fEnergy4OfMts.push_back(fEin2D);
    fCos4OfMts.push_back(fCos3D);
    fCosPdf4OfMts.push_back(fPdf3D);
    fCosCdf4OfMts.push_back(fCdf3D);
    fEin2D.clear();
    fCos3D.clear();
    fPdf3D.clear();
    fCdf3D.clear();
    fMtLct4Cos.clear();
    fMtNumbers4Cos.clear();
    // energy distributions
    fMt5Values.push_back(fMtNumbers5);
    fEnergy5OfMts.push_back(fEin2D5);
    fFraction5OfMts.push_back(fFrac2D5);
    fEnergyOut5OfMts.push_back(fEne3D5);
    fEnergyPdf5OfMts.push_back(fPdf3D5);
    fEnergyCdf5OfMts.push_back(fCdf3D5);
    fMtNumbers5.clear();
    fEin2D5.clear();
    fEne3D5.clear();
    fPdf3D5.clear();
    fCdf3D5.clear();
    fFrac2D5.clear();  
    // angle-energy distributions
    fCos6OfMts.push_back(fCos3d);
    fCosin6Pdf.push_back(fCosinpdf3d);
    fCosin6Cdf.push_back(fCosincdf3d);
    fEnergy6OfMts.push_back(fEin2d);
    fEnergyOut6OfMts.push_back(fEout4de);
    fEnergyPdf6OfMts.push_back(fPdf4de);
    fEnergyCdf6OfMts.push_back(fCdf4de);
    fMt6Values.push_back(fMtNumbers6);
    fMt6Neutron.push_back(fMtNumNeutron);
    //  fMt4Lct.push_back(fMtLct);
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
    // photon energy distribution
    fMt5ValuesPhoton.push_back(fMtNumbers15);
    fEnergy5OfMtsPhoton.push_back(fEin2D15);
    fFraction5OfMtsPhoton.push_back(fFrac2D15);
    fEnergyOut5OfMtsPhoton.push_back(fEin3D15);
    fEnergyPdf5OfMtsPhoton.push_back(fPdf3D15);
    fEnergyCdf5OfMtsPhoton.push_back(fCdf3D15);
    fMtNumbers15.clear();
    fEin2D15.clear();
    fEin3D15.clear();
    fPdf3D15.clear();
    fCdf3D15.clear();
    fFrac2D15.clear();
    
  }
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::FixupTotal(std::vector<double> &x1)
{
  for (int i = 0, SigmaOfMtsSize = fSigmaOfMts.size(); i != SigmaOfMtsSize; ++i) {
    int size    = fSigmaOfMts[i].size() / 2;
    int flagLoc = -1;
    for (int k = 0, x1Size = x1.size(); k != x1Size; ++k) {
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
int TNudyEndfRecoPoint::GetElementId(int ielemId)
{
for (int i = 0, ElementIdSize = fElementId.size(); i != ElementIdSize; ++i) {
  if (ielemId == fElementId[i]) {
    ielemId = i;
    break;
  }
}
return ielemId;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetSigmaTotal(int ielemId, double energyK)
{
  ielemId = GetElementId(ielemId);
  if (fMtValues[ielemId].size() <= 0) return 0;
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
  if (min == (int)fSigUniT[ielemId].size()) return fSigUniT[ielemId][min];
  if (min > (int)fSigUniT[ielemId].size()) return 0;
  return fSigUniT[ielemId][min] +
  (fSigUniT[ielemId][min + 1] - fSigUniT[ielemId][min]) * (energyK - fEneUni[ielemId][min]) /
  (fEneUni[ielemId][min + 1] - fEneUni[ielemId][min]);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetSigmaPartial(int ielemId, int mt, double energyK)
{
  ielemId = GetElementId(ielemId);
  int i = -1;
  if (fMtValues[ielemId].size() <= 0) return 0;
  for (int l = 0, MtValuesSize = fMtValues[ielemId].size(); l != MtValuesSize; ++l) {
    if (fMtValues[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  if (i < 0) return 0;
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
  if (minp <= 0) return 0;
  if (minp >= (int)fSigUniOfMt[ielemId][i].size()) return 0;
  return fSigUniOfMt[ielemId][i][minp] +
  (fSigUniOfMt[ielemId][i][minp + 1] - fSigUniOfMt[ielemId][i][minp]) * (energyK - fEneUni[ielemId][min]) /
  (fEneUni[ielemId][min + 1] - fEneUni[ielemId][min]);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetCos4(int ielemId, int mt, double energyK)
{
  ielemId = GetElementId(ielemId);
  int i = -1;
  if (fMt4Values[ielemId].size() <= 0) return GetCos6(ielemId, mt, energyK);
  for (int l = 0, Mt4Size = fMt4Values[ielemId].size(); l != Mt4Size; ++l) {
    if (fMt4Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  if (i < 0) return GetCos6(ielemId, mt, energyK);
  int min = 0;
  int max = fEnergy4OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= fEnergy4OfMts[ielemId][i][min])
    min = 0;
  else if (energyK >= fEnergy4OfMts[ielemId][i][max])
    min = max -1;
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
  double rnd1              = fRng.uniform();
  double rnd2              = fRng.uniform();
  if (rnd2 < fraction) min = min + 1;
  if (min > max) min = max;
  int k                    = 0;
  int size = fCosCdf4OfMts[ielemId][i][min].size();
  for (int j = 1; j < size; j++) {
    if (rnd1 <= fCosCdf4OfMts[ielemId][i][min][j]) {
      k                    = j - 1;
      if (k >= size - 2) k = size - 2;
      break;
    }
  }
  double plk = (fCosPdf4OfMts[ielemId][i][min][k + 1] - fCosPdf4OfMts[ielemId][i][min][k]) /
  (fCos4OfMts[ielemId][i][min][k + 1] - fCos4OfMts[ielemId][i][min][k]);
  double plk2 = fCosPdf4OfMts[ielemId][i][min][k] * fCosPdf4OfMts[ielemId][i][min][k];
  double plsq = plk2 + 2 * plk * (rnd1 - fCosCdf4OfMts[ielemId][i][min][k]);
  double Ang  = 0;
  if (plk != 0 && plsq > 0) {
    Ang = fCos4OfMts[ielemId][i][min][k] + (sqrt(std::fabs(plsq)) - fCosPdf4OfMts[ielemId][i][min][k]) / plk;
  } else {
    Ang = 2 * rnd1 - 1;
  }
  return Ang;
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetCos4Lct(int ielemId, int mt)
{
  ielemId = GetElementId(ielemId);
  int i = 0;
  for (int l = 0, Mt4Size = fMt4Values[ielemId].size(); l != Mt4Size; ++l) {
    if (fMt4Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  return fMt4Lct[ielemId][i];
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetEnergy5(int ielemId, int mt, double energyK)
{
  ielemId = GetElementId(ielemId);
  int i = -1;
  if (fMt5Values[ielemId].size() <= 0) return GetEnergy6(ielemId, mt, energyK);
  for (int l = 0, Mt5Size = fMt5Values[ielemId].size(); l != Mt5Size; ++l) {
    if (fMt5Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  if (i < 0) return GetEnergy6(ielemId, mt, energyK);
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
  if (min < 0) min = 0;
  double fraction = (energyK - fEnergy5OfMts[ielemId][i][min]) /
  (fEnergy5OfMts[ielemId][i][min + 1] - fEnergy5OfMts[ielemId][i][min]);
  double rnd1              = fRng.uniform();
  double rnd2              = fRng.uniform();
  if (rnd2 < fraction) min = min + 1;
  if (min > max) min = max;
  int k = 0;
  int size = fEnergyCdf5OfMts[ielemId][i][min].size();
  for (int j = 0, ESize = fEnergyPdf5OfMts[ielemId][i][min].size(); j != ESize; ++j) {
    if (rnd1 <= fEnergyCdf5OfMts[ielemId][i][min][j]) {
      k                    = j;
      if (k >= size - 1) k = size - 1;
      break;
    }
  }
  double plk = (fEnergyPdf5OfMts[ielemId][i][min][k] - fEnergyPdf5OfMts[ielemId][i][min][k - 1]) /
  (fEnergyOut5OfMts[ielemId][i][min][k] - fEnergyOut5OfMts[ielemId][i][min][k - 1]);
  double plk2 = fEnergyPdf5OfMts[ielemId][i][min][k - 1] * fEnergyPdf5OfMts[ielemId][i][min][k - 1];
  
  double edes = 0;
  if (plk != 0)
    edes = fEnergyOut5OfMts[ielemId][i][min][k - 1] +
    (sqrt(plk2 + 2 * plk * (rnd1 - fEnergyCdf5OfMts[ielemId][i][min][k - 1])) -
    fEnergyPdf5OfMts[ielemId][i][min][k - 1]) /
    plk;
  if (min == max) return edes;
  double emin = fEnergyOut5OfMts[ielemId][i][min][1] +
  fraction * (fEnergyOut5OfMts[ielemId][i][min + 1][1] - fEnergyOut5OfMts[ielemId][i][min][1]);
  double emax = fEnergyOut5OfMts[ielemId][i][min][size - 1] +
  fraction * (fEnergyOut5OfMts[ielemId][i][min + 1][fEnergyCdf5OfMts[ielemId][i][min + 1].size() - 1] -
  fEnergyOut5OfMts[ielemId][i][min][size - 1]);
  edes = fEnergyOut5OfMts[ielemId][i][min][1] +
  (edes - emin) * (fEnergyOut5OfMts[ielemId][i][min][size - 1] - fEnergyOut5OfMts[ielemId][i][min][1]) /
  (emax - emin);
  return edes;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetDelayedFraction(int ielemId, int mt, double energyK)
{
  ielemId = GetElementId(ielemId);
  int i = -1;
  for (int l = 0, Mt5Size = fMt5Values[ielemId].size(); l != Mt5Size; ++l) {
    if (fMt5Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  if (i < 0) return -99;
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
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetQValue(int ielemId, int mt)
{
  ielemId = GetElementId(ielemId);
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
  ielemId = GetElementId(ielemId);
  if (fMt4.size() <= 0) return -99;
  int size = fMt4[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (fMt4[ielemId][2 * i + 1] == mt) {
      return fMt4[ielemId][2 * i];
    }
  }
  return -99;
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetMt5(int ielemId, int mt)
{
  ielemId = GetElementId(ielemId);
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
double TNudyEndfRecoPoint::GetMt6(int ielemId, int mt)
{
  ielemId = GetElementId(ielemId);
  if (fMt6.size() <= 0) return -99;
  int size = fMt6[ielemId].size() / 2;
  for (int i = 0; i < size; i++) {
    if (fMt6[ielemId][2 * i + 1] == mt) {
      return fMt6[ielemId][2 * i];
    }
  }
  return -99;
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetMt6Neutron(int ielemId, int mt)
{
  ielemId = GetElementId(ielemId);
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
int TNudyEndfRecoPoint::GetAd6(int ielemId, int mt)
{
  ielemId = GetElementId(ielemId);
  int size = fAD6[ielemId].size();
  for (int i = 0; i < size; i++) {
    if (fMt6Values[ielemId][i] == mt) {
      return fAD6[ielemId][i];
    }
  }
  return 0;
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetZd6(int ielemId, int mt)
{
  ielemId = GetElementId(ielemId);
  int size = fZD6[ielemId].size();
  for (int i = 0; i < size; i++) {
    if (fMt6Values[ielemId][i] == mt) {
      return fZD6[ielemId][i];
    }
  }
  return 0;
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetLaw6(int ielemId, int mt)
{
  ielemId = GetElementId(ielemId);
  int size = fLaw6[ielemId].size();
  for (int i = 0; i < size; i++) {
    if (fMt6Values[ielemId][i] == mt) {
      return fLaw6[ielemId][i];
    }
  }
  return -1;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetCos6(int ielemId, int mt, double energyK)
{
  int i = -1;
  if (fMt6Values[ielemId].size() <= 0) return -99;
  for (int l = 0, Mt6Size = fMt6Values[ielemId].size(); l != Mt6Size; ++l) {
    if (fMt6Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  if (i < 0) return -99;
  int min = 0;
  int max = fEnergy6OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= fEnergy6OfMts[ielemId][i][min])
    min = 0;
  else if (energyK >= fEnergy6OfMts[ielemId][i][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEnergy6OfMts[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  double fraction = (energyK - fEnergy6OfMts[ielemId][i][min]) /
  (fEnergy6OfMts[ielemId][i][min + 1] - fEnergy6OfMts[ielemId][i][min]);
  double rnd1              = fRng.uniform();
  double rnd2              = fRng.uniform();
  if (rnd2 < fraction) min = min + 1;
  if (min > max) min = max;
  double Ang               = 0;
  int k    = 0;
  int size = fCos6OfMts[ielemId][i][min].size();
  for (int j = 0; j != size; ++j) {
    if (rnd1 < fCosin6Cdf[ielemId][i][min][j]) {
      k                    = j - 1;
      if (k >= size - 2) k = size - 2;
      break;
    }
  }
  double plk = (fCosin6Pdf[ielemId][i][min][k + 1] - fCosin6Pdf[ielemId][i][min][k]) /
  (fCos6OfMts[ielemId][i][min][k + 1] - fCos6OfMts[ielemId][i][min][k]);
  double plk2 = fCosin6Pdf[ielemId][i][min][k] * fCosin6Pdf[ielemId][i][min][k];
  double plsq = plk2 + 2 * plk * (rnd1 - fCosin6Cdf[ielemId][i][min][k]);
  if (plk != 0 && plsq > 0) {
    Ang = fCos6OfMts[ielemId][i][min][k] + (sqrt(std::fabs(plsq)) - fCosin6Pdf[ielemId][i][min][k]) / plk;
  } else {
    Ang = 2 * rnd1 - 1;
  }
  return Ang;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetEnergy6(int ielemId, int mt, double energyK)
{
  int i = -1;
  if (fMt6Values[ielemId].size() <= 0) return -99;
  for (int l = 0, Mt6Size = fMt6Values[ielemId].size(); l != Mt6Size; ++l) {
    if (fMt6Values[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  if (i < 0) return -99;
  int min = 0;
  int max = fEnergy6OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= fEnergy6OfMts[ielemId][i][min])
    min = 0;
  else if (energyK >= fEnergy6OfMts[ielemId][i][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEnergy6OfMts[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  double fraction = (energyK - fEnergy6OfMts[ielemId][i][min]) /
  (fEnergy6OfMts[ielemId][i][min + 1] - fEnergy6OfMts[ielemId][i][min]);
  double rnd1              = fRng.uniform();
  double rnd2              = fRng.uniform();
  double rnd3              = fRng.uniform();
  if (rnd2 < fraction) min = min + 1;
  if (min > max) min = max;
  int k    = 0;
  int size = fCos6OfMts[ielemId][i][min].size();
  for (int j = 0; j < size; j++) {
    if (rnd3 < fCosin6Cdf[ielemId][i][min][j]) {
      k                    = j - 1;
      if (k >= size - 2) k = size - 2;
      break;
    }
  }
  int m = 0;
  size  = fEnergyCdf6OfMts[ielemId][i][min][k].size();
  for (int j = 1, ESize = fEnergyPdf6OfMts[ielemId][i][min][k].size(); j != ESize; ++j) {
    if (rnd1 <= fEnergyCdf6OfMts[ielemId][i][min][k][j]) {
      m                    = j - 1;
      if (m >= size - 2) m = size - 2;
      break;
    }
  }
  double plk = (fEnergyPdf6OfMts[ielemId][i][min][k][m + 1] - fEnergyPdf6OfMts[ielemId][i][min][k][m]) /
  (fEnergyOut6OfMts[ielemId][i][min][k][m + 1] - fEnergyOut6OfMts[ielemId][i][min][k][m]);
  double plk2 = fEnergyPdf6OfMts[ielemId][i][min][k][m] * fEnergyPdf6OfMts[ielemId][i][min][k][m];
  // std::cout <<"plk "<< plk <<" plk2 "<< plk2  << std::endl;
  double edes = energyK;
  if (plk != 0) {
    edes = fEnergyOut6OfMts[ielemId][i][min][k][m] +
    (sqrt(plk2 + 2 * plk * (rnd1 - fEnergyCdf6OfMts[ielemId][i][min][k][m])) -
    fEnergyPdf6OfMts[ielemId][i][min][k][m]) /
    plk;
  }
  return edes;
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetFisYield(int ielemId, double energyK)
{
  ielemId = GetElementId(ielemId);
  int min = 0;
  int max = fEinfId[ielemId].size() - 2;
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
  double rnd1              = fRng.uniform();
  double rnd2              = fRng.uniform();
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
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetNuTotal(int ielemId, double energyK)
{
  ielemId = GetElementId(ielemId);
  return fRecoNuPh->GetNuTotal(ielemId, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetNuPrompt(int ielemId, double energyK)
{
  ielemId = GetElementId(ielemId);
  return fRecoNuPh->GetNuPrompt(ielemId, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetNuDelayed(int ielemId, double energyK)
{
  ielemId = GetElementId(ielemId);
  return fRecoNuPh->GetNuDelayed(ielemId, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetFissHeat(int ielemId, double energyK)
{
  ielemId = GetElementId(ielemId);
  return fRecoNuPh->GetFissHeat(ielemId, energyK);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetLambdaD(int ielemId, int time)
{
  ielemId = GetElementId(ielemId);
  return fRecoNuPh->GetLambdaD(ielemId, time);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetCos4Photon(int ielemId, int mt, double energyK)
{
  ielemId = GetElementId(ielemId);
  int i = -1;
  for (int l = 0, Mt4ValuesSize = fMt4ValuesPhoton[ielemId].size(); l != Mt4ValuesSize; ++l) {
    if (fMt4ValuesPhoton[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  int min = 0;
  int max = fEnergy4OfMtsPhoton[ielemId][i].size() - 2;
  int mid = 0;
  if (energyK <= fEnergy4OfMtsPhoton[ielemId][i][min])
    min = 0;
  else if (energyK >= fEnergy4OfMtsPhoton[ielemId][i][max])
    min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEnergy4OfMtsPhoton[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  double fraction = (energyK - fEnergy4OfMtsPhoton[ielemId][i][min]) / 
  (fEnergy4OfMtsPhoton[ielemId][i][min + 1] - fEnergy4OfMtsPhoton[ielemId][i][min]);
  double rnd1              = fRng.uniform();
  double rnd2              = fRng.uniform();
  if (rnd2 < fraction) min = min + 1;
  int k                    = 0;
  int size = fCosCdf4OfMtsPhoton[ielemId][i][min].size();
  for (int j = 1; j < size; j++) {
    if (rnd1 <= fCosCdf4OfMtsPhoton[ielemId][i][min][j]) {
      k                    = j - 1;
      if (k >= size - 2) k = size - 2;
      break;
    }
  }
  double plk = (fCosPdf4OfMtsPhoton[ielemId][i][min][k + 1] - fCosPdf4OfMtsPhoton[ielemId][i][min][k]) /
  (fCos4OfMtsPhoton[ielemId][i][min][k + 1] - fCos4OfMtsPhoton[ielemId][i][min][k]);
  double plk2       = fCosPdf4OfMtsPhoton[ielemId][i][min][k] * fCosPdf4OfMtsPhoton[ielemId][i][min][k];
  double plsq       = plk2 + 2 * plk * (rnd1 - fCosCdf4OfMtsPhoton[ielemId][i][min][k]);
  double Ang        = 0;
  if (plk == 0) Ang = fCos4OfMtsPhoton[ielemId][i][min][k];
  if (plk != 0 && plsq > 0) {
    Ang = fCos4OfMtsPhoton[ielemId][i][min][k] + (sqrt(std::fabs(plsq)) - 
    fCosPdf4OfMtsPhoton[ielemId][i][min][k]) / plk;
  } else {
    Ang = 2 * rnd1 - 1;
  }
  return Ang;
}
//________________________________________________________________________________________________________________
int TNudyEndfRecoPoint::GetCos4LctPhoton(int ielemId, int mt)
{
  ielemId = GetElementId(ielemId);
  int i = 0;
  for (int l = 0, Mt4ValuesSize = fMt4ValuesPhoton[ielemId].size(); l != Mt4ValuesSize; ++l) {
    if (fMt4ValuesPhoton[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  return fMt4LctPhoton[ielemId][i];
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetEnergy5Photon(int ielemId, int mt, double energyK)
{
  int i = -1;
  for (int l = 0, lmax = fMt5ValuesPhoton[ielemId].size(); l != lmax; ++l) {
    if (fMt5ValuesPhoton[ielemId][l] == mt) {
      i = l;
      break;
    }
  }
  if (i < 0) return -99;
  int min = 0;
  int max = fEnergy5OfMtsPhoton[ielemId][i].size() - 2;
  int mid = 0;
  if (energyK <= fEnergy5OfMtsPhoton[ielemId][i][min])
    min = 0;
  else if (energyK >= fEnergy5OfMtsPhoton[ielemId][i][max])
    min = max;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < fEnergy5OfMtsPhoton[ielemId][i][mid])
        max = mid;
      else
        min = mid;
    }
  }
  if (min < 0) min = 0;
  double fraction = (energyK - fEnergy5OfMtsPhoton[ielemId][i][min]) /
  (fEnergy5OfMtsPhoton[ielemId][i][min + 1] - fEnergy5OfMtsPhoton[ielemId][i][min]);
  double rnd1              = fRng.uniform();
  double rnd2              = fRng.uniform();
  if (rnd2 < fraction) min = min + 1;
  int k = 0;
  int size = fEnergyCdf5OfMtsPhoton[ielemId][i][min].size();
  for (int j = 0, jmax = fEnergyPdf5OfMtsPhoton[ielemId][i][min].size(); j != jmax; ++j) {
    if (rnd1 <= fEnergyCdf5OfMtsPhoton[ielemId][i][min][j]) {
      k                    = j;
      if (k >= size - 1) k = size - 1;
      break;
    }
  }
  double plk = (fEnergyPdf5OfMtsPhoton[ielemId][i][min][k] - fEnergyPdf5OfMtsPhoton[ielemId][i][min][k - 1]) /
  (fEnergyOut5OfMtsPhoton[ielemId][i][min][k] - fEnergyOut5OfMtsPhoton[ielemId][i][min][k - 1]);
  double plk2 = fEnergyPdf5OfMtsPhoton[ielemId][i][min][k - 1] * fEnergyPdf5OfMtsPhoton[ielemId][i][min][k - 1];
  
  double edes = 0;
  if (plk != 0)
    edes = fEnergyOut5OfMtsPhoton[ielemId][i][min][k - 1] +
    (sqrt(plk2 + 2 * plk * (rnd1 - fEnergyCdf5OfMtsPhoton[ielemId][i][min][k - 1])) -
    fEnergyPdf5OfMtsPhoton[ielemId][i][min][k - 1]) /
    plk;
  double emin = fEnergyOut5OfMtsPhoton[ielemId][i][min][1] +
  fraction * (fEnergyOut5OfMtsPhoton[ielemId][i][min + 1][1] - fEnergyOut5OfMtsPhoton[ielemId][i][min][1]);
  double emax = fEnergyOut5OfMtsPhoton[ielemId][i][min][size - 1] +
  fraction * (fEnergyOut5OfMtsPhoton[ielemId][i][min + 1][fEnergyCdf5OfMtsPhoton[ielemId][i][min + 1].size() - 1] -
  fEnergyOut5OfMtsPhoton[ielemId][i][min][size - 1]);
  edes = fEnergyOut5OfMtsPhoton[ielemId][i][min][1] +
  (edes - emin) * (fEnergyOut5OfMtsPhoton[ielemId][i][min][size - 1] - fEnergyOut5OfMtsPhoton[ielemId][i][min][1]) /
  (emax - emin);
  return edes;
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ModifyTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x)
{
  int NP = x1.size();
  secTab1->SetCont(0,x,0,0,1,NP);
  secTab1->SetNBT(NP,0);
  secTab1->SetINT(2,0);
  int np = 0;
  for (int j = 0; j < (NP - 1) / 3 + 1; ++j) {
    for (int k = 0; k < 3; ++k) {
      secTab1->SetX(x1[3 * j + k], np);
      secTab1->SetY(x2[3 * j + k], np);
      if (++np >= NP)break;
    }
  }
}  
// -------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ProcessTab2(Nudy::TNudyEndfTab2 *tab2, int &NR,int &NP,rowint &fnbt, rowint &fint)
{
  NR   = tab2->GetN1();
  NP   = tab2->GetN2();
  for (int cr = 0; cr < NR; cr++) {
    fnbt.push_back(tab2->GetNBT(cr));
    fint.push_back(tab2->GetINT(cr));
  }
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, 
                                     rowint &fint,rowd &x1, rowd &x2)
{  
  NR = tab1->GetNR();
  NP = tab1->GetNP();
  for (int i = 0; i < NR; i++) {
    fnbt.push_back(tab1->GetNBT(i));
    fint.push_back(tab1->GetINT(i));
  }
  for (int crs = 0; crs < NP; crs++) {
    x1.push_back(tab1->GetX(crs));
    x2.push_back(tab1->GetY(crs));
  }
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::CreateTab2(TNudyEndfTab2 *secTab2, int &NE)
{
  secTab2->SetCont(0,0,0,0,1,NE);
  secTab2->SetNBT(NE,0);
  secTab2->SetINT(2,0);
}
// -------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::CreateTab1(TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x)
{
  int NP  = x1.size();
  secTab1->SetCont(0,x,0,0,1,NP);
  secTab1->SetContMF(mtf[0],mtf[1],mtf[2]);
  secTab1->SetNBT(NP,0);
  secTab1->SetINT(2,0);
  int np = 0;
  for (int j = 0; j < (NP - 1) / 3 + 1; ++j) {
    for (int k = 0; k < 3; ++k) {
      secTab1->SetX(x1[3 * j + k], np);
      secTab1->SetY(x2[3 * j + k], np);
      if (++np >= NP)break;
    }
  }
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
