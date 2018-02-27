// 	Selection of secondary particle energy, angle
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: June 22, 2016
#include <iostream>
//#include "Geant/TNudyEndfDoppler.h"
//#include "Geant/TNudyEndfAng.h"
//#include "Geant/TNudyEndfEnergy.h"
//#include "Geant/TNudyEndfEnergyAng.h"
//#include "Geant/TNudyEndfFissionYield.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfRecoPoint.h"
#include "Geant/TNudySampling.h"
#ifdef USE_ROOT
#include "TRandom3.h"
#endif
#include "TCanvas.h"
#include "TFile.h"
#ifdef USE_ROOT
ClassImp(TNudySampling)
#endif

    TNudySampling::TNudySampling()
{
}
//------------------------------------------------------------------------------------------------------
TNudySampling::TNudySampling(Particle *particle, TNudyEndfRecoPoint *recoPoint)
{
  kineticE = particle[elemId].energy;
  events   = 1;
  fRnd     = new TRandom3(0);
  div_t divr;
  TFile *f = new TFile("test.root", "recreate");
  f->cd();
  h     = new TH2D("h", "", 100, -1, 1, 1000, 1E-13, kineticE / 1E9);
  h1    = new TH1D("h1", "", 100, -1, 1);
  h2    = new TH1D("h2", "", 1000, 0, kineticE / 1E9);
  fissA = new TH1D("fissA", "", 243, 0, 242);
  for (int i = 0; i < 10; i++) {
    fissA1[i] = new TH1D(Form("fissA1[%d]", i), "", 1000, 0, 20);
    hist[i]   = new TH2D(Form("hist[%d]", i), "", 1000, 0, 20, 10, 0, 10);
  }
  // std::cout <<"sampling sigmaTotal "<< recoPoint->GetSigmaTotal(elemId,kineticE) << std::endl;
  // std::cout <<"sampling sigmaPartial total "<< recoPoint->GetSigmaPartial(0,0,20) << std::endl;
  // std::cout <<"sampling sigmaPartial elstic "<< recoPoint->GetSigmaPartial(0,1,20) << std::endl;
  // determining reaction type from element;
  // double density = 1;
  // double charge = 1;
  // double avg = 6.022E23;
  // double ro = avg * density / mass;
  // int enemax = recoPoint->eneUni[elemId].size();
  // std::cout <<" target mass "<< particle[elemId].mass <<"  "<< particle[elemId].charge <<"  "<< kineticE <<
  // std::endl;
  // kineticE = fRnd->Uniform(1) * recoPoint->eneUni[elemId][enemax-1];
  for (unsigned int crsp = 0; crsp < recoPoint->MtValues[elemId].size(); crsp++) {
    crs.push_back(recoPoint->GetSigmaPartial(elemId, crsp, kineticE) / recoPoint->GetSigmaTotal(elemId, kineticE));
    // std::cout<<recoPoint->MtValues[elemId][crsp]<<"  "<< crs[crsp] <<"  "<<
    // recoPoint->GetSigmaPartial(elemId,crsp,kineticE) <<"  "<< recoPoint->GetSigmaTotal(elemId,kineticE) << std::endl;
  }
  do {
    double sum1 = 0;
    std::cout << counter << " bef " << recoPoint->MtValues[elemId].size() << std::endl;
    double rnd1 = fRnd->Uniform(1);
    // std::cout << counter <<" aft "<<rnd1<< std::endl;
    for (unsigned int crsp = 0; crsp < recoPoint->MtValues[elemId].size(); crsp++) {
      sum1 += crs[crsp];
      // std::cout << crs[crsp] <<"  "<< sum1 <<"  "<< rnd1 << std::endl;
      if (rnd1 <= sum1) {
        isel = crsp;
        break;
      }
    }
    // std::cout <<"isel "<< isel << std::endl;
    MT  = recoPoint->MtValues[elemId][isel];
    MF4 = recoPoint->GetMt4(elemId, MT);
    MF5 = recoPoint->GetMt5(elemId, MT);
    MF6 = recoPoint->GetMt6(elemId, MT);
    // std::cout <<" MT "<< MT<<"   "<< MF4 <<"  "<< MF5 <<"  "<< MF6 << std::endl;
    // selection of the data from file 4 5 and 6 for angle and energy calculations
    LCT = recoPoint->GetCos4Lct(elemId, MT);
    switch (MT) {
    case 2: { // elastic
      residueA     = particle[elemId].mass;
      residueZ     = particle[elemId].charge;
      cosCM        = recoPoint->GetCos4(elemId, MT, kineticE);
      cosLab       = TNudyCore::Instance()->cmToLabElasticCosT(cosCM, particle[elemId].mass);
      secEnergyLab = TNudyCore::Instance()->cmToLabElasticE(kineticE, cosCM, particle[elemId].mass);
    } break;
    case 11: // 2nd
      residueA = particle[elemId].mass - 3;
      residueZ = particle[elemId].charge - 1;
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 16: // 2n
      residueA = particle[elemId].mass - 1;
      residueZ = particle[elemId].charge;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // std::cout << cosLab <<"  "<<secEnergyLab << std::endl;
      FillHisto(cosLab, secEnergyLab);
      break;
    case 17: // 3n
      residueA = particle[elemId].mass - 2;
      residueZ = particle[elemId].charge;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 18: // fission
    {
      // double kineticRand = kineticE * fRnd->Uniform(1);
      double zaf = recoPoint->GetFisYield(elemId, kineticE);
      // std::cout<<" sample z "<< std::endl;
      divr     = div(zaf, 1000);
      residueZ = divr.quot;
      residueA = divr.rem;
      // fissA->Fill(residueA);
      // fissA->Fill(particle[elemId].mass - residueA);
      // double fissHeat = recoPoint->GetFissHeat(elemId, kineticRand);
      double nut = recoPoint->GetNuTotal(elemId, kineticE);
      // double nup = recoPoint->GetNuPrompt(elemId, kineticRand);
      // double nud = recoPoint->GetNuDelayed(elemId, kineticRand);
      int nu                              = (int)nut;
      if (fRnd->Uniform(1) < nut - nu) nu = nu + 1;
      // std::cout<<" sample heat "<< fissHeat/1E6 <<"  "<< nut <<"  "<< nup <<"  "<< nud <<"  "<< kineticRand <<
      // std::endl;
      // ene[ecounter] = kineticRand/1E6;
      // nu1[ecounter] = nut;
      // nu2[ecounter] = nup;
      // nu3[ecounter] = nud;
      // x[ecounter] = fissHeat/1E6;
      // ecounter ++;
      // hist[0]->Fill(kineticRand/1E6,nut);
      // hist[1]->Fill(kineticRand/1E6,nup);
      // hist[2]->Fill(kineticRand/1E6,nud);
      // hist[3]->Fill(kineticRand/1E6,fissHeat/1E6);
      /*
      int itime = 0;
      int mttime = 1000;
      do
      {
        double lambda = recoPoint->GetLambdaD(elemId, itime);
        double frac = recoPoint->GetDelayedFraction(elemId, mttime, kineticE);
        double deleyedE = recoPoint->GetEnergy5(elemId, mttime, kineticE);
              double frac0 = frac * exp(-lambda * 1);
        fissA1[0]->Fill(deleyedE/1E6, frac0);
              double frac1 = frac * exp(-lambda * 60);
        fissA1[1]->Fill(deleyedE/1E6, frac1);
              double frac2 = frac * exp(-lambda * 100);
        fissA1[2]->Fill(deleyedE/1E6, frac2);
              double frac3 = frac * exp(-lambda * 600);
        fissA1[3]->Fill(deleyedE/1E6, frac3);
              double frac4 = frac * exp(-lambda * 3600);
        fissA1[4]->Fill(deleyedE/1E6, frac4);
              double frac5 = frac * exp(-lambda * 7200);
        fissA1[5]->Fill(deleyedE/1E6, frac5);
        //std::cout<<1/lambda <<"  "<<frac <<"  "<<deleyedE << std::endl;
        mttime++;
        itime++;
      }while(recoPoint->GetLambdaD(elemId, itime) > 0);
      */
      for (int nui = 0; nui < nu; nui++) {
        secEnergyLab = recoPoint->GetEnergy5(elemId, MT, kineticE);
        cosLab       = 2 * fRnd->Uniform(1) - 1;
        fissA1[6]->Fill(secEnergyLab / 1E6);
        // FillHisto(cosLab, secEnergyLab);
      }
    } break;
    case 22: // n+alpha
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 23: // n+3a
      residueA = particle[elemId].mass - 12;
      residueZ = particle[elemId].charge - 6;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 24: // 2n+a
      residueA = particle[elemId].mass - 5;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 25: // 3n+a
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 28: // n+p
      residueA = particle[elemId].mass - 1;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 29: // n +2a
      residueA = particle[elemId].mass - 8;
      residueZ = particle[elemId].charge - 4;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 30: // 2n+2a
      residueA = particle[elemId].mass - 9;
      residueZ = particle[elemId].charge - 4;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 32: // n+d
      residueA = particle[elemId].mass - 2;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 33: // n+t
      residueA = particle[elemId].mass - 3;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 34: // n + He3
      residueA = particle[elemId].mass - 3;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 35: // n+d+2a
      residueA = particle[elemId].mass - 10;
      residueZ = particle[elemId].charge - 5;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 36: // n+t+2a
      residueA = particle[elemId].mass - 11;
      residueZ = particle[elemId].charge - 5;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 37: // 4n
      residueA = particle[elemId].mass - 3;
      residueZ = particle[elemId].charge;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 41: // 2n+p
      residueA = particle[elemId].mass - 2;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 42: // 3n+p
      residueA = particle[elemId].mass - 3;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 44: // n+2p
      residueA = particle[elemId].mass - 2;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 45: // n+p+a
      residueA = particle[elemId].mass - 5;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 50:
    case 51:
    case 52:
    case 53:
    case 54:
    case 55:
    case 56:
    case 57:
    case 58:
    case 59:
    case 60:
    case 61:
    case 62:
    case 63:
    case 64:
    case 65:
    case 66:
    case 67:
    case 68:
    case 69:
    case 70:
    case 71:
    case 72:
    case 73:
    case 74:
    case 75:
    case 76:
    case 77:
    case 78:
    case 79:
    case 80:
    case 81:
    case 82:
    case 83:
    case 84:
    case 85:
    case 86:
    case 87:
    case 88:
    case 89:
    case 90:
    case 91:
      residueA = particle[elemId].mass;
      residueZ = particle[elemId].charge;
      GetSecParameter(particle, recoPoint);
      // std::cout << MT <<"  "<< MF4 <<"  "<< MF5 <<"  "<< MF6 << std::endl;
      // fissA1[0]->Fill(secEnergyLab);
      // cosCM = recoPoint->GetCos4(elemId, MT, kineticE);
      // secEnergyCM = recoPoint->GetEnergy5(elemId, MT, kineticE);
      // secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(secEnergyCM, kineticE, cosCM, particle[elemId].mass);
      // cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
      // particle[elemId].mass);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 102: // capture
      residueA = particle[elemId].mass + 1;
      residueZ = particle[elemId].charge;
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 103: // p
      residueA = particle[elemId].mass;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 104: // d
      residueA = particle[elemId].mass - 1;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 105: // t
      residueA = particle[elemId].mass - 2;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      FillHisto(cosLab, secEnergyLab);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 106: // He3
      residueA = particle[elemId].mass - 2;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 107: // alpha
      residueA = particle[elemId].mass - 3;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 108: // 2a
      residueA = particle[elemId].mass - 7;
      residueZ = particle[elemId].charge - 4;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 109: // 3a
      residueA = particle[elemId].mass - 11;
      residueZ = particle[elemId].charge - 6;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 111: // 2p
      residueA = particle[elemId].mass - 1;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 112: // p+a
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 113: // t+2a
      residueA = particle[elemId].mass - 10;
      residueZ = particle[elemId].charge - 5;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 114: // d+2a
      residueA = particle[elemId].mass - 9;
      residueZ = particle[elemId].charge - 5;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 115: // p+d
      residueA = particle[elemId].mass - 2;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 116: // p+t
      residueA = particle[elemId].mass - 3;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 117: // d+a
      residueA = particle[elemId].mass - 5;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 152: // 5n
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 153: // 6n
      residueA = particle[elemId].mass - 5;
      residueZ = particle[elemId].charge;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 154: // 2n+t
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 155: // t+a
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 156: // 4n+p
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 157: // 3n+d
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 158: // n+d+a
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 159: // 2n+p+a
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 160: // 7n
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 161: // 8n
      residueA = particle[elemId].mass - 7;
      residueZ = particle[elemId].charge;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 162: // 5np
      residueA = particle[elemId].mass - 5;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 163: // 6np
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 164: // 7np
      residueA = particle[elemId].mass - 7;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 165: // 4n+a
      residueA = particle[elemId].mass - 7;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 166: // 5na
      residueA = particle[elemId].mass - 8;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 167: // 6na
      residueA = particle[elemId].mass - 9;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 168: // 7na
      residueA = particle[elemId].mass - 10;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 169: // 4nd
      residueA = particle[elemId].mass - 5;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 170: // 5nd
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 171: // 6nd
      residueA = particle[elemId].mass - 7;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 172: // 3nt
      residueA = particle[elemId].mass - 5;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 173: // 4nt
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 174: // 5nt
      residueA = particle[elemId].mass - 7;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 175: // 6nt
      residueA = particle[elemId].mass - 8;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 176: // 2n+He3
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 177: // 3n + He3
      residueA = particle[elemId].mass - 5;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 178: // 4n +He3
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 179: // 3n2p
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 180: // 3n2a
      residueA = particle[elemId].mass - 10;
      residueZ = particle[elemId].charge - 4;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 181: // 3npa
      residueA = particle[elemId].mass - 7;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 182: // dt
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 183: // npd
      residueA = particle[elemId].mass - 3;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 184: // npt
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 185: // ndt
      residueA = particle[elemId].mass - 5;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 186: // npHe3
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 187: // ndHe3
      residueA = particle[elemId].mass - 5;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 188: // ntHe3
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 189: // nta
      residueA = particle[elemId].mass - 7;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 190: // 2n2p
      residueA = particle[elemId].mass - 3;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 191: // pHe3
      residueA = particle[elemId].mass - 3;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 192: // dHe3
      residueA = particle[elemId].mass - 4;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 193: // aHe3
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge - 4;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 194: // 4n2p
      residueA = particle[elemId].mass - 5;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 195: // 4n2a
      residueA = particle[elemId].mass - 11;
      residueZ = particle[elemId].charge - 4;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 196: // 4npa
      residueA = particle[elemId].mass - 8;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 197: // 3p
      residueA = particle[elemId].mass - 2;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 198: // n3p
      residueA = particle[elemId].mass - 3;
      residueZ = particle[elemId].charge - 3;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 199: // 3n2pa
      residueA = particle[elemId].mass - 8;
      residueZ = particle[elemId].charge - 4;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 200: // 5n2p
      residueA = particle[elemId].mass - 6;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      break;
    }
    // double cosT = recoPoint->GetCos4(elemId, MT, kineticE);
    // double secEnergy = recoPoint->GetEnergy5(elemId, MT, kineticE);
    // std::cout<<"mass = "<< particle[elemId].mass << std::endl;
    // if(MT==2){
    // std::cout <<secEnergyLab/1E9 <<"  "<< cosLab <<  std::endl;
    //}
    crs.clear();
    counter++;
  } while (counter < events);
  fissA->SetLineWidth(3);
  fissA->GetXaxis()->SetTitle("A");
  fissA->GetXaxis()->SetLabelFont(22);
  fissA->GetXaxis()->SetTitleFont(22);
  fissA->GetYaxis()->SetTitle("Yield");
  fissA->GetYaxis()->SetLabelFont(22);
  fissA->GetYaxis()->SetTitleFont(22);
  fissA->GetYaxis()->SetLabelSize(0.035);
  fissA->GetXaxis()->SetLabelSize(0.035);
  h1->SetLineWidth(3);
  h1->GetXaxis()->SetTitle("Cos(\\theta)");
  h1->GetXaxis()->SetLabelFont(22);
  h1->GetXaxis()->SetTitleFont(22);
  h1->GetYaxis()->SetTitle("Counts");
  h1->GetYaxis()->SetLabelFont(22);
  h1->GetYaxis()->SetTitleFont(22);
  h1->GetYaxis()->SetLabelSize(0.035);
  h1->GetXaxis()->SetLabelSize(0.035);
  h1->Draw("P");
  h2->SetLineWidth(3);
  h2->GetXaxis()->SetTitle("Energy, GeV");
  h2->GetXaxis()->SetLabelFont(22);
  h2->GetXaxis()->SetTitleFont(22);
  h2->GetYaxis()->SetTitle("Counts");
  h2->GetYaxis()->SetLabelFont(22);
  h2->GetYaxis()->SetTitleFont(22);
  h2->GetYaxis()->SetLabelSize(0.035);
  h2->GetXaxis()->SetLabelSize(0.035);
  h->SetLineWidth(3);
  h->GetXaxis()->SetTitle("Cos(\\theta)");
  h->GetXaxis()->SetLabelFont(22);
  h->GetXaxis()->SetTitleFont(22);
  h->GetYaxis()->SetTitle("Energy, GeV");
  h->GetYaxis()->SetLabelFont(22);
  h->GetYaxis()->SetTitleFont(22);
  h->GetYaxis()->SetLabelSize(0.035);
  h->GetXaxis()->SetLabelSize(0.035);
  h->Draw("colz");
  h1->Draw("P");
  gr[0] = new TGraph(ecounter, ene, nu1);
  gr[1] = new TGraph(ecounter, ene, nu2);
  gr[2] = new TGraph(ecounter, ene, nu3);
  gr[3] = new TGraph(ecounter, ene, x);
  gr[0]->Draw();
  gr[0]->Write("nut");
  gr[1]->Draw();
  gr[1]->Write("nup");
  gr[2]->Draw();
  gr[2]->Write("nud");
  gr[3]->Draw();
  gr[3]->Write("Fission_Heat");
  /*
  gr1 = new TGraph (ecounter, y , x);
  gr1->GetXaxis()->SetTitle("Angle");
  gr1->GetYaxis()->SetTitle("Energy, GeV");
  gr1->GetXaxis()->SetTitleSize(0.05);
  gr1->GetYaxis()->SetTitleSize(0.05);
  gr1->GetYaxis()->SetTitleOffset(0.5);
  gr1->GetXaxis()->SetTitleOffset(0.85);
  gr1->GetXaxis()->SetTitleSize(0.05);
  gr1->Draw();
  gr1->Write("mygraph");
  */
  f->Write();
  f->Close();
  // std::cout<< std::endl;
}
//------------------------------------------------------------------------------------------------------
void TNudySampling::GetSecParameter(Particle *particle, TNudyEndfRecoPoint *recoPoint)
{
  // std::cout<<"MF "<< MF << " MT "<< MT << " ID "<< elemId << std::endl;
  if (MF4 == 4 && MF5 == 5) {
    cosCM        = recoPoint->GetCos4(elemId, MT, kineticE);
    secEnergyCM  = recoPoint->GetEnergy5(elemId, MT, kineticE);
    secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(secEnergyCM, kineticE, cosCM, particle[elemId].mass);
    if (LCT == 2) {
      cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
                                                           particle[elemId].mass);
    } else if (LCT == 1) {
      cosLab       = cosCM;
      secEnergyLab = secEnergyCM;
    }
  } else if (MF4 == 4 && MF5 == -1) {
    cosCM      = recoPoint->GetCos4(elemId, MT, kineticE);
    double fm1 = 1;
    double fm2 = particle[elemId].mass;
    double fm4 = residueA;
    double fm3 = particle[elemId].mass - residueA;
    // double fqval =  recoPoint->GetQValue(elemId, MT) ;
    double fsi  = fm1 + fm2;
    double fdi  = fm1 - fm2;
    double fsf  = fm3 + fm4;
    double fdf  = fm3 - fm4;
    double fs   = fsi * fsi + 2 * fm2 * kineticE;
    double fki2 = (fs - fsi * fsi) * (fs - fdi * fdi) / (4 * fs);
    double fkf2 = (fs - fsf * fsf) * (fs - fdf * fdf) / (4 * fs);
    double fz2  = sqrt(fki2 + fm2 * fm2);
    double fz3  = sqrt(fkf2 + fm3 * fm3);
    // double fz4 = sqrt(fkf2 + fm4 * fm4);
    double fe3 = (fz2 * fz3 - sqrt(fki2 * fkf2) / fm2) - fm3;
    // double fe4 = (fz2 * fz4 - sqrt(fki2 * fkf2)/fm2) - fm4;
    // double fp12 = kineticE * kineticE + 2 * fm1 * kineticE;
    // double fp32 = fe3 * fe3 + 2 * fm3 * fe3;
    // double fp42 = fe4 * fe4 + 2 * fm4 * fe4;
    // double fcos3 = ((kineticE + fsi)*(fe3 + fm3)- fz3 * sqrt(fs))/sqrt(fp12*fp32);
    // double fcos4 = ((kineticE + fsi)*(fe4 + fm4)- fz4 * sqrt(fs))/sqrt(fp12*fp42);

    // std::cout<<"cos CM = "<< cosCM << std::endl;
    // secEnergyCM = recoPoint->GetEnergy5(elemId, MT, kineticE);
    secEnergyLab = fe3;
    if (LCT == 2) {
      cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
                                                           particle[elemId].mass);
    } else if (LCT == 1) {
      cosLab = cosCM;
    }

  } else if (MF4 == 99 && MF5 == 5) {
    double fm1 = 1;
    double fm2 = particle[elemId].mass;
    double fm4 = residueA;
    double fm3 = particle[elemId].mass - residueA;
    // double fqval =  recoPoint->GetQValue(elemId, MT) ;
    double fsi  = fm1 + fm2;
    double fdi  = fm1 - fm2;
    double fsf  = fm3 + fm4;
    double fdf  = fm3 - fm4;
    double fs   = fsi * fsi + 2 * fm2 * kineticE;
    double fki2 = (fs - fsi * fsi) * (fs - fdi * fdi) / (4 * fs);
    double fkf2 = (fs - fsf * fsf) * (fs - fdf * fdf) / (4 * fs);
    double fz2  = sqrt(fki2 + fm2 * fm2);
    double fz3  = sqrt(fkf2 + fm3 * fm3);
    // double fz4 = sqrt(fkf2 + fm4 * fm4);
    double fe3 = (fz2 * fz3 - sqrt(fki2 * fkf2) / fm2) - fm3;
    // double fe4 = (fz2 * fz4 - sqrt(fki2 * fkf2)/fm2) - fm4;
    double fp12 = kineticE * kineticE + 2 * fm1 * kineticE;
    double fp32 = fe3 * fe3 + 2 * fm3 * fe3;
    // double fp42 = fe4 * fe4 + 2 * fm4 * fe4;
    cosLab = ((kineticE + fsi) * (fe3 + fm3) - fz3 * sqrt(fs)) / sqrt(fp12 * fp32);
    // double fcos4 = ((kineticE + fsi)*(fe4 + fm4)- fz4 * sqrt(fs))/sqrt(fp12*fp42);
    secEnergyLab = fe3;
  } else if (MF4 == 99 && MF5 == -1 && MF6 == 6) {
    int law = recoPoint->GetLaw6(elemId, MT);
    // std::cout<<"law "<< law << std::endl;
    switch (law) {
    case 2: {
      cosCM         = recoPoint->GetCos64(elemId, MT, kineticE);
      double fm1    = 1;
      double fm2    = particle[elemId].mass;
      double fm3    = particle[elemId].mass - residueA;
      double fm4    = residueA;
      double fqval  = recoPoint->GetQValue(elemId, MT);
      double fEt    = kineticE + fqval;
      double fm1m2  = fm1 + fm2;
      double fm3m4  = fm3 + fm4;
      double fA1    = fm1 * fm4 * (kineticE / fEt) / (fm1m2 * fm3m4);
      double fB1    = fm1 * fm3 * (kineticE / fEt) / (fm1m2 * fm3m4);
      double fC1    = fm2 * fm3 * (1 + fm1 * fqval * fm4 / (fEt * fm1)) / (fm1m2 * fm3m4);
      double fD1    = fm2 * fm4 * (1 + fm1 * fqval * fm4 / (fEt * fm1)) / (fm1m2 * fm3m4);
      secEnergyLab  = fEt * (fB1 + fD1 + 2 * sqrt(fA1 * fC1) * cosCM);
      double sinLab = fD1 * cosCM * fEt / secEnergyLab;
      cosLab        = sqrt(1 - sinLab * sinLab);
      // std::cout<< cosLab <<"  "<<secEnergyLab << std::endl;
    } break;
    case 3:
      break;
    case 1:
    case 6:
      cosCM = recoPoint->GetCos6(elemId, MT, kineticE);
      // std::cout<<"cos "<< cosCM << std::endl;
      secEnergyCM = recoPoint->GetEnergy6(elemId, MT, kineticE);
      // std::cout<<"E "<< secEnergyCM<< std::endl;
      secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(secEnergyCM, kineticE, cosCM, particle[elemId].mass);
      cosLab       = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
                                                           particle[elemId].mass);
      break;
    case 7:
      cosLab = recoPoint->GetCos6(elemId, MT, kineticE);
      // std::cout<< cosLab << std::endl;
      secEnergyLab = recoPoint->GetEnergy6(elemId, MT, kineticE);
      // std::cout<< secEnergyLab << std::endl;
      break;
    }
  }
}
//------------------------------------------------------------------------------------------------------
void TNudySampling::FillHisto(double icosLab, double isecEnergyLab)
{
  h->Fill(icosLab, isecEnergyLab / 1E9);
  h1->Fill(icosLab);
  h2->Fill(isecEnergyLab / 1E9);
  // x[ecounter] = isecEnergyLab/1E9 ;
  // y[ecounter] = icosLab ;
  // if(events<=1000000)ecounter ++ ;
}
//__________________________________________________________________________________________________________________
TNudySampling::~TNudySampling()
{
}
