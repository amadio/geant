// 	This class is reconstructing Fission yields
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 24, 2016

#include "TList.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfPhProd.h"
using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfPhProd)
#include "TRandom3.h"
#endif

    TNudyEndfPhProd::TNudyEndfPhProd()
{
}

//______________________________________________________________________________
TNudyEndfPhProd::TNudyEndfPhProd(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    // double ZA   = sec->GetC1();
    // double AWR  = sec->GetC2();
    // div_t divr;
    // int MT = sec->GetMT();
    // int LE = sec->GetL1();
    int NK = sec->GetN1();
    //       std::cout << " LO = " << LO <<" NK "<< NK << std::endl;
    TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
    if (NK == 1) {
      egk1.push_back(tab1->GetC1()); // photon energy for LP=0 or 1 or Binding Energy for LP=2. For a continuous photon
                                     // energy distribution, EGk ≡ 0.0 should be used.
      esk1.push_back(
          tab1->GetC2()); // energy of the level from which the photon originates, 0 for unknown and continuum
      lpk1.push_back(tab1->GetL1()); // photon energy if == 0,1 Egp = eg +awr/(1+awr) if ==2
      lfk1.push_back(tab1->GetL2()); // tabulated in file 15 if ==1, discrete if == 2
      nrk1.push_back(tab1->GetN1()); // regions
      npk1.push_back(tab1->GetN2()); // points
      for (int cr = 0; cr < tab1->GetN1(); cr++) {
        nbt1.push_back(tab1->GetNBT(cr));
        int1.push_back(tab1->GetINT(cr));
      }
      for (int crs = 0; crs < tab1->GetN2(); crs++) {
        eintk1.push_back(tab1->GetX(crs));
        sigPh1.push_back(tab1->GetY(crs));
        // 	      std::cout << " neutron energy =  "<<  tab1->GetX(crs)<<"  "<< tab1->GetY(crs) << std::endl;
      }
    } else if (NK > 1) {
      egk1.push_back(tab1->GetC1()); // photon energy for LP=0 or 1 or Binding Energy for LP=2. For a continuous photon
                                     // energy distribution, EGk ≡ 0.0 should be used.
      esk1.push_back(
          tab1->GetC2()); // energy of the level from which the photon originates, 0 for unknown and continuum
      lpk1.push_back(tab1->GetL1()); // photon energy if == 0,1 Egp = eg +awr/(1+awr) if ==2
      lfk1.push_back(tab1->GetL2()); // tabulated in file 15 if ==1, discrete if == 2
      nrk1.push_back(tab1->GetN1()); // regions
      npk1.push_back(tab1->GetN2()); // points
      for (int cr = 0; cr < tab1->GetN1(); cr++) {
        nbt1.push_back(tab1->GetNBT(cr));
        int1.push_back(tab1->GetINT(cr));
      }
      for (int crs = 0; crs < tab1->GetN2(); crs++) {
        eintk1.push_back(tab1->GetX(crs));
        sigPh1.push_back(tab1->GetY(crs));
        //  	      std::cout << " neutron energy =  "<<  tab1->GetX(crs)<<"  "<< tab1->GetY(crs) << std::endl;
      }
      for (int nki = 0; nki < NK; nki++) {
        TNudyEndfTab1 *tab2 = (TNudyEndfTab1 *)recIter.Next();
        esk1.push_back(tab2->GetC1());
        egk1.push_back(tab2->GetC2());
        lpk1.push_back(tab2->GetL1());
        lfk1.push_back(tab2->GetL2());
        nrk1.push_back(tab2->GetN1());
        npk1.push_back(tab2->GetN2());
        for (int cr = 0; cr < tab2->GetN1(); cr++) {
          nbt1.push_back(tab2->GetNBT(cr));
          int1.push_back(tab2->GetINT(cr));
        }
        for (int crs = 0; crs < tab2->GetN2(); crs++) {
          eintk1.push_back(tab2->GetX(crs));
          sigPh1.push_back(tab2->GetY(crs));
          // if (lpk1[lpk1.size() - 1] == 2) egk1[ egk1.size() - 1 ] = egk1[ egk1.size() - 1 ] + AWR * tab1->GetX(crs)/
          // (1 + AWR );
          // std::cout << " photon energy =  "<< tab1->GetC1() << "  " << tab1->GetX(crs)<<"  "<< tab1->GetY(crs) <<"
          // "<< egk1[ egk1.size() - 1 ] << std::endl;
          //  	      std::cout << " neutron energy =  "<<  tab2->GetX(crs)<<"  "<< tab2->GetY(crs) << std::endl;
        }
      }
    }
  }
}

TNudyEndfPhProd::~TNudyEndfPhProd()
{
}
/*
double TNudyEndfPhProd::GetFisYield(int ielemId, double energyK)
{
  fRnd = new TRandom3(0);
  // std::cout<<"element "<< fEinfId.size() << std::endl;
  // std::cout<<"energies "<< fEinfId[ielemId].size() << std::endl;
  return fRnd;
}
*/
