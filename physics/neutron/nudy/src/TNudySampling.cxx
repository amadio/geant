// 	Selection of secondary particle energy, angle
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: June 22, 2016
#include <iostream>
#include "Geant/TNudyCore.h"
#include "Geant/TNudyEndfRecoPoint.h"
#include "Geant/TNudySampling.h"
#include "TCanvas.h"
#include "TFile.h"

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudySampling)
#endif

TNudySampling::TNudySampling()
{
}
//------------------------------------------------------------------------------------------------------
TNudySampling::TNudySampling(TParticleTest *particle, TNudyEndfRecoPoint *recoPoint)
{
  kineticE = particle[elemId].energy;
  std::cout << "sampling start \t" << elemId << "\t" << kineticE << std::endl;
  events = 1000000;
  div_t divr;
  int jkl  = 0;
  TFile *f = new TFile("test.root", "recreate");
  f->cd();
  h     = new TH2D("h", "", 100, -1, 1, 100, 0, 20E6);
  h1    = new TH1D("h1", "", 101, -1, 1);
  h2    = new TH1D("h2", "", 100, 0, 20E6);
  fissA = new TH1D("fissA", "", 243, 0, 242);
  for (int i = 0; i < 21; i++) {
    fissA1[i] = new TH1D(Form("fissA1[%d]", i), "", 100, 0, 20E6);
    hist[i]   = new TH2D(Form("hist[%d]", i), "", 1000, 0, 20, 10, 0, 10);
  }
  std::cout << "sampling sigmaTotal " << recoPoint->GetSigmaTotal(elemId, kineticE) << std::endl;
  // std::cout <<"sampling sigmaPartial elstic "<< recoPoint->GetSigmaPartial(0,2,20) << std::endl;
  // determining reaction type from element;
  // double density = 1;
  // double charge = 1;
  // double avg = 6.022E23;
  // double ro = avg * density / mass;
  // int enemax = recoPoint->eneUni[elemId].size();
  std::cout << " target mass " << particle[elemId].mass << "  " << particle[elemId].charge << "  " << kineticE
            << std::endl;
  // kineticE = fRng.uniform() * recoPoint->eneUni[elemId][enemax-1];
  for (unsigned int crsp = 0; crsp < recoPoint->fMtValues[elemId].size(); crsp++) {
    int mtValue = recoPoint->fMtValues[elemId][crsp];
    std::cout << recoPoint->fMtValues[elemId][crsp] << "  \t" << recoPoint->GetSigmaTotal(elemId, kineticE) << "  \t"
    << recoPoint->GetSigmaPartial(elemId, mtValue, kineticE) << std::endl;
    crs.push_back(recoPoint->GetSigmaPartial(elemId, mtValue, kineticE) / recoPoint->GetSigmaTotal(elemId, kineticE));
    //     std::cout<<recoPoint->fMtValues[elemId][crsp]<<"  "<< crs[crsp] <<"  "<<
    //     recoPoint->GetSigmaPartial(elemId,crsp,kineticE) <<"  "<< recoPoint->GetSigmaTotal(elemId,kineticE) <<
    //     std::endl;
  }
  //  exit(1);
  do {
    kineticE    = particle[elemId].energy;
    double sum1 = 0;
    // std::cout << counter << " bef " << recoPoint->fMtValues[elemId].size() << std::endl;
    double rnd1 = fRng.uniform();
    // std::cout << counter <<" aft "<<rnd1<< std::endl;
    for (unsigned int crsp = 0; crsp < recoPoint->fMtValues[elemId].size(); crsp++) {
      sum1 += crs[crsp];
      // std::cout << crs[crsp] <<"  "<< sum1 <<"  "<< rnd1 << std::endl;
      if (rnd1 <= sum1) {
        isel = crsp;
        break;
      }
    }
    // std::cout <<"isel "<< isel << std::endl;
    MT  = recoPoint->fMtValues[elemId][isel];
    MF4 = recoPoint->GetMt4(elemId, MT);
    MF5 = recoPoint->GetMt5(elemId, MT);
    MF6 = recoPoint->GetMt6(elemId, MT);
    // int MT6N = recoPoint->GetMt6Neutron(elemId, MT);
    // std::cout << " bef "<< counter << " MT " << MT<<"   "<< MF4 <<"  "<< MF5 <<"  "<< MF6 <<"  "<< kineticE <<
    // std::endl;
    // selection of the data from file 4 5 and 6 for angle and energy calculations
    LCT = recoPoint->GetCos4Lct(elemId, MT);
    switch (MT) {
    case 2: { // elastic
      residueA = particle[elemId].mass;
      residueZ = particle[elemId].charge;
      cosCM    = recoPoint->GetCos4(elemId, MT, kineticE);
      //      std::cout << cosCM <<"\t"<< LCT << std::endl;
      cosLab       = TNudyCore::Instance()->CmToLabElasticCosT(cosCM, particle[elemId].mass);
      secEnergyLab = TNudyCore::Instance()->CmToLabElasticE(kineticE, cosCM, particle[elemId].mass);
      //      FillHisto(acos(cosCM)*180/3.14159, secEnergyLab);
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
      // FillHisto(cosLab, secEnergyLab);
      jkl = (int)((1 - cosLab) * 10);
      fissA1[jkl]->Fill(secEnergyLab);
      fissA1[20]->Fill(secEnergyLab);
      // std::cout << jkl <<"  "<< cosLab <<"  "<<secEnergyLab << std::endl;
      //       kineticE = kineticE - secEnergyLab ;
      //       GetSecParameter(particle, recoPoint);
      //       FillHisto(cosLab, secEnergyLab);
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
      // double kineticRand = kineticE * fRng.uniform();
      double zaf = recoPoint->GetFisYield(elemId, kineticE);
      // std::cout<<" sample z "<< std::endl;
      divr     = div(zaf, 10000);
      residueZ = divr.quot;
      divr     = div(divr.rem, 10);
      residueA = divr.quot;
      // std::cout <<"mass "<< residueA <<"  "<< residueZ <<"  "<< zaf << std::endl;
      fissA->Fill(residueA);
      fissA->Fill(particle[elemId].mass - residueA);
      // double fissHeat = recoPoint->GetFissHeat(elemId, kineticRand);
      double nut = recoPoint->GetNuTotal(elemId, kineticE);
      // double nup = recoPoint->GetNuPrompt(elemId, kineticRand);
      // double nud = recoPoint->GetNuDelayed(elemId, kineticRand);
      int nu                              = (int)nut;
      if (fRng.uniform() < nut - nu) nu = nu + 1;
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
        cosLab       = 2 * fRng.uniform() - 1;
        // fissA1[6]->Fill(secEnergyLab / 1E6);
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
      // secEnergyLab = TNudyCore::Instance()->CmToLabInelasticE(secEnergyCM, kineticE, cosCM, particle[elemId].mass);
      // cosLab = TNudyCore::Instance()->CmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
      // particle[elemId].mass);
      if (MT == 91) FillHisto(cosLab, secEnergyLab);
      break;
    case 102: // capture
      residueA = particle[elemId].mass + 1;
      residueZ = particle[elemId].charge;
      fissA->Fill(residueA);
      GetSecParameter(particle, recoPoint);
      std::cout << cosLab << "\t" << secEnergyLab << "\t" << std::endl;
      //      FillHisto(cosLab, secEnergyLab);
      //      FillHisto(acos(cosCM)*180/3.14159, secEnergyLab);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 103: // p
    case 600:
    case 601:
    case 602:
    case 603:
    case 604:
    case 605:
    case 606:
    case 607:
    case 608:
    case 609:
    case 610:
    case 611:
    case 612:
    case 613:
    case 614:
    case 615:
    case 616:
    case 617:
    case 618:
    case 619:
    case 620:
    case 621:
    case 622:
    case 623:
    case 624:
    case 625:
    case 626:
    case 627:
    case 628:
    case 629:
    case 630:
    case 631:
    case 632:
    case 633:
    case 634:
    case 635:
    case 636:
    case 637:
    case 638:
    case 639:
    case 640:
    case 641:
    case 642:
    case 643:
    case 644:
    case 645:
    case 646:
    case 647:
    case 648:
    case 649:
      residueA = particle[elemId].mass;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 104: // d
    case 650:
    case 651:
    case 652:
    case 653:
    case 654:
    case 655:
    case 656:
    case 657:
    case 658:
    case 659:
    case 660:
    case 661:
    case 662:
    case 663:
    case 664:
    case 665:
    case 666:
    case 667:
    case 668:
    case 669:
    case 670:
    case 671:
    case 672:
    case 673:
    case 674:
    case 675:
    case 676:
    case 677:
    case 678:
    case 679:
    case 680:
    case 681:
    case 682:
    case 683:
    case 684:
    case 685:
    case 686:
    case 687:
    case 688:
    case 689:
    case 690:
    case 691:
    case 692:
    case 693:
    case 694:
    case 695:
    case 696:
    case 697:
    case 698:
    case 699:
      residueA = particle[elemId].mass - 1;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 105: // t
    case 700:
    case 701:
    case 702:
    case 703:
    case 704:
    case 705:
    case 706:
    case 707:
    case 708:
    case 709:
    case 710:
    case 711:
    case 712:
    case 713:
    case 714:
    case 715:
    case 716:
    case 717:
    case 718:
    case 719:
    case 720:
    case 721:
    case 722:
    case 723:
    case 724:
    case 725:
    case 726:
    case 727:
    case 728:
    case 729:
    case 730:
    case 731:
    case 732:
    case 733:
    case 734:
    case 735:
    case 736:
    case 737:
    case 738:
    case 739:
    case 740:
    case 741:
    case 742:
    case 743:
    case 744:
    case 745:
    case 746:
    case 747:
    case 748:
    case 749:
      residueA = particle[elemId].mass - 2;
      residueZ = particle[elemId].charge - 1;
      GetSecParameter(particle, recoPoint);
      // FillHisto(cosLab, secEnergyLab);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 106: // He3
    case 750:
    case 751:
    case 752:
    case 753:
    case 754:
    case 755:
    case 756:
    case 757:
    case 758:
    case 759:
    case 760:
    case 761:
    case 762:
    case 763:
    case 764:
    case 765:
    case 766:
    case 767:
    case 768:
    case 769:
    case 770:
    case 771:
    case 772:
    case 773:
    case 774:
    case 775:
    case 776:
    case 777:
    case 778:
    case 779:
    case 780:
    case 781:
    case 782:
    case 783:
    case 784:
    case 785:
    case 786:
    case 787:
    case 788:
    case 789:
    case 790:
    case 791:
    case 792:
    case 793:
    case 794:
    case 795:
    case 796:
    case 797:
    case 798:
    case 799:
      residueA = particle[elemId].mass - 2;
      residueZ = particle[elemId].charge - 2;
      GetSecParameter(particle, recoPoint);
      fissA->Fill(residueA);
      // FillHisto(cosLab, secEnergyLab);
      break;
    case 107: // alpha
    case 800:
    case 801:
    case 802:
    case 803:
    case 804:
    case 805:
    case 806:
    case 807:
    case 808:
    case 809:
    case 810:
    case 811:
    case 812:
    case 813:
    case 814:
    case 815:
    case 816:
    case 817:
    case 818:
    case 819:
    case 820:
    case 821:
    case 822:
    case 823:
    case 824:
    case 825:
    case 826:
    case 827:
    case 828:
    case 829:
    case 830:
    case 831:
    case 832:
    case 833:
    case 834:
    case 835:
    case 836:
    case 837:
    case 838:
    case 839:
    case 840:
    case 841:
    case 842:
    case 843:
    case 844:
    case 845:
    case 846:
    case 847:
    case 848:
    case 849:
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
void TNudySampling::GetSecParameter(TParticleTest *particle, TNudyEndfRecoPoint *recoPoint)
{
  std::cout << "MF " << MF4 << " MF5 " << MF5 << " MF6 " << MF6 << " MT " << MT << " LCT " << LCT << std::endl;
  if (MF4 == 4 && MF5 == 5) {
    cosCM        = recoPoint->GetCos4(elemId, MT, kineticE);
    secEnergyCM  = recoPoint->GetEnergy5(elemId, MT, kineticE);
    secEnergyLab = TNudyCore::Instance()->CmToLabInelasticE(secEnergyCM, kineticE, cosCM, particle[elemId].mass);
    if (LCT == 2) {
      cosLab = TNudyCore::Instance()->CmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
                                                           particle[elemId].mass);
    } else if (LCT == 1) {
      cosLab       = cosCM;
      secEnergyLab = secEnergyCM;
    }
    std::cout << "cosLab " << cosLab << " secEnergyLab " << secEnergyLab << std::endl;
  } else if (MF4 == 4 && MF5 == -1) {
    cosCM      = recoPoint->GetCos4(elemId, MT, kineticE);
    double fm1 = 1E6;
    double fm2 = particle[elemId].mass * 1E6;
    double fm4 = residueA * 1E6;
    double fm3 = particle[elemId].mass * 1E6 - residueA * 1E6;
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

    // KinematicNonRel (particle, recoPoint) ;
    secEnergyLab = fe3;
    // std::cout<<"MT \t"<< MT <<"  \t" << cosLab  <<"   \t"<< cosCM  <<"   \t"<<  secEnergyLab << std::endl;
    // std::cout<<"cos CM = "<< cosCM << std::endl;
    // secEnergyCM = recoPoint->GetEnergy5(elemId, MT, kineticE);
    if (LCT == 2) {
      cosLab = TNudyCore::Instance()->CmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
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
    // int zi 	= recoPoint->GetZd6(elemId, MT);
    // int ai 	= recoPoint->GetAd6(elemId, MT);
    // std::cout<<"law "<< law <<" Zi \t" << zi <<" Ai \t" << ai << std::endl;
    switch (law) {
    case 2: {
//      cosCM      = recoPoint->GetCos64(elemId, MT, kineticE);
      double fm1 = 1;
      // double fm2    = particle[elemId].mass;
      double fm3   = particle[elemId].mass - residueA;
      double fm4   = residueA;
      double fqval = recoPoint->GetQValue(elemId, MT);
      double fg    = sqrt((fm1 * fm4 * kineticE) / (fm3 * (fm3 + fm4) * fqval + fm3 * (fm3 + fm4 - fm1) * kineticE));
      cosLab       = (fg + cosCM) / sqrt(1 + fg * fg + 2 * fg * cosCM);
      double fr    = sqrt(fm1 * fm4 * kineticE) * cosLab / (fm3 + fm4);
      double fs    = (kineticE * (fm3 - fm1) + fm3 * fqval) / (fm3 + fm4);
      double feb   = fr + sqrt(fr * fr + fs);
      secEnergyLab = feb * feb;

      //       std::cout<< cosLab  <<"   \t"<< cosCM  <<"   \t"<<  secEnergyLab << std::endl;
      //       KinematicNonRel (particle, recoPoint) ;
      //       std::cout<< cosLab  <<"   \t"<< cosCM  <<"   \t"<<  secEnergyLab << std::endl;
      //      KinematicRel (particle, recoPoint, cosLab, secEnergyLab) ;
      if (fm4 == 0) KinematicGama(particle, recoPoint);
      //       std::cout<< cosLab  <<"   \t"<< cosCM  <<"   \t"<<  secEnergyLab << std::endl;
    } break;
    case 3:
      break;
    case 1:
    case 6:
      cosCM = recoPoint->GetCos6(elemId, MT, kineticE);
      // std::cout<<"cos "<< cosCM << std::endl;
      secEnergyCM = recoPoint->GetEnergy6(elemId, MT, kineticE);
      // std::cout<<"E "<< secEnergyCM<< std::endl;
      secEnergyLab = TNudyCore::Instance()->CmToLabInelasticE(secEnergyCM, kineticE, cosCM, particle[elemId].mass);
      cosLab       = TNudyCore::Instance()->CmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
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
double TNudySampling::KinematicNonRel(TParticleTest * /*particle*/, TNudyEndfRecoPoint * /*recoPoint*/)
{
  /*
    double a1 	=  particle[elemId].mass ;
    double a2 	=  particle[elemId].mass - residueA ;
    double fqval  = recoPoint->GetQValue(elemId, MT);
    double cosCM 	= recoPoint->GetCos64(elemId, MT, kineticE);
    double beta 	= sqrt((a1* ( a1 + 1- a2 )/a2) * (1 + (1+a1)*fqval/(a1 * kineticE))) ;
    // double gama 	= (a2 * beta)/ (a1 + 1- a2);
    // double e3bye1 = a2 * beta * beta / ( a1 * a1 ) ;
    // double e1 	= (a1/(1+a1)) * (a1/(1+a1)) * kineticE ;
    // double mu3 	= cosCM ;
    // double e4bye1 = (a2 * e3bye1)/ (a1 + 1- a2) ;
    // double mu4 	= -cosCM ;
    cosLab	= (1 + beta *cosCM )/sqrt(beta * beta + 1 + 2 * beta * cosCM ) ;
    // double w4 	= (1 - gama * cosCM) / sqrt(1 + gama * gama  - 2 * gama * cosCM ) ;
    secEnergyLab 	= kineticE *(a2/((1 + a1)*(1 + a1)) * ( beta * beta + 1 + 2 * beta * cosCM )) ;
    double EnergyR= kineticE * ( a1 + 1- a2 ) * (1 + gama * gama  - 2 * gama * cosCM )/((1 + a1)*(1 + a1)) ;
  */
  return 0;
}
//------------------------------------------------------------------------------------------------------
double TNudySampling::KinematicRel(TParticleTest * /*particle*/, TNudyEndfRecoPoint * /*recoPoint*/)
{
  /*
    double a1 	=  1E6 ;
    double a2 	=  particle[elemId].mass*1E6 ;
    double a3 	=  residueA*1E6 ;
    double a4 	=  recoPoint->GetQValue(elemId, MT);
    double si	=  a1 + a2 ;
    double di	=  a1 - a2 ;
    double sf	=  a3 + a4 ;
    double df	=  a3 - a4 ;
    // double fqval  = recoPoint->GetQValue(elemId, MT);
    double cosCM 	= recoPoint->GetCos64(elemId, MT, kineticE);
    double s	= si * si + 2 * a2 * kineticE ;
  //   std::cout<<"q value \t" << fqval <<" cosCM \t"<< cosCM <<" s \t"<< s << std::endl;
    double ki2 	= (s - si * si) * (s - di * di)/(4 * s) ;
    double kf2 	= (s - sf * sf) * (s - df * df)/(4 * s) ;
    double z2 	= sqrt (ki2 + a2 * a2) ;
    double z3 	= sqrt (kf2 + a3 * a3) ;
    //double z4 	= sqrt (kf2 + a4 * a4) ;
  //   std::cout<<"ki2 \t" << ki2 <<" kf2 \t"<< kf2 <<" z2 \t"<< z2 <<" z3 \t"<< z3 <<" z4 \t"<< z4 << std::endl;
    secEnergyLab	= (z2 * z3 - sqrt(ki2 * kf2 * cosCM))/(a2) - a3;
    // double EnergyR= (z2 * z4 - sqrt(ki2 * kf2 * cosCM))/(a2) - a4;
  //   std::cout<<"secEnergyLab \t" << secEnergyLab <<" z2 * z3 \t"<< z2 * z3 <<" sqrt(ki2 * kf2 * cosCM) \t"<< sqrt(ki2
  * kf2 * cosCM) <<" a3 \t"<< a3 <<" z4 \t"<< z4 << std::endl;
    double p12	= kineticE * kineticE + 2 * a1 * kineticE ;
    double p32	= secEnergyLab * secEnergyLab + 2 * a3 * secEnergyLab ;
    // double p42	= EnergyR * EnergyR + 2 * a4 * EnergyR ;
    cosLab	= (kineticE + si) * (secEnergyLab + a3) -z3 * sqrt (s) / (p12 * p32);
    double w4 	= (kineticE + si) * (EnergyR + a4) -z4 * sqrt (s) / (p12 * p42);
  */
  return 0;
}
//------------------------------------------------------------------------------------------------------
double TNudySampling::KinematicGama(TParticleTest * /*particle*/, TNudyEndfRecoPoint * /*recoPoint*/)
{
  /*
    double a1 	=  1E6 ;
    double a2 	=  particle[elemId].mass*1E6 ;
    double a3 	=  residueA*1E6 ;
    // double a4 	=  recoPoint->GetQValue(elemId, MT);
    double fqval  =  recoPoint->GetQValue(elemId, MT);
    double cosCM 	=  recoPoint->GetCos64(elemId, MT, kineticE);
    double beta	=  sqrt(kineticE *(kineticE + 2*a1))/(a1 + a2 + kineticE);
    cosLab	=  (cosCM + beta)/(1 + beta * cosCM) ;
    secEnergyLab	=  (fqval * (a1 + a2 + a3)/2 + a2 * kineticE )/(a1 + a2 + kineticE - cosLab * sqrt(kineticE *(kineticE
    + 2 * a1)));
    //std::cout<<"secEnergyLab \t" << secEnergyLab <<" coslab \t"<< cosLab << std::endl;
  */
  return 0;
}
//------------------------------------------------------------------------------------------------------
void TNudySampling::FillHisto(double icosLab, double isecEnergyLab)
{
  // std::cout << icosLab <<"  "<< isecEnergyLab << std::endl;
  h->Fill(icosLab, isecEnergyLab);
  h1->Fill(icosLab);
  h2->Fill(isecEnergyLab);
  // x[ecounter] = isecEnergyLab/1E9 ;
  // y[ecounter] = icosLab ;
  // if(events<=1000000)ecounter ++ ;
}
//__________________________________________________________________________________________________________________
TNudySampling::~TNudySampling()
{
}
