// 	This class is reconstructing probability tables for Energy distribution 
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 22, 2016

#include "TNudyEndfEnergyAng.h"
class TNudyEndfFile;
class TNudyEndfSec;
class TNudyEndfCont;
class TNudyEndfTab1;
class TNudyEndfTab2;
class TNudyEndfList;

TNudyEndfEnergyAng::TNudyEndfEnergyAng(){}

//______________________________________________________________________________
TNudyEndfEnergyAng::TNudyEndfEnergyAng(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  std::vector<double>ein, cdf, pdf;
  std::vector<std::vector<double> >cdf2d, pdf2d;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    int NK = sec->GetN1();
    for (int k = 0; k < NK; k++) {
      TIter recIter(sec->GetRecords());
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      div_t divr; 
      AWRI = sec->GetC2();
      int MT = sec->GetMT();
      MtNumbers.push_back(MT);
      //      int LCT = sec->GetL2();
      divr = div(sec->GetC1(),1000);
      double ZA = divr.quot;
      double AA = divr.rem;
      double NA1 = AA -ZA;
      int ZAP =tab1->GetC1(); 
      double AWP =tab1->GetC2(); 
      int LIP =tab1->GetL1(); 
      int LAW = tab1->GetL2();
      NR = tab1->GetN1();
      NP = tab1->GetN2();
      divr =div(ZAP,1000);
      double particleA = divr.rem;
      double particleZ = divr.quot;
      std::cout<<"ZAP = "<< ZAP <<" AWP "<<AWP <<" LIP "<< LIP <<" parZ "<< particleZ <<" parA "<< particleA << std::endl;
      std::cout<<"LAW = "<< LAW <<" MT "<<MT <<" NR "<< NR <<" NP "<< NP << std::endl;
      if(LAW == 1){
	int ND,NA,NW,NEP,LANG,LEP;
	NBT1         = new int[NR]();
	INT1        = new int[NR]();
	fE1_file6    = new double[NP]();
	fp11_file6   = new double[NP]();
	for(int cr=0; cr < NR ; cr ++){
	  NBT1[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1_file6[crs]  = tab1->GetX(crs);
	  fp11_file6[crs] = tab1->GetY(crs);
//    std::cout<<"energy "<< fE1_file6[crs] <<" multip "<< fp11_file6[crs] << std::endl;
	}
	TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
        LANG  = tab2->GetL1();
        LEP  = tab2->GetL2();
        NR2  = tab2->GetN1();
        NE2  = tab2->GetN2();
        NBT2  = new int[NR2]();
        INT2  = new int[NR2]();
	for(int cr=0; cr < NR2 ; cr ++){
	  NBT2[cr] = tab2->GetNBT(cr);
	  INT2[cr] = tab2->GetINT(cr);
	}
//      std::cout<<"LANG "<< LANG <<" LEP "<< LEP <<" NE2 "<< NE2 << std::endl;
	for (int lis = 0; lis < NE2; lis++){
          TNudyEndfList *header = (TNudyEndfList *)recIter.Next();
	  ein.push_back(header->GetC2());
	  ND   = header->GetL1();
	  NA   = header->GetL2();
	  NW   = header->GetN1();
	  NEP  = header->GetN2();
          std::cout <<lis<<"  "<<ein[lis] <<" ND "<< ND <<" NA "<< NA <<" NEP "<< NEP << std::endl;
          if(LANG == 2 && NA==1){
	    TArrayD edes6(NEP),f06(NEP),r6(NEP);
	    for (int lis1 = 0; lis1 < NEP; lis1++){
	      edes6[lis1] = header->GetLIST(lis1*3 + 0);
	      f06[lis1] = header->GetLIST(lis1*3 + 1);
	      r6[lis1] = header->GetLIST(lis1*3 + 2);
 	      //std::cout<<lis1<<"  "<<edes6[lis1] << std::endl;
	    }
	    double epsa = ein[lis] *AWRI/(1. + AWRI);
	    double AC = 1. + AA, ZC = ZA , NC = AC - ZC;
	    double AB = AC - particleA, ZB = ZA - particleZ, NB = AB - ZB;
	    double Ia = 0.0;
	    double Ib = 0.0;
	    double Sa = 15.68*(AC - AA) - 28.07*((NC-ZC)*(NC-ZC)/AC - (NA1-ZA)*(NA1-ZA)/AA)-
			  -18.56*(pow(AC,2./3.)-pow(AA,2./3.)) + 33.22*((NC-ZC)*(NC-ZC)/pow(AC,4./3) - (NA1-ZA)*(NA1-ZA)/pow(AA,4./3))
			  -0.717*(ZC*ZC/pow(AC,1./3.)-ZA*ZA/pow(AA,1./3.)) + 1.211*(ZC*ZC/AC-ZA*ZA/AA) - Ia; 
	    double Sb = 15.68*(AC - AB) - 28.07*((NC-ZC)*(NC-ZC)/AC - (NB-ZB)*(NB-ZB)/AB)-
			  -18.56*(pow(AC,2./3.)-pow(AB,2./3.)) + 33.22*((NC-ZC)*(NC-ZC)/pow(AC,4./3) - (NB-ZB)*(NB-ZB)/pow(AB,4./3))
			  -0.717*(ZC*ZC/pow(AC,1./3.)-ZB*ZB/pow(AB,1./3.)) + 1.211*(ZC*ZC/AC-ZB*ZB/AB) - Ib; 
	    double C1 = 0.04, C2 = 1.8E-6, C3=6.7E-7, Et1=130.0, Et3 = 41.0;
	    double Mn1 = 1;//, Mp = 1, Md = 1,Mt = 1, Mhe3 = 1, Malpha = 0;
	    double mn = 0.5;//, mp = 1, md = 1, mt = 1, mhe3 = 1, malpha = 2;

	    int k1 = 0;
	    do
	    {
	      double x = -1. + k1*0.05;
	      double eoutlab, xlab, sum = 0; 
	      for (int j = 0; j < edes6.GetSize(); j++) {
		double epsb = edes6[j] *(particleA + AWRI)/(AWRI+1. - AWP);
		double ea = epsa + Sa, eb = epsb + Sb;
		double R1 = std::min(ea,Et1), R3 = std::min(eb,Et3);
		double X1 = R1*eb/ea, X3 = R3*eb/ea; 
		double a0 = C1 * X1 + C2 * X1*X1*X1 + C3 * Mn1 * mn * pow(X3,4);
		double prob = (a0*f06[j]/2/sinh(a0))*(cosh(a0*x) + r6[j]*sinh(a0*x));
		eoutlab = edes6[j] + ein[lis]*AWP/((1.+ AWRI)*(1.+ AWRI)) + 2.*sqrt(edes6[j]*ein[lis])*x/(1.+ AWRI);
		xlab = sqrt(edes6[j]/eoutlab)*x + sqrt(ein[lis]/eoutlab)*sqrt(AWP)/(1.+ AWRI);
		if(xlab <= 1 && xlab >= -1){
		  energyFile5.push_back(eoutlab);
		  energyPdfFile5.push_back(prob);
		  cosFile4.push_back(xlab);
		  cosPdfFile4.push_back(prob);
		  sum += prob;
//	          std::cout<< j <<"  "<< edes6[j]<<"  "<< eoutlab <<"  "<< prob <<"  "<< x <<"  "<< acos(xlab)*180/3.14159 << std::endl;
		}
	      }
	      k1++;
	    }while(k1<41);
	    //std::cout<< edes6[j] <<"  "<< sum << std::endl;
	    /*
	    for(int l=0; l<20; l++){
	      recursionLinearFile4(lis, cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
	    }
	    Sort(cosFile4, cosPdfFile4);
	    cdfGenerateT(cosFile4, cosPdfFile4);
	    for(int i = 0; i < cosFile4.size(); i++){
	      pdf.push_back(cosFile4[i]);
	      pdf.push_back(cosPdfFile4[i]);
	      cdf.push_back(cosFile4[i]);
	      cdf.push_back(cosCdfFile4[i]);	  
	    }
	    pdf2d.push_back(pdf);
	    cdf2d.push_back(cdf);
	    cosFile4.clear();	
	    cosPdfFile4.clear();	
	    cosCdfFile4.clear();	
	    pdf.clear();
	    cdf.clear();
	    */
	}else if (LANG == 2 && NA==2){
//	    std::cout<<"NA = "<< NA <<std::endl;
	  TArrayD edes6(NEP),f06(NEP),r6(NEP),a6(NEP);
	  for (int lis1 = 0; lis1 < NEP; lis1++){
	    edes6[lis1] = header->GetLIST(lis1*4 + 0);
	    f06[lis1] = header->GetLIST(lis1*4 + 1);
	    r6[lis1] = header->GetLIST(lis1*4 + 2);
	    a6[lis1] = header->GetLIST(lis1*4 + 3);
//	  std::cout<<lis1<<"  "<<edes6[lis1]<<"  "<<f06[lis1] <<"  "<<r6[lis1]<<"  "<<a6[lis1]<< std::endl;
	  }
	}
	
      }      
      energy5OfMts.push_back(ein);
//      energyPdf5OfMts.push_back(pdf2d);
//      energyCdf5OfMts.push_back(cdf2d);
      ein.clear();
//      pdf2d.clear();
//      cdf2d.clear();
      }else if (LAW == 2){
	int LANG,NW,NL;
	NBT1         = new int[NR]();
	INT1        = new int[NR]();
	fE1_file6    = new double[NP]();
	fp11_file6   = new double[NP]();
	
	for(int cr=0; cr < NR ; cr ++){
	  NBT1[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1_file6[crs]  = tab1->GetX(crs);
	  fp11_file6[crs] = tab1->GetY(crs);
        //std::cout<<"energy "<< fE1_file6[crs] <<" multip "<< fp11_file6[crs] << std::endl;
	}
	TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
        NR2  = tab2->GetN1();
        NE2  = tab2->GetN2();
        NBT2  = new int[NR2]();
        INT2  = new int[NR2]();
	for(int cr=0; cr < NR2 ; cr ++){
	  NBT2[cr] = tab2->GetNBT(cr);
	  INT2[cr] = tab2->GetINT(cr);
	}
        //std::cout<<"LANG "<< LANG <<" LEP "<< LEP <<" NE2 "<< NE2 << std::endl;
	for (int lis = 0; lis < NE2; lis++){
          TNudyEndfList *header = (TNudyEndfList *)recIter.Next();
	  ein.push_back(header->GetC2());
	  LANG   = header->GetL1();
	  NW   = header->GetN1();
	  NL   = header->GetN2();
          //std::cout<<"energy " <<ein[lis] <<" LANG "<< LANG << std::endl;
          if(LANG == 0){
	    lCoef = new TArrayD[NE2];
	    lCoef[lis].Set(NL);
//	    std::cout<<"energy "<< ein[lis] << std::endl; 
	    for (int j = 0; j < NL; j++){
	      lCoef[lis].SetAt(header->GetLIST(j), j);
//	    std::cout<<"legendre coef "<< lCoef[lis].At(j) << std::endl; 
	    }
	  
//        for (int i = 0; i < ein.GetSize(); i++) {
//          printf("Ein = %e\n", ein[i]);
	    int k1 = 0;double fme =0.0;
	    
            do
	    {
	      fme =1.0;
	      double x = -1. + k1*0.05;
	      for (int j = 0; j < lCoef[lis].GetSize(); j++) {
		double leg = ROOT::Math::legendre((unsigned int)(j+1), x);
		fme += 0.5*(2.*(j+1) + 1.)*lCoef[lis].At(j)*leg;
//            printf("a%d = %e leg= %e\n", j, lCoef[i].At(j),leg);
	      }
	      cosFile4.push_back(x);
	      cosPdfFile4.push_back(fme);
//            printf("%e %e\n", x, fme);
	      k1++;
	    }while(k1<41);
	    
	    for(int l=0; l<40; l++){
	      recursionLinearFile4(lis, cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
	    }
	    Sort(cosFile4, cosPdfFile4);
	    cdfGenerateT(cosFile4, cosPdfFile4);
	    for(unsigned long i = 0; i < cosFile4.size(); i++){
	      pdf.push_back(cosFile4[i]);
	      pdf.push_back(cosPdfFile4[i]);
	      cdf.push_back(cosFile4[i]);
	      cdf.push_back(cosCdfFile4[i]);	  
	    }
	    pdf2d.push_back(pdf);
	    cdf2d.push_back(cdf);
	    cosFile4.clear();	
	    cosPdfFile4.clear();	
	    cosCdfFile4.clear();	
	    pdf.clear();
	    cdf.clear();
	    delete[] lCoef;
	  }else if (LANG >0){
	    for (int i = 0; i < NL; i++) {
	      cosFile4.push_back(header->GetLIST(2*i+0));
	      cosPdfFile4.push_back(header->GetLIST(2*i+1));
	    }
	    Sort(cosFile4, cosPdfFile4);
	    cdfGenerateT(cosFile4, cosPdfFile4);
	    for(unsigned long i = 0; i < cosFile4.size(); i++){
	      pdf.push_back(cosFile4[i]);
	      pdf.push_back(cosPdfFile4[i]);
	      cdf.push_back(cosFile4[i]);
	      cdf.push_back(cosCdfFile4[i]);	  
	    }
	    pdf2d.push_back(pdf);
	    cdf2d.push_back(cdf);
	    cosFile4.clear();	
	    cosPdfFile4.clear();	
	    cosCdfFile4.clear();	
	    pdf.clear();
	    cdf.clear();
	  }
	}
	energy5OfMts.push_back(ein);
	cosPdf4OfMts.push_back(pdf2d);
	cosCdf4OfMts.push_back(cdf2d);
	ein.clear();
	pdf2d.clear();
	cdf2d.clear();
      }else if(LAW == 3) {
	NBT1         = new int[NR]();
	INT1        = new int[NR]();
	fE1_file6    = new double[NP]();
	fp11_file6   = new double[NP]();
	
	for(int cr=0; cr < NR ; cr ++){
	  NBT1[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1_file6[crs]  = tab1->GetX(crs);
	  fp11_file6[crs] = tab1->GetY(crs);
//    std::cout<<"energy1 "<< fE1_file6[crs] <<" multip1 "<< fp11_file6[crs] << std::endl;
	}
	TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
	double ZAPR =tab11->GetC1(); 
	//	double AWPR =tab11->GetC2(); 
	//int LIPR =tab11->GetL1(); 
	//int LAW = tab11->GetL2();
	NR = tab11->GetN1();
	NP = tab11->GetN2();
	divr =div(ZAPR,1000);
	//double ResidueA = divr.rem;
	//double ResidueZ = divr.quot;
	NBT2         = new int[NR]();
	INT2        = new int[NR]();
	fE2_file6    = new double[NP]();
	fp12_file6   = new double[NP]();
//std::cout << particleA <<" <- pA  pZ-> "<< particleZ << std::endl;	
//std::cout << ResidueA <<" <- rA  rZ-> "<< ResidueZ << std::endl;	
	for(int cr=0; cr < NR ; cr ++){
	  NBT2[cr] = tab11->GetNBT(cr);
	  INT2[cr] = tab11->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE2_file6[crs]  = tab11->GetX(crs);
	  fp12_file6[crs] = tab11->GetY(crs);
//    std::cout<<"energy "<< fE2_file6[crs] <<" multip "<< fp12_file6[crs] <<"  "<< LIPR << std::endl;
	}
//*****************************************************************************************	
      }else if(LAW == 4){
	
      }else if(LAW == 5){
	
      }else if(LAW == 6){
	int NPSX; double APSX;
	NBT1         = new int[NR]();
	INT1        = new int[NR]();
	fE1_file6    = new double[NP]();
	fp11_file6   = new double[NP]();
	
	for(int cr=0; cr < NR ; cr ++){
	  NBT1[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1_file6[crs]  = tab1->GetX(crs);
	  fp11_file6[crs] = tab1->GetY(crs);
//    std::cout<<"energy "<< fE1_file6[crs] <<" multip "<< fp11_file6[crs] << std::endl;
	}
	TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
	APSX   = header->GetC1();
	NPSX   = header->GetN2();
	double energy =-(1. + AWRI)*QValue[sec->GetMT()]/AWRI, eout =1E-5;
        //std::cout<<" begin incident energy "<< energy <<"  "<< QValue[sec->GetMT()] << std::endl;
	energy += 0.001*energy;
	do
	{
	  ein.push_back(energy);
	  double Ea = AWRI * energy/(1. + AWRI) + QValue[sec->GetMT()];
	  double EiMax = (APSX - AWP ) * Ea /APSX ;
	  //double EStar = energy/(AWP + AWRI);
	  double C[3] = {4./(PI*EiMax*EiMax), 105./(32. * pow(EiMax, 3.5)), 256./(14.*PI*pow(EiMax,5))};
          //std::cout<<energy<<"  "<<fE1_file6[NP-1]<<"  "<<EiMax<<std::endl;
          //std::cout<<"c[0] "<< C[0] <<" c[1] "<< C[1] <<" c[2] "<< C[2] << std::endl;
          //std::cout<<"APSX "<< APSX<<" NPSX "<< NPSX  <<" qval "<<QValue[sec->GetMT()] << std::endl;
	  eout =1E-5;
	  do
	  {
	    double ppeCm = 0.0;
	    ppeCm = C[NPSX-3]*sqrt(eout)*pow((EiMax - eout),1.5*NPSX-4);
	    //std:: cout<<"  eout "<< eout << "  " << ppeCm <<"  "<< energy  << std::endl;
	    //double eLab = cmToLabInelasticE(eout, energy, cosT, AWRI);
	    //double labCos = cmToLabInelasticCosT(eLab, eout, energy, cosT, AWRI);
	    //double bracket = EiMax - (EStar + eLab - 2 * labCos * sqrt(EStar*eLab));
	    //if(bracket >0 )ppeLab = C[NPSX-3]*sqrt(eLab)*pow(bracket,1.5*NPSX - 4);
	    std::cout<<" E = "<<energy <<"  "<<eout <<"  "<< ppeCm << std::endl;
	    energyFile5.push_back(eout);
	    energyPdfFile5.push_back(ppeCm);
	    eout *= 2;
	    }while(eout < EiMax); 
	    Sort(energyFile5, energyPdfFile5);
	    cdfGenerateE(energyFile5, energyPdfFile5);  
	    for(unsigned long i = 0; i < energyFile5.size(); i++){
	      if(energyPdfFile5[i] > 1E-15){
		pdf.push_back(energyFile5[i]);
		pdf.push_back(energyPdfFile5[i]);
		cdf.push_back(energyFile5[i]);
		cdf.push_back(energyCdfFile5[i]);
	      }
	    }
	    pdf2d.push_back(pdf);
	    cdf2d.push_back(cdf);
	    energyFile5.clear();	
	    energyPdfFile5.clear();	
	    energyCdfFile5.clear();	
	    pdf.clear();
	    cdf.clear();
	    energy *= 2;
	}while(energy < fE1_file6[NP-1]);	   
	// angle is isotropic in the CM System for this law
      }else if(LAW == 7){
	NBT1        = new int[NR]();
	INT1        = new int[NR]();
	fE1_file6    = new double[NP]();
	fp11_file6   = new double[NP]();
	
	for(int cr=0; cr < NR ; cr ++){
	  NBT1[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	  }
	for(int crs=0; crs < NP ; crs ++){
	  fE1_file6[crs]  = tab1->GetX(crs);
	  fp11_file6[crs] = tab1->GetY(crs);
//    std::cout<<"energy "<< fE1_file6[crs] <<" multip "<< fp11_file6[crs] << std::endl;
	}

	TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
	NR2 = tab2->GetN1();
	NE2  = tab2->GetN2();
	NBT2        = new int[NR2]();
	INT2        = new int[NR2]();
	for(int cr=0; cr < NR2 ; cr ++){
	  NBT2[cr] = tab2->GetNBT(cr);
	  INT2[cr] = tab2->GetINT(cr);
	}
	for(int cr1=0; cr1 < NE2 ; cr1 ++){
	  TNudyEndfTab2 *tab3 = (TNudyEndfTab2 *)recIter.Next();
	  ein.push_back(tab3->GetC2());
	  int NRM  = tab3->GetN1();
	  int NMU  = tab3->GetN2();
	  NBT3        = new int[NRM]();
	  INT3        = new int[NRM]();
	  for(int cr=0; cr < NRM ; cr ++){
	    NBT3[cr] = tab3->GetNBT(cr);
	    INT3[cr] = tab3->GetINT(cr);
	  }
	  TArrayD cosin(NMU);
	  for(int cr2=0; cr2 < NMU ; cr2 ++){
	    TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
	    cosin[cr2] = tab12->GetC2();
//	  std::cout<<"cosine "<< cosin[cr2] <<"  "<< ein[cr1] << std::endl;
	    NR3 = tab12->GetN1();
	    int NEP  = tab12->GetN2();
	    NBT3        = new int[NR3]();
	    INT3        = new int[NR3]();
	    for(int cr=0; cr < NR3 ; cr ++){
	      NBT3[cr] = tab12->GetNBT(cr);
	      INT3[cr] = tab12->GetINT(cr);
	    }
	  for(int crs=0; crs < NEP ; crs ++){
            energyFile5.push_back(tab12->GetX(crs));
            energyPdfFile5.push_back(tab12->GetY(crs));
	  }
	  cdfGenerateE(energyFile5, energyPdfFile5);
	  for (unsigned long i = 0; i< energyFile5.size(); i++){
//   std::cout<<"energyFile5 "<<std::setw(12)<< energyFile5[i] <<"   "<<std::setw(12)<<energyPdfFile5[i]<<"   "<<std::setw(12)<<energyCdfFile5[i]<< std::endl;
	  }
	  energyFile5.clear();	
	  energyPdfFile5.clear();	
	  energyCdfFile5.clear();	
	 }	
	}
      }
      energy5OfMts.push_back(ein);
      energyPdf5OfMts.push_back(pdf2d);
      energyCdf5OfMts.push_back(cdf2d);
      ein.clear();
      pdf2d.clear();
      cdf2d.clear();
    }  
  }
  Mt5Values.push_back(MtNumbers);
  MtNumbers.clear();
  /*  
  for(unsigned long i = 0; i < energy5OfMts.size(); i++){
      std::cout <<" mt "<<Mt5Values[0][i]<<" size "<< energy5OfMts[i].size() << std::endl;
    for(unsigned long j = 0; j < energy5OfMts[i].size(); j++){
      std::cout << energy5OfMts[i][j] << std::endl;
      for(unsigned long k =0; k < energyPdf5OfMts[i][j].size()/2; k++){
	std::cout << energyPdf5OfMts[i][j][2*k] <<"  "<< energyPdf5OfMts[i][j][2*k + 1] <<"  "<< energyCdf5OfMts[i][j][2*k + 1]<< std::endl;
      //for(unsigned long k =0; k < cosPdf4OfMts[i][j].size()/2; k++){
	//std::cout << cosPdf4OfMts[i][j][2*k] <<"  "<< cosPdf4OfMts[i][j][2*k + 1] <<"  "<< cosCdf4OfMts[i][j][2*k + 1]<< std::endl;
      }
    }
  }
  */
}

TNudyEndfEnergyAng::~TNudyEndfEnergyAng(){}

double TNudyEndfEnergyAng::ThinningDuplicate(std::vector<double> &x1){
  int size = x1.size();
  int size1 = x1.size();
  if(size>2){
    for(int i = 0; i< size1 - 1; i++){
      if(x1[i+1] == x1[i]) x1.erase(x1.begin()+i+1);
    size1 = x1.size();
    }
  }
  size1 = x1.size();
  if(size == size1)return 0;
  ThinningDuplicate(x1);
  return 0;
}
double TNudyEndfEnergyAng::ThinningDuplicate(std::vector<double> &x1,std::vector<double> &x2){
  int size = x1.size();
  int size1 = x1.size();
  if(size>2){
    for(int i = 0; i< size1 - 1; i++){
      if(x1[i+1] == x1[i]){ 
	x1.erase(x1.begin()+i+1);
	x2.erase(x2.begin()+i+1);
      }
      size1 = x1.size();
    }
  }
  size1 = x1.size();
  if(size == size1)return 0;
  ThinningDuplicate(x1,x2);
  return 0; 
}

void TNudyEndfEnergyAng::cdfGenerateE(std::vector<double> &x1,std::vector<double> &x2){
  double energy5Cdf = 0.0;
  int size = x1.size() ;
  for(int cr=0; cr < size ; cr ++){
    energy5Cdf += x2[cr];
  }
  double df = 0.0;
  for(int cr=0; cr < size ; cr ++){
    if(energy5Cdf > 0.0)x2[cr] = x2[cr]/energy5Cdf;
    df += x2[cr];
    energyCdfFile5.push_back(df);
//        std::cout <<"cos = "<< energyFile5[cr] <<"  "<< energyPdfFile5[cr] <<"  "<< energyCdfFile5[cr]  << std::endl;
  }
}

double TNudyEndfEnergyAng::recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
  pdf = TNudyCore::Instance()->Interpolate(NBT3, INT3, NR3, fE3_file5, fp93_file5, NE2, mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1)*(mid - x1)/(x2 - x1);
  if(fabs((pdf - pdfmid1)/pdfmid1)<=1E-3){
    return 0;
  }
  energyFile5.push_back(mid); 
  energyPdfFile5.push_back(pdf);
  recursionLinearFile5Prob(x1, mid, pdf1, pdf);
  recursionLinearFile5Prob(mid, x2, pdf, pdf2);
  return 0;  
}

double TNudyEndfEnergyAng::recursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
  double gx=0.0;	       
  double pe = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE1_file5,fp91_file5,NP,energy);
  double thetae = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,fp92_file5,NE, energy);
  if(thetae > 0.0)gx = TNudyCore::Instance()->Interpolate(NBT3,INT3,NR3,fE3_file5,fp93_file5,NE2, mid/thetae);
  pdf = 0.0;
  pdf = pe * gx;
  double pdfmid1 = pdf1 + (pdf2 - pdf1)*(mid - x1)/(x2 - x1);
  if(fabs((pdf - pdfmid1)/pdfmid1)<=1E-3){
    return 0;
  }
  energyFile5.push_back(mid); 
  energyPdfFile5.push_back(pdf);
  recursionLinearFile5GenEva(x1, mid, pdf1, pdf, energy);
  recursionLinearFile5GenEva(mid, x2, pdf, pdf2, energy);
  return 0;  
}

double TNudyEndfEnergyAng::recursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
  double pe = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE1_file5,fp91_file5,NP,energy);
  double I = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,INorm,NE, energy);
  double thetae = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,fp92_file5,NE, energy);
  pdf = 0.0;
  if(I > 0.0)pdf = pe * sqrt(mid) * exp (- mid / thetae)/I;
  double pdfmid1 = pdf1 + (pdf2 - pdf1)*(mid - x1)/(x2 - x1);
  if(fabs((pdf - pdfmid1)/pdfmid1)<=1E-3){
    return 0;
  }
  energyFile5.push_back(mid); 
  energyPdfFile5.push_back(pdf);
  recursionLinearFile5Maxwell(x1, mid, pdf1, pdf, energy);
  recursionLinearFile5Maxwell(mid, x2, pdf, pdf2, energy);
  return 0;  
}

double TNudyEndfEnergyAng::recursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
  double pe = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE1_file5,fp91_file5,NP,energy);
  double I = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,INorm,NE, energy);
  double ae = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,fp92_file5,NE, energy);
  double be = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE3_file5,fp93_file5,NE2, energy);
  pdf = 0.0;
  if(I > 0.0)pdf = pe * sinh(sqrt(be*mid)) * exp (-mid/ae)/I;
  double pdfmid1 = pdf1 + (pdf2 - pdf1)*(mid - x1)/(x2 - x1);
  if(fabs((pdf - pdfmid1)/pdfmid1)<=1E-3){
    return 0;
  }
  energyFile5.push_back(mid); 
  energyPdfFile5.push_back(pdf);
  recursionLinearFile5Watt(x1, mid, pdf1, pdf, energy);
  recursionLinearFile5Watt(mid, x2, pdf, pdf2, energy);
  return 0;
}
void TNudyEndfEnergyAng::Sort(std::vector<double>& x1, std::vector<double>& x2){
  std::multimap<double, double>map;
  std::multimap<double, double>::iterator i;
  int size = x1.size();
  for(int p = 0; p< size; p++)
    map.insert(std::make_pair(x1[p],x2[p]));
    int p1=0;
    for(i=map.begin(); i!=map.end();i++){
      x1[p1]=i->first;
      x2[p1]=i->second;
      p1++;
    }
}

double TNudyEndfEnergyAng::recursionLinearFile4(int i, double x1, double x2, double pdf1, double pdf2){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
//  std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< pdf1 <<"  "<< pdf2 << std::endl;
  for (int j = 0; j < lCoef[i].GetSize(); j++) {
    double leg = ROOT::Math::legendre((unsigned int)(j+1), mid);
    pdf += 0.5*(2.*(j+1) + 1.)*lCoef[i].At(j)*leg;
  }
//  std::cout<< mid <<" mid  "<< pdf1 <<std::endl;
//  double sigmid1 = TNudyCore::Instance()->LinearInterpolation(x1,pdf1,x2,pdf2,mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1)*(mid - x1)/(x2 - x1);
//  std::cout << mid <<" linear "<< pdf <<"  "<< pdfmid1 << std::endl;
  
  if(fabs((pdf - pdfmid1)/pdfmid1)<=1E-3){
    return 0;
  }
  cosFile4.push_back(mid); 
  cosPdfFile4.push_back(pdf);
  recursionLinearFile4(i, x1, mid, pdf1, pdf);
  recursionLinearFile4(i, mid, x2, pdf, pdf2);
  return 0;  
}

void TNudyEndfEnergyAng::cdfGenerateT(std::vector<double> &x1,std::vector<double> &x2){
  double cos4Cdf = 0.0;
  int size = x1.size();
  for(int cr=0; cr < size ; cr ++){
    cos4Cdf += x2[cr];
  }
  double df = 0.0;
  for(int cr=0; cr < size ; cr ++){
    if(cos4Cdf > 0.0)x2[cr] = x2[cr]/cos4Cdf;
    df += x2[cr];
    cosCdfFile4.push_back(df);
  }
}
