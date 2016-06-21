// 	This class is reconstructing probability tables for Energy distribution 
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 22, 2016

#include "TNudyEndfEnergyAng.h"

#ifdef USE_ROOT
ClassImp(TNudyEndfEnergyAng)
#endif

TNudyEndfEnergyAng::TNudyEndfEnergyAng(){}

//______________________________________________________________________________
TNudyEndfEnergyAng::TNudyEndfEnergyAng(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
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
	for(int cr=0; cr < NR ; cr ++){
	  nbt1.push_back (tab1->GetNBT(cr));
	  int1.push_back (tab1->GetINT(cr));
	}
	MtLct.push_back(1);
	for(int crs=0; crs < NP ; crs ++){
	  fE1.push_back (tab1->GetX(crs));
	  fP1.push_back (tab1->GetY(crs));
//    std::cout<<"energy "<< fE1_file6[crs] <<" multip "<< fp11_file6[crs] << std::endl;
	}
	TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
        LANG  = tab2->GetL1();
        LEP  = tab2->GetL2();
        nr2  = tab2->GetN1();
        np2  = tab2->GetN2();
//      std::cout<<"LANG "<< LANG <<" LEP "<< LEP <<" NE2 "<< NE2 << std::endl;
	for (int lis = 0; lis < np2; lis++){
          TNudyEndfList *header = (TNudyEndfList *)recIter.Next();
	  ein.push_back(header->GetC2());
	  ND   = header->GetL1();
	  NA   = header->GetL2();
	  NW   = header->GetN1();
	  NEP  = header->GetN2();
          std::cout <<lis<<"  "<<ein[lis] <<" ND "<< ND <<" NA "<< NA <<" NEP "<< NEP << std::endl;
          if(LANG == 2 && NA==1){
	    for (int lis1 = 0; lis1 < NEP; lis1++){
	      edes6.push_back (header->GetLIST(lis1*3 + 0));
	      f06.push_back (header->GetLIST(lis1*3 + 1));
	      r6.push_back (header->GetLIST(lis1*3 + 2));
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
	      for (unsigned long j = 0; j < edes6.size(); j++) {
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
	    TNudyCore::Instance()->cdfGenerateT(energyFile5, energyPdfFile5, energyCdfFile5);
	    TNudyCore::Instance()->cdfGenerateT(cosFile4, cosPdfFile4, cosCdfFile4);
	    for(unsigned long i = 0; i < cosFile4.size(); i++){
	      pdf.push_back(cosFile4[i]);
	      pdf.push_back(cosPdfFile4[i]);
	      pdf.push_back(energyFile5[i]);
	      pdf.push_back(energyPdfFile5[i]);	      
	      cdf.push_back(cosFile4[i]);
	      cdf.push_back(cosCdfFile4[i]);	  
	      cdf.push_back(energyFile5[i]);
	      cdf.push_back(energyCdfFile5[i]);
	    }
	    pdf2d.push_back(pdf);
	    cdf2d.push_back(cdf);
	    cosFile4.clear(); cosPdfFile4.clear(); cosCdfFile4.clear();	
	    pdf.clear(); cdf.clear();
	    edes6.clear(); f06.clear(); r6.clear();
	    nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
	}else if (LANG == 2 && NA==2){
//	    std::cout<<"NA = "<< NA <<std::endl;
	  for (int lis1 = 0; lis1 < NEP; lis1++){
	    edes6.push_back (header->GetLIST(lis1*4 + 0));
	    f06.push_back (header->GetLIST(lis1*4 + 1));
	    r6.push_back (header->GetLIST(lis1*4 + 2));
	    a6.push_back (header->GetLIST(lis1*4 + 3));
//	  std::cout<<lis1<<"  "<<edes6[lis1]<<"  "<<f06[lis1] <<"  "<<r6[lis1]<<"  "<<a6[lis1]<< std::endl;
	  }
	    edes6.clear(); f06.clear(); r6.clear(); a6.clear();
	}
      }      
      energy5OfMts.push_back(ein);
      energyPdf5OfMts.push_back(pdf2d);
      energyCdf5OfMts.push_back(cdf2d);
      ein.clear();
      pdf2d.clear();
      cdf2d.clear();
      }else if (LAW == 2){
	MtLct.push_back(2);
	int LANG,NW,NL;
	for(int cr=0; cr < NR ; cr ++){
	  nbt1.push_back (tab1->GetNBT(cr));
	  int1.push_back (tab1->GetINT(cr));
	}
	TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
        np2  = tab2->GetN2();
        //std::cout<<"LANG "<< LANG <<" LEP "<< LEP <<" NE2 "<< NE2 << std::endl;
	for (int lis = 0; lis < np2; lis++){
          TNudyEndfList *header = (TNudyEndfList *)recIter.Next();
	  ein.push_back(header->GetC2());
	  LANG   = header->GetL1();
	  NW   = header->GetN1();
	  NL   = header->GetN2();
	  NP = NL;
          //std::cout<<"energy " <<ein[lis] <<" LANG "<< LANG << std::endl;
          if(LANG == 0){
//	    std::cout<<"energy "<< ein[lis] << std::endl; 
	    for (int j = 0; j < NL; j++){
	      lCoef1.push_back (header->GetLIST(j));
	    }
	    lCoef.push_back(lCoef1);

	    int k1 = 0;double fme =0.0;
            do
	    {
	      fme =1.0;
	      double x = -1. + k1*0.05;
	      for (unsigned long j = 0; j < lCoef[lis].size(); j++) {
		double leg = ROOT::Math::legendre(j+1, x);
		fme += 0.5*(2.*(j+1) + 1.)*lCoef[lis][j]*leg;
	      }
	      cosFile4.push_back(x);
	      cosPdfFile4.push_back(fme);
	      k1++;
	    }while(k1<41);
	    
	    for(int l=0; l<40; l++){
	      recursionLinearFile4(lis, cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
	    }
	    TNudyCore::Instance()->Sort(cosFile4, cosPdfFile4);
	    TNudyCore::Instance()->cdfGenerateT(cosFile4, cosPdfFile4, cosCdfFile4);
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
	    lCoef1.clear();
	  }else if (LANG >0){
	    for (int i = 0; i < NL; i++) {
	      cosFile4.push_back(header->GetLIST(2*i+0));
	      cosPdfFile4.push_back(header->GetLIST(2*i+1));
	    }
	    int size = cosFile4.size();
	    for(int l = 0; l < size - 1; l++){
	      recursionLinearProb(cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
	    }	    
	    TNudyCore::Instance()->Sort(cosFile4, cosPdfFile4);
	    TNudyCore::Instance()->cdfGenerateT(cosFile4, cosPdfFile4, cosCdfFile4);
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
	energyPdf5OfMts.push_back(pdf2d);
	energyCdf5OfMts.push_back(cdf2d);
	ein.clear();
	pdf2d.clear();
	cdf2d.clear();
	lCoef.clear();
	nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
      }else if(LAW == 3) {
	for(int cr=0; cr < NR ; cr ++){
	  nbt1.push_back (tab1->GetNBT(cr));
	  int1.push_back (tab1->GetINT(cr));
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1.push_back (tab1->GetX(crs));
	  fP1.push_back (tab1->GetY(crs));
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
	for(int cr=0; cr < NR ; cr ++){
	  nbt2.push_back (tab11->GetNBT(cr));
	  int2.push_back (tab11->GetINT(cr));
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE2.push_back (tab11->GetX(crs));
	  fP2.push_back (tab11->GetY(crs));
//    std::cout<<"energy "<< fE2_file6[crs] <<" multip "<< fp12_file6[crs] <<"  "<< LIPR << std::endl;
	}
	nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
	nbt2.clear(); int2.clear(); fE2.clear(); fP2.clear();
//*****************************************************************************************	
      }else if(LAW == 4){
	
      }else if(LAW == 5){
	
      }else if(LAW == 6){
	MtLct.push_back(2);
	int NPSX; double APSX;
	
	for(int cr=0; cr < NR ; cr ++){
	  nbt1.push_back (tab1->GetNBT(cr));
	  int1.push_back (tab1->GetINT(cr));
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1.push_back (tab1->GetX(crs));
	  fP1.push_back (tab1->GetY(crs));
//    std::cout<<"energy "<< fE1_file6[crs] <<" multip "<< fp11_file6[crs] << std::endl;
	}
	TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
	APSX   = header->GetC1();
	NPSX   = header->GetN2();
	double energy = std::fabs(-(1. + AWRI)*QValue[sec->GetMT()]/AWRI), eout =1E-5;
        //std::cout<<" begin incident energy "<< energy <<"  "<< QValue[sec->GetMT()] << std::endl;
	energy += 0.001*(energy);
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
	    //std::cout<<" E = "<<energy <<"  "<<eout <<"  "<< ppeCm <<"  "<< QValue[sec->GetMT()] << std::endl;
	    energyFile5.push_back(eout);
	    energyPdfFile5.push_back(ppeCm);
	    eout *= 2;
	    }while(eout < EiMax); 
	    TNudyCore::Instance()->Sort(energyFile5, energyPdfFile5);
	    TNudyCore::Instance()->cdfGenerateT(energyFile5, energyPdfFile5, energyCdfFile5);  
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
	}while(energy < fE1[NP-1]);	   
	energy5OfMts.push_back(ein);
	energyPdf5OfMts.push_back(pdf2d);
	energyCdf5OfMts.push_back(cdf2d);
	ein.clear();
	pdf2d.clear();
	cdf2d.clear();
	nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
	// angle is isotropic in the CM System for this law
      }else if(LAW == 7){
	MtLct.push_back(2);
	for(int cr=0; cr < NR ; cr ++){
	  nbt1.push_back (tab1->GetNBT(cr));
	  int1.push_back (tab1->GetINT(cr));
	  }
	for(int crs=0; crs < NP ; crs ++){
	  fE1.push_back (tab1->GetX(crs));
	  fP1.push_back (tab1->GetY(crs));
//    std::cout<<"energy "<< fE1_file6[crs] <<" multip "<< fp11_file6[crs] << std::endl;
	}

	TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
	NR2 = tab2->GetN1();
	NE2  = tab2->GetN2();
	for(int cr=0; cr < NR2 ; cr ++){
	  nbt2.push_back (tab2->GetNBT(cr));
	  int2.push_back (tab2->GetINT(cr));
	}
	for(int cr1=0; cr1 < NE2 ; cr1 ++){
	  TNudyEndfTab2 *tab3 = (TNudyEndfTab2 *)recIter.Next();
	  ein.push_back(tab3->GetC2());
	  int NMU  = tab3->GetN2();
	  TArrayD cosin(NMU);
	  for(int cr2=0; cr2 < NMU ; cr2 ++){
	    TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
	    cosin[cr2] = tab12->GetC2();
//	  std::cout<<"cosine "<< cosin[cr2] <<"  "<< ein[cr1] << std::endl;
	    nr3 = tab12->GetN1();
	    np3 = tab12->GetN2();
	    for(int cr=0; cr < nr3 ; cr ++){
	      nbt3.push_back (tab12->GetNBT(cr));
	      int3.push_back (tab12->GetINT(cr));
	    }
	  for(int crs=0; crs < np3 ; crs ++){
            energyFile5.push_back(tab12->GetX(crs));
            energyPdfFile5.push_back(tab12->GetY(crs));
	  }
	  for(int cr=0; cr < np3 - 1 ; cr ++){
	    recursionLinearFile5Prob(energyFile5[cr], energyFile5[cr+1], energyPdfFile5[cr], energyPdfFile5[cr+1]);
	  }
	  TNudyCore::Instance()->Sort(energyFile5, energyPdfFile5);
	  TNudyCore::Instance()->ThinningDuplicate(energyFile5, energyPdfFile5);
	  TNudyCore::Instance()->cdfGenerateT(energyFile5, energyPdfFile5, energyCdfFile5);
	  for(unsigned long i = 0; i < energyFile5.size(); i++){
	    pdf.push_back(energyFile5[i]);
	    pdf.push_back(energyPdfFile5[i]);
	    cdf.push_back(energyFile5[i]);
	    cdf.push_back(energyCdfFile5[i]);
	    std::cout << energyFile5[i] << "  "<< energyPdfFile5[i] <<"  "<< energyCdfFile5[i] << std::endl;
	  }
	  pdf2d.push_back(pdf);
	  cdf2d.push_back(cdf);
	  energyFile5.clear();	
	  energyPdfFile5.clear();	
	  energyCdfFile5.clear();	
	  pdf.clear();
	  cdf.clear();
	  nbt3.clear();
	  int3.clear();
	 }	
	}
      energy5OfMts.push_back(ein);
      energyPdf5OfMts.push_back(pdf2d);
      energyCdf5OfMts.push_back(cdf2d);
      ein.clear();
      pdf2d.clear();
      cdf2d.clear();
      nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
      nbt2.clear(); int2.clear(); 
      nbt3.clear(); int3.clear(); 
      }
    }  
  }
  Mt5Values.push_back(MtNumbers);
  MtNumbers.clear();
  Mt4Lct.push_back(MtLct);
  MtLct.clear();
 // /*  
  for(unsigned long i = 0; i < energy5OfMts.size(); i++){
      std::cout <<" mt "<<Mt5Values[0][i]<<" size "<< energy5OfMts[i].size() << std::endl;
    for(unsigned long j = 0; j < energy5OfMts[i].size(); j++){
      std::cout << energy5OfMts[i][j] <<"  "<< energyPdf5OfMts[i][j].size() << std::endl;
      for(unsigned long k =0; k < energyPdf5OfMts[i][j].size()/2; k++){
	std::cout << energyPdf5OfMts[i][j][2*k] <<"  "<< energyPdf5OfMts[i][j][2*k + 1] <<"  "<< energyCdf5OfMts[i][j][2*k + 1]<< std::endl;
      //for(unsigned long k =0; k < cosPdf4OfMts[i][j].size()/2; k++){
	//std::cout << cosPdf4OfMts[i][j][2*k] <<"  "<< cosPdf4OfMts[i][j][2*k + 1] <<"  "<< cosCdf4OfMts[i][j][2*k + 1]<< std::endl;
      }
    }
  }
 // */
}

TNudyEndfEnergyAng::~TNudyEndfEnergyAng(){}

double TNudyEndfEnergyAng::recursionLinearFile4(int i, double x1, double x2, double pdf1, double pdf2){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
//  std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< pdf1 <<"  "<< pdf2 << std::endl;
  for (unsigned long j = 0; j < lCoef[i].size(); j++) {
    double leg = ROOT::Math::legendre(j+1, mid);
    pdf += 0.5*(2.*(j+1) + 1.)*lCoef[i][j]*leg;
  }
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
double TNudyEndfEnergyAng::recursionLinearProb(double x1, double x2, double pdf1, double pdf2){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
  //std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< pdf1 <<"  "<< pdf2 << std::endl;
  pdf = TNudyCore::Instance()->Interpolate(nbt1,int1,NR,cosFile4,cosPdfFile4,NP,mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1)*(mid - x1)/(x2 - x1);
//  std::cout << mid <<" linear "<< pdf <<"  "<< pdfmid1 << std::endl;
  
  if(fabs((pdf - pdfmid1)/pdfmid1)<=1E-3){
    return 0;
  }
  cosFile4.push_back(mid); 
  cosPdfFile4.push_back(pdf);
  recursionLinearProb(x1, mid, pdf1, pdf);
  recursionLinearProb(mid, x2, pdf, pdf2);
  return 0;  
}

double TNudyEndfEnergyAng::recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
  pdf = TNudyCore::Instance()->Interpolate(nbt3, int3, nr3, energyFile5, energyPdfFile5, np3, mid);
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

