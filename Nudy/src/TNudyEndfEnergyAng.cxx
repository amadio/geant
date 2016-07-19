// 	This class is reconstructing probability tables for Energy distribution 
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 22, 2016

#include "TList.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfSec.h"
#include "TNudyEndfTab1.h"
#include "TNudyEndfTab2.h"
#include "TNudyEndfList.h"
#include "TNudyCore.h"
#include "TNudyEndfEnergyAng.h"
#include "Math/SpecFuncMathMore.h"

#ifdef USE_ROOT
ClassImp(TNudyEndfEnergyAng)
#include "TRandom3.h"
#endif

TNudyEndfEnergyAng::TNudyEndfEnergyAng(){}

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
      //std::cout<<k<<" NK "<<NK<< std::endl;
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      div_t divr; 
      AWRI = sec->GetC2();
      int MT = sec->GetMT();
      //      int LCT = sec->GetL2();
      divr = div(sec->GetC1(),1000);
      double ZA = divr.quot;
      double AA = divr.rem;
      double NA1 = AA -ZA;
      int ZAP =tab1->GetC1(); 
      double AWP =tab1->GetC2(); 
      //int LIP =tab1->GetL1(); 
      int LAW = tab1->GetL2();
      if(LAW == 3 || LAW == 4 || LAW == 0)continue;
      law.push_back(LAW);
      NR = tab1->GetN1();
      NP = tab1->GetN2();
      divr =div(ZAP,1000);
      int particleA = divr.rem;
      int particleZ = divr.quot;
      zd.push_back(particleZ);
      ad.push_back(particleA);
      MtNumbers6.push_back(MT);
      //std::cout<<"LAW = "<< LAW <<" MT "<<MT <<" NP "<< NP <<" ZAP = "<< ZAP <<" AWP "<<AWP <<" parZ "<< particleZ <<" parA "<< particleA << std::endl;
      if(LAW == 2){
        MtNumbers4.push_back(MT);
	int LANG,NL;
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
	  //NW   = header->GetN1();
	  NL   = header->GetN2();
	  NP = NL;
	  //std::cout<<"energy " <<ein[lis] <<" LANG "<< LANG << std::endl;
	  if(LANG == 0){
	    //std::cout<<"energy "<< ein[lis] << std::endl; 
	    for (int j = 0; j < NL; j++){
	      lCoef1.push_back (header->GetLIST(j));
	    }
	    lCoef.push_back(lCoef1);

	    int k1 = 0;double fme =0.0;
	    do
	    {
	      fme =0.5;
	      double x = -1. + k1*0.02;
	      for (unsigned long j = 1; j < lCoef[lis].size(); j++) {
		double leg = ROOT::Math::legendre(j, x);
		fme += 0.5*(2.*j + 1.)*lCoef[lis][j]*leg;
	      }
	      cosFile4.push_back(x);
	      cosPdfFile4.push_back(fme);
	      k1++;
	    }while(k1<101);
	    
	    //for(int l = 0; l < 100; l++){
	      //recursionLinearFile4(lis, cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
	    //}
	    TNudyCore::Instance()->Sort(cosFile4, cosPdfFile4);
	    TNudyCore::Instance()->cdfGenerateT(cosFile4, cosPdfFile4, cosCdfFile4);
	    for(unsigned long i = 0; i < cosFile4.size(); i++){
	      cosc.push_back(cosFile4[i]);
	      pdfc.push_back(cosPdfFile4[i]);
	      cdfc.push_back(cosCdfFile4[i]);	  
	      //std::cout << cosFile4[i] << "  "<< cosPdfFile4[i] <<"  "<< cosCdfFile4[i] << std::endl;
	    }
	    cos2dc.push_back(cosc);
	    pdf2dc.push_back(pdfc);
	    cdf2dc.push_back(cdfc);
	    cosFile4.clear();	
	    cosPdfFile4.clear();	
	    cosCdfFile4.clear();	
	    cosc.clear();
	    pdfc.clear();
	    cdfc.clear();
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
	      cosc.push_back(cosFile4[i]);
	      pdfc.push_back(cosPdfFile4[i]);
	      cdfc.push_back(cosCdfFile4[i]);	  
	    }
	    cos2dc.push_back(cosc);
	    pdf2dc.push_back(pdfc);
	    cdf2dc.push_back(cdfc);
	    cosFile4.clear();	
	    cosPdfFile4.clear();	
	    cosCdfFile4.clear();	
	    cosc.clear();
	    pdfc.clear();
	    cdfc.clear();
	  }
	}
	ein2dc.push_back(ein);
	cos3dc.push_back(pdf2dc);
	pdf3dc.push_back(pdf2dc);
	cdf3dc.push_back(cdf2dc);
	ein.clear();
	cos2dc.clear();
	pdf2dc.clear();
	cdf2dc.clear();
	lCoef.clear();
	nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();  
	continue;
      }
      MtNumbers.push_back(MT);
      //std::cout<<"ZAP = "<< ZAP <<" AWP "<<AWP <<" parZ "<< particleZ <<" parA "<< particleA << std::endl;
//---------------------------------------------  
      switch(LAW){
	case 1:
	  {
	    int NA,NEP,LANG;
	    for(int cr=0; cr < NR ; cr ++){
	      nbt1.push_back (tab1->GetNBT(cr));
	      int1.push_back (tab1->GetINT(cr));
	    }
	    MtLct.push_back(1);
	    for(int crs=0; crs < NP ; crs ++){
	      fE1.push_back (tab1->GetX(crs));
	      fP1.push_back (tab1->GetY(crs));
	    }
	    TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
	    LANG  = tab2->GetL1();
	    //LEP  = tab2->GetL2();
	    nr2  = tab2->GetN1();
	    np2  = tab2->GetN2();
            //std::cout<<"LANG "<< LANG <<"  "<< NA<<"  "<< MT<< std::endl;
	    for (int lis = 0; lis < np2; lis++){
	      TNudyEndfList *header = (TNudyEndfList *)recIter.Next();
	      ein.push_back(header->GetC2());
	      double sumein =0;
	      //ND   = header->GetL1();
	      NA   = header->GetL2();
	      NEP  = header->GetN2();
	      //std::cout <<"ein  "<<header->GetC2() <<" NA "<< NA <<" lang "<< LANG << std::endl;
	      if(LANG == 2 && NA <2){
		for (int lis1 = 0; lis1 < NEP; lis1++){
		  edes6.push_back (header->GetLIST(lis1*3 + 0));
		  f06.push_back (header->GetLIST(lis1*3 + 1));
		  r6.push_back (header->GetLIST(lis1*3 + 2));
		  //std::cout<<lis1<<"  "<<edes6[lis1] << std::endl;
		}
		double epsa = header->GetC2() *AWRI/(1. + AWRI);
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
		  double x = -1. + k1*0.02;
		  //std::cout<<" k1 "<< k1 <<" cos "<< y << std::endl;
		  double sumprob = 0; 
		  for (unsigned long j = 0; j < edes6.size(); j++) {
		    //std::cout<<"j "<< j <<" cosine "<< x <<" eoutcm "<< edes6[j] << std::endl;
		    double epsb = edes6[j] *(AWP + AWRI)/(AWRI+1. - AWP);
		    double ea = epsa + Sa, eb = epsb + Sb;
		    double R1 = std::min(ea,Et1), R3 = std::min(eb,Et3);
		    double X1 = R1*eb/ea, X3 = R3*eb/ea; 
		    double a0 = C1 * X1 + C2 * X1*X1*X1 + C3 * Mn1 * mn * pow(X3,4);
		    double prob = (0.5*a0*f06[j]/sinh(a0))*(cosh(a0*x) + r6[j]*sinh(a0*x));
		    //std::cout<<"flaglab "<< flaglab <<"  "<< xlab <<"   "<< energyFile5.size() << std::endl;
		    //if(prob > 1E-15){
		      energyFile5.push_back(edes6[j]);
		      energyPdfFile5.push_back(prob);
		      sumprob += prob ;
		      //std::cout<< eoutlab<<"  "<< prob <<"  "<<xlab<<"  "<< y <<  std::endl;
		    //}
		  }
		  //std::cout<<" energies "<< energyFile5.size() << std::endl;
		  //for(unsigned long i = 0; i < energyFile5.size(); i++){
		    //std::cout << energyFile5[i] << "  "<< energyPdfFile5[i] <<" k1 "<< k1 << std::endl;
		  //}
		  double fsum = 0.0;
		  for(unsigned int cr=1; cr < energyFile5.size() ; cr ++){
		    fsum += 0.5 * (energyPdfFile5[cr] + energyPdfFile5[cr - 1]) * (energyFile5[cr] - energyFile5[cr - 1]);
		  }
		  //if(fsum > 1E-20){
		    sumein += fsum;
		    cosin.push_back(x);
		    cosinpdf.push_back(fsum);
		    //std::cout<<"fsum "<< fsum <<" x "<< x << std::endl;
		    TNudyCore::Instance()->Sort(energyFile5, energyPdfFile5);
		    TNudyCore::Instance()->ThinningDuplicate(energyFile5, energyPdfFile5);
		    TNudyCore::Instance()->cdfGenerateT(energyFile5, energyPdfFile5, energyCdfFile5);
		    for(unsigned long i = 0; i < energyFile5.size(); i++){
		      eoute.push_back(energyFile5[i]);
		      pdfe.push_back(energyPdfFile5[i]);
		      cdfe.push_back(energyCdfFile5[i]);
		      //std::cout <<"energy pdf "<< energyFile5[i] << "  "<< energyPdfFile5[i] <<"  "<< energyCdfFile5[i] << std::endl;
		    }
		    eout2de.push_back(eoute);
		    pdf2de.push_back(pdfe);
		    cdf2de.push_back(cdfe);
		  //}
		  energyFile5.clear();	
		  energyPdfFile5.clear();	
		  energyCdfFile5.clear();	
		  eoute.clear();
		  pdfe.clear();
		  cdfe.clear();
		  k1++;
		}while(k1<101);
		TNudyCore::Instance()->cdfGenerateT(cosin, cosinpdf, cosincdf);
         	//std::cout<<"cos size "<<  cosin.size() <<std::endl;
		//for(unsigned long i = 0; i < cosin.size(); i++){
		//std::cout << cosin[i] << " cospdf "<< cosinpdf[i] <<"  "<< cosincdf[i] << std::endl;
		//}
		//if(sumein > 1E-50){
		  cos2d.push_back(cosin);
		  cosinpdf2d.push_back(cosinpdf);
		  cosincdf2d.push_back(cosincdf);
		  eout3de.push_back(eout2de);
		  pdf3de.push_back(pdf2de);
		  cdf3de.push_back(cdf2de);
		  eout2de.clear();
		  pdf2de.clear();
		  cdf2de.clear();
		  cosin.clear();
		  cosinpdf.clear();
		  cosincdf.clear();
		//}
		//if(sumein < 1E-50 && ein.size()>0)ein.erase(ein.begin()+ ein.size()-1);
		//std::cout << "ein end "<< ein.size() << std::endl;
		edes6.clear(); f06.clear(); r6.clear();
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
	    } else if (LANG == 1){
		for (int lis1 = 0; lis1 < NEP; lis1++){
		  edes6.push_back (header->GetLIST(lis1*(NA+2) + 0));
		  for (int j = 0; j < NA + 1; j++){
		    lCoef1.push_back (header->GetLIST(lis1*(NA+2) + j + 1));
		    //std::cout <<  header->GetLIST(lis1*(NA+2) + j + 1) << std::endl;
		  }
		  lCoef.push_back(lCoef1);
		  lCoef1.clear();
		}
		int k1 = 0; double fme = 0;
		do
		{
		  double sumprob = 0; 
		  double x = -1. + k1*0.01;
		  for (unsigned long m = 0; m < edes6.size(); m++) {
		  fme = 0.5;
		    for (unsigned long j = 0; j < lCoef[m].size(); j++) {
		      double leg = ROOT::Math::legendre(j, x);
		      fme += 0.5*(2.*j + 1.)*lCoef[m][j]*leg;
		      //std::cout<<leg<<" leg "<<lCoef[m][j] <<"  "<<lCoef[m].size() <<std::endl;
		    }
		    //if(fme > 1E-15){
		      energyFile5.push_back(edes6[m]);
		      energyPdfFile5.push_back(fme);
		      sumprob += fme;
		      //std::cout<<"eout "<< edes6[m] <<"  "<< fme << std::endl;
		    //}
		  }
		    double fsum = 0.0;
		    for(unsigned int cr=1; cr < energyFile5.size() ; cr ++){
		      fsum += 0.5 * (energyPdfFile5[cr] + energyPdfFile5[cr - 1]) * (energyFile5[cr] - energyFile5[cr - 1]);
		    }
		  //if(fsum > 1E-20){
		    sumein += fsum;
		    cosin.push_back(x);
		    cosinpdf.push_back(fsum);
		    TNudyCore::Instance()->Sort(energyFile5, energyPdfFile5);
		    TNudyCore::Instance()->ThinningDuplicate(energyFile5, energyPdfFile5);
		    TNudyCore::Instance()->cdfGenerateT(energyFile5, energyPdfFile5, energyCdfFile5);
		    for(unsigned long i = 0; i < energyFile5.size(); i++){
		      eoute.push_back(energyFile5[i]);
		      pdfe.push_back(energyPdfFile5[i]);
		      cdfe.push_back(energyCdfFile5[i]);
		      //std::cout <<sumprob<<" energy "<< energyFile5[i] << "  "<< energyPdfFile5[i] <<"  "<< energyCdfFile5[i] << std::endl;
		    }
		    eout2de.push_back(eoute);
		    pdf2de.push_back(pdfe);
		    cdf2de.push_back(cdfe);
		  //}
		  energyFile5.clear();	
		  energyPdfFile5.clear();	
		  energyCdfFile5.clear();	
		  eoute.clear();
		  pdfe.clear();
		  cdfe.clear();
		  k1++;
		}while(k1<201);
		//if(sumein > 1E-50){
		  TNudyCore::Instance()->cdfGenerateT(cosin, cosinpdf, cosincdf);
		  //for(unsigned long i = 0; i < cosin.size(); i++){
		    //std::cout << cosin[i] << " cospdf "<< cosinpdf[i] <<"  "<< cosincdf[i] << std::endl;
		  //}
		  cos2d.push_back(cosin);
		  cosinpdf2d.push_back(cosinpdf);
		  cosincdf2d.push_back(cosincdf);
		  eout3de.push_back(eout2de);
		  pdf3de.push_back(pdf2de);
		  cdf3de.push_back(cdf2de);
		  eout2de.clear();
		  pdf2de.clear();
		  cdf2de.clear();
		  cosin.clear();
		  cosinpdf.clear();
		  cosincdf.clear();
		//}
		edes6.clear(); 
		lCoef.clear();
		//if(sumein < 1E-50 && ein.size()>0)ein.erase(ein.begin()+ ein.size()-1);
		//std::cout << "ein "<< ein.size() << std::endl;
	    }
	  } 
	    cos3d.push_back(cos2d);
	    cosinpdf3d.push_back(cosinpdf2d);
	    cosincdf3d.push_back(cosincdf2d);
	    ein2d.push_back(ein);
	    eout4de.push_back(eout3de);
	    pdf4de.push_back(pdf3de);
	    cdf4de.push_back(cdf3de);
	    ein.clear();
	    eout3de.clear();
	    pdf3de.clear();
	    cdf3de.clear();
	    cos2d.clear();
	    cosinpdf2d.clear();
	    cosincdf2d.clear();
	    nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
    //---------------------------------------------      
	  }break;
	case 2:
	  {
	    MtLct.push_back(2);
	    int LANG,NL;
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
	      //NW   = header->GetN1();
	      NL   = header->GetN2();
	      NP = NL;
	      //std::cout<<"energy " <<ein[lis] <<" LANG "<< LANG << std::endl;
	      if(LANG == 0){
		//std::cout<<"energy "<< ein[lis] << std::endl; 
		for (int j = 0; j < NL; j++){
		  lCoef1.push_back (header->GetLIST(j));
		}
		lCoef.push_back(lCoef1);

		int k1 = 0;double fme =0.0;
		do
		{
		  fme =0.5;
		  double x = -1. + k1*0.02;
		  for (unsigned long j = 1; j < lCoef[lis].size(); j++) {
		    double leg = ROOT::Math::legendre(j, x);
		    fme += 0.5*(2.*j + 1.)*lCoef[lis][j]*leg;
		  }
		  cosFile4.push_back(x);
		  cosPdfFile4.push_back(fme);
		  k1++;
		}while(k1<101);
		
		//for(int l = 0; l < 100; l++){
		  //recursionLinearFile4(lis, cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
		//}
		TNudyCore::Instance()->Sort(cosFile4, cosPdfFile4);
		TNudyCore::Instance()->cdfGenerateT(cosFile4, cosPdfFile4, cosCdfFile4);
		for(unsigned long i = 0; i < cosFile4.size(); i++){
		  cosc.push_back(cosFile4[i]);
		  pdfc.push_back(cosPdfFile4[i]);
		  cdfc.push_back(cosCdfFile4[i]);	  
		  //std::cout << cosFile4[i] << "  "<< cosPdfFile4[i] <<"  "<< cosCdfFile4[i] << std::endl;
		}
		cos2dc.push_back(cosc);
		pdf2dc.push_back(pdfc);
		cdf2dc.push_back(cdfc);
		cosFile4.clear();	
		cosPdfFile4.clear();	
		cosCdfFile4.clear();	
		cosc.clear();
		pdfc.clear();
		cdfc.clear();
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
		  cosc.push_back(cosFile4[i]);
		  pdfc.push_back(cosPdfFile4[i]);
		  cdfc.push_back(cosCdfFile4[i]);	  
		}
		cos2dc.push_back(cosc);
		pdf2dc.push_back(pdfc);
		cdf2dc.push_back(cdfc);
		cosFile4.clear();	
		cosPdfFile4.clear();	
		cosCdfFile4.clear();	
		cosc.clear();
		pdfc.clear();
		cdfc.clear();
	      }
	    }
	    ein2d.push_back(ein);
	    cos3dc.push_back(pdf2dc);
	    pdf3dc.push_back(pdf2dc);
	    cdf3dc.push_back(cdf2dc);
	    ein.clear();
	    cos2dc.clear();
	    pdf2dc.clear();
	    cdf2dc.clear();
	    lCoef.clear();
	    nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
    //---------------------------------------------      
	  }break;
	case 3:
	case 4:
	  {
	    for(int cr=0; cr < NR ; cr ++){
	      nbt1.push_back (tab1->GetNBT(cr));
	      int1.push_back (tab1->GetINT(cr));
	    }
	    for(int crs=0; crs < NP ; crs ++){
	      fE1.push_back (tab1->GetX(crs));
	      fP1.push_back (tab1->GetY(crs));
	      ein.push_back(tab1->GetX(crs));
	      //std::cout<<"energy1 "<< fE1[crs] <<" multip1 "<< fP1[crs] << std::endl;
	    }
	    nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
	    ein2d.push_back(ein);
	    ein.clear();
  //---------------------------------------------      
	  }break;
	case 5:
	  {
	  
  //---------------------------------------------      
	  }break;
	case 6:
	  {
	    MtLct.push_back(2);
	    int NPSX; double APSX;
	    
	    for(int cr=0; cr < NR ; cr ++){
	      nbt1.push_back (tab1->GetNBT(cr));
	      int1.push_back (tab1->GetINT(cr));
	    }
	    for(int crs=0; crs < NP ; crs ++){
	      fE1.push_back (tab1->GetX(crs));
	      fP1.push_back (tab1->GetY(crs));
	    }
	    TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
	    APSX   = header->GetC1();
	    NPSX   = header->GetN2();
	    double energy = std::fabs(-(1. + AWRI)*iQValue[sec->GetMT()]/AWRI), eout =1E-5;
	    //std::cout<<" begin incident energy "<< energy <<"  "<< iQValue[sec->GetMT()] << std::endl;
	    energy += 0.001*(energy);
	    do
	    {
	      ein.push_back(energy);
	      double Ea = AWRI * energy/(1. + AWRI) + iQValue[sec->GetMT()];
	      double EiMax = (APSX - AWP ) * Ea /APSX ;
	      //double EStar = energy/(AWP + AWRI);
	      double C[3] = {4./(PI*EiMax*EiMax), 105./(32. * pow(EiMax, 3.5)), 256./(14.*PI*pow(EiMax,5))};
	      //std::cout<<energy<<"  "<<fE1[NP-1]<<"  "<<EiMax<<std::endl;
	      //std::cout<<"c[0] "<< C[0] <<" c[1] "<< C[1] <<" c[2] "<< C[2] << std::endl;
	      //std::cout<<"APSX "<< APSX<<" NPSX "<< NPSX  <<" qval "<<iQValue[sec->GetMT()] << std::endl;
	      int k1 = 0;
	      double rn1 = 0, rn2 = 0, rn3 = 0 ,rn4 = 0 ,rn5 = 0 ,rn6 = 0, rn7 = 0, rn8 = 0, rn9 = 0;
	      double p = 0, xp = 0, yp = 0, tp = 0;
	      do
	      {
		double cosCM = -1 + k1 * 0.02;
		cosin.push_back(cosCM);
		int iei = 0;
		do
		{
		  do
		  {
		    rn1 = fRnd->Uniform(1);
		    rn2 = fRnd->Uniform(1);
		  }while(rn1 * rn1 + rn2 * rn2 > 1);
		  do
		  {
		    rn3 = fRnd->Uniform(1);
		    rn4 = fRnd->Uniform(1);
		  }while(rn3 * rn3 + rn4 * rn4 > 1);
		  rn5 = fRnd->Uniform(1);
		  rn6 = fRnd->Uniform(1);
		  rn7 = fRnd->Uniform(1);
		  rn8 = fRnd->Uniform(1);
		  if(NPSX == 3)p = rn5;
		  if(NPSX == 4)p = rn5 * rn6;
		  if(NPSX == 5)p = rn5 * rn6 * rn7 * rn8;
		  rn9 = fRnd->Uniform(1);
		  xp = -rn1*log(rn1 * rn1 + rn2 * rn2)/(rn1 * rn1 + rn2 * rn2) - log(rn9);
		  yp = -rn3*log(rn3 * rn3 + rn4 * rn4)/(rn3 * rn3 + rn4 * rn4) - log(p);
		  tp = xp/(xp + yp);
		  eout = tp * EiMax;
		  double ppeCm = 0.0;
		  ppeCm = C[NPSX-3]*sqrt(eout)*pow((EiMax - eout),1.5*NPSX-4);
		  energyFile5.push_back(eout);
		  energyPdfFile5.push_back(ppeCm);
		  iei++;
		  }while(iei < 300);
		  TNudyCore::Instance()->Sort(energyFile5, energyPdfFile5);
		  cosinpdf.push_back(0.5);
		  TNudyCore::Instance()->cdfGenerateT(energyFile5, energyPdfFile5, energyCdfFile5);
		  for(unsigned long i = 0; i < energyFile5.size(); i++){
		    eoute.push_back(energyFile5[i]);
		    pdfe.push_back(energyPdfFile5[i]);
		    cdfe.push_back(energyCdfFile5[i]);
		    //std::cout << energyFile5[i] << "  "<< energyPdfFile5[i] << std::endl;
		  }
		  eout2de.push_back(eoute);
		  pdf2de.push_back(pdfe);
		  cdf2de.push_back(cdfe);
		  energyFile5.clear();	
		  energyPdfFile5.clear();	
		  energyCdfFile5.clear();	
		  eoute.clear();
		  pdfe.clear();
		  cdfe.clear();
		  k1++;
	      }while(k1 < 101);
	      TNudyCore::Instance()->cdfGenerateT(cosin, cosinpdf, cosincdf);
	      cos2d.push_back(cosin);
	      cosinpdf2d.push_back(cosinpdf);
	      cosincdf2d.push_back(cosincdf);
	      eout3de.push_back(eout2de);
	      pdf3de.push_back(pdf2de);
	      cdf3de.push_back(cdf2de);
	      eout2de.clear();
	      pdf2de.clear();
	      cdf2de.clear();
	      cosin.clear();
	      cosinpdf.clear();
	      cosincdf.clear();
	      energy *= 2;
	    }while(energy < fE1[NP-1] + iQValue[sec->GetMT()]);
	    cos3d.push_back(cos2d);
	    cosinpdf3d.push_back(cosinpdf2d);
	    cosincdf3d.push_back(cosincdf2d);
	    ein2d.push_back(ein);
	    eout4de.push_back(eout3de);
	    pdf4de.push_back(pdf3de);
	    cdf4de.push_back(cdf3de);
	    ein.clear();
	    eout3de.clear();
	    pdf3de.clear();
	    cdf3de.clear();
	    cos2d.clear();
	    cosinpdf2d.clear();
	    cosincdf2d.clear();
	    nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
	    // angle is isotropic in the CM System for this law
    //---------------------------------------------      
	  }break;
	case 7:
	  {
	    MtLct.push_back(1);
	    for(int cr=0; cr < NR ; cr ++){
	      nbt1.push_back (tab1->GetNBT(cr));
	      int1.push_back (tab1->GetINT(cr));
	      }
	    for(int crs=0; crs < NP ; crs ++){
	      fE1.push_back (tab1->GetX(crs));
	      fP1.push_back (tab1->GetY(crs));
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
	      //std::cout <<"energy "<< tab3->GetC2() << std::endl;
	      int NMU  = tab3->GetN2();
	      for(int cr2=0; cr2 < NMU ; cr2 ++){
		TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
		cosin.push_back (tab12->GetC2());
		//std::cout<<"cosine "<< cosin[cr2] <<"  "<< ein[cr1] << std::endl;
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
		double fsum = 0.0;
		for(int cr=1; cr < np3 ; cr ++){
		  fsum += 0.5 * (energyPdfFile5[cr] + energyPdfFile5[cr - 1]) * (energyFile5[cr] - energyFile5[cr - 1]);
		}
		cosinpdf.push_back(fsum);
		for(int cr=0; cr < np3 - 1 ; cr ++){
		  recursionLinearFile5Prob(energyFile5[cr], energyFile5[cr+1], energyPdfFile5[cr], energyPdfFile5[cr+1]);
		}
		TNudyCore::Instance()->Sort(energyFile5, energyPdfFile5);
		TNudyCore::Instance()->ThinningDuplicate(energyFile5, energyPdfFile5);
		TNudyCore::Instance()->cdfGenerateT(energyFile5, energyPdfFile5, energyCdfFile5);
		for(unsigned long i = 0; i < energyFile5.size(); i++){
		  eoute.push_back(energyFile5[i]);
		  pdfe.push_back(energyPdfFile5[i]);
		  cdfe.push_back(energyCdfFile5[i]);
		  //std::cout << energyFile5[i] << "  "<< energyPdfFile5[i] <<"  "<< energyCdfFile5[i] << std::endl;
		}
		eout2de.push_back(eoute);
		pdf2de.push_back(pdfe);
		cdf2de.push_back(cdfe);
		energyFile5.clear();	
		energyPdfFile5.clear();	
		energyCdfFile5.clear();	
		eoute.clear();
		pdfe.clear();
		cdfe.clear();
		nbt3.clear();
		int3.clear();
	      }
	      TNudyCore::Instance()->cdfGenerateT(cosin, cosinpdf, cosincdf);
	      //for(unsigned long i = 0; i < cosin.size(); i++){
		//std::cout << cosin[i] << "  "<< cosinpdf[i] <<"  "<< cosincdf[i] << std::endl;
	      //}
	      cos2d.push_back(cosin);
	      cosinpdf2d.push_back(cosinpdf);
	      cosincdf2d.push_back(cosincdf);
	      eout3de.push_back(eout2de);
	      pdf3de.push_back(pdf2de);
	      cdf3de.push_back(cdf2de);
	      eout2de.clear();
	      pdf2de.clear();
	      cdf2de.clear();
	      cosin.clear();
	      cosinpdf.clear();
	      cosincdf.clear();
	    }
	    cos3d.push_back(cos2d);
	    cosinpdf3d.push_back(cosinpdf2d);
	    cosincdf3d.push_back(cosincdf2d);
	    ein2d.push_back(ein);
	    eout4de.push_back(eout3de);
	    pdf4de.push_back(pdf3de);
	    cdf4de.push_back(cdf3de);
	    ein.clear();
	    eout3de.clear();
	    pdf3de.clear();
	    cdf3de.clear();
	    cos2d.clear();
	    cosinpdf2d.clear();
	    cosincdf2d.clear();
	    nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
	    nbt2.clear(); int2.clear(); 
	    nbt3.clear(); int3.clear(); 
	  }break;
      }
    }
  }
  cos6OfMts.push_back(cos3d);
  cosin6Pdf.push_back(cosinpdf3d);
  cosin6Cdf.push_back(cosincdf3d);
  energy5OfMts.push_back(ein2d);
  energyOut6OfMts.push_back(eout4de);
  energyPdf6OfMts.push_back(pdf4de);
  energyCdf6OfMts.push_back(cdf4de);
  Mt6Values.push_back(MtNumbers6);
  Mt5Values.push_back(MtNumbers);
  Mt4Lct.push_back(MtLct);
  Law6.push_back(law);
  ZD6.push_back(zd);
  AD6.push_back(ad);
  law.clear();
  MtNumbers.clear();
  MtNumbers6.clear();
  zd.clear();
  ad.clear();
  MtLct.clear();
  ein2d.clear();
  eout4de.clear();
  pdf4de.clear();
  cdf4de.clear();
  cos3d.clear();
  Mt4Values.push_back(MtNumbers4);
  energy4OfMts.push_back(ein2dc);
  cos4OfMts.push_back(cos3dc);
  cosPdf4OfMts.push_back(pdf3dc);
  cosCdf4OfMts.push_back(cdf3dc);
  std::vector<int>().swap(MtNumbers4);
  ein2dc.clear();
  cos3dc.clear();
  pdf3dc.clear();
  cdf3dc.clear();
  
//  std::cout << Mt4Values[0].size() << std::endl;
//  std::cout << Mt5Values[0].size() << std::endl;
//  std::cout << Mt6Values[0].size() << std::endl;
//  std::cout << Law6[0].size() << std::endl;
//  for(unsigned int i = 0; i< Mt6Values[0].size(); i++ ){
//  std::cout<<"MT = "<< Mt6Values[0][i]<<"  "<<Mt6Values[0].size()<<"  "<< Law6[0][i] <<"  "<< ZD6[0][i] <<"  "<< AD6[0][i] << std::endl;
//  }
//     for(unsigned i = 0; i < energy5OfMts[0].size(); i++)
//       std::cout<<energy5OfMts[0][i].size()<< std::endl;
//  }
  //for(unsigned int i = 0; i< Mt4Values[0].size(); i++ ){
  //std::cout<<"MT = "<< Mt4Values[0][i]<<"  "<<Mt4Values[0].size()<< std::endl;
  //   for(unsigned j = 0; j < energy4OfMts[0][i].size(); j++)
  //     std::cout<<energy4OfMts[0][i][j]<<"  "<<energy4OfMts[0][i].size()<< std::endl;
  //}
}
//------------------------------------------------------------------------------------------------------

TNudyEndfEnergyAng::~TNudyEndfEnergyAng(){
  cosFile4.shrink_to_fit(); 
  energyFile5.shrink_to_fit();
  cosPdfFile4.shrink_to_fit(); 
  energyPdfFile5.shrink_to_fit();
  cosCdfFile4.shrink_to_fit(); 
  energyCdfFile5.shrink_to_fit();
  edes6.shrink_to_fit();
  f06.shrink_to_fit();
  r6.shrink_to_fit();
  a6.shrink_to_fit();
  fE1.shrink_to_fit(); 
  fP1.shrink_to_fit(); 
  fE2.shrink_to_fit(); 
  fP2.shrink_to_fit(); 
  fE3.shrink_to_fit(); 
  fP3.shrink_to_fit(); 
  INorm.shrink_to_fit();
  law.shrink_to_fit();
  MtNumbers.shrink_to_fit();
  zd.shrink_to_fit();
  ad.shrink_to_fit();
  MtLct.shrink_to_fit();
  nbt1.shrink_to_fit();
  int1.shrink_to_fit();
  nbt2.shrink_to_fit();
  int2.shrink_to_fit();
  nbt3.shrink_to_fit();
  int3.shrink_to_fit();
  ein.shrink_to_fit();
  cosc.shrink_to_fit();
  cdfc.shrink_to_fit();
  pdfc.shrink_to_fit();
  lCoef1.shrink_to_fit();
  cosin.shrink_to_fit();
  cosinpdf.shrink_to_fit();
  cosincdf.shrink_to_fit();
  eoute.shrink_to_fit();
  cdfe.shrink_to_fit();
  pdfe.shrink_to_fit();
  cos2d.shrink_to_fit();
  cosinpdf2d.shrink_to_fit();
  cosincdf2d.shrink_to_fit();
  cos2dc.shrink_to_fit();
  pdf2dc.shrink_to_fit();
  cdf2dc.shrink_to_fit();
  lCoef.shrink_to_fit();
  ein2d.shrink_to_fit();
  cos3d.shrink_to_fit();
  cosinpdf3d.shrink_to_fit();
  cosincdf3d.shrink_to_fit();
  cos3dc.shrink_to_fit();
  pdf3dc.shrink_to_fit();
  cdf3dc.shrink_to_fit();
  eout2de.shrink_to_fit();
  pdf2de.shrink_to_fit();
  cdf2de.shrink_to_fit();
  eout3de.shrink_to_fit();
  pdf3de.shrink_to_fit();
  cdf3de.shrink_to_fit();
  eout4de.shrink_to_fit();
  pdf4de.shrink_to_fit();
  cdf4de.shrink_to_fit();
  Mt4Values.shrink_to_fit();
}
void TNudyEndfEnergyAng::Law2OnlyAngle(){
}
//------------------------------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2){
  double pdf =1.0;
  double mid     = 0.5 * ( x1 + x2 );
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
//------------------------------------------------------------------------------------------------------
int TNudyEndfEnergyAng::GetAd6( int ielemId, int mt ){
  int size = AD6[ielemId].size();
  for (int i = 0; i < size; i++){
    if(Mt6Values [ ielemId ][ i ] == mt){
      return  AD6 [ ielemId ][ i ] ;
    }
  }
  return 0;
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfEnergyAng::GetZd6( int ielemId, int mt ){
  int size = ZD6[ielemId].size();
  for (int i = 0; i < size; i++){
    if(Mt6Values [ ielemId ][ i ] == mt){
      return  ZD6 [ ielemId ][ i ] ;
    }
  }
  return 0;
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfEnergyAng::GetLaw6( int ielemId, int mt ){
  int size = Law6[ielemId].size();
  for (int i = 0; i < size; i++){
    if(Mt6Values [ ielemId ][ i ] == mt){
      return  Law6 [ ielemId ][ i ] ;
    }
  }
  return -1;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::GetCos64(int ielemId, int mt, double energyK){
  fRnd = new TRandom3(0);
  int i = -1;
  for(unsigned int l =0; l < Mt4Values[ielemId].size(); l++){
    if(Mt4Values[ielemId][l] == mt){
      i = l;
      break;
    }
  }
  int min = 0;
  int max = energy4OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= energy4OfMts[ielemId][i][min])min = 0;
  else if (energyK >= energy4OfMts[ielemId][i][max]) min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < energy4OfMts[ielemId][i][mid]) max = mid;
      else min = mid;
    }
  }
  double fraction = (energyK - energy4OfMts[ielemId][i][min])/
                    (energy4OfMts[ielemId][i][min+1] - energy4OfMts[ielemId][i][min]);
		    //std::cout <<" fraction "<< fraction <<"  "<< energyK <<"  "<< energy4OfMts[ielemId][i][min] << std::endl;
  double rnd1 = fRnd->Uniform(1);
  double rnd2 = fRnd->Uniform(1);
  if(rnd2 < fraction)min = min + 1;
  int k =0;
  //std::cout<<" pdf size "<< cosPdf4OfMts[ielemId][i][min].size()/2 << std::endl;
  int size = cosCdf4OfMts[ielemId][i][min].size();
  for(int j = 1; j < size; j++){
    //std::cout<<"cdf "<< cosCdf4OfMts[ielemId][i][min][2 * j ] <<"  "<< cosCdf4OfMts[ielemId][i][min][2 * j + 1] << std::endl;
    if(rnd1 <= cosCdf4OfMts[ielemId][i][min][j]){
      k = j - 1 ;
      if(k >= size - 2) k = size - 2;
      break;
    }
  }
  //for(int j = 0; j < size; j++){
  //  std::cout <<cos4OfMts[ielemId][i][min][j] <<"  "<<cosPdf4OfMts[ielemId][i][min][j] << std::endl;
  //}
  
  //std::cout<< k <<"  "<<cos4OfMts[ielemId][i][min][k]<<"  "<<cosPdf4OfMts[ielemId][i][min][k] <<"  "<<cosCdf4OfMts[ielemId][i][min][ k] << std::endl;
  //std::cout<<" pdf "<<k<<"  "<< cosPdf4OfMts[ielemId][i][min][2 * k + 3]<<"  "<< cosPdf4OfMts[ielemId][i][min][2 * k + 1] << std::endl;
  //std::cout<<" cos "<< cosPdf4OfMts[ielemId][i][min][2 * k + 2]<<"  "<< cosPdf4OfMts[ielemId][i][min][2 * k ] << std::endl;
  double plk = (cosPdf4OfMts[ielemId][i][min][k + 1] - cosPdf4OfMts[ielemId][i][min][k])/
               (cos4OfMts[ielemId][i][min][k + 1] - cos4OfMts[ielemId][i][min][k]);
  double plk2 =  cosPdf4OfMts[ielemId][i][min][k] * cosPdf4OfMts[ielemId][i][min][k];
  double plsq = plk2 + 2 * plk *(rnd1 - cosCdf4OfMts[ielemId][i][min][k]);
  double Ang = 0 ;
  if(plk !=0 && plsq > 0){Ang = cos4OfMts[ielemId][i][min][k] + (sqrt(std::fabs(plsq)) -
                     cosPdf4OfMts[ielemId][i][min][k])/plk;
  }else {
  //std::cout <<"uniform angle "<< mt << std::endl;
    Ang = 2 * rnd1 - 1; 
  }
  return Ang ;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::GetCos6(int ielemId, int mt, double energyK){
  fRnd = new TRandom3(0);
  int i = -1;
  for(unsigned int l =0; l < Mt5Values[ielemId].size(); l++){
    //std::cout<< Mt5Values[ielemId][l] << std::endl;
    if(Mt5Values[ielemId][l] == mt){
      i = l;
      break;
    }
  }
  //std::cout<<" i "<<i <<std::endl;
  int min = 0;
  int max = energy5OfMts[ielemId][i].size() - 1;
  //std::cout<<"mt "<< mt <<"  "<< i <<"  "<<law <<"  "<< energy5OfMts[ielemId][i].size()<<"  "<<energy5OfMts[ielemId][i][max] << std::endl;   
  int mid = 0;
  if (energyK <= energy5OfMts[ielemId][i][min])min = 0;
  else if (energyK >= energy5OfMts[ielemId][i][max]) min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < energy5OfMts[ielemId][i][mid]) max = mid;
      else min = mid;
    }
  }
  double fraction = (energyK - energy5OfMts[ielemId][i][min])/
                    (energy5OfMts[ielemId][i][min+1] - energy5OfMts[ielemId][i][min]);
		    //std::cout <<" fraction "<< fraction <<"  "<< energyK <<"  "<< energy4OfMts[ielemId][i][min] << std::endl;
  double rnd1 = fRnd->Uniform(1);
  double rnd2 = fRnd->Uniform(1);
  if(rnd2 < fraction)min = min + 1;
  double Ang = 0 ;
  //std::cout<<"energy bin "<< min <<"   "<< energy5OfMts[ielemId][i][min] << std::endl; 
      int k = 0;
      int size = cos6OfMts[ielemId][i][min].size(); 
      //std::cout<<"cosine size "<<  size << std::endl; 
      for(unsigned int j = 0; j < cos6OfMts[ielemId][i][min].size(); j++){
	if(rnd1 < cosin6Cdf[ielemId][i][min][j]){
	  k = j - 1 ;
	  if(k >= size - 2) k = size - 2;
	  break;
	}
      }
      //std::cout<<"cosine bin "<< k <<"  "<< cos6OfMts[ielemId][i][min][k] << std::endl; 
      double plk = (cosin6Pdf[ielemId][i][min][k + 1] - cosin6Pdf[ielemId][i][min][k])/
		  (cos6OfMts[ielemId][i][min][k + 1] - cos6OfMts[ielemId][i][min][k]);
      double plk2 =  cosin6Pdf[ielemId][i][min][k] * cosin6Pdf[ielemId][i][min][k];
      double plsq = plk2 + 2 * plk *(rnd1 - cosin6Cdf[ielemId][i][min][k]);
      //std::cout <<"plk "<< plk <<" plk2 "<< plk2 <<"  "<< k << std::endl;
      if(plk !=0 && plsq > 0){Ang = cos6OfMts[ielemId][i][min][k] + (sqrt(std::fabs(plsq)) -
			cosin6Pdf[ielemId][i][min][k])/plk;
      //std::cout<< Ang <<" first "<< rnd1 << std::endl;
      }else {
	Ang = 2 * rnd1 - 1; 
      //std::cout<< Ang <<" sec  "<< rnd1 << std::endl;
      }
      //std::cout<< Ang << std::endl;
      return Ang ;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfEnergyAng::GetEnergy6(int ielemId, int mt, double energyK){
  fRnd = new TRandom3(0);
  int i = -1;
  for(unsigned int l =0; l < Mt5Values[ielemId].size(); l++){
    if(Mt5Values[ielemId][l] == mt){
      i = l;
      break;
    }
  }
  //std::cout<<"mt "<< mt <<"  "<< energyK << std::endl; 
  int min = 0;
  int max = energy5OfMts[ielemId][i].size() - 1;
  int mid = 0;
  if (energyK <= energy5OfMts[ielemId][i][min])min = 0;
  else if (energyK >= energy5OfMts[ielemId][i][max]) min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < energy5OfMts[ielemId][i][mid]) max = mid;
      else min = mid;
    }
  }
  double fraction = (energyK - energy5OfMts[ielemId][i][min])/
                    (energy5OfMts[ielemId][i][min+1] - energy5OfMts[ielemId][i][min]);
  double rnd1 = fRnd->Uniform(1);
  double rnd2 = fRnd->Uniform(1);
  double rnd3 = fRnd->Uniform(1);
  if(rnd2 < fraction)min = min + 1;
  //std::cout<<"energy bin "<< min <<"   "<< energy5OfMts[ielemId][i][min] << std::endl;
    int k = 0;
    int size = cos6OfMts[ielemId][i][min].size();
    for(int j = 0; j < size ; j++){
      if(rnd3 < cosin6Cdf[ielemId][i][min][j]){
	k = j - 1 ;
	if(k >= size - 2) k = size - 2;
	break;
      }
    }
    //std::cout<<"cosine bin "<< k <<"  "<< cos6OfMts[ielemId][i][min][k] << std::endl; 
    int m =0;
    size = energyCdf6OfMts[ielemId][i][min][k].size(); 
    for(unsigned int j = 1; j < energyPdf6OfMts[ielemId][i][min][k].size(); j++){
      if(rnd1 <= energyCdf6OfMts[ielemId][i][min][k][j]){
	m = j - 1;
	if(m >= size - 2) m = size - 2;
	break;
      }
    }
    //std::cout<<"energy bin "<< m <<"  "<< energyOut6OfMts[ielemId][i][min][k][m] << std::endl;     
    double plk = (energyPdf6OfMts[ielemId][i][min][k][m + 1] - energyPdf6OfMts[ielemId][i][min][k][m])/
		(energyOut6OfMts[ielemId][i][min][k][m + 1] - energyOut6OfMts[ielemId][i][min][k][m]);
    double plk2 =  energyPdf6OfMts[ielemId][i][min][k][m]* energyPdf6OfMts[ielemId][i][min][k][m];
    //std::cout <<"plk "<< plk <<" plk2 "<< plk2  << std::endl;    
    double edes = 0;
    if(plk != 0)edes = energyOut6OfMts[ielemId][i][min][k][m] + (sqrt(plk2 + 2 * plk *(rnd1 - energyCdf6OfMts[ielemId][i][min][k][m])) -
		  energyPdf6OfMts[ielemId][i][min][k][m])/plk;
    return edes;
}

  
