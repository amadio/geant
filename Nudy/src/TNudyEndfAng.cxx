// 	This class is reconstructing probability tables for Angular distribution 
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 22, 2016

#include "TNudyEndfAng.h"

#ifdef USE_ROOT
ClassImp(TNudyEndfAng)
#endif

TNudyEndfAng::TNudyEndfAng(){}

//______________________________________________________________________________
TNudyEndfAng::TNudyEndfAng(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    //if (sec->GetMT() == (int)fReaction) {
    TIter recIter(sec->GetRecords());
    TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
    int MT = sec->GetMT();
    MtNumbers.push_back(MT);
    int LTT = sec->GetL2();
    int LI = header->GetL1();
    int LCT = header->GetL2();
    MtLct.push_back(LCT);
    //	printf("LCT = %d LTT = %d LI = %d\n",LCT, LTT, LI);
    //Legendre polynomial coefficients
    if (LTT == 1 && LI == 0) { 
      TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0; i < tab2->GetN2(); i++) {
	TNudyEndfList *tab = (TNudyEndfList *)recIter.Next();
	ein.push_back(tab->GetC2());
	for (int j = 0; j < tab->GetNPL(); j++){
	  lCoef1.push_back (tab->GetLIST(j));
	}
	lCoef.push_back(lCoef1);
	lCoef1.clear();
      }
      for (unsigned long i = 0; i < ein.size(); i++) {
        //printf("Ein = %e\n", ein[i]);
	int k1 = 0;
	double fme =0.0;
	do
	{
	  fme =1.0;
	  double x = -1. + k1*0.05;
          for (unsigned long j = 0; j < lCoef[i].size(); j++) {
	    double leg = ROOT::Math::legendre(j+1, x);
	    fme += 0.5*(2.*(j+1) + 1.)*lCoef[i][j]*leg;
          //printf("a%d = %e leg= %e\n", j, lCoef[i].At(j),leg);
          }
	  cosFile4.push_back(x);
	  cosPdfFile4.push_back(fme);
          //printf("%e %e\n", x, fme);
          k1++;
	}while(k1<41);
	for(int l=0; l<40; l++){
	  recursionLinearLeg(i, cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
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
      energy4OfMts.push_back(ein);
      cosPdf4OfMts.push_back(pdf2d);
      cosCdf4OfMts.push_back(cdf2d);
      ein.clear();
      pdf2d.clear();
      cdf2d.clear();
      lCoef.clear();
      //Tabulated probability tables
    } else if (LTT == 2 && LI == 0) {
      TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0; i < tab2->GetN2(); i++) {
	TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
	ein.push_back(tab->GetC2());
	nr = tab->GetNR();
	np = tab->GetNP();      
	for(int i = 0; i < tab->GetNR(); i++ ){
	  nbt1.push_back(tab->GetNBT(i));
	  int1.push_back(tab->GetINT(i));
	}
	for (int j = 0; j < tab->GetNP(); j++) {
	  cosFile4.push_back(tab->GetX(j));
	  cosPdfFile4.push_back(tab->GetY(j));
	}
	for(int l=0; l < tab->GetNP() - 1; l++){
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
	nbt1.clear();
	int1.clear();
      }
      energy4OfMts.push_back(ein);
      cosPdf4OfMts.push_back(pdf2d);
      cosCdf4OfMts.push_back(cdf2d);
      ein.clear();
      pdf2d.clear();
      cdf2d.clear();
      // Low energy given by legendre polynomial and high energy by tabulated probability tables
    }else if (LTT == 3 && LI == 0) {
      TNudyEndfTab2 *lowE = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0; i < lowE->GetN2(); i++) {
	TNudyEndfList *tab = (TNudyEndfList *)recIter.Next();
	ein.push_back(tab->GetC2());
	for (int j = 0; j < tab->GetNPL(); j++){
	  lCoef1.push_back(tab->GetLIST(j));
	}
	lCoef.push_back(lCoef1);
	lCoef1.clear();	
      }
      for (unsigned long i = 0; i < ein.size(); i++) {
//          printf("Ein = %e\n", ein[i]);
	int k1 = 0;double fme =0.0;
	do
	{
	  fme =1.0;
	  double x = -1. + k1*0.05;
          for (unsigned long j = 0; j < lCoef[i].size(); j++) {
	    double leg = ROOT::Math::legendre(j+1, x);
	    fme += 0.5*(2.*(j+1) + 1.)*lCoef[i][j]*leg;
//            printf("a%d = %e leg= %e\n", j, lCoef[i].At(j),leg);
          }
	  cosFile4.push_back(x);
	  cosPdfFile4.push_back(fme);
//            printf("%e %e\n", x, fme);
          k1++;
	}while(k1<41);
	    
	for(int l=0; l<40; l++){
	  recursionLinearLeg(i, cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
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
      } lCoef.clear();
      TNudyEndfTab2 *highE = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0; i < highE->GetN2(); i++) {
        TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
        ein.push_back(tab->GetC2());
	nr = tab->GetNR();
	np = tab->GetNP();      
	for(int i = 0; i < tab->GetNR(); i++ ){
	  nbt1.push_back(tab->GetNBT(i));
	  int1.push_back(tab->GetINT(i));
	}
        for (int j = 0; j < tab->GetNP(); j++) {
          cosFile4.push_back(tab->GetX(j));
          cosPdfFile4.push_back(tab->GetY(j));
        }
	for(int l=0; l < tab->GetNP() - 1; l++){
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
	nbt1.clear();
	int1.clear();
      }
      energy4OfMts.push_back(ein);
      cosPdf4OfMts.push_back(pdf2d);
      cosCdf4OfMts.push_back(cdf2d);
      ein.clear();
      pdf2d.clear();
      cdf2d.clear();
    } else if (LTT == 0 && LI == 1) {
      ein.push_back(1E-14);
      ein.push_back(1.5E8);
      for(int j = 0; j < 2; j++){
	cosFile4.push_back(1);
	cosPdfFile4.push_back(0.5);
	cosFile4.push_back(-1);
	cosPdfFile4.push_back(0.5);
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
      energy4OfMts.push_back(ein);
      cosPdf4OfMts.push_back(pdf2d);
      cosCdf4OfMts.push_back(cdf2d);
      ein.clear();
      pdf2d.clear();
      cdf2d.clear();
    }
      // Low energy given by legendre polynomial and high energy by tabulated probability tables
  }//end while loop
  Mt4Values.push_back(MtNumbers);
  Mt4Lct.push_back(MtLct);
  std::vector<int>().swap(MtNumbers);
  std::vector<int>().swap(MtLct);
  ///*
  for(unsigned long i = 0; i < energy4OfMts.size(); i++){
    std::cout <<" mt "<<Mt4Values[0][i]<<" size "<< energy4OfMts[i].size() << std::endl;
    for(unsigned long j = 0; j < energy4OfMts[i].size(); j++){
      std::cout << energy4OfMts[i][j] << std::endl;
      for(unsigned long k =0; k < cosPdf4OfMts[i][j].size()/2; k++){
	std::cout << cosPdf4OfMts[i][j][2*k] <<"  "<< cosPdf4OfMts[i][j][2*k + 1] <<"  "<< cosCdf4OfMts[i][j][2*k + 1]<< std::endl;
      }
    }
  }
  //*/
}//end class

TNudyEndfAng::~TNudyEndfAng(){}

double TNudyEndfAng::recursionLinearLeg(int i, double x1, double x2, double pdf1, double pdf2){
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
  recursionLinearLeg(i, x1, mid, pdf1, pdf);
  recursionLinearLeg(i, mid, x2, pdf, pdf2);
  return 0;  
}
double TNudyEndfAng::recursionLinearProb(double x1, double x2, double pdf1, double pdf2){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
  //std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< pdf1 <<"  "<< pdf2 << std::endl;
  pdf = TNudyCore::Instance()->Interpolate(nbt1,int1,nr,cosFile4,cosPdfFile4,np,mid);
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
