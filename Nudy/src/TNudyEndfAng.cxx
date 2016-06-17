// 	This class is reconstructing probability tables for Angular distribution 
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 22, 2016

#include "TNudyEndfAng.h"
class TNudyEndfFile;
class TNudyEndfSec;
class TNudyEndfCont;
class TNudyEndfTab1;
class TNudyEndfTab2;
class TNudyEndfList;

TNudyEndfAng::TNudyEndfAng(){}

//______________________________________________________________________________
TNudyEndfAng::TNudyEndfAng(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  std::vector<double>ein,cdf,pdf;
  std::vector<std::vector<double> >pdf2d,cdf2d;
  int size1, size2, sizecos;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    //if (sec->GetMT() == (int)fReaction) {
    TIter recIter(sec->GetRecords());
    TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
    int MT = sec->GetMT();
    MtNumbers.push_back(MT);
    int LTT = sec->GetL2();
    int LI = header->GetL1();
    int LCT = header->GetL2();
    //printf("LTT = %d LCT = %d LI = %d\n",LTT, LCT, LI);
    //Legendre polynomial coefficients
    if (LTT == 1 && LI == 0) { 
      TNudyEndfTab2 *subheader = (TNudyEndfTab2 *)recIter.Next();
      lCoef = new TArrayD[subheader->GetN2()];
      for (int i = 0; i < subheader->GetN2(); i++) {
	TNudyEndfList *tab = (TNudyEndfList *)recIter.Next();
	lCoef[i].Set(tab->GetNPL());
	ein.push_back(tab->GetC2());
	for (int j = 0; j < tab->GetNPL(); j++){
	  lCoef[i].SetAt(tab->GetLIST(j), j);
	}
      }
      size1 = ein.size();
      for (int i = 0; i < size1; i++) {
        //printf("Ein = %e\n", ein[i]);
	int k1 = 0;
	double fme =0.0;
	do
	{
	  fme =1.0;
	  double x = -1. + k1*0.05;
          for (int j = 0; j < lCoef[i].GetSize(); j++) {
	    double leg = ROOT::Math::legendre((unsigned int)(j+1), x);
	    fme += 0.5*(2.*(j+1) + 1.)*lCoef[i].At(j)*leg;
          //printf("a%d = %e leg= %e\n", j, lCoef[i].At(j),leg);
          }
	  cosFile4.push_back(x);
	  cosPdfFile4.push_back(fme);
          //printf("%e %e\n", x, fme);
          k1++;
	}while(k1<41);
	for(int l=0; l<40; l++){
	  recursionLinearFile4(i, cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
	}
	Sort(cosFile4, cosPdfFile4);
	cdfGenerateT(cosFile4, cosPdfFile4);
	sizecos = cosFile4.size();
	for(int i = 0; i < sizecos; i++){
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
      delete[] lCoef;
      //Tabulated probability tables
    } else if (LTT == 2 && LI == 0) {
      TNudyEndfTab2 *prob = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0; i < prob->GetN2(); i++) {
	TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
	ein.push_back(tab->GetC2());
	for (int j = 0; j < tab->GetNP(); j++) {
	  cosFile4.push_back(tab->GetX(j));
	  cosPdfFile4.push_back(tab->GetY(j));
	}
	cdfGenerateT(cosFile4, cosPdfFile4);
	sizecos = cosFile4.size();
	for(int i = 0; i < sizecos; i++){
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
      // Low energy given by legendre polynomial and high energy by tabulated probability tables
    }else if (LTT == 3 && LI == 0) {
      TNudyEndfTab2 *lowE = (TNudyEndfTab2 *)recIter.Next();
      lCoef = new TArrayD[lowE->GetN2()];
      for (int i = 0; i < lowE->GetN2(); i++) {
	TNudyEndfList *tab = (TNudyEndfList *)recIter.Next();
	lCoef[i].Set(tab->GetNPL());
	ein.push_back(tab->GetC2());
	for (int j = 0; j < tab->GetNPL(); j++){
	  lCoef[i].SetAt(tab->GetLIST(j), j);
	}
      }
      size2 = ein.size();
      for (int i = 0; i < size2; i++) {
//          printf("Ein = %e\n", ein[i]);
	int k1 = 0;double fme =0.0;
	do
	{
	  fme =1.0;
	  double x = -1. + k1*0.05;
          for (int j = 0; j < lCoef[i].GetSize(); j++) {
	    double leg = ROOT::Math::legendre((unsigned int)(j+1), x);
	    fme += 0.5*(2.*(j+1) + 1.)*lCoef[i].At(j)*leg;
//            printf("a%d = %e leg= %e\n", j, lCoef[i].At(j),leg);
          }
	  cosFile4.push_back(x);
	  cosPdfFile4.push_back(fme);
//            printf("%e %e\n", x, fme);
          k1++;
	}while(k1<41);
	    
	for(int l=0; l<40; l++){
	  recursionLinearFile4(i, cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
	}
	Sort(cosFile4, cosPdfFile4);
	cdfGenerateT(cosFile4, cosPdfFile4);
	sizecos = cosFile4.size();
	for(int i = 0; i < sizecos; i++){
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
      }delete[] lCoef;
      TNudyEndfTab2 *highE = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0; i < highE->GetN2(); i++) {
        TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
        ein.push_back(tab->GetC2());
        for (int j = 0; j < tab->GetNP(); j++) {
          cosFile4.push_back(tab->GetX(j));
          cosPdfFile4.push_back(tab->GetY(j));
        }
	cdfGenerateT(cosFile4, cosPdfFile4);
	sizecos = cosFile4.size();
	for(int i = 0; i < sizecos; i++){
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
    } else if (LTT == 0 && LI == 1) {
      ein.push_back(1E-14);
      ein.push_back(1.5E8);
      for(int j = 0; j < 2; j++){
	cosFile4.push_back(1);
	cosPdfFile4.push_back(0.0);
	cosFile4.push_back(-1);
	cosPdfFile4.push_back(1.0);
	cdfGenerateT(cosFile4, cosPdfFile4);
	sizecos = cosFile4.size();
	for(int i = 0; i < sizecos; i++){
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
  std::vector<int>().swap(MtNumbers);
  ///*
  size1 = energy4OfMts.size();
  for(int i = 0; i < size1; i++){
    std::cout <<" mt "<<Mt4Values[0][i]<<" size "<< energy4OfMts[i].size() << std::endl;
    size2 = energy4OfMts[i].size();
    for(int j = 0; j < size2; j++){
      std::cout << energy4OfMts[i][j] << std::endl;
      int size3 = cosPdf4OfMts[i][j].size()/2;
      for(int k =0; k < size3; k++){
	std::cout << cosPdf4OfMts[i][j][2*k] <<"  "<< cosPdf4OfMts[i][j][2*k + 1] <<"  "<< cosCdf4OfMts[i][j][2*k + 1]<< std::endl;
      }
    }
  }
  //*/
}//end function

TNudyEndfAng::~TNudyEndfAng(){}

double TNudyEndfAng::recursionLinearFile4(int i, double x1, double x2, double pdf1, double pdf2){
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

void TNudyEndfAng::cdfGenerateT(std::vector<double> &x1,std::vector<double> &x2){
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
void TNudyEndfAng::Sort(std::vector<double>& x1, std::vector<double>& x2){
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
