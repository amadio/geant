// 	This class is reconstructing probability tables for Angular distribution 
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 22, 2016

#include "TList.h"
#include "TNudyEndfAng.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfCont.h"
#include "TNudyEndfTab1.h"
#include "TNudyEndfTab2.h"
#include "TNudyEndfList.h"
#include "TNudyCore.h"
#include "Math/SpecFuncMathMore.h"

#ifdef USE_ROOT
ClassImp(TNudyEndfAng)
#include "TRandom3.h"
#endif

TNudyEndfAng::TNudyEndfAng () : TNudyEndfRecoPoint(){}
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
    	//printf("LCT = %d LTT = %d LI = %d\n",LCT, LTT, LI);
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
	  double x = -1. + k1*0.02;
          for (unsigned long j = 0; j < lCoef[i].size(); j++) {
	    double leg = ROOT::Math::legendre(j+1, x);
	    fme += 0.5*(2.*(j+1) + 1.)*lCoef[i][j]*leg;
          //printf("a%d = %e leg= %e\n", j, lCoef[i].At(j),leg);
          }
          if(fme > 0.0) {
	    cosFile4.push_back(x);
	    cosPdfFile4.push_back(fme);
	  }
          //printf("%e %e\n", x, fme);
          k1++;
	}while(k1<101);
	for(int l=0; l < 100; l++){
	  recursionLinearLeg(i, cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
	}
	fillPdf1d();
      }
      fillPdf2d();
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
	fillPdf1d();
	nbt1.clear();
	int1.clear();
      }
      fillPdf2d();
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
	  double x = -1. + k1*0.02;
          for (unsigned long j = 0; j < lCoef[i].size(); j++) {
	    double leg = ROOT::Math::legendre(j+1, x);
	    fme += 0.5*(2.*(j+1) + 1.)*lCoef[i][j]*leg;
//            printf("a%d = %e leg= %e\n", j, lCoef[i].At(j),leg);
          }
          if(fme > 0.0) {
	    cosFile4.push_back(x);
	    cosPdfFile4.push_back(fme);
	  }
//            printf("%e %e\n", x, fme);
          k1++;
	}while(k1<101);
	    
	for(int l=0; l < 100; l++){
	  recursionLinearLeg(i, cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
	}
	fillPdf1d();
      } lCoef.clear();
      TNudyEndfTab2 *highE = (TNudyEndfTab2 *)recIter.Next();
      for (int i = 0; i < highE->GetN2(); i++) {
        TNudyEndfTab1 *tab = (TNudyEndfTab1 *)recIter.Next();
        ein.push_back(tab->GetC2());
	//std::cout <<"energy "<< ein[ein.size() - 1] << std::endl;
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
        if(np > 2){
	  for(int l=0; l < tab->GetNP() - 1; l++){
	    recursionLinearProb(cosFile4[l], cosFile4[l+1], cosPdfFile4[l], cosPdfFile4[l+1]); 
	  }
	}
	fillPdf1d();
	nbt1.clear();
	int1.clear();
      }
      fillPdf2d();
    } else if (LTT == 0 && LI == 1) {
      ein.push_back(1E-14);
      ein.push_back(1.5E8);
      for(int j = 0; j < 2; j++){
	cosFile4.push_back(1);
	cosPdfFile4.push_back(0.5);
	cosFile4.push_back(-1);
	cosPdfFile4.push_back(0.5);
	fillPdf1d();
      }
      fillPdf2d();
    }
      // Low energy given by legendre polynomial and high energy by tabulated probability tables
  }//end while loop
  Mt4Values.push_back(MtNumbers);
  Mt4Lct.push_back(MtLct);
  energy4OfMts.push_back(ein2d);
  cos4OfMts.push_back(cos3d);
  cosPdf4OfMts.push_back(pdf3d);
  cosCdf4OfMts.push_back(cdf3d);
  std::vector<int>().swap(MtNumbers);
  std::vector<int>().swap(MtLct);
  ein2d.clear();
  cos3d.clear();
  pdf3d.clear();
  cdf3d.clear();
  /*
  std::cout<<"energies "<< energy4OfMts.size() << std::endl;
  for(unsigned long i = 0; i < energy4OfMts.size(); i++){
    std::cout <<" mt "<<Mt4Values[0][i]<<" size "<< energy4OfMts[i].size() << std::endl;
    for(unsigned long j = 0; j < energy4OfMts[i].size(); j++){
      std::cout << energy4OfMts[i][j] << std::endl;
      for(unsigned long k =0; k < cosPdf4OfMts[i][j].size()/2; k++){
	std::cout << cosPdf4OfMts[i][j][2*k] <<"  "<< cosPdf4OfMts[i][j][2*k + 1] <<"  "<< cosCdf4OfMts[i][j][2*k + 1]<< std::endl;
      }
    }
  }
  */
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
  //std::cout <<x1<<"  " <<x2<<"  "<< pdf1 <<"  "<< pdf2 << std::endl;
  //std::cout <<mid<<"  "<< fabs((pdf - pdfmid1)/pdf) <<"  "<< pdf <<"  "<< pdfmid1 << std::endl;
  if(pdf > 0 && fabs((pdf - pdfmid1)/pdf) <= 2E-4){
    return 0;
  }
  cosFile4.push_back(mid); 
  cosPdfFile4.push_back(pdf);
  recursionLinearLeg(i, x1, mid, pdf1, pdf);
  recursionLinearLeg(i, mid, x2, pdf, pdf2);
  return 0;  
}
//--------------------------------------------------------------------------------------
double TNudyEndfAng::recursionLinearProb(double x1, double x2, double pdf1, double pdf2){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
  //std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< pdf1 <<"  "<< pdf2 << std::endl;
  pdf = TNudyCore::Instance()->Interpolate(nbt1,int1,nr,cosFile4,cosPdfFile4,np,mid);
  double pdfmid1 = pdf1 + (pdf2 - pdf1)*(mid - x1)/(x2 - x1);
//  std::cout << mid <<" linear "<< pdf <<"  "<< pdfmid1 << std::endl;
  
  if(pdf > 0 && fabs((pdf - pdfmid1)/pdf) <= 2E-4){
    return 0;
  }
  cosFile4.push_back(mid); 
  cosPdfFile4.push_back(pdf);
  recursionLinearProb(x1, mid, pdf1, pdf);
  recursionLinearProb(mid, x2, pdf, pdf2);
  return 0;  
}
//-----------------------------------------------------------------------------------------
void TNudyEndfAng::fillPdf1d() {
  TNudyCore::Instance()->Sort(cosFile4, cosPdfFile4);
  TNudyCore::Instance()->cdfGenerateT(cosFile4, cosPdfFile4, cosCdfFile4);
  for(unsigned long i = 0; i < cosFile4.size(); i++){
    cos4.push_back(cosFile4[i]);
    pdf.push_back(cosPdfFile4[i]);
    cdf.push_back(cosCdfFile4[i]);
    //std::cout<<cosFile4[i] <<"  "<<cosPdfFile4[i] <<"  "<< cosCdfFile4[i] << std::endl;
  }
  cos2d.push_back(cos4);
  pdf2d.push_back(pdf);
  cdf2d.push_back(cdf);
  cosFile4.clear();	
  cosPdfFile4.clear();	
  cosCdfFile4.clear();	
  cos4.clear();
  pdf.clear();
  cdf.clear();
}
//--------------------------------------------------------------------------------------
void TNudyEndfAng::fillPdf2d() {
  ein2d.push_back(ein);
  cos3d.push_back(cos2d);
  pdf3d.push_back(pdf2d);
  cdf3d.push_back(cdf2d);
  ein.clear();
  cos2d.clear();
  pdf2d.clear();
  cdf2d.clear();
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfAng::GetCos4(int ielemId, int mt, double energyK){
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
//________________________________________________________________________________________________________________
int TNudyEndfAng::GetCos4Lct(int ielemId, int mt){
  int i =0;
  for(unsigned int l =0; l < Mt4Values[ielemId].size(); l++){
    if(Mt4Values[ielemId][l] == mt){
      i = l;
      break;
    }
  }
  return Mt4Lct[ielemId][i];  
}
