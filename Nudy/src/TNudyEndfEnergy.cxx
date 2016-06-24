// 	This class is reconstructing probability tables for Energy distribution 
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 22, 2016

#include "TList.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfTab1.h"
#include "TNudyEndfTab2.h"
#include "TNudyCore.h"
#include "TNudyEndfEnergy.h"
#include "Math/SpecFuncMathMore.h"
#include "TMath.h"

TNudyEndfEnergy::TNudyEndfEnergy(){}

//______________________________________________________________________________
TNudyEndfEnergy::TNudyEndfEnergy(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    for (int k = 0; k < sec->GetN1(); k++) {
      TIter recIter(sec->GetRecords());
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      int MT = sec->GetMT();
      MtNumbers.push_back(MT);
      int LF = tab1->GetL2();
      std::cout<<" LF = "<< LF <<" MT "<< MT  << std::endl;
      NR = tab1->GetN1();
      NP = tab1->GetN2();
//****************************************************************************	
      // arbitrary tabulated function
      if(LF==1){
	for(int cr=0; cr < NR ; cr ++){
	  nbt1.push_back (tab1->GetNBT(cr));
	  int1.push_back (tab1->GetINT(cr));
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1.push_back (tab1->GetX(crs));
	  fP1.push_back (tab1->GetY(crs));
	}
	TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
	nr2 = tab2->GetN1();
	np2  = tab2->GetN2();
	
	for(int cr=0; cr < nr2 ; cr ++){
	  nbt2.push_back (tab2->GetNBT(cr));
	  int2.push_back (tab2->GetINT(cr));
	}
	for(int cr=0; cr < np2 ; cr ++){
	  TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
          ein.push_back(tab12->GetC2());
	  nr3 = tab12->GetNR();
	  np3 = tab12->GetNP();      
	  for (int i = 0; i < nr3 ; i++ ) {
	    nbt3.push_back(tab12->GetNBT(i));
	    int3.push_back(tab12->GetINT(i));
	  }
	  for(int crs=0; crs < np3 ; crs ++){
	    energyFile5.push_back(tab12->GetX(crs));
	    energyPdfFile5.push_back(tab12->GetY(crs));
	  }
	  for(int cr=0; cr < np3 - 1 ; cr ++){
	    recursionLinearFile5Prob(energyFile5[cr], energyFile5[cr+1], energyPdfFile5[cr], energyPdfFile5[cr+1]);
	  }
	  fillPdf1d();
	  nbt3.clear();
	  int3.clear();
	}
	nbt1.clear();
	int1.clear();
	nbt2.clear();
	int2.clear();
	fE1.clear();
	fP1.clear();
//****************************************************************************	
      // general evaporation spectrum
      }else if(LF==5){
	double u = tab1->GetC1();
	for(int cr=0; cr < NR ; cr ++){
	  nbt1.push_back (tab1->GetNBT(cr));
	  int1.push_back (tab1->GetINT(cr));
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1.push_back (tab1->GetX(crs));
	  fP1.push_back (tab1->GetY(crs));
	}
        TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
	nr2 = tab11->GetN1();
	np2  = tab11->GetN2();
	for(int cr=0; cr < nr2 ; cr ++){
	  nbt2.push_back (tab11->GetNBT(cr));
	  int2.push_back (tab11->GetINT(cr));
	}
	for(int crs=0; crs < np2 ; crs ++){
          fE2.push_back (tab11->GetX(crs));
          fP2.push_back (tab11->GetY(crs));
	  ein.push_back(tab11->GetX(crs));
	}
        TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
	nr3  = tab12->GetN1();
	np3  = tab12->GetN2();
	for(int cr=0; cr < nr3 ; cr ++){
	  nbt3.push_back (tab12->GetNBT(cr));
	  int3.push_back (tab12->GetINT(cr));
	}
	for(int crs=0; crs < np3 ; crs ++){
          fE3.push_back (tab12->GetX(crs));
          fP3.push_back (tab12->GetY(crs));
	}
	double energy, eout;
	for(int i = 0; i < np2; i++){
	  energy = fE2[i];
	  eout =1E-5;
	  do
	  {
	    double gx=0.0;	       
	    double pe = TNudyCore::Instance()->Interpolate(nbt1,int1,NR,fE1,fP1,NP,energy);
	    double thetae = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,fP2,np2, energy);
	    if(thetae > 0.0)gx = TNudyCore::Instance()->Interpolate(nbt3,int3,nr3,fE3,fP3,np3, eout/thetae);
	    double ppe = 0.0;
	    ppe = pe * gx;
	    if(ppe > 0.0) {
	      energyFile5.push_back(eout);
	      energyPdfFile5.push_back(ppe);
	    }
            //std:: cout<<"i = "<< i <<" pe "<< pe <<" thetae "<< thetae <<"  prob "<< ppe <<"  ein "<< energy <<"  eout "<< eout <<"  "<< u << std::endl;
	    eout *= 2;
	  }while(eout < energy - u);
	  int size = energyFile5.size();
	  for(int cr=0; cr < size - 1 ; cr ++){
	    recursionLinearFile5GenEva(energyFile5[cr], energyFile5[cr+1], energyPdfFile5[cr], energyPdfFile5[cr+1],energy);
	  }
	  fillPdf1d();
	}
	nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
	nbt2.clear(); int2.clear(); fE2.clear(); fP2.clear();
	nbt3.clear(); int3.clear(); fE3.clear(); fP3.clear();
//****************************************************************************
      // simple maxwellian fission spectrum
      }else if(LF==7){
	double u = tab1->GetC1();
	for(int cr=0; cr < NR ; cr ++){
	  nbt1.push_back (tab1->GetNBT(cr));
	  int1.push_back (tab1->GetINT(cr));
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1.push_back (tab1->GetX(crs));
	  fP1.push_back (tab1->GetY(crs));
	}
        TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
	nr2 = tab11->GetN1();
	np2  = tab11->GetN2();
	for(int cr=0; cr < nr2 ; cr ++){
	  nbt2.push_back (tab11->GetNBT(cr));
	  int2.push_back (tab11->GetINT(cr));
	}
	for(int crs=0; crs < np2 ; crs ++){
          fE2.push_back (tab11->GetX(crs));
          fP2.push_back (tab11->GetY(crs));
	  double sq = (fE2[crs] - u)/fP2[crs];
	  double sqt = sqrt(sq);
          INorm.push_back (pow(fP2[crs], 1.5)*(0.5*1.7724538529055*erf(sqt)- sqt*exp(-sq)));
	  ein.push_back(tab11->GetX(crs));
	}
	double energy, eout;
	for(int i = 0; i < np2; i++){
	  energy = fE2[i];
	  eout =1E-5;
	  do
	  {
	    double pe = TNudyCore::Instance()->Interpolate(nbt1,int1,NR,fE1,fP1,NP,energy);
	    double I = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,INorm,np2, energy);
	    double thetae = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,fP2,np2, energy);
	    double ppe = 0.0;
	    if(I > 0.0)ppe = pe * sqrt(eout) * exp (-eout/thetae)/I;
	    if(ppe > 0.0){
	      energyFile5.push_back(eout);
	      energyPdfFile5.push_back(ppe);
	    }
//    std:: cout<<"I = "<< I <<" pe "<< pe <<" thetae "<< thetae <<"  prob "<< ppe <<"  ein "<< energy <<"  eout "<< eout << std::endl;
	    eout *= 2;
	  }while(eout < energy - u);
	  int size = energyFile5.size();
	  for(int cr=0; cr < size - 1 ; cr ++){
	    recursionLinearFile5Maxwell(energyFile5[cr], energyFile5[cr+1], energyPdfFile5[cr], energyPdfFile5[cr+1],energy);
	  }
	  fillPdf1d();
	}	   
	nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
	nbt2.clear(); int2.clear(); fE2.clear(); fP2.clear();INorm.clear();
///////////////////////////////////////////////////////////////////////////////////
      // evaporation spectrum
      }else if(LF==9){
	double u = tab1->GetC1();
	for(int cr=0; cr < NR ; cr ++){
	  nbt1.push_back (tab1->GetNBT(cr));
	  int1.push_back (tab1->GetINT(cr));
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1.push_back (tab1->GetX(crs));
	  fP1.push_back (tab1->GetY(crs));
	}
        TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
	nr2 = tab11->GetN1();
	np2  = tab11->GetN2();
	for(int cr=0; cr < nr2 ; cr ++){
	  nbt2.push_back (tab11->GetNBT(cr));
	  int2.push_back (tab11->GetINT(cr));
	}
	for(int crs=0; crs < np2 ; crs ++){
          fE2.push_back (tab11->GetX(crs));
          fP2.push_back (tab11->GetY(crs));
	  INorm.push_back (fP2[crs]*fP2[crs]*(1. - exp(-(fE2[crs] - u)/fP2[crs])*
		      (1.0 + (fE2[crs] - u)/fP2[crs])));
	  ein.push_back(tab11->GetX(crs));
	}
	double energy, eout;
	for(int i = 0; i < np2; i++){
	  energy = fE2[i];
	  eout =1E-5;
	  do
	  {
	    double pe = TNudyCore::Instance()->Interpolate(nbt1,int1,NR,fE1,fP1,NP,energy);
	    double I = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,INorm,np2, energy);
	    double thetae = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,fP2,np2, energy);
	    double ppe = 0.0;
	    if(I > 0.0)ppe = pe * eout * exp (-eout/thetae)/I;
	    if(ppe > 0.0) {
	      energyFile5.push_back(eout);
	      energyPdfFile5.push_back(ppe);
	    }
  //    std:: cout<<"I = "<< I <<" pe "<< pe <<" thetae "<< thetae <<"  prob "<< ppe <<"  ein "<< energy <<"  eout "<< eout << std::endl;
	    eout *= 2;
	  }while(eout < u);
	  int size = energyFile5.size();
	  for(int cr=0; cr < size - 1 ; cr ++){
	    recursionLinearFile5Maxwell(energyFile5[cr], energyFile5[cr+1], energyPdfFile5[cr], energyPdfFile5[cr+1],energy);
	  }
	  fillPdf1d();
	}
      //std:: cout<< std::endl;
	nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
	nbt2.clear(); int2.clear(); fE2.clear(); fP2.clear();INorm.clear();
/////////////////////////////////////////////////////////////////////////////////////////
      // energy dependent watt spectrum
      }else if(LF==11){
	double u = tab1->GetC1();
	
	for(int cr=0; cr < NR ; cr ++){
	  nbt1.push_back (tab1->GetNBT(cr));
	  int1.push_back (tab1->GetINT(cr));
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1.push_back (tab1->GetX(crs));
	  fP1.push_back (tab1->GetY(crs));
	}
	TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
	nr2 = tab11->GetN1();
	np2  = tab11->GetN2();
	for(int cr=0; cr < nr2 ; cr ++){
	  nbt2.push_back (tab11->GetNBT(cr));
	  int2.push_back (tab11->GetINT(cr));
	}

	for(int crs=0; crs < np2 ; crs ++){
          fE2.push_back (tab11->GetX(crs));
          fP2.push_back (tab11->GetY(crs));
	  ein.push_back(tab11->GetX(crs));
	}
	TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
	nr3  = tab12->GetN1();
	np3  = tab12->GetN2();
	for(int cr=0; cr < nr3 ; cr ++){
	  nbt3.push_back (tab12->GetNBT(cr));
	  int3.push_back (tab12->GetINT(cr));
	}
	for(int crs=0; crs < np3 ; crs ++){
          fE3.push_back (tab12->GetX(crs));
          fP3.push_back (tab12->GetY(crs));
	      
	  double a = fP2[crs];
	  double b = fP3[crs];
	  double eua = (fE3[crs] - u )/a;
	    
	  INorm.push_back (0.5*sqrt(0.25*PI*a*a*a*b)*exp(0.25*a*b)*(erf(sqrt(eua)- sqrt(0.25*a*b)) + 
		      erf(sqrt(eua) + sqrt(0.25*a*b))) - a*exp(-eua)*sinh(sqrt(b*(fE3[crs] - u))));
	}
	double energy, eout;
	for(int i = 0; i < np2; i++){
	  energy = fE2[i];
	  eout =1E-5;
	  do
	  {
	    double pe = TNudyCore::Instance()->Interpolate(nbt1,int1,NR,fE1,fP1,NP,energy);
	    double I = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,INorm,np2, energy);
	    double ae = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,fP2,np2, energy);
	    double be = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE3,fP3,np3, energy);
	    double ppe = 0.0;
	    if(I > 0.0)ppe = pe * sinh(sqrt(be*eout)) * exp (-eout/ae)/I;
	    if(ppe > 0.0) {
	      energyFile5.push_back(eout);
	      energyPdfFile5.push_back(ppe);
	    }
  //    std:: cout<<"I = "<< I <<" pe "<< pe <<" ae "<< ae <<"  prob "<< ppe <<"  ein "<< energy <<"  eout "<< eout << std::endl;
	    eout *= 2;
	  }while(eout < energy - u);	   
	  int size = energyFile5.size();
	  for(int cr=0; cr < size - 1 ; cr ++){
	    recursionLinearFile5Watt(energyFile5[cr], energyFile5[cr+1], energyPdfFile5[cr], energyPdfFile5[cr+1],energy);
	  }
	  fillPdf1d();
	}
	nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
	nbt2.clear(); int2.clear(); fE2.clear(); fP2.clear();
	nbt3.clear(); int3.clear(); fE3.clear(); fP3.clear();INorm.clear();
      }else if(LF==12){
	for(int cr=0; cr < NR ; cr ++){
	  nbt1.push_back (tab1->GetNBT(cr));
	  int1.push_back (tab1->GetINT(cr));
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1.push_back (tab1->GetX(crs));
	  fP1.push_back (tab1->GetY(crs));
	}
	TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
	double efl = tab11->GetC1();
	double efh = tab11->GetC2();
	  
	nr2 = tab11->GetN1();
	np2  = tab11->GetN2();
	for(int cr=0; cr < nr2 ; cr ++){
	  nbt2.push_back (tab11->GetNBT(cr));
	  int2.push_back (tab11->GetINT(cr));
	}
	for(int crs=0; crs < np2 ; crs ++){
          fE2.push_back (tab11->GetX(crs));
          fP2.push_back (tab11->GetY(crs));
	  ein.push_back(tab11->GetX(crs));
	}
	  double energy, eout;
	  for(int i = 0; i < np2; i++){
	    energy = fE2[i];
	    eout =1E-5;
	  do
	  {
	    double pe = TNudyCore::Instance()->Interpolate(nbt1,int1,NR,fE1,fP1,NP,energy);
	    double tm = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,fP2,np2, energy);
	    double u1l = (sqrt(eout) - sqrt(efl))*(sqrt(eout) - sqrt(efl))/tm;
	    double u2l = (sqrt(eout) + sqrt(efl))*(sqrt(eout) + sqrt(efl))/tm;
	    double u1h = (sqrt(eout) - sqrt(efh))*(sqrt(eout) - sqrt(efh))/tm;
	    double u2h = (sqrt(eout) + sqrt(efh))*(sqrt(eout) + sqrt(efh))/tm;
  //std::cout<<" u1l "<< u1l <<" u2l "<< u2l <<" u1h "<< u1h <<" u2h "<< u2h << std::endl;
	    double e1ul = ROOT::Math::expint(u1l);
	    double e2ul = ROOT::Math::expint(u2l);
	    double e1uh = ROOT::Math::expint(u1h);
	    double e2uh = ROOT::Math::expint(u2h);
  //std::cout<<" e1ul "<< e1ul <<" e2ul "<< e2ul <<" e1uh "<< e1uh <<" e2uh "<< e2uh << std::endl;
	    double a1 = 1.5;
	    double gamau1l = TMath::Gamma(a1,u1l);
	    double gamau2l = TMath::Gamma(a1,u2l);
	    double gamau1h = TMath::Gamma(a1,u1h);
	    double gamau2h = TMath::Gamma(a1,u2h);
  //std::cout<<" gamau2l "<< gamau2l <<" gamau1l "<< gamau1l <<" gamau2h "<< gamau2h <<" gamau1h "<< gamau1h << std::endl;
	    double gl = (1./(3*sqrt(efl*tm)))*(pow(u2l,1.5)*e2ul - pow(u1l,1.5)*e1ul + gamau2l - gamau1l );
	    double gh = (1./(3*sqrt(efh*tm)))*(pow(u2h,1.5)*e2uh - pow(u1h,1.5)*e1uh + gamau2h - gamau1h );
  /*
  double alp = sqrt(tm);
  double bet = sqrt(efl);
  double a1 = (sqrt(eout) + bet) * (sqrt(eout) + bet)/tm;
  double b1 = (sqrt(30E+6) + bet) * (sqrt(3E7) + bet)/tm;
  double a2 = (sqrt(eout) - bet) * (sqrt(eout) - bet)/tm;
  double b2 = (sqrt(30E+6) - bet) * (sqrt(3E7) - bet)/tm;
  */
	    double ppe = 0.5 * pe * (gl + gh);
	    energyFile5.push_back(eout);
	    energyPdfFile5.push_back(ppe);
  //    std:: cout<< eout <<"  "<< ppe <<"  "<< energy <<"  "<< eout << std::endl;
  //    std:: cout<<"I = "<< I <<" pe "<< pe <<" ae "<< ae <<"  prob "<< ppe <<"  ein "<< energy <<"  eout "<< eout << std::endl;
	    eout *= 2;
	  }while(eout < fE1[NP-1]);
	  fillPdf1d();
	}
	nbt1.clear(); int1.clear(); fE1.clear(); fP1.clear();
	nbt2.clear(); int2.clear(); fE2.clear(); fP2.clear();
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
  for(unsigned long i = 0; i < energy5OfMts.size() ; i++){
      std::cout <<" mt "<<Mt5Values[0][i]<<" size "<< energy5OfMts[i].size() << std::endl;
    for(unsigned long j = 0; j < energy5OfMts[i].size(); j++){
      std::cout << energy5OfMts[i][j] << std::endl;
      for(unsigned long k =0; k < energyPdf5OfMts[i][j].size()/2; k++){
	std::cout << energyPdf5OfMts[i][j][2*k] <<"  "<< energyPdf5OfMts[i][j][2*k + 1] <<"  "<< energyCdf5OfMts[i][j][2*k + 1]<< std::endl;
      }
    }
  }
  */
}

TNudyEndfEnergy::~TNudyEndfEnergy(){}

double TNudyEndfEnergy::recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2){
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

double TNudyEndfEnergy::recursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
  double gx=0.0;	       
  double pe = TNudyCore::Instance()->Interpolate(nbt1,int1,NR,fE1,fP1,NP,energy);
  double thetae = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,fP2,np2, energy);
  if(thetae > 0.0)gx = TNudyCore::Instance()->Interpolate(nbt3,int3,nr3,fE3,fP3,np3, mid/thetae);
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

double TNudyEndfEnergy::recursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
  double pe = TNudyCore::Instance()->Interpolate(nbt1,int1,NR,fE1,fP1,NP,energy);
  double I = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,INorm,np2, energy);
  double thetae = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,fP2,np2, energy);
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

double TNudyEndfEnergy::recursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy){
  double pdf =1.0;
  double mid     = 0.5 * (x1 + x2);
  if((pdf1==0.0 && pdf2 ==0.0) || x1==x2)return 0;
  double pe = TNudyCore::Instance()->Interpolate(nbt1,int1,NR,fE1,fP1,NP,energy);
  double I = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,INorm,np2, energy);
  double ae = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE2,fP2,np2, energy);
  double be = TNudyCore::Instance()->Interpolate(nbt2,int2,nr2,fE3,fP3,np3, energy);
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
void TNudyEndfEnergy::fillPdf1d(){
  TNudyCore::Instance()->Sort(energyFile5, energyPdfFile5);
  TNudyCore::Instance()->ThinningDuplicate(energyFile5, energyPdfFile5);	  
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
}