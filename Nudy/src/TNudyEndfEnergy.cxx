// 	This class is reconstructing probability tables for Energy distribution 
// 	of the secondatries
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 22, 2016

#include "TNudyEndfEnergy.h"
class TNudyEndfFile;
class TNudyEndfSec;
class TNudyEndfCont;
class TNudyEndfTab1;
class TNudyEndfTab2;
class TNudyEndfList;

TNudyEndfEnergy::TNudyEndfEnergy(){}

//______________________________________________________________________________
TNudyEndfEnergy::TNudyEndfEnergy(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  int size1, size2, size3;
  std::vector<double>ein, cdf, pdf;
  std::vector<std::vector<double> >cdf2d, pdf2d;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    for (int k = 0; k < sec->GetN1(); k++) {
      TIter recIter(sec->GetRecords());
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      int MT = sec->GetMT();
      MtNumbers.push_back(MT);
      int LF = tab1->GetL2();
//      std::cout<<" LF = "<< LF <<" MT "<< MT  << std::endl;
      NR = tab1->GetN1();
      NP = tab1->GetN2();
//****************************************************************************	
      // arbitrary tabulated function
      if(LF==1){
	NBT         = new int[NR]();
	INT1        = new int[NR]();
	fE1_file5    = new double[NP]();
	fp91_file5   = new double[NP]();
	for(int cr=0; cr < NR ; cr ++){
	  NBT[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1_file5[crs]  = tab1->GetX(crs);
	  fp91_file5[crs] = tab1->GetY(crs);
	}
	TNudyEndfTab2 *tab2 = (TNudyEndfTab2 *)recIter.Next();
	NR2 = tab2->GetN1();
	NE  = tab2->GetN2();
	NBT2         = new int[NR2]();
	INT2        = new int[NR2]();
	for(int cr=0; cr < NR2 ; cr ++){
	  NBT2[cr] = tab2->GetNBT(cr);
	  INT2[cr] = tab2->GetINT(cr);
	}
        lCoef = new TArrayD[NE];
	for(int cr=0; cr < NE ; cr ++){
	  TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
	  NR3 = tab12->GetN1();
	  NE2  = tab12->GetN2();
	  ein.push_back(tab12->GetC2());
	  NBT3         = new int[NR3]();
	  INT3        = new int[NR3]();
	  fE3_file5    = new double[NE2]();
	  fp93_file5   = new double[NE2]();
	  for(int cr=0; cr < NR3 ; cr ++){
	    NBT3[cr] = tab12->GetNBT(cr);
	    INT3[cr] = tab12->GetINT(cr);
	  }
	  for(int crs=0; crs < NE2 ; crs ++){
	    fE3_file5[crs]  = tab12->GetX(crs);
	    fp93_file5[crs] = tab12->GetY(crs);
	    energyFile5.push_back(tab12->GetX(crs));
	    energyPdfFile5.push_back(tab12->GetY(crs));
	  }
	  for(int cr=0; cr < NE2 - 1 ; cr ++){
	    recursionLinearFile5Prob(fE3_file5[cr], fE3_file5[cr+1], fp93_file5[cr], fp93_file5[cr+1]);
	  }
	  Sort(energyFile5, energyPdfFile5);
	  ThinningDuplicate(energyFile5, energyPdfFile5);
	  cdfGenerateE(energyFile5, energyPdfFile5);
	  size1 = energyFile5.size();	  
	  for(int i = 0; i < size1; i++){
	    pdf.push_back(energyFile5[i]);
	    pdf.push_back(energyPdfFile5[i]);
	    cdf.push_back(energyFile5[i]);
	    cdf.push_back(energyCdfFile5[i]);	  
	  }
	  pdf2d.push_back(pdf);
	  cdf2d.push_back(cdf);
	  energyFile5.clear();	
	  energyPdfFile5.clear();	
	  energyCdfFile5.clear();	
	  pdf.clear();
	  cdf.clear();
	}
	energy5OfMts.push_back(ein);
	energyPdf5OfMts.push_back(pdf2d);
	energyCdf5OfMts.push_back(cdf2d);
	ein.clear();
	pdf2d.clear();
	cdf2d.clear();
//****************************************************************************	
      // general evaporation spectrum
      }else if(LF==5){
	double u = tab1->GetC1();
	NBT         = new int[NR]();
	INT1        = new int[NR]();
	fE1_file5    = new double[NP]();
	fp91_file5   = new double[NP]();
	for(int cr=0; cr < NR ; cr ++){
	  NBT[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1_file5[crs]  = tab1->GetX(crs);
	  fp91_file5[crs] = tab1->GetY(crs);
	}
        TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
	NR2 = tab11->GetN1();
	NE  = tab11->GetN2();
	NBT2         = new int[NR2]();
	INT2        = new int[NR2]();
	fE2_file5    = new double[NE]();
	fp92_file5   = new double[NE]();
	for(int cr=0; cr < NR2 ; cr ++){
	  NBT2[cr] = tab11->GetNBT(cr);
	  INT2[cr] = tab11->GetINT(cr);
	}
	for(int crs=0; crs < NE ; crs ++){
          fE2_file5[crs]  = tab11->GetX(crs);
          fp92_file5[crs] = tab11->GetY(crs);
	  ein.push_back(tab11->GetX(crs));
	}
        TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
	NR3 = tab12->GetN1();
	NE2  = tab12->GetN2();
	NBT3         = new int[NR3]();
	INT3        = new int[NR3]();
	fE3_file5    = new double[NE2]();
	fp93_file5   = new double[NE2]();
	for(int cr=0; cr < NR3 ; cr ++){
	  NBT3[cr] = tab12->GetNBT(cr);
	  INT3[cr] = tab12->GetINT(cr);
	}
	for(int crs=0; crs < NE2 ; crs ++){
          fE3_file5[crs]  = tab12->GetX(crs);
          fp93_file5[crs] = tab12->GetY(crs);
	}
	double energy, eout;
	for(int i = 0; i < NE; i++){
	  energy = fE2_file5[i];
	  eout =1E-5;
	  do
	  {
	    double gx=0.0;	       
	    double pe = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE1_file5,fp91_file5,NP,energy);
	    double thetae = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,fp92_file5,NE, energy);
	    if(thetae > 0.0)gx = TNudyCore::Instance()->Interpolate(NBT3,INT3,NR3,fE3_file5,fp93_file5,NE2, eout/thetae);
	    double ppe = 0.0;
	    ppe = pe * gx;
	    energyFile5.push_back(eout);
	    energyPdfFile5.push_back(ppe);
 //   std:: cout<<"k = "<< k <<" pe "<< pe <<" thetae "<< thetae <<"  prob "<< ppe <<"  ein "<< energy <<"  eout "<< eout << std::endl;
	    eout *= 2;
	  }while(eout < energy - u);
	  int size = energyFile5.size();
	  for(int cr=0; cr < size - 1 ; cr ++){
	    recursionLinearFile5GenEva(energyFile5[cr], energyFile5[cr+1], energyPdfFile5[cr], energyPdfFile5[cr+1],energy);
	  }
	  Sort(energyFile5, energyPdfFile5);
	  ThinningDuplicate(energyFile5, energyPdfFile5);
	  cdfGenerateE(energyFile5, energyPdfFile5);
	  size2 = energyFile5.size();	  
	  for(int i = 0; i < size2; i++){
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
//****************************************************************************
      // simple maxwellian fission spectrum
      }else if(LF==7){
	double u = tab1->GetC1();
	NBT         = new int[NR]();
	INT1        = new int[NR]();
	fE1_file5    = new double[NP]();
	fp91_file5   = new double[NP]();
	for(int cr=0; cr < NR ; cr ++){
	  NBT[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1_file5[crs]  = tab1->GetX(crs);
	  fp91_file5[crs] = tab1->GetY(crs);
	}
        TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
	NR2 = tab11->GetN1();
	NE  = tab11->GetN2();
	NBT2         = new int[NR2]();
	INT2        = new int[NR2]();
	fE2_file5    = new double[NE]();
	fp92_file5   = new double[NE]();
	INorm        = new double[NE]();
	for(int cr=0; cr < NR2 ; cr ++){
	  NBT2[cr] = tab11->GetNBT(cr);
	  INT2[cr] = tab11->GetINT(cr);
	}
	for(int crs=0; crs < NE ; crs ++){
          fE2_file5[crs]  = tab11->GetX(crs);
          fp92_file5[crs] = tab11->GetY(crs);
	  double sq = (fE2_file5[crs] - u)/fp92_file5[crs];
	  double sqt = sqrt(sq);
          INorm[crs] = pow(fp92_file5[crs], 1.5)*(0.5*1.7724538529055*erf(sqt)- sqt*exp(-sq));
	  ein.push_back(tab11->GetX(crs));
	}
	double energy, eout;
	for(int i = 0; i < NE; i++){
	  energy = fE2_file5[i];
	  eout =1E-5;
	  do
	  {
	    double pe = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE1_file5,fp91_file5,NP,energy);
	    double I = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,INorm,NE, energy);
	    double thetae = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,fp92_file5,NE, energy);
	    double ppe = 0.0;
	    if(I > 0.0)ppe = pe * sqrt(eout) * exp (-eout/thetae)/I;
	    energyFile5.push_back(eout);
	    energyPdfFile5.push_back(ppe);
//    std:: cout<<"I = "<< I <<" pe "<< pe <<" thetae "<< thetae <<"  prob "<< ppe <<"  ein "<< energy <<"  eout "<< eout << std::endl;
	    eout *= 2;
	  }while(eout < energy - u);
	  int size = energyFile5.size();
	  for(int cr=0; cr < size - 1 ; cr ++){
	    recursionLinearFile5Maxwell(energyFile5[cr], energyFile5[cr+1], energyPdfFile5[cr], energyPdfFile5[cr+1],energy);
	  }
	  Sort(energyFile5, energyPdfFile5);
	  ThinningDuplicate(energyFile5, energyPdfFile5);	  
	  cdfGenerateE(energyFile5, energyPdfFile5);
	  size3 = energyFile5.size();
	  for(int i = 0; i < size3; i++){
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
///////////////////////////////////////////////////////////////////////////////////
      // evaporation spectrum
      }else if(LF==9){
	double u = tab1->GetC1();
	NBT         = new int[NR]();
	INT1        = new int[NR]();
	fE1_file5    = new double[NP]();
	fp91_file5   = new double[NP]();
	for(int cr=0; cr < NR ; cr ++){
	  NBT[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1_file5[crs]  = tab1->GetX(crs);
	  fp91_file5[crs] = tab1->GetY(crs);
	}
	TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
	NR2 = tab11->GetN1();
	NE  = tab11->GetN2();
	NBT2         = new int[NR2]();
	INT2        = new int[NR2]();
	fE2_file5    = new double[NE]();
	fp92_file5   = new double[NE]();
	INorm        = new double[NE]();
      
	for(int cr=0; cr < NR2 ; cr ++){
	  NBT2[cr] = tab11->GetNBT(cr);
	  INT2[cr] = tab11->GetINT(cr);
	}
	for(int crs=0; crs < NE ; crs ++){
	  fE2_file5[crs]  = tab11->GetX(crs);
	  fp92_file5[crs] = tab11->GetY(crs);
	  INorm[crs] = fp92_file5[crs]*fp92_file5[crs]*(1. - exp(-(fE2_file5[crs] - u)/fp92_file5[crs])*
		      (1.0 + (fE2_file5[crs] - u)/fp92_file5[crs]));
	  ein.push_back(tab11->GetX(crs));
	}
	double energy, eout;
	for(int i = 0; i < NE; i++){
	  energy = fE2_file5[i];
	  eout =1E-5;
	  do
	  {
	    double pe = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE1_file5,fp91_file5,NP,energy);
	    double I = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,INorm,NE, energy);
	    double thetae = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,fp92_file5,NE, energy);
	    double ppe = 0.0;
	    if(I > 0.0)ppe = pe * eout * exp (-eout/thetae)/I;
	    energyFile5.push_back(eout);
	    energyPdfFile5.push_back(ppe);
  //    std:: cout<<"I = "<< I <<" pe "<< pe <<" thetae "<< thetae <<"  prob "<< ppe <<"  ein "<< energy <<"  eout "<< eout << std::endl;
	    eout *= 2;
	  }while(eout < u);
	  int size = energyFile5.size();
	  for(int cr=0; cr < size - 1 ; cr ++){
	    recursionLinearFile5Maxwell(energyFile5[cr], energyFile5[cr+1], energyPdfFile5[cr], energyPdfFile5[cr+1],energy);
	  }
	  Sort(energyFile5, energyPdfFile5);
	  ThinningDuplicate(energyFile5, energyPdfFile5);	  
	  cdfGenerateE(energyFile5, energyPdfFile5); 
	  size1 = energyFile5.size();
	  for(int i = 0; i < size1; i++){
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
      //std:: cout<< std::endl;
/////////////////////////////////////////////////////////////////////////////////////////
      // energy dependent watt spectrum
      }else if(LF==11){
	double u = tab1->GetC1();
	NBT         = new int[NR]();
	INT1        = new int[NR]();
	fE1_file5    = new double[NP]();
	fp91_file5   = new double[NP]();
	
	for(int cr=0; cr < NR ; cr ++){
	  NBT[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1_file5[crs]  = tab1->GetX(crs);
	  fp91_file5[crs] = tab1->GetY(crs);
	}
	TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
	NR2 = tab11->GetN1();
	NE  = tab11->GetN2();
	NBT2         = new int[NR2]();
	INT2        = new int[NR2]();
	fE2_file5    = new double[NE]();
	fp92_file5   = new double[NE]();
	INorm        = new double[NE]();
	
	for(int cr=0; cr < NR2 ; cr ++){
	  NBT2[cr] = tab11->GetNBT(cr);
	  INT2[cr] = tab11->GetINT(cr);
	}
	for(int crs=0; crs < NE ; crs ++){
	  fE2_file5[crs]  = tab11->GetX(crs);
	  fp92_file5[crs] = tab11->GetY(crs);
	  ein.push_back(tab11->GetX(crs));
	}
	TNudyEndfTab1 *tab12 = (TNudyEndfTab1 *)recIter.Next();
	NR3 = tab12->GetN1();
	NE2  = tab12->GetN2();
	NBT3         = new int[NR3]();
	INT3        = new int[NR3]();
	fE3_file5    = new double[NE2]();
	fp93_file5   = new double[NE2]();
	
	for(int cr=0; cr < NR3 ; cr ++){
	  NBT3[cr] = tab12->GetNBT(cr);
	  INT3[cr] = tab12->GetINT(cr);
	}
	for(int crs=0; crs < NE2 ; crs ++){
	  fE3_file5[crs]  = tab12->GetX(crs);
	  fp93_file5[crs] = tab12->GetY(crs);
	      
	  double a = fp92_file5[crs];
	  double b = fp93_file5[crs];
	  double eua = (fE3_file5[crs] - u )/a;
	    
	  INorm[crs] = 0.5*sqrt(0.25*PI*a*a*a*b)*exp(0.25*a*b)*(erf(sqrt(eua)- sqrt(0.25*a*b)) + 
		      erf(sqrt(eua) + sqrt(0.25*a*b))) - a*exp(-eua)*sinh(sqrt(b*(fE3_file5[crs] - u)));
	}
	double energy, eout;
	for(int i = 0; i < NE; i++){
	  energy = fE2_file5[i];
	  eout =1E-5;
	  do
	  {
	    double pe = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE1_file5,fp91_file5,NP,energy);
	    double I = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,INorm,NE, energy);
	    double ae = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,fp92_file5,NE, energy);
	    double be = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE3_file5,fp93_file5,NE2, energy);
	    double ppe = 0.0;
	    if(I > 0.0)ppe = pe * sinh(sqrt(be*eout)) * exp (-eout/ae)/I;
	    energyFile5.push_back(eout);
	    energyPdfFile5.push_back(ppe);
  //    std:: cout<<"I = "<< I <<" pe "<< pe <<" ae "<< ae <<"  prob "<< ppe <<"  ein "<< energy <<"  eout "<< eout << std::endl;
	    eout *= 2;
	  }while(eout < energy - u);	   
	  int size = energyFile5.size();
	  for(int cr=0; cr < size - 1 ; cr ++){
	    recursionLinearFile5Watt(energyFile5[cr], energyFile5[cr+1], energyPdfFile5[cr], energyPdfFile5[cr+1],energy);
	  }
	  Sort(energyFile5, energyPdfFile5);
	  ThinningDuplicate(energyFile5, energyPdfFile5);	  
	  cdfGenerateE(energyFile5, energyPdfFile5);  
	  size2 = energyFile5.size();
	  for(int i = 0; i < size2; i++){
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
      }else if(LF==12){
	NBT         = new int[NR]();
	INT1        = new int[NR]();
	fE1_file5    = new double[NP]();
	fp91_file5   = new double[NP]();
	
	for(int cr=0; cr < NR ; cr ++){
	  NBT[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE1_file5[crs]  = tab1->GetX(crs);
	  fp91_file5[crs] = tab1->GetY(crs);
	}
	TNudyEndfTab1 *tab11 = (TNudyEndfTab1 *)recIter.Next();
	double efl = tab11->GetC1();
	double efh = tab11->GetC2();
	NR2 = tab11->GetN1();
	NE  = tab11->GetN2();
	NBT2         = new int[NR2]();
	INT2        = new int[NR2]();
	fE2_file5    = new double[NE]();
	fp92_file5   = new double[NE]();
	  
	for(int cr=0; cr < NR2 ; cr ++){
	  NBT2[cr] = tab11->GetNBT(cr);
	  INT2[cr] = tab11->GetINT(cr);
	}
	for(int crs=0; crs < NE ; crs ++){
	  fE2_file5[crs]  = tab11->GetX(crs);
	  fp92_file5[crs] = tab11->GetY(crs);
	  ein.push_back(tab11->GetX(crs));
	}
	  double energy, eout;
	  for(int i = 0; i < NE; i++){
	    energy = fE2_file5[i];
	    eout =1E-5;
	  do
	  {
	    double pe = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE1_file5,fp91_file5,NP,energy);
	    double tm = TNudyCore::Instance()->Interpolate(NBT2,INT2,NR2,fE2_file5,fp92_file5,NE, energy);
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
	  }while(eout < fE1_file5[NP-1]);   
	  cdfGenerateE(energyFile5, energyPdfFile5);  
	  size3 = energyFile5.size();
	  for(int i = 0; i < size3 ; i++){
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
      }
      energy5OfMts.push_back(ein);
      energyPdf5OfMts.push_back(pdf2d);
      energyCdf5OfMts.push_back(cdf2d);
      ein.clear();
      pdf2d.clear();
      cdf2d.clear();
      if (NBT) {    delete[] NBT;    NBT = 0;  }
      if (INT1) {    delete[] INT1;    INT1 = 0;  }
      if (fE1_file5) {    delete[] fE1_file5;    fE1_file5 = 0;  }
      if (fp91_file5) {    delete[] fp91_file5;    fp91_file5 = 0;  }
      if (NBT2) {    delete[] NBT2;    NBT2 = 0;  }
      if (INT2) {    delete[] INT2;    INT2 = 0;  }
      if (fE2_file5) {    delete[] fE2_file5;    fE2_file5 = 0;  }
      if (fp92_file5) {    delete[] fp92_file5;    fp92_file5 = 0;  }
      if (NBT3) {    delete[] NBT3;    NBT3 = 0;  }
      if (INT3) {    delete[] INT3;    INT3 = 0;  }
      if (fE3_file5) {    delete[] fE3_file5;    fE3_file5 = 0;  }
      if (fp93_file5) {    delete[] fp93_file5;    fp93_file5 = 0;  }
    }
  }
  Mt5Values.push_back(MtNumbers);
  MtNumbers.clear();
  size1 =  energy5OfMts.size();
  for(int i = 0; i < size1 ; i++){
    size2 =  energy5OfMts[i].size();
      std::cout <<" mt "<<Mt5Values[0][i]<<" size "<< energy5OfMts[i].size() << std::endl;
    for(int j = 0; j < size2; j++){
      size3 =  energyPdf5OfMts[i][j].size()/2;
      std::cout << energy5OfMts[i][j] << std::endl;
      for(int k =0; k < size3; k++){
	std::cout << energyPdf5OfMts[i][j][2*k] <<"  "<< energyPdf5OfMts[i][j][2*k + 1] <<"  "<< energyCdf5OfMts[i][j][2*k + 1]<< std::endl;
      }
    }
  }
}

TNudyEndfEnergy::~TNudyEndfEnergy(){}

double TNudyEndfEnergy::ThinningDuplicate(std::vector<double> &x1){
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
double TNudyEndfEnergy::ThinningDuplicate(std::vector<double> &x1,std::vector<double> &x2){
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

void TNudyEndfEnergy::cdfGenerateE(std::vector<double> &x1,std::vector<double> &x2){
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

double TNudyEndfEnergy::recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2){
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

double TNudyEndfEnergy::recursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy){
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

double TNudyEndfEnergy::recursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy){
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

double TNudyEndfEnergy::recursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy){
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
void TNudyEndfEnergy::Sort(std::vector<double>& x1, std::vector<double>& x2){
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
