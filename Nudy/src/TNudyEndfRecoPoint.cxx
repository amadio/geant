// This class is reconstructing ENDF data and creating probability tables for the angle
// and energy distributions of the secondatries
// Author: Dr. Harphool Kumawat
// Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// date of creation: March 22, 2016

#include "TNudyEndfDoppler.h"
#include "TNudyEndfAng.h"
#include "TNudyEndfRecoPoint.h"

TNudyEndfRecoPoint::TNudyEndfRecoPoint(): TObject(){}

double TNudyEndfRecoPoint::recursionLinearPh(double x1, double x2, double sig1, double sig2){
  double siga;
  double mid     = 0.5 * (x1 + x2);
  if((sig1==0.0 && sig2 ==0.0) || x1==x2 || x1 < 1E-5 || x2 < 1E-5)return 0;
//  std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< sig1 <<"  "<< sig2 << std::endl;
  siga = TNudyCore::Instance()->Interpolate(NBT1, INT1, NR, fE_file1,fph_file1, NP, mid);
//  std::cout<< mid <<" mid  "<< siga1 <<std::endl;
//  double sigmid1 = TNudyCore::Instance()->LinearInterpolation(x1,sig1,x2,sig2,mid);
  double sigmid1 = sig1 + (sig2 - sig1)*(mid - x1)/(x2 - x1);
//  std::cout << mid <<" linear "<< siga <<"  "<< sigmid1 << std::endl;
  
  if(fabs((siga - sigmid1)/sigmid1)<=sigDiff){
    return 0;
  }
  einFile1.push_back(mid); 
  phFile1.push_back(siga);
  recursionLinearPh(x1, mid, sig1, siga);
  recursionLinearPh(mid, x2, siga, sig2);
  return 0;  
}

double TNudyEndfRecoPoint::recursionLinearNu(double x1, double x2, double sig1, double sig2){
  double siga;
  double mid     = 0.5 * (x1 + x2);
  if((sig1==0.0 && sig2 ==0.0) || x1==x2 || x1 < 1E-5 || x2 < 1E-5)return 0;
//  std::cout <<" beg   "<< x1 <<"  "<< x2 <<"  "<< sig1 <<"  "<< sig2 << std::endl;
  siga = TNudyCore::Instance()->Interpolate(NBT1, INT1, NR, fE_file1,fnu_file1, NP, mid);
//  std::cout<< mid <<" mid  "<< siga1 <<std::endl;
//  double sigmid1 = TNudyCore::Instance()->LinearInterpolation(x1,sig1,x2,sig2,mid);
  double sigmid1 = sig1 + (sig2 - sig1)*(mid - x1)/(x2 - x1);
//  std::cout << mid <<" linear "<< siga <<"  "<< sigmid1 << std::endl;
  
  if(fabs((siga - sigmid1)/sigmid1)<=sigDiff){
    return 0;
  }
  einFile1.push_back(mid); 
  nuFile1.push_back(siga);
  recursionLinearNu(x1, mid, sig1, siga);
  recursionLinearNu(mid, x2, siga, sig2);
  return 0;  
}	

void TNudyEndfRecoPoint::ReadFile1(TNudyEndfFile *file) {
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
//    std::cout <<"file  "<<sec->GetMT() <<std::endl;
    TIter recIter(sec->GetRecords());
    //double ZA   = sec->GetC1();
    //double AWR  = sec->GetC2();
    int MT = sec->GetMT();
    if(MT == 452){// Total fission neutron multiplicity
      int LNU = sec->GetL2();
      if(LNU == 1){
//      std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
	TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
	TArrayD cnc(list->GetN1());
//	      std::cout<<list->GetNPL()<<"  "<< list->GetN1() << std::endl;
        for (int j = 0; j < list->GetN1(); j++){
	  cnc[j] = list->GetLIST(j);
//	      std::cout<<cnc[j]<<" cnc "<< list->GetN1() << std::endl;
        }
	double ein = 1.0;
	do
	{
	double  nun = 0;
	for(int i = 0; i < list->GetN1(); i++){
	  nun += cnc[i]*pow(ein,i); 
	  }
	eintFile1.push_back(ein);
	nutFile1.push_back(nun);
	ein *= 2;
	}while(ein < 21E8);
	for(int cr=0; cr < eintFile1.size() ; cr ++){
//        std::cout <<"nu = "<< eintFile1[cr] <<"  "<< nutFile1[cr]  << std::endl;
	}	
	     
      } else {
//	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
      NR = tab1->GetN1();
      NP = tab1->GetN2();
      NBT1        = new int[NR]();
      INT1        = new int[NR]();
      fE_file1    = new double[NP]();
      fnu_file1 = new double[NP]();
//	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
      for(int cr=0; cr < NR ; cr ++){
	NBT1[cr] = tab1->GetNBT(cr);
	INT1[cr] = tab1->GetINT(cr);
      }
      for(int crs=0; crs < NP ; crs ++){
	fE_file1[crs]  = tab1->GetX(crs);
	fnu_file1[crs] = tab1->GetY(crs);
	eintFile1.push_back(fE_file1[crs]);
	nutFile1.push_back(fnu_file1[crs]);
//    std::cout<<fE_file1[crs]<<"  "<<fnu_file1[crs]<< std::endl;
      }
      for(int cr=0; cr < NP-1 ; cr ++){
	recursionLinearNu(fE_file1[cr], fE_file1[cr+1], fnu_file1[cr], fnu_file1[cr+1]);
      }
      Sort(eintFile1, nutFile1);
      for(int cr=0; cr < eintFile1.size() ; cr ++){
//        std::cout <<"nu = "<< einFile1[cr] <<"  "<< nuFile1[cr]  << std::endl;
      }	
//eintFile1.clear();	
//nutFile1.clear();	
      if (NBT1) { delete[] NBT1;  NBT1 = 0; }
      if (INT1) { delete[] INT1;  INT1 = 0;}
      if (fE_file1) { delete[] fE_file1; fE_file1 = 0; }
      if (fnu_file1) { delete[] fnu_file1; fnu_file1 = 0; }
    }
    } else if(MT == 455){// delayed neutron multiplicity
      int LDG   = sec->GetL1();
      int LNU   = sec->GetL2();
//	 std::cout<<" LNU "<< LNU <<"  LDG "<< LDG << std::endl;
      if(LNU == 1 && LDG == 0){
	
      }else if(LNU == 1 && LDG == 1){
	
      }else if(LNU == 2 && LDG == 0){
	TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
	int NNF = list->GetNPL();
	TArrayD nui(NNF);
	for(int i = 0; i < NNF; i++){
	  nui[i] = list->GetLIST(i); 
//	 std::cout<<" lambda "<< nui[i] << std::endl;
	}
        TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
        NR = tab1->GetN1();
        NP = tab1->GetN2();
	NBT1        = new int[NR]();
	INT1        = new int[NR]();
	fE_file1    = new double[NP]();
	fnu_file1 = new double[NP]();
//	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
	for(int cr=0; cr < NR ; cr ++){
	  NBT1[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE_file1[crs]  = tab1->GetX(crs);
	  fnu_file1[crs] = tab1->GetY(crs);
	  einFile1.push_back(fE_file1[crs]);
	  nuFile1.push_back(fnu_file1[crs]);
//    std::cout<<fE_file1[crs]<<"  "<<fnu_file1[crs]<< std::endl;
	}
	for(int cr=0; cr < NP-1 ; cr ++){
	  recursionLinearNu(fE_file1[cr], fE_file1[cr+1], fnu_file1[cr], fnu_file1[cr+1]);
	}
	Sort(einFile1, nuFile1);
	for(int cr=0; cr < einFile1.size() ; cr ++){
//        std::cout <<"nu = "<< einFile1[cr] <<"  "<< nuFile1[cr]  << std::endl;
	}	
	einFile1.clear();	
	nuFile1.clear();	
	
      }else if(LNU == 2 && LDG == 1){
	
	
      }      
    }else if(MT == 456){//prompt neutron multiplicity
	
      int LNU = sec->GetL2();
      if(LNU == 1){
//      std::cout<<"prompt nu = "<< ZA <<" LNU "<< LNU << std::endl;
	TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
	//double nup = list->GetLIST(0);
      }else {
//	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
	TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
	NR = tab1->GetN1();
	NP = tab1->GetN2();
	NBT1        = new int[NR]();
	INT1        = new int[NR]();
	fE_file1    = new double[NP]();
	fnu_file1 = new double[NP]();
//	std::cout<<"ZA = "<< ZA <<" LNU "<< LNU << std::endl;
	for(int cr=0; cr < NR ; cr ++){
	  NBT1[cr] = tab1->GetNBT(cr);
	  INT1[cr] = tab1->GetINT(cr);
	}
	for(int crs=0; crs < NP ; crs ++){
	  fE_file1[crs]  = tab1->GetX(crs);
	  fnu_file1[crs] = tab1->GetY(crs);
	  einFile1.push_back(fE_file1[crs]);
	  nuFile1.push_back(fnu_file1[crs]);
//    std::cout<<fE_file1[crs]<<"  "<<fnu_file1[crs]<< std::endl;
	}
	for(int cr=0; cr < NP-1 ; cr ++){
	  recursionLinearNu(fE_file1[cr], fE_file1[cr+1], fnu_file1[cr], fnu_file1[cr+1]);
	}
	Sort(einFile1, nuFile1);
	for(int cr=0; cr < einFile1.size() ; cr ++){
//        std::cout <<"nu = "<< einFile1[cr] <<"  "<< nuFile1[cr]  << std::endl;
	}	
	einFile1.clear();	
	nuFile1.clear();	
	if (NBT1) { delete[] NBT1;  NBT1 = 0; }
	if (INT1) { delete[] INT1;  INT1 = 0;}
	if (fE_file1) { delete[] fE_file1; fE_file1 = 0; }
	if (fnu_file1) { delete[] fnu_file1; fnu_file1 = 0; }
      }
    }else if(MT == 458){
      TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
      int NPLY = list->GetL2();
//	  std::cout<< "NPLY "<< NPLY << std::endl;
      if(NPLY == 0){
	double EFR = list->GetLIST(0);
	double ENP = list->GetLIST(2);
	double END = list->GetLIST(4);
	double EGP = list->GetLIST(6);
	double EGD = list->GetLIST(8);
	double EB  = list->GetLIST(10);
	double ENU = list->GetLIST(12);
	double ER  = list->GetLIST(14);
	//double ET  = list->GetLIST(16);
	double ein = 1.0;
	do
	{
	  double  EFis =0;
	  EFis = EFR + ENP + END + EGP + EGD + EB + ENU; 
	  EFis -= 0.100*ein;//nuetrino energy dependence
	  EFis -= 0.075*ein;//delayed gamma energy dependence
	  EFis -= 0.075*ein;//delayed beta energy dependence
//	      for(int i = 0; i < eintFile1.size(); i++)
//		std::cout<< eintFile1[i] <<"  "<<nutFile1[i] << std::endl; 
	  int n0 = 0;
	  int max = eintFile1.size() - 1;
	  int mid = 0;
	  if (ein <= eintFile1[n0]){
	    n0 = 0;
	  }else if (ein >= eintFile1[max]){
	    n0 = max - 1;
	  }else {
	    while (max - n0 > 1) {
	      mid = (n0 + max) / 2;
	      if (ein < eintFile1[mid])
		max = mid;
	      else
		n0 = mid;
//      std::cout<<"n0 "<< n0 <<" max "<< max << std::endl;
	      }
	  }
//  std::cout<< eintFile1[n0] <<std::endl;
	  double nue = TNudyCore::Instance()->LinearInterpolation(eintFile1[n0],nutFile1[n0],eintFile1[n0+1],nutFile1[n0+1], ein);
	  double nu0 = nutFile1[0];
//   std::cout<<"n0 "<< n0 <<" nu0 "<< nu0 << std::endl;
	  EFis -= -1.307*ein + 8.07*1E6*(nue - nu0);//prompt neutron energy dependence
//      std::cout<<"ein "<< ein <<"  "<< EFis <<std::endl;
	  ein *= 10;
	  }while(ein < 21E8);
//	std::cout<< "NPLY "<< NPLY <<" ER "<< ER <<" ET "<< ET << std::endl;
      }else {
	TArrayD c0(9*(NPLY+1)),c1(9*(NPLY+1));
	for(int i = 0; i < 9*(NPLY+1); i++){
	  c0[i] = list->GetLIST(i);
//	     std::cout <<"c0  "<< c0[i] << std::endl; 
	}
	for(int i = 9*(NPLY+1); i < 18*(NPLY+1); i++){
	  c1[i-9*(NPLY+1)] = list->GetLIST(i);
//	     std::cout <<"c1  "<< c1[i-9*(NPLY+1)] << std::endl;
	}
	double ein = 1.0;
	do
	{
	  double  EFis =0;
	  for(int i = 0; i < 9*NPLY/2 - 2; i++){
	    EFis += c0[i*2] + c1[i*2]*ein; 
//      std::cout<<"ein "<< ein <<"  "<< EFis <<std::endl;
	  }
	  ein *= 10;
	}while(ein < 21E8);
      }
	    
      }else if(MT == 460){
	
	int LO  = sec->GetL1();
	int NG  = sec->GetN1();
      
//      std::cout<<" Lo "<< LO << std::endl;
	if(LO == 1){
	  for(int ng = 0; ng < NG; ng++)  {
	    TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)recIter.Next();
	    NR = tab1->GetN1();
	    NP = tab1->GetN2();
	    NBT1        = new int[NR]();
	    INT1        = new int[NR]();
	    fE_file1    = new double[NP]();
	    fph_file1 = new double[NP]();
	    for(int cr=0; cr < NR ; cr ++){
	      NBT1[cr] = tab1->GetNBT(cr);
	      INT1[cr] = tab1->GetINT(cr);
  //	std::cout<<"NR = "<< NR <<" NP "<< NP << " NBT1 "<<NBT1[cr]<< " INT1 "<<INT1[cr] <<std::endl;
	    }
	    for(int crs=0; crs < NP ; crs ++){
	      fE_file1[crs]  = tab1->GetX(crs);
	      fph_file1[crs] = tab1->GetY(crs);
	      einFile1.push_back(fE_file1[crs]);
	      phFile1.push_back(fph_file1[crs]);
  //    std::cout<<fE_file1[crs]<<"  "<<fph_file1[crs]<< std::endl;
	    }
	    for(int cr=0; cr < NP-1 ; cr ++){
	      recursionLinearPh(fE_file1[cr], fE_file1[cr+1], fph_file1[cr], fph_file1[cr+1]);
	    }
	    Sort(einFile1, phFile1);
	    for(int cr=0; cr < einFile1.size() ; cr ++){
  //        std::cout <<"nu = "<< einFile1[cr] <<"  "<< nuFile1[cr]  << std::endl;
	    }
  //	std::cout<<" linear size "<< einFile1.size() << std::endl;
	    einFile1.clear();	
	    phFile1.clear();	
	    if (NBT1) { delete[] NBT1;  NBT1 = 0; }
	    if (INT1) { delete[] INT1;  INT1 = 0;}
	    if (fE_file1) { delete[] fE_file1; fE_file1 = 0; }
	    if (fph_file1) { delete[] fph_file1; fph_file1 = 0; }
	  }	 
	  
	}else{
          TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
	  int NNF = list->GetN1();
//	std::cout<< "Photon NNF "<< NNF <<" LO "<< LO << std::endl;
	  TArrayD lambda(NNF);
	  for(int i = 0; i < NNF; i++){
	    lambda[i] = list->GetLIST(i);
//	     std::cout <<"lambda  "<< lambda[i] << std::endl; 
	  }
	}
     }
  }
}
//______________________________________________________________________________
void TNudyEndfRecoPoint::ReadFile2(TNudyEndfFile *file) {
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  NLS = 0;double gjdeno = 1; LSSF=0;
  std::vector<double> rowI;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    for (int k = 0; k < sec->GetN1(); k++) {
      TIter recIter(sec->GetRecords());
      ZAI   = sec->GetC1();
      AWRI  = sec->GetC2();
      NIS   = sec->GetN1();
      TNudyEndfCont *cont1 = (TNudyEndfCont *)recIter.Next();
      ABN = cont1->GetC2();
      LFW = cont1->GetL2(); 
      NER = cont1->GetN1();
      std::cout<< "NER "<< NER <<" AWRI "<< AWRI << std::endl;
      for(int j = 0; j < NER; j++){
	TNudyEndfCont *cont2 = (TNudyEndfCont *)recIter.Next();
	eLo = cont2->GetC1();  
	eHi = cont2->GetC2();  
	LRU = cont2->GetL1(); 
	LRF = cont2->GetL2();  
	NRO = cont2->GetN1();  
	NAPS = cont2->GetN2();
	std::cout<<" LRU "<<LRU<<" LRF "<<LRF<< std::endl;
      std::cout<<"eLo "<<eLo<<" eHi  "<<eHi<<" NLS "<<NLS<<" LRU "<<LRU<<" LRF "<<LRF<< std::endl;
      
	if(LRU==2 && LRF==1 && LFW == 1){
	}else{
	  TNudyEndfCont *cont3 = (TNudyEndfCont *)recIter.Next();
	  SPI = cont3->GetC1(); 
	  AP = cont3->GetC2();
	  if(LRU==2)LSSF = cont3->GetL1();
	  //if(LRF==3)int LAD = cont3->GetL1();
	  if(LRU==1)NLS = cont3->GetN1();
	  NLS2 = cont3->GetN1();
	  NRS.resize(NLS, 0);
	  gjdeno = 4.0 * SPI + 2.0;
//       exit(1);
	}
	switch(LRU)
	{
	  case 0:
	  {
	  }break;  
	  case 1:
	  {//resolved resonance region
	  switch(LRF)
	  {
	    case 1: //single level resonance region
	    case 2://multi level resonance region
	    case 3://RM resonance region
	    case 4://Adler - Adler resonance region
	    {
	      flagResolve =1;
	      if(j==0)eLo1 = eLo; 
	      eHi1 = eHi;
	      for (int lseries = 0; lseries < NLS; lseries++) { 
		TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
		AWRI          = list->GetC1();  
		QX            = list->GetC2();
		if(QX == 0.0)APL[lseries] = AP;
		if(LRF==3 && QX >0.0)APL[lseries] = list->GetC2(); 
		l.push_back(list->GetL1()); 
		LRX           = list->GetL2();
		NRS[lseries]  = list->GetN2();
		A = Mn * AWRI; 
		Z = (int)(ZA - A)/1000.0;  
		rad_a = 0.08 + 0.123 * pow(1.00866491578*AWRI, (1./3.));
		if(AP==0.0)AP=rad_a;
//	std::cout << APL[lseries] <<"   "<< AP <<"  "<< rad_a <<std::endl;
		factor_k = kconst*(AWRI/(AWRI+1.0));
		JMIN = (std::fabs(SPI - l[lseries]) - 0.5);
		JMAX = (SPI + l[lseries] + 0.5);
		int nrsl0 = 0;
		nrsl0 = (!lseries) ? 0 :nrsl0 + NRS[lseries-1];  // last NRS value
		double jtemp[NRS[lseries]];
//      std::cout<<nrsl0<<"  "<< jtemp[NRS[lseries]] <<"  "<< NRS[lseries] << std::endl;
		if(LRF == 1 || LRF == 2){
		  for (int ii = 0; ii < NRS[lseries]; ii++) {
		  cueMat = nrsl0 + ii;
		  Er.push_back(list->GetLIST(ii*6+0));
		  J.push_back(list->GetLIST(ii*6+1));//J value
		  Gamma_r.push_back(list->GetLIST(ii*6+2));//total width
		  Gamma_n.push_back(list->GetLIST(ii*6+3));//neutron width
		  Gamma_g.push_back(list->GetLIST(ii*6+4));//gamma width
		  Gamma_f.push_back(list->GetLIST(ii*6+5));//fission width	
		  GJ.push_back((2.0 * std::fabs(J[cueMat])+1.0)/gjdeno); //(2J+1)/2(2I+1)
		  PhiEr.push_back(calcPene(std::fabs(list->GetLIST(ii*6+0)), lseries));
		  ShiftEr.push_back(calcShift(std::fabs(list->GetLIST(ii*6+0)), lseries));
		  Gamma_n[nrsl0 + ii] = Gamma_n[nrsl0 + ii]/PhiEr[nrsl0 + ii];
		  jtemp[ii] = list->GetLIST(ii*6+1);
//       std::cout<<Er[nrsl0 + ii]<<"  "<<PhiEr[nrsl0 + ii]<< " "<< ShiftEr[cueMat]<<"  "<<cueMat<<"  "<<nrsl0<<std::endl;
//	std::cout<<Gamma_r[ii]<<" Gamma_n  "<<Gamma_n[ii]<<" Gamma_g "<<Gamma_g[ii]<<" Gamma_f "<<Gamma_f[ii]<<std::endl;
//	std::cout<<Er[ii]<<" cueMat  "<<cueMat<<" J "<<J[ii]<<" gjdeno "<<gjdeno<<std::endl;
		}
		double gjfound = 0.0, jmis = -99, temp,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9;
		for(int isort = 0; isort < NRS[lseries]; isort++){
		  for(int isort1 = 1; isort1 < NRS[lseries]; isort1++){
		    if(jtemp[isort1] < jtemp[isort1-1]){
		      temp = jtemp[isort1];
		      jtemp[isort1] = jtemp[isort1-1];
		      jtemp[isort1-1] = temp;
	  
		      temp1 = Er[nrsl0 + isort1];
		      Er[nrsl0 + isort1] = Er[nrsl0 + isort1-1];
		      Er[nrsl0 + isort1-1] = temp1;
	  
		      temp2               = J[nrsl0 + isort1];
		      J[nrsl0 + isort1]   = J[nrsl0 + isort1-1];
		      J[nrsl0 + isort1-1] = temp2;
	  
		      temp3                     = Gamma_r[nrsl0 + isort1];
		      Gamma_r[nrsl0 + isort1]   = Gamma_r[nrsl0 + isort1-1];
		      Gamma_r[nrsl0 + isort1-1] = temp3;
	  
		      temp4                     = Gamma_n[nrsl0 + isort1];
		      Gamma_n[nrsl0 + isort1]   = Gamma_n[nrsl0 + isort1-1];
		      Gamma_n[nrsl0 + isort1-1] = temp4;
	  
		      temp5                     = Gamma_g[nrsl0 + isort1];
		      Gamma_g[nrsl0 + isort1]   = Gamma_g[nrsl0 + isort1-1];
		      Gamma_g[nrsl0 + isort1-1] = temp5;
		      
		      temp6                     = Gamma_f[nrsl0 + isort1];
		      Gamma_f[nrsl0 + isort1]   = Gamma_f[nrsl0 + isort1-1];
		      Gamma_f[nrsl0 + isort1-1] = temp6;
		      
		      temp7                     = GJ[nrsl0 + isort1];
		      GJ[nrsl0 + isort1]        = GJ[nrsl0 + isort1-1];
		      GJ[nrsl0 + isort1-1]      = temp7;
		      
		      temp8                     = PhiEr[nrsl0 + isort1];
		      PhiEr[nrsl0 + isort1]     = PhiEr[nrsl0 + isort1-1];
		      PhiEr[nrsl0 + isort1-1]   = temp8;
		      
		      temp9                     = ShiftEr[nrsl0 + isort1];
		      ShiftEr[nrsl0 + isort1]   = ShiftEr[nrsl0 + isort1-1];
		      ShiftEr[nrsl0 + isort1-1] = temp9;
		    }
		  }
		}
		int ju = 0;
		NJValue[l[lseries]] = 0; MisGj[l[lseries]] = 0;
		for(int j1 = 0; j1 < NRS[lseries]; j1++ ){
		  if(jtemp[j1] != jmis){
		    MissingJ[l[lseries]][ju] = jtemp[j1];
		    NJValue[l[lseries]] += 1 ; 
		    jmis = jtemp[j1];
		    gjfound += (2 *fabs(jtemp[j1]) + 1 )/gjdeno;	
		    ju += 1;
		  }
		}
		MisGj[l[lseries]] = 2*l[lseries] +1 - gjfound;
		} else if (LRF == 3){
		  for (int ii = 0; ii < NRS[lseries]; ii++) {
		    cueMat = nrsl0 + ii;
		    Er.push_back(list->GetLIST(ii*6+0));
		    J.push_back(list->GetLIST(ii*6+1));//J value
		    Gamma_n.push_back(list->GetLIST(ii*6+2));//total width
		    Gamma_g.push_back(list->GetLIST(ii*6+3));//neutron width
		    Gamma_fa.push_back(list->GetLIST(ii*6+4));//gamma width
		    Gamma_fb.push_back(list->GetLIST(ii*6+5));//fission width	
		    Gamma_fasq.push_back(sqrt(0.5*std::fabs(list->GetLIST(ii*6+4))));//gamma width
		    Gamma_fbsq.push_back(sqrt(0.5*std::fabs(list->GetLIST(ii*6+5))));//fission width
		    if(Gamma_fa[cueMat]<0.0)Gamma_fasq[cueMat] = - Gamma_fasq[cueMat];
		    if(Gamma_fb[cueMat]<0.0)Gamma_fbsq[cueMat] = - Gamma_fbsq[cueMat];
		    GJ.push_back((2.0 * std::fabs(J[cueMat])+1.0)/gjdeno); //(2J+1)/2(2I+1)
		    jtemp[ii] = list->GetLIST(ii*6+1);
		    PhiEr.push_back(calcPene(std::fabs(list->GetLIST(ii*6+0)), lseries));
		    ShiftEr.push_back(calcShift(std::fabs(list->GetLIST(ii*6+0)), lseries));
		    Gamma_n[nrsl0 + ii] = Gamma_n[nrsl0 + ii]/PhiEr[nrsl0 + ii];
    //       std::cout<<lseries<<"  "<<Er[nrsl0 + ii]<<"  "<<PhiEr[nrsl0 + ii]<< " "<< ShiftEr[cueMat]<<"  "<<cueMat<<"  "<<nrsl0<<std::endl;
    //	       std::cout<< GJ[ii] <<" J " << std::fabs(J[cueMat]) <<"   "<< Er[ii] << std::endl;
		  }
    //	std::cout << nrsl0 << std::endl;
		  double gjfound = 0.0, jmis = -99, temp,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11;
		  for(int isort = 0; isort < NRS[lseries]; isort++){
		    for(int isort1 = 1; isort1 < NRS[lseries]; isort1++){
		      if(jtemp[isort1] < jtemp[isort1-1]){
			temp = jtemp[isort1];
			jtemp[isort1] = jtemp[isort1-1];
			jtemp[isort1-1] = temp;
			    
			temp1 = Er[nrsl0 + isort1];
			Er[nrsl0 + isort1]   = Er[nrsl0 + isort1-1];
			  Er[nrsl0 + isort1-1] = temp1;
			
			temp2               = J[nrsl0 + isort1];
			J[nrsl0 + isort1]   = J[nrsl0 + isort1-1];
			J[nrsl0 + isort1-1] = temp2;
			
			temp3                      = Gamma_fa[nrsl0 + isort1];
			Gamma_fa[nrsl0 + isort1]   = Gamma_fa[nrsl0 + isort1-1];
			Gamma_fa[nrsl0 + isort1-1] = temp3;
			
			temp4                      = Gamma_n[nrsl0 + isort1];
			Gamma_n[nrsl0 + isort1]    = Gamma_n[nrsl0 + isort1-1];
			Gamma_n[nrsl0 + isort1-1]  = temp4;
			
			temp5                      = Gamma_g[nrsl0 + isort1];
			Gamma_g[nrsl0 + isort1]    = Gamma_g[nrsl0 + isort1-1];
			Gamma_g[nrsl0 + isort1-1]  = temp5;
			
			temp6                      = Gamma_fb[nrsl0 + isort1];
			Gamma_fb[nrsl0 + isort1]   = Gamma_fb[nrsl0 + isort1-1];
			Gamma_fb[nrsl0 + isort1-1] = temp6;
			
			temp7                      = GJ[nrsl0 + isort1];
			GJ[nrsl0 + isort1]         = GJ[nrsl0 + isort1-1];
			GJ[nrsl0 + isort1-1]       = temp7;
			
			temp8                      = PhiEr[nrsl0 + isort1];
			PhiEr[nrsl0 + isort1]      = PhiEr[nrsl0 + isort1-1];
			PhiEr[nrsl0 + isort1-1]    = temp8;
			
			temp9                      = ShiftEr[nrsl0 + isort1];
			ShiftEr[nrsl0 + isort1]    = ShiftEr[nrsl0 + isort1-1];
			ShiftEr[nrsl0 + isort1-1]  = temp9;
			
			temp10                       = Gamma_fasq[nrsl0 + isort1];
			Gamma_fasq[nrsl0 + isort1]   = Gamma_fasq[nrsl0 + isort1-1];
			Gamma_fasq[nrsl0 + isort1-1] = temp10;
			
			temp11                       = Gamma_fbsq[nrsl0 + isort1];
			Gamma_fbsq[nrsl0 + isort1]   = Gamma_fbsq[nrsl0 + isort1-1];
			Gamma_fbsq[nrsl0 + isort1-1] = temp11;
		      }
		    }
		  }
		  int ju = 0;
		  NJValue[l[lseries]] = 0; 
		  MisGj[l[lseries]] = 0; 
		  for(int j1 = 0; j1 < NRS[lseries]; j1++ ){
		    if(jtemp[j1] != jmis){
		      MissingJ[l[lseries]][ju] = jtemp[j1];
		      NJValue[l[lseries]] += 1 ; 
		      jmis = jtemp[j1];
		      gjfound += (2 *fabs(jtemp[j1]) + 1 )/gjdeno;	
		      ju += 1;
		    }
		  }
		  MisGj[l[lseries]] = 2*l[lseries] +1 - gjfound;
    //	      std::cout << MisGj[l[lseries]] << std::endl;
		} else if (LRF == 4){//Adler -Adler resonance region
		  for (int ii = 0; ii < NRS[lseries]; ii++) { // line 6 onwards data
		    cueMat = nrsl0 + ii;
		    at1.push_back(list->GetLIST(ii*6+0));  
		    at2.push_back(list->GetLIST(ii*6+1));  
		    at3.push_back(list->GetLIST(ii*6+2));  
		    at4.push_back(list->GetLIST(ii*6+3));  
		    bt1.push_back(list->GetLIST(ii*6+4));  
		    bt2.push_back(list->GetLIST(ii*6+5));  
		  }
		  TNudyEndfList *list1 = (TNudyEndfList *)recIter.Next();
		  for (int ii = 0; ii < list1->GetN2(); ii++) { 
		    det1.push_back(list1->GetLIST(ii*6+0));  
		    dwt1.push_back(list1->GetLIST(ii*6+1));  
		    grt1.push_back(list1->GetLIST(ii*6+2));  
		    git1.push_back(list1->GetLIST(ii*6+3));  
		    def1.push_back(list1->GetLIST(ii*6+4));  
		    dwf1.push_back(list1->GetLIST(ii*6+5));  
		    grf1.push_back(list1->GetLIST(ii*6+6+0));  
		    gif1.push_back(list1->GetLIST(ii*6+6+1));  
		    dec1.push_back(list1->GetLIST(ii*6+6+2));  
		    dwc1.push_back(list1->GetLIST(ii*6+6+3));  
		    grc1.push_back(list1->GetLIST(ii*6+6+4));  
		    gic1.push_back(list1->GetLIST(ii*6+6+5));  
		  }	      
		}
	      }// loop for L values
	    }break;
	    case 7:
	    { //R-Matrix
	    }break;
	  }//break;
//	  std::cout <<" linear lru "<< LRU <<" LRF "<< LRF << std::endl;
	  Linearize(j);
// clearing vectors for different energy regions  
	  NRS.clear();
	  l.clear();
	  Er.clear();
	  J.clear();
	  Gamma_r.clear();
	  Gamma_n.clear();
	  Gamma_g.clear();
	  Gamma_f.clear();
	  GJ.clear();
	  PhiEr.clear();
	  ShiftEr.clear();
	  Gamma_fa.clear();
	  Gamma_fasq.clear();
	  Gamma_fb.clear();
	  Gamma_fbsq.clear();
	  at1.clear();
	  at2.clear();
	  at3.clear();
	  at4.clear();
	  bt1.clear();
	  bt2.clear();
	  det1.clear();
	  dwt1.clear();
	  grt1.clear();
	  git1.clear();
	  def1.clear();
	  dwf1.clear();
	  grf1.clear();
	  gif1.clear();
	  dec1.clear();
	  dwc1.clear();
	  grc1.clear();
	  gic1.clear();
	  
	}break;
	case 2:
	{//unresolved resonance region
	flagUnResolve =1;
	std::cout<<" LRF "<< LRF <<" LFW "<< LFW <<" LSSF "<< LSSF << std::endl;
	switch(LRF)
	{
	case 2:
	{
	  for (int lseries = 0; lseries < NLS2; lseries++) {  
	    TNudyEndfList *list2 = (TNudyEndfList *)recIter.Next();
	    AWRI = list2->GetC1();  
	    l.push_back(list2->GetL1());
	    NJS = list2->GetN1();
	    NRJ.push_back(NJS);
	    APL[NLS+lseries] = AP;
//	std::cout<<lseries<<"  "<< AWRI <<"  "<< list2->GetN2() << std::endl; 
//	NRS[NLS+lseries] = list2->GetN2();
	    for (int jseries = 0; jseries < NJS; jseries++) {  
	      TNudyEndfList *list3 = (TNudyEndfList *)recIter.Next();
	      JSM.push_back(list3->GetN2()) ;
	      rad_a = 0.08 + 0.123 * pow(1.00866491578*AWRI, (1./3.));
	      if(AP==0.0)AP=rad_a;
	      factor_k = kconst*(AWRI/(AWRI+1.0));
	      JMIN = list3->GetC1(); 
	      INT = list3->GetL1();
	      amux.push_back(list3->GetLIST(2));//total width
	      amun.push_back(list3->GetLIST(3));//neutron width
	      amug.push_back(list3->GetLIST(4));//gamma width
	      amuf.push_back(list3->GetLIST(5));//fission width	
	      GJ.push_back((2.0 * JMIN +1.0)/gjdeno); //(2J+1)/2(2I+1)
	      for (int ii = 0; ii < list3->GetN2(); ii++) { // line 6 onwards data
		Es.push_back(list3->GetLIST((ii+1)*6+0));
		D.push_back(list3->GetLIST((ii+1)*6+1));//J value
		GX.push_back(list3->GetLIST((ii+1)*6+2));//total width
		GNO.push_back(list3->GetLIST((ii+1)*6+3));//neutron width
		GG.push_back(list3->GetLIST((ii+1)*6+4));//gamma width
		GF.push_back(list3->GetLIST((ii+1)*6+5));//fission width	
	      }	   
	    }//complete NJS reading
	  }//complete NLS2 reading
	  eLo2 = eLo; eHi2 = eHi;
	  std::cout << " eLo2 "<< eLo2 <<" eHi2 "<< eHi2 << std::endl;
	  if(LSSF==0)Linearize(j);
	  l.clear();
	  NRS.clear();
	  NRJ.clear();
	  JSM.clear();
	  amux.clear();
	  amun.clear();
	  amug.clear();
	  amuf.clear();
	  GJ.clear();
	  Es.clear();
	  D.clear();
	  GX.clear();
	  GNO.clear();
	  GG.clear();
	  GF.clear();
	}break;
	case 1:
	{
	switch(LFW)
	{
	case 0:
	{
	  for (int lseries = 0; lseries < NLS2; lseries++) {  
	    TNudyEndfList *list2 = (TNudyEndfList *)recIter.Next();
	    AWRI = list2->GetC1();  
	    l.push_back(list2->GetL1());
	    NJS = list2->GetN2();
	    NRJ.push_back(NJS);
	    APL[NLS+lseries] = AP;
//	std::cout<<" NJS "<< NJS <<std::endl;
	    for (int ii = 0; ii < NJS; ii++) { // line 6 onwards data
	      JSM.push_back(1) ;
	      rad_a = 0.08 + 0.123 * pow(1.00866491578*AWRI, (1./3.));
	      if(AP==0.0)AP=rad_a;
	      factor_k = kconst*(AWRI/(AWRI+1.0));
	      INT = 1;
	      Es.push_back(eLo);
	      D.push_back(list2->GetLIST((ii)*6+0));//Average level spacing for resonances with spin J
	      amun.push_back(list2->GetLIST((ii)*6+2));//degree of freedom neutron width
	      GNO.push_back(list2->GetLIST((ii)*6+3));//neutron width
	      GG.push_back(list2->GetLIST((ii)*6+4));//gamma width	
	      GF.push_back(0.0);//fission width	
	      GX.push_back(0.0);//fission width	
	      amug.push_back(0.0);//degree of freedom gamma width
	      amuf.push_back(0.0);//degree of freedom fission width	
	      amux.push_back(0.0);//degree of freedom competitive width
	      GJ.push_back((2.0 * list2->GetLIST((ii)*6+1) +1.0)/gjdeno); //(2J+1)/2(2I+1)
//	std::cout<<ii<<" D "<< NJS <<"  "<< list2->GetLIST((ii)*6+0)<<"  "<< D.size() << std::endl;
	      }	   
	    }
	    eLo2 = eLo; eHi2 = eHi;
	    if(LSSF==0)Linearize(j);
	    l.clear();
	    NRS.clear();
	    NRJ.clear();
	    JSM.clear();
	    amux.clear();
	    amun.clear();
	    amug.clear();
	    amuf.clear();
	    GJ.clear();
	    Es.clear();
	    D.clear();
	    GX.clear();
	    GNO.clear();
	    GG.clear();
	    GF.clear();
	  }break;
	  case 1:
	  {
	    TNudyEndfList *list = (TNudyEndfList *)recIter.Next();
	    SPI = list->GetC1(); 
	    AP = list->GetC2();
	    LSSF = list->GetL1();
	    NLS2 = list->GetN2(); 
	    NE = list->GetN1();
	    NRS.resize(NLS, 0);
	    gjdeno = 4.0 * SPI + 2.0;
	    for (int lseries = 0; lseries < NLS2; lseries++) {  
	      TNudyEndfCont *list2 = (TNudyEndfCont *)recIter.Next();
	      AWRI = list2->GetC1();  
	      l.push_back(list2->GetL1());
	      NJS = list2->GetN1();
	      NRJ.push_back(NJS);
	      APL[NLS+lseries] = AP;
//	 std::cout<<"njs "<< NJS << std::endl;
	      for (int jseries = 0; jseries < NJS; jseries++) {  
		TNudyEndfList *list3 = (TNudyEndfList *)recIter.Next();
		JSM.push_back(NE) ;
		rad_a = 0.08 + 0.123 * pow(1.00866491578*AWRI, (1./3.));
		if(AP==0.0)AP=rad_a;
		factor_k = kconst*(AWRI/(AWRI+1.0));
		INT = 1;
		amux.push_back(0.0);//total width
		amun.push_back(list3->GetLIST(2));//neutron width
		amug.push_back(0.0);//gamma width
		amuf.push_back(list3->GetL2());//fission width	
		GJ.push_back((2.0 * list3->GetLIST((1)) +1.0)/gjdeno); //(2J+1)/2(2I+1)
		for (int ii = 0; ii < NE; ii++) { // line 6 onwards data
		  Es.push_back(list->GetLIST(ii));
		  GF.push_back(list3->GetLIST(6+ii));//fission width	
		  GG.push_back(list3->GetLIST(4));//gamma width
		  GX.push_back(0.0);//gamma width
		  D.push_back(list3->GetLIST(0));//J value
		  GNO.push_back(list3->GetLIST(3));//neutron width
//	  std::cout<<list3->GetLIST(0)<<"  "<<list3->GetL2()<<"  "<<list3->GetLIST(6+ii)<<"  "<<list->GetLIST(ii)<<std::endl;
		}	   
	      }//complete NJS reading
	    }
	    eLo2 = eLo; eHi2 = eHi;
	    if(LSSF==0)Linearize(j);
	    l.clear();
	    NRS.clear();
	    NRJ.clear();
	    JSM.clear();
	    amux.clear();
	    amun.clear();
	    amug.clear();
	    amuf.clear();
	    GJ.clear();
	    Es.clear();
	    D.clear();
	    GX.clear();
	    GNO.clear();
	    GG.clear();
	    GF.clear();		
	  }break;
	  }
	  }break;
	  }
	  }break;
	}
      }
    }
  }
}
double TNudyEndfRecoPoint::recursionLinearFile3(double x1, double x2, double sig1, double sig2){
  double siga;int flag = 0;
  double mid     = 0.5 * (x1 + x2);
//  std::cout <<"beg "<< x1 <<"  "<< x2 <<"  "<< mid <<"  "<< sig1 <<"  "<< sig2 <<std::endl;
  if(sig1==0.0 && sig2 ==0.0)return 0;
  if(x1==x2 || x1 < 1E-5 || x2 < 1E-5){
    flag = 1;
    return 0;
  }
  siga = TNudyCore::Instance()->Interpolate(NBT, INT1, NR, fE_file3,fXsec_file3, NP, mid);
  double sigmid1 = sig1 + (sig2 - sig1)*(mid - x1)/(x2 - x1);
//  std::cout << x1 <<"  "<< x2 <<"  "<< mid <<"  "<< sig1 <<"  "<< sig2 <<"  "<< siga<<"  "<<sigmid1<<"  "<<std::fabs((siga - sigmid1)/siga)/*<<" NR "<< NR <<" NBT "<< NBT[0] <<" INT1 "<< INT1[0]*/<< std::endl;
  if(siga==0.0)return 0;
  if(std::fabs((siga - sigmid1)/siga)<=sigDiff || x1==x2){
    return 0;
  } else{
    eLinearFile3.push_back(mid); 
    xLinearFile3.push_back(siga);
  }
//  std::cout<<" linear " << mid <<"  "<<x1<< "  "<<x2 << "  "<< siga <<"  "<< sigmid1 <<"  "<< (siga - sigmid1)/siga<< std::endl;
  recursionLinearFile3(x1, mid, sig1, siga);
  recursionLinearFile3(mid, x2, siga, sig2);
  return 0;  
}


void TNudyEndfRecoPoint::ReadFile3(TNudyEndfFile *file) {
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  std::vector<double>eneExtra;
  if(LRF==0)eHi=0.0;
  int MtNumber = 0;
  double eneLow = 1E-5;
  do
    {
      eneExtra.push_back(eneLow);
      eneLow *= 1.15;
    }while(eneLow < 1E3);
  
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    std::cout <<"MT  "<<sec->GetMT() <<std::endl;
    int MT = sec->GetMT();
    MtNumbers.push_back(MT);
    MtNumber += 1;
    TIter recIter(sec->GetRecords());
    TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
    //double ZA   = sec->GetC1();
    double AWRI  = sec->GetC2();
    //double AWR  = sec->GetC2();
    //double QM   = header->GetC1();
    double QI   = header->GetC2();
    QValue[sec->GetMT()] = QI;
//      std::cout<< "file3 "<< QValue[sec->GetMT()] << std::endl; 
    //int LR = header->GetL2();
    NR = header->GetN1();
    NP = header->GetN2();
//std::cout<<sec->GetMT() <<"  "<< LR <<"  "<< NR <<"  "<< NP << std::endl;      
    TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
    NBT         = new int[NR]();
    INT1         = new int[NR]();
    fE_file3    = new double[NP]();
    fXsec_file3 = new double[NP]();

    for(int cr=0; cr < NR ; cr ++){
      NBT[cr] = tab1->GetNBT(cr);
      INT1[cr] = tab1->GetINT(cr);
//      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
//   std::cout << tab1->GetNBT(cr) <<" NBT "<< tab1->GetINT(cr)<< std::endl;
//   std::cout << tab1->GetINT(cr)<< std::endl;
    }
    for(int crs=0; crs < NP ; crs ++){
      fE_file3[crs]    = 0.0;
      fXsec_file3[crs] = 0.0;
    }
    for(int crs=0; crs < NP ; crs ++){
      eneTemp.push_back(tab1->GetX(crs));
      sigTemp.push_back(tab1->GetY(crs));
      if(eneTemp[crs]==eneTemp[crs-1])eneTemp[crs] += 0.001;//adding small difference for boundary points
      if(eneTemp[crs]==eneTemp[crs-2])eneTemp[crs] += 0.002;//adding small difference for double boundary points
//	std::cout<<"hi "<< npp <<"  "<< fE_file3[npp]<< std::endl;
    }
    Sort(eneTemp, sigTemp);
    ThinningDuplicate(eneTemp, sigTemp);
    NP = eneTemp.size();
    for(int crs=0; crs < NP ; crs ++){
      fE_file3[crs]    = eneTemp[crs];
      fXsec_file3[crs] = sigTemp[crs];
      eLinearFile3.push_back(fE_file3[crs]);
      xLinearFile3.push_back(fXsec_file3[crs]);
    }
//	std::cout << NBT[0] <<" "<< NBT[1] <<"  "<< INT1[0] <<"  "<< INT1[1] <<"  "<< eLinearFile3.size() << std::endl;
// Linearization of file 3 data
    for(int cr=0; cr < NP-1 ; cr ++){
//	  std::cout<<"cr = "<<cr<<"  "<<eLinearFile3[cr] <<"  "<< xLinearFile3[cr] << std::endl;
      recursionLinearFile3(fE_file3[cr], fE_file3[cr+1], fXsec_file3[cr], fXsec_file3[cr+1]);
    }
    for(int ie = 0; ie < eneExtra.size(); ie++){
      double sigExtra = TNudyCore::Instance()->Interpolate(NBT, INT1, NR, fE_file3,fXsec_file3, NP, eneExtra[ie]);
      if(sigExtra > 0.0){
	eLinearFile3.push_back(eneExtra[ie]);
	xLinearFile3.push_back(sigExtra);
      }
//	  std::cout<<"extra= "<<MT<<"  "<<std::setprecision(12)<<eneExtra[ie] <<"  "<<std::setprecision(12)<< sigExtra << std::endl;
    }
    Sort(eLinearFile3, xLinearFile3);
//    for(int cr=0; cr < eLinearFile3.size() ; cr ++){
//      std::cout<<"1Mt = "<<MT<<"  "<<eLinearFile3[cr] <<"  "<< xLinearFile3[cr] << std::endl;
//    }    
    ThinningDuplicate(eLinearFile3, xLinearFile3);
	
//    for(int cr=0; cr < eLinearFile3.size() ; cr ++){
//      std::cout<<"2Mt = "<<MT<<"  "<<eLinearFile3[cr] <<"  "<< xLinearFile3[cr] << std::endl;
//    }    
//Filling of array to interpolate and add with file 2 data if it is given    
    if (NBT) { delete[] NBT;    NBT = 0;  }
    if (INT1) {    delete[] INT1;    INT1 = 0;  }
    if (fE_file3) {    delete[] fE_file3;    fE_file3 = 0;  }
    if (fXsec_file3) {    delete[] fXsec_file3;    fXsec_file3 = 0;  }
    fE_file3    = new double[eLinearFile3.size()]();
    fXsec_file3 = new double[eLinearFile3.size()]();
    for(int cr=0; cr < eLinearFile3.size() ; cr ++){
      fE_file3[cr] = eLinearFile3[cr]; 
      fXsec_file3[cr] = xLinearFile3[cr];
    }
    NBT     = new int[1]();
    INT1    = new int[1]();
    NBT[0]  = eLinearFile3.size(); 
    INT1[0] = 2;
/////////////////////////////////////////////////////////////////////////    
//	std::cout<<"LSSF "<< LSSF <<" LRU "<< LRU <<" eLo1 "<<eLo1<<" eHi "<< eHi << std::endl;
//	  for(int cr=0; cr < eLinearFile3.size() ; cr ++)
//	  std::cout<<"Mt = "<<MT<<"  "<<std::setprecision(12)<<eLinearFile3[cr] <<"  "<<std::setprecision(12)<< xLinearFile3[cr] << std::endl;

// resolved range is always added in file 3
    if(LRU != 0){
      if(flagResolve != 0){
	switch(MT)
	{
	  case 2:
	  {
	  int linsize = eLinElastic.size();
	  int size1 = eLinearFile3.size();
	  for(int cr1=0; cr1 < linsize ; cr1 ++){
	    if(eLinElastic[cr1] >=eLo1 && eLinElastic[cr1] <=eHi1){
	      double crs = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE_file3,fXsec_file3,size1,eLinElastic[cr1]);
	      xLinElastic[cr1] += crs;
	    }
	  }
	  }break;
	  case 18:
	  {
	  int linsize = eLinFission.size();
	  int size1 = eLinearFile3.size();
	  for(int cr1=0; cr1 < linsize ; cr1 ++){
	    if(eLinFission[cr1] >=eLo1 && eLinFission[cr1] <=eHi1){
	      double crs = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE_file3,fXsec_file3,size1,eLinFission[cr1]);
	      xLinFission[cr1] += crs;
	    }
	  }
	  }break;
	  case 102:
	  {
	  int linsize = eLinCapture.size();
	  int size1 = eLinearFile3.size();
	  for(int cr1=0; cr1 < linsize ; cr1 ++){
	    if(eLinCapture[cr1] >=eLo1 && eLinCapture[cr1] <=eHi1){
	      double crs = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE_file3,fXsec_file3,size1,eLinCapture[cr1]);
	      xLinCapture[cr1] += crs;
	    }
	  }
	  }break;
	}
      }
//unresolved resonance region is added if LSSF = 0
      if(flagUnResolve != 0){
	switch(LSSF)
	{
	  case 0:
	  {
	    switch(MT)
	    {
	    case 2:
	    {
	      int linsize = eLinElastic.size();
	      int size1 = eLinearFile3.size();
	      for(int cr1=0; cr1 < linsize ; cr1 ++){
		if(eLinElastic[cr1] >eLo2 && eLinElastic[cr1] <=eHi2){
		  double crs = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE_file3,fXsec_file3,size1,eLinElastic[cr1]);
		  xLinElastic[cr1] += crs;
		}
	      }
	    }break;
	    case 18:
	    {
	      int linsize = eLinFission.size();
	      int size1 = eLinearFile3.size();
	      for(int cr1=0; cr1 < linsize ; cr1 ++){
		if(eLinFission[cr1] >eLo2 && eLinFission[cr1] <=eHi2){
		  double crs = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE_file3,fXsec_file3,size1,eLinFission[cr1]);
		  xLinFission[cr1] += crs;
		}
	      }
	    }break;
	    case 102:
	    {
	      int linsize = eLinCapture.size();
	      int size1 = eLinearFile3.size();
	      for(int cr1=0; cr1 < linsize ; cr1 ++){
		if(eLinCapture[cr1] >eLo2 && eLinCapture[cr1] <=eHi2){
		  double crs = TNudyCore::Instance()->Interpolate(NBT,INT1,NR,fE_file3,fXsec_file3,size1,eLinCapture[cr1]);
		  xLinCapture[cr1] += crs;
		}
	      }
	    }break;
	  }
	}break;
	case 1:
	{
	  switch(MT)
	  {
	    case 2:
	    {
	      int size1 = eLinearFile3.size();
	      for(int cr1=0; cr1 < size1 ; cr1 ++){
		if((eLinearFile3[cr1] >eLo2 && eLinearFile3[cr1] <=eHi2) || eLinearFile3[cr1] >eHi1){
		  eLinElastic.push_back(eLinearFile3[cr1]);
		  xLinElastic.push_back(xLinearFile3[cr1]);
		}
	      }
	    }break;
	    case 18:
	    {
	      int size1 = eLinearFile3.size();
	      for(int cr1=0; cr1 < size1 ; cr1 ++){
		if((eLinearFile3[cr1] >eLo2 && eLinearFile3[cr1] <=eHi2) || eLinearFile3[cr1] >eHi1){
		  eLinFission.push_back(eLinearFile3[cr1]);
		  xLinFission.push_back(xLinearFile3[cr1]);
		}
	      }
	    }break;
	    case 102:
	    {
	      int size1 = eLinearFile3.size();
	      for(int cr1=0; cr1 < size1 ; cr1 ++){
		if((eLinearFile3[cr1] >eLo2 && eLinearFile3[cr1] <=eHi2) || eLinearFile3[cr1] >eHi1){
		  eLinCapture.push_back(eLinearFile3[cr1]);
		  xLinCapture.push_back(xLinearFile3[cr1]);
		}
	      }
	    }break;
	  }	  
	}break;
      }
    }
//This adds additional points to the elastic, fission and capture below resolved range and above uresolved range
     switch(MT)
     {
       case 2:
	{
	  if(flagResolve != 0){
	    for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	      if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] <=eLo1){
		eLinElastic.push_back(eLinearFile3[cr]);
		xLinElastic.push_back(xLinearFile3[cr]);
	      }
	    }
	  }else if(flagResolve == 0 && flagUnResolve != 0){
	    for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	      if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] <=eLo2){
		eLinElastic.push_back(eLinearFile3[cr]);
		xLinElastic.push_back(xLinearFile3[cr]);
	      }
	    }
	  }
	  if(flagUnResolve != 0){
	    for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	      if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] >eHi2){
		eLinElastic.push_back(eLinearFile3[cr]);
		xLinElastic.push_back(xLinearFile3[cr]);
	      }
	    }
	  }else if(flagUnResolve == 0){
	    for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	      if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] >eHi1){
		eLinElastic.push_back(eLinearFile3[cr]);
		xLinElastic.push_back(xLinearFile3[cr]);
	      }
	    }
	  }
	}break;
       case 18:
	{
	  if(flagResolve != 0){
	    for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	      if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] <=eLo1){
		eLinFission.push_back(eLinearFile3[cr]);
		xLinFission.push_back(xLinearFile3[cr]);
	      }
	    }
	  }else if(flagResolve == 0 && flagUnResolve != 0){
	    for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	      if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] <=eLo2){
		eLinFission.push_back(eLinearFile3[cr]);
		xLinFission.push_back(xLinearFile3[cr]);
	      }
	    }
	  }
	if(flagUnResolve != 0){
	  for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	    if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] >eHi2){
	      eLinFission.push_back(eLinearFile3[cr]);
	      xLinFission.push_back(xLinearFile3[cr]);
	    }
	  }
	}else if(flagUnResolve == 0){
	  for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	    if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] >eHi1){
	      eLinFission.push_back(eLinearFile3[cr]);
	      xLinFission.push_back(xLinearFile3[cr]);
	    }
	  }
	}
      }break;
      case 102:
      {
	if(flagResolve != 0){
	  for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	    if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] <=eLo1){
	      eLinCapture.push_back(eLinearFile3[cr]);
	      xLinCapture.push_back(xLinearFile3[cr]);
	    }
	   }
	}else if(flagResolve == 0 && flagUnResolve != 0){
	  for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	    if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] <=eLo2){
	      eLinCapture.push_back(eLinearFile3[cr]);
	      xLinCapture.push_back(xLinearFile3[cr]);
	    }
	  }
	}
	if(flagUnResolve != 0){
	  for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	    if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] >eHi2){
	      eLinCapture.push_back(eLinearFile3[cr]);
	      xLinCapture.push_back(xLinearFile3[cr]);
	    }
	  }
	}else if(flagUnResolve == 0){
	  for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	    if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] >eHi1){
	      eLinCapture.push_back(eLinearFile3[cr]);
	      xLinCapture.push_back(xLinearFile3[cr]);
	    }
	  }
	}
      }break;
    }
  }else{// no resonance parameters are given
    switch(MT)
    {
      case 2:
	{
	  for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	    eLinElastic.push_back(eLinearFile3[cr]);
	    xLinElastic.push_back(xLinearFile3[cr]);
	  }
	}break;
      case 18:
	{
	  for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	    eLinFission.push_back(eLinearFile3[cr]);
	    xLinFission.push_back(xLinearFile3[cr]);
	  }
	}break;
      case 102:
	{
	  for(int cr=0; cr < eLinearFile3.size() ; cr ++){
	    eLinCapture.push_back(eLinearFile3[cr]);
	    xLinCapture.push_back(xLinearFile3[cr]);
	  }
	}break;
      }
    }    
  
// check the sum of cross sections with total cross-section given in the MT=1 + resonance
    if(MT==1){
      for(int cr=0; cr<eLinearFile3.size();cr++){
	energy.push_back(eLinearFile3[cr]);
	sigmaT.push_back(xLinearFile3[cr]);
      }
    }

//      int size = eLinearFile3.size();
      
//      for(int i = 0; i < eLinearFile3.size(); i++){
//	if(xLinearFile3[i] < 1E-20){
//	  xLinearFile3.erase(xLinearFile3.begin()+i);
//	  eLinearFile3.erase(eLinearFile3.begin()+i);
//	}
//      }

      if(MT == 2) {
//	energyMts.insert(std::end(energyMts), std::begin(eLinElastic), std::end(eLinElastic));
	sigmaMts.insert(std::end(sigmaMts), std::begin(eLinElastic), std::end(eLinElastic));
	sigmaMts.insert(std::end(sigmaMts), std::begin(xLinElastic), std::end(xLinElastic));
	energyUni.insert(std::end(energyUni), std::begin(eLinElastic), std::end(eLinElastic));
      } else if(MT == 18) {
//	energyMts.insert(std::end(energyMts), std::begin(eLinFission), std::end(eLinFission));
	sigmaMts.insert(std::end(sigmaMts), std::begin(eLinFission), std::end(eLinFission));
	sigmaMts.insert(std::end(sigmaMts), std::begin(xLinFission), std::end(xLinFission));
	energyUni.insert(std::end(energyUni), std::begin(eLinFission), std::end(eLinFission));
      } else if(MT == 102) {
//	energyMts.insert(std::end(energyMts), std::begin(eLinCapture), std::end(eLinCapture));
	sigmaMts.insert(std::end(sigmaMts), std::begin(eLinCapture), std::end(eLinCapture));
	sigmaMts.insert(std::end(sigmaMts), std::begin(xLinCapture), std::end(xLinCapture));
	energyUni.insert(std::end(energyUni), std::begin(eLinCapture), std::end(eLinCapture));
      } else if (MT != 2 && MT != 18 && MT != 102) {
//	energyMts.insert(std::end(energyMts), std::begin(eLinearFile3), std::end(eLinearFile3));
	sigmaMts.insert(std::end(sigmaMts), std::begin(eLinearFile3), std::end(eLinearFile3));
	sigmaMts.insert(std::end(sigmaMts), std::begin(xLinearFile3), std::end(xLinearFile3));
	energyUni.insert(std::end(energyUni), std::begin(eLinearFile3), std::end(eLinearFile3));
      }
//    if ( MT==2 )energyUni.insert(it,eLinElastic.begin(),eLinElastic.end());
//    if ( MT ==18 )energyUni.insert(it,eLinFission.begin(),eLinFission.end());
//    if ( MT == 102 )energyUni.insert(it,eLinCapture.begin(),eLinCapture.end());
    if(MT != 1 && MT != 2 && MT != 3 && MT != 4 && MT != 27 && MT != 18 && MT != 19 && MT != 20 && MT != 21 && MT != 38 && MT != 101 && MT != 102 && MT < 250){
//  fixupTotal(eLinearFile3, xLinearFile3);
    }
    eLinearFile3.clear();	
    xLinearFile3.clear();
    eneTemp.clear();	
    sigTemp.clear();	

    if (NBT) { delete[] NBT;    NBT = 0;  }
    if (INT1) {    delete[] INT1;    INT1 = 0;  }
    if (fE_file3) {    delete[] fE_file3;    fE_file3 = 0;  }
    if (fXsec_file3) {    delete[] fXsec_file3;    fXsec_file3 = 0;  }
    
//    for(int i = 0; i < sigmaMts.size(); i++){
//      std::cout<<sigmaMts[i]<<std::endl;
//    }
    
    sigmaOfMts.push_back(sigmaMts);
    sigmaMts.clear();
//    std::vector<double>().swap(sigmaMts);
 
  }
  MtValues.push_back(MtNumbers);
  std::vector<int>().swap(MtNumbers);
//  sigmaOfMts.push_back(energyMts);
//  sigmaOfMts.push_back(sigmaMts);
//  std::vector<int>().swap(energyMts);
//  std::vector<double>().swap(sigmaMts);
  NoOfElements += 1;
}

void TNudyEndfRecoPoint::cdfGenerateE(std::vector<double> &x1,std::vector<double> &x2){
  double energy5Cdf = 0.0;
  for(int cr=0; cr < x1.size() ; cr ++){
    energy5Cdf += x2[cr];
  }
  double df = 0.0;
  for(int cr=0; cr < x1.size() ; cr ++){
    if(energy5Cdf > 0.0)x2[cr] = x2[cr]/energy5Cdf;
    df += x2[cr];
    energyCdfFile5.push_back(df);
//        std::cout <<"cos = "<< energyFile5[cr] <<"  "<< energyPdfFile5[cr] <<"  "<< energyCdfFile5[cr]  << std::endl;
  }
}

double TNudyEndfRecoPoint::recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2){
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

double TNudyEndfRecoPoint::recursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy){
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

double TNudyEndfRecoPoint::recursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy){
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

double TNudyEndfRecoPoint::recursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy){
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
//______________________________________________________________________________
void TNudyEndfRecoPoint::ReadFile5(TNudyEndfFile *file) {
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
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
	  for(int i = 0; i < energyFile5.size(); i++){
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
	  for(int i = 0; i < energyFile5.size(); i++){
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
	  for(int i = 0; i < energyFile5.size(); i++){
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
	  for(int i = 0; i < energyFile5.size(); i++){
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
	  for(int i = 0; i < energyFile5.size(); i++){
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
	  for(int i = 0; i < energyFile5.size(); i++){
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
  /*  
  for(int i = 0; i < energy5OfMts.size(); i++){
      std::cout <<" mt "<<Mt5Values[0][i]<<" size "<< energy5OfMts[i].size() << std::endl;
    for(int j = 0; j < energy5OfMts[i].size(); j++){
      std::cout << energy5OfMts[i][j] << std::endl;
      for(int k =0; k < energyPdf5OfMts[i][j].size()/2; k++){
	std::cout << energyPdf5OfMts[i][j][2*k] <<"  "<< energyPdf5OfMts[i][j][2*k + 1] <<"  "<< energyCdf5OfMts[i][j][2*k + 1]<< std::endl;
      }
    }
  }
  */
}

//______________________________________________________________________________
void TNudyEndfRecoPoint::ReadFile6(TNudyEndfFile *file) {
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
      int LCT = sec->GetL2();
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
	    double Mn1 = 1, Mp = 1, Md = 1,Mt = 1, Mhe3 = 1, Malpha = 0;
	    double mn = 0.5, mp = 1, md = 1, mt = 1, mhe3 = 1, malpha = 2;

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
	    delete[] lCoef;
	  }else if (LANG >0){
	    for (int i = 0; i < NL; i++) {
	      cosFile4.push_back(header->GetLIST(2*i+0));
	      cosPdfFile4.push_back(header->GetLIST(2*i+1));
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
	double AWPR =tab11->GetC2(); 
	int LIPR =tab11->GetL1(); 
	int LAW = tab11->GetL2();
	NR = tab11->GetN1();
	NP = tab11->GetN2();
	divr =div(ZAPR,1000);
	double ResidueA = divr.rem;
	double ResidueZ = divr.quot;
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
	    for(int i = 0; i < energyFile5.size(); i++){
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
	  for (int i = 0; i< energyFile5.size(); i++){
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
  for(int i = 0; i < energy5OfMts.size(); i++){
      std::cout <<" mt "<<Mt5Values[0][i]<<" size "<< energy5OfMts[i].size() << std::endl;
    for(int j = 0; j < energy5OfMts[i].size(); j++){
      std::cout << energy5OfMts[i][j] << std::endl;
      for(int k =0; k < energyPdf5OfMts[i][j].size()/2; k++){
	std::cout << energyPdf5OfMts[i][j][2*k] <<"  "<< energyPdf5OfMts[i][j][2*k + 1] <<"  "<< energyCdf5OfMts[i][j][2*k + 1]<< std::endl;
      //for(int k =0; k < cosPdf4OfMts[i][j].size()/2; k++){
	//std::cout << cosPdf4OfMts[i][j][2*k] <<"  "<< cosPdf4OfMts[i][j][2*k + 1] <<"  "<< cosCdf4OfMts[i][j][2*k + 1]<< std::endl;
      }
    }
  }
  */
}
  double TNudyEndfRecoPoint::cmToLabElasticE(double inE, double cmCos, double awr){
    return inE * ((1 + awr * awr + 2 * awr * cmCos)/((1 + awr) * (1 + awr)));
  }
  double TNudyEndfRecoPoint::cmToLabElasticCosT(double cmCos, double awr){
    double sint = sqrt(1 - cmCos * cmCos);
    return atan(awr * sint/(awr * cmCos + 1));
  }
  double TNudyEndfRecoPoint::cmToLabInelasticE(double cmEOut, double inE, double cmCos, double awr){
    double mass = awr + 1;
    return cmEOut + (inE + 2 * cmCos * mass * sqrt(inE * cmEOut))/(mass*mass);
  }
  double TNudyEndfRecoPoint::cmToLabInelasticCosT(double labEOut, double cmEOut, double inE, double cmCos, double awr){
    double mass = awr + 1;
    return cmCos * sqrt(cmEOut / labEOut) + sqrt (inE / labEOut) / mass;
  }

void TNudyEndfRecoPoint::ReadFile8(TNudyEndfFile *file) { //fission product yield data
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
//    std::cout <<"file  "<<sec->GetMT() <<std::endl;
    TIter recIter(sec->GetRecords());
    //double ZA   = sec->GetC1();
    //double AWR  = sec->GetC2();
    int MT = sec->GetMT();
    int LE = sec->GetL1();
    if(MT == 454){// Neutron induced independent fission yield
      TArrayD ein(LE);
      zafp = new TArrayD[LE];
      fps = new TArrayD[LE];
      yi = new TArrayD[LE];
      dyi = new TArrayD[LE];
      for (int i =0; i < LE; i++){
        TNudyEndfList *list1 = (TNudyEndfList *)recIter.Next();
	ein[i]   = list1->GetC1();
	//int NN   = list1->GetN1();
	int NFP  = list1->GetN2();
//      std::cout<<"energy i " <<ein[i] << std::endl;
	zafp[i].Set(NFP);
	fps[i].Set(NFP);
	yi[i].Set(NFP);
	dyi[i].Set(NFP);
//	    std::cout<<"energy "<< ein[lis] << std::endl; 
          for (int j = 0; j < NFP; j++){
            zafp[i].SetAt(list1->GetLIST(4*j+0), j);
            fps[i].SetAt(list1->GetLIST(4*j+1), j);
            yi[i].SetAt(list1->GetLIST(4*j+2), j);
            dyi[i].SetAt(list1->GetLIST(4*j+3), j);
//	    std::cout<<"fission yield i "<< zafp[i].At(j) << std::endl; 
	  }
	}
      }else if(MT == 459){
          TArrayD einc(LE);
	  zafpc = new TArrayD[LE];
          fpsc = new TArrayD[LE];
	  yc = new TArrayD[LE];
	  dyc = new TArrayD[LE];
	  for (int i =0; i < LE; i++){
	    TNudyEndfList *list1 = (TNudyEndfList *)recIter.Next();
	    einc[i]   = list1->GetC1();
	    //int NN   = list1->GetN1();
	    int NFP  = list1->GetN2();
//      std::cout<<"energy c " <<einc[i] << std::endl;
	    zafpc[i].Set(NFP);
	    fpsc[i].Set(NFP);
	    yc[i].Set(NFP);
	    dyc[i].Set(NFP);
//	    std::cout<<"energy "<< ein[lis] << std::endl; 
	    for (int j = 0; j < NFP; j++){
	      zafpc[i].SetAt(list1->GetLIST(4*j+0), j);
	      fpsc[i].SetAt(list1->GetLIST(4*j+1), j);
              yc[i].SetAt(list1->GetLIST(4*j+2), j);
	      dyc[i].SetAt(list1->GetLIST(4*j+3), j);
//	    std::cout<<"fission yield c "<< zafpc[i].At(j) << std::endl; 
	    }
	  }
      }
   }
}

double TNudyEndfRecoPoint::K_wnum(double x) {
  double k; 
  k = factor_k * sqrt(std::fabs(x));
  return k;
}

double TNudyEndfRecoPoint::GetRho(double x1, int lVal) {
  if (!NAPS && !NRO) return factor_k * sqrt(std::fabs(x1)) * rad_a;   	// Use rad_a in the penetrabilities Pl and shift factors Sl , and AP in the hard-sphere phase shifts φl 
  if (!NAPS && NRO) return factor_k * sqrt(std::fabs(x1)) *  rad_a;    	// Use rad_a in the penetrabilities Pl and shift factors Sl,AP(E) in the hard-sphere phase shifts φl
  if ( NAPS && !NRO) return factor_k * sqrt(std::fabs(x1)) * APL[lVal];      	// use AP in the penetrabilities and shift factor as well as in the phase shifts
  if ( NAPS && NRO) return factor_k * sqrt(std::fabs(x1)) * APL[lVal];       	// read AP(E) and use it in all three places, Pl , Sl , φl
  return 0;  
}
 
double TNudyEndfRecoPoint::GetRhoC(double x, int isDiff, int lVal) {
  if (!isDiff) return factor_k * sqrt(std::fabs(x)) * APL[lVal];    		//read AP(E) and use it in all three places, Pl , Sl , φl
  if (isDiff) return factor_k * sqrt(std::fabs(x)) * APL[lVal];                	//read AP(E) and use it in all three places, Pl , Sl , φl
  return 0;  
} 

double TNudyEndfRecoPoint::calcPhi(double x1, int l) {
  double  x = GetRhoC(x1, 0, l);
  switch (l) 
  {
    case 0:
    return x;
    case 1:
    return (x - atan(x));
    case 2:
    return (x - atan(3.0 * x/(3.0 - x2(x))));
    case 3:
    return (x - atan(x * (15.0 - x*x)/(15.0 - 6.0 * x2(x))));
    case 4:
    return (x - atan(x * (105.0 - 10.0 * x2(x))/(105.0 - 45.0 * x2(x) + x4(x))));
  }
  return 0;
}

double TNudyEndfRecoPoint::calcShift(double x1, int l) {
  double  x = GetRho(x1, l);
  switch (l) 
  {
    case 0:
    return 0.0;
    case 1:
    return (-1.0 / (1.0 + x*x));
    case 2:
    return (-(18.0 + 3.0 * x2(x)) / Fac2(x));
    case 3:
    return (-(675.0 + 90.0 * x2(x) + 6.0 * x4(x))/(x6(x) + 6.0 * x4(x) + 45.0 * x2(x) + 225.0));
    case 4:
    return (-(44100.0 + 4725.0 * x2(x) + 270.0 * x4(x) + 10.0 * x6(x))/(x8(x) + 10.0 * x6(x) + 135.0 * x4(x) + 1575.0 * x2(x) + 11025.0));
  }
  return 0;
}

double TNudyEndfRecoPoint::calcPene(double x1, int l) {
  double x = GetRho(x1, l);
  switch (l) 
  {
    case 0:
    return x;
    case 1:
    return (x*x*x)/(1.0 + x*x);
    case 2:
    return (x5(x)/Fac2(x));
    case 3:
    return ((x * x6(x))/(225.0 + 45.0 * x2(x) + 6.0 * x4(x) + x6(x)));
    case 4:
    return ((x * x8(x))/(x8(x) + 10.0 * x6(x) + 135.0 * x4(x) + 1575.0 * x2(x) + 11025.0));
  }
  return 0;
}


double TNudyEndfRecoPoint::GetERP(double x, int r, int lVal) {
  double er = Er[r];
  if(lVal==0)return er;
  return er + (ShiftEr[r] - calcShift(x, lVal))/(2.0 * PhiEr[r]) * Gamma_nrE(std::fabs(er), r, lVal);
}

double TNudyEndfRecoPoint::Gamma_nrE(double x, int ii, int lval) {
  return  calcPene(x, lval) * Gamma_n[ii];
}

double TNudyEndfRecoPoint::Gamma_xrE(int ii, int lrx) {
  if (!lrx) {
    return 0;
  } else {
    return Gamma_r[ii] - (Gamma_n[ii] + Gamma_g[ii] + Gamma_f[ii]);
  }
}

double TNudyEndfRecoPoint::Gamma_rE(double x, int ii, int lval, int lrx) {
  return Gamma_nrE(x, ii, lval) + Gamma_xrE(ii, lrx) + Gamma_g[ii] + Gamma_f[ii];
}
int TNudyEndfRecoPoint::widthFluctuation(double GNR, double gx, double gg, double gf,int jval ){
  double XX[10][5]={
  {3.0013465E-03,1.3219203E-02,1.0004488E-03,1.3219203E-02,1.0E+0},
  {7.8592886E-02,7.2349624E-02,2.6197629E-02,7.2349624E-02,0.0E+0},
  {4.3282415E-01,1.9089473E-01,1.4427472E-01,1.9089473E-01,0.0E+0},
  {1.3345267E+00,3.9528842E-01,4.4484223E-01,3.9528842E-01,0.0E+0},
  {3.0481846E+00,7.4083443E-01,1.0160615E+00,7.4083443E-01,0.0E+0},
  {5.8263198E+00,1.3498293E+00,1.9421066E+00,1.3498293E+00,0.0E+0},
  {9.9452656E+00,2.5297983E+00,3.3150885E+00,2.5297983E+00,0.0E+0},
  {1.5782128E+01,5.2384894E+00,5.2607092E+00,5.2384894E+00,0.0E+0},
  {2.3996824E+01,1.3821772E+01,7.9989414E+00,1.3821772E+01,0.0E+0},
  {3.6216208E+01,7.5647525E+01,1.2072069E+01,7.5647525E+01,0.0E+0}
  };
  
  double WW[10][5]={
  {1.1120413E-01,3.3773418E-02,3.3376214E-04,1.7623788E-03,1.0E+0},
  {2.3546798E-01,7.9932171E-02,1.8506108E-02,2.1517749E-02,0.0E+0},
  {2.8440987E-01,1.2835937E-01,1.2309946E-01,8.0979849E-02,0.0E+0},
  {2.2419127E-01,1.7652616E-01,2.9918923E-01,1.8797998E-01,0.0E+0},
  {1.0967668E-01,2.1347043E-01,3.3431475E-01,3.0156335E-01,0.0E+0},
  {3.0493789E-02,2.1154965E-01,1.7766657E-01,2.9616091E-01,0.0E+0},
  {4.2930874E-03,1.3365186E-01,4.2695894E-02,1.0775649E-01,0.0E+0},
  {2.5827047E-04,2.2630659E-02,4.0760575E-03,2.5171914E-03,0.0E+0},
  {4.9031965E-06,1.6313638E-05,1.1766115E-04,8.9630388E-10,0.0E+0},
  {1.4079206E-08,0.0000000E+00,5.0989546E-07,0.0000000E+00,0.0E+0}
  };
  
  RN = 0.0;
  RG = 0.0;
  RF = 0.0; 
  int mux = (int)amux[jval];
  int mun = (int)amun[jval];
  int muf = (int)amuf[jval];
//     GNR = GNO[i];
//     std::cout<<GNR<<std::endl;
  if(GNR <= 0.0 || gg <= 0.0 ){return 0;}
  if(mun < 1.0 || mun >4.0)mun=5.0;
  if(mun < 1.0 || mun >4.0)mun=5.0;
  if(mux < 1.0 || mux >4.0)mux=5.0;
  if(gf < 0.0 ){return 0;}
  if(gf == 0.0 && gx == 0.0){
    for(int j=0; j<10; j++){
      double XJ = XX[j][mun-1];
      double fak = XJ * WW[j][mun-1]/(GNR *XJ + gg);
      RN += fak *XJ;
      RG += fak;
    }//for loop
  }//if
  if(gf == 0.0 && gx > 0.0){
    for(int j=0; j<10; j++){
      double XJ = XX[j][mun-1];
      double WJXJ = WW[j][mun-1]*XJ;
      double EFFJ = GNR*XJ + gg;
      for(int k=0; k<10; k++){
	double XK = XX[k][mux-1];
	double fk = WW[k][mux-1]*WJXJ/(EFFJ + gx*XK);
	RN += XJ*fk;
	RG += fk;
	RX += XK*fk;
      }//2 for loop
    }//1 for loop
  }//if
  if(gf > 0.0 && gx == 0.0){
    for(int j=0; j<10; j++){
      double XJ = XX[j][mun-1];
      double WJXJ = WW[j][mun-1]*XJ;
      double EFFJ = GNR*XJ + gg;
      for(int k=0; k<10; k++){
	double fk = WW[k][muf-1]*WJXJ/(EFFJ + gf*XX[k][muf-1]);
	RN += XJ*fk;
	RG += fk;
	RF += XX[k][muf-1]*fk;
      }//2 for loop
    }//1 for loop
  }//if
  if(gf > 0.0 && gx > 0.0){
//      std::cout<<GNR<<" GNR "<<gg<<" gg "<<gf<<" gf "<<gx<<std::endl;
    
    for(int j=0; j<10; j++){
      double XJ = XX[j][mun-1];
      double WJXJ = WW[j][mun-1]*XJ;
      double EFFJ = GNR*XJ + gg;
      for(int k=0; k<10; k++){
	double XK = XX[k][muf-1];
	double WKWJXJ = WW[k][muf-1]*WJXJ;
	double EFFJK = gf*XK + EFFJ;
	for(int l=0; l<10; l++){
	  double XL = XX[l][mux-1];
	  double fk = WW[l][mux-1]*WKWJXJ/(EFFJK + gx*XL);
	  RN += XJ*fk;
	  RG += fk;
	  RX += fk*XL;
	  RF += XK*fk;
	}//3 for loop
      }//2 for loop
    }//1 for loop
  }//if
  return 0;
}
void TNudyEndfRecoPoint::InverseMatrix() {
  double A[3][3], C[3][3], C2[3][3], C3[3][3], AI[3][3], AIB[3][3];
  for(int i =0; i <3; i++){
    for(int j =0; j <3; j++){
      A[i][j] = R[i][j];
    }
    A[i][i] = A[i][i] + 1.0;
  }
  double det1 = A[1][1] * A[2][2] - A[1][2]*A[2][1];
  double det2 = A[1][2] * A[2][0] - A[1][0]*A[2][2];
  double det3 = A[1][0] * A[2][1] - A[1][1]*A[2][0];
  double det = A[0][0] * det1 + A[0][1]*det2 + A[0][2]*det3;
  AI[0][0] = det1/det;
  AI[1][0] = det2/det;
  AI[2][0] = det3/det;
  AI[0][1] = AI[1][0];
  AI[1][1] = (A[0][0]*A[2][2] - A[2][0] *A[0][2])/det; 
  AI[2][1] = (A[0][1]*A[2][0] - A[0][0] *A[2][1])/det; 
  AI[0][2] = AI[2][0];
  AI[1][2] = AI[2][1]; 
  AI[2][2] = (A[0][0]*A[1][1] - A[0][1] *A[1][0])/det;
  
  for(int i =0; i <3; i++){
    for(int j =0; j <3; j++){
      double sum1 = 0.0;
      for(int k =0; k <3; k++){
	sum1 += AI[i][k]*S[k][j];
      }
      AIB[i][j] = sum1;
    }
  }
  for(int i =0; i <3; i++){
    for(int j =0; j <3; j++){
      double sum1 = 0.0;
      for(int k =0; k <3; k++){
	sum1 += S[i][k]*AIB[k][j];
      }
      C[i][j] = A[i][j] + sum1;
      C2[i][j] = R[i][j] + sum1;
    }
  }
  det1 = C[1][1] * C[2][2] - C[1][2]*C[2][1];
  det2 = C[1][2] * C[2][0] - C[1][0]*C[2][2];
  det3 = C[1][0] * C[2][1] - C[1][1]*C[2][0];
  det = C[0][0] * det1 + C[0][1]*det2 + C[0][2]*det3;
  C3[0][0] = det1/det;
  C3[1][0] = det2/det;
  C3[2][0] = det3/det;
  C3[0][1] = C3[1][0];
  C3[1][1] = (C[0][0]*C[2][2] - C[2][0] *C[0][2])/det; 
  C3[2][1] = (C[0][1]*C[2][0] - C[0][0] *C[2][1])/det; 
  C3[0][2] = C3[2][0];
  C3[1][2] = C3[2][1]; 
  C3[2][2] = (C[0][0]*C[1][1] - C[0][1] *C[1][0])/det;
  for(int i =0; i <3; i++){
    for(int j =0; j <3; j++){
      double sum1 = 0.0;
      for(int k =0; k <3; k++){
	sum1 += C3[i][k]*C2[k][j];
      }
      C[i][j] = - sum1;
      RI[i][j] = - sum1;
    }
    C[i][i] = C[i][i] + 1.0;
  }
  for(int i =0; i <3; i++){
    for(int j =0; j <3; j++){
      double sum1 = 0.0;
      for(int k =0; k <3; k++){
	sum1 += AIB[i][k]*C[k][j];
      }
      SI[i][j] = - sum1;
    }
  }
}
void TNudyEndfRecoPoint::GetSigmaRMP(double x, double &siga, double &sigb, double &sigc) {
  std::vector<double> row;
  int nrsl,nrs=0,jVal,nrstemp;
  double fac2;
  double phil, sfil, cfil, sfil2,s2fil;
  double deno=0.0;
  double ERP, GNR;
  double sumCrsElastic = 0.0, sumCrsCapture = 0.0, sumCrsFission = 0.0 ; 
  double GR;
  double reduced_n, reduced_g, reduced_fa, reduced_fb;
  x = std::fabs(x);
  double pibyk2 = PI/x2(K_wnum(x));
  double rt11=0.0, st11=0.0;
  double sigrm1=0.0, sigab1=0.0, sigf1=0.0, sigel=0.0;
  double gj;
  
  for (int lcnt = 0; lcnt < NLS; lcnt++) {
    fac2 = 0.0;
    int lval = l[lcnt];
//    std::cout<<"L= "<<lval<<"  "<<APL[lval]<<std::endl;
    nrsl = NRS[lval];
    phil = calcPhi(x, l[lcnt]);
    sfil = sin(phil);
    cfil = cos(phil);
    sfil2 = sfil*sfil;
    s2fil = 2. * sfil *cfil ;
    fac2 = 4.0 * MisGj[l[lcnt]] * sfil2* pibyk2;
    
//    std::cout<<" "<< GetRho(x, 0, lval) <<"  "<< GetRhoC(x, 0, lval) <<"  "<< rad_a <<"  "<< AP << std::endl;
    sigrm1=0.0; sigel=0.0, sigab1=0.0, sigf1=0.0;
	     
//    for (int jVal = 0; jVal < Jnum; jVal++) {
//    for (int jVal = 0; jVal < NJValue[l[lcnt]]; jVal++) {
    jVal = 0; nrstemp = nrs; int nrend = nrs + nrsl;
    do
    {
      rt11 = 0.0; st11 = 0.0;
      for(int i =0; i <3; i++){
        for(int j =0; j <3; j++){
	  R[i][j]=0.0;
	  S[i][j]=0.0;
	  RI[i][j]=0.0;
	  SI[i][j]=0.0;
        }
      }
//        std::cout<<" L= "<<l[lcnt]<<" E "<< x<<" jVal "<<jVal<<"  "<< MissingJ[l[lcnt]][jVal] <<std::endl;
      for (int r = nrs; r < nrend; r++) {
//       std::cout<<" jVal "<<jVal<<"  "<<  J[r] <<"  "<< MissingJ[l[lcnt]][jVal] <<"  "<< nrs << "  "<<NJValue[l[lcnt]]<<std::endl;
        if(MissingJ[l[lcnt]][jVal] != (J[r])){
	  nrs=r;
	  break;
	}
//        std::cout<<r<< "  "<< J[r]<<"   "<< nrs<< "  "<<Gamma_g[r]<< "  "<< nrsl<<  std::endl;
//	{
//	  continue;}	
	// Gamma_n,r term
	reduced_n = 0.5 * Gamma_nrE(x, r, l[lcnt]);//calcPene(x, l[lcnt]) * std::fabs(Gamma_n[r])/calcPene(std::fabs(Er[r]), l[lcnt]);  // 
	reduced_g = 0.5 * std::fabs(Gamma_g[r]);  // 
//	GR = reduced_n + reduced_g + reduced_fa + reduced_fb;
	ERP = Er[r];
//	ERP = GetERP(x, r,  lVal, 0);
	GR = reduced_g;
	GNR = reduced_n;
	double DE     = ERP - x;
	deno   = DE * DE + GR * GR;
        double de2    = DE/deno;
	double gamma2 = GR/deno;
	if(LFW==0){
	  rt11 += gamma2*GNR;
	  st11 += de2*GNR;
        }else{
	  reduced_fa = Gamma_fasq[r];  // 
	  reduced_fb = Gamma_fbsq[r];  //
	  double GR2 = sqrt(GNR);
	  double GFA = reduced_fa * reduced_fa;
	  double GFB = reduced_fb * reduced_fb;
	  double GFAN = reduced_fa * GR2;
	  double GFBN = reduced_fb * GR2;
	  double GFAB = reduced_fa * reduced_fb;
	  R[0][0] += gamma2*GNR;
	  R[0][1] += gamma2*GFAN;
	  R[0][2] += gamma2*GFBN;
	  R[1][1] += gamma2*GFA;
	  R[1][2] += gamma2*GFAB;
	  R[2][2] += gamma2*GFB;
	  S[0][0] += de2*GNR;
	  S[0][1] += de2*GFAN;
	  S[0][2] += de2*GFBN;
	  S[1][1] += de2*GFA;
	  S[1][2] += de2*GFAB;
	  S[2][2] += de2*GFB;
	  R[1][0] = R[0][1];
	  S[1][0] = S[0][1];
	  R[2][0] = R[0][2];
	  S[2][0] = S[0][2];
	  R[2][1] = R[1][2];
	  S[2][1] = S[1][2];
	}
      }
      if(LFW==0){ 
	gj = (2.*std::fabs(MissingJ[l[lcnt]][jVal]) + 1.)/(4.*SPI + 2.0);
        double det = (rt11 + 1.0) * (rt11 + 1.0) + st11 * st11;
//	    std::cout<<st11<<"  "<<rt11<<std::endl;
	double SI11 = - st11/det;
        double RI11 = -(rt11*(rt11 + 1.0) + st11*st11)/det;
	sigrm1 += -4. *gj*(RI11 +(RI11 * RI11 + SI11 * SI11));   
	sigel += gj * ((2.0*sfil2 + 2.*RI11)*(2.0*sfil2 + 2.*RI11)+(s2fil + 2.0*SI11)*(s2fil + 2.0*SI11));
	sigf1 = 0.0;
      }else{
	InverseMatrix();
	gj = (2.*std::fabs(MissingJ[l[lcnt]][jVal]) + 1.)/(4.*SPI + 2.0);
	double SI11 = SI[0][0];
	double RI11 = RI[0][0];
	double SI12 = SI[0][1];
	double RI12 = RI[0][1];
	double SI13 = SI[0][2];
	double RI13 = RI[0][2];
	sigab1 += -4.0*gj * (RI11 + (RI11*RI11 + SI11*SI11));
	sigel += gj * ((2.0*sfil2 + 2.*RI11)*(2.0*sfil2 + 2.*RI11)+(s2fil + 2.0*SI11)*(s2fil + 2.0*SI11));
	sigf1 += 4.0*gj*(RI12*RI12 + RI13*RI13 + SI12*SI12 + SI13*SI13);
	sigrm1 = sigab1 - sigf1;
      }
      jVal += 1;
    }while(jVal < NJValue[l[lcnt]]);
    nrs = nrstemp;
    nrs += nrsl;
    sumCrsElastic +=  pibyk2 * (sigel) + fac2;
    sumCrsCapture += sigrm1 * pibyk2;
    sumCrsFission += sigf1 * pibyk2;
///	   std::cout<<" l= "<< lVal<< "  "<< pibyk2 * (sigel)<<"   " <<fac1 <<"  "<< 2.0*sfil2 <<"  "<< fac3 << "  "<<pibyk2<< std::endl;
//   std::cout<<phil<<"  "<<philp<<"  "<<std::setprecision(8)<< x << "  "<<std::setprecision(10)<<sumCrsElastic << "  "<<std::setprecision(10)<<sumCrsCapture<< "  "<<std::setprecision(10)<<sumCrsFission << std::endl;
  }
//         energy.push_back(x);
//	  sigma.push_back(sumCrsElastic);
  siga = sumCrsElastic;	    
  sigb = sumCrsCapture;	    
  sigc = sumCrsFission;	    
}
void TNudyEndfRecoPoint::GetSigma(int lrfa, double x, double &siga, double &sigb, double &sigc) {
  std::vector<double> row;
  int nrsl,nrs=0,lVal=0, jVal,nrstemp;
  double fac1=0.0, fac2 = 0.0;
  double phil, sfil, cfil, sfil2, s2fil;
  double deno=0.0;
  double ERP, ERPP1, GNR;
  double fachNumElastic=0.0 , fachNumCapture=0.0 ,  fachNumFission=0.0;
  double facElastic = 0.0 , facCapture = 0.0 ,  facFission = 0.0;
  double sumElastic = 0.0, sumCapture = 0.0 , sumFission = 0.0;
  double sumCrsElastic = 0.0, sumCrsCapture = 0.0, sumCrsFission = 0.0 ; 
  double GNRS, GR, GS;
  double fachNumElasticG ,fachNumElasticH , fachNumElasticM=0.0;
  double facElasticM = 0.0 , sumElasticM = 0.0 ;
  x = std::fabs(x);
  double pibyk2 = PI/x2(K_wnum(x));
  switch(lrfa)
  {
// SLBW--------
    case 1:
      switch(LRU)
      {
	case 1:
	  {
	    for (int lcnt = 0; lcnt < NLS; lcnt++) {
//     if(lcnt>0)continue;
	      sumElastic = 0.0; sumCapture = 0.0; sumFission = 0.0;
	      sumElasticM = 0.0; 
	      lVal = l[lcnt];
	      nrsl = NRS[lVal];
	      phil = calcPhi(x, lVal);
	      sfil = sin(phil);
	      cfil = cos(phil);
	      sfil2 = sfil*sfil;
	      //cfil2 = cfil*cfil;
	      s2fil = 2. * sfil *cfil ;
	      //c2fil = 1. - 2.*sfil2;
	      fac1 = 4.0 * (2.0 * lVal + 1.0) * sfil2;
	      fac2 = 4.0 * MisGj[l[lcnt]] * sfil2* pibyk2;
//     for (int jVal = 0; jVal <= Jnum; jVal++) {
//    for (int jVal = 0; jVal < NJValue[l[lcnt]]; jVal++) {
	      jVal = 0; nrstemp = nrs;  int nrend = nrs + nrsl;  
	      do
	      {
		facElastic =0.0; facCapture = 0.0; facFission = 0.0;
		for (int r = nrs; r < nrend; r++) {	
		  if(MissingJ[l[lcnt]][jVal] != (J[r])){
		    nrs = r; 
		    break;
		  }
		GNR = Gamma_nrE(x, r, lVal);  // Gamma_nr(E)	
		//GXR = Gamma_xrE(x, r, LRX);
		GR  = Gamma_rE(x, r, lVal, LRX);
		ERP = GetERP(x, r,  lVal);
		//ERPP = GetERP(x, r, lVal, 1);
		deno = (x-ERP) * (x-ERP) + 0.25 * GR * GR;
		fachNumElastic =  GNR*((GNR - 2.0 * GR * sfil2 + 2.0 * (x - ERP) * s2fil)/deno);
		fachNumCapture =  (GNR *Gamma_g[r] /deno);
		fachNumFission =  (GNR *Gamma_f[r] /deno);
		facElastic += GJ[r] * fachNumElastic;
		facCapture += GJ[r] * fachNumCapture;
		facFission += GJ[r] * fachNumFission;
//	  std::cout<< r <<"  "<< GJ[r] << "  "<< MissingJ[l[lcnt]][jVal]<< "  "<< NJValue[l[lcnt]]<< "  "<< x<<"  "<< Er[r] << std::endl;
	      }   
	      sumElastic += facElastic; 
	      sumCapture += facCapture; 
	      sumFission += facFission; 
	      jVal += 1;
	      }while(jVal < NJValue[l[lcnt]]);
	      nrs = nrstemp;
	      nrs += nrsl;
	      sumCrsElastic += pibyk2 * (fac1 + sumElastic) + fac2;
	      sumCrsCapture += pibyk2 * (sumCapture);
	      sumCrsFission += pibyk2 * (sumFission);
//   std::cout<<phil<<"  "<<std::setprecision(8)<< x << "  "<<std::setprecision(10)<<sumCrsElastic << "  "<<std::setprecision(10)<<sumCrsCapture<< "  "<<std::setprecision(10)<<sumCrsFission << std::endl;
	    }
//	   std::cout<< "E  "<< x << "  "<<sumCrsElastic <<"  "<< sumCrsCapture <<"  "<< sumCrsFission << std::endl;
//          energy.push_back(x);
//	  sigma.push_back(sumCrsElastic);
	    siga = sumCrsElastic;	    
	    sigb = sumCrsCapture;	    
	    sigc = sumCrsFission;	
	  
	  }break;
	case 2:
	  switch(LFW)
	  {
	    case 0:
	      {
//	    std::cout<<"case 2 nls2 "<< NLS2 <<std::endl;
		int nrsj = 0; nrs = 0; 
		for (int lcnt = 0; lcnt < NLS2; lcnt++) {
		  sumElastic = 0.0; sumCapture = 0.0; sumFission = 0.0;
		  sumElasticM = 0.0; 
		  lVal = l[lcnt];
	      //    lVal = l[NLS + lcnt];
		  phil = calcPhi(x, lVal);
		  sfil = sin(phil);
		  cfil = cos(phil);
		  sfil2 = sfil *sfil ;
		  //cfil2 = cfil * cfil;
		  s2fil = 2. * sfil *cfil ;
		  //c2fil = 1. - 2.*sfil2;
		  fac1 = 4.0 * (2.0 * lVal + 1.0) * sfil2;
		  NJS = NRJ[lcnt];
//	    std::cout<<"case 2 njs "<< NJS <<std::endl;
		  for (int jVal = 0; jVal < NJS; jVal++) {
		    nrsl = JSM[nrsj];
		    facElastic =0.0; facCapture = 0.0; facFission = 0.0; facElasticM = 0.0;
		    int min = nrs;
		    int max = nrs + nrsl - 1;
		    int mid = 0;
		    if (x <= Es[min])min = 0;
		    else if (x >= Es[max]) min = max - 1;
		    else {
		    while (max - min > 1) {
		      mid = (min + max) / 2;
		      if (x < Es[mid])
			max = mid;
		      else
			min = mid;
		      }
		    }
		    double gnr = GNO[min];
		    double gx = GX[min];
		    double gg = GG[min];
		    double gf = GF[min];
		    double dr = D[min];
       
		    GNR = gnr * calcPene(x, lVal)*sqrt(x)/GetRho(x, lVal);
		    widthFluctuation(GNR, gx, gg, gf, nrsj);
		    double temp1 =2.*PI*GNR * GJ[nrsj]/dr;
		    facElastic =  amun[nrsj]*(GNR*RN-2.*sfil2)*temp1;
		    facCapture =   amun[nrsj] * gg*temp1*RG;
		    facFission =   amun[nrsj] *  gf*temp1*RF;

		    sumElastic += facElastic; 
		    sumCapture += facCapture; 
		    sumFission += facFission; 
		    nrs += nrsl;
		    nrsj +=1;
		  }
		  sumCrsElastic += pibyk2 * (fac1 + sumElastic);
		  sumCrsCapture += pibyk2 * (sumCapture);
		  sumCrsFission += pibyk2 * (sumFission);
		}
//       std::cout<< x <<"  "<< sumCrsElastic <<"  "<< sumCrsCapture <<"  "<< sumCrsFission << std::endl;
		siga = sumCrsElastic;	    
		sigb = sumCrsCapture;	    
		sigc = sumCrsFission;	    
//          energy.push_back(x);
//	  sigma.push_back(sumCrsCapture);
	    
	      }break;
	    case 1:
	      {
//	    std::cout<<"case 2 nls2 "<< NLS2 <<std::endl;
	      int nrsj = 0; nrs = 0; 
	      for (int lcnt = 0; lcnt < NLS2; lcnt++) {
		sumElastic = 0.0; sumCapture = 0.0; sumFission = 0.0;
		sumElasticM = 0.0; 
		lVal = l[lcnt];
	    //   lVal = l[NLS + lcnt];
		phil = calcPhi(x, lVal);
		sfil = sin(phil);
		cfil = cos(phil);
		sfil2 = sfil *sfil ;
		//cfil2 = cfil * cfil;
		s2fil = 2. * sfil *cfil ;
		//c2fil = 1. - 2.*sfil2;
		fac1 = 4.0 * (2.0 * lVal + 1.0) * sfil2;
		NJS = NRJ[lcnt];
//	    std::cout<<"case 2 njs "<< NJS <<std::endl;
		for (int jVal = 0; jVal < NJS; jVal++) {
		  nrsl = JSM[nrsj];
		  facElastic =0.0; facCapture = 0.0; facFission = 0.0; facElasticM = 0.0;
		  int min = nrs;
		  int max = nrs + nrsl - 1;
		  int mid = 0;
		  if (x <= Es[min])min = 0;
		  else if (x >= Es[max]) min = max - 1;
		  else {
		  while (max - min > 1) {
		    mid = (min + max) / 2;
		    if (x < Es[mid])
		      max = mid;
		    else
		      min = mid;
		    }
		  }
//  std::cout<<"min "<<min<<"  "<<Es[min]<<"  "<<Es[min+1]<<" nrsl "<< nrsl <<" nrs "<< nrs <<std::endl;
		  double gnr = GNO[min];
		  double gx = GX[min];
		  double gg = GG[min];
		  double gf = GF[min];
		  double dr = D[min];
		  GNR = gnr * calcPene(x, lVal)*sqrt(x)/GetRho(x, lVal);
		  widthFluctuation(GNR, gx, gg, gf, nrsj);
		  double temp1 =2.*PI*GNR * GJ[nrsj]/dr;
		  facElastic =  amun[nrsj]*(GNR*RN-2.*sfil2)*temp1;
		  facCapture =   amun[nrsj] * gg*temp1*RG;
		  facFission =   amun[nrsj] *  gf*temp1*RF;

		  sumElastic += facElastic; 
		  sumCapture += facCapture; 
		  sumFission += facFission; 
		  nrs += nrsl;
		  nrsj +=1;
		}
		sumCrsElastic += pibyk2 * (fac1 + sumElastic);
		sumCrsCapture += pibyk2 * (sumCapture);
		sumCrsFission += pibyk2 * (sumFission);
	      }
//       std::cout<< x <<"  "<< sumCrsElastic <<"  "<< sumCrsCapture <<"  "<< sumCrsFission << std::endl;
	      siga = sumCrsElastic;	    
	      sigb = sumCrsCapture;	    
	      sigc = sumCrsFission;	    
//          energy.push_back(x);
//	  sigma.push_back(sumCrsCapture);
	    
	  }break;
	}
	break;
      }break;
// MLBW--------  
	  case 2:
	    {
	      switch(LRU)
	      {
		case 1:
		  {
		    for (int lcnt = 0; lcnt < NLS; lcnt++) {
		      sumElastic = 0.0; sumCapture = 0.0; sumFission = 0.0;
		      sumElasticM = 0.0; 
		      lVal = l[lcnt];
		  //    std::cout<<APL[l[lcnt]]<<std::endl;
		      nrsl = NRS[lVal];
		      phil = calcPhi(x, lVal);
		  //    std::cout<<"phil "<< phil <<" rho "<<GetRho(x, 0, lVal)/sqrt(x)<<" pE "<< calcPene(x, lVal) << std::endl;
		      sfil = sin(phil);
		      cfil = cos(phil);
		      sfil2 = sfil *sfil ;
		      //cfil2 = cfil * cfil;
		      s2fil = 2. * sfil *cfil ;
		      //c2fil = 1. - 2.*sfil2;
		      fac1 = 4.0 * (2.0 * lVal + 1.0) * sfil2;
		      fac2 = 4.0 * MisGj[l[lcnt]] * sfil2* pibyk2;
//     std::cout<<"beg "<< lcnt <<" E "<< x <<" nrsl "<< nrsl <<  std::endl;
//    std::cout<<"phi "<< phil <<" shift "<<calcShift(x, lVal)<<" p "<< calcPene(x, lVal) <<" E "<< x << std::endl;
//     for (int jVal = 0; jVal < Jnum; jVal++) {
//    for (int jVal = 0; jVal < NJValue[l[lcnt]]; jVal++) {
		      jVal = 0; nrstemp = nrs;   int nrend = nrs + nrsl; 
		      do
		      {
//      std::cout<<std::endl;
			facElastic =0.0; facCapture = 0.0; facFission = 0.0; facElasticM = 0.0;
			for (int r = nrs; r < nrend; r++) {	
			  fachNumElastic=0.0; fachNumCapture=0.0;  fachNumFission=0.0; fachNumElasticM = 0.0;
			  if(MissingJ[l[lcnt]][jVal] != (J[r])){
			    nrs = r; 
			    break;
			  }
			  GNR = Gamma_nrE(x, r, lVal);  // Gamma_nr(E)	
			  //GXR = Gamma_xrE(x, r, LRX);
			  GR  = Gamma_rE(x, r, lVal, LRX);
			  ERP = (GetERP(x, r,  lVal));
			  //ERPP = (GetERP(x, r, lVal, 0));
			  deno = (x - ERP) * (x - ERP) + 0.25 * GR *GR;
		  //	std::cout<<GNR<<"  "<<Gamma_g[r]<<"  "<<Gamma_f[r]<<"  "<<ERP<<"  "<<lVal<<" E "<<x<<std::endl; 
			  fachNumElastic =  GNR*((GNR - 2.0 * GR * sfil2 + 2.0 * (x - ERP) * s2fil)/deno);
			  fachNumCapture =  (GNR *Gamma_g[r] /deno);
			  fachNumFission =  (GNR *Gamma_f[r] /deno);
			  facElastic +=  GJ[r] * fachNumElastic;
			  facCapture +=  GJ[r] * fachNumCapture;
			  facFission +=  GJ[r] * fachNumFission;
			  fachNumElasticG = 0.0; fachNumElasticH = 0.0;
			  for (int rs = nrs; rs < nrend; rs++) {
		  //        if(jstart != J[rs])continue;	
			    if(MissingJ[l[lcnt]][jVal] != J[rs])continue;	
			    GNRS = Gamma_nrE(x, rs, lVal);
			    GS  = Gamma_rE(x, rs, lVal, LRX);
			    ERPP1 = (GetERP(x, rs,  lVal));
			    if(r != rs){
			      fachNumElasticG +=  0.5 * GNRS * GNR * (GR + GS)/((ERP -ERPP1) * (ERP -ERPP1) +0.25* (GR + GS) * (GR + GS));
			      fachNumElasticH +=  GNRS * GNR * (ERP -ERPP1)/((ERP -ERPP1) * (ERP -ERPP1) + 0.25 * (GR + GS) * (GR + GS));
			    }
			  }
			  fachNumElasticM = (fachNumElasticG * GR + 2.* fachNumElasticH * (x - ERP)) / ((x - ERP) * (x - ERP) + 0.25 * GR * GR); 
			  facElasticM += GJ[r] * fachNumElasticM;
//	 std::cout<<"facCapture "<<facCapture<<"  "<< GJ[r] * fachNumCapture <<"  "<< x <<" Er "<< Er[r] <<" " << GJ[r]<< std::endl;
			} 
			sumElastic += facElastic; 
			sumCapture += facCapture; 
			sumFission += facFission; 
			sumElasticM +=  facElasticM; 
//	 std::cout<<"sumCapture "<< sumCapture <<" facCapture "<< facCapture <<" E "<< x <<" jVal "<< jVal<<" lVal "<< lVal << std::endl;
			jVal += 1;
		      }while(jVal < NJValue[l[lcnt]]);
		      nrs = nrstemp;
		      nrs  += nrsl;
		      sumCrsElastic += pibyk2 * (fac1 + sumElastic +sumElasticM) + fac2;
		      sumCrsCapture += pibyk2 * (sumCapture);
		      sumCrsFission += pibyk2 * (sumFission);
//std::cout<<pibyk2<<"  "<< x << "  "<<pibyk2 * (fac1 + sumElastic +sumElasticM) <<"  "<< pibyk2 * (sumCapture) <<"  "<< pibyk2 * (sumFission) << std::endl;
		    }
//          energy.push_back(x);
//	  sigma.push_back(sumCrsElastic);
		    siga = sumCrsElastic;	    
		    sigb = sumCrsCapture;	    
		    sigc = sumCrsFission;	    
//std::cout<<"enegy "<< x << "  "<<sumCrsElasticM <<"  "<< sumCrsCapture <<"  "<< sumCrsFission << std::endl;
		  }break;
		case 2:
		  {
//	    std::cout<<"case 2 nls2 "<< NLS2 <<std::endl;
		    int nrsj = 0; nrs = 0; 
		    for (int lcnt = 0; lcnt < NLS2; lcnt++) {
		      sumElastic = 0.0; sumCapture = 0.0; sumFission = 0.0;
		      sumElasticM = 0.0; 
		      lVal = l[lcnt];
		  //    lVal = l[NLS + lcnt];
		      phil = calcPhi(x, lVal);
		      sfil = sin(phil);
		      cfil = cos(phil);
		      sfil2 = sfil *sfil ;
		      //cfil2 = cfil * cfil;
		      s2fil = 2. * sfil *cfil ;
		      //c2fil = 1. - 2.*sfil2;
		      fac1 = 4.0 * (2.0 * lVal + 1.0) * sfil2;
		      NJS = NRJ[lcnt];
//		      std::cout<<"case 2 njs "<< NJS <<" E= "<< x << std::endl;
		      for (int jVal = 0; jVal < NJS; jVal++) {
		      
			nrsl = JSM[nrsj];
			facElastic =0.0; facCapture = 0.0; facFission = 0.0; facElasticM = 0.0; 
			int min = nrs;
			int max = nrs + nrsl - 1;
			int mid = 0;
			if (x <= Es[min])min = 0;
			else if (x >= Es[max]) min = max - 1;
			else {
			  while (max - min > 1) {
			    mid = (min + max) / 2;
			    if (x < Es[mid])
			      max = mid;
			    else
			      min = mid;
			  }
			}
			  double gnr,gx,gg,gf,dr;
			  gnr = GNO[min];
			  gx = GX[min];
			  gg = GG[min];
			  gf = GF[min];
			  dr = D[min];
			  GNR = gnr * calcPene(x, lVal)*sqrt(x)/GetRho(x, lVal);
			  widthFluctuation(GNR, gx, gg, gf, nrsj);
			  double temp1 =2.*PI*GNR * GJ[nrsj]/dr;
			  facElastic =  amun[nrsj]*(GNR*RN-2.*sfil2)*temp1;
			  facCapture =   amun[nrsj] * gg*temp1*RG;
			  facFission =   amun[nrsj] *  gf*temp1*RF;

			  sumElastic += facElastic; 
			  sumCapture += facCapture; 
			  sumFission += facFission;
			nrs += nrsl;
			nrsj +=1;
		      }
		      sumCrsElastic += pibyk2 * (fac1 + sumElastic);
		      sumCrsCapture += pibyk2 * (sumCapture);
		      sumCrsFission += pibyk2 * (sumFission);
		    }
//       std::cout<< x <<"  "<< sumCrsElastic <<"  "<< sumCrsCapture <<"  "<< sumCrsFission << std::endl;
		    siga = sumCrsElastic;	    
		    sigb = sumCrsCapture;	    
		    sigc = sumCrsFission;	    
//          energy.push_back(x);
//	  sigma.push_back(sumCrsCapture);
		}break;
	    }break;
	  }break;
// Reich-Moore  
	  case 3:
	  {
	    GetSigmaRMP(x, siga,sigb,sigc);
	  }break;
// Adler-Adler	  
	  case 4:
	  {
	  }break;
// R-Matrix  
	  case 7:
	  {
	  }break;
	}
}
double TNudyEndfRecoPoint::backCrsAdler(double x, int nx){
  double crs1=0.0, crs2 = 0.0, backCrs = 0.0, pibyk2;
  double phil, sfil, sfil2, s2fil, c2fil, deno;
  phil = calcPhi(x, 0);
  sfil = sin(phil);
  sfil2 = x2(sfil);
  s2fil = sin(2.0 * phil);
  c2fil = cos(2.0 * phil);

  pibyk2  = PI/x2(K_wnum(x));
  backCrs = sqrt(x)*pibyk2/(at1[nx] + at2[nx]/x + at3[nx]/x2(x) + at4[nx]/(x*x*x) + bt1[nx] * x + bt2[nx] * x2(x));
  crs1    = 4.0 * pibyk2 * sfil2;
  switch(nx)
  {
    case 0:
      {
	for (int nrs = 0; nrs < totalAdler; nrs++){
	  deno    = (det1[nrs] - x)*(det1[nrs] - x) +  dwt1[nrs] *  dwt1[nrs] ;
	  crs2    = sqrt (x) * ( dwt1[nrs] * (grt1[nrs] * c2fil + git1[nrs] * s2fil) + (det1[nrs] - x) * (git1[nrs] * c2fil - grt1[nrs] * s2fil) )/deno;
	}
	crsAdler[nx] = crs1 + crs2 + backCrs;
      }break;
    case 1:
      {
	for (int nrs = 0; nrs < totalAdler; nrs++){
	  deno    = (def1[nrs] - x)*(def1[nrs] - x) +  dwf1[nrs] *  dwf1[nrs] ;
	  crs2    = sqrt (x) * ( dwf1[nrs] * grf1[nrs]  + (def1[nrs] - x) * gif1[nrs] )/deno;
	}
	crsAdler[nx] = crs2 + backCrs;
      }break;
    case 2:
      {
	for (int nrs = 0; nrs < totalAdler; nrs++){
	  deno    = (dec1[nrs] - x)*(dec1[nrs] - x) +  dwc1[nrs] *  dwc1[nrs] ;
	  crs2    = sqrt (x) * ( dwc1[nrs] * grc1[nrs]  + (dec1[nrs] - x) * gic1[nrs] )/deno;
	}
	crsAdler[nx] = crs2 + backCrs;
      }break;
  }
  return 0;
}

double TNudyEndfRecoPoint::recursionLinear(double x1, double x2, double sig1, double sig2, double sig3, double sig4, double sig5, double sig6){
  double siga, sigb, sigc, elRatio= sigDiff-0.1*sigDiff, capRatio= sigDiff-0.1*sigDiff, fisRatio= sigDiff-0.1*sigDiff;
  double mid     = 0.5 * (x1 + x2);
  if((sig1==0.0 && sig2 ==0.0 )|| x1==x2 || (x1 < 1E-5 || x2 < 1E-5))return 0;
//  std::cout<<" beg " << mid <<"  "<< x1 <<"  "<< x2 <<"  "<< sig1 <<"  "<< sig2 <<"  "<< sig3 <<"  "<< sig4<<"  "<< sig5 <<"  "<< sig6<< std::endl;
  GetSigma(LRF, mid, siga, sigb, sigc);
//  double sigmid1 = TNudyCore::Instance()->LinearInterpolation(x1,sig1,x2,sig2,mid);
  double sigmid1 = sig1 + (sig2 - sig1)*(mid - x1)/(x2 - x1);
  double sigmid2 = sig3 + (sig4 - sig3)*(mid - x1)/(x2 - x1);
  double sigmid3 = sig5 + (sig6 - sig5)*(mid - x1)/(x2 - x1);
  if(siga>0)elRatio = std::fabs((siga - sigmid1)/siga);
  if(sigb>0)capRatio = std::fabs((sigb - sigmid2)/sigb);
  if(sigc>0)fisRatio = std::fabs((sigc - sigmid3)/sigc);
//  if(elRatio <= sigDiff){
  if(elRatio <= sigDiff  && capRatio <= sigDiff && fisRatio <= sigDiff){
  return 0;
  }else{
//    if(elRatio > sigDiff){ 
      eLinElastic.push_back(mid);
      xLinElastic.push_back(siga);
//  std::cout<<" el1 " << mid <<"  "<< x1 <<"  "<< x2 <<"  "<< siga <<"  "<< sigb <<"  "<< sigc <<"  "<< elRatio <<"  "<< capRatio<<"  "<<fisRatio << std::endl;
//    }	  
//    if(capRatio > sigDiff){ 
      eLinCapture.push_back(mid);
      xLinCapture.push_back(sigb);
//  std::cout<<" cap " << mid <<"  "<< x1 <<"  "<< x2 <<"  "<< siga <<"  "<< sigb <<"  "<< sigc <<"  "<< elRatio <<"  "<< capRatio<<"  "<<fisRatio << std::endl;
//    }	  
//    if(fisRatio > sigDiff){ 
      eLinFission.push_back(mid);
      xLinFission.push_back(sigc);
//  std::cout<<" fis " << mid <<"  "<< x1 <<"  "<< x2 <<"  "<< siga <<"  "<< sigb <<"  "<< sigc <<"  "<< elRatio <<"  "<< capRatio<<"  "<<fisRatio << std::endl;
//    }
  }
  recursionLinear(x1, mid, sig1, siga, sig3, sigb, sig5, sigc);
  recursionLinear(mid, x2, siga, sig2, sigb, sig4, sigc, sig6);
  return 0;
}
double TNudyEndfRecoPoint::recursionLinear(double x1, double x2, double sig1, double sig2){
  double siga, sigb, sigc, elRatio;//= sigDiff-0.1*sigDiff;
  double mid     = 0.5 * (x1 + x2);
  if((sig1==0.0 && sig2 ==0.0 )|| x1==x2 || (x1 < 1E-5 || x2 < 1E-5))return 0;
//  std::cout<<" x1 " << mid <<"  "<< x1 <<"  "<< x2 <<"  "<< sig1 <<"  "<< sig2 << std::endl;
  GetSigma(LRF, mid, siga, sigb, sigc);
//  double sigmid1 = TNudyCore::Instance()->LinearInterpolation(x1,sig1,x2,sig2,mid);
  double sigmid1 = sig1 + (sig2 - sig1)*(mid - x1)/(x2 - x1);
  elRatio=std::fabs((siga - sigmid1)/siga);
//  std::cout<<" mid " << mid <<"  "<< x1 <<"  "<< x2 <<"  "<< elRatio <<"  "<< capRatio<<"  "<<fisRatio << std::endl;
  if(elRatio <= sigDiff){
    return 0;
  }else{
//  std::cout<<" el1 " << mid <<"  "<< x1 <<"  "<< x2 <<"  "<< siga <<"  "<< sigb <<"  "<< sigc <<"  "<< elRatio << std::endl;
    eLinElastic.push_back(mid);
    xLinElastic.push_back(siga);
    eLinCapture.push_back(mid);
    xLinCapture.push_back(sigb);
    eLinFission.push_back(mid);
    xLinFission.push_back(sigc);
  }
  recursionLinear(x1, mid, sig1, siga);
  recursionLinear(mid, x2, siga, sig2);
  return 0;
}

void TNudyEndfRecoPoint::additionalSigma(int LRF, double x){
  double siga, sigb, sigc;
  GetSigma(LRF, x, siga, sigb, sigc);
//    std::cout<<x<<"  "<<siga<<"  "<<sigb<<"  "<<sigc<<std::endl;
  eLinElastic.push_back(x);
  eLinCapture.push_back(x);
  eLinFission.push_back(x);
  xLinElastic.push_back(siga);
  xLinCapture.push_back(sigb);
  xLinFission.push_back(sigc);
}
void TNudyEndfRecoPoint::Sort(std::vector<double>& x1, std::vector<double>& x2){
  std::multimap<double, double>map;
  std::multimap<double, double>::iterator i;
  for(int p = 0; p< x1.size(); p++)
    map.insert(std::make_pair(x1[p],x2[p]));
    int p1=0;
    for(i=map.begin(); i!=map.end();i++){
      x1[p1]=i->first;
      x2[p1]=i->second;
      p1++;
    }
}
 void TNudyEndfRecoPoint::recoPlusBroad(int flagNer){
  int nrs = 0, nrsl;
  std::cout<<" reco called LRU "<< LRU <<" LSSF "<< LSSF << std::endl;
  if (LRU==1){
    for (int j = 0; j < NLS; j++) {
      nrsl = NRS[l[j]]; 
      for (int i = nrs; i < nrs + NRS[l[j]]; i++) {
	if (Er[i] < eLo || Er[i] > eHi1){
	  continue;
	}
//	      std::cout<<"L = "<< j <<" r " << i <<" energy = "<< Er[i] <<" nrs "<< nrs << " nrsl "<< nrsl <<" Er size "<< Er.size() << std::endl;
	additionalSigma(LRF, Er[i]);
      }
      nrs += nrsl;
    }
    double del = 0.001*eHi1;
    additionalSigma(LRF, eHi1-del);
//    std::cout << eLinElastic.size() <<std::endl;
    if(flagNer == 0){
      double eneLow = 1.0E-5;
      do
      {
	if (eneLow>=eLo)additionalSigma(LRF, eneLow);
	  eneLow *= 1.15;
      }while(eneLow < 1.0E3 && eneLow < eHi1 );
    } 
//    std::cout << eLinElastic.size() <<std::endl;
    Sort(eLinElastic, xLinElastic);
    Sort(eLinCapture, xLinCapture);
    Sort(eLinFission, xLinFission);
    ThinningDuplicate(eLinElastic, xLinElastic);
    ThinningDuplicate(eLinCapture, xLinCapture);
    ThinningDuplicate(eLinFission, xLinFission);
    int nvectorend = eLinElastic.size();
//	std::cout<<"eLo1 "<< eLo1 <<" eHi1 "<< eHi1 <<"   "<< eLinElastic[nvectorend-1] << std::endl;
    for(int ju = intLinLru1; ju < nvectorend-1 ; ju++){
//	  std::cout<<"ju  "<<ju<<"  "<<eLinElastic[ju]<<"  "<<eLinElastic[ju+1]<<"  "<< xLinElastic[ju]<<"  "<<xLinElastic[ju+1]<<"  "<< xLinCapture[ju]<<"  "<< xLinCapture[ju+1]<<"  "<< xLinFission[ju]<<"  "<< xLinFission[ju+1]<< std::endl;
      recursionLinear(eLinElastic[ju], eLinElastic[ju+1], xLinElastic[ju], xLinElastic[ju+1], xLinCapture[ju], xLinCapture[ju+1], xLinFission[ju], xLinFission[ju+1]);
      intLinLru1 = eLinElastic.size();
    }
  }else{
  intLinLru1 = eLinElastic.size();
  /*  
  int unresoPoint = (int)(Es[1] - Es[0]);
  double es = eLo2 ,esdel = Es[0]/unresoPoint/10;
  double del = 0.001;
  additionalSigma(LRF, es+del);
  double del1 = 0.1;
  es += del1;
  do
  { 
    additionalSigma(LRF, es);
    es += esdel;
  }while(es <= Es[1]);
  additionalSigma(LRF, eHi2);
  */
  for(int l = 0; l<JSM[0]; l++){
    additionalSigma(LRF, Es[l]); 
  }
  
    //  for (int i =0;i<eLinCapture.size(); i++)
    //  std::cout<<"k1 "<<eLinCapture[i]<<"  "<<xLinCapture[i]<<std::endl;
//	std::cout<<"before "<< eLinElastic[intLinLru1] <<"  "<<eLinElastic[eLinElastic.size()-1]<<"  "<<intLinLru1<<"  "<<eLinElastic.size()<<std::endl;
//  recursionLinear(eLinElastic[intLinLru1], eLinElastic[eLinElastic.size()-1], xLinElastic[intLinLru1], xLinElastic[eLinElastic.size()-1]);
    //  for (int i =0;i<eLinCapture.size(); i++)
    //  std::cout<<"k2 "<<eLinCapture[i]<<"  "<<xLinCapture[i]<<std::endl;
    //    std::cout<<"vector size LRF"<< LRF <<" = "<<eLinElastic.size()<<std::endl;
    //    xSort();
    //      dopplerBroadning(t11, t12); 
  }
}
//************  Uniform Linearization ***************
void TNudyEndfRecoPoint::Linearize(int flagNer) {
  double a = 0.0;
  switch(LRF)
  {
    case 1:
    case 2:
    case 3:
    {
     recoPlusBroad(flagNer);
    }break;
    case 4:
    {
    for (int i = 1; i < (int)eHi; i++) {
      a = (double)i;
      if(eHi>100000.0)a = (double)i*10.;
      if(eHi<10000.0)a = (double)i/10.;
      if(eHi<5000.0)a = (double)i/100.;
      a += eLo;
      crsAdler[0] = backCrsAdler(a,0);
      crsAdler[1] = backCrsAdler(a,1);
      crsAdler[2] = backCrsAdler(a,2);
      crsAdler[3] = crsAdler[0] - crsAdler[1] - crsAdler[2];
    }
    }break;
  }
}
//******************** Get Resonant PArameter and cross-section data from ENDF file **********
void TNudyEndfRecoPoint::GetData(const char *rENDF) {
  //TVNudyModel *model;
  //TGeoElementRN *mat;
  std::vector<double> rowI;
//  if (!gGeoManager) gGeoManager = new TGeoManager("rENDF Nudy Manager","");
//  TNudyManager *nudy = TNudyManager::Instance();
//  nudy->OpenDatabase("endf",rENDF);
//  nudy->AddEndfLibrary("endf",fENDF);
//  nudy->LoadLibrary("endf","endf");
  
//  TVNudyModel *model = nudy->GetModel(13,27,0,2,300.0,"neutron");
//  std::cout<< "end of get model " << std::endl;

  TFile *rEND = TFile::Open(rENDF);
  if (!rEND || rEND->IsZombie()) printf("Error: TFile :: Cannot open file %s\n", rENDF);
  TKey *rkey = (TKey*)rEND->GetListOfKeys()->First();
  TNudyEndfTape *rENDFVol = (TNudyEndfTape*)rkey->ReadObj();
  TNudyEndfMat *tMat = 0; 
  TList *mats = (TList*)rENDFVol->GetMats(); 
  
  // I am checking for Char_t* as Root Browser is showing as TString
  //std::cout << " I am checking whether fname and fTitle crashes or not for Char_t *...\n";
  //Char_t *theKeyName = (Char_t*)rENDFVol->GetName();
  //Char_t *theTitle   = (Char_t*)rENDFVol->GetTitle();
  //std::cout << " Key Name :" << theKeyName << "\nTitle: " << theTitle << std::endl;
  //std::cout << "-------------------------------------------------------\n";
  //std::cout << " I am checking whether fname and fTitle crashes or not for TString ...\n";
  //TString strKeyName = rENDFVol->GetName();
  //TString strTitle   = rENDFVol->GetTitle();
  //std::cout << " Using TString fName: " << strKeyName << "\nTitle: " << strTitle << std::endl;
 //std::cout << "-------------------------------------------------------\n";
  
  int nmats = mats->GetEntries();
  for (int iMat = 0; iMat < nmats; iMat++) {
    tMat = (TNudyEndfMat*)mats->At(iMat);
//    std::cout<<" MAT number "<< tMat->GetMAT() <<"  "<< nmats << std::endl;
    TIter iter(tMat->GetFiles());
    TNudyEndfFile *file;
    while ((file = (TNudyEndfFile *)iter.Next())) {
    //std::cout <<" mf "<< file->GetMF() << std::endl;
    // Read File data into class structures
      switch (file->GetMF()) 
      {
        case 1:
	  ReadFile1(file);
	  break;
        case 2:
	  ReadFile2(file);
//      std::cout<<" file2 finished "<< std::endl;
	  if(LRU !=0){
	    Sort(eLinElastic, xLinElastic);
	    Sort(eLinCapture, xLinCapture);
	    Sort(eLinFission, xLinFission);
	    ThinningDuplicate(eLinElastic, xLinElastic);
	    ThinningDuplicate(eLinCapture, xLinCapture);
	    ThinningDuplicate(eLinFission, xLinFission);
	    std::cout<<"Elastic points  "<< eLinElastic.size() << std::endl;
	    std::cout<<"Capture points  "<< eLinCapture.size() << std::endl;
	    std::cout<<"Fission points  "<< eLinFission.size() << std::endl;
//	    for(int cr1=0; cr1 < eLinElastic.size() ; cr1 ++){
//	      std::cout<<"Elastic  "<<eLinElastic[cr1]<<"  "<< xLinElastic[cr1] << std::endl;
//	    }
	  }
	  std::cout<<"file 2 OK "<< std::endl;
	  break;
        case 3:
	  ReadFile3(file);
///*	  
	  std::cout<<"file 3 OK "<<std::endl;
	  //std::cout<<"Total energy points "<< energyUni.size()<<std::endl;
	  std::sort(energyUni.begin(), energyUni.end());//Unionization of energy grid for cross-section 
	  //for(int j=0;j<energyUni.size();j++){
	  //  std::cout<<energyUni[j]<<std::endl;
	  //}
	  ThinningDuplicate(energyUni);
	  std::cout<<"Total energy points "<< energyUni.size()<<std::endl;
//      if(LRU==0){
	  std::cout<<"NoOfElements "<< MtValues.size() <<" MT values "<< MtValues[0].size() << std::endl;
	  for (int i = 0; i < MtValues.size(); i++){
	    for (int j = 0; j < MtValues[i].size(); j++){
	      std::cout <<"MT value "<< MtValues[i][j] << std::endl;
	    }
	  }
	  std::cout<<"NoOfsigma "<< sigmaOfMts.size() <<" No of points of 1, 2 "<< sigmaOfMts[0].size() <<"  "<< sigmaOfMts[1].size()<< std::endl;
	  for (int i = 0; i < sigmaOfMts.size(); i++){
	    int size = sigmaOfMts[i].size()/2;
//	    std::cout<<"NoOfsigma before "<< size << std::endl;
	    for (int k = 0; k < energyUni.size(); k++){
	      int min = 0;
	      int max = size - 1;
	      int mid = 0;
	      if (energyUni[k] <= sigmaOfMts[i][min])min = 0;
	      else if (energyUni[k] >= sigmaOfMts[i][max]) min = max - 1;
	      else {
		while (max - min > 1) {
		  mid = (min + max) / 2;
		  if (energyUni[k] < sigmaOfMts[i][mid])
		    max = mid;
		  else
		    min = mid;
		  }
	      }
//	      std::cout<<"energy "<< energyUni[k] << std::endl;
	      if(energyUni[k] == sigmaOfMts[i][min] && sigmaOfMts[i][size + min] > 1E-20){
		sigmaMts.push_back(energyUni[k]);
		sigmaMts.push_back(sigmaOfMts[i][size + min]);
//		std::cout<<"energy =  "<< energyUni[k] << std::endl;
	      }else{
		double sigmaAdd = sigmaOfMts[i][size + min] + (sigmaOfMts[i][size + min + 1] - sigmaOfMts[i][size + min]) 
				* (energyUni[k]-sigmaOfMts[i][min])/(sigmaOfMts[i][min + 1] - sigmaOfMts[i][min]); //linear interpolation
		if(sigmaAdd > 1E-20){
		  sigmaMts.push_back(energyUni[k]);
		  sigmaMts.push_back(sigmaAdd);
//		    std::cout <<"energy!= "<< energyUni[k]  << std::endl;
		}
//		std::cout <<"energy "<< sigmaOfMts[i][j] <<" sigma "<<sigmaOfMts[i][size + j] << std::endl;
	      }
	      if(sigmaMts.size() == 2){
//		std::cout <<" first location "<< k << std::endl;
		energyLocationMts.push_back(k);
	      }
	    }
//	    std::cout<<"NoOfsigma after "<< sigmaMts.size() << std::endl;
//	    for(int l = 0 ; l < sigmaMts.size(); l++){
//	      std::cout << l <<"   "<< sigmaMts[l] << std::endl;
//	    }
	    sigmaUniOfMts.push_back(sigmaMts);
	    sigmaMts.clear();
	  }
	  sigmaUniTotal.resize(energyUni.size());
	  for (int i = 0; i < sigmaUniOfMts.size(); i++){
	    int size = sigmaUniOfMts[i].size();
//	    std::cout<<"NoOfsigma after add "<< size/2 <<" location "<< energyLocationMts[i] << " MT = "<< MtValues[0][i] << std::endl;
	    int MT = MtValues[0][i];
	    if(MT != 1 && MT != 3 && MT != 4 && MT != 27 && MT != 19 && MT != 20 && MT != 21 && MT != 38 && MT != 101 && MT < 250){
	      for (int j = 0; j < size/2; j++){
		sigmaUniTotal[energyLocationMts[i] + j] += sigmaUniOfMts[i][2*j + 1];
//		std::cout <<"MT= "<< MT <<"  "<< energyUni[energyLocationMts[i]+j] << "  "<< sigmaUniOfMts[i][2*j] <<"  "<< sigmaUniOfMts[i][2*j + 1]<< std::endl;
	      }
	    }
	  }
	  sigma.clear();
	  //std::cout<<" AWRI "<< AWRI <<"   "<< GetAWR() << std::endl;
	  if(energyUni.size()>0){
	    doppler = new TNudyEndfDoppler(AWRI,0.0, 293.6, energyUni, sigmaUniTotal);
	    //dopplerBroadning(0.0, 293.6, energyUni, sigmaUniTotal);
	    std::cout<<"Total Doppler done "<<doppler->sigma.size()<<std::endl;
	    dopplerBroad=1;
	    sigmaUniTotal.clear();
	    for(int j=0; j< energyUni.size(); j++){
	      sigmaUniTotal.push_back(doppler->sigma[j]);
	    }
	    doppler->sigma.clear();
	  }
	  outtotal << energyUni.size() << std::endl;
	  for(int j=0;j<energyUni.size();j++){
	    outtotal << energyUni[j] << "  " << sigmaUniTotal[j] << std::endl;
	  }
	    std::cout<<"file 3 OK "<<std::endl;
	    sigma.clear();
  //      dopplerBroadning(0.0, 293.6,  energy, sigmaT);
  //      for(int j =0; j< energy.size(); j++)
  //	if(fabs((sigmaT[j]-sigma[j])*100/sigmaT[j]>1.0))std::cout <<"Total "<< energy[j] <<"  "<<sigmaT[j]<<"  "<< sigma[j]<<"   "<<(sigmaT[j]-sigma[j])*100/sigmaT[j]<< std::endl;
  //      energy.clear(); sigmaT.clear();sigma.clear();
  // std::cout<<LRU <<"  "<< LSSF<<std::endl;     
	    std::cout<<"sorting before Doppler begins "<<std::endl;
	    Sort(eLinElastic, xLinElastic);
	    Sort(eLinCapture, xLinCapture);
	    Sort(eLinFission, xLinFission);
	    std::cout<<"sorting before Doppler done "<<std::endl;
  //           if(LRU !=0){
  //           Thinning(eLinElastic, xLinElastic);
  //           Thinning(eLinCapture, xLinCapture);
  //           Thinning(eLinFission, xLinFission);
  //	  }
  //	  for(int j=0; j< eLinElastic.size(); j++){
  //	    std::cout<<eLinElastic[j]<<"  "<<xLinElastic[j]<<std::endl;
  //	  }
	    std::cout<<"Elstic Doppler begins "<<std::endl;
	    if(eLinElastic.size()>0){
	      doppler = new TNudyEndfDoppler(AWRI,0.0, 293.6, eLinElastic, xLinElastic);
	      //dopplerBroadning(0.0, 293.6, eLinElastic, xLinElastic);
	      std::cout<<"Elstic Doppler done "<<std::endl;
	      dopplerBroad=1;
	      for(int j=0; j< eLinElastic.size(); j++){
  //	     energy.push_back(eLinElastic[j]);
		xBroadElastic.push_back(doppler->sigma[j]);}
		doppler->sigma.clear();
  //	   for(int j=0; j< eLinElastic.size(); j++)
  //	     std::cout<<eLinElastic[j]<<"  "<<xBroadElastic[j]<<std::endl;
		Thinning(eLinElastic, xBroadElastic);
	      }
	      std::cout<<"Capture Doppler begins "<<std::endl;
	      if(eLinCapture.size()>0){
		doppler = new TNudyEndfDoppler(AWRI,0.0, 293.6, eLinCapture, xLinCapture);
		//dopplerBroadning(0.0, 293.6, eLinCapture, xLinCapture);
		std::cout<<"Capture Doppler done "<<std::endl;
		for(int j=0; j< eLinCapture.size(); j++){
  //	     energy.push_back(eLinCapture[j]);
		  xBroadCapture.push_back(doppler->sigma[j]);
		}
		doppler->sigma.clear();
		Thinning(eLinCapture, xBroadCapture);
	      }
	      std::cout<<"Fission Doppler begins "<<std::endl;
	      if(eLinFission.size()>0){
		doppler = new TNudyEndfDoppler(AWRI,0.0, 293.6, eLinFission, xLinFission);
		//dopplerBroadning(0.0, 293.6, eLinFission, xLinFission);
		std::cout<<"Fission Doppler done "<<std::endl;
		for(int j=0; j< eLinFission.size(); j++){
  //	     energy.push_back(eLinFission[j]);
		  xBroadFission.push_back(doppler->sigma[j]);
		}
		doppler->sigma.clear();
		Thinning(eLinFission, xBroadFission);
	      }
	      dopplerBroad=0;
  //      std::sort(energy.begin(),energy.end());
	    out << eLinElastic.size() << std::endl;
	    for(int j =0; j< eLinElastic.size(); j++)
	      out << eLinElastic[j] <<"  "<< xBroadElastic[j] << std::endl;
	      out << eLinCapture.size() << std::endl;
	    for(int j =0; j< eLinCapture.size(); j++)
	      out << eLinCapture[j] <<"  "<< xBroadCapture[j] << std::endl;
	      out << eLinFission.size() << std::endl;
	    for(int j =0; j< eLinFission.size(); j++)
	      out << eLinFission[j] <<"  "<< xBroadFission[j] << std::endl;
//*/     
//    }
	  break;
        case 4:
	  recoAng = new TNudyEndfAng(file);
	  break;
        case 5:
	  ReadFile5(file);
	  break;
        case 6:
	  ReadFile6(file);
	  break;
        case 8:
//      ReadFile8(file);
	  break;
      }
    }

//     model->ReadFile(tMat);
    //if (!strcmp(tMat->GetName(),"90-Th-232")) {       //94-Pu-238  90-Th-232
    //if (tMat->GetMAT()==9434) {
    if (tMat->GetMAT()==9434) {
      std::cout << " Found Material " << tMat->GetName() << " < MAT # " << tMat->GetMAT() << ">" << std::endl;
      break;
    }
  }
}	 
double TNudyEndfRecoPoint::GetSigmaTotal(double energyK){
  int min = 0;
  int max = energyUni.size() - 1;
  int mid = 0;
  if (energyK <= energyUni[min])min = 0;
  else if (energyK >= energyUni[max]) min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < energyUni[mid]) max = mid;
      else min = mid;
    }
  }
  return sigmaUniTotal[min] + (sigmaUniTotal[min + 1] - sigmaUniTotal[min])*
	 (energyK - energyUni[min])/(energyUni[min + 1] - energyUni[min]);
}
double TNudyEndfRecoPoint::ThinningDuplicate(std::vector<double> &x1){
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
double TNudyEndfRecoPoint::ThinningDuplicate(std::vector<double> &x1,std::vector<double> &x2){
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

double TNudyEndfRecoPoint::Thinning(std::vector<double> &x1, std::vector<double> &x2){
  int size = x1.size();
  int size1 = x1.size();
  if(size<=0)return 0;
  for(int i = 0; i< size1 - 2; i++){
    if(x1[i]<1E3 && dopplerBroad==0)continue;
    if(x1[i] > eHi && LRU != 0 && eHi > 0)continue;
  //  double sigmid1 = TNudyCore::Instance()->LinearInterpolation(x1,sig1,x2,sig2,mid);
    double sigmid1 = x2[i] + (x2[i+2] - x2[i])*(x1[i+1] - x1[i])/(x1[i+2] - x1[i]);
//  std::cout<<" mid " << mid <<"  "<< siga <<"  "<< sigmid1 << std::endl;
  
    if(fabs((x2[i+1] - sigmid1)/sigmid1)<=sigDiff){
      x1.erase(x1.begin()+i+1);
      x2.erase(x2.begin()+i+1);
    } 
    if(x2[i]<0.0){
      x1.erase(x1.begin()+i);
      x2.erase(x2.begin()+i);
    }  
  size1 = x1.size();
  }
  size1 = x1.size();
  if(size == size1)return 0;
  Thinning(x1, x2);
  return 0;
}

void TNudyEndfRecoPoint::fixupTotal(std::vector<double> &x1, std::vector<double> &x2){
  int size = energy.size();
  int size1 = x1.size();
  double sigmid1=0.0;
  if (fE_file3) {    delete[] fE_file3;    fE_file3 = 0;  }
  if (fXsec_file3) {    delete[] fXsec_file3;    fXsec_file3 = 0;  }
  fE_file3    = new double[x1.size()]();
  fXsec_file3 = new double[x1.size()]();
  for(int cr=0; cr < size1 ; cr ++){
    fE_file3[cr] = x1[cr]; 
    fXsec_file3[cr] = x2[cr];
  }
  for(int i = 0; i< size; i++){
    int index = TNudyCore::Instance()->BinarySearch(fE_file3, size1, energy[i]);
//    std::cout<<"index "<< index <<"  "<<std::setprecision(10)<< fE_file3[index] <<"  "<<std::setprecision(10)<< fE_file3[index+1]<< std::endl;
    if(energy[i]>=fE_file3[0]){
      sigmid1 = TNudyCore::Instance()->LinearInterpolation(fE_file3[index],fXsec_file3[index],fE_file3[index+1],fXsec_file3[index+1],energy[i]);
// std::cout<<i<<"  "<<std::setprecision(10)<< energy[i]<<"  "<<std::setprecision(10)<<sigmid1<<"  "<<std::setprecision(10)<< sigma[i] << std::endl;
      sigma[i] += sigmid1;
// std::cout<<i<<"  "<<std::setprecision(10)<< energy[i]<<"  "<<std::setprecision(10)<<sigmid1<<"  "<<std::setprecision(10)<< sigma[i] << std::endl;
    }
  }
  if (fE_file3) {    delete[] fE_file3;    fE_file3 = 0;  }
  if (fXsec_file3) {    delete[] fXsec_file3;    fXsec_file3 = 0;  }
}
TNudyEndfRecoPoint::~TNudyEndfRecoPoint(){}
