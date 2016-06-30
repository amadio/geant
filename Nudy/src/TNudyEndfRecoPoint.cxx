// This class is reconstructing ENDF data and creating probability tables for the angle
// and energy distributions of the secondatries
// Author: Dr. Harphool Kumawat
// Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// date of creation: March 22, 2016

#include "TList.h"
#include "TFile.h"
#include "TKey.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfTape.h"
#include "TNudyEndfList.h"
#include "TNudyEndfTab1.h"
#include "TNudyEndfMat.h"
#include "TNudyEndfNuPh.h"
#include "TNudyEndfDoppler.h"
#include "TNudyEndfAng.h"
#include "TNudyEndfEnergy.h"
#include "TNudyEndfEnergyAng.h"
#include "TNudyEndfFissionYield.h"
#include "TNudyCore.h"
#include "TNudyEndfRecoPoint.h"

#ifdef USE_ROOT
ClassImp(TNudyEndfRecoPoint)
#endif

#ifdef USE_ROOT
#include "TRandom.h"
#endif

TNudyEndfRecoPoint::TNudyEndfRecoPoint() : elemId(0),rENDF(),sigDiff(0){}
TNudyEndfRecoPoint::TNudyEndfRecoPoint(int ielemId, const char *irENDF,double isigDiff)
:elemId(ielemId),
rENDF(irENDF),
sigDiff(isigDiff)
{
  GetData(ielemId,irENDF,isigDiff);
}

//______________________________________________________________________________
void TNudyEndfRecoPoint::ReadFile2(TNudyEndfFile *file) {
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  NLS = 0;double gjdeno = 1; LSSF=0;
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
		    gjfound += (2 *std::fabs(jtemp[j1]) + 1 )/gjdeno;	
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
		      gjfound += (2 *std::fabs(jtemp[j1]) + 1 )/gjdeno;	
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
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::recursionLinearFile3(double x1, 
						double x2, 
						double sig1, 
						double sig2, 
						std::vector<double> ene, 
						std::vector<double> sig){
  double siga;
  double mid     = 0.5 * (x1 + x2);
//  std::cout <<"beg "<< x1 <<"  "<< x2 <<"  "<< mid <<"  "<< sig1 <<"  "<< sig2 <<std::endl;
  if(sig1==0.0 && sig2 ==0.0)return 0;
  if(x1==x2 || x1 < 1E-5 || x2 < 1E-5){
    return 0;
  }
  siga = TNudyCore::Instance()->Interpolate(nbt1, int1, NR, ene,sig, NP, mid);
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
  recursionLinearFile3(x1, mid, sig1, siga, ene, sig);
  recursionLinearFile3(mid, x2, siga, sig2, ene, sig);
  return 0;  
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::addFile3Resonance(double &x1, 
					     double &x2, 
				std::vector<double> &x3, 
				std::vector<double> &x4){
  int linsize = x3.size();
  int size1 = eLinearFile3.size();
  for(int cr1=0; cr1 < linsize ; cr1 ++){
    if(x3[cr1] >= x1 && x3[cr1] <= x2){
      double crs = TNudyCore::Instance()->Interpolate(nbt1,int1,NR,eLinearFile3,xLinearFile3,size1,x3[cr1]);
      x4[cr1] += crs;
    }
  }
  return 0;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::insertFile3(std::vector<double> &x1, 
				       std::vector<double> &x2){
  int size1 = eLinearFile3.size();
  for(int cr1=0; cr1 < size1 ; cr1 ++){
    if((eLinearFile3[cr1] >eLo2 && eLinearFile3[cr1] <=eHi2) || eLinearFile3[cr1] >eHi1){
      x1.push_back(eLinearFile3[cr1]);
      x2.push_back(xLinearFile3[cr1]);
    }
  }
  return 0;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::insertFile3High(std::vector<double> &x1, 
					   std::vector<double> &x2){
  if(flagResolve != 0){
    for(unsigned long cr=0; cr < eLinearFile3.size() ; cr ++){
      if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] <=eLo1){
	x1.push_back(eLinearFile3[cr]);
	x2.push_back(xLinearFile3[cr]);
      }
    }
  }else if(flagResolve == 0 && flagUnResolve != 0){
    for(unsigned long cr=0; cr < eLinearFile3.size() ; cr ++){
      if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] <=eLo2){
	x1.push_back(eLinearFile3[cr]);
	x2.push_back(xLinearFile3[cr]);
      }
    }
  }
  if(flagUnResolve != 0){
    for(unsigned long cr=0; cr < eLinearFile3.size() ; cr ++){
      if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] >eHi2){
	x1.push_back(eLinearFile3[cr]);
	x2.push_back(xLinearFile3[cr]);
      }
    }
  }else if(flagUnResolve == 0){
    for(unsigned long cr=0; cr < eLinearFile3.size() ; cr ++){
      if(xLinearFile3[cr]>0.0 && eLinearFile3[cr] >eHi1){
	x1.push_back(eLinearFile3[cr]);
	x2.push_back(xLinearFile3[cr]);
      }
    }
  }
  return 0;
}

//------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::ReadFile3(TNudyEndfFile *file) {
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  std::vector<double>eneExtra;
  if(LRF==0)eHi=0.0;
  double eneLow = 1E-5;
  do
    {
      eneExtra.push_back(eneLow);
      eneLow *= 1.15;
    }while(eneLow < 1E3);
  
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    std::cout <<"MT  "<<sec->GetMT() <<std::endl;
    int MT = sec->GetMT();
    if(MT != 1 && MT != 3 && MT != 4 && MT != 27 && MT != 19 && MT != 20 && MT != 21 && MT != 38 && MT != 101 && MT < 250){
      MtNumbers.push_back(MT);
    TIter recIter(sec->GetRecords());
    TNudyEndfCont *header = (TNudyEndfCont *)recIter.Next();
    //double ZA   = sec->GetC1();
    //    double AWRI  = sec->GetC2();
    //double AWR  = sec->GetC2();
    //double QM   = header->GetC1();
    double QI   = header->GetC2();
    QValue[sec->GetMT()] = QI;
    SetQValue(QI,MT);
    //  std::cout<< "file3 "<< QValue[sec->GetMT()] <<"  "<< GetQValue(MT) << std::endl; 
      
    //int LR = header->GetL2();
    NR = header->GetN1();
    NP = header->GetN2();
//std::cout<<sec->GetMT() <<"  "<< LR <<"  "<< NR <<"  "<< NP << std::endl;      
    TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
    for(int cr=0; cr < NR ; cr ++){
      nbt1.push_back (tab1->GetNBT(cr));
      int1.push_back (tab1->GetINT(cr));
//      TNudyEndfTab1 *tab1 = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
//   std::cout << tab1->GetNBT(cr) <<" NBT "<< tab1->GetINT(cr)<< std::endl;
//   std::cout << tab1->GetINT(cr)<< std::endl;
    }
    for(int crs=0; crs < NP ; crs ++){
      eneTemp.push_back(tab1->GetX(crs));
      sigTemp.push_back(tab1->GetY(crs));
      if(eneTemp[crs]==eneTemp[crs-1])eneTemp[crs] += 0.001;//adding small difference for boundary points
      if(eneTemp[crs]==eneTemp[crs-2])eneTemp[crs] += 0.002;//adding small difference for double boundary points
//	std::cout<<"hi "<< npp <<"  "<< fE_file3[npp]<< std::endl;
    }
    TNudyCore::Instance()->Sort(eneTemp, sigTemp);
    TNudyCore::Instance()->ThinningDuplicate(eneTemp, sigTemp);
    NP = eneTemp.size();
    for(int crs=0; crs < NP ; crs ++){
      eLinearFile3.push_back(eneTemp[crs]);
      xLinearFile3.push_back(sigTemp[crs]);
    }
//	std::cout << NBT[0] <<" "<< NBT[1] <<"  "<< INT1[0] <<"  "<< INT1[1] <<"  "<< eLinearFile3.size() << std::endl;
// Linearization of file 3 data
    for(int cr=0; cr < NP - 1 ; cr ++){
//	  std::cout<<"cr = "<<cr<<"  "<<eLinearFile3[cr] <<"  "<< xLinearFile3[cr] << std::endl;
      recursionLinearFile3(eLinearFile3[cr], eLinearFile3[cr+1], xLinearFile3[cr], xLinearFile3[cr+1], eLinearFile3, xLinearFile3);
    }
    for(unsigned long ie = 0; ie < eneExtra.size(); ie++){
      double sigExtra = TNudyCore::Instance()->Interpolate(nbt1, int1, NR, eLinearFile3,xLinearFile3, NP, eneExtra[ie]);
      if(sigExtra > 0.0){
	eLinearFile3.push_back(eneExtra[ie]);
	xLinearFile3.push_back(sigExtra);
      }
//	  std::cout<<"extra= "<<MT<<"  "<<std::setprecision(12)<<eneExtra[ie] <<"  "<<std::setprecision(12)<< sigExtra << std::endl;
    }
    TNudyCore::Instance()->Sort(eLinearFile3, xLinearFile3);
//    for(int cr=0; cr < eLinearFile3.size() ; cr ++){
//      std::cout<<"1Mt = "<<MT<<"  "<<eLinearFile3[cr] <<"  "<< xLinearFile3[cr] << std::endl;
//    }    
    TNudyCore::Instance()->ThinningDuplicate(eLinearFile3, xLinearFile3);
	
//    for(int cr=0; cr < eLinearFile3.size() ; cr ++){
//      std::cout<<"2Mt = "<<MT<<"  "<<eLinearFile3[cr] <<"  "<< xLinearFile3[cr] << std::endl;
//    }    
//Filling of array to interpolate and add with file 2 data if it is given    
    nbt1.clear(); int1.clear();
    
    nbt1.push_back(eLinearFile3.size()); 
    int1.push_back(2);
    NR = 1;
/////////////////////////////////////////////////////////////////////////    
//	std::cout<<"LSSF "<< LSSF <<" LRU "<< LRU <<" eLo1 "<<eLo1<<" eHi "<< eHi << std::endl;
//	  for(unsigned long cr=0; cr < eLinearFile3.size() ; cr ++)
//	  std::cout<<"Mt = "<<MT<<"  "<<std::setprecision(12)<<eLinearFile3[cr] <<"  "<<std::setprecision(12)<< xLinearFile3[cr] << std::endl;

// resolved range is always added in file 3
    if(LRU != 0){
      if(flagResolve != 0){
	switch(MT)
	{
	  case 2:
	  {
	    addFile3Resonance(eLo1, eHi1, eLinElastic, xLinElastic);  
	  }break;
	  case 18:
	  {
	    addFile3Resonance(eLo1, eHi1, eLinFission, xLinFission);  
	  }break;
	  case 102:
	  {
	    addFile3Resonance(eLo1, eHi1, eLinCapture, xLinCapture);  
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
	      addFile3Resonance(eLo2, eHi2, eLinElastic, xLinElastic);  
	    }break;
	    case 18:
	    {
	      addFile3Resonance(eLo2, eHi2, eLinFission, xLinFission);
	    }break;
	    case 102:
	    {
	      addFile3Resonance(eLo2, eHi2, eLinCapture, xLinCapture);
	    }break;
	  }
	}break;
	case 1:
	{
	  switch(MT)
	  {
	    case 2:
	    {
	      insertFile3(eLinElastic, xLinElastic); 
	    }break;
	    case 18:
	    {
	      insertFile3(eLinFission, xLinFission); 
	    }break;
	    case 102:
	    {
	      insertFile3(eLinCapture, xLinCapture); 
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
	  insertFile3High(eLinElastic, xLinElastic); 
	}break;
       case 18:
	{
	  insertFile3High(eLinFission, xLinFission); 
      }break;
      case 102:
	{
	  insertFile3High(eLinCapture, xLinCapture); 
	}break;
      }
  }else{// no resonance parameters are given
    switch(MT)
    {
      case 2:
	{
	  for(unsigned long cr=0; cr < eLinearFile3.size() ; cr ++){
	    eLinElastic.push_back(eLinearFile3[cr]);
	    xLinElastic.push_back(xLinearFile3[cr]);
	  }
	}break;
      case 18:
	{
	  for(unsigned long cr=0; cr < eLinearFile3.size() ; cr ++){
	    eLinFission.push_back(eLinearFile3[cr]);
	    xLinFission.push_back(xLinearFile3[cr]);
	  }
	}break;
      case 102:
	{
	  for(unsigned long cr=0; cr < eLinearFile3.size() ; cr ++){
	    eLinCapture.push_back(eLinearFile3[cr]);
	    xLinCapture.push_back(xLinearFile3[cr]);
	  }
	}break;
      }
    }    
  
// check the sum of cross sections with total cross-section given in the MT=1 + resonance
//      int size = eLinearFile3.size();
      
//      for(unsigned long i = 0; i < eLinearFile3.size(); i++){
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
    eLinearFile3.clear();	
    xLinearFile3.clear();
    eneTemp.clear();	
    sigTemp.clear();	
    nbt1.clear(); int1.clear();
    
//    for(unsigned long i = 0; i < sigmaMts.size(); i++){
//      std::cout<<sigmaMts[i]<<std::endl;
//    }
    
    sigmaOfMts.push_back(sigmaMts);
    sigmaMts.clear();
    }
  }
  MtValues.push_back(MtNumbers);
  MtNumbers.clear();
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::K_wnum(double x) {
  double k; 
  k = factor_k * sqrt(std::fabs(x));
  return k;
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetRho(double x1, int lVal) {
  if (!NAPS && !NRO) return factor_k * sqrt(std::fabs(x1)) * rad_a;   	// Use rad_a in the penetrabilities Pl and shift factors Sl , and AP in the hard-sphere phase shifts φl 
  if (!NAPS && NRO) return factor_k * sqrt(std::fabs(x1)) *  rad_a;    	// Use rad_a in the penetrabilities Pl and shift factors Sl,AP(E) in the hard-sphere phase shifts φl
  if ( NAPS && !NRO) return factor_k * sqrt(std::fabs(x1)) * APL[lVal];      	// use AP in the penetrabilities and shift factor as well as in the phase shifts
  if ( NAPS && NRO) return factor_k * sqrt(std::fabs(x1)) * APL[lVal];       	// read AP(E) and use it in all three places, Pl , Sl , φl
  return 0;  
}
 
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetRhoC(double x, int isDiff, int lVal) {
  if (!isDiff) return factor_k * sqrt(std::fabs(x)) * APL[lVal];    		//read AP(E) and use it in all three places, Pl , Sl , φl
  if (isDiff) return factor_k * sqrt(std::fabs(x)) * APL[lVal];                	//read AP(E) and use it in all three places, Pl , Sl , φl
  return 0;  
} 

//------------------------------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------------------------------
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


//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetERP(double x, int r, int lVal) {
  double er = Er[r];
  if(lVal==0)return er;
  return er + (ShiftEr[r] - calcShift(x, lVal))/(2.0 * PhiEr[r]) * Gamma_nrE(std::fabs(er), r, lVal);
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::Gamma_nrE(double x, int ii, int lval) {
  return  calcPene(x, lval) * Gamma_n[ii];
}

//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::Gamma_xrE(int ii, int lrx) {
  if (!lrx) {
    return 0;
  } else {
    return Gamma_r[ii] - (Gamma_n[ii] + Gamma_g[ii] + Gamma_f[ii]);
  }
}
//------------------------------------------------------------------------------------------------------

double TNudyEndfRecoPoint::Gamma_rE(double x, int ii, int lval, int lrx) {
  return Gamma_nrE(x, ii, lval) + Gamma_xrE(ii, lrx) + Gamma_g[ii] + Gamma_f[ii];
}
//------------------------------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------------------------------
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
  }
  siga = sumCrsElastic;	    
  sigb = sumCrsCapture;	    
  sigc = sumCrsFission;	    
}
//------------------------------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::recursionLinear(double x1, double x2, double sig1, double sig2, double sig3, double sig4, double sig5, double sig6){
  double siga, sigb, sigc, elRatio= sigDiff-0.1*sigDiff, capRatio= sigDiff-0.1*sigDiff, fisRatio= sigDiff-0.1*sigDiff;
  double mid     = 0.5 * (x1 + x2);
  if((sig1==0.0 && sig2 ==0.0 )|| x1==x2 || (x1 < 1E-5 || x2 < 1E-5))return 0;
  GetSigma(LRF, mid, siga, sigb, sigc);
  double sigmid1 = sig1 + (sig2 - sig1)*(mid - x1)/(x2 - x1);
  double sigmid2 = sig3 + (sig4 - sig3)*(mid - x1)/(x2 - x1);
  double sigmid3 = sig5 + (sig6 - sig5)*(mid - x1)/(x2 - x1);
  if(siga>0)elRatio = std::fabs((siga - sigmid1)/siga);
  if(sigb>0)capRatio = std::fabs((sigb - sigmid2)/sigb);
  if(sigc>0)fisRatio = std::fabs((sigc - sigmid3)/sigc);
  if(elRatio <= sigDiff  && capRatio <= sigDiff && fisRatio <= sigDiff){
  return 0;
  }else{
  //std::cout<<elRatio <<"  "<<capRatio<<"  "<< fisRatio <<"  "<<std::setprecision(6)<< mid <<"  "<< x1 <<"  "<< x2 << std::endl;
    //if(elRatio > sigDiff){ 
      eLinElastic.push_back(mid);
      xLinElastic.push_back(siga);
//  std::cout<<" el1 " << mid <<"  "<< x1 <<"  "<< x2 <<"  "<< siga <<"  "<< sigb <<"  "<< sigc <<"  "<< elRatio <<"  "<< capRatio<<"  "<<fisRatio << std::endl;
    //}	  
    //if(capRatio > sigDiff){ 
      eLinCapture.push_back(mid);
      xLinCapture.push_back(sigb);
//  std::cout<<" cap " << mid <<"  "<< x1 <<"  "<< x2 <<"  "<< siga <<"  "<< sigb <<"  "<< sigc <<"  "<< elRatio <<"  "<< capRatio<<"  "<<fisRatio << std::endl;
    //}	  
    //if(fisRatio > sigDiff){ 
      eLinFission.push_back(mid);
      xLinFission.push_back(sigc);
//  std::cout<<" fis " << mid <<"  "<< x1 <<"  "<< x2 <<"  "<< siga <<"  "<< sigb <<"  "<< sigc <<"  "<< elRatio <<"  "<< capRatio<<"  "<<fisRatio << std::endl;
    //}
  }
  recursionLinear(x1, mid, sig1, siga, sig3, sigb, sig5, sigc);
  recursionLinear(mid, x2, siga, sig2, sigb, sig4, sigc, sig6);
  return 0;
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::recursionLinear(double x1, double x2, double sig1, double sig2){
  double siga, sigb, sigc, elRatio; //= sigDiff-0.1*sigDiff;
  double mid     = 0.5 * (x1 + x2);
  if((sig1==0.0 && sig2 ==0.0 )|| x1==x2 || (x1 < 1E-5 || x2 < 1E-5))return 0;
  GetSigma(LRF, mid, siga, sigb, sigc);
  double sigmid1 = sig1 + (sig2 - sig1)*(mid - x1)/(x2 - x1);
  elRatio=std::fabs((siga - sigmid1)/siga);
  if(elRatio <= sigDiff){
    return 0;
  }else{
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

//------------------------------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------------------------------
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
    TNudyCore::Instance()->Sort(eLinElastic, xLinElastic);
    TNudyCore::Instance()->Sort(eLinCapture, xLinCapture);
    TNudyCore::Instance()->Sort(eLinFission, xLinFission);
    TNudyCore::Instance()->ThinningDuplicate(eLinElastic, xLinElastic);
    TNudyCore::Instance()->ThinningDuplicate(eLinCapture, xLinCapture);
    TNudyCore::Instance()->ThinningDuplicate(eLinFission, xLinFission);
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
//------------------------------------------------------------------------------------------------------
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
void TNudyEndfRecoPoint::GetData(int ielemId, const char *rENDF, double isigDiff) {
  elemId = ielemId;
  sigDiff = isigDiff;
  TFile *rEND = TFile::Open(rENDF);
  if (!rEND || rEND->IsZombie()) printf("Error: TFile :: Cannot open file %s\n", rENDF);
  TKey *rkey = (TKey*)rEND->GetListOfKeys()->First();
  TNudyEndfTape *rENDFVol = (TNudyEndfTape*)rkey->ReadObj();
  TNudyEndfMat *tMat = 0; 
  TList *mats = (TList*)rENDFVol->GetMats();
  int nmats = mats->GetEntries();
  for (int iMat = 0; iMat < nmats; iMat++) {
    tMat = (TNudyEndfMat*)mats->At(iMat);
  for (int inf = 0; inf < tMat->GetNXC(); inf++){
    if(tMat->GetMFn(inf)>3){
      MtNumbers.push_back(tMat->GetMFn(inf));
      MtNumbers.push_back(tMat->GetMTn(inf));
    }
  }
  Mt456.push_back(MtNumbers);
  std::vector<int>().swap(MtNumbers);
//    std::cout<<" MAT number "<< tMat->GetMAT() <<"  "<< nmats << std::endl;
    TIter iter(tMat->GetFiles());
    TNudyEndfFile *file;
    while ((file = (TNudyEndfFile *)iter.Next())) {
    std::cout <<" mf "<< file->GetMF() << std::endl;
    // Read File data into class structures
      switch (file->GetMF()) 
      {
        case 1:
	  recoNuPh = new TNudyEndfNuPh(file);
	  std::cout<<"file 1 OK "<< std::endl;
	  break;
        case 2:
	  ReadFile2(file);
          //      std::cout<<" file2 finished "<< std::endl;
	  if(LRU !=0){
	    TNudyCore::Instance()->Sort(eLinElastic, xLinElastic);
	    TNudyCore::Instance()->Sort(eLinCapture, xLinCapture);
	    TNudyCore::Instance()->Sort(eLinFission, xLinFission);
	    TNudyCore::Instance()->ThinningDuplicate(eLinElastic, xLinElastic);
	    TNudyCore::Instance()->ThinningDuplicate(eLinCapture, xLinCapture);
	    TNudyCore::Instance()->ThinningDuplicate(eLinFission, xLinFission);
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
	  sigma.clear();
	  std::cout<<"file 3 OK "<<std::endl;
	  fixupTotal(energyUni, sigmaUniTotal);
	  eneUni.push_back(energyUni);
	  sigUniT.push_back(sigmaUniTotal);
	  std::cout<<"Union Total OK "<<  energyUni.size()  <<std::endl;
	  energyUni.clear();
	  sigmaUniTotal.clear();
	  sigma.clear();
  // std::cout<<LRU <<"  "<< LSSF<<std::endl; 
	  std::cout<<"before Doppler begins "<<std::endl;
	  broadSigma(eLinElastic, xLinElastic, xBroadElastic);
	  broadSigma(eLinCapture, xLinCapture, xBroadCapture);
	  broadSigma(eLinFission, xLinFission, xBroadFission);
	  std::cout<<"Doppler done "<<std::endl;
	  dopplerBroad=0;

	  out << eLinElastic.size() << std::endl;
	  for(unsigned long j =0; j< eLinElastic.size(); j++)
	    out << eLinElastic[j] <<"  "<< xBroadElastic[j] << std::endl;
	    eLinElastic.clear(); xLinElastic.clear(); xBroadElastic.clear();
	    out << eLinCapture.size() << std::endl;
	  for(unsigned long j =0; j< eLinCapture.size(); j++)
	    out << eLinCapture[j] <<"  "<< xBroadCapture[j] << std::endl;
	    eLinCapture.clear(); xLinCapture.clear(); xBroadCapture.clear();
	    out << eLinFission.size() << std::endl;
	  for(unsigned long j =0; j< eLinFission.size(); j++)
	    out << eLinFission[j] <<"  "<< xBroadFission[j] << std::endl;
	    eLinFission.clear(); xLinFission.clear(); xBroadFission.clear();
//*/     
//    }
	  break;
        case 4:
	  recoAng = new TNudyEndfAng(file);
	  break;
        case 5:
	  recoEnergy = new TNudyEndfEnergy(file);
	  break;
        case 6:
	  recoEnergyAng = new TNudyEndfEnergyAng(file,QValue);
	  break;
        case 8:
	  recoFissY = new TNudyEndfFissionYield(file);
	  break;
      }
    }
  }
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::broadSigma(std::vector<double> &x1, 
				    std::vector<double> &x2, 
				    std::vector<double> &x3){
  TNudyCore::Instance()->Sort(x1, x2);
  if(x1.size()>0){
    doppler = new TNudyEndfDoppler(AWRI,0.0, 293.6, x1, x2);
    dopplerBroad=1;
    for(unsigned long j=0; j< x1.size(); j++){
      x3.push_back(doppler->sigma[j]);
    }
    doppler->sigma.clear();
    //Thinning(x1, x3);
  }
}
//------------------------------------------------------------------------------------------------------
void TNudyEndfRecoPoint::fixupTotal(std::vector<double> &x1, 
				    std::vector<double> &x2){
  
  std::sort(x1.begin(), x1.end());//Unionization of energy grid for cross-section 
  TNudyCore::Instance()->ThinningDuplicate(x1);
  std::cout<<"Total energy points "<< x1.size()<<std::endl;
  //for (unsigned int i = 0; i < x1.size(); i++)
    //std::cout << std::setprecision(12) << x1[i] << std::endl;
  //std::cout<<"NoOfElements "<< MtValues.size() << std::endl;
  // for (unsigned long i = 0; i < MtValues.size(); i++){
  //   for (unsigned long j = 0; j < MtValues[i].size(); j++){
  //     std::cout <<"MT value "<< MtValues[i][j] << std::endl;
  //   }
  // }
   std::cout<<"NoOfsigma "<< sigmaOfMts.size() << std::endl;
  for (unsigned long i = 0; i < sigmaOfMts.size(); i++){
    int size = sigmaOfMts[i].size()/2;
    for (unsigned long k = 0; k < x1.size(); k++){
      int min = 0;
      int max = size - 1;
      int mid = 0;
      if (x1[k] <= sigmaOfMts[i][min])min = 0;
      else if (x1[k] >= sigmaOfMts[i][max]) min = max - 1;
      else {
	while (max - min > 1) {
	  mid = (min + max) / 2;
	  if (x1[k] < sigmaOfMts[i][mid])
	    max = mid;
	  else
	    min = mid;
	}
      }
      if(x1[k] == sigmaOfMts[i][min] && sigmaOfMts[i][size + min] > 1E-20){
	eneTemp.push_back(x1[k]);
	sigTemp.push_back(sigmaOfMts[i][size + min]);
      }else{
 	double sigmaAdd = sigmaOfMts[i][size + min] + (sigmaOfMts[i][size + min + 1] - sigmaOfMts[i][size + min]) 
			* (x1[k]-sigmaOfMts[i][min])/(sigmaOfMts[i][min + 1] - sigmaOfMts[i][min]); //linear interpolation
	if(sigmaAdd > 1E-20){
	  eneTemp.push_back(x1[k]);
	  sigTemp.push_back(sigmaAdd);
	}
      }
      if(eneTemp.size() == 1){
	energyLocationMts.push_back(k);
      }
    }
    broadSigma(eneTemp, sigTemp, sigma);
    sigmaUniOfMts.push_back(sigma);
    eneTemp.clear(); sigTemp.clear(); sigma.clear();
  }
  x2.resize(x1.size());
  sigmaOfMts.clear();
  ///*
  for (unsigned long i = 0; i < sigmaUniOfMts.size(); i++){
    int size = sigmaUniOfMts[i].size();
      //std::cout<<"size "<<  energyLocationMts.size  << std::endl;
    for (int j = 0; j < size; j++){
      x2[energyLocationMts[i] + j] += sigmaUniOfMts[i][j];
	//std::cout <<"MT= "<< MtValues[elemId][i] <<"  "<< x1[energyLocationMts[i]+j] << "  "<< sigmaUniOfMts[i][j]<< std::endl;
    }
  }
 // */
 sigUniOfMt.push_back(sigmaUniOfMts);
 energyLocMtId.push_back(energyLocationMts);
 sigmaUniOfMts.clear();
 sigma.clear();
 energyLocationMts.clear();
  //std::cout<<" AWRI "<< AWRI <<"   "<< GetAWR() << std::endl;
  outtotal << x1.size() << std::endl;
  for(unsigned long j = 0; j < x1.size(); j++){
    outtotal << x1[j] << "  " << x2[j] << std::endl;
  }
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetSigmaTotal(int ielemId, double energyK){
  int min = 0;
  int max = eneUni[ielemId].size() - 1;
  int mid = 0;
  if (energyK <= eneUni[ielemId][min])min = 0;
  else if (energyK >= eneUni[ielemId][max]) min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < eneUni[ielemId][mid]) max = mid;
      else min = mid;
    }
  }
  return sigUniT[ielemId][min] + (sigUniT[ielemId][min + 1] - sigUniT[ielemId][min])*
	 (energyK - eneUni[ielemId][min])/(eneUni[ielemId][min + 1] - eneUni[ielemId][min]);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetSigmaPartial(int ielemId, int i, double energyK){
  int min = 0;
  int max = eneUni[ielemId].size() - 1;
  int mid = 0;
  if (energyK <= eneUni[ielemId][min])min = 0;
  else if (energyK >= eneUni[ielemId][max]) min = max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (energyK < eneUni[ielemId][mid]) max = mid;
      else min = mid;
    }
  }
  //std::cout<<min <<"  "<< max <<"  "<<energyK<<"  "<< eneUni[ielemId][min] << std::endl;
  return sigUniOfMt[ielemId][i][energyLocMtId[ielemId][i]+min] + (sigUniOfMt[ielemId][i][energyLocMtId[ielemId][i]+min + 1]
         - sigUniOfMt[ielemId][i][energyLocMtId[ielemId][i]+min])*
	 (energyK - eneUni[ielemId][min])/(eneUni[ielemId][min + 1] - eneUni[ielemId][min]);
}

//------------------------------------------------------------------------------------------------------
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
  
    if(std::fabs((x2[i+1] - sigmid1)/sigmid1) <= sigDiff){
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

//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetCos6(int ielemId, int mt, double energyK){
  return recoEnergyAng->GetCos4(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetEnergy6(int ielemId, int mt, double energyK){
  return recoEnergyAng->GetEnergy5(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetCos4(int ielemId, int mt, double energyK){
  return recoAng->GetCos4(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetEnergy5(int ielemId, int mt, double energyK){
  return recoEnergy->GetEnergy5(ielemId, mt, energyK);
}
//------------------------------------------------------------------------------------------------------
int TNudyEndfRecoPoint::GetLaw6(int ielemId, int mt){
  return recoEnergyAng->GetLaw6(ielemId, mt);
}
//-------------------------------------------------------------------------------------------------------
double TNudyEndfRecoPoint::GetMt456(int ielemId, int mt){
  int size = Mt456[ielemId].size()/2;
  for (int i = 0; i < size; i++){
    if(Mt456[ielemId][ 2 * i + 1] == mt){
      return  Mt456[ielemId][ 2 * i ] ;
    }
  }
  return 0;
}
TNudyEndfRecoPoint::~TNudyEndfRecoPoint(){}
