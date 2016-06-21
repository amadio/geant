// 	This class is reconstructing Fission yields  
// 	Author: Dr. Harphool Kumawat
// 	Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// 	date of creation: March 22, 2016

#include "TList.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfList.h"
#include "TNudyEndfFissionYield.h"

TNudyEndfFissionYield::TNudyEndfFissionYield(){}

//______________________________________________________________________________
TNudyEndfFissionYield::TNudyEndfFissionYield(TNudyEndfFile *file)
{
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    TIter recIter(sec->GetRecords());
    //double ZA   = sec->GetC1();
    //double AWR  = sec->GetC2();
    int MT = sec->GetMT();
    int LE = sec->GetL1();
    if(MT == 454){// Neutron induced independent fission yield
      for (int i =0; i < LE; i++){
        TNudyEndfList *list1 = (TNudyEndfList *)recIter.Next();
	ein.push_back ( list1->GetC1() );
	//int NN   = list1->GetN1();
	int NFP  = list1->GetN2();
//      std::cout<<"energy i " <<ein[i] << std::endl;
	for (int j = 0; j < NFP; j++){
	  zafp1.push_back (list1->GetLIST(4*j+0));
	  fps1.push_back(list1->GetLIST(4*j+1));
	  yi1.push_back(list1->GetLIST(4*j+2));
	  dyi1.push_back(list1->GetLIST(4*j+3));
	}
	zafp.push_back (zafp1);
	fps.push_back(fps1);
	yi.push_back(yi1);
	dyi.push_back(dyi1);
	zafp1.clear();
	fps1.clear();
	yi1.clear();
	dyi1.clear();	  
      }
    }else if(MT == 459){// Neutron induced cummulative fission yield
	for (int i =0; i < LE; i++){
	  TNudyEndfList *list1 = (TNudyEndfList *)recIter.Next();
	  einc.push_back (list1->GetC1());
	  //int NN   = list1->GetN1();
	  int NFP  = list1->GetN2();
//          std::cout<<"energy c " <<einc[i] << std::endl;
	  for (int j = 0; j < NFP; j++){
	    zafpc1.push_back(list1->GetLIST(4*j+0));
	    fpsc1.push_back(list1->GetLIST(4*j+1));
	    yc1.push_back(list1->GetLIST(4*j+2));
	    dyc1.push_back(list1->GetLIST(4*j+3));
//	    std::cout<<"fission yield c "<< zafpc[i].At(j) << std::endl; 
	  }
	  zafpc.push_back (zafpc1);
	  fpsc.push_back(fpsc1);
	  yc.push_back(yc1);
	  dyc.push_back(dyc1);
          zafpc1.clear();
	  fpsc1.clear();
	  yc1.clear();
	  dyc1.clear();	  	    
	}
      }
   }
}

TNudyEndfFissionYield::~TNudyEndfFissionYield(){}
