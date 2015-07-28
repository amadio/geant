/*
   This is the main class supporting an ENDF section in R-ENDF format

*/

// @(#)root/meta:$Id: TNuEndf.h 29000 2009-06-15 13:53:52Z rdm $
// Author: F.Carminati 02/05/09

#include <string.h>

#include <TString.h>
#include <TNudyEndfTab2.h>

ClassImp(TNudyEndfTab2)

//_______________________________________________________________________________
TNudyEndfTab2::TNudyEndfTab2() :
  TNudyEndfCont(),
  fNBT(NULL),
  fINT(NULL)
{
  //
  // Default constructor
  //
}

//_______________________________________________________________________________
TNudyEndfTab2::TNudyEndfTab2(Double_t c1, Double_t c2,
			     Int_t l1, Int_t l2, Int_t n1, Int_t n2) :
  TNudyEndfCont(c1, c2, l1, l2, n1, n2),
  fNBT(new Int_t[n1]),
  fINT(new Int_t[n1])
{
  //
  // Standard constructor
  //
}

//______________________________________________________________________________
TNudyEndfTab2::~TNudyEndfTab2()
{
  //printf("Deleting Tab2\n");
  SafeDelete(fNBT);
  SafeDelete(fINT);
}

//_______________________________________________________________________________
void TNudyEndfTab2::SetCont(Double_t c1, Double_t c2,
			    Int_t l1, Int_t l2, Int_t n1, Int_t n2)
{
  TNudyEndfCont::SetCont(c1, c2, l1, l2, n1, n2);
  delete [] fNBT;
  delete [] fINT;
  fNBT=new Int_t[n1];
  fINT=new Int_t[n1];
}

void TNudyEndfTab2::DumpENDF(Int_t mat, Int_t mf, Int_t mt, Int_t& ns,Int_t flags = 1)
{
  //Print Tab2 CONT Record
  Char_t s1[14],s2[14];
  F2F(fC1,s1); F2F(fC2,s2);
  printf("%11s%11s%11d%11d%11d%11d", s1,s2, fL1,fL2, fN1,fN2);
  printf("%4d%2d%3d%5d", mat, mf, mt, ns);
  if (ns < 99999)
    ns++;
  else
    ns =1;
  if(flags)
    printf("  ---CONT TAB2\n");
  else
    printf("\n");
    
  for(Int_t i=0; i<GetNR(); i++) {                      //print NBT(N) INT(N)
    if(i%3==0 && i!=0)  {
      printf("%4d%2d%3d%5d", mat, mf, mt, ns);
      if (ns < 99999)
	ns++;
      else
	ns =1;
      if(flags)
	printf("  ---NBT(%d,%d,%d) TAB2\n",i-2,i-1,i);
      else
	printf("\n");
    }
    printf("%11d%11d", GetNBT(i), GetINT(i));
  }
  //Pad blank columns
  if(3-(GetNR()%3) < 3){
    for(Int_t i = 0; i < 3-(GetNR()%3);i++) {
      printf("%22s"," ");
    }
  }
  printf("%4d%2d%3d%5d", mat, mf, mt, ns);
  if (ns < 99999)
    ns++;
  else
    ns =1;
  if(flags)
    printf("  ---NBT(%d,%d,%d) TAB2\n",(GetNR()-(GetNR()%3))+1,(GetNR()-(GetNR()%3))+2,(GetNR()-(GetNR()%3))+3);
  else
    printf("\n");
}

