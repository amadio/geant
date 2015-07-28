/*
   This is the main class supporting an ENDF section in R-ENDF format

*/

// @(#)root/meta:$Id: TNuEndf.h 29000 2009-06-15 13:53:52Z rdm $
// Author: F.Carminati 02/05/09

/*
   This is the main class supporting an ENDF material in R-ENDF format

*/

#include <string.h>

#include <TString.h>
#include <TNudyEndfCont.h>

ClassImp(TNudyEndfCont)

//_______________________________________________________________________________
TNudyEndfCont::TNudyEndfCont() :
  fC1(0),
  fC2(0),
  fL1(0),
  fL2(0),
  fN1(0),
  fN2(0)
{
  //
  // Default constructor
  //
}

//_______________________________________________________________________________
TNudyEndfCont::TNudyEndfCont(Double_t c1, Double_t c2,
			     Int_t l1, Int_t l2, Int_t n1, Int_t n2) :
  fC1(c1),
  fC2(c2),
  fL1(l1),
  fL2(l2),
  fN1(n1),
  fN2(n2)
{
  //
  // Standard constructor
  //
}


//_______________________________________________________________________________
void TNudyEndfCont::SetCont(Double_t c1, Double_t c2,
			    Int_t l1, Int_t l2, Int_t n1, Int_t n2)
{
  fC1=c1;
  fC2=c2;
  fL1=l1;
  fL2=l2;
  fN1=n1;
  fN2=n2;
}

//
// Dump Data to screen in ENDF format
//______________________________________________________________________________
void TNudyEndfCont::DumpENDF(Int_t mat,Int_t mf, Int_t mt,Int_t& ns, Int_t flags = 1)
{
  Char_t s1[14],s2[14];
  F2F(fC1,s1); F2F(fC2,s2);
  printf("%11s%11s%11d%11d%11d%11d", s1,s2, fL1,fL2, fN1,fN2);
  printf("%4d%2d%3d%5d", mat, mf, mt, ns);
  if(flags)
    printf("  ---CONT\n");
  else
    printf("\n");
  if (ns < 99999)
    ns++;
  else
    ns =1;
}

//
//Float_t to Fortran style string, dim of s should be 14 or bigger
//______________________________________________________________________________
Char_t * TNudyEndfCont::F2F(Double_t f, char s[])
{
  snprintf(s,14,"%12.5E",f);
  if(s[10]!='0') {
    s[8]=s[9];
    s[9]=s[10];
      s[10]=s[11];
  } else {
    snprintf(s,14,"%13.6E",f);
    s[9]= s[10];
      s[10]= s[12];
  }
  s[11]= '\0';
   return s;
}
