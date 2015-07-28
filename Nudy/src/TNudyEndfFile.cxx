/*
   This igs the main class supporting an ENDF file in R-ENDF format

*/

// @(#)root/meta:$Id: TNuEndf.h 29000 2009-06-15 13:53:52Z rdm $
// Author: F.Carminati 02/05/09

/*
   This is the main class supporting an ENDF material in R-ENDF format

*/

#include <string.h>

#include <TString.h>
#include <TNudyEndfFile.h>

ClassImp(TNudyEndfFile)

//_______________________________________________________________________________
TNudyEndfFile::TNudyEndfFile() :
  fMAT(0),
  fMF(0)
{
  //
  // Default constructor
  //
  strcpy(fName,"");
  fSecs = new TList();
}

//_______________________________________________________________________________
TNudyEndfFile::TNudyEndfFile(Int_t mat, Int_t mf) :
  fMAT(mat),
  fMF(mf)
{
  //
  // Standard constructor
  //
  snprintf(fName,8,"%04d-%02d",mat,mf);
  fSecs = new TList();
}

//______________________________________________________________________________
TNudyEndfFile::~TNudyEndfFile()
{
  //  printf("Destroying File %s\n", fName);
  fSecs->Delete();
  SafeDelete(fSecs);
}

//_______________________________________________________________________________
TNudyEndfSec* TNudyEndfFile::GetSec(Int_t MT)
{
	for(Int_t i=0;i<=this->fSecs->LastIndex();i++)
	{
		TNudyEndfSec *thisSec = (TNudyEndfSec*)this->fSecs->At(i);
		if(thisSec->GetMT()==MT)
			return thisSec;
	}
	Error("TNudyEndfFile::GetSec(Int_t)","Could not find section %d on tape",MT);
	return NULL;
}


//_______________________________________________________________________________
void TNudyEndfFile::DumpENDF(Int_t flags = 1)
{
	//Sections
	for(Int_t i=0;i<=fSecs->LastIndex();i++)
	{
		TNudyEndfSec *sec = (TNudyEndfSec*)fSecs->At(i);
		sec->DumpENDF(flags);
	}
	//FEND
	//	cout<<setw(66)<<" "<<setw(4)<<fMAT<<setw(2)<<"0"<<setw(3)<<"0"<<setw(5)<<"0"<<endl;
	printf("%66s%4d%2d%3d%5d"," 0.000000+0 0.000000+0          0          0          0          0",fMAT,0,0,0);
	if(flags)
	  printf("  ---FEND\n");
	else
	  printf("\n");
}



