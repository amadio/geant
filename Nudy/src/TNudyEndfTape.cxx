/*
   This is the main class to read a file in ENDF format and write a file
   in R-ENDF format

   fca 2-mai-2010
*/

#include <string.h>

#include <Riostream.h>
#include <TNudyEndfTape.h>
#include <TNudyEndfMat.h>

ClassImp(TNudyEndfTape)

//_______________________________________________________________________________
TNudyEndfTape::TNudyEndfTape() :
  fLogLev(0){
  //
  // Standard constructor
  //
  fName[0]='\0';
  fMats = new TList();
}

//_______________________________________________________________________________
TNudyEndfTape::TNudyEndfTape(const Char_t *name, UChar_t loglev) :
  fLogLev(loglev)
{
  strncpy(fName,name,80);
  std::cout << "Creating ENDF Tape:" << std::endl << name << std::endl;
  fMats= new TList();

};

//______________________________________________________________________________
TNudyEndfTape::~TNudyEndfTape()
{
  //  printf("Destroying TAPE %s\n",fName);
  fMats->Delete();
  SafeDelete(fMats);
}

//_______________________________________________________________________________
void TNudyEndfTape::DumpENDF(Int_t flags = 1)
{
	//Name of the tape
	printf("%80s\n",fName);
	//Materials
	for(Int_t i=0;i<=fMats->LastIndex();i++)
	{
		TNudyEndfMat *mat = (TNudyEndfMat*)fMats->At(i);
		mat->DumpENDF(flags);
	}
	//TEND
	printf("%66s"," 0.000000+0 0.000000+0          0          0          0          0");
	printf("%4d%2d%3d%5d",-1,0,0,0);
	if(flags)
	  printf("  ---TEND\n");
	else
	  printf("\n");
}

//_______________________________________________________________________________
void TNudyEndfTape::AddMat(TNudyEndfMat* mat) 
{
  fMats->Add(mat);
}

//_______________________________________________________________________________
TNudyEndfMat* TNudyEndfTape::GetMAT(Int_t MAT)
{
	for(Int_t i=0;i<=this->GetMats()->LastIndex();i++)
	{
		TNudyEndfMat *thisMat = (TNudyEndfMat*)this->GetMats()->At(i);
		if(thisMat->GetMAT()==MAT)
			return thisMat;
	}
	Error("TNudyEndfMat::GetMAT(Int_t)","Could not find material %d on tape",MAT);
	return NULL;

}

//_______________________________________________________________________________
TNudyEndfMat* TNudyEndfTape::GetMAT(Int_t Z, Int_t A)
{
	Int_t ZA = 1000*Z + A;
	for(Int_t i=0;i<=this->GetMats()->LastIndex();i++)
	{
		TNudyEndfMat *thisMat = (TNudyEndfMat*)this->GetMats()->At(i);
		if(thisMat->GetZA()==ZA)
			return thisMat;
	}
	return NULL;
}
TNudyEndfMat* TNudyEndfTape::GetMAT(TString name)
{
	for(Int_t i=0;i<=this->GetMats()->LastIndex();i++)
	{
		TNudyEndfMat *thisMat = (TNudyEndfMat*)this->GetMats()->At(i);
		if(name.CompareTo(thisMat->GetName()) == 0)
			return thisMat;
	}
	return NULL;
}


