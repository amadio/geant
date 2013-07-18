#include <TFinState.h>

Int_t TFinState::fVerbose=0;

ClassImp(TFinState)

//_________________________________________________________________________
TFinState::TFinState():
   fNFstates(0),
   fNsecs(0),
   fPID(0),
   fNpart(0),
   fWeight(0),
   fKerma(0),
   fMom(0)
{
}

//_________________________________________________________________________
TFinState::TFinState(Int_t nfstates, const Float_t weight[], const Float_t kerma[], 
		     const Int_t npart[], const Float_t (*mom)[3], const Int_t pid[]):
   fNFstates(nfstates),
   fNsecs(0),
   fPID(0),
   fNpart(new Int_t[fNFstates]),
   fWeight(new Float_t[fNFstates]),
   fKerma(new Float_t[fNFstates]),
   fMom(0)
{
   memcpy(fNpart,npart,fNFstates*sizeof(Int_t));
   memcpy(fWeight,weight,fNFstates*sizeof(Float_t));
   memcpy(fKerma,kerma,fNFstates*sizeof(Float_t));

   fNsecs = 0;
   for(Int_t j=0; j<fNFstates; ++j) fNsecs+=fNpart[j];

   fMom = new Float_t[fNsecs][3];
   memcpy(fMom,mom,3*fNsecs*sizeof(Float_t));
   fPID = new Int_t[fNsecs];
   memcpy(fPID,pid,fNsecs*sizeof(Int_t));
}

//_________________________________________________________________________
TFinState::~TFinState()
{
   delete [] fMom;
   delete [] fPID;
}

//_________________________________________________________________________
TFinState& TFinState::operator=(const TFinState& right)
{
   return *this;
}

//_________________________________________________________________________
Bool_t TFinState::SetFinState(Int_t nfstates, const Float_t weight[],
			      const Float_t kerma[], const Int_t npart[],
			      const Float_t (*mom)[3], const Int_t pid[])
{
   fNFstates = nfstates;

   delete [] fNpart;
   fNpart = new Int_t[fNFstates];
   memcpy(fNpart,npart,fNFstates*sizeof(Int_t));
   delete [] fWeight;
   fWeight = new Float_t[fNFstates];
   memcpy(fWeight,weight,fNFstates*sizeof(Float_t));
   delete [] fKerma;
   fKerma = new Float_t[fNFstates];
   memcpy(fKerma,kerma,fNFstates*sizeof(Float_t));
   
   fNsecs = 0;
   for(Int_t j=0; j<fNFstates; ++j) fNsecs+=fNpart[j];

   delete [] fMom;
   fMom = new Float_t[fNsecs][3];
   memcpy(fMom,mom,3*fNsecs*sizeof(Float_t));
   delete [] fPID;
   fPID = new Int_t[fNsecs];
   memcpy(fPID,pid,fNsecs*sizeof(Int_t));
   return kTRUE;
}

