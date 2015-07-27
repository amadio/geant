#include "TNudyLibrary.h"
#include <TFile.h>
ClassImp(TNudyLibrary)

  //______________________________________________________________________________
TNudyLibrary::TNudyLibrary()
{
  //Create TNudy Library class 
  printf("Making %s\n",GetName());
  //Create new sublibrary Hashtable;
  fSubLib = new THashTable();
  fCurSubLib = NULL;
}
TNudyLibrary::~TNudyLibrary()
{
  printf("Deleting %s\n",GetName());
  fSubLib->Delete();
  delete fSubLib;
  fSubLib = 0;
}
//_______________________________________________________________________________
void TNudyLibrary::ListModels(){
  printf("---- Library %s \n",GetName());
  if(fCurSubLib)
    fCurSubLib->ListModels();
  else{
    TIter subLibIter(fSubLib);
    TNudySubLibrary* subLib = NULL;
    while((subLib = (TNudySubLibrary*)subLibIter.Next())){
      subLib->ListModels();
    }
  }
}
//_______________________________________________________________________________
TNudyLibrary::TNudyLibrary(const char* name,const char* title)
{
  printf("Making %s\n",GetName());
  SetName(name);
  SetTitle(title);  
  //Create new sublibrary Hashtable;
  fSubLib = new THashTable();
  fCurSubLib = NULL;
}

//_______________________________________________________________________________
void TNudyLibrary::ReadTape(TNudyEndfTape *tape){
  //Function to read and process ROOT ENDF tape
  int objcount = 0;
  //Create Material Iterator
  TIter iter(tape->GetMats());
  TNudyEndfMat *mat;
  Info("TNudyLibrary::ReadTape","Reading RENDF Tape");
  THashList* pList = (THashList*)TNudyCore::Instance()->GetParticleList()->Clone();

  //Get Iterator for list of particles
  TIter pIter(pList);
  TParticlePDG *particle;
  //Store mass of interested particle right now only neutron
  //**Make more general
  double neutron = TNudyCore::Instance()->GetParticlePDG(kNeutron)->Mass();
  //Iterate through all materials
  while((mat = (TNudyEndfMat*)iter.Next())){
    //Store relative mass to neutron of incident particle
    double awi = mat->GetAWI();
    //    printf("In Material %s\n",mat->GetName());
    //Iterate through all particles
    pIter.Reset();
    while((particle = (TParticlePDG*)pIter.Next())){
      //If particles relative masses are equal material is read into corresponding sub-library and that it is not a anti particle
      //Anti-particles have negative pdgcodes
      if(TMath::AreEqualAbs(awi,particle->Mass()/neutron,0.0000001) == kTRUE && particle->PdgCode() > 0){
        //Create a new TParticle of the right type
        //Create of new Sublibrary of TParticle part
        TNudySubLibrary* subLib = NULL;
        if(!(subLib = (TNudySubLibrary*)fSubLib->FindObject(particle->GetName()))){

          subLib = new TNudySubLibrary(particle);
          fSubLib->Add(subLib);
          printf("Making new sub library %s\n",subLib->GetName());
        }
        //If folder does not exist create it
        //	printf("Material %s\n",mat->GetName());
        //	gDirectory->pwd();
        if(!(gDirectory->FindKey(subLib->GetName()))) {
          gDirectory->mkdir(subLib->GetName());
          //Go to Sub-Library folder	
        }
        gDirectory->cd(subLib->GetName());  
        //Read and Process Material 
        subLib->ReadMat(mat);
        objcount++;
        //Return to root directory of Library
        gDirectory->cd("..");       
      } 
    }
  }
  //  pList->Delete();
  delete pList;
  printf("%d Models written in Library %s\n",objcount,GetName());
}

//_______________________________________________________________________________
TNudySubLibrary* TNudyLibrary::AddSubLib(TParticlePDG *particle){
  if(!particle) return NULL;
  TNudySubLibrary *newSubLib;
  if((newSubLib=(TNudySubLibrary*)fSubLib->FindObject(particle->GetName())))
    return newSubLib;
  newSubLib = new TNudySubLibrary(particle);
  fSubLib->Add(newSubLib);
  return SetSubLib(particle);
}

//_______________________________________________________________________________
Bool_t TNudyLibrary::IsHandled(TParticlePDG *particle, TGeoElementRN *targets, ULong_t temp){
  return (kTRUE);
}
