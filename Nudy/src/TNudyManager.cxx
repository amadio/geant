#include "TNudyManager.h"
#include <TROOT.h>

ClassImp(TNudyManager)

TNudyManager* TNudyManager::fgInstance = 0;

//______________________________________________________________________________
TNudyManager::~TNudyManager(){
  //Destructor
  fNudyDB->Delete();
  delete fNudyDB;
  fLibrary->Delete();
  delete fLibrary;
  fLibrary = NULL;
  fCurNudyDB = NULL;
  fCurLibrary = NULL;
  if(fResult) delete fResult;
  fResult = NULL;
  delete fCore;
  gROOT->GetListOfSpecials()->Remove(this);
  fgInstance = NULL;
}

//______________________________________________________________________________
TNudyManager::TNudyManager() : TNamed("NudyManager","Manager of Nudy ENDF Framework"){
  fCore = TNudyCore::Instance();
  fNudyDB = new THashTable();
  fLibrary = new THashTable();
  fResult = NULL;
  fCurNudyDB = NULL;
  fCurLibrary = NULL;
  if(fgInstance){
    Warning("TNudyManager","object already instantiated");
  } else {
    fgInstance = this;
    gROOT->GetListOfSpecials()->Add(this);
  }
}

//______________________________________________________________________________
TNudyManager* TNudyManager::Instance(){
  if(!fgInstance)
    fgInstance = new TNudyManager();
  return fgInstance;
}

//_______________________________________________________________________________
void TNudyManager::DumpTape(const char* rendf,const int debug){
	TFile *file = TFile::Open(rendf,"OLD");
	if(!file)
		Fatal("ctor","Could not open RENDF file %s",rendf);
	// Read to Tape Identifier
	TNudyEndfTape *tape = new TNudyEndfTape();
	tape->Read(file->GetListOfKeys()->First()->GetName());
	tape->DumpENDF(debug);
	file->Close();
	delete tape;
	delete file;
}

//______________________________________________________________________________
void TNudyManager::ProcessTape(const char* endf, const char* rendf){
  TNudyENDF *newRendf = new TNudyENDF(endf,rendf,"RECREATE",0);

  if(newRendf){
    newRendf->Process();
  }else{
    Fatal("ProcessTape","Error processing tape !");
  }
  delete newRendf;
}

//______________________________________________________________________________
void TNudyManager::AddEndfLibrary(const char* name, const char* endf){
  ProcessTape(endf,Form("%s.rendf.root",endf));
  AddLibrary(name,Form("%s.rendf.root",endf));
}

//______________________________________________________________________________
TNudyDB* TNudyManager::OpenDatabase(const char* name,const char* file){
  TNudyDB *newDB = new TNudyDB(name,name,file);
  fNudyDB->Add(newDB);
  fCurNudyDB = newDB;
  return fCurNudyDB;
}

//______________________________________________________________________________
TNudyDB* TNudyManager::SetDatabase(const char* name){
  TNudyDB* db = (TNudyDB*)fNudyDB->FindObject(name);
  if(db)
    fCurNudyDB = db;
  return db;
}
TNudyDB* TNudyManager::SetDatabase(const TNudyDB* name){
  TNudyDB* db = (TNudyDB*)fNudyDB->FindObject(name);
  if(db)
    fCurNudyDB = db;
  return db;
}

//______________________________________________________________________________
int TNudyManager::CloseDatabase(const char* name){
  if(name){
    TNudyDB* db = (TNudyDB*)fNudyDB->FindObject(name);
    if(db){
      delete db;
      return 1;
    }
    return 0;
  }
  else{
    if(fCurNudyDB){
      delete fCurNudyDB;
      fCurNudyDB = NULL;
      return 1;
    }
  }
  return 0;
}

//______________________________________________________________________________
TNudyLibrary* TNudyManager::LoadLibrary(const char* memLibName, const char* diskLibName,const char *sublib ,TGeoElementRN *mat ,Reaction_t reac, ULong_t temp ){
  if(!fCurNudyDB) return NULL;
  TFile *dbFile=fCurNudyDB->GetDBFile(); 
  dbFile->cd();
  gDirectory->cd("/");
  if(!dbFile->GetDirectory(diskLibName))
    return NULL;
  TNudyLibrary* newLib;
  if(!(newLib=GetLibrary(memLibName)))
    newLib = new TNudyLibrary(memLibName, Form("%s library",memLibName));
  dbFile->cd(diskLibName);
  if(sublib){
    if(!gDirectory->GetDirectory(sublib))
       return NULL;
    TNudySubLibrary *newSubLib;
    if(!(newSubLib=newLib->GetSubLib(TNudyCore::Instance()->GetParticlePDG(sublib))))
      newSubLib = newLib->AddSubLib(TNudyCore::Instance()->GetParticlePDG(sublib));
    gDirectory->cd(newSubLib->GetName());
    TList *models = gDirectory->GetListOfKeys();
    TIter iter(models);
    TObject *model;
    printf("Models");
    while((model = iter.Next())){
      if(TNudyCore::Instance()->IsMaterial(mat,model->GetName()) && TNudyCore::Instance()->IsReaction(reac,model->GetName()) && TNudyCore::Instance()->IsTemperature(temp,model->GetName())){
    	Info("Match","Model %s read",model->GetName());
        TVNudyModel* newModel = (TVNudyModel*)gDirectory->GetObjectUnchecked(model->GetName());
        newSubLib->AddModel(newModel);
      }
      else{
    	  Info("No Match","Model %s not read",model->GetName());
      }
    }
    gDirectory->cd("..");
  }
  else{
    //    printf("Get all sublibraries\n");
    TList *subdirs =  gDirectory->GetListOfKeys();
    subdirs->Print();
    TIter dirIter(subdirs);
    TObject *dir;
    while((dir = dirIter.Next())){
      //printf("Check sublibs\n");
      TNudySubLibrary *newSubLib;
      TParticlePDG* subdir = TNudyCore::Instance()->GetParticlePDG(dir->GetName());
      if(!subdir)
        Error("LoadModel","Particle %s is not handled",dir->GetName());
      if(!(newSubLib=newLib->GetSubLib(subdir))){
        //printf("Creating new Sublib\n");
        newSubLib = newLib->AddSubLib(subdir);
        //	printf("Created new Sublib\n");
      }
      //      printf("Changing dir to %d",newSubLib);
      gDirectory->cd(newSubLib->GetName());
      printf("Listof keys\n");
      TList *models = gDirectory->GetListOfKeys();
      TIter iter(models);
      models->Print();
      TObject *model;
      //      printf("Models\n");
      while((model = iter.Next())){
	TVNudyModel* newModel = (TVNudyModel*)gDirectory->GetObjectUnchecked(model->GetName());
	//	Info("LoadLibrary","Adding %s",model->GetName());
        newSubLib->AddModel(newModel);
      } 
      gDirectory->cd("..");
    }
  }
  fLibrary->Add(newLib);
  return SetLibrary(newLib);
}

//______________________________________________________________________________
void TNudyManager::ListModels(){
	if(fCurLibrary)
		fCurLibrary->ListModels();
	else{
		TIter libIter(fLibrary);
		TNudyLibrary *lib = NULL;
		while((lib = (TNudyLibrary*)libIter.Next())){
			lib->ListModels();
		}
	}

}
TVNudyModel* TNudyManager::GetModel(const int a, const int z, const int iso, const int reaction, const ULong_t temp, const char *particleName){
	if(fResult) delete fResult;
	fResult = new TBtree();
	TVNudyModel *model = NULL;
	TGeoElementRN *tar = TNudyCore::Instance()->GetMaterial(a,z,iso);
	if(particleName)
		model = fCurLibrary->GetSubLib(TNudyCore::Instance()->GetParticlePDG(particleName))->GetModel(tar,(Reaction_t)reaction,temp);
	else
		model = fCurLibrary->GetSubLib()->GetModel(tar,(Reaction_t)reaction,temp);
	fResult->Add(model);
	return model;
}

//______________________________________________________________________________
TBtree* TNudyManager::GetAllModels(const int a,const int z,const int iso, const int reaction,const ULong_t temp, const char *particleName) {
	if(fResult) delete fResult;
	fResult = new TBtree();
	TGeoElementRN *tar = TNudyCore::Instance()->GetMaterial(a,z,iso);
	TParticlePDG *proj = TNudyCore::Instance()->GetParticlePDG(particleName);
/** need to complete*/
	fResult->Add(fCurLibrary->GetSubLib(proj)->GetModel(tar,(Reaction_t)reaction,temp));
	return fResult;
}
