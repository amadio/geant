#include "TNudySubLibrary.h"
#include "TParticlePDG.h"
#include "TBtree.h"
#include "TVNudyModel.h"
#include "TNudyCore.h"
#include "TNudyEndfMat.h"
#include "TROOT.h"
#include "TFile.h"

#ifdef USE_ROOT
ClassImp(TNudySubLibrary)
#endif

//______________________________________________________________________________
TNudySubLibrary::TNudySubLibrary() {
  // Default Constructor to create a new SubLibrary
  printf("Making SubLibrary %s\n", GetName());
  // BTree to store data for fast processing
  fIndex = new TBtree();
  fBuffer = NULL;
  fProjectile = NULL;
}

//______________________________________________________________________________
TNudySubLibrary::TNudySubLibrary(TParticlePDG *projectile) {
  // Constructor to create a new SubLibrary for TParticle projectile
  SetName(projectile->GetName());
  printf("Making SubLibrary %s\n", GetName());
  SetTitle(Form("%s sub library", projectile->GetName()));
  fProjectile = projectile;
  // BTree to store data for fast processing
  fIndex = new TBtree();
  fBuffer = NULL;
}

//______________________________________________________________________________
void TNudySubLibrary::AddModel(TVNudyModel *model) {
  // Every Model should be unique
  if (!(fIndex->FindObject(model->GetName()))) {
    fIndex->Add(model);
  } else {
    //    Error("AddModel","Duplicate Model %s",model->GetTitle());
  }
}

//______________________________________________________________________________
void TNudySubLibrary::ListModels() {
  TIter modelIter(fIndex);
  TVNudyModel *model;
  printf("---- SubLibrary %s Model List\n\n", GetName());
  while ((model = (TVNudyModel *)modelIter.Next())) {
    printf(" %10s | A=%-3d Z=%-3d ISO=%-1d Temperature=%ldK  Reaction=%s | %s\n ", model->GetMaterialName(),
           model->GetA(), model->GetZ(), model->GetISO(), model->GetTemp(),
           TNudyCore::Instance()->ExpandReaction(model->GetReaction()), model->GetTitle());
  }
}

//______________________________________________________________________________
TNudySubLibrary::~TNudySubLibrary() {
  // Destructor for TNudySubLibrary
  printf("Deleting SubLibrary %s\n", GetName());
  delete fProjectile;
  fIndex->Delete();
  delete fIndex;
  if (fBuffer)
    delete fBuffer;
  fBuffer = NULL;
}

//______________________________________________________________________________
void TNudySubLibrary::ReadMat(TNudyEndfMat *material) {
  // Function to Read and Process a Material
  TList reactions;
  TIter iter(material->GetFiles());
  TNudyEndfFile *file = NULL;
  TNudyEndfSec *sec = NULL;

  // Loop through all sections and create a list of all reactions that occur
  while ((file = (TNudyEndfFile *)iter.Next())) {
    TIter secIter(file->GetSections());
    while ((sec = (TNudyEndfSec *)secIter.Next())) {
      if (!reactions.FindObject(TString(Form("%d", sec->GetMT())))) {
        reactions.Add(new TNamed(TString(Form("%d", sec->GetMT())), TString("")));
      }
    }
  }
  TIter rIter(&reactions);
  TNamed *obj = NULL;
  int ZA = material->GetZA();
  // Get Default Isotopes
  if (ZA % 1000 == 0) {
    ZA = ZA + (int)TNudyCore::Instance()->GetElementTable()->GetElement(ZA / 1000)->A();
  }
  TGeoElementRN *mat = TNudyCore::Instance()->GetMaterial(ZA * 10 + (int)material->GetLISO());
  if (!mat) {
    Warning("ReadMat", "Material %s MAT %d ISO %d does not exist. Check the ENDF file.", material->GetName(),
            material->GetMAT(), material->GetLISO());
    return;
  }
  //  printf("Material Obtained %s\n",mat->GetName());
  while ((obj = (TNamed *)rIter.Next())) {
    unsigned long temp = material->GetTEMP();
    Reaction_t reac = (Reaction_t)TString(obj->GetName()).Atoi();
    if (temp == 0) {
      Error("ReadMat", "Data for MAT=%d is not Doppler Broadened - Temperature not set\n", material->GetMAT());
    }
    if (!(gDirectory->FindKey(TNudyCore::Instance()->GetKey(mat, reac, temp)))) {
      if (mat->ENDFCode() == 120240)
        Info("ReadMat", "Reading materials %s MAT - %d %s into Sublibrary Tape %s",
             TNudyCore::Instance()->GetKey(mat, reac, temp), material->GetMAT(), material->GetName(), GetName());
      TVNudyModel *newModel = new TVNudyModel(mat, reac, temp, fProjectile, material);
      //      fIndex->Add(newModel);
      //      printf("New Model Created %s\n",newModel->GetName());
      newModel->Write();
      delete newModel;
    } else {
      Warning("ReadFile", "Model %s - %s MAT - %d %s already present in file %s ",
              TNudyCore::Instance()->GetKey(mat, reac, temp), mat->GetName(), material->GetMAT(), material->GetName(),
              gFile->GetName());
      // TVNudyModel *x = (TVNudyModel*)gDirectory->Read(TNudyCore::Instance()->GetKey(mat,reac,temp));
      // x->Print();
      // reactions.Print();
    }
  }
  reactions.Delete();
}
TVNudyModel *TNudySubLibrary::GetModel(const TGeoElementRN *mat, const Reaction_t reac, const unsigned long temp) {
  if (fBuffer)
    delete fBuffer;
  fBuffer = new TBtree();
  if (!mat)
    return NULL;
  TVNudyModel *model = (TVNudyModel *)fIndex->FindObject(TNudyCore::Instance()->GetKey(mat, reac, temp));
  if (model)
    fBuffer->Add(model);
  return model;
}

TBtree *TNudySubLibrary::GetAllModels(const TGeoElementRN *mat, const Reaction_t reac, const unsigned long temp) {
  if (fBuffer)
    delete fBuffer;
  fBuffer = new TBtree();
  TIter modelIter(fIndex);
  TVNudyModel *model;
  while ((model = (TVNudyModel *)modelIter.Next())) {
    if (TNudyCore::Instance()->IsMaterial(mat, model->GetName()) &&
        TNudyCore::Instance()->IsReaction(reac, model->GetName()) &&
        TNudyCore::Instance()->IsTemperature(temp, model->GetName())) {
      fBuffer->Add(model);
    }
  }
  return fBuffer;
}
