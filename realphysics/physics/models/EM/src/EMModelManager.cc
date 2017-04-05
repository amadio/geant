
#include "EMModelManager.h"

#include "EMModel.h"
#include "Particle.h"

#include <iostream>


namespace geantphysics {

EMModelManager::EMModelManager() : fNumModels(0) {}

// it owns the EMModels so it will be responsible to delete them
EMModelManager::~EMModelManager() {
  for (unsigned long i=0; i<fModelList.size(); ++i) {
    delete fModelList[i];
  }
}


void EMModelManager::Initialise(std::vector<bool> &listofisactiveregions, const std::string &procname,
                                const std::vector<Particle*> &particlelist) {
  int numRegions = listofisactiveregions.size();
  for (int i=0; i<fNumModels; ++i) {
    // set the default active regions i.e. those where the process is active
    (fModelList[i]->GetListActiveRegions()).resize(numRegions);
    for (unsigned long j=0; j<listofisactiveregions.size(); ++j) {
      (fModelList[i]->GetListActiveRegions())[j] = listofisactiveregions[j];
    }
    // check the list of reagions where the user inactivated the model
    // two consistency check only with warnings:
    // 1. user tries to inactive model for region with region index that is > maximum region index.
    // 2. user tries to inactive model for region where the process itself is not active
    std::vector<int> userInact = fModelList[i]->GetListOfUserRequestedInActiveRegions();
    for (unsigned long j=0; j<userInact.size(); ++j) {
      if (userInact[j]>numRegions) {
        std::string partNames = "   \n";
        for (unsigned long p=0; p<particlelist.size(); ++p)
          partNames += particlelist[p]->GetName() + "  ";
        partNames += "\n";
        std::cerr<<" *** WARNING: EMModelManager::Initialise() \n"
                 <<"  User requested inactivation of model Name = " << fModelList[i]->GetName() << "\n"
                 <<"  For process with Name = " << procname << " that is assigned for particle(s):" << partNames
                 <<"  in region index = " << userInact[j] << " is ignored since max region index = "
                 << numRegions-1 <<" !"
                 << std::endl;
        continue;
      }
      if ((fModelList[i]->GetListActiveRegions())[userInact[j]]) {
        (fModelList[i]->GetListActiveRegions())[userInact[j]] = false;
      } else {
        std::string partNames = "   \n";
        for (unsigned long p=0; p<particlelist.size(); ++p)
          partNames += particlelist[p]->GetName() + "  ";
        partNames += "\n";
        std::cerr<<" *** WARNING: EMModelManager::Initialise() \n"
                 <<"  User requested inactivation of model Name = " << fModelList[i]->GetName() << "\n"
                 <<"  For process with Name = " << procname << " that is assigned for particle(s):" << partNames
                 <<"  in region index = " << userInact[j] << " is ignored since the process itself is inactive in"
                 <<" that region. "
                 << std::endl;
      }
    }
  }
  // set up model lists per regions
  // clear
  for (unsigned long j=0; j<fModelListsPerRegions.size(); ++j) {
    fModelListsPerRegions[j].clear();
  }
  fModelListsPerRegions.resize(numRegions);
  for (unsigned long i=0; i<fModelList.size(); ++i) {
    EMModel *emModel = fModelList[i];
    for (int j=0; j<numRegions; ++j) {
      if (emModel->IsActiveRegion(j)) {
        fModelListsPerRegions[j].push_back(emModel);
      }
    }
  }
  // initialise the models; they will be initialised only for their active regions
  for (int i=0; i<fNumModels; ++i) {
    fModelList[i]->Initialize();
  }
}


int EMModelManager::AddModel(EMModel *mod) {
    mod->SetIndex(fModelList.size());
    fModelList.push_back(mod);
    ++fNumModels;
    return fNumModels-1;
}


EMModel* EMModelManager::SelectModel(double ekin, int regionindx) {
  // get models for the given region
  const std::vector<EMModel*> models = GetModelListInRegion(regionindx);
  int numModels = models.size();
  EMModel *selectedEMModel = nullptr;
  if (numModels>0) {
    int modelIndx = 0;
    if (numModels>1) {
      modelIndx = numModels;
      do {--modelIndx;} while (modelIndx && ekin<=models[modelIndx]->GetLowEnergyUsageLimit());
      //modelIndx = numModels-1;
      //for (; modelIndx && ekin<=models[modelIndx]->GetLowEnergyUsageLimit(); --modelIndx){}
    }
    selectedEMModel = models[modelIndx];
  } else {
    std::cerr<< "  **** ERROR EMModelManager::SelectModel() \n"
             << "     No model in region so discrete interaction should not happen!"
             <<std::endl;
    exit(-1);
  }
  return selectedEMModel;
}


} // namespace geantphysics
