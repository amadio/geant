

#include "EMModel.h"

#include "EMElementSelector.h"

#include "Material.h"
#include "MaterialCuts.h"
#include "MaterialProperties.h"
#include "Element.h"

#include "Particle.h"

#include "PhysicsParameters.h"
#include "EMElementSelector.h"

#include <cmath>
#include <iostream>


namespace geantphysics {

EMModel::EMModel(const std::string& name) : fName(name), fIndex(-1), fPhysicsParameters(nullptr) {
  // it is just to make all EMModel by default active in region zero only for testing individual models with MaterialCuts
  // having 0 region index;
  // it will be cleared and reset by the EMModelManager of the process when the model is used in a process
  //fListActiveRegions.resize(1,true);
  fIsElementSelectorsPerMaterial  = false; // by default, the derived model must tell this
  fLowestSecondaryEnergy = 0.0;            // set to zero energy by defult i.e. no secondary generation limit
  fIsUseSamplingTables   = true;           // by default we use sampling tables
}

EMModel::~EMModel() {
  // delete all EMElementSelector-s if there is any
  ClearElementSelectros();
}


void  EMModel::Initialize(){
  for (unsigned long i=0; i<fListActiveRegions.size(); ++i) {
    if (fListActiveRegions[i]) {
      fPhysicsParameters = PhysicsParameters::GetPhysicsParametersForRegion(i);
      break;
    }
  }
}


double EMModel::ComputeDEDX(const MaterialCuts * /*matcut*/, double /*kinenergy*/, const Particle * /*particle*/,
                            bool /*istotal*/) { return 0.; }


// [u,v,w] is the direction given in the scattering frame i.e. taking [u1,u2,u3] direction equal to z-direction of
// this scattering frame
void EMModel::RotateToLabFrame(double &u, double &v, double &w, double u1, double u2, double u3) {
  double up = u1*u1 + u2*u2;
  if (up>0.) {
    up = std::sqrt(up);
    double px = u;
    double py = v;
    double pz = w;
    u = (u1*u3*px - u2*py)/up + u1*pz;
    v = (u2*u3*px + u1*py)/up + u2*pz;
    w =    -up*px +             u3*pz;
  } else if (u3<0.) {       // phi=0  teta=pi
    u = -u;
    w = -w;
  }
}


// only if InitialiseElementSelectors had been called from the derived class other
int EMModel::SampleTargetElementIndex(const MaterialCuts* matcut, double ekin, double rndm) {
  int elIndx = 0;
  int num    = matcut->GetMaterial()->GetNumberOfElements();
  if (num>1) {
    EMElementSelector *els = nullptr;
    // check if there is any element selector
    if (fElementSelectors.size()==0) {
      // error
      std::cerr<<"    **** ERROR EMModel::SampleTargetElementIndex() \n"
               <<"         Size of element selector = 0. Probably forgot to call \n"
               <<"         EMModel::InitialiseElementSelectors() at the end of the Initialise() \n"
               <<"         method of the derived class name = "<<GetName()<<" \n"
               <<std::endl;
    } else if (fIsElementSelectorsPerMaterial) {
      els = fElementSelectors[matcut->GetMaterial()->GetIndex()];
    } else {
      els = fElementSelectors[matcut->GetIndex()];
    }
    // just for the dev. to make sure
    if (els==nullptr) {
      std::cerr<<" ---- nullptr the EMElementSelector : emmodel"<< GetName() << "  matcut = "<< matcut->GetMaterial()->GetName()<< "  matcut indx = "<<matcut->GetIndex() << "  matindx = "<< matcut->GetMaterial()->GetIndex()<<" \n"; exit(-1);
    }
    elIndx = els->SampleTargetElement(ekin,rndm);
  }
  return elIndx;
}


void  EMModel::InitialiseElementSelectors(EMModel *emmodel, const Particle *part, bool ispermaterial) {
   std::cerr<<"  ELSelector for model = "<<GetName()<<std::endl;
   ClearElementSelectros();
   int numMatCuts       = MaterialCuts::GetTheMaterialCutsTable().size();
   int numMaterials     = Material::GetTheMaterialTable().size();
   if (ispermaterial) {
     fIsElementSelectorsPerMaterial = true;
     // the container will store nullptr at material indices that do not appear in any of the regions where this model
     // is active or the material is built up from a single element
     fElementSelectors.resize(numMaterials,nullptr);
     for (int i=0; i<numMatCuts; ++i) {
       const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[i];
       const Material     *mat    = matCut->GetMaterial();
       // if the model is active in region where this MaterialCuts is and there are more than 1 elements in the material
       // and the element selector for the corresponding material has not built yet (note: the model might active in
       // more than one region and the same material appears in more of these regions)
       if (fListActiveRegions[matCut->GetRegionIndex()] && mat->GetNumberOfElements()>1 && !fElementSelectors[mat->GetIndex()]) {
         // build an element selector for this material and insert it into fElementSelectors vector with index = to
         // the index of the material
         double emin = emmodel->MinimumPrimaryEnergy(matCut,part);
         if (emin<fLowEnergyUsageLimit) {
           emin = fLowEnergyUsageLimit;
         }
         double emax = emin*10.;
         if (emax<fHighEnergyUsageLimit) {
           emax = fHighEnergyUsageLimit;
         }
         if (emax>emin) {
           int binsPerDecade = fPhysicsParameters->GetNumLambdaTableBinsPerDecade();
           int numElems      = mat->GetNumberOfElements();
           EMElementSelector  *elSelector = new EMElementSelector(emmodel, emin, emax, binsPerDecade, numElems);
           elSelector->Build(matCut, part);
           fElementSelectors[mat->GetIndex()] = elSelector;
         }
       }
     }
   } else {
     fIsElementSelectorsPerMaterial = false;
     fElementSelectors.resize(numMatCuts,nullptr);
     for (int i=0; i<numMatCuts; ++i) {
       const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[i];
       // if the model is active in region where this MaterialCuts is and there are more than 1 elements in the material
       if (fListActiveRegions[matCut->GetRegionIndex()] && matCut->GetMaterial()->GetNumberOfElements()>1) {
         // build an element selector for this MaterialCuts and insert it into fElementSelectors vector with index = to
         // the index of the MaterialCuts
         double emin = emmodel->MinimumPrimaryEnergy(matCut,part);
         if (emin<fLowEnergyUsageLimit) {
           emin = fLowEnergyUsageLimit;
         }
         double emax = emin*10.;
         if (emax<fHighEnergyUsageLimit) {
           emax = fHighEnergyUsageLimit;
         }
         if (emax>emin) {
           int binsPerDecade = fPhysicsParameters->GetNumLambdaTableBinsPerDecade();
           int numElems      = matCut->GetMaterial()->GetNumberOfElements();
           EMElementSelector  *elSelector = new EMElementSelector(emmodel, emin, emax, binsPerDecade, numElems);
           elSelector->Build(matCut, part);
           fElementSelectors[matCut->GetIndex()] = elSelector;
         }
       }
     }
   }
   std::cerr<<"   ----> DONE!"<<std::endl;
}


void EMModel::ClearElementSelectros() {
  for (unsigned long i=0; i<fElementSelectors.size(); ++i) {
    if (fElementSelectors[i]) {
      delete fElementSelectors[i];
    }
  }
  fElementSelectors.clear();
}


} // namespace geantphysics
