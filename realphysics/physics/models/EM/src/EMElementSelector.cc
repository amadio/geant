

#include "EMElementSelector.h"
// from material
#include "Types.h"


#include "EMModel.h"

#include "MaterialCuts.h"
#include "Material.h"
#include "MaterialProperties.h"
#include "Element.h"

#include <cmath>

namespace geantphysics {

EMElementSelector::EMElementSelector(EMModel *emmodel, double emin, double emax, int binsperdecade, int numelement)
: fNumEnergyBins(0), fMinEnergy(emin), fMaxEnergy(emax), fLogMinEnergy(-1.), fEnergyILDelta(-1.), fEnergyGrid(nullptr),
  fEMModel(emmodel) {
  InitializeEnergyGrid(binsperdecade);
  fProbsPerElements.resize(numelement,nullptr);
  for (int iel=0; iel<numelement; ++iel) {
    fProbsPerElements[iel] = new double[fNumEnergyBins]();
  }
}


EMElementSelector::~EMElementSelector() {
  if (fEnergyGrid) {
    delete [] fEnergyGrid;
  }
  for (unsigned long i=0; i<fProbsPerElements.size(); ++i) {
    if (fProbsPerElements[i]) {
      delete [] fProbsPerElements[i];
    }
  }
  fProbsPerElements.clear();
}


int EMElementSelector::SampleTargetElement(double ekin, double rndm) {
  //int elIndx = -1;
  unsigned long elIndx = 0;
  if (ekin>=fMinEnergy && ekin<=fMaxEnergy) {
    double logE     = std::log(ekin);
    int    lowEIndx = (int) ((logE-fLogMinEnergy)*fEnergyILDelta);
    // we might put it under verbose build since
    // protection against very small numerical uncertainties
//      if (lowEIndx>0 && ekin<fEnergyGrid[lowEIndx]) {
//        --lowEIndx;
//      } else if (ekin>fEnergyGrid[lowEIndx+1]) {
//        ++lowEIndx;
//      }
    if (lowEIndx>=fNumEnergyBins-1) --lowEIndx;
    // linear interpolation on log E scale
    double factor = (ekin-fEnergyGrid[lowEIndx])/(fEnergyGrid[lowEIndx+1]-fEnergyGrid[lowEIndx]);
    for (elIndx=0; elIndx<fProbsPerElements.size()-1; ++elIndx) {
      double val =  fProbsPerElements[elIndx][lowEIndx]
                  + (fProbsPerElements[elIndx][lowEIndx+1]-fProbsPerElements[elIndx][lowEIndx])*factor;
      if (rndm<=val) {
        break;
      }
    }
  }
  return elIndx;
}

void EMElementSelector::Build(const MaterialCuts *matcut, const Particle *part) {
  const Vector_t<Element*> elemVect       = matcut->GetMaterial()->GetElementVector();
  const double *theAtomicNumDensityVector = matcut->GetMaterial()->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElems = elemVect.size();
  for (int i=0; i<fNumEnergyBins; ++i) {
    double ekin = fEnergyGrid[i];
    double sum = 0.0;
    // fill per element cumulatives; will be normalized below
    for (int iel=0; iel<numElems; ++iel) {
      double xSecPerElem  = fEMModel->ComputeXSectionPerAtom(elemVect[iel], matcut, ekin, part);
      sum += xSecPerElem*theAtomicNumDensityVector[iel];
      fProbsPerElements[iel][i] = sum;
    }
  }
  // correct possible problems at min/max energy limits: all models might gives back zero x-section at min energies
  // if the min energy corresponds to production energy threshold
  if (fProbsPerElements[numElems-1][0]==0.0) {
    for (int i=0; i<numElems; ++i) {
      fProbsPerElements[i][0] = fProbsPerElements[i][1];
    }
  }
  if (fProbsPerElements[numElems-1][fNumEnergyBins-1]==0.0) {
    for (int i=0; i<numElems; ++i) {
      fProbsPerElements[i][fNumEnergyBins-1] = fProbsPerElements[i][fNumEnergyBins-2];
    }
  }
  // normalization i.e. final step to prepare cumulatives per elements
  for (int i=0; i<fNumEnergyBins; ++i) {
    double norm = fProbsPerElements[numElems-1][i];
    if (norm>0.0) {
      for (int iel=0; iel<numElems; ++iel) {
        fProbsPerElements[iel][i] /= norm;
      }
    }
  }
}


void EMElementSelector::InitializeEnergyGrid(int binsperdecade) {
  static const double invlog106 = 1.0/(6*std::log(10.));
  fNumEnergyBins = (int)(binsperdecade*std::log(fMaxEnergy/fMinEnergy)*invlog106);
  if (fNumEnergyBins<3) {
    fNumEnergyBins = 3;
  }
  ++fNumEnergyBins;
  if (!fEnergyGrid) {
    delete [] fEnergyGrid;
    fEnergyGrid = nullptr;
  }
  fEnergyGrid                      = new double[fNumEnergyBins]();
  fLogMinEnergy                    = std::log(fMinEnergy);
  double delta                     = std::log(fMaxEnergy/fMinEnergy)/(fNumEnergyBins-1.0);
  fEnergyILDelta                   = 1.0/delta;
  fEnergyGrid[0]                   = fMinEnergy;
  fEnergyGrid[fNumEnergyBins-1] = fMaxEnergy;
  for (int i=1; i<fNumEnergyBins-1; ++i) {
    fEnergyGrid[i] = std::exp(fLogMinEnergy+i*delta);
  }
}


} // namespace geantphysics
