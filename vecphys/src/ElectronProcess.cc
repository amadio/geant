#include "GUAliasSampler.h"
#include "GUAliasTable.h"

#include "ElectronProcess.h"

// from VecGeom
#include "materials/Material.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_CUDA_HOST
ElectronProcess::ElectronProcess(Random_t *states, int tid) : EmProcess<ElectronProcess>(states, tid)
{
  fNumberOfProcess = 2;
  fIonisation = 0;
  fBremsstrahlung = 0;

  fLogEnergyLowerBound = log(fEnergyLowerBound);
  fInverseLogEnergyBin = fNumberOfEnergyBin / (log(fEnergyUpperBound) - fLogEnergyLowerBound);

  Initialization();
}

VECCORE_CUDA_HOST
ElectronProcess::~ElectronProcess()
{
  delete fIonisation;
  delete fBremsstrahlung;

  for (int i = 0; i < fNumberOfMaterialBin; ++i) {
    free(fElectronCrossSectionData[i]);
  }
}

VECCORE_CUDA_HOST void ElectronProcess::Initialization()
{
  fIonisation = new IonisationMoller(0, -1);
  fBremsstrahlung = new BremSeltzerBerger(0, -1);

  fNumberOfMaterialBin = (vecgeom::Material::GetMaterials()).size();
  fElectronCrossSectionData =
      (ElectronCrossSectionData **)malloc(sizeof(ElectronCrossSectionData *) * fNumberOfMaterialBin);

  for (int i = 0; i < fNumberOfMaterialBin; ++i) {
    fElectronCrossSectionData[i] =
        (ElectronCrossSectionData *)malloc(sizeof(ElectronCrossSectionData) * fNumberOfEnergyBin);
  }

  // initialize table
  for (int i = 0; i < fNumberOfMaterialBin; ++i) {
    for (int j = 0; j < fNumberOfEnergyBin; ++j) {
      for (int k = 0; k < fNumberOfProcess; ++k)
        fElectronCrossSectionData[i][j].fAlias[k] = 0;
      fElectronCrossSectionData[i][j].fSigma = 0.0;
      for (int k = 0; k < fNumberOfProcess - 1; ++k)
        fElectronCrossSectionData[i][j].fWeight[k] = 0.0;
    }
  }

  BuildCrossSectionTable();
  // PrintCrossSectionTable();
}

VECCORE_CUDA_HOST void ElectronProcess::BuildCrossSectionTable()
{
  // Get the global material table
  //  const std::vector<geant::Material*> &mtable = geant::Material::GetTheMaterialTable();
  // Get the material DB
  std::vector<vecgeom::Material *> &mtable = vecgeom::Material::GetMaterials();

  double cross[3] = {
      0., 0., 0.,
  };
  double energy = 0;

  double logBinInterval = (log(fEnergyUpperBound) - fLogEnergyLowerBound) / fNumberOfEnergyBin;

  for (int i = 0; i < fNumberOfMaterialBin; ++i) {
    for (int j = 0; j < fNumberOfEnergyBin; ++j) {
      energy = exp(fLogEnergyLowerBound + logBinInterval * (j + 0.5));

      cross[0] = fIonisation->G4CrossSectionPerVolume((mtable)[i], energy);
      cross[1] = fBremsstrahlung->G4CrossSectionPerVolume((mtable)[i], energy);
      // add-msc

      // fill cross section information (total and weights)
      fElectronCrossSectionData[i][j].fSigma = cross[0] + cross[1];

      if (fElectronCrossSectionData[i][j].fSigma != 0.0) {
        fElectronCrossSectionData[i][j].fWeight[0] = cross[0] / (fElectronCrossSectionData[i][j].fSigma);
        fElectronCrossSectionData[i][j].fWeight[1] = cross[1] / (fElectronCrossSectionData[i][j].fSigma);
      }

      // alias table
      int *a = (int *)malloc(fNumberOfProcess * sizeof(int));
      double *ap = (double *)malloc(fNumberOfProcess * sizeof(double));

      double cp = 1. / fNumberOfProcess;

      // pdf for alias
      double pdf[3] = {fElectronCrossSectionData[i][j].fWeight[0], fElectronCrossSectionData[i][j].fWeight[1],
                       1.0 - fElectronCrossSectionData[i][j].fWeight[0] - fElectronCrossSectionData[i][j].fWeight[1]};

      for (int k = 0; k < fNumberOfProcess; ++k) {
        a[k] = -1;
        ap[k] = pdf[k];
      }

      // O(n) iterations
      int iter = fNumberOfProcess;

      do {
        int donor = 0;
        int recip = 0;

        // A very simple search algorithm
        for (int k = donor; k < fNumberOfProcess; ++k) {
          if (ap[k] >= cp) {
            donor = k;
            break;
          }
        }

        for (int k = recip; k < fNumberOfProcess; ++k) {
          if (ap[k] >= 0.0 && ap[k] < cp) {
            recip = k;
            break;
          }
        }

        // alias and non-alias probability
        fElectronCrossSectionData[i][j].fAlias[recip] = donor;

        // update pdf
        ap[donor] = ap[donor] - (cp - ap[recip]);
        ap[recip] = -1.0;
        --iter;

      } while (iter > 0);

      free(a);
      free(ap);
    }
  }
}

VECCORE_CUDA_HOST void ElectronProcess::PrintCrossSectionTable()
{
  printf("Number of Material Bins = %d\n", fNumberOfMaterialBin);
  printf("Number of Energy Bins = %d\n", fNumberOfEnergyBin);

  for (int i = 0; i < fNumberOfMaterialBin; ++i) {
    for (int j = 0; j < fNumberOfEnergyBin; ++j) {
      printf("[M=%d][E=%d] = [ %g %g %g %d %d %d]\n", i, j, fElectronCrossSectionData[i][j].fSigma,
             fElectronCrossSectionData[i][j].fWeight[0], fElectronCrossSectionData[i][j].fWeight[1],
             fElectronCrossSectionData[i][j].fAlias[0], fElectronCrossSectionData[i][j].fAlias[1],
             fElectronCrossSectionData[i][j].fAlias[2]);
    }
  }
}

} // end namespace impl
} // end namespace vecphys
