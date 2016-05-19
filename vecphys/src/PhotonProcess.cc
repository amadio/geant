#include "GUAliasSampler.h"
#include "GUAliasTable.h"

#include "PhotonProcess.h"
#include "GUG4TypeDef.h"

//from VecGeom
#include "materials/Material.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_CUDA_HOST
PhotonProcess::PhotonProcess(Random_t* states, int tid) 
  : EmProcess<PhotonProcess>(states,tid)
{
  fNumberOfProcess = 3;
  fCompton = 0;
  fConversion = 0;
  fPhotoElectron = 0;
  fLogEnergyLowerBound = log(fEnergyLowerBound);
  fInverseLogEnergyBin = fNumberOfEnergyBin/(log(fEnergyUpperBound)-fLogEnergyLowerBound );

  Initialization();
}

VECCORE_CUDA_HOST
PhotonProcess::~PhotonProcess() 
{
  delete fCompton;
  delete fConversion;
  delete fPhotoElectron;

  for (int i = 0 ; i < fNumberOfMaterialBin ; ++i) {
    free(fPhotonCrossSectionData[i]);
  }
}

VECCORE_CUDA_HOST void 
PhotonProcess::Initialization()
{

  fCompton = new ComptonKleinNishina(0,-1);
  fConversion = new ConversionBetheHeitler(0,-1);
  fPhotoElectron = new PhotoElectronSauterGavrila(0,-1);

  fNumberOfMaterialBin = (vecgeom::Material::GetMaterials()).size();

  fPhotonCrossSectionData = (PhotonCrossSectionData**) malloc(sizeof(PhotonCrossSectionData*)*fNumberOfMaterialBin);
  for (int i = 0 ; i < fNumberOfMaterialBin ; ++i) {
    fPhotonCrossSectionData[i] = (PhotonCrossSectionData*) malloc(sizeof(PhotonCrossSectionData)*fNumberOfEnergyBin);
  }

  //initialize table
  for (int i = 0 ; i < fNumberOfMaterialBin ; ++i) {
    for (int j = 0 ; j < fNumberOfEnergyBin ; ++j) {
      for (int k = 0 ; k < 3 ; ++k) fPhotonCrossSectionData[i][j].alias[k] = 0;
      fPhotonCrossSectionData[i][j].sigma = 0.0;
      for (int k = 0 ; k < 2 ; ++k) fPhotonCrossSectionData[i][j].weight[k] = 0.0;
    }
  }
  
  BuildCrossSectionTable();

  // PrintCrossSectionTable();
  
}

VECCORE_CUDA_HOST void 
PhotonProcess::BuildCrossSectionTable()
{
  //Get the material DB
  std::vector<vecgeom::Material*>& mtable = vecgeom::Material::GetMaterials();

  double cross[3] = {0.,0.,0.,};
  double energy = 0;

  double logBinInterval = (log(fEnergyUpperBound)-fLogEnergyLowerBound)/fNumberOfEnergyBin;

  for (int i = 0 ; i < fNumberOfMaterialBin ; ++i) {
    for (int j = 0 ; j < fNumberOfEnergyBin ; ++j) {
      energy = exp(fLogEnergyLowerBound + logBinInterval*(j+0.5));

      cross[0] = fCompton->G4CrossSectionPerVolume((mtable)[i],energy);
      cross[1] = fConversion->G4CrossSectionPerVolume((mtable)[i],energy);
      cross[2] = fPhotoElectron->G4CrossSectionPerVolume((mtable)[i],energy);

      //fill cross section information (total and weights)
      fPhotonCrossSectionData[i][j].sigma = cross[0] + cross[1] + cross[2];

      if(fPhotonCrossSectionData[i][j].sigma !=0.0) {
        fPhotonCrossSectionData[i][j].weight[0] = cross[0]/(fPhotonCrossSectionData[i][j].sigma);
        fPhotonCrossSectionData[i][j].weight[1] = cross[1]/(fPhotonCrossSectionData[i][j].sigma);
      }

      //alias table
      int *a     = (int*)   malloc(fNumberOfProcess*sizeof(int)); 
      double *ap = (double*)malloc(fNumberOfProcess*sizeof(double)); 

      const double cp = 1./fNumberOfProcess;

      //copy and initialize
      //      double pdf[fNumberOfProcess] = {fPhotonCrossSectionData[i][j].weight[0],
      double pdf[3] = {fPhotonCrossSectionData[i][j].weight[0],
                       fPhotonCrossSectionData[i][j].weight[1],
                       1.0-fPhotonCrossSectionData[i][j].weight[0]-fPhotonCrossSectionData[i][j].weight[1]};

      for(int k = 0; k < fNumberOfProcess ; ++k) {
        a[k] = -1;
        ap[k] = pdf[k];
      }

      //O(n) iterations
      int iter = fNumberOfProcess;
  
      do {
        int donor = 0;
        int recip = 0;

        // A very simple search algorithm
        for(int k = donor; k < fNumberOfProcess ; ++k) {
           if(ap[k] >= cp) {
              donor = k;
              break;
           }
        }

        for(int k = recip; k < fNumberOfProcess ; ++k) {
           if(ap[k] >= 0.0 && ap[k] < cp) {
              recip = k;
              break;
           }
        }

        //alias and non-alias probability
        fPhotonCrossSectionData[i][j].alias[recip] = donor;

        //update pdf 
        ap[donor] = ap[donor] - (cp-ap[recip]);
        ap[recip] = -1.0;
        --iter;

      }
      while (iter > 0);

      free(a);
      free(ap);
    }
  }
}

VECCORE_CUDA_HOST void 
PhotonProcess::PrintCrossSectionTable()
{
  printf("Number of Material Bins = %d\n",fNumberOfMaterialBin);
  printf("Number of Energy Bins = %d\n",fNumberOfEnergyBin); 

  for (int i = 0 ; i < fNumberOfMaterialBin ; ++i) {
    for (int j = 0 ; j < fNumberOfEnergyBin ; ++j) {
      printf("[M=%d][E=%d] = [ %g %g %g %d %d %d ]\n",i,j,
             fPhotonCrossSectionData[i][j].sigma,
	     fPhotonCrossSectionData[i][j].weight[0],
	     fPhotonCrossSectionData[i][j].weight[1],
	     fPhotonCrossSectionData[i][j].alias[0],
	     fPhotonCrossSectionData[i][j].alias[1],
	     fPhotonCrossSectionData[i][j].alias[2]
      );
    }
  }
}

} // end namespace impl
} // end namespace vecphys
