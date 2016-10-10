
#include "CutConverter.h"

#include "Material.h"
#include "MaterialProperties.h"
#include "Element.h"

#include "GLIntegral.h"
#include "Spline.h"

#include <cmath>
#include <iostream>

namespace geantphysics {

CutConverter::CutConverter(int particleindx, int numebins, double mincutenergy, double maxcutenergy)
: fParticleIndx(particleindx), fMaxZet(-1), fNumEBins(numebins), fMinCutEnergy(mincutenergy), fMaxCutEnergy(maxcutenergy),
  fEnergyGrid(nullptr), fLengthVector(nullptr), fElossOrAbsXsecTable(nullptr) , fMaxLengthIndx(-1) {}


CutConverter::~CutConverter() {
  //clean up if there is anything to clean
  if (!fElossOrAbsXsecTable) {
    for (int i=0; i<fMaxZet; ++i) {
      if (!fElossOrAbsXsecTable[i]) {
        delete [] fElossOrAbsXsecTable[i];
        fElossOrAbsXsecTable[i] = nullptr;
      }
    }
    delete [] fElossOrAbsXsecTable;
    fElossOrAbsXsecTable = nullptr;
  }
  if (!fEnergyGrid) {
    delete [] fEnergyGrid;
    fEnergyGrid = nullptr;
  }
  if (!fLengthVector) {
    delete fLengthVector;
    fLengthVector = nullptr;
  }
}


// can be called after setting any new value for fNumEBins, fMinEnergy, fMaxEnergy
void CutConverter::Initialise() {
  //clean up if there is anything to clean
  if (!fElossOrAbsXsecTable) {
    for (int i=0; i<fMaxZet; ++i) {
      if (!fElossOrAbsXsecTable[i]) {
        delete [] fElossOrAbsXsecTable[i];
        fElossOrAbsXsecTable[i] = nullptr;
      }
    }
    delete [] fElossOrAbsXsecTable;
    fElossOrAbsXsecTable = nullptr;
  }
  if (!fEnergyGrid) {
    delete [] fEnergyGrid;
    fEnergyGrid = nullptr;
  }
  if (!fLengthVector) {
    delete fLengthVector;
    fLengthVector = nullptr;
  }
  // determine maximum Z number
  const std::vector<Element*> elemTable = Element::GetTheElementTable();
  int   numElements = elemTable.size();
  double maxZet = 0.0;
  for (int i=0; i<numElements; ++i) {
    double zet = elemTable[i]->GetZ();
    if (zet>maxZet)
      maxZet = zet;
  }
  fMaxZet = std::lrint(maxZet);
  // allocate space for atomic energyloss or absorption x-section data vectors
  fElossOrAbsXsecTable = new double*[fMaxZet]();
  for (int i=0; i<numElements; ++i) {
    int izet = std::lrint(elemTable[i]->GetZ());
    if (!fElossOrAbsXsecTable[izet-1])
      fElossOrAbsXsecTable[izet-1] = new double[fNumEBins]();
  }
  // allocate space for the length vector that will be used to compute the macroscopic length (range or absorption)
  fLengthVector = new double[fNumEBins]();
  // allocate space and generate the energy grid: logspacing between fMinCutEnergy-fMaxCutEnergy with fNumEBins bin
  double delta = std::log(fMaxCutEnergy/fMinCutEnergy)/(fNumEBins-1);
  double base  = std::log(fMinCutEnergy)/delta;
  fEnergyGrid  = new double[fNumEBins]();
  int i = 0;
  for (; i<fNumEBins-1; ++i)
    fEnergyGrid[i] = std::exp((base+i)*delta);
  fEnergyGrid[i] = fMaxCutEnergy;
}


// converts production cut given in length/energy to energy/legth
double CutConverter::Convert(const Material *mat, double cut, bool isfromlength) {
  fMaxLengthIndx = -1;
  BuildLengthVector(mat);
  if (!isfromlength)
   return ConvertKineticEnergyToLength(cut); // NOTE: ONLY INFORMAL i.e. NO lenght-energy-length conversion guaranted
  // go for length - to - energy conversation
  double cutInEnergy = ConvertLengthToKineticEnergy(cut);
  // low energy correction for e-/e+
  const double lowen = 30.0*geant::keV;
  const double tune  = 0.0025*geant::g/geant::cm2;
  if ((fParticleIndx==1 || fParticleIndx==2) && cutInEnergy<lowen) {
    cutInEnergy /= (1.0+(1.0-cutInEnergy/lowen)*tune/(cut*mat->GetDensity()));
  }
  // check limits
  if (cutInEnergy<fMinCutEnergy) {
    cutInEnergy = fMinCutEnergy;
  } else if (cutInEnergy>fMaxCutEnergy) {
    cutInEnergy = fMaxCutEnergy;
  }
  return cutInEnergy;
}


// Builds the energy-length function
void CutConverter::BuildLengthVector(const Material *mat) {
  const std::vector<Element*> elemVect    = mat->GetElementVector();
  const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  int   numElements = elemVect.size();
  // fill lossvect with the material dE/dx
  double *lossvect = new double[fNumEBins]();
  for (int iener=0; iener<fNumEBins; ++iener) {
    double eloss = 0.0;
    for (int ielem=0; ielem<numElements; ++ielem) {
      int izet     = std::lrint(elemVect[ielem]->GetZ());
      double *vect = fElossOrAbsXsecTable[izet-1];
      eloss       += theAtomicNumDensityVector[ielem]*vect[iener];
    }
    lossvect[iener] = eloss;
  }
  // integrate 1/dE/dx to get the range like in Geant4; //NOTE: that it is not correct in Geant4 neithr here because
  // the integral between [0,E_1] is missing!!! We keep it like this to guarantee consistency with Geant4 since it is
  // only used when production cut is given in length.
  double deltae = std::log(fMaxCutEnergy/fMinCutEnergy);
  deltae /=(fNumEBins-1);
  fLengthVector[0] = fEnergyGrid[0]/lossvect[0]*deltae;
  double sum = 0.5*fEnergyGrid[0]/lossvect[0];
  for (int i=1; i<fNumEBins; ++i) {
    double q = fEnergyGrid[i]/lossvect[i];
    sum += q;
    fLengthVector[i] = (sum-0.5*q)*deltae;
  }
//  for (int i=0; i<fNumEBins; ++i)
//    std::cerr<<fEnergyGrid[i]/geant::MeV<<" "<<fLengthVector[i]/geant::mm<<std::endl;
}


// must be called only once
void CutConverter::BuildElossOrAbsXsecTable() {
//  std::cerr<<"  Table is building : max Z now is = "<<fMaxZet<<std::endl;
  for (int i=0; i<fMaxZet; ++i) {
    if (fElossOrAbsXsecTable[i]) {
      double zet   = i+1.0;
      double *vect = fElossOrAbsXsecTable[i];
      for (int iekin=0; iekin<fNumEBins; ++iekin) {
        vect[iekin] = ComputeELossOrAbsXsecPerAtom(zet, fEnergyGrid[iekin]);
      }
    }
  }
}


double CutConverter::ConvertLengthToKineticEnergy(double cutlenght) {
/*
// This is a more compact and accurate way that results in a more exact (compared to geant4) lenght-energy-length
// conversions
std::cerr<<" --- "<<fMaxLengthIndx<<std::endl;
  // if it is photon(i.e. fMaxLengthIndx was set) and the cut length is higher than the maximum photoabsorption lenght
  // then return with the maximum production energy threshold
  if (fMaxLengthIndx>-1 && cutlenght>=fLengthVector[fMaxLengthIndx])
    return fMaxCutEnergy;
  if (cutlenght<=fLengthVector[0])
    return fMinCutEnergy;

  int numdata = fNumEBins;
  // check if it was photon and set only the monotonic lower part for the spline
  if (fMaxLengthIndx>-1)
    numdata = fMaxLengthIndx+1;
  // set a spline interpolator on length-energy
  Spline *sp = new Spline(fLengthVector, fEnergyGrid, numdata, true);
  double res = sp->GetValueAt(cutlenght);
  delete sp;
  return res;
*/
 //
 // NOTE: This is the less efficient and accurate way a la Geant4. We keep this to guarantee consistency with Geant4
 //       since it is used only if the production cut is given in length.
 //
 // THIS the Geant4 like way to do it i.e. using bisection.
  //  find max. range and the corresponding energy (rmax,Tmax)
  const double eps =0.01;
  double rmax      = -1.0;

  double T1 = fMinCutEnergy;
  double r1 = fLengthVector[0];
  double T2 = fMaxCutEnergy;

  // check theCutInLength < r1
  if (cutlenght<=r1)
    return T1;

  int i = 0;
  for (; i<fNumEBins; ++i) {
    double T = fEnergyGrid[i];
    double r = fLengthVector[i];
    if (r>rmax)
      rmax=r;
    if (r<cutlenght) {
      T1 = T;
      r1 = r;
    } else if (r>cutlenght) {
      T2 = T;
      break;
    }
  }
  // check if cut in length is smaller than range/radioationlength max
  if (cutlenght>=rmax)
    return  fMaxCutEnergy;
  // convert range to energy
  Spline *sp = new Spline(fEnergyGrid, fLengthVector, fNumEBins, true);
  double T3 = std::sqrt(T1*T2);
  double r3 = sp->GetValueAt(T3);
  const int maxloop = 1000;
  for (int itr=0; itr<maxloop; ++itr) {
    if (std::abs(1.0-r3/cutlenght)<eps)
      break;
    if (cutlenght<=r3) {
      T2 = T3;
    } else {
      T1 = T3;
    }
    T3 = std::sqrt(T1*T2);
    r3 = sp->GetValueAt(T3);
  }
  delete sp;
  return T3;
}


double CutConverter::ConvertKineticEnergyToLength(double cutenergy) {
//  if (fMaxLengthIndx>-1 && cutenergy>=fEnergyGrid[fMaxLengthIndx])
//    return 1.01*fLengthVector[fMaxLengthIndx];
  if (cutenergy<=fMinCutEnergy)
    return fLengthVector[0];
  if (cutenergy>=fMaxCutEnergy)
    return fLengthVector[fNumEBins-1];

  int istart  = 0;
  int numdata = fNumEBins;
  // for photons
  if (fMaxLengthIndx>-1) {
    if (cutenergy<fEnergyGrid[fMaxLengthIndx]) {
      numdata = fMaxLengthIndx+1;
    } else {
      istart  = fMaxLengthIndx;
      numdata = fNumEBins-fMaxLengthIndx;
    }
  }

  Spline *sp = new Spline(&(fEnergyGrid[istart]), &(fLengthVector[istart]), numdata, true);
  double cutInLenght = sp->GetValueAt(cutenergy);
  delete sp;
  return cutInLenght;
}

}  // namespace geantphysics
