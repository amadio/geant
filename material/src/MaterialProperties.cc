
#include "MaterialProperties.h"

#include "Types.h"

#include "PhysicalConstants.h"

#include "Material.h"
#include "Element.h"
#include "NISTMaterialData.h"
#include "DensityEffectData.h"

#include <cmath>
#include <limits>
#include <iomanip>

namespace geantphysics {
inline namespace GEANT_IMPL_NAMESPACE {


//
// CTR
MaterialProperties::MaterialProperties(Material *mat) : fMaterial(mat) {
  InitialiseMembers();
  ComputeBasicMaterialParameters();
  ComputeIonizationParameters();
  ComputeDensityEffectParameters();
  ComputeRadiationLength();
}


MaterialProperties::~MaterialProperties() {
  if (fNumOfAtomsPerVolVect)
    delete [] fNumOfAtomsPerVolVect;
}


void MaterialProperties::InitialiseMembers() {
  fNumOfAtomsPerVolVect      = nullptr;

  fTotalNumOfAtomsPerVol     = 0.0;
  fTotalNumOfElectronsPerVol = 0.0;
  fZeff                      = 0.0;
  fMeanExcitationEnergy      = 0.0;

  fPlasmaEnergy                   = 0.0;
  fDensityEffectParameterC        = 0.0;
  fDensityEffectParameterX0       = 0.0;
  fDensityEffectParameterX1       = 0.0;
  fDensityEffectParameterA        = 0.0;
  fDensityEffectParameterM        = 0.0;
  fDensityEffectParameterDelta0   = 0.0;

  fRadiationLength           = 0.0;


}


void MaterialProperties::ComputeBasicMaterialParameters() {
  using geant::kAvogadro;
  // get number of elements in the material and the density
  int     numElems = fMaterial->GetNumberOfElements();
  double  density  = fMaterial->GetDensity();

  // initialise corresponding members
  fTotalNumOfAtomsPerVol     = 0.0;
  fTotalNumOfElectronsPerVol = 0.0;
  fZeff                      = 0.0;
  if (fNumOfAtomsPerVolVect)
    delete fNumOfAtomsPerVolVect;
  fNumOfAtomsPerVolVect = new double[numElems];

  // get some element composition realted data from the material
  const Vector_t<Element*> elemVector = fMaterial->GetElementVector();
  const double *massFractionVect      = fMaterial->GetMassFractionVector();
  for (int i=0; i<numElems; ++i) {
    double zeff = elemVector[i]->GetZ();
    double aeff = elemVector[i]->GetA();  // atomic mass in internal [weight/amount of substance]
    fNumOfAtomsPerVolVect[i]    = kAvogadro*density*massFractionVect[i]/aeff;
    fTotalNumOfAtomsPerVol     += fNumOfAtomsPerVolVect[i];
    fTotalNumOfElectronsPerVol += fNumOfAtomsPerVolVect[i]*zeff;
  }
  fZeff = fTotalNumOfElectronsPerVol/fTotalNumOfAtomsPerVol;
}


void MaterialProperties::ComputeIonizationParameters() {
  // init the mean excitation energy
  fMeanExcitationEnergy = 0.0;
  // try to get the mean exc. energy for the material from the NISTMaterialData
  // NOTE: it will work only if we can identify this material in the NISTMaterialData
  // by its name
  int indx = NISTMaterialData::Instance().FindNISTMaterialDataIndex(fMaterial->GetName());
  if (indx>-1) {  // material could be identified by its name in the NISTMaterialData database
    fMeanExcitationEnergy = NISTMaterialData::Instance().GetMeanExcitationEnergy(indx);
  }
  // if material could not be found in the NISTMaterialData database or its
  // mean excitation energy was not defined then compute it based on the elemental
  // mean excitation energies.
  if (fMeanExcitationEnergy<=0.0) {
    fMeanExcitationEnergy = 0.0;
    // get the max Z(atomic number) that we have elemental data for in the NISTMaterialData database
    int maxZ     = NISTMaterialData::Instance().GetNumberOfElementalNISTMaterialData();
    // gent number of elements, number of atoms per unit volume vector
    int numElems = fMaterial->GetNumberOfElements();
    const Vector_t<Element*> elemVect = fMaterial->GetElementVector();
    for (int i=0; i<numElems; ++i) {
      double zeff = elemVect[i]->GetZ();
      int    z    = std::lrint(zeff);
      // get mean excitation energy for the given element
      if (z>maxZ)
        z = maxZ;
      double meanExcEnerZ = NISTMaterialData::Instance().GetMeanExcitationEnergy(z-1);
      fMeanExcitationEnergy += fNumOfAtomsPerVolVect[i]*zeff*std::log(meanExcEnerZ);
    }
    fMeanExcitationEnergy = std::exp(fMeanExcitationEnergy/fTotalNumOfElectronsPerVol);
  }
}


void MaterialProperties::ComputeDensityEffectParameters() {
  using geant::kPi;
  using geant::kClassicElectronRadius;
  using geant::kHBarPlanckCLightSquare;
  using geant::eV;
  using geant::kSTPPressure;
  using geant::kSTPTemperature;

  const double twolog10 = 2.0*std::log(10.0);
  // density and/or non-tandard gas corrector factor
  double corrFactor = 0.0;

  MaterialState matState = fMaterial->GetMaterialState();

  // get the material name and check if we take density effect parameters from the DensityEffectData
  // we take the parameters only if the material is a NIST material and we can find data in the DB
  int indx = DensityEffectData::Instance().FindDensityEffectDataIndex(fMaterial->GetName());

  // here Geant4 checks if 1. it is an elemental material i.e. has only one elemnt
  // and if yes it will take density effect data for that NIST_MAT_"Symbol" if the
  // material state does atch as well
  // Moreover it will correct the parameters nased on the denisties.
  // - Imean should be the same as in the NIST_MAT_! (by default it is ok but can be set to be different!)
  // - electrn denisty should be the same i.e. molar mass and density should be the same
  //   (density is corrected here but with the asumption that the molar mass is the same!)
  // Nontheleass let's do it like in G4 to keep consistency
  int numElems = fMaterial->GetNumberOfElements();
  if (indx<0 && numElems==1) {
    int z0 = std::lrint(fMaterial->GetElementVector()[0]->GetZ());
    indx = DensityEffectData::Instance().GetElementalIndex(z0, fMaterial->GetMaterialState()); // elements trat with index=1 in DensityeffectData
    //indx will be -1 if there is no macth
    // set corrector factor if the material was found with the correct state
    if (indx>0) { // simple NIST material DenistyEffectData was found with the correct state
                   // we have the corresponding NISTMaterialData as well so get the denisty from there
      // get the density from the orgiginal NIST DB
      double orgDensity = NISTMaterialData::Instance().GetDensity(z0-1);
      corrFactor = std::log(fMaterial->GetDensity()/orgDensity);
      // the found(below) density effect parameters will be accepted only if
      // the correction is less then 1.0 otherwise the parameters might be very different
      // so cancel the index and don't use the parameters
      if (std::abs(corrFactor)>1.0)
        indx = -1;
    }
  }

  // check if material was found in the DensityEffectData database (IMPORTANT:: Mean excitation energy
  // should be the NIST one i.e. material should be the NIST one !)
  if (indx>=0) {
    fPlasmaEnergy                  = DensityEffectData::Instance().GetPlasmaEnergy(indx);
    // we always recompute the Cbar parameter
    // fDensityEffectParameterC       = 2.0*std::log(fMeanExcitationEnergy/fPlasmaEnergy) + 1;
    fDensityEffectParameterC       = DensityEffectData::Instance().GetParameterC(indx);
    fDensityEffectParameterX0      = DensityEffectData::Instance().GetParameterX0(indx);
    fDensityEffectParameterX1      = DensityEffectData::Instance().GetParameterX1(indx);
    fDensityEffectParameterA       = DensityEffectData::Instance().GetParameterA(indx);
    fDensityEffectParameterM       = DensityEffectData::Instance().GetParameterM(indx);
    fDensityEffectParameterDelta0  = DensityEffectData::Instance().GetParameterDelta0(indx);

    // apply density correction: note, that if it was a real nist material then corrFactor = 0 at this point
    // and if it was not a nist material but single element material (case 2. that indx>0)
    fDensityEffectParameterC  -= corrFactor; // cbar_\eta = cbar - \ln(10)log10(\eta)
    fDensityEffectParameterX0 -= corrFactor/twolog10;  // X_{1\eta} = X_1 - 0.5log10(\eta)
    fDensityEffectParameterX1 -= corrFactor/twolog10;  // X_{0\eta} = X_0 - 0.5log10(\eta)
  } else {
    // we compute the density effect parameters based on:
    // Sternheimer, R. M., and R. F. Peierls. "General expression for the density
    // effect for the ionization loss of charged particles." Physical Review B 3.11 (1971): 3681.

    // fDensityEffectParameterDelta0 i.e. \f$ \delta(X=X_0) = 0 \f$ will remain
    // zero in this case because Sternheimer
    // used \f$ \delta(X=X_0) = 0 \f$ not only in the case of insulators and gases but
    // also in the case of conductors. According to the equation
    // \f$ \delta(X) = 2\ln(10)X + a(X_1-X_0)^m + C\f$ if \f$ \delta(X=X_0) = 0 \f$
    // then \f$ 0=\delta(X=X_0) = 2\ln(10)X_0 + a(X_1-X_0)^m + C \f$ that gives
    // \f$ a = (-C - 2\ln(10)X_0)/(X_1-X_0)^m \f$ always if \f$ \delta(X=X_0) = 0 \f$
    // THIS is also true for the case if we take the parameters from the DensityEffectData
    // and the taken fDensityEffectParameterDelta0 = 0

    // the plasma energy of the material considered as a free electron gas
    constexpr double dum0 = 4.0*kPi*kClassicElectronRadius*kHBarPlanckCLightSquare;
    fPlasmaEnergy  = std::sqrt(dum0*fTotalNumOfElectronsPerVol); // in internal [energy] units
    // paramater \f$ -C = 2\ln{I/(h\nu_p)} +1 \f$
    fDensityEffectParameterC = 2.0*std::log(fMeanExcitationEnergy/fPlasmaEnergy) + 1.0;
    fDensityEffectParameterM = 3.0;
    // other parameters depends on the state of the material
             // for solids and liquids
    if (matState==MaterialState::kStateSolid || matState==MaterialState::kStateLiquid) {
      if (fMeanExcitationEnergy<100.0*eV) { // Imean<100 eV
        fDensityEffectParameterX1 = 2.0;
        fDensityEffectParameterX0 = 0.2;
        if (fDensityEffectParameterC>=3.681)
          fDensityEffectParameterX0 = 0.326*fDensityEffectParameterC - 1.0;
      } else {                              // Imean>=100 eV
        fDensityEffectParameterX1 = 3.0;
        fDensityEffectParameterX0 = 0.2;
        if (fDensityEffectParameterC>=5.215)
          fDensityEffectParameterX0 = 0.326*fDensityEffectParameterC - 1.5;
      }
      // Liquid Hydrogen is special [Table 9 in Seltzer Berger : Evaluation of the collision
      // stopping power of elements and compounds for electrons and positrons(1981)]
      int z0 = std::lrint(fMaterial->GetElementVector()[0]->GetZ());
      if (1==fMaterial->GetNumberOfElements() && 1==z0) {
        fDensityEffectParameterX0 = 0.425;
        fDensityEffectParameterX1 = 2.0;
        fDensityEffectParameterM  = 5.949;
      }
      // correction should be applied here as well if denisty is not equal to the
      // the "normal" density i.e. to density that was used to derive the parameters
      // but this is exactly what we don't know here
    } else if (matState==MaterialState::kStateGas) { // for gas
      if (fDensityEffectParameterC<10.0) {
        fDensityEffectParameterX0 = 1.6;
        fDensityEffectParameterX1 = 4.0;
      } else if (fDensityEffectParameterC<10.5) {
        fDensityEffectParameterX0 = 1.7;
        fDensityEffectParameterX1 = 4.0;
      } else if (fDensityEffectParameterC<11.0) {
        fDensityEffectParameterX0 = 1.8;
        fDensityEffectParameterX1 = 4.0;
      } else if (fDensityEffectParameterC<11.5) {
        fDensityEffectParameterX0 = 1.9;
        fDensityEffectParameterX1 = 4.0;
      } else if (fDensityEffectParameterC<12.25) {
        fDensityEffectParameterX0 = 2.0;
        fDensityEffectParameterX1 = 4.0;
      } else if (fDensityEffectParameterC<13.804) {
        fDensityEffectParameterX0 = 2.0;
        fDensityEffectParameterX1 = 5.0;
      } else {
        fDensityEffectParameterX0 = 0.326*fDensityEffectParameterC - 2.5;
        fDensityEffectParameterX1 = 5.0;
      }
      // Hydrogen gas is special [Table 9 in Seltzer Berger : Evaluation of the collision
      // stopping power of elements and compounds for electrons and positrons(1981)]
      int z0 = std::lrint(fMaterial->GetElementVector()[0]->GetZ());
      if (1==fMaterial->GetNumberOfElements() && 1==z0) {
        fDensityEffectParameterX0 = 1.837;
        fDensityEffectParameterX1 = 3.0;
        fDensityEffectParameterM  = 4.754;
      }
      // Helium gas is special [Table 9 in Seltzer Berger : Evaluation of the collision
      // stopping power of elements and compounds for electrons and positrons(1981)]
      if (1==fMaterial->GetNumberOfElements() && 2==z0) {
        fDensityEffectParameterX0 = 2.191;
        fDensityEffectParameterX1 = 3.0;
        fDensityEffectParameterM  = 3.297;
      }
    }  // end gas
  }

  // correction for gas at non standard conditions (parameters were derived under
  // 0 Celsius and 1 Atmosphere that we call here STP conditions; the correction
  // factor is detemined by current_gas_density/what_would_be_this_gas_density_at_stp)
  if (matState==MaterialState::kStateGas) {
    //double curDensity     = fMaterial->GetDensity();
    double curTemperature = fMaterial->GetTemperature();
    double curPressure    = fMaterial->GetPressure();
    //double stpDenisty     = curDensity*curTemperature*kSTPPressure/(curPressure*kSTPTemperature);
    corrFactor             = std::log(curPressure*kSTPTemperature/(curTemperature*kSTPPressure));

    fDensityEffectParameterC  -= corrFactor; // cbar_\eta = cbar - \ln(10)log10(\eta)
    fDensityEffectParameterX0 -= corrFactor/twolog10;  // X_{1\eta} = X_1 - 0.5log10(\eta)
    fDensityEffectParameterX1 -= corrFactor/twolog10;  // X_{0\eta} = X_0 - 0.5log10(\eta)
  }

  // set parameter a (works both case i.e. if we compute the parameters or take them
  // because if we take them and delta(X_0)==0 we just recompute parameterA here)
  // if \f$\delta(X=X_0)=0\f$ then parameter 'a' must staisfy the \f$ 0 = 2\ln(10)X + a(X_1-X_0)^m + C\f$ at X=X_0
  if (fDensityEffectParameterDelta0==0.0) {
    double dum0 = std::pow(fDensityEffectParameterX1-fDensityEffectParameterX0,fDensityEffectParameterM);
    fDensityEffectParameterA = (fDensityEffectParameterC-twolog10*fDensityEffectParameterX0)/dum0;
  }
}


double MaterialProperties::GetDensityEffectFunctionValue(const double atx) {
  const double twolog10 = 2.0*std::log(10.0);
  double delta = 0.0;

  if (atx<fDensityEffectParameterX0) { // below lower limit i.e. either zero or for conductors might be a small constant
    // conductor and delta converge to non-zero at low momentum/mass values
    if (fDensityEffectParameterDelta0>0.0) // delta(X) = delta(X_0) 10^{2(X-X_0)} = delta(X_0) exp[2ln(10)(X-X_0)]
      delta = fDensityEffectParameterDelta0 * std::exp(twolog10*(atx-fDensityEffectParameterX0));
  } else if (atx>=fDensityEffectParameterX1) {
    delta = twolog10*atx - fDensityEffectParameterC;
  } else {
    delta = twolog10*atx + fDensityEffectParameterA*std::exp(fDensityEffectParameterM*std::log(fDensityEffectParameterX1-atx))
           - fDensityEffectParameterC;
  }
  return delta;
}


// based on Tsai complete screening
void MaterialProperties::ComputeRadiationLength() {
  // use these elastic and inelatic form factors for light elements instead of TFM
  // under the complete screening approximation
  // Tsai Table.B2.
  constexpr double FelLowZet  [] = {0.0, 5.310, 4.790, 4.740, 4.710, 4.680, 4.620, 4.570};
  constexpr double FinelLowZet[] = {0.0, 6.144, 5.621, 5.805, 5.924, 6.012, 5.891, 5.788};
  // up the constant factor
  constexpr double factor = 4.0*geant::kFineStructConst*geant::kClassicElectronRadius*geant::kClassicElectronRadius;
  // constant factors for L_el and L_inel under TFM and complete screening i.e. Tsai
  const double factorLel   = std::log(184.1499);
  const double factorLinel = std::log(1193.923);

  // the element composition data of the material
  const Vector_t<Element*> theElements = fMaterial->GetElementVector();
  int   numElems = theElements.size();

  double invRadLength = 0.0;
  for (int i=0; i<numElems; ++i) {
    // Coulomb correction from Davis,Bethe,Maximom PRL 1954 Eqs.(36-38)
    double zet = theElements[i]->GetZ();
    int izet   = std::lrint(zet);
    double mu  = zet*geant::kFineStructConst;
    double mu2 = mu*mu;
    double mu4 = mu2*mu2;
    double mu6 = mu2*mu4;
    double coulombCorrection = mu2*(1.0/(1.0+mu2) + 0.20206 - 0.0369*mu2 + 0.0083*mu4 - 0.002*mu6);
    // compute radlength
    double logZper3 = std::log(zet)/3.0;
    double dum0 =0.0;
    if (zet<=3) { // it should be < 5 in order to be consistent with the relativistic brems. models
      dum0 = zet*zet*(FelLowZet[izet]-coulombCorrection)+zet*FinelLowZet[izet];//+1/18.*zet*(zet+1.0);
    } else {
      dum0 = zet*zet*(factorLel-logZper3-coulombCorrection)+zet*(factorLinel-2.0*logZper3);//+1/18.*zet*(zet+1.0);
    }
    invRadLength += factor*fNumOfAtomsPerVolVect[i]*dum0;
  }
  fRadiationLength = (invRadLength <= 0.0 ? std::numeric_limits<double>::max() : 1.0/invRadLength);
}


//
// Printouts
std::ostream& operator<<(std::ostream& flux, const MaterialProperties* matprop) {
  using geant::eV;
  using geant::cm;
  using geant::g;
  using geant::mole;

  std::ios::fmtflags mode = flux.flags();
  flux.setf(std::ios::fixed,std::ios::floatfield);
  long prec = flux.precision(3);

  flux
    << "   Material properties for material : " << matprop->fMaterial->GetName() << "\n"
    << "     ---> Radiation length          = "    << std::setw(8) << std::setprecision(6) <<  matprop->fRadiationLength/cm << " [cm]\n"
    << "     ---> I(mean excitation neergy) = "    << std::setw(8) << std::setprecision(6) <<  matprop->fMeanExcitationEnergy/eV << " [eV]\n"
    << "     ---> Plasma energy             = "    << std::setw(8) << std::setprecision(6) <<  matprop->fPlasmaEnergy/eV << " [eV]\n"
    << "     ---> Density effect parameters : \n"
    << "            Parameter     Cbar  = "    << std::setw(8) << std::setprecision(6) <<  matprop->fDensityEffectParameterC << "\n"
    << "            Parameter       X0  = "    << std::setw(8) << std::setprecision(6) <<  matprop->fDensityEffectParameterX0 << "\n"
    << "            Parameter       X1  = "    << std::setw(8) << std::setprecision(6) <<  matprop->fDensityEffectParameterX1 << "\n"
    << "            Parameter        a  = "    << std::setw(8) << std::setprecision(6) <<  matprop->fDensityEffectParameterA << "\n"
    << "            Parameter        m  = "    << std::setw(8) << std::setprecision(6) <<  matprop->fDensityEffectParameterM << "\n"
    << "            Parameter  delta_0  = "    << std::setw(8) << std::setprecision(6) <<  matprop->fDensityEffectParameterDelta0 << "\n";

  flux.precision(prec);
  flux.setf(mode,std::ios::floatfield);
  return flux;
}


 std::ostream& operator<<(std::ostream& flux, const MaterialProperties& matprop) {
  flux << &matprop;
  return flux;
}

}  // inline namespace GEANT_IMPL_NAMESPACE {
}  // namespace geantphysics
