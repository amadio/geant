
#include "GSPWACorrections.h"

#include "PhysicalConstants.h"

#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

#include "Material.h"
#include "MaterialProperties.h"
#include "Element.h"
#include "MaterialCuts.h"

#include <iostream>
#include <fstream>
#include <cmath>


namespace geantphysics {

const std::string GSPWACorrections::gElemSymbols[] = {"H","He","Li","Be","B" ,
 "C" ,"N" ,"O" ,"F" ,"Ne","Na","Mg","Al","Si","P" , "S","Cl","Ar","K" ,"Ca","Sc",
 "Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb",
 "Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,
 "Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
 "Yb","Lu","Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
 "Rn","Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf"};

GSPWACorrections::GSPWACorrections(bool iselectron) : fIsElectron(iselectron) {
  // init grids related data member values
  fMaxEkin        = geant::kElectronMassC2*(1./std::sqrt(1.-gMaxBeta2)-1.);
  fLogMinEkin     = std::log(gMinEkin);
  fInvLogDelEkin  = (gNumEkin-gNumBeta2)/std::log(gMidEkin/gMinEkin);
  double pt2      = gMidEkin*(gMidEkin+2.0*geant::kElectronMassC2);
  fMinBeta2       = pt2/(pt2+geant::kElectronMassC2*geant::kElectronMassC2);
  fInvDelBeta2    = (gNumBeta2-1.)/(gMaxBeta2-fMinBeta2);
}


GSPWACorrections::~GSPWACorrections() {
  ClearDataPerElement();
  ClearDataPerMaterial();
}


void  GSPWACorrections::GetPWACorrectionFactors(double logekin, double beta2, int matindx, double &corToScr,
                                                double &corToQ1, double &corToG2PerG1) {
  int    ekinIndxLow = 0;
  double remRfaction = 0.;
  if (beta2>=gMaxBeta2) {
    ekinIndxLow = gNumEkin - 1;
    // remRfaction = -1.
  } else if (beta2>=fMinBeta2) {  // linear interpolation on \beta^2
    remRfaction   = (beta2 - fMinBeta2) * fInvDelBeta2;
    ekinIndxLow   = (int)remRfaction;
    remRfaction  -= ekinIndxLow;
    ekinIndxLow  += (gNumEkin - gNumBeta2);
  } else if (logekin>=fLogMinEkin) {
    remRfaction   = (logekin - fLogMinEkin) * fInvLogDelEkin;
    ekinIndxLow   = (int)remRfaction;
    remRfaction  -= ekinIndxLow;
  } // the defaults otherwise i.e. use the lowest energy values when ekin is smaller than the minum ekin
  //
  DataPerMaterial *data = fDataPerMaterial[matindx];
  corToScr      = data->fCorScreening[ekinIndxLow];
  corToQ1       = data->fCorFirstMoment[ekinIndxLow];
  corToG2PerG1  = data->fCorSecondMoment[ekinIndxLow];
  if (remRfaction>0.) {
    corToScr      += remRfaction*(data->fCorScreening[ekinIndxLow+1]    - data->fCorScreening[ekinIndxLow]);
    corToQ1       += remRfaction*(data->fCorFirstMoment[ekinIndxLow+1]  - data->fCorFirstMoment[ekinIndxLow]);
    corToG2PerG1  += remRfaction*(data->fCorSecondMoment[ekinIndxLow+1] - data->fCorSecondMoment[ekinIndxLow]);
  }
}


void  GSPWACorrections::Initialise(const std::vector<bool>& activeregionv) {
  // load PWA correction data for each elements that belongs to materials that are used in the detector
  InitDataPerElement(activeregionv);
  // clear  PWA correction data per material
  ClearDataPerMaterial();
  // initialise PWA correction data for the materials that are used in the detector
  InitDataPerMaterials(activeregionv);
}


void GSPWACorrections::InitDataPerElement(const std::vector<bool>& activeregionv) {
  // do it only once
  if (fDataPerElement.size()<gMaxZet+1) {
    fDataPerElement.resize(gMaxZet+1,nullptr);
  }
  // loop over all material-cuts, for those that are in region where the model is active get the list of elements and
  // load data from file if the corresponding data has not been loaded yet
  const std::vector<MaterialCuts*> &matCutTable = MaterialCuts::GetTheMaterialCutsTable();
  size_t numMatCuts = matCutTable.size();
  for (size_t imc=0; imc<numMatCuts; ++imc) {
    const MaterialCuts *matCut = matCutTable[imc];
    if (!activeregionv[matCut->GetRegionIndex()]) {
      continue;
    }
    const Vector_t<Element*> &elemVect = matCut->GetMaterial()->GetElementVector();
    //
    size_t numElems = elemVect.size();
    for (size_t ielem=0; ielem<numElems; ++ielem) {
      const Element *elem = elemVect[ielem];
      int izet = std::lrint(elem->GetZ());
      if (izet>gMaxZet) {
        izet = gMaxZet;
      }
      if (!fDataPerElement[izet]) {
        LoadDataElement(elem);
      }
    }
  }
}


void GSPWACorrections::InitDataPerMaterials(const std::vector<bool>& activeregionv) {
  // prepare size of the container
  size_t numMaterials = Material::GetNumberOMaterias();
  if (fDataPerMaterial.size()!=numMaterials) {
    fDataPerMaterial.resize(numMaterials);
  }
  // init. Mott-correction data for the Materials that are in regions where the model is active
  const std::vector<MaterialCuts*> &matCutTable = MaterialCuts::GetTheMaterialCutsTable();
  size_t numMatCuts = matCutTable.size();
  for (size_t imc=0; imc<numMatCuts; ++imc) {
    const MaterialCuts *matCut = matCutTable[imc];
    if (!activeregionv[matCut->GetRegionIndex()]) {
      continue;
    }
    const Material* mat = matCut->GetMaterial();
    if (!fDataPerMaterial[mat->GetIndex()]) {
      InitDataMaterial(mat);
    }
  }
}


// it's called only if data has not been loaded for this element yet
void GSPWACorrections::LoadDataElement(const Element *elem) {
  // allocate memory
  int izet = std::lrint(elem->GetZ());
  if (izet>gMaxZet) {
    izet = gMaxZet;
  }
  // load data from file
  char* tmppath = getenv("GEANT_PHYSICS_DATA");
  if (!tmppath) {
    std::cerr<<"******   ERROR in GSPWACorrections::LoadDataElement() \n"
             <<"         GEANT_PHYSICS_DATA is not defined! Set the GEANT_PHYSICS_DATA\n"
             <<"         environmental variable to the location of Geant data directory!\n"
             <<std::endl;
    exit(1);
  }
  std::string path(tmppath);
  if (fIsElectron) {
    path += "/msc_GS/PWACor/el/";
  } else {
    path += "/msc_GS/PWACor/pos/";
  }
  std::string   fname = path+"cf_"+gElemSymbols[izet-1];
  std::ifstream infile(fname,std::ios::in);
  if (!infile.is_open()) {
    std::cerr<<"******   ERROR in GSPWACorrections::LoadDataElement() \n"
             <<"         "<< fname << " could not be found!\n"
             <<std::endl;
    exit(1);
  }
  // allocate data structure
  DataPerMaterial *perElem = new DataPerMaterial();
  perElem->fCorScreening.resize(gNumEkin,0.0);
  perElem->fCorFirstMoment.resize(gNumEkin,0.0);
  perElem->fCorSecondMoment.resize(gNumEkin,0.0);
  fDataPerElement[izet]  = perElem;
  double dum0;
  for (int iek=0; iek<gNumEkin; ++iek) {
    infile >> dum0;
    infile >> perElem->fCorScreening[iek];
    infile >> perElem->fCorFirstMoment[iek];
    infile >> perElem->fCorSecondMoment[iek];
  }
  infile.close();
}


void GSPWACorrections::InitDataMaterial(const Material *mat) {
  constexpr double const1   = 7821.6;      // [cm2/g]
  constexpr double const2   = 0.1569;      // [cm2 MeV2 / g]
  constexpr double finstrc2 = 5.325135453E-5; // fine-structure const. square

  double constFactor        = geant::kElectronMassC2*geant::kFineStructConst/0.88534;
  constFactor              *= constFactor;  // (mc^2)^2\alpha^2/( C_{TF}^2)
  // allocate memory
  DataPerMaterial *perMat   = new DataPerMaterial();
  perMat->fCorScreening.resize(gNumEkin,0.0);
  perMat->fCorFirstMoment.resize(gNumEkin,0.0);
  perMat->fCorSecondMoment.resize(gNumEkin,0.0);
  fDataPerMaterial[mat->GetIndex()] = perMat;
  //
  const Vector_t<Element*> &elemVect = mat->GetElementVector();
  const double *nbAtomsPerVolVect = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  const double  totNbAtomsPerVol  = mat->GetMaterialProperties()->GetTotalNumOfAtomsPerVol();
  int numElems = elemVect.size();
  // 1. Compute material dependent part of Moliere's b_c \chi_c^2
  //    (with \xi=1 (i.e. total sub-threshold scattering power correction)
  double moliereBc  = 0.0;
  double moliereXc2 = 0.0;
  double zs         = 0.0;
  double ze         = 0.0;
  double zx         = 0.0;
  double sa         = 0.0;
  double xi         = 1.0;
  for (int ielem=0; ielem<numElems; ++ielem) {
    double zet = elemVect[ielem]->GetZ();
    if (zet>gMaxZet) {
      zet = (double)gMaxZet;
    }
    double iwa  = elemVect[ielem]->GetA()*geant::mole/geant::g;  // [g/mole]
    double ipz  = nbAtomsPerVolVect[ielem]/totNbAtomsPerVol;
    double dum  = ipz*zet*(zet+xi);
    zs           += dum;
    ze           += dum*(-2.0/3.0)*std::log(zet);
    zx           += dum*std::log(1.0+3.34*finstrc2*zet*zet);
    sa           += ipz*iwa;
  }
  double density = mat->GetDensity()*geant::cm3/geant::g; // [g/cm3]
  //
  moliereBc  = const1*density*zs/sa*std::exp(ze/zs)/std::exp(zx/zs);  //[1/cm]
  moliereXc2 = const2*density*zs/sa;  // [MeV2/cm]
  // change to Geant4 internal units of 1/length and energ2/length
  moliereBc  *= 1.0/geant::cm;
  moliereXc2 *= geant::MeV*geant::MeV/geant::cm;
  //
  // 2. loop over the kinetic energy grid
  for (int iek=0; iek<gNumEkin; ++iek) {
    // 2./a. set current kinetic energy and pt2 value
      double ekin = std::exp(fLogMinEkin+iek/fInvLogDelEkin);
      double pt2  = ekin*(ekin+2.0*geant::kElectronMassC2);
      if (ekin>gMidEkin) {
        double b2   = fMinBeta2+(iek-(gNumEkin-gNumBeta2))/fInvDelBeta2;
        ekin = geant::kElectronMassC2*(1./std::sqrt(1.-b2)-1.);
        pt2  = ekin*(ekin+2.0*geant::kElectronMassC2);
      }
    // 2./b. loop over the elements at the current kinetic energy point
    for (int ielem=0; ielem<numElems; ++ielem) {
      const Element *elem = elemVect[ielem];
      double zet  = elem->GetZ();
      if (zet>gMaxZet) {
        zet = (double)gMaxZet;
      }
      int izet = std::lrint(zet);
      // loaded PWA corrections for the current element
      DataPerMaterial *perElem  = fDataPerElement[izet];
      //
      // xi should be one i.e. z(z+1) since total sub-threshold scattering power correction
      double nZZPlus1  = nbAtomsPerVolVect[ielem]*zet*(zet+1.0)/totNbAtomsPerVol;
      double Z23       = std::pow(zet,2./3.);
      //
      // 2./b./(i) Add the 3 PWA correction factors
      double mcScrCF = perElem->fCorScreening[iek];     // \kappa_i[1.13+3.76(\alpha Z_i)^2] with \kappa_i=scr_mc/scr_sr
      // compute the screening parameter correction factor (Z_i contribution to the material)
      // src_{mc} = C \exp\left[ \frac{ \sum_i n_i Z_i(Z_i+1)\ln[Z_{i}^{2/3}\kappa_i(1.13+3.76(\alpha Z_i)^2)] } {\sum_i n_i Z_i(Z_i+1)}
      // with C = \frac{(mc^2)^\alpha^2} {4(pc)^2 C_{TF}^2} = constFactor/(4*(pc)^2)
      // here we compute the \sum_i n_i Z_i(Z_i+1)\ln[Z_{i}^{2/3}\kappa_i(1.13+3.76(\alpha Z_i)^2)] part
      perMat->fCorScreening[iek] += nZZPlus1*std::log(Z23*mcScrCF);
      // compute the corrected screening parameter for the current Z_i and E_{kin}
      // src(Z_i)_{mc} = \frac{(mc^2)^\alpha^2 Z_i^{2/3}} {4(pc)^2 C_{TF}^2} \kappa_i[1.13+3.76(\alpha Z_i)^2]
      mcScrCF *= constFactor*Z23/(4.*pt2);
      // compute first moment correction factor
      // q1_{mc} = \frac{ \sum_i n_i Z_i(Z_i+1) A_i  B_i } {\sum_i n_i Z_i(Z_i+1)} \frac{1}{C}
      // where:
      // A_i(src(Z_i)_{mc}) = [\ln(1+1/src(Z_i)_{mc}) - 1/(1+src(Z_i)_{mc})]; where \sigma(Z_i)_{tr1}^(sr) = A_i(src(Z_i)_{mc}) [2\pi r_0 Z_i mc^2/(pc)\beta]^2
      // B_i = \beta_i \gamma_i with beta_i(Z_i) = \sigma(Z_i)_{tr1}^(PWA)/\sigma(Z_i,src(Z_i)_{mc})_{tr1}^(sr)
      // and \gamma_i = \sigma(Z_i)_{el}^(MC-DCS)/\sigma(Z_i,src(Z_i)_{mc})_{el}^(sr)
      // C(src_{mc}) = [\ln(1+1/src_{mc}) - 1/(1+src_{mc})]; where \sigma_{tr1}^(sr) = C(src_{mc}) [2\pi r_0 Z_i mc^2/(pc)\beta]^2
      // A_i x B_i is stored in file per e-/e+, E_{kin} and Z_i
      // here we compute the \sum_i n_i Z_i(Z_i+1) A_i  B_i part
      perMat->fCorFirstMoment[iek] += nZZPlus1*(std::log(1.+1./mcScrCF)-1./(1.+mcScrCF))*perElem->fCorFirstMoment[iek];
      // compute the second moment correction factor
      // [G2/G1]_{mc} = \frac{ \sum_i n_i Z_i(Z_i+1) A_i } {\sum_i n_i Z_i(Z_i+1)} \frac{1}{C}
      // with A_i(Z_i) = G2(Z_i)^{PWA}/G1(Z_i)^{PWA} and C=G2(Z_i,scr_{mc})^{sr}/G1(Z_i,scr_{mc})^{sr}}
      // here we compute the \sum_i n_i Z_i(Z_i+1) A_i part
      perMat->fCorSecondMoment[iek] += nZZPlus1*perElem->fCorSecondMoment[iek];
      //
      // 2./b./(ii) When the last element has been added:
      if (ielem==numElems-1) {
        //
        // 1. the remaining part of the sreening correction and divide the corrected screening par. with Moliere's one:
        //    (Moliere screening parameter = moliereXc2/(4(pc)^2 moliereBc) )
        double dumScr   = std::exp(perMat->fCorScreening[iek]/zs);
        perMat->fCorScreening[iek] = constFactor*dumScr*moliereBc/moliereXc2;
        //
        // 2. the remaining part of the first moment correction and divide by the one computed by using the corrected
        //    screening parameter (= (mc^2)^\alpha^2/(4(pc)^2C_{TF}^2) dumScr
        double scrCorTed = constFactor*dumScr/(4.*pt2);
        double dum0      = std::log(1.+1./scrCorTed);
        perMat->fCorFirstMoment[iek] = perMat->fCorFirstMoment[iek]/(zs*(dum0-1./(1.+scrCorTed)));
        //
        // 3. the remaining part of the second moment correction and divide by the one computed by using the corrected
        //    screening parameter
        double G2PerG1   =  3.*(1.+scrCorTed)*((1.+2.*scrCorTed)*dum0-2.)/((1.+scrCorTed)*dum0-1.);
        perMat->fCorSecondMoment[iek] = perMat->fCorSecondMoment[iek]/(zs*G2PerG1);
      }
    }
  }
}



void GSPWACorrections::ClearDataPerElement() {
  for (size_t i=0; i<fDataPerElement.size(); ++i) {
    if (fDataPerElement[i]) {
      fDataPerElement[i]->fCorScreening.clear();
      fDataPerElement[i]->fCorFirstMoment.clear();
      fDataPerElement[i]->fCorSecondMoment.clear();
      delete fDataPerElement[i];
    }
  }
  fDataPerElement.clear();
}


void GSPWACorrections::ClearDataPerMaterial() {
  for (size_t i=0; i<fDataPerMaterial.size(); ++i) {
    if (fDataPerMaterial[i]) {
      fDataPerMaterial[i]->fCorScreening.clear();
      fDataPerMaterial[i]->fCorFirstMoment.clear();
      fDataPerMaterial[i]->fCorSecondMoment.clear();
      delete fDataPerMaterial[i];
    }
  }
  fDataPerMaterial.clear();
}

} // namespace geantphysics
