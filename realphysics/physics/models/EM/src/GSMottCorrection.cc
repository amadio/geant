
#include "GSMottCorrection.h"

#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

#include "Material.h"
#include "MaterialProperties.h"
#include "Element.h"
#include "MaterialCuts.h"

// only for rng
#include "GeantTaskData.h"

#include <iostream>
#include <fstream>
#include <cmath>

namespace geantphysics {

const std::string GSMottCorrection::gElemSymbols[] = {"H","He","Li","Be","B" ,
 "C" ,"N" ,"O" ,"F" ,"Ne","Na","Mg","Al","Si","P" , "S","Cl","Ar","K" ,"Ca","Sc",
 "Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb",
 "Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,
 "Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
 "Yb","Lu","Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
 "Rn","Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf"};

GSMottCorrection::GSMottCorrection(bool iselectron) : fIsElectron(iselectron) {
  // init grids related data member values
  fMaxEkin        = geant::kElectronMassC2*(1./std::sqrt(1.-gMaxBeta2)-1.);
  fLogMinEkin     = std::log(gMinEkin);
  fInvLogDelEkin  = (gNumEkin-gNumBeta2)/std::log(gMidEkin/gMinEkin);
  double pt2      = gMidEkin*(gMidEkin+2.0*geant::kElectronMassC2);
  fMinBeta2       = pt2/(pt2+geant::kElectronMassC2*geant::kElectronMassC2);
  fInvDelBeta2    = (gNumBeta2-1.)/(gMaxBeta2-fMinBeta2);
  fInvDelDelta    = (gNumDelta-1.)/gMaxDelta;
  fInvDelAngle    = gNumAngle-1.;
}


GSMottCorrection::~GSMottCorrection() {
  ClearMCDataPerElement();
  ClearMCDataPerMaterial();
}


void GSMottCorrection::GetMottCorrectionFactors(double logekin, double beta2, int matindx, double &mcToScr,
                                                double &mcToQ1, double &mcToG2PerG1) {
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
  DataPerEkin *perEkinLow  = fMCDataPerMaterial[matindx]->fDataPerEkin[ekinIndxLow];
  mcToScr      = perEkinLow->fMCScreening;
  mcToQ1       = perEkinLow->fMCFirstMoment;
  mcToG2PerG1  = perEkinLow->fMCSecondMoment;
  if (remRfaction>0.) {
    DataPerEkin *perEkinHigh = fMCDataPerMaterial[matindx]->fDataPerEkin[ekinIndxLow+1];
    mcToScr      += remRfaction*(perEkinHigh->fMCScreening    - perEkinLow->fMCScreening);
    mcToQ1       += remRfaction*(perEkinHigh->fMCFirstMoment  - perEkinLow->fMCFirstMoment);
    mcToG2PerG1  += remRfaction*(perEkinHigh->fMCSecondMoment - perEkinLow->fMCSecondMoment);
  }
}


// accept cost if rndm [0,1] < return value
double GSMottCorrection::GetMottRejectionValue(double logekin, double beta2, double q1, double cost, int matindx,
                                               int &ekindx, int &deltindx, Geant::GeantTaskData* td) {
  double val   = 1.0;
  double delta = q1/(0.5+q1);
  // check if converged to 1 for all angles => accept cost
  if (delta>=gMaxDelta) {
    return val;
  }
  //
  // check if kinetic energy index needs to be determined
  if (ekindx<0) {
    int    ekinIndxLow  = 0;
    double probIndxHigh = 0.;  // will be the prob. of taking the ekinIndxLow+1 bin
    if (beta2>gMaxBeta2) {
      ekinIndxLow = gNumEkin - 1;
      // probIndxHigh = -1.
    } else if (beta2>=fMinBeta2) {    // linear interpolation on \beta^2
      probIndxHigh  = (beta2 - fMinBeta2) * fInvDelBeta2;
      ekinIndxLow   = (int)probIndxHigh;
      probIndxHigh -= ekinIndxLow;
      ekinIndxLow  += (gNumEkin - gNumBeta2);
    } else if (logekin>fLogMinEkin) { // linear interpolation on \ln(E_{kin})
      probIndxHigh  = (logekin - fLogMinEkin) * fInvLogDelEkin;
      ekinIndxLow   = (int)probIndxHigh;
      probIndxHigh -= ekinIndxLow;
    } // the defaults otherwise i.e. use the lowest energy values when ekin is smaller than the minum ekin
    //
    // check if need to take the higher ekin index
    if (td->fRndm->uniform()<probIndxHigh) {
      ++ekinIndxLow;
    }
    // set kinetic energy grid index
    ekindx = ekinIndxLow;
  }
  // check if delta value index needs to be determined (note: in case of single scattering deltindx will be set to 0 by
  // by the caller but the ekindx will be -1: kinetic energy index is not known but the delta index is known)
  if (deltindx<0) {
    // note: delta is for sure < gMaxDelta at this point ( and minimum delta value is 0)
    double probIndxHigh = delta*fInvDelDelta;  // will be the prob. of taking the deltIndxLow+1 bin
    int    deltIndxLow  = (int)probIndxHigh;
    probIndxHigh       -= deltIndxLow;
    // check if need to take the higher delta index
    if (td->fRndm->uniform()<probIndxHigh) {
      ++deltIndxLow;
    }
    // set the delta value grid index
    deltindx = deltIndxLow;
  }
  //
  // get the corresponding distribution
  DataPerDelta *perDelta  = fMCDataPerMaterial[matindx]->fDataPerEkin[ekindx]->fDataPerDelta[deltindx];
  //
  // determine lower index of the angular bin
  double ang         = std::sqrt(0.5*(1.-cost)); // sin(0.5\theta) in [0,1]
  double remRfaction = ang*fInvDelAngle;
  int    angIndx     = (int)remRfaction;
  remRfaction       -= angIndx;
  if (angIndx<gNumAngle-2) { // normal case: linear interpolation
    val        = remRfaction*(perDelta->fRejFuntion[angIndx+1]-perDelta->fRejFuntion[angIndx]) + perDelta->fRejFuntion[angIndx];
  } else {   // last bin
    double dum = ang-1.+1./fInvDelAngle;
    val        = perDelta->fSA + dum*(perDelta->fSB + dum*(perDelta->fSC + dum*perDelta->fSD));
  }
  return val;
}


void GSMottCorrection::Initialise(const std::vector<bool>& activeregionv) {
  // load Mott-correction data for each elements that belongs to materials that are used in the detector
  InitMCDataPerElement(activeregionv);
  // clrea Mott-correction data per material
  ClearMCDataPerMaterial();
  // initialise Mott-correction data for the materials that are used in the detector
  InitMCDataPerMaterials(activeregionv);
}


void GSMottCorrection::InitMCDataPerElement(const std::vector<bool>& activeregionv) {
  // do it only once
  if (fMCDataPerElement.size()<gMaxZet+1) {
    fMCDataPerElement.resize(gMaxZet+1,nullptr);
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
      if (!fMCDataPerElement[izet]) {
        LoadMCDataElement(elem);
      }
    }
  }
}


void GSMottCorrection::InitMCDataPerMaterials(const std::vector<bool>& activeregionv) {
  // prepare size of the container
  size_t numMaterials = Material::GetNumberOMaterias();
  if (fMCDataPerMaterial.size()!=numMaterials) {
    fMCDataPerMaterial.resize(numMaterials);
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
    if (!fMCDataPerMaterial[mat->GetIndex()]) {
      InitMCDataMaterial(mat);
    }
  }
}


// it's called only if data has not been loaded for this element yet
void GSMottCorrection::LoadMCDataElement(const Element *elem) {
  // allocate memory
  int izet = std::lrint(elem->GetZ());
  if (izet>gMaxZet) {
    izet = gMaxZet;
  }
  DataPerMaterial *perElem = new DataPerMaterial();
  AllocateDataPerMaterial(perElem);
  fMCDataPerElement[izet]  = perElem;
  //
  // load data from file
  char* tmppath = getenv("GEANT_PHYSICS_DATA");
  if (!tmppath) {
    std::cerr<<"******   ERROR in GSMottCorrection::LoadMCDataElement() \n"
             <<"         GEANT_PHYSICS_DATA is not defined! Set the GEANT_PHYSICS_DATA\n"
             <<"         environmental variable to the location of Geant data directory!\n"
             <<std::endl;
    exit(1);
  }
  std::string path(tmppath);
  if (fIsElectron) {
    path += "/msc_GS/MottCor/el/";
  } else {
    path += "/msc_GS/MottCor/pos/";
  }
  std::string fname = path+"rej_"+gElemSymbols[izet-1];
  std::ifstream infile(fname,std::ios::in);
  if (!infile.is_open()) {
    std::cerr<<"******   ERROR in GSMottCorrection::LoadMCDataElement() \n"
             <<"         "<< fname << " could not be found!\n"
             <<std::endl;
    exit(1);
  }
  for (int iek=0; iek<gNumEkin; ++iek) {
    DataPerEkin *perEkin = perElem->fDataPerEkin[iek];
    // 1. get the 3 Mott-correction factors for the current kinetic energy
    infile >> perEkin->fMCScreening;
    infile >> perEkin->fMCFirstMoment;
    infile >> perEkin->fMCSecondMoment;
    // 2. load each data per delta:
    for (int idel=0; idel<gNumDelta; ++idel) {
      DataPerDelta *perDelta = perEkin->fDataPerDelta[idel];
      // 2./a.  : first the rejection function values
      for (int iang=0; iang<gNumAngle; ++iang) {
        infile >> perDelta->fRejFuntion[iang];
      }
      // 2./b. : then the 4 spline parameter for the last bin
      infile >> perDelta->fSA;
      infile >> perDelta->fSB;
      infile >> perDelta->fSC;
      infile >> perDelta->fSD;
    }
  }
  infile.close();
}



void GSMottCorrection::InitMCDataMaterial(const Material *mat) {
  constexpr double const1   = 7821.6;      // [cm2/g]
  constexpr double const2   = 0.1569;      // [cm2 MeV2 / g]
  constexpr double finstrc2 = 5.325135453E-5; // fine-structure const. square

  double constFactor        = geant::kElectronMassC2*geant::kFineStructConst/0.88534;
  constFactor              *= constFactor;  // (mc^2)^2\alpha^2/( C_{TF}^2)
  // allocate memory
  DataPerMaterial *perMat   = new DataPerMaterial();
  AllocateDataPerMaterial(perMat);
  fMCDataPerMaterial[mat->GetIndex()] = perMat;
  //
  const Vector_t<Element*> &elemVect = mat->GetElementVector();
  const double *nbAtomsPerVolVect = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
  const double  totNbAtomsPerVol  = mat->GetMaterialProperties()->GetTotalNumOfAtomsPerVol();
  int numElems = elemVect.size();
  //
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
    zs         += dum;
    ze         += dum*(-2.0/3.0)*std::log(zet);
    zx         += dum*std::log(1.0+3.34*finstrc2*zet*zet);
    sa         += ipz*iwa;
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
      int izet         = std::lrint(zet);
      // xi should be one i.e. z(z+1) since total sub-threshold scattering power correction
      double nZZPlus1  = nbAtomsPerVolVect[ielem]*zet*(zet+1.0)/totNbAtomsPerVol;
      double Z23       = std::pow(zet,2./3.);
      //
      DataPerEkin *perElemPerEkin  = fMCDataPerElement[izet]->fDataPerEkin[iek];
      DataPerEkin *perMatPerEkin   = perMat->fDataPerEkin[iek];
      //
      // 2./b./(i) Add the 3 Mott-correction factors
      double mcScrCF = perElemPerEkin->fMCScreening;     // \kappa_i[1.13+3.76(\alpha Z_i)^2] with \kappa_i=scr_mc/scr_sr
      // compute the screening parameter correction factor (Z_i contribution to the material)
      // src_{mc} = C \exp\left[ \frac{ \sum_i n_i Z_i(Z_i+1)\ln[Z_{i}^{2/3}\kappa_i(1.13+3.76(\alpha Z_i)^2)] } {\sum_i n_i Z_i(Z_i+1)}
      // with C = \frac{(mc^2)^\alpha^2} {4(pc)^2 C_{TF}^2} = constFactor/(4*(pc)^2)
      // here we compute the \sum_i n_i Z_i(Z_i+1)\ln[Z_{i}^{2/3}\kappa_i(1.13+3.76(\alpha Z_i)^2)] part
      perMatPerEkin->fMCScreening += nZZPlus1*std::log(Z23*mcScrCF);
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
      perMatPerEkin->fMCFirstMoment += nZZPlus1*(std::log(1.+1./mcScrCF)-1./(1.+mcScrCF))*perElemPerEkin->fMCFirstMoment;
      // compute the second moment correction factor
      // [G2/G1]_{mc} = \frac{ \sum_i n_i Z_i(Z_i+1) A_i } {\sum_i n_i Z_i(Z_i+1)} \frac{1}{C}
      // with A_i(Z_i) = G2(Z_i)^{PWA}/G1(Z_i)^{PWA} and C=G2(Z_i,scr_{mc})^{sr}/G1(Z_i,scr_{mc})^{sr}}
      // here we compute the \sum_i n_i Z_i(Z_i+1) A_i part
      perMatPerEkin->fMCSecondMoment += nZZPlus1*perElemPerEkin->fMCSecondMoment;
      //
      // 2./b./(ii) Go for the rejection funtion part
      // I. loop over delta values
      for (int idel=0; idel<gNumDelta; ++idel) {
        DataPerDelta *perMatPerDelta  = perMatPerEkin->fDataPerDelta[idel];
        DataPerDelta *perElemPerDelta = perElemPerEkin->fDataPerDelta[idel];
        // I./a. loop over angles (i.e. the \sin(0.5\theta) values) and add the rejection function
        for (int iang=0; iang<gNumAngle; ++iang) {
          perMatPerDelta->fRejFuntion[iang] += nZZPlus1*perElemPerDelta->fRejFuntion[iang];
        }
        // I./b. get the last bin spline parameters and add them (a+bx+cx^2+dx^3)
        perMatPerDelta->fSA += nZZPlus1*perElemPerDelta->fSA;
        perMatPerDelta->fSB += nZZPlus1*perElemPerDelta->fSB;
        perMatPerDelta->fSC += nZZPlus1*perElemPerDelta->fSC;
        perMatPerDelta->fSD += nZZPlus1*perElemPerDelta->fSD;
      }
      //
      // 2./b./(iii) When the last element has been added:
      if (ielem==numElems-1) {
        //
        // 1. the remaining part of the sreening correction and divide the corrected screening par. with Moliere's one:
        //    (Moliere screening parameter = moliereXc2/(4(pc)^2 moliereBc) )
        double dumScr   = std::exp(perMatPerEkin->fMCScreening/zs);
        perMatPerEkin->fMCScreening = constFactor*dumScr*moliereBc/moliereXc2;
        //
        // 2. the remaining part of the first moment correction and divide by the one computed by using the corrected
        //    screening parameter (= (mc^2)^\alpha^2/(4(pc)^2C_{TF}^2) dumScr
        double scrCorTed = constFactor*dumScr/(4.*pt2);
        double dum0      = std::log(1.+1./scrCorTed);
        perMatPerEkin->fMCFirstMoment = perMatPerEkin->fMCFirstMoment/(zs*(dum0-1./(1.-scrCorTed)));
        //
        // 3. the remaining part of the second moment correction and divide by the one computed by using the corrected
        //    screening parameter
        double G2PerG1   =  3.*(1.+scrCorTed)*((1.+2.*scrCorTed)*dum0-2.)/((1.+scrCorTed)*dum0-1.);
        perMatPerEkin->fMCSecondMoment = perMatPerEkin->fMCSecondMoment/(zs*G2PerG1);
        //
        // 4. scale the maximum of the rejection function to unity and correct the last bin spline parameters as well
        // I. loop over delta values
        for (int idel=0; idel<gNumDelta; ++idel) {
          DataPerDelta *perMatPerDelta  = perMatPerEkin->fDataPerDelta[idel];
          double maxVal = -1.;
          // II. llop over angles
          for (int iang=0; iang<gNumAngle; ++iang) {
            if (perMatPerDelta->fRejFuntion[iang]>maxVal)
              maxVal = perMatPerDelta->fRejFuntion[iang];
          }
          for (int iang=0; iang<gNumAngle; ++iang) {
            perMatPerDelta->fRejFuntion[iang] /=maxVal;
          }
          perMatPerDelta->fSA /= maxVal;
          perMatPerDelta->fSB /= maxVal;
          perMatPerDelta->fSC /= maxVal;
          perMatPerDelta->fSD /= maxVal;
        }
      }
    }
  }
}


void GSMottCorrection::AllocateDataPerMaterial(DataPerMaterial *data) {
  data->fDataPerEkin = new DataPerEkin*[gNumEkin]();
  for (int iek=0; iek<gNumEkin; ++iek) {
    DataPerEkin *perEkin   = new DataPerEkin();
    perEkin->fDataPerDelta = new DataPerDelta*[gNumDelta]();
    for (int idel=0; idel<gNumDelta; ++idel) {
      DataPerDelta *perDelta       = new DataPerDelta();
      perDelta->fRejFuntion        = new double[gNumAngle]();
      perEkin->fDataPerDelta[idel] = perDelta;
    }
    data->fDataPerEkin[iek] = perEkin;
  }
}

void GSMottCorrection::DeAllocateDataPerMaterial(DataPerMaterial *data) {
  for (int iek=0; iek<gNumEkin; ++iek) {
    DataPerEkin *perEkin = data->fDataPerEkin[iek]; //new DataPerEkin();
    for (int idel=0; idel<gNumDelta; ++idel) {
      DataPerDelta *perDelta = perEkin->fDataPerDelta[idel];
      delete [] perDelta->fRejFuntion;
      delete perDelta;
    }
    delete [] perEkin->fDataPerDelta;
    delete perEkin;
  }
  delete [] data->fDataPerEkin;
}


void GSMottCorrection::ClearMCDataPerElement() {
  for (size_t i=0; i<fMCDataPerElement.size(); ++i) {
    if (fMCDataPerElement[i]) {
      DeAllocateDataPerMaterial(fMCDataPerElement[i]);
      delete fMCDataPerElement[i];
    }
  }
  fMCDataPerElement.clear();
}

void GSMottCorrection::ClearMCDataPerMaterial() {
  for (size_t i=0; i<fMCDataPerMaterial.size(); ++i) {
    if (fMCDataPerMaterial[i]) {
      DeAllocateDataPerMaterial(fMCDataPerMaterial[i]);
      delete fMCDataPerMaterial[i];
    }
  }
  fMCDataPerMaterial.clear();
}

} // namespace geantphysics
