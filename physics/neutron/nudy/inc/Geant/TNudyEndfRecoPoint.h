//===-- Nudy/TNudyEndfRecoPoint.h - Instruction class definition -------*- C++ -*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \brief The class reads doppler broadened cross-sections from the root file
///  and gives all the corss-sections/ secondary parameters forthe outgoing particles
/// \class TNudyEndfRecoPoint
/// \author H. Kumawat
/// \date March 2016
//===----------------------------------------------------------------------===//
#ifndef __TNudyEndfRecoPoint__
#define __TNudyEndfRecoPoint__

#include <vector>
#include <fstream>
#include "Geant/RngWrapper.h"
namespace Nudy {
  class TNudyEndfFile;
  class TNudyEndfList;
  class TNudyEndfTab1;
  class TNudyEndfTab2;
}

class TList;

namespace NudyPhysics {
  class TNudyEndfNuPh;
  class TNudyEndfEnergyAng;
  class TNudyEndfPhYield;
  class TNudyEndfPhProd;
  class TNudyEndfPhAng;
  class TNudyEndfPhEnergy;
}

#ifdef USE_ROOT
#include "Rtypes.h"
#endif

namespace NudyPhysics {
  
  typedef std::vector<double> rowd;
  typedef std::vector<int> rowint;
  typedef std::vector<rowint> matrixint;
  typedef std::vector<rowd> matrixd2;
  typedef std::vector<std::vector<rowd>> matrixd3;
  typedef std::vector<std::vector<std::vector<rowd>>> matrixd4;
  typedef std::vector<std::vector<std::vector<std::vector<rowd>>>> matrixd5;
  
  class TNudyEndfRecoPoint {
    
  public:
    TNudyEndfRecoPoint();
    /// \brief Default constructure
    TNudyEndfRecoPoint(int ielemId, const char *irENDF);
    /// \brief constructure to be called for any interface class to get data
    virtual ~TNudyEndfRecoPoint();
    void GetData(int elemid, const char *irENDF);
    /// \brief main function to process the data
    double GetSigmaTotal(int elemid, double energyK);
    /// \brief getting total XSec.
    double GetSigmaPartial(int elemid, int i, double energyK);
    /// \brief getting partial cross-section
    double GetCos4(int elemid, int mt, double energyK);
    /// \brief getting cosine from file 4
    int GetCos4Lct(int elemid, int mt);
    /// \brief getting CM vs. LAB flag
    double GetEnergy5(int elemid, int mt, double energyK);
    /// \brief getting secondary neutron energy from file 5
    double GetDelayedFraction(int ielemId, int mt, double energyK);
    /// \brief getting delayed fraction for each group of the emitter family
    int GetLaw6(int ielemId, int mt);
    /// \brief getting LAW of the energy distribution in file 6
    int GetZd6(int ielemId, int mt);
    /// \brief getting Z value
    int GetAd6(int ielemId, int mt);
    /// \brief getting A value
    int GetMt6Neutron(int ielemId, int mt);
    /// \brief getting MT values for neutron only
    //double GetCos64(int elemid, int mt, double energyK);
    /// \brief getting cosine from file 6 of type 4
    double GetCos6(int elemid, int mt, double energyK);
    /// \brief getting cosine from file 6
    double GetEnergy6(int elemid, int mt, double energyK);
    /// \brief getting energy from file 6
    double GetCos4Photon(int elemid, int mt, double energyK);
    /// \brief getting cosine from file 4 for photon
    int GetCos4LctPhoton(int elemid, int mt);
    /// \brief getting CM vs. LAB flag for photon
    double GetEnergy5Photon(int elemid, int mt, double energyK);
    /// \brief getting secondary photon energy from file 5
    virtual double GetQValue(int elemid, int mt);
    /// \brief getting q value for the reaction specially for the excited state emission
    virtual double GetMt4(int elemid, int mt);
    /// \brief getting MT values for file 4 for which angular distributions are given
    virtual double GetMt5(int elemid, int mt);
    /// \brief getting MT values for which energy distributions are given in file 5
    virtual double GetMt6(int elemid, int mt);
    /// \brief getting MT values for which angle-energy correlated distributions are given in file 6
    virtual double GetNuTotal(int elemid, double energyK);
    /// \brief getting total fission neutron multiplicity
    virtual double GetNuPrompt(int elemid, double energyK);
    /// \brief getting prompt fission neutron multiplicity
    virtual double GetNuDelayed(int elemid, double energyK);
    /// \brief getting delayed fission neutron multiplicity
    virtual double GetFissHeat(int elemid, double energyK);
    /// \brief getting total fission Heat
    double GetFisYield(int elemid, double energyK);
    /// \brief getting fission yield
    virtual double GetLambdaD(int elemid, int time);
    /// \brief getting lambda of delayed neutron emitter family
    matrixint fMtValues;
    /// \brief MT values for which cross-section/ heating values are given
    rowint fElementId;
    /// \brief z*1000+A of the element in the material list of the simulation
    int GetElementId(int elemid);
    /// \brief getting isotope index for stored vector to retrive data
  protected:
    int fElemId;
    /// \brief serial no of the element in the material list of the simulation
    const char *rENDF;
    /// \brief root file name
    matrixd2 fEneUni, fSigUniT;
    /// \brief unionization of energy and total cross-section
    matrixd3 fSigUniOfMt;
    /// \brief Xsec for each reaction after unionization of energy for all reactions
    matrixint fEnergyLocMtId;
    /// \brief MT wise starting energy locations for all cross-sections
    matrixd4 fCos4OfMts;
    /// \brief cosine from file 4 for each reaction
    matrixd4 fCosPdf4OfMts;
    /// \brief pdf from file 4 for each reaction
    matrixd4 fCosCdf4OfMts;
    /// \brief cdf from file 4 for each reaction
    matrixd3 fEnergy4OfMts;
    /// \brief incident energy in file 4 for each reaction
    matrixint fMt4Values;
    ///  \brief MT values for which angular distributions are given in file 4
    matrixint fMt4Lct;
    /// \brief CM and Lab flag for angular distributions as given in file 4
    matrixd4 fCos4OfMtsPhoton;
    /// \brief 4-D vector: photon cosine in file 4 for each given reaction and element
    matrixd4 fCosPdf4OfMtsPhoton;
    /// \brief cosine photon pdf in file 4 for each given reaction and element
    matrixd4 fCosCdf4OfMtsPhoton;
    /// \brief cosine photon cdf from file 4 for each given reaction and element
    matrixd3 fEnergy4OfMtsPhoton;
    /// \brief incident neutron energy in file 4 for which cosine angles are given
    matrixint fMt4ValuesPhoton;
    /// \brief MT values for which angular distributions of photons are given in file 4
    matrixint fMt4LctPhoton;
    /// CM and Lab flag for angular distributions of photons as given in file 4
    matrixd4 fEnergyOut5OfMts, fEnergyOut5OfMtsPhoton;
    /// \brief energy from file 5 for each reaction
    matrixd4 fEnergyPdf5OfMts, fEnergyPdf5OfMtsPhoton;
    /// \brief pdf from file 5 for each reaction
    matrixd4 fEnergyCdf5OfMts, fEnergyCdf5OfMtsPhoton;
    /// \brief cdf from file 5 for each reaction
    matrixd4 fCos6OfMts;
    /// \brief cosine 6 for each reaction
    matrixd4 fCosin6Pdf, fCosin6Cdf;
    /// \brief pdf cdf file 6 for each reaction
    matrixd5 fEnergyOut6OfMts;
    /// \brief energy from file 6 for each reaction
    matrixd5 fEnergyPdf6OfMts;
    /// \brief pdf from file 6 for each reaction
    matrixd5 fEnergyCdf6OfMts;
    /// \brief cdf from file 6 for each reaction
    matrixd3 fEnergy5OfMts, fEnergy5OfMtsPhoton, fEnergy6OfMts;
    /// \brief incident energy in file 5 for each reaction
    matrixd3 fFraction5OfMts, fFraction5OfMtsPhoton;
    /// \brief fraction for incident energy in file 5 for each reaction
    matrixint fMt5Values, fMt5ValuesPhoton;
    /// \brief MT values for which energy distributions are given in file 5
    matrixint fLaw6;
    /// \brief law 6 for angular-energy distributions are given in file 6
    matrixint fZD6, fAD6;
    /// \brief Z, A of law 6 distributions
    matrixd2 fEint, fNut;
    /// \brief total incident energy and nu
    matrixd2 fEinp, fNup;
    /// \brief prompt incident energy and nu
    matrixd2 fEind, fNud, fLambdaD;
    /// \brief delayed incident energy, nu and lambda
    matrixd2 fEinFissHeat, fFissHeat;
    /// \brief fission incident energy and heat
    matrixd2 fEinfId, fQvalue;
    /// \brief incident energy for fission yield and q-value
    matrixd3 fZafId, fPdfYieldId, fCdfYieldId;
    /// \brief za and yield fission
    matrixint fMt6Neutron, fMt6Photon, fMt6Charge;
    /// \brief MT numbers if it is for neutron photon or charge particle
    matrixint fMt6Values;
    /// \brief MT numbers in file 6
    double AWRI;
    /// \brief Mass in units of the neutron
    rowd fEinc;
    /// \brief incident energy for the cummulative yield
    matrixd2 fZafp, fFps, fZafpc, fFpsc, fYi, fCyi, fDyi, fYc, fDyc;
    /// \brief charge, mass, yield (independent and cummulative and error)
    rowd fZafp1, fFps1, fZafpc1, fFpsc1, fYi1, fCyi1, fDyi1, fYc1, fDyc1;
    /// \brief charge, mass, yield (independent and cummulative)
    int fNK, fNI;
    /// \brief Total no. of discrete + continuous photon distributions, discrete distributions
    
  private:
    void ReadFile2(Nudy::TNudyEndfFile *file);
    /// \brief function to read file 2 (resonance parameters)
    void ReadFile3(Nudy::TNudyEndfFile *file);
    /// \brief function to read file 3 (cross-section)
    void FixupTotal(rowd &x1);
    /// \brief making total cross-section
    void ReadFile4(Nudy::TNudyEndfFile *file);
    /// \brief function to read file 4 (angular distribution)
    void ReadFile5(Nudy::TNudyEndfFile *file);
    /// \brief function to read file 5 (energy distribution)
    void ReadFile6(Nudy::TNudyEndfFile *file);
    /// \brief function to read file 6 (angle-energy distribution)
    void ReadFile8(Nudy::TNudyEndfFile *file);
    /// \brief function to read file 8 (fission yield)
    void ReadFile14(Nudy::TNudyEndfFile *file);
    /// \brief function to read file 14 (photon angular distributions)
    void ReadFile15(Nudy::TNudyEndfFile *file);
    /// \brief function to read file 15 (photon energy distributions)
    void ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, rowint &fint,rowd &x1, rowd &x2);
    /// \brief process tab1 entry
    void ProcessTab2(Nudy::TNudyEndfTab2 *tab2, int &NR,int &NP,rowint &fnbt, rowint &fint); 
    /// \brief process tab2 entry 
    void ModifyTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x);
    /// \brief modify Tab1 record after linearization
    void CreateTab2(Nudy::TNudyEndfTab2 *secTab2, int &NE);
    /// \brief create tab2 record in file5 (MF = 5) to make LF=1 structure
    void CreateTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x);
    /// \brief create tab1 record in file5 (MF = 5) to make LF=1 structure
    void FillPdf1D(rowd &x1, rowd &x2, rowd &x3, matrixd2 &x4, matrixd2 &x5, matrixd2 &x6);
    /// \brief filling 1 dimentional pdf  distribution
    void FillPdf2D();
    /// \brief filling 2 dimentional pdf  distribution
    
    int mtf[3];
    /// \brief MAT, MT, MF parameters
    int fFlagRead = -1;
    /// \brief flag for reading charge particle and photon production cross-sections
    double fQValue[999];
    /// \brief q-value for all the reactions
    int fNR, fNP;
    /// \brief standard ENDF parameters for range and interpolation
    matrixint fMt4, fMt5, fMt6;
    /// \brief MT values for which angular, energy/ angular-energy distributions are given in file 4, 5, 6
    matrixd2 fSigmaOfMts;
    /// \brief sigma for each reaction
    matrixd2 fSigmaUniOfMts;
    /// \brief sigma for each reaction after unionization of energy
    rowint fEnergyLocationMts;
    /// \breif MT wise starting energy for cross-section
    rowint fMtNumbers, fMtNum4, fMtNum5, fMtNum6;
    /// \breif MT numbers
    rowd fEnergyMts, fSigmaMts, fQvalueTemp;
    /// \brief MT numbers for sigma in file3
    rowd fELinearFile3, fXLinearFile3;
    /// \breif energy and Xsec
    rowd fEneTemp, fSigTemp;
    /// \breif temporary vectors to store energy and sigma
    rowd fE1, fP1, fE2, fP2, fE3, fP3, INorm;
    /// \brief standard file 5 parameters from ENDF manual
    rowint fNbt1, fInt1;
    /// \brief standard ENDF interpolation parameter \cite ENDF Manual
    rowint fNbt2, fInt2;
    /// \brief standard ENDF interpolation parameter \cite ENDF Manual
    int fNr2, fNp2;
    /// \brief standard ENDF parameters for no .of regions and points for interpolation
    rowint fNbt3, fInt3;
    /// \brief standard ENDF interpolation parameter \cite ENDF Manual
    int fNr3, fNp3;
    /// \brief standard ENDF parameters for no .of regions and points for interpolation
    rowd fEnergyFile5, fEnergyPdfFile5, fEnergyCdfFile5;
    /// \brief energy, pdf and cdf for energy distribution for the file 5 with integrated angular distribution
    rowd fEin, fEneE, fCos, fCdf, fPdf, fCosIn, fCosInPdf, fCosInCdf;
    /// \brief Temp. energy, cosine, cdf and pdf for angle energy distribution
    matrixd2 fEin2D5, fFrac2D5;
    /// \brief Temp. energy and fractional contribution for energy distribution in file 5
    matrixd2 fEne2D, fFrac2D, fCdf2D, fPdf2D, fEin2D;
    /// \brief Temp. energy, cdf and pdf for energy distribution for the file 5
    matrixd3 fEne3D5, fCdf3D5, fPdf3D5;
    /// \brief Temp. energy, cdf and pdf for energy distribution for the file 5
    matrixd3 fEne3D, fCdf3D, fPdf3D;
    /// \brief Temp. energy, cdf and pdf for energy distribution for the file 4
    matrixd2 fCos2D;
    /// \brief Temp. variables for cosine angular distribution
    matrixd3 fCos3D;
    /// \brief Temp. variables for cosine angular distribution
    rowd fCosFile4, fCosPdfFile4, fCosCdfFile4;
    /// \brief cosine, pdf and cdf for angular distribution
    rowd fYmulti;
    /// \brief yield of the reaction in file 6
    rowint fMtLct, fMtLct4Cos;
    /// \brief temp. LCT numbers (flag for angular distribution in cm or lab system)
    rowint fZd, fAd;
    /// \brief z and A
    rowint fLaw;
    /// \brief law6 numbers for endf file 6
    rowint fMtNumbers6, fMtNumbers4;
    /// \brief temp. MT numbers
    rowint fMtNumNeutron, fMtNumPhoton, fMtNumCharge;
    /// \brief MT numbers either for neutron, photon or charge particles
    matrixd2 fCos2d, fCosinpdf2d, fCosincdf2d, fEin2d;
    /// \brief temp. energy, pdf and cdf for energy distribution
    matrixd3 fCos3d, fCosinpdf3d, fCosincdf3d;
    /// \brief temp. energy, pdf and cdf for energy distribution
    rowd fEoute, fCdfe, fPdfe;
    /// \brief temp. energy, pdf and cdf for energy distribution
    matrixd2 fEout2de, fPdf2de, fCdf2de;
    /// \brief temp. energy, pdf and cdf for energy distribution
    matrixd3 fEout3de, fPdf3de, fCdf3de;
    /// \brief temp. energy, pdf and cdf for energy distribution
    matrixd4 fEout4de, fPdf4de, fCdf4de;
    /// \brief energy, pdf and cdf for energy distribution
    rowint fMtNumbers4Cos, fMtNumbers5;
    /// \brief MT no. for cosine distribution from file 4, 6 and 5
    rowint fMtNumbers15;
    /// \brief MT no. for photon energy distribution from file 6 and 15
    matrixd2 fEin2D15, fFrac2D15;
    /// \brief Temp. energy and fractional contribution for energy distribution in file 5
    matrixd3 fEin3D15, fCdf3D15, fPdf3D15;
    /// \brief Temp. energy, cdf and pdf for energy distribution for the file 15
    
    NudyPhysics::TNudyEndfNuPh *fRecoNuPh;
    NudyPhysics::TNudyEndfPhYield *fRecoPhYield;
    NudyPhysics::TNudyEndfPhProd *fRecoPhProd;
    geant::RngWrapper fRng;
    #ifdef USE_ROOT
    ClassDef(TNudyEndfRecoPoint, 1) // class for an ENDF reconstruction
    #endif
  };
}
#endif
