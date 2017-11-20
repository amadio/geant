#include "SauterGavrilaPhotoElectricModel.h"

#include "PhysicalConstants.h"
#include "Material.h"
#include "Element.h"
#include "MaterialProperties.h"

#include "MaterialCuts.h"

#include "Spline.h"
#include "GLIntegral.h"
#include "AliasTable.h"
#include "XSectionsVector.h"

#include "PhysicsParameters.h"
#include "Gamma.h"
#include "Electron.h"
#include "LightTrack.h"
#include "PhysicsData.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

// from geantV
#include "Geant/Typedefs.h"
#include "GeantTaskData.h"


using namespace std;
namespace geantphysics {
    
    std::vector<double>*  SauterGavrilaPhotoElectricModel::fParamHigh[] = {nullptr};
    std::vector<double>*  SauterGavrilaPhotoElectricModel::fParamLow[] = {nullptr};
    
    int                   SauterGavrilaPhotoElectricModel::fNShells[] = {0};
    int                   SauterGavrilaPhotoElectricModel::fNShellsUsed[] = {0};
    
    Material*             SauterGavrilaPhotoElectricModel::fWater = nullptr;
    double                SauterGavrilaPhotoElectricModel::fWaterEnergyLimit = 0.0;
    
    XSectionsVector** SauterGavrilaPhotoElectricModel::fShellVector[]= {nullptr};
    XSectionsVector * SauterGavrilaPhotoElectricModel::fLECSVector[]= {nullptr};
    XSectionsVector * SauterGavrilaPhotoElectricModel::fCSVector[]= {nullptr};
    bool* SauterGavrilaPhotoElectricModel::fCrossSection = nullptr;
    bool* SauterGavrilaPhotoElectricModel::fCrossSectionLE = nullptr;

    SauterGavrilaPhotoElectricModel::SauterGavrilaPhotoElectricModel(const std::string &modelname, bool aliasActive)
    : EMModel(modelname){
        
        SetUseSamplingTables(aliasActive);
        fMinPrimEnergy           =  1.e-12*geant::eV;  // Minimum of the gamma kinetic energy grid, used to sample the photoelectron direction
        fMaxPrimEnergy           =  100*geant::MeV;    // Maximum of the gamma kinetic energy grid (after this threshold the e- is considered to follow the same direction as the incident gamma)
        
        fPrimEnLMin              = 0.;       // will be set in InitSamplingTables if needed
        fPrimEnILDelta           = 0.;       // will be set in InitSamplingTables if needed
        fSamplingPrimEnergies    = nullptr;  // will be set in InitSamplingTables if needed
        fLSamplingPrimEnergies   = nullptr;  // will be set in InitSamplingTables if needed
        
        fAliasData               = nullptr;  // will be set in InitSamplingTables if needed
        fAliasSampler            = nullptr;
        
        fCrossSection            = nullptr;
        fCrossSectionLE          = nullptr;
        
    }
    
    SauterGavrilaPhotoElectricModel::~SauterGavrilaPhotoElectricModel() {
        
        //CLEANING fParamHigh and fParamLow
        for(int i=0; i<gMaxSizeData; ++i) {
            delete fParamHigh[i];
            delete fParamLow[i];
            fParamHigh[i] = 0;
            fParamLow[i] = 0;
        }
        
        if (fSamplingPrimEnergies)
            delete [] fSamplingPrimEnergies;
        if (fLSamplingPrimEnergies)
            delete [] fLSamplingPrimEnergies;
        
        if (fAliasData) {
            for (int i=0; i<fNumSamplingPrimEnergies; ++i) {
                if (fAliasData[i]) {
                    delete [] fAliasData[i]->fXdata;
                    delete [] fAliasData[i]->fYdata;
                    delete [] fAliasData[i]->fAliasW;
                    delete [] fAliasData[i]->fAliasIndx;
                    delete fAliasData[i];
                }
            }
            delete [] fAliasData;
        }
        
        if (fAliasSampler)
            delete fAliasSampler;
        
    }
    
    void SauterGavrilaPhotoElectricModel::Initialize() {
        EMModel::Initialize();
        fSecondaryInternalCode    = Electron::Definition()->GetInternalCode();
        InitializeModel();
    }
    
    void SauterGavrilaPhotoElectricModel::InitializeModel() {
        
        //ALLOCATION fCrossSection
        if (fCrossSection) {
            delete [] fCrossSection;
            fCrossSection = nullptr;
        }
        
        //ALLOCATION fCrossSectionLE
        if (fCrossSectionLE) {
            delete [] fCrossSectionLE;
            fCrossSectionLE = nullptr;
            
        }

        fCrossSection   = new bool [gMaxSizeData];
        fCrossSectionLE = new bool [gMaxSizeData];
        
        for (int i=0; i<gMaxSizeData; ++i) {
            fCrossSection[i]   = false;
            fCrossSectionLE[i] = false;
        }
        
        fVerboseLevel=1;
        LoadData();
        if(GetUseSamplingTables())
            InitSamplingTables();
        
    }
    
    void SauterGavrilaPhotoElectricModel::SetVerboseLevel(int lev){
        fVerboseLevel=lev;
    }
    
    void SauterGavrilaPhotoElectricModel::LoadData()
    {
        int numMatCuts = MaterialCuts::GetTheMaterialCutsTable().size();
        // get list of active region
        std::vector<bool> isActiveInRegion = GetListActiveRegions();
        for (int i=0; i<numMatCuts; ++i) {
            const MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[i];
            // if this MaterialCuts belongs to a region where this model is active:
            if (isActiveInRegion[matCut->GetRegionIndex()]) {
                // get the list of elements
                const Vector_t<Element*> &theElements =  matCut->GetMaterial()->GetElementVector();
                int numElems = theElements.size();
                for (int j=0; j<numElems; ++j) {
                    double zet = theElements[j]->GetZ();
                    int elementIndx = std::lrint(zet);
                    //std::cout<<"\n***Calling ReadData on: "<<elementIndx<<" ....\n";
                    ReadData(elementIndx);
                }
            }
        }
    }
    
    void SauterGavrilaPhotoElectricModel::ReadData(int Z)
    {
        //bool debug= false;
        using geant::MeV;
        using geant::barn;
        
        if (fVerboseLevel > 1)
        {
            std::cout << "Calling ReadData() of SauterGavrilaPhotoElectricModel"
            << std::endl;
        }
        
        if(  (fCrossSection[Z]) && ( (fCrossSectionLE[Z] && Z>2) || (!fCrossSectionLE[Z] && Z<3)) )
        {//std::cout<<"Data loaded before!\n";
            return;}
        //else std::cout<<"Data not loaded before!\n";
        
        // get the path to the main physics data directory
        char *path = std::getenv("GEANT_PHYSICS_DATA");
        if (!path) {
            std::cerr<<"******   ERROR in SauterGavrilaPhotoElectricModel::ReadData() \n"
            <<"         GEANT_PHYSICS_DATA is not defined! Set the GEANT_PHYSICS_DATA\n"
            <<"         environment variable to the location of Geant data directory!\n"
            <<std::endl;
            exit(1);
        }
        
        std::ostringstream ost;
        ost << path << "/livermore/phot_epics2014/pe-cs-" << Z <<".dat";
        std::ifstream fin(ost.str().c_str());
        if( !fin.is_open()){
            std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost.str().c_str()
            << "> is not opened!" << std::endl;
            
            return;
        } else {
            
            fCrossSection[Z] =true;
            if(fVerboseLevel > 2) { std::cout << "File " << ost.str().c_str()
                << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;}
        
            fCSVector[Z]= new XSectionsVector;
            fin >> fCSVector[Z]->fEdgeMin >> fCSVector[Z]->fEdgeMax >> fCSVector[Z]->fNumberOfNodes;
            
            int siz=0;
            fin >> siz;
            if (fin.fail() || siz<=0)
            {
                std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost.str().c_str()
                << "> is not opened!" << std::endl;
            }
            
            fCSVector[Z]->fBinVector.reserve(siz);
            fCSVector[Z]->fDataVector.reserve(siz);
            
            
            double vBin, vData;
            
            for(int i = 0; i < siz ; i++)
            {
                vBin = 0.;
                vData= 0.;
                fin >> vBin >> vData;
                
                if (fin.fail())  {
                    std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost.str().c_str()
                    << "> is not opened!" << std::endl;
                }
                
                fCSVector[Z]->fBinVector.push_back(vBin*MeV);
                fCSVector[Z]->fDataVector.push_back(vData*barn);
            }
            
            // to remove any inconsistency
            fCSVector[Z]->fNumberOfNodes = siz;
            fCSVector[Z]->fEdgeMin = fCSVector[Z]->fBinVector[0];
            fCSVector[Z]->fEdgeMax = fCSVector[Z]->fBinVector[fCSVector[Z]->fNumberOfNodes-1];
            
            
            //Use spline interpolator for Cross-sections vector
            fCSVector[Z]->fSplineInt= new Spline((fCSVector[Z]->fBinVector.data()),(fCSVector[Z]->fDataVector.data()),fCSVector[Z]->fNumberOfNodes);
            fin.close();
        }
        
        // read high-energy fit parameters
        fParamHigh[Z] = new std::vector<double>;
        
        int n1 = 0;
        int n2 = 0;
        double x;
        std::ostringstream ost1;
        ost1 << path << "/livermore/phot_epics2014/pe-high-" << Z <<".dat";
        std::ifstream fin1(ost1.str().c_str());
        if( !fin1.is_open()) {
            std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost1.str().c_str()
            << "> is not opened!" << std::endl;
            return;
        } else {
            if(fVerboseLevel > 2) {
                std::cout << "File " << ost1.str().c_str()
                << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
            }
            fin1 >> n1;
            if(fin1.fail()) { return; }
            if(0 > n1 || n1 >= INT_MAX) { n1 = 0; }
            
            fin1 >> n2;
            if(fin1.fail()) { return; }
            if(0 > n2 || n2 >= INT_MAX) { n2 = 0; }
            
            fin1 >> x;
            if(fin1.fail()) { return; }
            
            fNShells[Z] = n1;
            fParamHigh[Z]->reserve(7*n1+1);
            fParamHigh[Z]->push_back(x*MeV);
            for(int i=0; i<n1; ++i) {
                for(int j=0; j<7; ++j) {
                    fin1 >> x;
                    if(0 == j) { x *= MeV; }
                    else       { x *= barn; }
                    fParamHigh[Z]->push_back(x);
                }
            }
            fin1.close();
        }
        
        //std::cout<<"pe-high parameterization for ["<<Z<<"], loaded\n";
        
        // read low-energy fit parameters
        fParamLow[Z] = new std::vector<double>;
        int n1_low = 0;
        int n2_low = 0;
        double x_low;
        std::ostringstream ost1_low;
        ost1_low << path << "/livermore/phot_epics2014/pe-low-" << Z <<".dat";
        std::ifstream fin1_low(ost1_low.str().c_str());
        if( !fin1_low.is_open()) {
            std::cerr << "SauterGavrilaPhotoElectricModel data file <" << ost1_low.str().c_str()
            << "> is not opened!" << std::endl;
            return;
        } else {
            if(fVerboseLevel > 2) {
                std::cout << "File " << ost1_low.str().c_str()
                << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
            }
            fin1_low >> n1_low;
            if(fin1_low.fail()) { return; }
            if(0 > n1_low || n1_low >= INT_MAX) { n1_low = 0; }
            
            fin1_low >> n2_low;
            if(fin1_low.fail()) { return; }
            if(0 > n2_low || n2_low >= INT_MAX) { n2_low = 0; }
            
            fin1_low >> x_low;
            if(fin1_low.fail()) { return; }
            
            fNShells[Z] = n1_low;
            fParamLow[Z]->reserve(7*n1_low+1);
            fParamLow[Z]->push_back(x_low*MeV);
            for(int i=0; i<n1_low; ++i) {
                for(int j=0; j<7; ++j) {
                    fin1_low >> x_low;
                    if(0 == j) { x_low *= MeV; }
                    else       { x_low *= barn; }
                    fParamLow[Z]->push_back(x_low);
                }
            }
            fin1_low.close();
        }
        
        //std::cout<<"pe-low parameterization for ["<<Z<<"], loaded\n";
        
        
        // there is a possibility to use only main shells
        if(gNShellLimit < n2)
        {
            n2 = gNShellLimit;
        }
        fNShellsUsed[Z]=n2;
        fShellVector[Z] = new XSectionsVector*[n2];
        for(int i=0; i<n2; i++)
            fShellVector[Z][i]=new XSectionsVector;
        
        for (int i=0 ; i<n2 ; i++)
        {
            fShellVector[Z][i]->fDataVector.clear();
            fShellVector[Z][i]->fBinVector.clear();
            fShellVector[Z][i]->fNumberOfNodes=0;
        }
        
        fNShellsUsed[Z] = n2;
        
        //If more than one shell is used -> Read sub-shells cross section data
        if(1 < n2) {
            std::ostringstream ost2;
            ost2 << path << "/livermore/phot_epics2014/pe-ss-cs-" << Z <<".dat";
            std::ifstream fin2(ost2.str().c_str());
            if( !fin2.is_open()) {
                std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost2.str().c_str()
                << "> is not opened!" << std::endl;
                return;
            } else {
                if(fVerboseLevel > 2) {
                    std::cout << "File " << ost2.str().c_str()
                    << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
                }
                
                int n3, n4;
                double y;
                
                for(int i=0; i<n2; ++i)
                {
                    fin2 >> x >> y >> n3 >> n4;
                    
                    fShellVector[Z][i]->fBinVector.clear();
                    fShellVector[Z][i]->fDataVector.clear();
                    fShellVector[Z][i]->fBinVector.reserve(n3);
                    fShellVector[Z][i]->fDataVector.reserve(n3);
                    fShellVector[Z][i]->fEdgeMin=x*MeV;
                    fShellVector[Z][i]->fEdgeMax=y*MeV;
                    if(fVerboseLevel>3) std::cout<<"fShellVector["<<Z<<"]["<<i<<"]->fEdgeMin: "<<fShellVector[Z][i]->fEdgeMin<<"\t fShellVector["<<Z<<"]["<<i<<"]->fEdgeMax: "<<fShellVector[Z][i]->fEdgeMax<<std::endl;
                    
                    for(int j=0; j<n3; ++j)
                    {
                        fin2 >> x >> y;
                        fShellVector[Z][i]->fBinVector.push_back(x*MeV);
                        fShellVector[Z][i]->fDataVector.push_back(y*barn);
                        fShellVector[Z][i]->fNumberOfNodes++;
                    }
                    fShellVector[Z][i]->fCompID=n4;
                }
                
                fin2.close();
            }
        }
        //std::cout<<"pe-ss-cs- cross sections for ["<<Z<<"], loaded\n";
        
        // no spline for photoeffect total x-section below K-shell
        if(1 < fNShells[Z]) {
            
            std::ostringstream ost3;
            ost3 << path << "/livermore/phot_epics2014/pe-le-cs-" << Z <<".dat";
            std::ifstream fin3(ost3.str().c_str());
            if( !fin3.is_open()) {
                std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost3.str().c_str()
                << "> is not opened!" << std::endl;
                return;
            } else {
                if(fVerboseLevel > 2) {
                    std::cout << "File " << ost3.str().c_str()
                    << " is opened by SauterGavrilaPhotoElectricModel" << std::endl;
                }
                
                fCrossSectionLE[Z] = true;
                fLECSVector[Z]= new XSectionsVector;
                
                fin3 >> fLECSVector[Z]->fEdgeMin >> fLECSVector[Z]->fEdgeMax >> fLECSVector[Z]->fNumberOfNodes;
                int siz=0;
                fin3 >> siz;
                if (fin3.fail() || siz<=0)
                {
                    std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost3.str().c_str()
                    << "> is not opened!" << std::endl;
                }
                
                fLECSVector[Z]->fBinVector.reserve(siz);
                fLECSVector[Z]->fDataVector.reserve(siz);
                
                double vBin, vData;
                
                for(int i = 0; i < siz ; i++)
                {
                    vBin = 0.;
                    vData= 0.;
                    fin3 >> vBin >> vData;
                    
                    //Scale vector
                    vBin  *= MeV;
                    vData *= barn;
                    
                    if (fin3.fail())  {
                        std::cerr<< "SauterGavrilaPhotoElectricModel data file <" << ost3.str().c_str()
                        << "> is not opened!" << std::endl;
                    }
                    
                    fLECSVector[Z]->fBinVector.push_back(vBin);
                    fLECSVector[Z]->fDataVector.push_back(vData);
                    
                }
                
                // to remove any inconsistency
                fLECSVector[Z]->fNumberOfNodes = siz;
                fLECSVector[Z]->fEdgeMin = fLECSVector[Z]->fBinVector[0];
                fLECSVector[Z]->fEdgeMax = fLECSVector[Z]->fBinVector[fLECSVector[Z]->fNumberOfNodes-1];
                
                fin3.close();
            }
            //std::cout<<"pe-le-cs- cross sections for ["<<Z<<"], loaded\n";
        }
    }
    
    
    //____________________
    //NB: cosTheta is supposed to contain the dirZ of the incoming photon
    void SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirection_Rejection(double gammaEnIn, double &cosTheta, Geant::GeantTaskData *td) const{
        
        //1) initialize energy-dependent variables
        // Variable naming according to Eq. (2.24) of Penelope Manual
        // (pag. 44)
        double gamma = 1.0 + gammaEnIn/geant::kElectronMassC2;
        double gamma2 = gamma*gamma;
        double beta = std::sqrt((gamma2-1.0)/gamma2);
            
        // ac corresponds to "A" of Eq. (2.31)
        //
        double ac = (1.0/beta) - 1.0;
        double a1 = 0.5*beta*gamma*(gamma-1.0)*(gamma-2.0);
        double a2 = ac + 2.0;
        // gtmax = maximum of the rejection function according to Eq. (2.28), obtained for tsam=0
        double gtmax = 2.0*(a1 + 1.0/ac);
            
        double tsam = 0;
        double gtr = 0;
            
        //2) sampling. Eq. (2.31) of Penelope Manual
        // tsam = 1-std::cos(theta)
        // gtr = rejection function according to Eq. (2.28)
        double rndArray[2];
        do {
          td->fRndm->uniform_array(2, rndArray);
          tsam = 2.0*ac * (2.0*rndArray[0] + a2*std::sqrt(rndArray[0])) / (a2*a2 - 4.0*rndArray[0]);
          gtr = (2.0 - tsam) * (a1 + 1.0/(ac+tsam));
        } while (rndArray[1]*gtmax > gtr);
        
        cosTheta = 1.0-tsam;
        return;
    }
    
    
    //____________________
    double SauterGavrilaPhotoElectricModel::SamplePhotoElectronDirection_Alias(double primekin, double r1, double r2, double r3) const{
        
        // determine primary energy lower grid point
        if (primekin > 100*geant::MeV) {
            return 1.;
        } else
        {
            //std::cout<<"::::SamplePhotoElectronDirection_Alias::::\n";
            double lGammaEnergy  = std::log(primekin);
            int gammaEnergyIndx  = (int) ((lGammaEnergy-fPrimEnLMin)*fPrimEnILDelta);
            //
            if (gammaEnergyIndx>=fNumSamplingPrimEnergies-1)
                gammaEnergyIndx = fNumSamplingPrimEnergies-2;
            //
            double pLowerGammaEner = (fLSamplingPrimEnergies[gammaEnergyIndx+1]-lGammaEnergy)*fPrimEnILDelta;
            if (r1>pLowerGammaEner) {
                ++gammaEnergyIndx;
            }
            // sample the outgoing electron cosTheta
            double ecosTheta = fAliasSampler->SampleLinear(fAliasData[gammaEnergyIndx]->fXdata, fAliasData[gammaEnergyIndx]->fYdata,fAliasData[gammaEnergyIndx]->fAliasW, fAliasData[gammaEnergyIndx]->fAliasIndx,fAliasData[gammaEnergyIndx]->fNumdata,r2,r3);
            
            //This have to be seen (Transformation)
            //double xsi=fAliasSampler->SampleLinear(fAliasData[gammaEnergyIndx]->fXdata, fAliasData[gammaEnergyIndx]->fYdata,fAliasData[gammaEnergyIndx]->fAliasW, fAliasData[gammaEnergyIndx]->fAliasIndx,fAliasData[gammaEnergyIndx]->fNumdata,r2,r3);
            
            //double ecosTheta= std::exp(xsi)-2;
            return ecosTheta;
        }
        
    }
    
    double SauterGavrilaPhotoElectricModel::ComputeXSectionPerAtom(const Element *elem, const MaterialCuts*, double kinenergy,
                                                                   const Particle*) {
        
        double xsec  = 0.0;
        if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
            return xsec;
        }
        // compute the parameterized atomic cross section: depends only on target Z and gamma energy.
        xsec = ComputeXSectionPerAtom(elem->GetZ(), kinenergy);
        return xsec;
    }
    
    
    double SauterGavrilaPhotoElectricModel::ComputeXSectionPerAtom(double zeta, double energy) const
    {
        
        int verboseLevel = 0;
        using geant::keV;
        using geant::MeV;
        using geant::barn;
        
        
        if (verboseLevel > 3) {
            std::cout << "G4LivermorePhotoElectricModel_new::ComputeCrossSectionPerAtom():"
            << " Z= " << zeta << "  R(keV)= " << energy/keV << std::endl;
        }
        double cs = 0.0;
        int Z = std::lrint(zeta);
        if(Z < 1 || Z >= gMaxSizeData) { return cs; }
        
        
        //EXTRA CHECK --
        //if( !fCrossSection[Z] && !fCrossSectionLE[Z] )
        //{
        //    ReadData(Z);
        //    if(!fCrossSectionLE[Z] && !fCrossSection[Z]) { return cs; }
        //}
        
        int idx = fNShells[Z]*7 - 5; //7: rows in the parameterization file; 5: number of parameters
        
        if (energy < (*(fParamHigh[Z]))[idx-1]) { energy = (*(fParamHigh[Z]))[idx-1];}
        
        double x1 = (MeV)/(energy);
        double x2 = x1*x1;
        double x3 = x2*x1;
        
        // (*) High energy parameterisation
        if(energy >= (*(fParamHigh[Z]))[0]) {
            double x4 = x2*x2;
            double x5 = x4*x1;
            
            cs = x1*((*(fParamHigh[Z]))[idx] + x1*(*(fParamHigh[Z]))[idx+1]
                     + x2*(*(fParamHigh[Z]))[idx+2] + x3*(*(fParamHigh[Z]))[idx+3]
                     + x4*(*(fParamHigh[Z]))[idx+4]+ x5*(*(fParamHigh[Z]))[idx+5]);
            
            
        }
        // (**) Low energy parameterisation
        else if(energy >= (*(fParamLow[Z]))[0]) {
            double x4 = x2*x2;
            double x5 = x4*x1;//this variable usage can probably be optimized
            cs = x1*((*(fParamLow[Z]))[idx] + x1*(*(fParamLow[Z]))[idx+1]
                     + x2*(*(fParamLow[Z]))[idx+2] + x3*(*(fParamLow[Z]))[idx+3]
                     + x4*(*(fParamLow[Z]))[idx+4]+ x5*(*(fParamLow[Z]))[idx+5]);
        }
        
        // (***) Tabulated values above k-shell ionization energy
        else if(energy >= (*(fParamHigh[Z]))[1]) {
            //this is an extra-check - if energy<fCSVector[Z]->edgeMin it should go to the next else (****)
            if(energy<fCSVector[Z]->fEdgeMin)
                cs=x3*fCSVector[Z]->fDataVector[0];
            else
            {
                size_t idx=0;
                idx=fCSVector[Z]->FindCSBinLocation(energy, idx);
                cs=x3*fCSVector[Z]->fSplineInt->GetValueAt(energy, idx);
            }
            
        }
        
        //(****) Tabulated values below k-shell ionization energy
        else
        {
            //this check is needed to have a constant cross-section(cs) value for energies below the lowest cs point.
            //in this way the gamma is always absorbed - instead of giving zero cs -
            if(energy<fLECSVector[Z]->fEdgeMin)
                cs=x3*fLECSVector[Z]->fDataVector[0];
            else
            {
                size_t idx=0;
                cs=x3*fLECSVector[Z]->GetValue(energy, idx);
            }
        }
        if (verboseLevel > 1) {
            std::cout << "LivermorePhotoElectricModel: E(keV)= " << energy/keV
            << " Z= " << Z << " cross(barn)= " << cs/barn << std::endl;
        }
        
        return cs;
        
    }
    
    double SauterGavrilaPhotoElectricModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                                       const Particle* /*particle*/) {
        
        double xsec = 0.0;
        if (kinenergy<GetLowEnergyUsageLimit() || kinenergy>GetHighEnergyUsageLimit()) {
            return xsec;
        }
        // compute the macroscopic cross section as the sum of the atomic cross sections weighted by the number of atoms in
        // in unit volume.
        const Material *mat =  matcut->GetMaterial();
        double       egamma = kinenergy;
        // we need the element composition of this material
        const Vector_t<Element*> &theElements    = mat->GetElementVector();
        const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
        int   numElems = theElements.size();
        for (int iel=0; iel<numElems; ++iel) {
            xsec += theAtomicNumDensityVector[iel]*ComputeXSectionPerAtom(theElements[iel]->GetZ(), egamma);
        }
        return xsec;
    }
    
    
    size_t SauterGavrilaPhotoElectricModel::SampleTargetElementIndex (const MaterialCuts *matCut, double gammaekin0, Geant::GeantTaskData *td) const
    {
        size_t index =0;
        std::vector<double> mxsec(20,0.);
        
        const Material *mat =  matCut->GetMaterial();
        const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
        
        const Vector_t<Element*> &theElements = mat->GetElementVector();
        size_t num    = matCut->GetMaterial()->GetNumberOfElements();
        if (num>mxsec.size()) {mxsec.resize(num,0.);}
        double cum =0.;
        for (size_t i=0; i<num; ++i) {
            //double xx=theAtomicNumDensityVector[i]* ComputeXSectionPerAtom(theElements[i]->GetZ(), gammaekin0);
            //calculate the Macroscopic Cross section
            cum += theAtomicNumDensityVector[i]* ComputeXSectionPerAtom(theElements[i]->GetZ(), gammaekin0);;
            //store directly the cumulative
            mxsec[i] = cum;
        }
        double rnd=cum * td->fRndm->uniform();
        //double cumxsec=mxsec[0];
        for(; index<num-1 && rnd>mxsec[index]; ++index) {/*cumxsec += mxsec[index+1];*/}
        return index;

//        int index=0;
//        //retrieve the elements vector
//        const Vector_t<Element*> &theElements = matCut->GetMaterial()->GetElementVector();
//        //retrieve the number of elements in the material
//        int num    = matCut->GetMaterial()->GetNumberOfElements();
//        if (num > 1)
//        {
//            double macxsec=ComputeMacroscopicXSection(matCut,gammaekin0, Gamma::Definition());
//            double rnd=macxsec * td->fRndm->uniform();
//            
//            const Material *mat =  matCut->GetMaterial();
//            const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
//            double cumxsec=0.;
//            for(; index<num-1; ++index)
//            {
//                double xsec = theAtomicNumDensityVector[index]* ComputeXSectionPerAtom(theElements[index]->GetZ(), gammaekin0);
//                cumxsec += xsec;
//                if (rnd <= cumxsec)
//                {
//                    break;
//                }
//            }
//        }
//        return index;
  
  }
    
    void SauterGavrilaPhotoElectricModel::TestSampleTargetElementIndex(const MaterialCuts *matcut, double energy, Geant::GeantTaskData *td) const{
        
        std::cout<<"testSampleTargetElementIndex\n";
        size_t index=0;
        double sum=0;
        //retrieve the elements vector
        const Vector_t<Element*> &theElements = matcut->GetMaterial()->GetElementVector();
        //retrieve the number of elements in the material
        int num    = matcut->GetMaterial()->GetNumberOfElements();
        double xsec[num];
        double xsecSampled[num];
        
        const Material *mat =  matcut->GetMaterial();
        const double* theAtomicNumDensityVector = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
        for (int i=0; i<num; i++)
        {
            xsec[i] = theAtomicNumDensityVector[i]* ComputeXSectionPerAtom(theElements[i]->GetZ(), energy);
            xsecSampled[i]=0.;
            sum+=xsec[i];
        }
        for (int i=0; i<1000000000; i++){
            index= SampleTargetElementIndex(matcut, energy, td);
            xsecSampled[index]++;
        }
        for (int i=0; i<num; i++)
        {
            xsec[i]/=sum;
            xsecSampled[i]/=1000000000;
        }
        
        char filename[521];
        sprintf(filename,"SampleTargetElementIndexTest_%s", (matcut->GetMaterial()->GetName()).c_str() );
        FILE *f     = fopen(filename,"w");
        
        for (int i=0; i<num; ++i) {
            fprintf(f,"%d\t%.8g\t%.8g\n",i,xsec[i],xsecSampled[i]);
        }
        fclose(f);
        
    }
    
    
    int SauterGavrilaPhotoElectricModel::SampleSecondaries(LightTrack &track,
                                                           Geant::GeantTaskData *td){
        
        using geant::MeV;
        double gammaekin0          = track.GetKinE();
        // check if kinetic energy is below fLowEnergyUsageLimit and do nothing if yes;
        // check if kinetic energy is above fHighEnergyUsageLimit and do nothing if yes;
        
        if (gammaekin0<GetLowEnergyUsageLimit() || gammaekin0>GetHighEnergyUsageLimit())
        {
            return 0; //numSecondaries is zero since the interaction is not happening
        }
        
        //interaction is possible so sample target element
        MaterialCuts *matCut = MaterialCuts::GetTheMaterialCutsTable()[track.GetMaterialCutCoupleIndex()];
        const Vector_t<Element*> &theElements = matCut->GetMaterial()->GetElementVector();
        
        size_t targetElemIndx = 0;
        if (theElements.size()>1) {
            //uncomment the following line to test SampleTargetElementIndex
            //testSampleTargetElementIndex(matCut, gammaekin0, td );
            targetElemIndx = SampleTargetElementIndex(matCut, gammaekin0, td);
        }
        double  zeta  = theElements[targetElemIndx]->GetZ();
        int     Z = std::lrint(zeta);
       
        // if element was not initialised, gamma should be absorbed
        if(!fCrossSectionLE[Z] && !fCrossSection[Z]) {
            track.SetEnergyDeposit(gammaekin0);
            //std::cout<<"Model not initialized, Exiting!\n";
            return 0;
        }

        //SAMPLING OF THE SHELL
        size_t shellIdx = 0;
        size_t nn = fNShellsUsed[Z];
    
        if(nn > 1)
        {
            // sample gamma energy
            //double *rndArray = td->fDblArray;
            double rand=td->fRndm->uniform();
            //td->fRndm->uniform_array(1, rndArray);
            
            // (*) High energy parameterisation
            if(gammaekin0 >= (*(fParamHigh[Z]))[0])
            {
                double x1 = (MeV)/gammaekin0;
                double x2 = x1*x1;
                double x3 = x2*x1;
                double x4 = x3*x1;
                double x5 = x4*x1;
                int idx   = nn*7 - 5;
                // when do sampling common factors are not taken into account
                // so cross section is not real
                double cs0 = rand*(     (*(fParamHigh[Z]))[idx]
                                          + x1*(*(fParamHigh[Z]))[idx+1]
                                          + x2*(*(fParamHigh[Z]))[idx+2]
                                          + x3*(*(fParamHigh[Z]))[idx+3]
                                          + x4*(*(fParamHigh[Z]))[idx+4]
                                          + x5*(*(fParamHigh[Z]))[idx+5]);
                
                for(shellIdx=0; shellIdx<nn; ++shellIdx)
                {
                    idx = shellIdx*7 + 2;
                    if(gammaekin0 > (*(fParamHigh[Z]))[idx-1])
                    {
                        double cs =
                        (*(fParamHigh[Z]))[idx]
                        + x1*(*(fParamHigh[Z]))[idx+1]
                        + x2*(*(fParamHigh[Z]))[idx+2]
                        + x3*(*(fParamHigh[Z]))[idx+3]
                        + x4*(*(fParamHigh[Z]))[idx+4]
                        + x5*(*(fParamHigh[Z]))[idx+5];
                        if(cs >= cs0) {break;}
                    }
                }
                if(shellIdx >= nn) { shellIdx = nn-1; }
        
            }
            // (**) Low energy parameterisation
            else if(gammaekin0 >= (*(fParamLow[Z]))[0])
            {
                double x1 = MeV/gammaekin0;
                double x2 = x1*x1;
                double x3 = x2*x1;
                double x4 = x3*x1;
                double x5 = x4*x1;
                int idx   = nn*7 - 5;
                double cs0 = rand*((*(fParamLow[Z]))[idx]
                                          + x1*(*(fParamLow[Z]))[idx+1]
                                          + x2*(*(fParamLow[Z]))[idx+2]
                                          + x3*(*(fParamLow[Z]))[idx+3]
                                          + x4*(*(fParamLow[Z]))[idx+4]
                                          + x5*(*(fParamLow[Z]))[idx+5]);
                for(shellIdx=0; shellIdx<nn; ++shellIdx)
                {
                    idx = shellIdx*7 + 2;
                    if(gammaekin0 > (*(fParamHigh[Z]))[idx-1])
                    {
                        double cs = (*(fParamLow[Z]))[idx] + x1*(*(fParamLow[Z]))[idx+1]
                        + x2*(*(fParamLow[Z]))[idx+2] + x3*(*(fParamLow[Z]))[idx+3]
                        + x4*(*(fParamLow[Z]))[idx+4]+ x5*(*(fParamLow[Z]))[idx+5];
                        if(cs >= cs0) { break; }
                    }
                }
                if(shellIdx >= nn) {shellIdx = nn-1;}
            }
            else
            {
                double cs = rand;
                size_t idx= 0;
                // (***) Tabulated values above k-shell ionization energy
                if(gammaekin0 >= (*(fParamHigh[Z]))[1]) {
                    //this is an extra-check - if energy<fCSVector[Z]->edgeMin it should go to the next else (****)
                    if(gammaekin0<fCSVector[Z]->fEdgeMin)
                        cs*=fCSVector[Z]->fDataVector[0];
                    else
                    {
                        idx=fCSVector[Z]->FindCSBinLocation(gammaekin0, idx);
                        cs*=fCSVector[Z]->fSplineInt->GetValueAt(gammaekin0, idx);
                    }
                }
                //(****) Tabulated values below k-shell ionization energy
                else
                {
                    //this check is needed to have a constant cross-section(cs) value for energies below the lowest cs point.
                    //in this way the gamma is always absorbed - instead of giving zero cs -
                    if(gammaekin0<fLECSVector[Z]->fEdgeMin)
                        cs*=fLECSVector[Z]->fDataVector[0];
                    else
                        //cs*=fLECSVector[Z]->GetValueAt(gammaekin0);
                        cs*=fLECSVector[Z]->GetValue(gammaekin0, idx);
                }
                //size_t j=0;
                for(size_t j=0; j<nn; ++j)
                {
                    
                    shellIdx=(size_t)fShellVector[Z][j]->fCompID;
                    if(gammaekin0 > (*(fParamLow[Z]))[7*shellIdx+1]) {
                        size_t idx=0;
                        //cs-=fShellVector[Z][j]->GetValueAt(gammaekin0);
                        cs-=fShellVector[Z][j]->GetValue(gammaekin0, idx);
                    }
                    
                    if(cs <= 0.0 || j+1 == nn)
                        break;
                }
            }
        }
        
        // END: SAMPLING OF THE SHELL
        
        //Retrieving ionized shell bindingEnergy
        double bindingEnergy = (*(fParamHigh[Z]))[shellIdx*7 + 1];
        
        if (gammaekin0 < bindingEnergy) {
          track.SetEnergyDeposit(gammaekin0);
          return 0; //numSecondaries
        }
        
        //since edep is equal to bindingenergy I get rid of it
        //double edep = bindingEnergy;
        
        //  deexcitation is MISSING for now
        /*
         DEEXCITATION
         */
        //
        
        
        // Create the secondary particle: the photoelectron
        double elecKineEnergy = gammaekin0 - bindingEnergy;
        double cosTheta = track.GetDirZ();
        double sinTheta = 0.0;
        double phi = 0.0;
        double eDirX1;
        double eDirY1;
        double eDirZ1;
        
        if (gammaekin0 <= 100*geant::MeV) {
            if (!GetUseSamplingTables()) {
              SamplePhotoElectronDirection_Rejection(gammaekin0, cosTheta, td);
            } else {
              double *rndArray = td->fDblArray;
              td->fRndm->uniform_array(3, rndArray);
              cosTheta=SamplePhotoElectronDirection_Alias(gammaekin0, rndArray[0], rndArray[1], rndArray[2]);
            }
            sinTheta = std::sqrt((1 - cosTheta)*(1 + cosTheta));
            double rnd = td->fRndm->uniform();
            phi     = geant::kTwoPi * rnd;
            
            // new photoelectron direction in the scattering frame
            eDirX1  = sinTheta*std::cos(phi);
            eDirY1  = sinTheta*std::sin(phi);
            eDirZ1  = cosTheta;
 
            // rotate new photoelectron direction to the lab frame:
            RotateToLabFrame(eDirX1, eDirY1, eDirZ1,track.GetDirX(), track.GetDirY(), track.GetDirZ());
            
        } else {
          eDirX1  = track.GetDirX();
          eDirY1  = track.GetDirY();
          eDirZ1  = track.GetDirZ();
        }
    
        // create the secondary particle i.e. the photoelectron
        //numSecondaries = 1;
        
        // current capacity of secondary track container
        int curSizeOfSecList = td->fPhysicsData->GetSizeListOfSecondaries();
        // currently used secondary tracks in the container
        int curNumUsedSecs   = td->fPhysicsData->GetNumUsedSecondaries();
        if (curSizeOfSecList-curNumUsedSecs<1) {
            td->fPhysicsData->SetSizeListOfSecondaries(2*curSizeOfSecList);
        }
        int secIndx = curNumUsedSecs;
        curNumUsedSecs +=1;
        td->fPhysicsData->SetNumUsedSecondaries(curNumUsedSecs);
        
        std::vector<LightTrack>& sectracks = td->fPhysicsData->GetListOfSecondaries();
        
        // this is known since it is a secondary track
        //  sectracks[secIndx].SetTrackStatus(LTrackStatus::kNew); // to kew
        
        sectracks[secIndx].SetDirX(eDirX1);
        sectracks[secIndx].SetDirY(eDirY1);
        sectracks[secIndx].SetDirZ(eDirZ1);
        sectracks[secIndx].SetKinE(elecKineEnergy);
        sectracks[secIndx].SetGVcode(fSecondaryInternalCode);  // electron GV code
        sectracks[secIndx].SetMass(geant::kElectronMassC2);
        sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); // parent GeantTrack index
        
        
        /*if(fabs(gammaekin0 - elecKineEnergy - esec - edep) > geant::eV) {
         std::cout << "### SauterGavrilaPhotoElectricModel dE(eV)= "
         << (gammaekin0 - elecKineEnergy - esec - edep)/geant::eV
         << "  shell= " << shellIdx
         << "  E(keV)= " << gammaekin0/geant::keV
         << "  Ebind(keV)= " << bindingEnergy/geant::keV
         << "  Ee(keV)= " << elecKineEnergy/geant::keV
         << "  Esec(keV)= " << esec/geant::keV
         << "  Edep(keV)= " << edep/geant::keV
         << std::endl;
         }*/
        
        
        //always kill primary photon
        track.SetTrackStatus(LTrackStatus::kKill);
        track.SetKinE(0.0);
        //edep is = bindingEnergy
        //if(edep > 0.0) {
        if (bindingEnergy > 0.0) {
          track.SetEnergyDeposit(bindingEnergy);
            //track.SetEnergyDeposit(edep);
        }
        // return with number of secondaries i.e. 1 photoelectron
        return 1;
    }
    
    
    void SauterGavrilaPhotoElectricModel::InitSamplingTables() {
        
        // set number of primary gamma energy grid points
        // keep the prev. value of primary energy grid points.
        int oldNumGridPoints = fNumSamplingPrimEnergies;
        fNumSamplingPrimEnergies = fNumSamplingPrimEnergiesPerDecade*std::lrint(std::log10(fMaxPrimEnergy/fMinPrimEnergy))+1;
        if (fNumSamplingPrimEnergies<2) {
            fNumSamplingPrimEnergies = 2;
        }
        
        // set up the initial gamma energy grid
        if (fSamplingPrimEnergies) {
            delete [] fSamplingPrimEnergies;
            delete [] fLSamplingPrimEnergies;
            fSamplingPrimEnergies  = nullptr;
            fLSamplingPrimEnergies = nullptr;
        }
        fSamplingPrimEnergies  = new double[fNumSamplingPrimEnergies];
        fLSamplingPrimEnergies = new double[fNumSamplingPrimEnergies];
        fPrimEnLMin    = std::log(fMinPrimEnergy);
        double delta   = std::log(fMaxPrimEnergy/fMinPrimEnergy)/(fNumSamplingPrimEnergies-1.0);
        fPrimEnILDelta = 1.0/delta;
        fSamplingPrimEnergies[0]  = fMinPrimEnergy;
        fLSamplingPrimEnergies[0] = fPrimEnLMin;
        fSamplingPrimEnergies[fNumSamplingPrimEnergies-1]  = fMaxPrimEnergy;
        fLSamplingPrimEnergies[fNumSamplingPrimEnergies-1] = std::log(fMaxPrimEnergy);
        for (int i=1; i<fNumSamplingPrimEnergies-1; ++i) {
            fLSamplingPrimEnergies[i] = fPrimEnLMin+i*delta;
            fSamplingPrimEnergies[i]  = std::exp(fPrimEnLMin+i*delta);
        }
        //
        // build the sampling tables at each primary gamma energy grid point.
        //
        // prepare the array that stores pointers to sampling data structures
        if (fAliasData) {
            for (int i=0; i<oldNumGridPoints; ++i) {
                if (fAliasData[i]) {
                    delete [] fAliasData[i]->fXdata;
                    delete [] fAliasData[i]->fYdata;
                    delete [] fAliasData[i]->fAliasW;
                    delete [] fAliasData[i]->fAliasIndx;
                    delete fAliasData[i];
                }
            }
            delete [] fAliasData;
        }
        // create new fAliasData array
        fAliasData = new LinAlias*[fNumSamplingPrimEnergies];
        for (int i=0; i<fNumSamplingPrimEnergies; ++i) {
            fAliasData[i] = nullptr;
        }
        // create one sampling data structure at each primary gamma energy grid point:
        // -first create an AliasTable object
        if (fAliasSampler) {
            delete fAliasSampler;
        }
        // -the prepare each table one-by-one
        for (int i=0; i<fNumSamplingPrimEnergies; ++i) {
            double egamma = fSamplingPrimEnergies[i];
            double tau  = egamma/geant::kElectronMassC2;
            BuildOneLinAlias(i,tau);
        }
    }
    
    //This method is calculating the differential cross section in the transformed variable xsi
    double SauterGavrilaPhotoElectricModel::CalculateDiffCrossSectionLog(double tau, double xsi) const
    {
        
        // Based on Geant4 : G4SauterGavrilaAngularDistribution
        // SauterGavrila approximation for K-shell, correct to the first \alphaZ order
        // input  : energy0  (incoming photon energy)
        // input  : cosTheta (cons(theta) of photo-electron)
        // output : dsigma   (differential cross section, K-shell only)
        
        //double tau = energy0 / geant::kElectronMassC2;
        
        //std::cout<<"CalculateDiffCrossSectionLog. tau: "<<tau<<" and xsi: "<<xsi<<std::endl;
        double cosTheta= std::exp(xsi) - 2;
        
        //gamma and beta: Lorentz factors of the photoelectron
        double gamma = tau + 1.0;
        double beta = std::sqrt(tau * (tau + 2.0)) / gamma;
        
        double z = 1 - beta * cosTheta;
        double z2 = z * z;
        double z4 = z2 * z2;
        double y = 1 - cosTheta * cosTheta; //sen^2(theta)
        
        double dsigmadcostheta = (y / z4) * (1 + 0.5 * gamma * (tau) * (gamma - 2) * z);
        double dsigmadxsi=dsigmadcostheta*(std::exp(xsi));
        //std::cout<<"dsigmadcostheta: "<<dsigmadcostheta<<" and dsigmadxsi: "<<dsigmadxsi<<std::endl;
        return dsigmadxsi;
        
    }
    
    
    double SauterGavrilaPhotoElectricModel::CalculateDiffCrossSection(double tau, double cosTheta)const
    {
        
        // Based on Geant4 : G4SauterGavrilaAngularDistribution
        // SauterGavrila approximation for K-shell, correct to the first \alphaZ order
        // input  : energy0  (incoming photon energy)
        // input  : cosTheta (cons(theta) of photo-electron)
        // output : dsigma   (differential cross section, K-shell only)
        
        //double tau = energy0 / geant::kElectronMassC2;
        
        //gamma and beta: Lorentz factors of the photoelectron
        double gamma = tau + 1.0;
        double beta = std::sqrt(tau * (tau + 2.0)) / gamma;
        
        double z = 1 - beta * cosTheta;
        double z2 = z * z;
        double z4 = z2 * z2;
        double y = 1 - cosTheta * cosTheta; //sen^2(theta)
        
        double dsigma = (y / z4) * (1 + 0.5 * gamma * (tau) * (gamma - 2) * z);
        return dsigma;
        
    }
    
    
    int SauterGavrilaPhotoElectricModel::PrepareLinAlias(double tau, std::vector<double> & x, std::vector<double> & y){
        
        int numpoints=40;
        int curNumData = 5;
        double maxErrorThreshold = gsingleTableErrorThreshold;
        
        x.resize(numpoints);
        y.resize(numpoints);
        
        //cosTheta variable between [-1, 1]
        //start with 5 points: thetaMin, thetaMin+0.1, 0, thetaMax-0.1, thetaMax
        double thetaMin=-1.;
        double thetaMax=1.;
        x[0] = thetaMin;
        x[1] = thetaMin+0.1;
        x[2] = 0.;
        x[3] = thetaMax-0.1;
        x[4] = thetaMax;
        y[0] = CalculateDiffCrossSection(tau, x[0]);
        y[1] = CalculateDiffCrossSection(tau, x[1]);
        y[2] = CalculateDiffCrossSection(tau, x[2]);
        y[3] = CalculateDiffCrossSection(tau, x[3]);
        y[4] = CalculateDiffCrossSection(tau, x[4]);
        double maxerr     = 1.0;
        
        // expand the data up to the required precision level
        while (curNumData<numpoints &&  maxerr>=maxErrorThreshold ) {
            // find the lower index of the bin, where we have the biggest linear interp. error
            maxerr     = 0.0; // value of the current maximum error
            double thexval    = 0.0;
            double theyval    = 0.0;
            int    maxerrindx = 0;   // the lower index of the corresponding bin
            for (int i=0; i<curNumData-1; ++i) {
                
                double xx    = 0.5*(x[i]+x[i+1]);    // mid x point
                double yy    = 0.5*(y[i]+y[i+1]);    // lin. interpolated pdf value at the mid point
                double val   = CalculateDiffCrossSection(tau, xx); // real pdf value at the mid point
                double err   = std::abs(1.-(yy/val));
                if (err>maxerr) {
                    maxerr     = err;
                    maxerrindx = i;
                    thexval    = xx;
                    theyval    = val;
                    
                }
                
            }
            // extend x,y data by putting a new real value at the mid point of the highest error bin
            // first shift all values to the right
            for (int j=curNumData; j>maxerrindx+1; --j) {
                x[j] = x[j-1];
                y[j] = y[j-1];
            }
            // fill x mid point
            x[maxerrindx+1] = thexval;
            y[maxerrindx+1] = theyval;
            
            // increase number of data
            ++curNumData;
            if(curNumData>=numpoints)
            {
                numpoints*=2;
                x.resize(numpoints);
                y.resize(numpoints);
                
            }
        }// end while
        
        x.resize(curNumData);
        y.resize(curNumData);
        return curNumData;
    }
    
    
    void SauterGavrilaPhotoElectricModel::BuildOneLinAlias(int indx, double tau){
        
        std::vector<double> x;
        std::vector<double> y;
        fNumSamplingAngles= PrepareLinAlias(tau, x, y);
        
        fAliasData[indx]              = new LinAlias();
        fAliasData[indx]->fNumdata    = fNumSamplingAngles;
        fAliasData[indx]->fXdata      = new double[fNumSamplingAngles];
        fAliasData[indx]->fYdata      = new double[fNumSamplingAngles];
        fAliasData[indx]->fAliasW     = new double[fNumSamplingAngles];
        fAliasData[indx]->fAliasIndx  = new    int[fNumSamplingAngles];
        
        //INTEGRAL CALCULATION : is not needed, the pdf can be not normalized
        //GLIntegral  *gl   = new GLIntegral(pointsForIntegral, -1., 1.0);
        //std::vector<double> glx  = gl->GetAbscissas();
        //std::vector<double> glw  = gl->GetWeights();
        
        //calculate the integral of the differential cross section
        //double integral=0.;
        //for(int i = 0 ; i < pointsForIntegral ; i++)
        //    integral+= (glw[i]* CalculateDiffCrossSection(tau, glx[i]));
        
        //for (int i = 0; i < (int)y.size(); ++i) {
        //    y[i]/=integral;
       // }
        
        //copy data
        for (int i = 0; i < (int) x.size(); i++)
        {
            fAliasData[indx]->fXdata[i]=x[i];
            fAliasData[indx]->fYdata[i]=y[i];
        }
        
        // prepare the alias data for this PDF(x,y)
        fAliasSampler->PreparLinearTable(fAliasData[indx]->fXdata, fAliasData[indx]->fYdata,fAliasData[indx]->fAliasW, fAliasData[indx]->fAliasIndx,fAliasData[indx]->fNumdata);
        
    }
    
    
}   // namespace geantphysics
