
#ifndef SAUTERGAVRILAPHOTOELECTRICMODEL_H
#define SAUTERGAVRILAPHOTOELECTRICMODEL_H

#include "EMModel.h"


// from geantV
#include "Geant/Config.h"
namespace Geant {
    inline namespace GEANT_IMPL_NAMESPACE {
        class GeantTaskData;
    }
}

namespace geantphysics {
    inline namespace GEANT_IMPL_NAMESPACE {
        class Material;
        class Element;
    }
}


#include <string>

namespace geantphysics {
    
    class MaterialCuts;
    class AliasTable;
    class Particle;
    class LightTrack;
    class XSectionsVector;
    
    /**
     * @brief   Photoelectric model for gamma.
     * @class   SauterGavrilaPhotoElectricModel
     * @author  M Bandieramonte
     * @date    June 2017
     *
     * PhotoElectric model for gamma based on SauterGavrila differential cross section for the calculation of the scattering angle of the secondary particle e- (photoelectron) and on the livermore/epics2014 cross-sections data. Depending on the input energy of the incident gamma, the model provides cross-sections based either on the interpolation of tabulated cross-sections data or on the parameterization of cross-sections data obtained through two fits in two different energy ranges.
     * \cite
     */
    
    class SauterGavrilaPhotoElectricModel : public EMModel {
    public:
        
        /**
         * @name Constructor, destructor:
         */
        //@{
        /**
         * @brief Constructor.
         *
         * @param[in] modelname   Name of the model.
         */
        SauterGavrilaPhotoElectricModel(const std::string &modelname = "SauterGavrilaPhotoElectric");
        /** @brief Destructor. */
        ~SauterGavrilaPhotoElectricModel();
        //@}
        
        
        /**
         * @name Implemented EMModel base class methods:
         */
        //@{
        virtual void Initialize();
        virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle);
        virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts*, double kinenergy, const Particle*);
        virtual int SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td);
        //@}
        
        
        /**
         * @name Model specific public methods:
         */
        //@{
        
        //@}
        
    private:
        
        /**
         * @name Model specific private methods.
         */
        //@{
        
        //---------------------------------------------
        //Initialise: Initialization
        /**
         * @brief Public method to Initialize the model and load the cross-section and parameterizations data.
         *
         */
        void InitializeModel();
        
        
        //---------------------------------------------
        //SampleTargetElementIndex: SampleTargetElementIndex
        /**
         * @brief Public method to Sample the target index
         *
         *
         * @param[in] matcut    material cut to retrieve the element composition of the material
         * @param[in] energy    primary particle kinetic energy (gamma)
         * @param[in] td        GeantTaskData needed to generate random numbers
         * @return              Index of the element sampled from the material composition to be the one involved in photoelectric effect
         */
        
        int SampleTargetElementIndex(const MaterialCuts *matcut, double energy, Geant::GeantTaskData *td);
        
        //---------------------------------------------
        //TestSampleTargetElementIndex: TestSampleTargetElementIndex
        /**
         * @brief Public method to Test the correct sampling of the target element index
         *
         *
         * @param[in] matcut    material cut to retrieve the element composition of the material
         * @param[in] energy    primary particle energy (gamma)
         * @param[in] td        GeantTaskData needed to generate random numbers
         * @return              Output file SampleTargetElementIndexTest_Z that contains the expected pdf and the sampled one
         */
        void TestSampleTargetElementIndex(const MaterialCuts *matcut, double energy, Geant::GeantTaskData *td);
        
        //---------------------------------------------
        //CalculateDiffCrossSection:: Differential cross section based on SauterGavrila distribution for k-shell
        /**
         * @brief Public method to calculate Differential cross section based on SauterGavrila distribution for k-shell/
         *
         *
         * @param[in] energy    primary particle kinetic energy (gamma)
         * @param[in] costheta  cosTheta of the secondary particle (photoelectron)
         * @return              differential cross section based on SauterGavrila distribution for k-shell
         */
        double CalculateDiffCrossSection(double energy, double costheta);
        
        
        //---------------------------------------------
        //CalculateDiffCrossSectionLog:: Differential cross section based on SauterGavrila distribution for k-shell, but using the transformed variable xsi
        /**
         * @brief Public method to calculate Differential cross section based on SauterGavrila distribution for k-shell  but using the transformed variable xsi/
         *
         *
         * @param[in] energy    primary particle kinetic energy (gamma)
         * @param[in] costheta  cosTheta of the secondary particle (photoelectron)
         * @return              differential cross section based on SauterGavrila distribution for k-shell
         */
        double CalculateDiffCrossSectionLog(double energy, double costheta);
        
        
        //---------------------------------------------
        //ComputeXSectionPerAtom Cross-section per atom
        /**
         * @brief Public method to compute the cross-section per atom for the photoelectric effect. The computation is performed in a different way depending on the Z of the element and on the energy of the incident gamma. There are two parameterizations (low and high energy) and two sets of tabulated cross-section data (below k-shell binding energy and above k-shell binding energies).
         * It will use:
         *   - le-cs data: for energies below k-shell binding energy - Linear interpolation of cs tabulated data (as in Geant4)
         *   - cs data   : for energies above k-shell binding energy - Spline interpolatio of cs tabulated data (as in Geant4)
         *   - low-param : for energies above the lowParamThreshold (i.e fParamLow[Z]))[0])
         *   - high param: for energies above the highParamThreshold (i.e fParamHigh[Z]))[0])
         *
         *
         * @param[in] zeta      Z of the element involved in the pe interaction
         * @param[in] energy    primary particle kinetic energy (gamma)
         * @return              cross-section per atom for the photoelectric effect
         */
        double ComputeXSectionPerAtom(double zeta, double energy);
        
        
        //---------------------------------------------
        //SamplePhotoElectronDirection_Alias: Sample PhotoElectron direction with Alias sampling
        /**
         * @brief Public method to sample PhotoElectron direction with Alias sampling. The tables are created in the
         *
         *
         * @param[in] energy    primary particle kinetic energy (gamma)
         * @param[in] r1        random number used for sampling
         * @param[in] r2        random number used for sampling
         * @param[in] r3        random number used for sampling
         * @return              cosTheta of the secondary particle (photoelectron e-)
         */
        double SamplePhotoElectronDirection_Alias(double energy,
                                                  double r1,
                                                  double r2,
                                                  double r3);
        
        
        //---------------------------------------------
        //SamplePhotoElectronDirection_Rejection: Sample PhotoElectron direction with rejection sampling
        /**
         * @brief Public method to sample PhotoElectron direction with rejection sampling.
         *
         *
         * @param[in]  energy       primary particle kinetic energy (gamma)
         * @param[in]  td           Geant::GeantTaskData used to generate random numbers
         * @param[out] sintheta     sin of the polar angle (theta) of the secondary particle (photoelectron e-)
         * @param[out] costheta     cos of the polar angle (theta) of the secondary particle (photoelectron e-)
         * @param[out] phi          azimuthal angle of the secondary particle (photoelectron e-)
         *
         */
        void SamplePhotoElectronDirection_Rejection(double energy,
                                                    double &sintheta,
                                                    double &costheta,
                                                    double &phi,
                                                    Geant::GeantTaskData *td);
        

        //---------------------------------------------
        //GetValue
        /**
         * @brief Public method to
         *
         *
         * @param[in]  energy       primary particle kinetic energy (gamma)
         * @param[in]  Z            Atomic number of the element
         * @param[in]  shellIdx
         * @return     value
         */
        
        inline double GetValue(double energy, int Z, size_t shellIdx){
            size_t bin = 0;

            size_t numberofnodes= fShellCrossSection[Z]->fCompLength[shellIdx];
            //to do: check this: it needs a value
            if(energy < fShellCrossSection[Z]->fCompBinVector[shellIdx][1]) return std::min(bin, numberofnodes-2);
            if(energy >= fShellCrossSection[Z]->fCompBinVector[shellIdx][numberofnodes-2]) return numberofnodes - 2;
            if(bin >= numberofnodes || energy < fShellCrossSection[Z]->fCompBinVector[shellIdx][bin] || energy > fShellCrossSection[Z]->fCompBinVector[shellIdx][bin+1])
                bin = std::lower_bound(fShellCrossSection[Z]->fCompBinVector[shellIdx].begin(), fShellCrossSection[Z]->fCompBinVector[shellIdx].end(), energy) - fShellCrossSection[Z]->fCompBinVector[shellIdx].begin() - 1;
            bin=std::min(bin, numberofnodes-2);
        
            //to do: check performance here
            return LinearInterpolation (energy,fShellCrossSection[Z]->fCompBinVector[shellIdx],fShellCrossSection[Z]->fCompDataVector[shellIdx], bin);
        }

        
        //---------------------------------------------
        //Linear interpolation
        inline double LinearInterpolation(double energy, std::vector<double>   &binvector, std::vector<double>   &datavector,  size_t idx) const
        {
            // Linear interpolation is used to get the interpolated value for lowEnergy cross sections (below K-shell binding energy).
            //Before this method is called it is ensured that the energy is inside the bin
            // 0 < idx < numberOfNodes-1
            return datavector[idx] +( datavector[idx + 1]-datavector[idx] ) * (energy - binvector[idx]) /( binvector[idx + 1]-binvector[idx] );
        }
        
        //@}
        
        // data members
    private:
        static const int gMaxSizeData = 102;        //Maximum number of Z elements
        static const int gNShellLimit= 100;         //Maximum number of shells per element - to do: check this value
        
        //to do: this could be seen as a struct with an array of XSectionsVector + the int* fCompID;
        struct ShellData{
            
            std::vector<double>* fCompBinVector;    // Bins for every shell of element Z
            std::vector<double>* fCompDataVector;   // Data for every shell of element Z
            int* fCompID;                           // Id for every shell of element Z
            size_t* fCompLength;                    // Total number of shells per element Z
        };
        
        static std::vector<double>*  fParamHigh[gMaxSizeData];  //High-energy parameterization data
        static std::vector<double>*  fParamLow[gMaxSizeData];   //Low-energy parameterization data
        
        void   LoadData();
        void   ReadData(int Z);
        void   InitSamplingTables();
        void   BuildOneLinAlias(int indxlalias, double gcut);
        
        int  fVerboseLevel;                         //Verbose level to control the printout
        //bool fDeexcitationActive;                 //True if deexitation is active - not used at the moment
        
        static ShellData  **fShellCrossSection;            //Several shells cross-sections data per Z
        
        XSectionsVector * fLECSVector[gMaxSizeData];  //one LE cross-section vector per Z
        XSectionsVector * fCSVector[gMaxSizeData];    //one !LE cross-section vector per Z
        
        //to do:  check the use of these members
        bool* fCrossSection;                        //true if we have CrossSections data (for energies above k-shell binding energy)
        bool* fCrossSectionLE;                      //true if we have Low-Energy CrossSections data (for energies below k-shell binding energy)
        
        static int                   fNShells[gMaxSizeData];
        static int                   fNShellsUsed[gMaxSizeData];
        static Material*             fWater;
        static double                fWaterEnergyLimit;
        
        int      fSecondaryInternalCode; // electron GV code set at initialization
        
        /**
         * @name Members to describe the discrete photon energy grid for sampling tables:
         *
         * These variables describe and define .......
         *
         *
         */
        //@{
        //---------------------------------------------
        /** @brief Number of gamma kinetic energy grid points in [SauterGavrilaPhotoElectricModel::fMinGammaEnergy,
         *        SauterGavrilaPhotoElectricModel::fMaxGammaEnergy].
         */
        int     fNumSamplingPrimEnergies;
        
        
        /** @brief Number of primary gamma kinetic energy grid points per decade. */
        int     fNumSamplingPrimEnergiesPerDecade;
        
        
        //---------------------------------------------
        /** @brief Number of emitted pe cosTheta in [-1, 1] or transformed emitted photoelectron angle related variable in [log(e-12),log(2)]. */
        int     fNumSamplingAngles;
        
        //---------------------------------------------
        /** @brief Minimum of the gamma kinetic energy grid (default 1.0 [keV] that is the minimum available.) */
        double  fMinPrimEnergy;
        
        //---------------------------------------------
        /** @brief Maximum of the gamma kinetic energy grid (default 10.0 [GeV] that is the maximum available.) */
        double  fMaxPrimEnergy;
        
        //---------------------------------------------
        /** @brief Logarithm of SauterGavrilaPhotoElectricModel::fMinGammaEnergy i.e. ln(fMinGammaEnergy) . */
        double  fPrimEnLMin;                         // log min gamma energy
        
        //---------------------------------------------
        /** @brief Inverse of the gamma kinetic energy grid delta i.e.
         *        ln[fMaxGammaEnergy/fMinGammaEnergy]/(fNumSamplingGammaEnergies-1)
         */
        double  fPrimEnILDelta;                      // 1 log delta gamma energy of the gamma energy grid
        
        //---------------------------------------------
        /** @brief The logarithmically spaced gamma kinetic energy grid.
         *
         *        Size of the array is SauterGavrilaPhotoElectricModel::fNumSamplingGammaEnergies points in the
         *        [SauterGavrilaPhotoElectricModel::fMinGammaEnergy, SauterGavrilaPhotoElectricModel::fMaxGammaEnergy] interval.
         */
        double *fSamplingPrimEnergies;   // the common gamma energy grid which we build sampling tables above
        
        //---------------------------------------------
        /** @brief The logarithm of SauterGavrilaPhotoElectricModel::fSamplingGammaEnergies grid.
         *
         *  Size of the array is SauterGavrilaPhotoElectricModel::fNumSamplingGammaEnergies points in the
         *  [ln(SauterGavrilaPhotoElectricModel::fMinGammaEnergy), ln(SauterGavrilaPhotoElectricModel::fMaxGammaEnergy)]
         *  interval.
         */
        double *fLSamplingPrimEnergies;            // log of sampling gamma energies
        //@}
        
        
        
        //---------------------------------------------
        /** @brief Internal data structure to store data for sampling the emitted photoelectron direction
         *
         *  This data structure is set up at initialisation
         *  over the gamma kinetic energy grid to sample the emitted photoelectron angle using a combination
         *  of Walker's alias sampling and liner approximation. At most we have as many data structure as
         *  SauterGavrilaPhotoElectricModel::fNumSamplingGammaEnergies
         *  and these data structure pointers are stored in the SauterGavrilaPhotoElectricModel::fAliasData linear array.
         */
        struct LinAlias{
            /** @brief Number of data points i.e. size of the arrays = SauterGavrilaPhotoElectricModel::fNumSamplingGammaEnergies. */
            int     fNumdata;
            
            /** @brief This must be the gamma energies - maybe beginning of the bin? (transformed or not..it depends). */
            double *fXdata;
            /** @brief The probability density function values (not necessarily normalised) over the photoelectron angle
             *        variable values.
             */
            double *fYdata;
            
            /** @brief The alias probabilities (not necessarily normalised) over the photoelectron angle
             *        variable values.
             */
            double *fAliasW;
            /** @brief The alias indices over the photon energy variable values. */
            int    *fAliasIndx; // alias indices
        };
        
        //---------------------------------------------
        /** @brief Linear array to store pointers to LinAlias data structures.
         *
         * The size is
         * SauterGavrilaPhotoElectricModel::fNumSamplingPrimEnergies.
         */
        LinAlias   **fAliasData;                   //alias data structure
        
        //---------------------------------------------
        /** @brief An alias sampler used at run-time sampling of the emitted photon energy. */
        AliasTable  *fAliasSampler;
        
    };
    
}      // namespace geantphysics

#endif // SAUTERGAVRILAPHOTOELECTRICMODEL_H
