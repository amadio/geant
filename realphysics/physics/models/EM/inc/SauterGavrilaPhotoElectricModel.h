
#ifndef SAUTERGAVRILAPHOTOELECTRICMODEL_H
#define SAUTERGAVRILAPHOTOELECTRICMODEL_H

#include "EMModel.h"
#include "Spline.h"

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
    
    /**
     * @brief   Photoelectric model for gamma.
     * @class   SauterGavrilaPhotoElectricModel
     * @author  M Bandieramonte
     * @date    June 2017
     *
     * PhotoElectric model for gamma based on SauterGavrila differential cross section
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
        
        //int   GetNumberOfPhotonEnergiesPerDecade()   const { return fNumSamplingPrimEnergiesPerDecade; }
        
        
        //@}
        
    private:
        /**
         * @name Model specific private methods.
         */
        //@{
        
        //---------------------------------------------
        //Initialise: Initialization
        /**
         * @brief Public method to Initialize the Model
         *
         */
        void InitializeModel();
        
        
        //---------------------------------------------
        //MySampleTargetElementIndex: MySampleTargetElementIndex
        /**
         * @brief Public method to Sample the target index
         *
         */

        int MySampleTargetElementIndex(const MaterialCuts *matcut, double energy, Geant::GeantTaskData *td);
        
        
        //---------------------------------------------
        //CalculateDiffCrossSection:: Differential cross section based on SauterGavrila distribution for k-shell
        /**
         * @brief Public method to calculate Differential cross section based on SauterGavrila distribution for k-shell/
         *
         *
         * @param[in] energy
         * @param[in] costheta
         * @return
         */
        double CalculateDiffCrossSection(double energy, double costheta);
        
        
        //---------------------------------------------
        //CalculateDiffCrossSection:: Differential cross section based on SauterGavrila distribution for k-shell, but using the transformed variable xsi
        /**
         * @brief Public method to calculate Differential cross section based on SauterGavrila distribution for k-shell  but using the transformed variable xsi/
         *
         *
         * @param[in] energy
         * @param[in] costheta
         * @return
         */
        double CalculateDiffCrossSectionLog(double energy, double costheta);
        
        
        //---------------------------------------------
        //ComputeXSectionPerAtom Cross-section per atom
        /**
         * @brief Public method to compute the cross-section per atom for the photoelectric effect
         *
         *
         * @param[in] zeta
         * @param[in] energy
         * @return
         */
        double ComputeXSectionPerAtom(double zeta, double energy);
        
        
        
        //---------------------------------------------
        //SamplePhotoElectronDirection_Alias: Sample PhotoElectron direction with Alias sampling
        /**
         * @brief Public method to sample PhotoElectron direction with Alias sampling
         *
         *
         * @param[in] gammaenin
         * @param[in] r1
         * @param[in] r2
         * @param[in] r3
         * @return
         */
        double SamplePhotoElectronDirection_Alias(double gammaenin,
                                                  double r1,
                                                  double r2,
                                                  double r3);
        
        
        //---------------------------------------------
        //SamplePhotoElectronDirection_Rejection: Sample PhotoElectron direction with rejection sampling
        /**
         * @brief Public method to sample PhotoElectron direction with rejection sampling
         *
         *
         * @param[in]  gammaekinin
         * @param[in]  td
         * @param[out] sintheta
         * @param[out] costheta
         * @param[out] phi
         * @return
         */
        void SamplePhotoElectronDirection_Rejection(double gammaekinin,
                                                    double &sintheta,
                                                    double &costheta,
                                                    double &phi,
                                                    Geant::GeantTaskData *td);
        
        //---------------------------------------------
        //Method to retrieve, given the energy of the incoming gamma, the corresponding bin index of the CrossSectionVector
        inline size_t FindCSBinLocation(double energy,  size_t index, size_t numberofnodes, std::vector<double>   binvector){
            
            size_t bin= index;
            if(energy < binvector[1]) {
                bin = 0;
            } else if(energy >= binvector[numberofnodes-2]) {
                bin = numberofnodes - 2;
            } else if(bin >= numberofnodes || energy < binvector[bin]
                      || energy > binvector[bin+1])
            {
                // Bin location proposed by K.Genser (FNAL) from G4
                bin = std::lower_bound(binvector.begin(), binvector.end(), energy) - binvector.begin() - 1;
            }
            size_t minV=std::min(bin, numberofnodes-2);
            return minV;
        }
        
        //---------------------------------------------
        //Linear interpolation
        inline double LinearInterpolation(double energy, std::vector<double>   binvector, std::vector<double>   datavector,  size_t idx) const
        {
            // Linear interpolation is used to get the interpolated value for lowEnergy cross sections (below K-shell binding energy).
            //Before this method is called it is ensured that the energy is inside the bin
            // 0 < idx < numberOfNodes-1
            //std::cout<<"LinearInterpolation for index: "<<idx<<" : "<<datavector[idx]<<" ---- "<<datavector[idx+1]<<" ---- "<<binvector[idx]<<" ---- "<<binvector[idx+1]<<" ---- energy: "<<energy<<"\n";
            return datavector[idx] +
            ( datavector[idx + 1]-datavector[idx] ) * (energy - binvector[idx]) /( binvector[idx + 1]-binvector[idx] );
        }
        
        //@}
        
        // data members
    private:
        static const int gMaxSizeData = 120;        //Maximum number of Z elements
        static const int gNShellLimit= 100;         //Maximum number of shells per element - controllare questo numero
        
        
        struct ShellData{
            
            std::vector<double>* fCompBinVector;    // Bins for every shell of element Z
            std::vector<double>* fCompDataVector;   // Data for every shell of element Z
            int* fCompID;                           // Id for every shell of element Z
            size_t* fCompLength;                    // Total number of shells per element Z
        };
        
        
        //Struct to handle tabulated cross-sections data -> one struct per element
        struct CrossSectionsVector{
            
            std::vector<double>   fBinVector;       //Cross sections bin vector (i.e. x coordinate)
            std::vector<double>   fDataVector;      //Cross sections data vector (i.e. y coordinate)
            
            size_t numberOfNodes;                   // Number of elements
            double edgeMin;                         // Energy of first point
            double edgeMax;                         // Energy of last point
            Spline     *sp;                         // Spline interpolator
            
        };
        
        
        static std::vector<double>*  fParamHigh[gMaxSizeData];  //High-energy parameterization data
        static std::vector<double>*  fParamLow[gMaxSizeData];   //Low-energy parameterization data
        
        
        void   LoadData();
        void   ReadData(int Z);
        void   InitSamplingTables();
        void   BuildOneLinAlias(int indxlalias, double gcut);
        
        int  fVerboseLevel;                         //Verbose level to control the printout
        //bool fDeexcitationActive;                 //True if deexitation is active - not used at the moment
        
        ShellData  **fShellCrossSection;            //Several shells cross-sections data per Z
        CrossSectionsVector ** fLECSVector;         //one LE cross-section struct per Z
        CrossSectionsVector ** fCSVector;           //one !LE cross-section struct per Z
        
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
