
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
     * @brief  This class implements a PhotoElectric model for gamma based on SauterGavrila differential cross section for the sampling of the 
     * scattering angle of the secondary particle (e-) and on the EPICS2014 cross-sections data.
     *
     * @class   SauterGavrilaPhotoElectricModel
     * @author  M Bandieramonte
     * @date    June 2017
     *
     * The photoelectric effect is the ejection of an electron from a material after a photon has been absorbed by that material.
     * Depending on the input energy of the incident gamma, the model provides cross-sections based either on the interpolation of tabulated
     * cross-sections data or on the parameterization of cross-sections data obtained through two fits in two different energy ranges (low and 
     * high). 
     * Note that for this process the lambda table is not built, since the cross-section is not a smooth function of the energy,
     * therefore in all calculations the cross section is used directly.
     * Final state sampling: The incident photon is absorbed and an electron is emitted with a direction that is calculated
     * using the Sauter-Gavrila distribution (for K-shell) \cite sautergavrila.
     * The electron kinetic energy is the difference between the incident photon energy and the binding energy of the electron before the
     * interaction. The sub-shell, from which the electron is emitted, is randomly selected according to the relative cross-sections of all
     * subshells, determined at the given energy \f$E_0\f$, by interpolating the evaluated cross-section data from the EPICS (Electron Photon
     * Interaction Cross Sections) v.2014 data bank \cite epics2014.
     * The interaction leaves the atom in an excited state; the deexcitation process is not implemented yet.
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
         * @param[in] modelname     Name of the model.
         * @param[in] aliasActive   Boolean true if we want to use Alias Sampling to sample the secondary particle direction. By default is set to false and 
         *                          model uses composition-rejection sampling
         */
        SauterGavrilaPhotoElectricModel(const std::string &modelname = "SauterGavrilaPhotoElectric", bool aliasActive = false);
        
        /** @brief Destructor. */
        ~SauterGavrilaPhotoElectricModel();
        //@}
        
        /**
         * @name Implemented EMModel base class methods:
         */
        //@{
        virtual void Initialize();
        
        //ComputeMacroscopicXSection
        /**
         *@brief Method to compute macroscopic cross section for a given MaterialCuts, Particle, kinetic energy.
         *
         * This method is called at initialization to build the lambda tables through the corresponding PhysicsProcess by the
         * PhysicsManagerPerParticle object. This method is implemented even if for photoelectric effect the table is not built 
         * to be used at run-time, but the cross-sections are calculated on the flight.
         *
         * @param[in] matcut	  Pointer to the MaterialCuts object in which the macroscopic cross section must be computed.
         * @param[in] kinenergy	  Kinetic energy of the Particle at which the macroscopic cross section must be computed.
         * @param[in] particle	  Pointer to the Particle object for which the macroscopic cross section must be computed.
         * @return                Macroscopic cross section computed by the given electromagnetic model in internal [1/length] units for the given 
         * Particle, MaterialCuts/Material and particle kinetic energy combination.
         */
        virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle);
        
        // Documentation from the EMModel base class method
        virtual double ComputeXSectionPerAtom(const Element *elem, const MaterialCuts*, double kinenergy, const Particle*);
        
        //SampleSecondaries
        /**
         * @brief Method for the final state sampling. 
         * Primary track properties are updated and secondary track(s) are generated according to the photoelectric effect interaction.
         * After photoelectric interaction the incident photon is absorbed and an electron is emitted. 
         * The electron kinetic energy is the difference between the incident photon energy and the binding energy of the electron 
         * before the interaction. The sub-shell, from which the electron is emitted, is randomly
         * selected according to the relative cross-sections of all subshells, determined at the given energy. The interaction leaves the atom in 
         * an excited state but the deexcitation of the atom is not simulated yet.
         *
         * @param[in,out] track     Primary track. At input, it stores the pre-interaction primary particle properties and
         *                          some information about the current material-cut couple (if needed). It is updated by the method and
         *                          it stores the post-interaction primary track properties at output. After the interaction the gamma is
         *                          absorbed and the status of the track is set to "kKill"
         * @param[in,out] td        Pointer to a Geant thread local data object. At output, its fPhysicsData object will store
         *                          the secondary tracks generated in the interaction.
         * @return                  Number of secondary tracks generated in the interaction.
         */
        virtual int SampleSecondaries(LightTrack &track, Geant::GeantTaskData *td);
        //@}
        
        
        /**
         * @name Model specific public methods:
         */
        //@{
        //SetVerboseLevel
        /**
         * @brief Public method to set the verbose level (fVerboseLevel attribute).*/
        void  SetVerboseLevel(int lev);
        int   GetNumberOfPhotonEnergiesPerDecade()   const { return fNumSamplingPrimEnergiesPerDecade; }
        // before initialisation
        void  SetNumberOfPhotonEnergiesPerDecade(int val)  { fNumSamplingPrimEnergiesPerDecade = val;  }
        //
        int   GetNumberOfDiscretePDFSamples()        const { return fNumSamplingAngles;              }
        // before initialisation
        void  SetNumberOfDiscretePDFSamples(int val)       { fNumSamplingAngles = val;               }
        
        //@}
        
    private:
        
        /**
         * @name Model specific private methods.
         */
        //@{
        
        //---------------------------------------------
        //InitializeModel
        /**
         * @brief Internal method to initialize the model and load and store the cross-section and parameterizations data.
         *
         */
        void InitializeModel();
        
        
        //---------------------------------------------
        //SampleTargetElementIndex
        /**
         * @brief Private method to Sample the target index of the element involved in the interaction.
         *
         *
         * @param[in] matcut    MaterialCuts to retrieve the element composition of the material.
         * @param[in] energy    Primary particle (gamma) kinetic energy.
         * @param[in] td        GeantTaskData needed to generate random numbers.
         * @return              Index of the element sampled from the material composition to be the one involved in photoelectric effect.
         */
        
        size_t SampleTargetElementIndex(const MaterialCuts *matcut, double energy, Geant::GeantTaskData *td);
        
        //---------------------------------------------
        //TestSampleTargetElementIndex
        /**
         * @brief Private method to Test the correct sampling of the target element index.
         * This is a unit test for SampleTargetElementIndex method. The methos is calculating the Macroscopic cross-section for all the materials 
         * active in the simulation and storing it in xsec vector. Then SampleTargetElementIndex is called 10^9 times and sampled index used to 
         * increment the corresponding entry of xsecSampled vector. The xsec and xsecSampled (normalized) vectors are then stored in an output file.
         *
         * @param[in] matcut    MaterialCuts to retrieve the element composition of the material.
         * @param[in] energy    primary particle (gamma) energy.
         * @param[in] td        GeantTaskData needed to generate random numbers.
         * @return              Output file SampleTargetElementIndexTest_Z that contains the expected pdf and the sampled one.
         */
        void TestSampleTargetElementIndex(const MaterialCuts *matcut, double energy, Geant::GeantTaskData *td);
        
        //---------------------------------------------
        //CalculateDiffCrossSection
        /**
         * @brief Private method to calculate Differential cross section based on SauterGavrila distribution for k-shell.
         *
         * SauterGavrila approximation for K-shell, correct to the first \f$\alpha Z \f$ order
         *
         * @param[in] energy    primary particle (gamma) kinetic energy.
         * @param[in] costheta  cosTheta of the secondary particle (photoelectron).
         * @return    dsigma    differential cross section based on SauterGavrila distribution (K-shell only).
         */
        double CalculateDiffCrossSection(double energy, double costheta);
        
        
        //---------------------------------------------
        //CalculateDiffCrossSectionLog
        /**
         * @brief Private method to calculate Differential cross section based on SauterGavrila distribution for k-shell  but using the transformed variable \f$xsi\f$.
         *
         *
         * @param[in] energy    primary particle (gamma) kinetic energy.
         * @param[in] costheta  cosTheta of the secondary particle (photoelectron).
         * @return              differential cross section based on SauterGavrila distribution for k-shell.
         */
        double CalculateDiffCrossSectionLog(double energy, double costheta);
        
        
        //---------------------------------------------
        //ComputeXSectionPerAtom Cross-section per atom
        /**
         * @brief Private method to compute the cross-section per atom for the photoelectric effect. The computation is performed in a different way depending on the Z of the element and on the energy of the incident gamma. There are two parameterizations (low and high energy) and two sets of tabulated cross-section data (below k-shell binding energy and above k-shell binding energy).
         *
         * The total photoelectric and single shell cross-sections are tabulated from threshold to a lowParamThreshold that starts from 5keV
         * and increases with the k-shell binding energy of the corresponding element. Above the lowParamThreshold, EPICS2014 cross sections 
         * \cite epics2014 are parameterized as following:
         *
         * \f[ 
         * \sigma(E)= \frac{a_1}{E}+\frac{a_2}{E^2} +\frac{a_3}{E^3} +\frac{a_4}{E^4} +\frac{a_5}{E^5}+\frac{a_6}{E^6}
         * \f]
         * This method will use the following data files:
         *   - le-cs data : for energies below k-shell binding energy - Linear interpolation of low-energy cs tabulated data (as done in Geant4).
         *   - cs data    : for energies above k-shell binding energy - Spline interpolation of cs tabulated data (as done in Geant4).
         *   - low-param  : for energies above the lowParamThreshold  (stored in fParamLow  [Z][0] ).
         *   - high param : for energies above the highParamThreshold (stored in fParamHigh [Z][0] ).
         *
         * To avoid tracking problems for very low-energy gamma, the photoelectric cross section is not equal to zero below the first ionisation
         * potential but it has a constant value (equal to corresponding value at the ionization potential). 
         * This means that all types of media are not transparent for gamma.
         *
         * @param[in] zeta      Z of the element involved in the pe interaction.
         * @param[in] energy    primary particle (gamma) kinetic energy.
         * @return              photoelectric effect cross-section per atom.
         */
        double ComputeXSectionPerAtom(double zeta, double energy);
        
        
        //---------------------------------------------
        //SamplePhotoElectronDirection_Alias
        /**
         * @brief Private method to sample PhotoElectron direction corresponding to a specific gamma kinetic energy, with Alias sampling.
         * The tables are created at initialization phase.
         *
         * @param[in] energy    primary particle (gamma) kinetic energy.
         * @param[in] r1        random number used for sampling.
         * @param[in] r2        random number used for sampling.
         * @param[in] r3        random number used for sampling.
         * @return              cosTheta of the secondary particle (photoelectron e-).
         */
        double SamplePhotoElectronDirection_Alias(double energy,
                                                  double r1,
                                                  double r2,
                                                  double r3);
        
        
        //---------------------------------------------
        //SamplePhotoElectronDirection_Rejection
        /**
         * @brief Private method to sample PhotoElectron direction with rejection sampling.
         *
         * The polar angle of the photoelectron is sampled from the Sauter-Gavrila distribution (for K-shell) \cite sautergavrila, which is correct only to zero order in \f$\alpha Z\f$ :
         *  \f[
         *    \frac{d\sigma}{d\cos\theta} \simeq \frac{{\sin^2\theta}}{(1-\beta\cos\theta)^4} \left\{1+\frac{1}{2}\gamma (\gamma -1)(\gamma -2) (1-\beta\cos\theta) \right\}
         *  \f]
         *
         * where \f$\beta\f$ and \f$\gamma\f$ are the Lorentz factors of the photoelectron. \f$\f$
         * Introducing the variable \f$\nu = 1 - cos \theta_e \f$, the angular distribution can be expressed as:
         *
         *  \f[
         *    p(\nu) = (2-\nu) [\frac{1}{A+\nu} + \frac{1}{2} \beta \gamma (\gamma -1)(\gamma -2)] \frac{\nu}{(A+\nu)^3}
         *  \f]
         *  apart from a normalisation constant. With
         *  \f[
         * \gamma = 1+ \frac{E_e}{m_e c^2}\text{,  }  A=\frac{1}{\beta}-1
         * \f]
         * where  \f$E_e\f$ is the electron energy, \f$m_e\f$ its rest mass and \f$\beta\f$ its velocity in units of the speed of light \f$c\f$.
         * Random sampling of \f$\nu\f$ from this distribution can be performed analytically. For more details have a look at the Penelope manual 
         * \cite salvat2006penelope.
         * Though the Sauter distribution, strictly speaking, is adequate only for ionisation of the K-shell by high-energy photons, in many
         * practical simulations it does not introduce appreciable errors in the description of any photoionisation event, irrespective of the
         * atomic shell or of the photon energy.
         * @param[in]  energy        Primary particle (gamma) kinetic energy.
         * @param[in]  td            Geant::GeantTaskData used to generate random numbers.
         * @param[out] costheta      Cosinus of the polar angle (theta) of the secondary particle (photoelectron e-).
         * 
         *
         */
        void SamplePhotoElectronDirection_Rejection(double energy,
                                                    double &costheta,
                                                    Geant::GeantTaskData *td);
        
        
        //---------------------------------------------
        //LoadData
        /**
         * @brief Internal method to load parameterization data.
         *
         *  Used at initialisation of the model to load parameterization data used to calculate cross-sections. It selects only the materials
         *  present in the list of active regions and for each of them calls the internal method ReadData to read and store the corresponding
         * parameterization files.
         *
         **/
        void   LoadData();
        
        //---------------------------------------------
        //ReadData
        /**
         * @brief Internal method to read parameterization data corresponding to the 'active' element with atomic number Z.
         *
         *  Used at initialisation of the model to read and store in the corresponding class members parameterization data 
         *  of the selected element with atomic number Z.
         *
         *  The files read are:
         *  - pe-le-cs-zeta.dat:    Low energy (below k-shell binding energy) cross-section data. They are stored in fLECSVector[Z].
         *  - pe-cs-zeta.dat:       Cross-section data above k-shell binding energy. They are stored in fCSVector[Z].
         *  - pe-low-zeta.dat:      Low-energy parameterization data. They are stored in fParamLow[Z].
         *  - pe-high-zeta.dat:     High-energy paramterization data. They are stored in fParamHigh[Z].
         *  - pe-ss-cs-zeta.dat:    Subshells cross-sections data. They are stored in fShellVector[Z].
         *
         *  @param[in]  zeta   Atomic number of the element for which to load the corresponding parameterization data.
         **/
        void   ReadData(int zeta);
        
        //---------------------------------------------
        //InitSamplingTables
        /**
         * @brief Internal method to build (post interaction) photoelectron angle sampling tables.
         *
         *  Used at initialisation of the model to prepare sampling tables over an initial gamma energy grid.
         *  The gamma energy grid is determined by the low/high energies at which the angle can be calculated and the
         *  fNumSamplingPrimEnergiesPerDecade variable. For energies higher than 100MeV the secondary is considered to go straight
         *  and being ejected from the atom with the same direction as the incident gamma.
         *  A sampling table, using Walker's discrete alias method combined with <em>linear approximation of the p.d.f.</em>, is 
         *  built at each gamma energy grid point using the BuildOneLinAlias() method. 
         *  These sampling tables are used at run-time to sample the secondary (photoelectron) direction.
         **/
        void   InitSamplingTables();
        
        //---------------------------------------------
        //BuildOneLinAlias
        /** @brief Internal method to build photoelectron angle sampling tables for a given initial gamma energy.
         *
         *  This method is used by the InitSamplingTables() method to build sampling table (using Walker's discrete alias
         *  method combined with <em>linear approximation of the p.d.f.</em>) at a given initial gamma energy. The method will
         *  call the CalculateDiffCrossSection internal method to compute the Differential cross section based on SauterGavrila 
         *  distribution for k-shell \cite sautergavrila.
         *
         *  @param[in]  indx   Index of the alias table data structure to build in the container.
         *  @param[in]  tau    Initial photon energy \f$ E_0 \f$ expressend in \f$e_m c^2\f$ units.
         */
        void   BuildOneLinAlias(int indx, double tau);
        
        
        
        /**
         * @brief Public method to prepare sampling table for discretized continuous distribution with combination of alias
         *        sampling and linear approximation of the p.d.f.  
         *        This method is preparing data structure for Alias Sampling using an adaptive binning, The number of points used for 
         *        binning is not fixed but varies from energy to energy to meet the precision requirements that is set to the value of
         *        gsingleTableErrorThreshold variable
         *
         *
         * @param[in,out] xdata      Array of discrete samples of the random variable between its minimum and maximum values.
         *
         * @param[in,out] ydata      Array of the (not necessarily normalised) p.d.f. at the discrete sample points of the
         *                           random variable given in xdata. It will be used at sampling as well. 
         *
         * @param[in]     tau        Initial photon energy \f$ E_0 \f$ expressend in \f$e_m c^2\f$ units.
         *
         **/
        
         int PrepareLinAlias(double tau, std::vector<double> & x, std::vector<double> & y);
        
        //@}
        
        // data members
    private:
        /** @brief Maximum number of Z elements. */
        static const int        gMaxSizeData                = 100;          //Maximum number of Z elements
        /** @brief Maximum number of shells per element. */
        static const int        gNShellLimit                = 100;          //Maximum number of shells per element
        /** @brief Maximum error introduced by the use of Alias sampling at each decade energy*/
        static constexpr double gsingleTableErrorThreshold   = 2.e-3;       //2 per mille error threshold
        //static const int        gpointsForIntegral           = 200;       //not used for the moment
        
        /** @brief Vector storing high-energy parameterization data. */
        static std::vector<double>*  fParamHigh[gMaxSizeData];   //High-energy parameterization data
        /** @brief Vector storing low-energy parameterization data. */
        static std::vector<double>*  fParamLow [gMaxSizeData];   //Low-energy parameterization data
        
        /** @brief Verbose level to control the printout. */
        int  fVerboseLevel;                         //Verbose level to control the printout
        //bool fDeexcitationActive;                 //True if deexitation is active - not used at the moment

        
        /** @brief Vector of pointers to XSectionsVector cross-sections. Several subshell XSectionsVector per Z. */
        XSectionsVector** fShellVector[gMaxSizeData];     //Several subshell cross-section vector per Z
        
        /** @brief Array of pointers to Low-energy XSectionsVector (one LE cross-section vector per Z).*/
        XSectionsVector * fLECSVector[gMaxSizeData];  //one LE cross-section vector per Z
        
        /** @brief Array of pointers to High-energy XSectionsVector vector (one !LE cross-section vector per Z). */
        XSectionsVector * fCSVector[gMaxSizeData];    //one !LE cross-section vector per Z
        
        //to do:  check the use of these members
        /** @brief Vector of booleans. fCrossSection[Z] is true if there are CrossSections data (for energies above k-shell binding energy) for element Z. */
        bool* fCrossSection;                        //true if there are CrossSections data (for energies above k-shell binding energy)
        
        /** @brief Vector of booleans. fCrossSectionLE[Z] is true if there are Low-Energy CrossSections data (for energies below k-shell binding energy) for element Z. */
        bool* fCrossSectionLE;                      //true if there are Low-Energy CrossSections data (for energies below k-shell binding energy)
        
        /** @brief Total number of shells per each element Z. */
        static int                   fNShells[gMaxSizeData];
        
        /** @brief Number of used shells per each element Z. */
        static int                   fNShellsUsed[gMaxSizeData];
        
        /** @brief Pointer to Material to handle water as a special case. */
        static Material*             fWater;
        
        /** @brief Water energy limit. */
        static double                fWaterEnergyLimit;
        
        /** @brief Secondary particle (e-) GeantV code, set at initialization. */
        int      fSecondaryInternalCode; // electron GV code set at initialization
        
        /**
         * @name Members to describe the discrete gamma energy grid for sampling tables that will be used to sample the photoelectron (e-) direction
         *
         * These variables describe and define the primary gamma kinetic energy grid on which we build sampling tables in the
         * InitSamplingTables() method for run-time sampling of the post interaction photoelectron (e-) direction.
         * The min/max of the table is set to fMinPrimEnergy and fMaxPrimEnergy limits and the number of discrete gamma energy points
         * will be determined by the value of SauterGavrilaPhotoElectricModel::fNumSamplingPrimEnergiesPerDecade variable. 
         * The default value can be changed by the user through the SetNumberOfPhotonEnergiesPerDecade()(to do: add this method)
         * public method (must be set before the initialisation of the model!).
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
        /** @brief Number of emitted photoelectron cosTheta considered in the range[-1, 1] or number of transformed emitted photoelectron angle related variable in [log(e-12),log(2)].*/
        int     fNumSamplingAngles;
        
        //---------------------------------------------
        /** @brief Minimum of the gamma kinetic energy grid (default 1.0 [eV]).*/
        double  fMinPrimEnergy;
        
        //---------------------------------------------
        /** @brief Maximum of the gamma kinetic energy grid (default 1.0 [GeV]). After this threshold the photoelectron is considered to go straight with the same direction as the incident photon */
        double  fMaxPrimEnergy;
        
        //---------------------------------------------
        /** @brief Logarithm of SauterGavrilaPhotoElectricModel::fMinGammaEnergy i.e. ln(fMinGammaEnergy). */
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
         *  of Walker's alias sampling and linear approximation. At most we have as many data structure as
         *  SauterGavrilaPhotoElectricModel::fNumSamplingGammaEnergies
         *  and these data structure pointers are stored in the SauterGavrilaPhotoElectricModel::fAliasData linear array.
         */
        struct LinAlias{
            /** @brief Number of data points i.e. size of the arrays = SauterGavrilaPhotoElectricModel::fNumSamplingGammaEnergies. */
            int     fNumdata;
            
            /** @brief Gamma discrete kinetic energies stored in the alias table -  (transformed or not..it depends). */
            double *fXdata;
            
            /** @brief The probability density function values (not necessarily normalised) over the photoelectron angle
             *        variable values. */
            double *fYdata;
            
            /** @brief The alias probabilities (not necessarily normalised) over the photoelectron angle
             *        variable values.
             */
            double *fAliasW;
            
            /** @brief The alias indices over the photon energy variable values. */
            int    *fAliasIndx; // Alias indices
        };
        
        //---------------------------------------------
        /** @brief Linear array to store pointers to LinAlias data structures.
         *
         * The size is
         * SauterGavrilaPhotoElectricModel::fNumSamplingPrimEnergies.
         */
        LinAlias   **fAliasData;                   //alias data structure
        
        //---------------------------------------------
        /** @brief An alias sampler used at run-time to sample the emitted e- (photoelectron) direction. */
        AliasTable  *fAliasSampler;
        
    };
    
}      // namespace geantphysics

#endif // SAUTERGAVRILAPHOTOELECTRICMODEL_H
