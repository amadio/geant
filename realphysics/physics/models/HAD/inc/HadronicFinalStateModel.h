//====================================================================================================================
// HadronicFinalStateModel.h - Geant-V Prototype
/**
 * @file HadronicFinalStateModel.h
 * @brief Base class for all Geant-V hadronic final-state models.
 *
 * Similarly to the Geant4 G4HadronicInteraction class, this is the abstract base class for any hadronic final-state
 * model, for one or more projectile particle types.
 *
 * @author Mihaly Novak & Alberto Ribon (Aug 2016)
 */
//====================================================================================================================

#ifndef HADRONIC_FINAL_STATE_MODEL
#define HADRONIC_FINAL_STATE_MODEL

#include <string>
#include <vector>

#include "GeantTaskData.h"


// Forward declarations
namespace geantphysics {
  inline namespace GEANT_IMPL_NAMESPACE {
    class Isotope;
    class Material;
    class Element;
  }
}

namespace geantphysics {

class LightTrack;

/** Hadronic final-state model types */
enum class HadronicModelType {  
             kNotDefined,        /** it does not fit any other type: default */
             kString,            /** high-energy hadronic string model */
             kRescattering,      /** re-scattering of the hadronic string model */
             kCascade,           /** intranuclear cascade model */
             kPreEquilibrium,    /** pre-equilibrium model */ 
             kDeexcitation,      /** de-excitation model */
             kRadioactiveDecay,  /** radioactive decay model */
             kParticleHP,        /** low-energy high-precision model */
             kLeptoNuclear,      /** lepto-nuclear model */
             kQuasiElastic,      /** quasi-elastic model */
             kElastic,           /** elastic model */
             kUserDefined        /** generic, user-defined hadronic model */
           };

/**
 * @brief Class HadronicFinalStateModel
 */
class HadronicFinalStateModel {
public:
  /** @brief HadronicFinalStateModel default constructor */
  HadronicFinalStateModel();

  /** @brief HadronicFinalStateModel complete constructor */
  HadronicFinalStateModel( const std::string name, const std::vector< int > &projectilecodevec,
                           const HadronicModelType type, const double minenergy, const double maxenergy,
                           const double mintargetz, const double maxtargetz, 
                           const double mintargetn, const double maxtargetn );

  /** @brief HadronicFinalStateModel copy constructor */
  HadronicFinalStateModel( const HadronicFinalStateModel &other );

  /** @brief Operator = */
  HadronicFinalStateModel& operator=( const HadronicFinalStateModel &other );

  /** @brief HadronicFinalStateModel destructor */
  virtual ~HadronicFinalStateModel();

  /** @brief Method that is called once only at initialization. 
   *  SIGNATURE (INPUT PARAMETERS AND RETURN TYPE) AND ACTION TO BE DEFINED EVENTUALLY LATER...
   */
  virtual void Initialize( /* Not yet defined */ );

  /** @brief Method that returns "true" if the model is applicable, "false" otherwise.
   *  
   *  @param projectilecode is the GV particle code of the projectile
   *  @param projectilekineticenergy is the projectile kinetic energy in GeV
   *  @param targetisotope is the pointer to the target isotope
   */
  virtual bool IsApplicable( const int projectilecode, const double projectilekineticenergy, 
                             const Isotope* targetisotope );

  /** @brief Main method of this class: sampling of the final state.
   *  
   *  TO-BE-COMPLETED
   *
   *  @param xxx
   */

  virtual int SampleFinalState(LightTrack &track, Isotope* targetisotope, Geant::GeantTaskData *td) = 0;

  //--- Getters ---

  /** Method that returns the vector of GV particle codes for the allowed projectiles */
  const std::vector< int >& GetProjectileCodeVec() const { return fProjectileCodeVec; }

  /** Method that returns the name of this hadronic cross section */
  std::string GetName() const { return fName; }

  /** Method that returns the type of this hadronic final-state model */
  HadronicModelType GetType() const { return fType; }

  /** Method that returns the minimum required projectile kinetic energy in GeV */
  double GetLowEnergyUsageLimit() const { return fLowEnergyUsageLimit; }

  /** Method that returns the maximum required projectile kinetic energy in GeV */
  double GetHighEnergyUsageLimit() const { return fHighEnergyUsageLimit; }

  /** Method that returns the minimum required target Z (atomic number) */
  int GetMinTargetZ() const { return fMinTargetZ; }

  /** Method that returns the maximum required target Z (atomic number) */
  int GetMaxTargetZ() const { return fMaxTargetZ; }

  /** Method that returns the minimum required target N (number of nucleons) */
  int GetMinTargetN() const { return fMinTargetN; }

  /** Method that returns the maximum required target N (number of nucleons) */
  int GetMaxTargetN() const { return fMaxTargetN; }

  //--- Setters ---

  /** Method that sets the GV particle codes of the allowed projectiles */
  void SetProjectileCodeVec( const std::vector< int > &projectileCodeVec ) { 
    fProjectileCodeVec.clear();
    for ( unsigned int i = 0; i < projectileCodeVec.size(); i++ ) {
      fProjectileCodeVec.push_back( projectileCodeVec[i] );
    } 
  }

  /** Method that sets the name of this model */
  void SetName( const std::string &name ) { fName = name; }

  /** Method that sets the type of this model */
  void SetType( const HadronicModelType type ) { fType = type; }

  /** Method that sets the minimum required projectile kinetic energy in GeV */
  void SetLowEnergyUsageLimit( const double minenergy ) { fLowEnergyUsageLimit = minenergy; }

  /** Method that sets the maximum required projectile kinetic energy in GeV */
  void SetHighEnergyUsageLimit( const double maxenergy ) { fHighEnergyUsageLimit = maxenergy; }

  /** Method that sets the minimum required target Z (atomic number) */
  void SetMinTargetZ( const int mintargetz ) { fMinTargetZ = mintargetz; }

  /** Method that sets the maximum required target Z (atomic number) */
  void SetMaxTargetZ( const int maxtargetz ) { fMaxTargetZ = maxtargetz; }

  /** Method that sets the minimum required target N (number of nucleons) */
  void SetMinTargetN( const int mintargetn ) { fMinTargetN = mintargetn; }

  /** Method that sets the maximum required target N (number of nucleons) */
  void SetMaxTargetN( const int maxtargetn ) { fMaxTargetN = maxtargetn; }

private:
  std::vector< int > fProjectileCodeVec;  /** Vector of GV particle codes for the allowed projectiles */
  std::string fName;                      /** Model name */
  HadronicModelType fType;                /** Type of this model */
  double fLowEnergyUsageLimit;                      /** Minimum required projectile kinetic energy in GeV */
  double fHighEnergyUsageLimit;                      /** Maximum required projectile kinetic energy in GeV */
  int fMinTargetZ;                        /** Minimum required target Z (atomic number) */
  int fMaxTargetZ;                        /** Maximum required target Z (atomic number) */
  int fMinTargetN;                        /** Minimum required target N (number of nucleons) */
  int fMaxTargetN;                        /** Maximum required target N (number of nucleons) */
};

}  // end of namespace geantphysics

#endif
