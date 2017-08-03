//====================================================================================================================
// HadronicProcess.h - Geant-V Prototype
/**
 * @file HadronicProcess.h
 * @brief Base class for all hadronic physics processes in the Geant-V prototype.
 *
 * Similarly to the Geant4 class G4HadronicProcess, the class HadronicProcess is the base class for all
 * hadronic physics processes. 
 *
 * @author Mihaly Novak & Alberto Ribon (Aug 2016)
 */
//====================================================================================================================

#ifndef HADRONIC_PROCESS
#define HADRONIC_PROCESS

#include <string>
#include <vector>
#include "PhysicsProcess.h"


namespace geantphysics {

// Forward declarations
class LightTrack;
class LightTrack_v;
class HadronicCrossSection;
class HadronicCrossSectionStore;
class HadronicFinalStateModelStore;
class HadronicFinalStateModel;
 
 inline namespace GEANT_IMPL_NAMESPACE {
   class Isotope;
   class Material;
   class Element;
 }
  

/** Hadronic process types */
enum class HadronicProcessType {  
             kNotDefined,        /** it does not fit any other type: default */
             kElastic,           /** elastic process */
             kQuasiElastic,      /** quasi-elastic process */
             kCapture,           /** nuclear capture process */
             kFission,           /** nuclear fission process */
             kInelastic,         /** any other inelastic process */
             kRadioactiveDecay,  /** radioactive decay process */
             kLeptoNuclear,      /** lepto-nuclear process */
             kUserDefined        /** generic, user-defined hadronic process */
           };

/**
 * @brief Class HadronicProcess
 */
class HadronicProcess : public PhysicsProcess {
public:
  /** @brief HadronicProcess default constructor */
  HadronicProcess();

  HadronicProcess( const std::string &name );

  /** @brief HadronicProcess complete constructor */
  HadronicProcess( const std::string &name, const std::vector< int > &particlecodevec, 
                   const HadronicProcessType type, const bool isatrest,
                   HadronicCrossSectionStore* xsecstore, HadronicFinalStateModelStore* modelstore );

  /** @brief HadronicProcess destructor */
  virtual ~HadronicProcess();

  // The methods below are those inherited from PhysicsProcess
 
  /** Method that returns "true" ("false") if the specified GV particle code is (not) accepted by this process */
  virtual bool IsApplicable( const LightTrack &/*track*/ ) const;

  virtual double ComputeMacroscopicXSection(const MaterialCuts * /*matcut*/, double /*kinenergy*/,
                                            const Particle * /*particle*/, double /*mass*/) const;

  /** Method that returns the atomic cross section, for an in-flight hadronic process */
  virtual double GetAtomicCrossSection( const int particlecode, const double particlekineticenergy, const double particlemass,
                                        const Element* targetelement, const Material* targetmaterial ) const;

  /** Method that returns the mean lifetime for an at-rest hadronic process */
  //  virtual double AverageLifetime( const LightTrack &track ) const;

  /** Method that returns the sampled target isotope for an in-flight hadronic process.
      Z and N are also set in the track. It should be called only by the PostStepDoIt method. */
  Isotope* SampleTarget( LightTrack &track ) const;

  /** Main method that produces the secondaries for an in-flight hadronic process */
  virtual int PostStepDoIt( LightTrack &track, Geant::GeantTaskData *td);

  /** Main method that sample the target isotope and produces the secondaries for an at-rest hadronic process */
  virtual void AtRestDoIt( LightTrack &track, Geant::GeantTaskData * td);

  /** Method to add model to the process **/
  void AddModel(HadronicFinalStateModel *model);

  /** Method to add cross sections to the process **/
  void AddCrossSection(HadronicCrossSection *xsection);

  //--- Getters ---

  /** Method that returns the vector of GV particle codes for the allowed particles */
  const std::vector< int >& GetParticleCodeVec() const { return fParticleCodeVec; }

  /** Method that returns the type of this hadronic process */
  HadronicProcessType GetType() const { return fType; }

  /** Method that returns the pointer to the hadronic cross section store */
  const HadronicCrossSectionStore* GetCrossSectionStore() const { return fXsecStore; }

  /** Method that returns the pointer to the hadronic final-state model store */
  const HadronicFinalStateModelStore* GetFinalStateModelStore() const { return fModelStore; }

  //--- Setters ---

  /** Method that sets the GV particle codes of the allowed particles */
  void SetParticleCodeVec( const std::vector< int > &particlecodevec ) { 
    fParticleCodeVec.clear();
    for ( size_t i = 0; i < particlecodevec.size(); i++ ) {
      fParticleCodeVec.push_back( particlecodevec[i] );
    } 
  }

  /** Method that sets the type of this process */
  void SetType( const HadronicProcessType type ) { fType = type; }

  /** Method that sets the pointer to the hadronic cross section store */
  void SetCrossSectionStoreStore( HadronicCrossSectionStore* xsecstore ) { fXsecStore = xsecstore; }

  /** Method that returns the pointer to the hadronic final-state model store */
  void SetFinalStateModelStoreStore( HadronicFinalStateModelStore* modelstore ) { fModelStore = modelstore; }
  
private:
  std::vector< int > fParticleCodeVec;        /** Vector of GV particle codes for the allowed particles */
  HadronicProcessType fType;                  /** Type of this hadronic process */
  HadronicCrossSectionStore* fXsecStore;      /** Set of hadronic cross sections for this process */
  HadronicFinalStateModelStore* fModelStore;  /** Set of hadronic final-state models for this process */
};

}  // end of namespace geant

#endif
