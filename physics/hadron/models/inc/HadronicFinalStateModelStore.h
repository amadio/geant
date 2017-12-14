//====================================================================================================================
// HadronicFinalStateModelStore.h - Geant-V Prototype
/**
 * @file HadronicFinalStateModelStore.h
 * @brief Class for all Geant-V hadronic final-state models.
 *
 * This class keeps a collection of hadronic final-state models - as pointers to objects derived from the abstract
 * base class HadronicFinalState.
 * The ordering in which these models are registered in the store does not matter.
 * Overlapping in the projectile kinetic energy between models is allowed, but the following restrictions:
 * -  not more than 2 models can overlap in the same interval;
 * -  two models cannot overlap in their full applicability range.
 * When two models overlap in an energy interval, at each interaction (i.e. call to the model's SampleFinalState
 * method) one model is chosen with the following rule:
 * -  a random number is thrown;
 * -  the probability to choose one of the two models is a linear function of the (projectile kinetic) energy,
 *    between 0.0 at the point where the overlaps starts, and 1.0 where the overlap ends.
 *
 * @author Mihaly Novak & Alberto Ribon (Aug 2016)
 */
//====================================================================================================================

#ifndef HADRONIC_FINAL_STATE_MODEL_STORE
#define HADRONIC_FINAL_STATE_MODEL_STORE

#include <string>
#include <vector>
#include "Geant/Config.h"

// Forward declarations

namespace geantphysics {
    inline namespace GEANT_IMPL_NAMESPACE {
    class Isotope;
    }
  }


namespace geantphysics {

  
class HadronicFinalStateModel;


/**
 * @brief Class HadronicFinalStateModelStore
 */
class HadronicFinalStateModelStore {
public:
  /** @brief HadronicFinalStateModelStore default constructor */
  HadronicFinalStateModelStore();

  /** @brief HadronicFinalStateModelStore complete constructor */
  HadronicFinalStateModelStore( const std::string name );

  /** @brief HadronicFinalStateModelStore copy constructor */
  HadronicFinalStateModelStore( const HadronicFinalStateModelStore &other );

  /** @brief Operator = */
  HadronicFinalStateModelStore& operator=( const HadronicFinalStateModelStore &other );

  /** @brief HadronicFinalStateModelStore destructor */
  ~HadronicFinalStateModelStore();

  /** @brief Method that is called once only at initialization. 
   *  SIGNATURE (INPUT PARAMETERS AND RETURN TYPE) AND ACTION TO BE DEFINED EVENTUALLY LATER...
   */
  void Initialize( /* Not yet defined */ );

  /** @brief Method that register a hadronic final-state model to the store. 
   *
   *  Note that the ordering in which hadronic final-state models are registered does not matter: ...
   *
   *  @param ptrhadfs is a pointer to a HadronicFinalState object
  */
  void RegisterHadronicFinalStateModel( HadronicFinalStateModel* ptrhadfs );

  /** @brief Method that returns the index of the chosen hadronic final-state model. If none, returns -1 .
   *  
   *  @param projectilecode is the GV particle code of the projectile
   *  @param projectilekineticenergy is the projectile kinetic energy in GeV
   *  @param targetisotope is the pointer to the target isotope
   */
  int GetIndexChosenFinalStateModel( const int projectilecode, const double projectilekineticenergy,
                                     const Isotope* targetisotope ) const;

  //--- Getters ---

  /** Method that returns the vector of HadronicFinalStates */
  const std::vector< HadronicFinalStateModel* >& GetHadronicFinalStateModelVec() const { return fHadFsVec; }

  /** Method that returns the name of this hadronic final-state model store */
  std::string GetName() const { return fName; }

  //--- Setters ---

  /** Method that sets the name of this hadronic cross section */
  void SetName( const std::string &name ) { fName = name; }

private:
  std::vector< HadronicFinalStateModel* > fHadFsVec;  /* Vector of hadronic final-state models */
  std::string fName;                                  /* Hadronic final-state model store name */
};

}  // end of namespace geant

#endif
