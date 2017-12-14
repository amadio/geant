//====================================================================================================================
// HadronicCrossSectionStore.h - Geant-V Prototype
/**
 * @file HadronicCrossSectionStore.h
 * @brief Class for all Geant-V hadronic process cross sections.
 *
 * Similarly to the Geant4 G4CrossSectionDataStore class, this class specifies a cross section 
 * (e.g. elastic or inelastic) for a hadronic process. More hadronic cross sections - as objects derived
 * from the abstract base class HadronicCrossSection - can be combined, with the convention (as in Geant4)
 * of "last-in-first-out", i.e. the cross sections are checked starting from the back of the list, and the
 * first one which is applicable will be used.
 *
 * Two types of cross sections are provided by this class:
 * - "microscopic" cross sections :
 *   o  per-isotope cross section
 *   o  per-element cross section, i.e. the sum of the per-isotope cross sections weighted with the relative
 *                  isotopic abundance composition of the element;
 * - "macroscopic" cross section, i.e. per-material cross section multiplied by the number of atoms
 *                per unit of volume of the target. This has the dimension of 1/length.
 *
 * These cross sections are used at initialization, when the lambda tables are built.
 * In the event simulation, the sampling of the target isotope is needed: this can be done on-the-fly, by
 * calling the "microscopic" cross sections, or, better, in a much faster way by using some special tables
 * built during the initialization. Currently the first method is used, but later the second one should be tried out.
 * 
 * @author Mihaly Novak & Alberto Ribon (Aug 2016)
 */
//====================================================================================================================

#ifndef HADRONIC_CROSS_SECTION_STORE
#define HADRONIC_CROSS_SECTION_STORE

#include <string>
#include <vector>

#include "Geant/Config.h"

namespace geantphysics {

// Forward declarations
class HadronicCrossSection;

    inline namespace GEANT_IMPL_NAMESPACE {
    class Isotope;
    class Material;
    class Element;
    }

/**
 * @brief Class HadronicCrossSectionStore
 */
class HadronicCrossSectionStore {
public:
  /** @brief HadronicCrossSectionStore default constructor */
  HadronicCrossSectionStore();

  /** @brief HadronicCrossSectionStore complete constructor */
  HadronicCrossSectionStore( const std::string name );

  /** @brief HadronicCrossSectionStore copy constructor */
  HadronicCrossSectionStore( const HadronicCrossSectionStore &other );

  /** @brief Operator = */
  HadronicCrossSectionStore& operator=( const HadronicCrossSectionStore &other );

  /** @brief HadronicCrossSectionStore destructor */
  ~HadronicCrossSectionStore();

  /** @brief Method that is called once only at initialization. 
   *  SIGNATURE (INPUT PARAMETERS AND RETURN TYPE) AND ACTION TO BE 
   *  DEFINED EVENTUALLY LATER...
   */
  void Initialize( /* Not yet defined */ );

  /** @brief Method that register a hadronic cross section to the store. 
   *
   *  Note that the ordering in which hadronic cross sections are registered is crucial: the later a cross section
   *  is registered, the higher is the priority of that cross sections (i.e. LIFO (Last In First Out) structure).
   *
   *  @param ptrhadxsec is a pointer to a HadronicCrossSection object
  */
  void RegisterHadronicCrossSection( HadronicCrossSection* ptrhadxsec );

  /** @brief Method that returns the index of the applicable cross section that has been registered last. 
   *         If none, returns -1.
   *  
   *  @param projectilecode is the GV particle code of the projectile
   *  @param projectilekineticenergy is the projectile kinetic energy in GeV
   *  @param targetelement is the pointer to the target element
   *  @param targetmaterial is the pointer to the target material (made of the single element)
   */
  int GetIndexFirstApplicableXsec( const int projectilecode, const double projectilekineticenergy,
                                   const Element* targetelement, const Material* targetmaterial );

  /** @brief Methods that return, respectively, the isotope, element, and macroscopic hadronic cross sections.
   *
   *  The first two, isotope and element hadronic cross sections are microscopic ones, with units of cm^2 .
   *  The third one is the macroscopic, material hadronic cross section, i.e. the cross section in the material
   *  (i.e. sum of the element cross sections weighted with their relative abundances) multiplied by the number
   *  of atoms per unit of volume; its unit is 1/cm.
   *
   *  @param projectilecode is the GV particle code of the projectile
   *  @param projectilekineticenergy is the projectile kinetic energy in GeV
   *  @param projectilemass is the projectile mass in GeV
   *  @param Z is the atomic number of the target
   *  @param N is the number of nucleons of the target
   */
  double GetIsotopeCrossSection( const int projectilecode, const double projectilekineticenergy, const double projectilemass,
				 const Isotope* targetisotope, const Element* targetelement, const Material* targetmaterial);

  double GetElementCrossSection( const int projectilecode, const double projectilekineticenergy, const double projectilemass,
                                 const Element* targetelement, const Material* targetmaterial );

  double GetMacroscopicCrossSection( const int projectilecode, const double projectilekineticenergy, const double projectilemass,
                                     const Material* targetmaterial );

  /** @brief Method that returns a target nucleus, as a pair of Z (atomic number) and N (number of nucleons). 
   *  
   *  To sample a nucleus in the target material, ...
   *
   *  @param projectilecode is the GV particle code of the projectile
   *  @param projectilekineticenergy is the projectile kinetic energy in GeV
   *  @param targetmaterial is the pointer to the target material
   */
  std::pair< int, int > SampleTarget( const int projectilecode, const double projectilekineticenergy, const double projectilemass,
                                      const Material* targetmaterial );

  //--- Getters ---

  /** Method that returns the vector of HadronicCrossSections */
  std::vector< HadronicCrossSection* >& GetHadronicCrossSectionVec() { return fHadXsecVec; }

  /** Method that returns the name of this hadronic cross section */
  std::string GetName() const { return fName; }

  //--- Setters ---

  /** Method that sets the name of this hadronic cross section */
  void SetName( const std::string &name ) { fName = name; }

private:
  std::vector< HadronicCrossSection* > fHadXsecVec;  /* Vector of hadronic cross sections */
  std::string fName;                                 /* Hadronic cross section store name */
};

}  // end of namespace geant

#endif
