//====================================================================================================================
// HadronicCrossSection.h - Geant-V Prototype
/**
 * @file HadronicCrossSection.h
 * @brief Base class for all Geant-V microscopic hadronic cross sections.
 *
 * This abstract class specifies one microscopic, i.e. per-isotope, hadronic cross section (e.g. elastic or inelastic) 
 * for one or more projectile particle types.
 * This is the (abstract) base class from which all microscopic hadronic cross sections must derived from.
 * This class does NOT cover the case of thermal neutron cross sections, where parameters like temperature or element
 * needs to be provided (a separate base class will be provided for that case).
 *
 * @author Mihaly Novak & Alberto Ribon (Aug 2016)
 */
//====================================================================================================================

#ifndef HADRONIC_CROSS_SECTION
#define HADRONIC_CROSS_SECTION

#include <string>
#include <vector>


namespace geantphysics {

/**
  * @brief Class HadronicCrossSection
  */
class HadronicCrossSection {

 public:
  /** @brief HadronicCrossSection default constructor */
  HadronicCrossSection();

  /** @brief HadronicCrossSection complete constructor */
  HadronicCrossSection( const std::string name, const std::vector< int > &projectilecodevec, 
                        const double minenergy, const double maxenergy,
                        const double mintargetz, const double maxtargetz, 
                        const double mintargetn, const double maxtargetn );

  /** @brief HadronicCrossSection copy constructor */
  HadronicCrossSection( const HadronicCrossSection &other );

  /** @brief Operator = */
  HadronicCrossSection& operator=( const HadronicCrossSection &other );

  /** @brief HadronicCrossSection destructor */
  virtual ~HadronicCrossSection();

  /** @brief Method that is called once only at initialization. 
   *  SIGNATURE (INPUT PARAMETERS AND RETURN TYPE) AND ACTION TO BE 
   *  DEFINED EVENTUALLY LATER...
   */
  virtual void Initialize( /* Not yet defined */ );

  /** @brief Method that returns "true" if the cross section is applicable, "false" otherwise.
   *  
   *  (Witek) I am not sure we will need that method, but I leave it for the time being
   *
   *  To be applicable, the projectile type should be one of the allowed values, and the projectile kinetic energy,
   *  target atomic number (Z) and target number of nucleons (N) should all be within the allowed ranges.
   *  In order to be able to treat also the special case of thermal neutrons, we could introduce the temperature 
   *  is an optional parameter, but mostly likely that will be handled by a separate cross-section base class.
   *  
   *  @param projectilecode is the GV particle code of the projectile
   *  @param projectilekineticenergy is the projectile kinetic energy in GeV
   *  @param Z is the atomic number 
   *  @param N is the number of nucleons
   */
  virtual bool IsApplicable( const int projectilecode, const double projectilekineticenergy,
                             const int Z, const int N);

  /** @brief Method that returns the isotope hadronic cross section (unit: cm^2).
   * 
   *  The GetIsotopeCrossSection is a pure virtual method, so a definition must be provided by each class
   *  inheriting from this abstract class.
   *
   *  The current decision is NOT to provide a method to return element cross-sections. This would be only needed
   *  for the case of thermal neutrons, and we think it should be handled by a separate cross-section base class. 
   *
   *  @param projectilecode is the GV particle code of the projectile
   *  @param projectilekineticenergy is the projectile kinetic energy in GeV
   *  @param projectilemass is the projectile mass in GeV
   *  @param Z is the atomic number of the target
   *  @param N is the number of nucleons of the target
   */
  virtual double GetIsotopeCrossSection( const int projectilecode, const double projectilekineticenergy,
                                         const double projectilemass, const int Z, const int N ) = 0;

  //--- Getters ---

  /** Method that returns the vector of GV particle codes for the allowed projectiles */
  const std::vector< int >& GetProjectileCodeVec() const { return fProjectileCodeVec; }

  /** Method that returns the name of this hadronic cross section */
  std::string GetName() const { return fName; }

  /** Method that returns the minimum required projectile kinetic energy in GeV */
  double GetMinEnergy() const { return fMinEnergy; }

  /** Method that returns the maximum required projectile kinetic energy in GeV */
  double GetMaxEnergy() const { return fMaxEnergy; }

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
    for ( unsigned long i = 0; i < projectileCodeVec.size(); i++ ) {
      fProjectileCodeVec.push_back( projectileCodeVec[i] );
    } 
  }

  /** Method that sets the name of this hadronic cross section */
  void SetName( const std::string &name ) { fName = name; }

  /** Method that sets the minimum required projectile kinetic energy in GeV */
  void SetMinEnergy( const double minenergy ) { fMinEnergy = minenergy; }

  /** Method that sets the maximum required projectile kinetic energy in GeV */
  void SetMaxEnergy( const double maxenergy ) { fMaxEnergy = maxenergy; }

  /** Method that sets the minimum required target Z (atomic number) */
  void SetMinTargetZ( const int mintargetz ) { fMinTargetZ = mintargetz; }

  /** Method that sets the maximum required target Z (atomic number) */
  void SetMaxTargetZ( const int maxtargetz ) { fMaxTargetZ = maxtargetz; }

  /** Method that sets the minimum required target N (number of nucleons) */
  void SetMinTargetN( const int mintargetn ) { fMinTargetN = mintargetn; }

  /** Method that sets the maximum required target N (number of nucleons) */
  void SetMaxTargetN( const int maxtargetn ) { fMaxTargetN = maxtargetn; }

private:
  std::vector< int > fProjectileCodeVec;  /* Vector of GV particle codes for the allowed projectiles */
  std::string fName;                      /* Hadronic cross section name */
  double fMinEnergy;                      /* Minimum required projectile kinetic energy in GeV */
  double fMaxEnergy;                      /* Maximum required projectile kinetic energy in GeV */
  int fMinTargetZ;                        /* Minimum required target Z (atomic number) */
  int fMaxTargetZ;                        /* Maximum required target Z (atomic number) */
  int fMinTargetN;                        /* Minimum required target N (number of nucleons) */
  int fMaxTargetN;                        /* Maximum required target N (number of nucleons) */
};

}  // end of namespace geantphysics

#endif
