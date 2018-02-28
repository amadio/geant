
#ifndef ELEMENTPROPERTIES_H
#define ELEMENTPROPERTIES_H

#include "Geant/Config.h"

namespace geantphysics {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief   Class to represent some frequently used element/atom related variables.
 * @class   ElementProperties
 * @author  M Novak, A Ribon
 * @date    May 2017
 */

class Element;

class ElementProperties {
public:
  /**
  * @name Constructors/destructors:
  */
  //@{
  /** @brief Constructor to create element properties object.
   *
   *  Initialisation of parameters are done at construction time.
   *
   *  @param[in] Pointer to the corresponding element object.
   */
  ElementProperties(Element *elem);
  /** @brief Destructor.
   */
  ~ElementProperties();
  //@}

  /**
  * @name Public methods to obtain simple elemnt properties:
  */
  //@{
  /** @brief \f$ Z^{1/3} \f$ */
  double GetZ13() const { return fZ13; }
  /** @brief \f$ Z^{2/3} \f$ */
  double GetZ23() const { return fZ23; }
  /** @brief \f$ \log(Z) \f$ */
  double GetLogZ() const { return fLogZ; }
  /** @brief \f$ \frac{1}{3} \log(Z) \f$ */
  double GetLogZ13() const { return fLogZ13; }
  /** @brief \f$ \frac{2}{3} \log(Z) \f$ */
  double GetLogZ23() const { return fLogZ23; }
  /** @brief Coulomb correction */
  double GetCoulombCorrection() const { return fCoulombCorrection; }
  //@}

private:
  /** @brief Internal method to compute some simple element properties. */
  void InitialiseMembers();

  /** @brief Internal method to compute Coulomb correction.
    *
    * param[in] z Atomic number.
    */
  void ComputeCoulombCorrection(double z);

private:
  /** @brief Pointer to the element that this element parameters belong to (doesn't own the object)*/
  Element *fElement;

  /** @brief \f$ Z^{1/3} \f$ */
  double fZ13;
  /** @brief \f$ Z^{2/3} \f$ */
  double fZ23;
  /** @brief \f$ \log(Z) \f$ */
  double fLogZ;
  /** @brief \f$ \frac{1}{3} \log(Z) \f$ */
  double fLogZ13;
  /** @brief \f$ \frac{2}{3} \log(Z) \f$ */
  double fLogZ23;
  /** @brief Coulomb correction */
  double fCoulombCorrection;
};

} // inline namespace GEANT_IMPL_NAMESPACE {
} // namespace geantphysics

#endif // ELEMENTPROPERTIES_H
