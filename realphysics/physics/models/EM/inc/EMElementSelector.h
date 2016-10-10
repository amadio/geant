
#ifndef EMELEMENTSELECTOR_H
#define EMELEMENTSELECTOR_H

#include <vector>

namespace geantphysics {

class EMModel;
class MaterialCuts;
class Particle;

/**
 * @brief   Utility class for EMModel-s to provide common easy way to build at initialization and use at run-time target
 *          element selection for a given discrete interaction described by the given EMModel.
 * @class   EMElementSelector
 * @author  M Novak, A Ribon
 * @date    august 2016
 *
 * Each physics models that derived from EMModel base class has the possibility to use this functionality. The derived
 * models needs to call explicitly the EMModel base class InitialiseElementSelectors method in their Initialize method
 * at the end after the model has already been initialised properly.
 *
 * Each EMElementSelector object will be specific for one EMModel, one Material/MaterialCuts and for a given Particle.
 * The collection of EMElementSelector-s (for a given EMModel, for a given Particle, for all MaterialCuts/Material)
 * will be created, build, stored at initialization and handled at run-time by the EMModel base class.
 */


class EMElementSelector {
public:
  /**
   * @brief Contructor to create element selector for a given EMModel for one Material/MaterialCuts.
   *
   * EMElementSelector objects are created and built in the EMModel::InitialiseElementSelectors() method.
   *
   * @param[in] emmodel       Pointer to the EMModel that this element selector belongs to.
   * @param[in] emin          Minimum of the kinetic energy grid over this element selector must be built.
   * @param[in] emax          Maximum of the kinetic energy grid over this element selector must be built.
   * @param[in] binsperdecade Number of energy bins per decade of the kinetic energy grid.
   * @param[in] numelement    Number of elements in the Material this element selector belongs to.
   */
  EMElementSelector(EMModel *emmodel, double emin, double emax, int binsperdecade, int numelement);

  /**
   * @brief DTR.
   */
 ~EMElementSelector();


 /**
  * @brief Method to build up i.e. initilize the element selector.
  *
  * The method is expceted to be called from EMModel::InitialiseElementSelectors() after the EMElementSelector object
  * has been created.
  *
  * @param[in] matcut  Pointer to the MaterialCuts object which(or which Material) this element selector belongs to.
  * param[in]  part    Pointer to the Particle object the EMModel, that this element selector belongs to, assigned to.
  */
 void    Build(const MaterialCuts *matcut, const Particle *part);

  /**
   * @brief Run-time method to sample one traget element.
   *
   * This method is expected to be called from the EMModel::SampleTargetElementIndex() method.
   * Todo: will return with the selected const Element* later!
   *
   * @param[in] ekin  Kinetic energy of the particle at wich target element must be sampled for discrete interaction.
   * @param[in] rndm  Random number uniformly distributed on [0,1).
   */
  int     SampleTargetElement(double ekin, double rndm);

//
// private methods
//
private:
  /**
   * @brief Utility method to set up the kinetic energy grid over this element selector needs to be built.
   *
   * @param[in]  binsperdecade Number of energy bins per decade.
   */
  void InitializeEnergyGrid(int binsperdecade);

//
// data members
//
private:
  // data to describe the energy grid of this element selector
  int         fNumEnergyBins;
  double      fMinEnergy;
  double      fMaxEnergy;
  double      fLogMinEnergy;
  double      fEnergyILDelta;
  double     *fEnergyGrid;      // owned; size is fNumEnergyBins

  EMModel    *fEMModel;        // NOT owned

  std::vector<double*>  fProbsPerElements; // normalized cumulative probabilities per target elements for each energy
                                           // grid point, for each element; the vector will store #elements,
                                           // fNumEnergyBins size double arrays; the class do own these arrays
};

}  // namespace geantphysics

#endif //  EMELEMENTSELECTOR_H
