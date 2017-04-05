
#ifndef EMMODELMANAGER_H
#define EMMODELMANAGER_H

#include <vector>
#include <string>

namespace geantphysics {

class EMModel;
class Particle;

/**
 * @brief   Class to handle EMModels that has been added to an EMPhysicsProcess.
 * @class   EMModelManager
 * @author  M Novak, A Ribon
 * @date    august 2016
 *
 * Each EMPhysicsProcess has an EMModelManager member. EMModels can be added and accessed through this EMModelManager
 * member of the EMPhysicsProcess. The EMModelManager object owns all the models that have been added to it.
 *
 * At the initialization of the EMModelManager the EMModelManager: (1) will set for each EMModel-s that have been added
 * to the corresponding EMPhysicsProcess the list of active regions: (a) first the default active regions which are
 * determined by the active regions of the EMPhysicsProcess; (b) on the top of this, the user requested inactive regions
 * are considered in case of each individual EMModel-s; (2) the EMModelManager will prepare the list of active EMModel-s
 * per region collections. (3) After each EMModel-s active regions list is set, the EMModelManager will initilise all
 * EMModel-s: each EMModel-s will be initilized only for reagions where they are active.
 */


class EMModelManager {
public:
  /**
   * @brief CTR
   */
  EMModelManager();

  /**
   * @brief DTR
   */
  ~EMModelManager();

  // must be called after all the models are added i.e. when the process is initialised
  // will init the models as well
  /**
   * @brief Method to initilize the model manager.
   *
   * Must be called after all the models are added. It is called from the EMPhysicsProcess Initialize() method.
   * - will set for each EMModel-s, that have been added to the corresponding EMPhysicsProcess, the list of active
   *   regions:
   *    - first the default active regions which are determined by the active regions of the EMPhysicsProcess;
   *    - on the top of this, the user requested inactive regions are considered in case of each individual EMModel-s;
   * - the EMModelManager will prepare the list of active EMModel-s per region collections.
   * - after each EMModel-s active regions list is set, the EMModelManager will initilise all EMModel-s: each EMModel-s
   *   will be initilized only for reagions where they are active.
   *
   * @param[in] listofisactiveregions  List of active region vector (that is the same as the corresponding
   *                                   EMPhysicsProcess list of active region).
   * @param[in] procname               Name of the EMPhysicsProcess to provide information in case of warnings.
   * @param[in] particlelist           List of particles the EMPhysicsProcess to provide information in case of warnings.
   */
  void Initialise(std::vector<bool> &listofisactiveregions, const std::string &procname,
                  const std::vector<Particle*> &particlelist);


  /**
   * @brief Method to add EMModel to the model manager.
   *
   * Expected to be called from the EMPhysicsProcess::AddModel() method. Each EMModel pointer that has been added to
   * this model manager are stored in a vector. The index of the model will be set to its index in this vector.
   *
   * @param[in] mod  Pointer to an EMModel object to be added to the model manager.
   * @return    The index of the EMModel pointer in the container.
   */
  int AddModel(EMModel *mod);

  // selects one EMModel that is active in the region and the given kinetic energy of
  // the particle is within the usage limits of the model. If the kinetic energy is lower than the lowest usgae limits
  // of the models it will return with the lowest energy model. So we always check in the model if we are within the
  // usage limits of the model in the final state sampling
  EMModel* SelectModel(double ekin, int regionindx);


  /**
   * @brief Public method to obtain the list of EMModels from the manager that are active in a given region.
   *
   * @param[in] regionindx  Index of the region in which the list of active models is required. (validity of the index
   *                        is not checked)
   * @return    List of EMModel-s that are active in the given region.
   */
  const std::vector<EMModel*>& GetModelListInRegion(int regionindx) const { return fModelListsPerRegions[regionindx]; }


  /**
   * @brief Public method to obtain the collections of active EMModel-s per region.
   *
   * The first dimention in the returned collection is the region index and the list of EMModel-s that are active in
   * that region is stored along the second dimension.
   *
   * @return    Collections of active EMModel-s per region.
   */
  const std::vector< std::vector< EMModel*> >& GetModelListPerRegions() const { return fModelListsPerRegions; }

  /**
   * @brief Public method to get the number of EMModel-s that have been added so far to thid mamanger.
   *
   * @return Number of EMModel-s that have been added so far to thid mamanger.
   */
  int GetNumberOfModels() const { return fNumModels; }


  /**
   * @brief Method to obtain one EMModel by its index.
   *
   * Each EMModel that has been added to the model manager is stored in a vector. The index of the model is set to its
   * index in this vector.
   *
   * @param[in]  modelindx  Index of the required EMModel.
   * @return     Pointer to the EMModel with the specified index.
   */
  EMModel* GetModelByIndex(int modelindx) { return fModelList[modelindx]; }


  /**
   * @brief Public method to obtain the full list of EMModel-s that have been added to this model manager so far.
   *
   * @return List of EMModel-s that have been added so far to this model manager.
   */
  const std::vector<EMModel*>& GetModelList() const { return fModelList; }


//
// data members
//
private:
  int                      fNumModels; /** number of models added to this manager */
  std::vector<EMModel*>    fModelList; /** list of models added to this manager; the manager owns the models; the model
                                           indices are set to the index of the model in this vector */
  std::vector< std::vector<EMModel*> >  fModelListsPerRegions;
//  std::vector<double>      fLowEnergyUsageLimits; /** the low energy usage limits for each model */

};

} // namespace geantphysics

#endif // EMMODELMANAGER_H
