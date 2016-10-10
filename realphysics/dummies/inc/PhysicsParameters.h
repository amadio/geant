#ifndef PHYSICSPARAMETERS_H
#define PHYSICSPARAMETERS_H

#include <vector>


namespace geantphysics {
/**
 * @brief   Class to store some physics parameters.
 * @class   PhysicsParameters
 * @author  M Novak, A Ribon
 * @date    july 2016
 *
 * Each PhysicsList will have a PhysicsParameters object. The object will be created in the PhysicsList ctr but owned by
 * the PhysicsParameters class. Each PhysicsParameters object will store the list of its active regions.
 */
class PhysicsParameters{
public:
  PhysicsParameters();
 ~PhysicsParameters() { /*nothing to do at the moment*/}

  // delete all PhysicsParameters objects
  static void Clear();

  std::vector<bool>& GetListActiveRegions() { return fListActiveRegions; }
  bool IsActiveRegion(int regionindx) const { return fListActiveRegions[regionindx]; }

  // list of all physics parameters objects created so far
  static const std::vector<PhysicsParameters*>& GetThePhysicsParametersTable() { return gThePhysicsParametersTable; }
  // get physics parameters object that is active in a given region; used only at initialization
  static const PhysicsParameters*  GetPhysicsParametersForRegion(int regionindx);


  // secondary production kinetic energy threshold related minimum/maximum values
  static void   SetMinAllowedGammaCutEnergy   (double val) { fMinAllowedGammaCutEnergy    = val; }
  static void   SetMinAllowedElectronCutEnergy(double val) { fMinAllowedElectronCutEnergy = val; }
  static void   SetMinAllowedPositronCutEnergy(double val) { fMinAllowedPositronCutEnergy = val; }
  static void   SetMaxAllowedGammaCutEnergy   (double val) { fMaxAllowedGammaCutEnergy    = val; }
  static void   SetMaxAllowedElectronCutEnergy(double val) { fMaxAllowedElectronCutEnergy = val; }
  static void   SetMaxAllowedPositronCutEnergy(double val) { fMaxAllowedPositronCutEnergy = val; }

  static double GetMinAllowedGammaCutEnergy()    { return fMinAllowedGammaCutEnergy;    }
  static double GetMinAllowedElectronCutEnergy() { return fMinAllowedElectronCutEnergy; }
  static double GetMinAllowedPositronCutEnergy() { return fMinAllowedPositronCutEnergy; }
  static double GetMaxAllowedGammaCutEnergy()    { return fMaxAllowedGammaCutEnergy;    }
  static double GetMaxAllowedElectronCutEnergy() { return fMaxAllowedElectronCutEnergy; }
  static double GetMaxAllowedPositronCutEnergy() { return fMaxAllowedPositronCutEnergy; }

  static void   SetDefaultGammaCutInLength(double val)     { fDefaultGammaCutInLength    = val; }
  static void   SetDefaultElectronCutInLength(double val)  { fDefaultElectronCutInLength = val; }
  static void   SetDefaultPositronCutInLength(double val)  { fDefaultPositronCutInLength = val; }

  static double GetDefaultGammaCutInLength()     { return fDefaultGammaCutInLength; }
  static double GetDefaultElectronCutInLength()  { return fDefaultElectronCutInLength; }
  static double GetDefaultPositronCutInLength()  { return fDefaultPositronCutInLength; }

  static void   SetDefaultGammaCutInEnergy(double val)     { fDefaultGammaCutInEnergy    = val; }
  static void   SetDefaultElectronCutInEnergy(double val)  { fDefaultElectronCutInEnergy = val; }
  static void   SetDefaultPositronCutInEnergy(double val)  { fDefaultPositronCutInEnergy = val; }

  static double GetDefaultGammaCutInEnergy()     { return fDefaultGammaCutInEnergy; }
  static double GetDefaultElectronCutInEnergy()  { return fDefaultElectronCutInEnergy; }
  static double GetDefaultPositronCutInEnergy()  { return fDefaultPositronCutInEnergy; }

  //
  // EM loss table related variables
  void    SetMinLossTableEnergy(double val);
  double  GetMinLossTableEnergy() const { return fMinLossTableEnergy; }

  void    SetMaxLossTableEnergy(double val);
  double  GetMaxLossTableEnergy() const { return fMaxLossTableEnergy; }

  void    SetNumLossTableBins(int val);
  int     GetNumLossTableBins() const { return fNumLossTableBins; }

  void    SetNumLossTableBinsPerDecade(int val);
  int     GetNumLossTableBinsPerDecade()   const { return fNumLossTableBinsPerDecade; }

  void    SetIsComputeCSDARange(bool val) { fIsComputeCSDARange=val; }
  bool    GetIsComputeCSDARange() const { return fIsComputeCSDARange; }

  //
  void    SetLowestElectronTrackingEnergy(double val);
  double  GetLowestElectronTrackingEnergy() const { return fLowestElectronTrackingEnergy; }

  void    SetLinearEnergyLossLimit(double val);
  double  GetLinearEnergyLossLimit() const { return fLinearEnergyLossLimit; }


  double  GetDRoverRange() const { return fDRoverRange; }
  double  GetFinalRange()  const { return fFinalRange;  }
  void    SetStepFunction(double roverrange, double finalrange);


  //
  // Lambda table related variables
  void    SetMinLambdaTableEnergy(double val);
  double  GetMinLambdaTableEnergy() const { return fMinLambdaTableEnergy; }

  void    SetMaxLambdaTableEnergy(double val);
  double  GetMaxLambdaTableEnergy() const { return fMaxLambdaTableEnergy; }

  void    SetNumLambdaTableBins(int val);
  int     GetNumLambdaTableBins() const { return fNumLambdaTableBins; }

  void    SetNumLambdaTableBinsPerDecade(int val);
  int     GetNumLambdaTableBinsPerDecade()   const { return fNumLambdaTableBinsPerDecade; }




  // printouts
  friend std::ostream& operator<<(std::ostream&, PhysicsParameters*);
  friend std::ostream& operator<<(std::ostream&, PhysicsParameters&);




private:
  // secondary production kinetic energy threshold related minimum/maximum values
  // These must be the same in each reagion
  static double   fDefaultGammaCutInLength;
  static double   fDefaultElectronCutInLength;
  static double   fDefaultPositronCutInLength;
  static double   fDefaultGammaCutInEnergy;
  static double   fDefaultElectronCutInEnergy;
  static double   fDefaultPositronCutInEnergy;

  static double   fMinAllowedGammaCutEnergy;
  static double   fMaxAllowedGammaCutEnergy;
  static double   fMinAllowedElectronCutEnergy;
  static double   fMaxAllowedElectronCutEnergy;
  static double   fMinAllowedPositronCutEnergy;
  static double   fMaxAllowedPositronCutEnergy;

  //
  // EM energy loss table related variables
  //
  /** minimum of the kinetic energy grid of the energy loss table; (G4->minKinEnergy) */
  double          fMinLossTableEnergy;
  /** maximum of the kinetic energy grid of the energy loss table; (G4->maxKinEnergy) */
  double          fMaxLossTableEnergy;
  /** number of energy bins in the loss table between  fMinLossTableEnergy and fMaxLossTableEnergy */
  int             fNumLossTableBins;
  /** number of energy bins in the loss table per decade */
  int             fNumLossTableBinsPerDecade;

  bool            fIsComputeCSDARange; //if true  unrestricted total (brem+ioni) range will also be computed

  //
  // Lambda tables related variables
  //
  /** minimum of the kinetic energy grid of the lambda table; (G4->minKinEnergy) */
  double          fMinLambdaTableEnergy;
  /** maximum of the kinetic energy grid of the lambda table; (G4->maxKinEnergy) */
  double          fMaxLambdaTableEnergy;
  /** number of energy bins in the lambda tables between  fMinLambdaTableEnergy and fMaxMaxTableEnergy */
  int             fNumLambdaTableBins;
  /** number of energy bins in the lambda tables per decade */
  int             fNumLambdaTableBinsPerDecade;

  // EM tracking cuts for charged particles: particles with kinetic energy lower than this value will be stopped in
  // the EnergyLoss process AlongStepDoIt.
  double          fLowestElectronTrackingEnergy;

  double          fLinearEnergyLossLimit;  // result of linear energy loss approximation is accepted in the along
                                           // step energy loss computation (ELossTable::GetMeanEnergyAfterAStep) if the
                                           // (linear energy loss)/(initial kinetic energy) <= fLinearEnergyLossLimit
  // continuous step limit function parameters for kEnergyLoss EMPhysicsProcess-es
  double          fFinalRange;
  double          fDRoverRange;

  std::vector<bool> fListActiveRegions;  /** is this PhysicsParameters active in the i-th region?;
                                             will be set by the PhysicsListManager */
  /** each created physics parameters will be registered in this table; the class DO OWN the objects  */
  static std::vector<PhysicsParameters*>  gThePhysicsParametersTable;
};

}  // namespace geantphysics

#endif // PHYSICSPARAMETERS_H
