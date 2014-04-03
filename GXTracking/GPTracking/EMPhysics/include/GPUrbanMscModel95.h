#ifndef GPUrbanMscModel95_h
#define GPUrbanMscModel95_h 1

//#include "G4VMscModel.hh"
//#include "G4MscStepLimitType.hh"

//class G4ParticleChangeForMSC;
//class G4SafetyHelper;
//class G4LossTableManager;
#include "GPVParticleChange.h"

#include "GPThreeVector.h"
#include "GXTrack.h"
#include "GPMaterial.h"
#include "GPPhysicsTable.h"

#include "GPRandom.h"

class GPUrbanMscModel95 //: public G4VMscModel
{

public:

  FQUALIFIER GPUrbanMscModel95(curandState* devStates, 
			       int threadId);

  FQUALIFIER ~GPUrbanMscModel95();

  //  void Initialise(const G4ParticleDefinition*, const G4DataVector&);
  FQUALIFIER void Initialise();

  FQUALIFIER void StartTracking();

  /* 
  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition* particle,
				      G4double KineticEnergy,
				      G4double AtomicNumber,
				      G4double AtomicWeight=0., 
				      G4double cut =0.,
				      G4double emax=DBL_MAX);
  */

  FQUALIFIER
  GPThreeVector SampleScattering(GPMaterial* material,
				 GPThreeVector direction, 
				 G4double safety);

  FQUALIFIER
  G4double ComputeTruePathLengthLimit(GPMaterial* material,
				      G4double kineticEnergy,
				      G4double xlambda,
				      G4double& currentMinimalStep);

  FQUALIFIER
  G4double ConvertTrueToGeom(GPMaterial* material,
			     G4double& tlength,
			     G4double& glength);

  FQUALIFIER 
  G4double ComputeGeomPathLength(GPMaterial* material);

  FQUALIFIER
  G4double ComputeTruePathLengthLimit2(GPMaterial* material,
				      G4double kineticEnergy,
				      G4double xlambda,
				      G4double* currentMinimalStep);

  FQUALIFIER
  G4double ConvertTrueToGeom2(GPMaterial* material,
			      G4double* tlength,
			      G4double* glength);

  FQUALIFIER 
    G4double ComputeGeomPathLength2(GPMaterial* material);

  FQUALIFIER 
  G4double ComputeTrueStepLength(G4double geomStepLength);

  FQUALIFIER 
  G4double ComputeTheta0(G4double truePathLength,
                         G4double KineticEnergy);


  //from base G4VMscModel

  FQUALIFIER 
  void SetSampleZ(G4bool val);

  FQUALIFIER
  G4double GetDEDX(G4double kinEnergy);

  FQUALIFIER
  G4double GetRange(GPMaterial* material,
		    G4double kinEnergy);

  FQUALIFIER
  G4double GetEnergy(GPMaterial* material,
		     G4double range);

  FQUALIFIER
  G4double GetTransportMeanFreePath(G4double ekin);

  FQUALIFIER
  void SetLambdaTable(GPPhysicsTable* lambdaTable);

  //G4VEmModel
  FQUALIFIER void SetParticleChange(GPVParticleChange* p);
  FQUALIFIER G4double HighEnergyLimit();
  FQUALIFIER G4double LowEnergyLimit();
  FQUALIFIER G4double HighEnergyActivationLimit();
  FQUALIFIER G4double LowEnergyActivationLimit();
  FQUALIFIER void SetHighEnergyLimit(G4double val);
  FQUALIFIER void SetLowEnergyLimit(G4double val);
  FQUALIFIER void SetActivationHighEnergyLimit(G4double val);
  FQUALIFIER void SetActivationLowEnergyLimit(G4double val);
  FQUALIFIER G4bool IsActive(G4double kinEnergy);

private:

  FQUALIFIER 
  G4double SimpleScattering(G4double xmeanth, G4double x2meanth);

  FQUALIFIER 
  G4double SampleCosineTheta(GPMaterial* material,
			     G4double trueStepLength,
			     G4double KineticEnergy);
  //			     G4double xlambda);

  FQUALIFIER 
  G4double SampleDisplacement();

  FQUALIFIER 
  G4double LatCorrelation();

  FQUALIFIER 
  void SetParticle();

  FQUALIFIER 
  void UpdateCache();


  //  hide assignment operator
  //  GPUrbanMscModel95 & operator=(const  GPUrbanMscModel95 &right);
  //  GPUrbanMscModel95(const  GPUrbanMscModel95&);

  //  const G4ParticleDefinition* particle;
  //  G4ParticleChangeForMSC*     fParticleChange;

  //  const G4MaterialCutsCouple* couple;
  //  G4LossTableManager*         theManager;

  //curand                  
  int fThreadId;
  curandState* fDevStates;  
  GPPhysicsTable* theLambdaTable;

  //from base G4VMscModel
  G4double skin;
  G4double dtrl;
  G4double facrange;
  G4double facsafety;
  G4bool latDisplasment;
  G4bool samplez;
  G4double dedx;
  G4double localtkin;
  G4double localrange;

  G4double mass;
  G4double charge;
  G4double masslimite,lambdalimit,fr;

  G4double taubig;
  G4double tausmall;
  G4double taulim;
  G4double currentTau;
  G4double tlimit;
  G4double tlimitmin;
  G4double tlimitminfix;
  G4double tgeom;

  G4double geombig;
  G4double geommin;
  G4double geomlimit;
  G4double skindepth;
  G4double smallstep;

  G4double presafety;

  G4double lambda0;
  G4double lambdaeff;
  G4double tPathLength;
  G4double zPathLength;
  G4double par1,par2,par3;

  G4double stepmin;

  G4double currentKinEnergy;
  G4double currentRange; 
  G4double rangeinit;
  G4double currentRadLength;

  G4int    currentMaterialIndex;

  G4double Zold;
  G4double Zeff,Z2,Z23,lnZ;
  G4double coeffth1,coeffth2;
  G4double coeffc1,coeffc2,coeffc3,coeffc4;
  G4double scr1ini,scr2ini,scr1,scr2;

  G4bool   firstStep;
  G4bool   inside;
  G4bool   insideskin;

  //G4VEmModel
  GPVParticleChange* pParticleChange;
  G4double        lowLimit;
  G4double        highLimit;
  G4double        eMinActive;
  G4double        eMaxActive;

};

#endif

