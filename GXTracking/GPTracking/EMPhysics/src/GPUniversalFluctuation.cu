#include "GPUniversalFluctuation.h"
#include "GPPhysicalConstants.h"
#include "GPSystemOfUnits.h"

//#include "Randomize.hh"
//#include "G4Poisson.hh"
//#include "G4Step.hh"
//#include "GPMaterial.h"
#include "GPIonisParamMat.h"
//#include "G4DynamicParticle.hh"
//#include "G4ParticleDefinition.hh"

#include "stdio.h"

FQUALIFIER 
GPUniversalFluctuation::GPUniversalFluctuation(curandState* devStates, 
					       int threadId)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

  lastMaterial = 0;
  minNumberInteractionsBohr = 10.0;
  theBohrBeta2 = 50.0*keV/proton_mass_c2;
  minLoss = 10.*eV;
  nmaxCont = 16.;
  rate = 0.55;
  fw = 4.;

  particleMass = chargeSquare = ipotFluct = electronDensity = f1Fluct = f2Fluct 
    = e1Fluct = e2Fluct = e1LogFluct = e2LogFluct = ipotLogFluct = e0 = esmall 
    = e1 = e2 = 0;

}

FQUALIFIER GPUniversalFluctuation::~GPUniversalFluctuation()
{}

//void GPUniversalFluctuation::InitialiseMe(const G4ParticleDefinition* part)
FQUALIFIER
void GPUniversalFluctuation::InitialiseMe()
{
  //  particle       = part;
  //  particleMass   = part->GetPDGMass();
  particleMass   = electron_mass_c2;
  //  G4double q     = part->GetPDGCharge()/eplus;
  G4double q     = electron_charge/eplus;
  chargeSquare   = q*q;
}

FQUALIFIER
G4double GPUniversalFluctuation::SampleFluctuations(GPMaterial* material,
						    G4double kineticEnergy,
						    G4double& tmax,
						    G4double& length,
						    G4double& averageLoss)
{
  //  printf("Inside SampleFluctuations\n");

  // Calculate actual loss from the mean loss.
  // The model used to get the fluctuations is essentially the same
  // as in Glandz in Geant3 (Cern program library W5013, phys332).
  // L. Urban et al. NIM A362, p.416 (1995) and Geant4 Physics Reference Manual

  // shortcut for very small loss or from a step nearly equal to the range
  // (out of validity of the model)
  //
  G4double meanLoss = averageLoss;
  // G4double tkin  = dp->GetKineticEnergy();
  G4double tkin  = kineticEnergy; 
  if (meanLoss < minLoss) { return meanLoss; }

  // if(!particle) { InitialiseMe(dp->GetDefinition()); }
  InitialiseMe();
  
  G4double tau   = tkin/particleMass;
  G4double gam   = tau + 1.0;
  G4double gam2  = gam*gam;
  G4double beta2 = tau*(tau + 2.0)/gam2;

  G4double loss(0.), siga(0.);
  
  // Gaussian regime
  // for heavy particles only and conditions
  // for Gauusian fluct. has been changed 
  //
  if ((particleMass > electron_mass_c2) &&
      (meanLoss >= minNumberInteractionsBohr*tmax))
  {
    G4double massrate = electron_mass_c2/particleMass ;
    G4double tmaxkine = 2.*electron_mass_c2*beta2*gam2/
                        (1.+massrate*(2.*gam+massrate)) ;
    if (tmaxkine <= 2.*tmax)   
    {
      electronDensity = GPMaterial_GetElectronDensity(material);
      siga = sqrt((1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
		  * electronDensity * chargeSquare);

     
      G4double sn = meanLoss/siga;
  
      // thick target case 
      if (sn >= 2.0) {

	G4double twomeanLoss = meanLoss + meanLoss;
	do {
	  // loss = G4RandGauss::shoot(meanLoss,siga);
	  loss = GPRandGauss(fDevStates,fThreadId,meanLoss,siga);
	} while (0.0 > loss || twomeanLoss < loss);

	// Gamma distribution
      } else {

	G4double neff = sn*sn;
	// loss = meanLoss*CLHEP::RandGamma::shoot(neff,1.0)/neff;
	loss = meanLoss*GPRandGamma(fDevStates,fThreadId,neff,1.0)/neff;
      }
      return loss;
    }
  }

  // Glandz regime : initialisation
  //
  if (material != lastMaterial) {

    GPIonisParamMat aIonisParamMat;
    GPIonisParamMat_Constructor(&aIonisParamMat,material);
    /*
    f1Fluct      = material->GetIonisation()->GetF1fluct();
    f2Fluct      = material->GetIonisation()->GetF2fluct();
    e1Fluct      = material->GetIonisation()->GetEnergy1fluct();
    e2Fluct      = material->GetIonisation()->GetEnergy2fluct();
    e1LogFluct   = material->GetIonisation()->GetLogEnergy1fluct();
    e2LogFluct   = material->GetIonisation()->GetLogEnergy2fluct();
    ipotFluct    = material->GetIonisation()->GetMeanExcitationEnergy();
    ipotLogFluct = material->GetIonisation()->GetLogMeanExcEnergy();
    e0           = material->GetIonisation()->GetEnergy0fluct();
    */
    f1Fluct      = GPIonisParamMat_GetF1fluct(&aIonisParamMat);
    f2Fluct      = GPIonisParamMat_GetF2fluct(&aIonisParamMat);
    e1Fluct      = GPIonisParamMat_GetEnergy1fluct(&aIonisParamMat);
    e2Fluct      = GPIonisParamMat_GetEnergy2fluct(&aIonisParamMat);
    e1LogFluct   = GPIonisParamMat_GetLogEnergy1fluct(&aIonisParamMat);
    e2LogFluct   = GPIonisParamMat_GetLogEnergy2fluct(&aIonisParamMat);
    ipotFluct    = GPIonisParamMat_GetMeanExcitationEnergy(&aIonisParamMat);
    ipotLogFluct = GPIonisParamMat_GetLogMeanExcEnergy(&aIonisParamMat);
    e0           = GPIonisParamMat_GetEnergy0fluct(&aIonisParamMat);

    esmall = 0.5*sqrt(e0*ipotFluct);  
    lastMaterial = material;   
  }

  // very small step or low-density material
  if(tmax <= e0) { return meanLoss; }

  G4double losstot = 0.;
  G4int    nstep = 1;
  if(meanLoss < 25.*ipotFluct)
    {
      // if(G4UniformRand() < 0.04*meanLoss/ipotFluct)
      if(GPUniformRand(fDevStates,fThreadId) < 0.04*meanLoss/ipotFluct)
	{ nstep = 1; }
      else
	{ 
	  nstep = 2;
	  meanLoss *= 0.5; 
	}
    }

  for (G4int istep=0; istep < nstep; ++istep) {
    
    loss = 0.;

    G4double a1 = 0. , a2 = 0., a3 = 0. ;

    if(tmax > ipotFluct) {
      G4double w2 = log(2.*electron_mass_c2*beta2*gam2)-beta2;

      if(w2 > ipotLogFluct)  {
	G4double C = meanLoss*(1.-rate)/(w2-ipotLogFluct);
	a1 = C*f1Fluct*(w2-e1LogFluct)/e1Fluct;
	if(w2 > e2LogFluct) {
	  a2 = C*f2Fluct*(w2-e2LogFluct)/e2Fluct;
	}
	if(a1 < nmaxCont) { 
	  //small energy loss
	  G4double sa1 = sqrt(a1);
	  // if(G4UniformRand() < exp(-sa1))
	  if( GPUniformRand(fDevStates,fThreadId) < exp(-sa1))
	    {
	      e1 = esmall;
	      a1 = meanLoss*(1.-rate)/e1;
	      a2 = 0.;
	      e2 = e2Fluct;
	    } 
	  else
	    {
	      a1 = sa1 ;    
	      e1 = sa1*e1Fluct;
	      e2 = e2Fluct;
	    }

	} else {
	    //not small energy loss
	    //correction to get better fwhm value
	    a1 /= fw;
	    e1 = fw*e1Fluct;
	    e2 = e2Fluct;
	}
      }   
    }

    G4double w1 = tmax/e0;
    if(tmax > e0) {
      a3 = rate*meanLoss*(tmax-e0)/(e0*tmax*log(w1));
    }
    //'nearly' Gaussian fluctuation if a1>nmaxCont&&a2>nmaxCont&&a3>nmaxCont  
    G4double emean = 0.;
    G4double sig2e = 0., sige = 0.;
    G4double p1 = 0., p2 = 0., p3 = 0.;
 
    // excitation of type 1
    if(a1 > nmaxCont)
      {
	emean += a1*e1;
	sig2e += a1*e1*e1;
      }
    else if(a1 > 0.)
      {
	// p1 = G4double(G4Poisson(a1));
	p1 = G4double(GPPoisson(fDevStates,fThreadId,a1));
	loss += p1*e1;
	if(p1 > 0.) {
	  // loss += (1.-2.*G4UniformRand())*e1;
	  loss += (1.-2.*GPUniformRand(fDevStates,fThreadId))*e1;
	}
      }

    // excitation of type 2
    if(a2 > nmaxCont)
      {
	emean += a2*e2;
	sig2e += a2*e2*e2;
      }
    else if(a2 > 0.)
      {
	// p2 = G4double(G4Poisson(a2));
	p2 = G4double(GPPoisson(fDevStates,fThreadId,a2));
	loss += p2*e2;
	if(p2 > 0.) 
	  // loss += (1.-2.*G4UniformRand())*e2;
	  loss += (1.-2.*GPUniformRand(fDevStates,fThreadId))*e2;
      }

    if(emean > 0.)
      {
	sige   = sqrt(sig2e);
	// loss += max(0.,G4RandGauss::shoot(emean,sige));
	loss += fmax(0.,GPRandGauss(fDevStates,fThreadId,emean,sige));
      }

    // ionisation 
    G4double lossc = 0.;
    if(a3 > 0.) {
      emean = 0.;
      sig2e = 0.;
      sige = 0.;
      p3 = a3;
      G4double alfa = 1.;
      if(a3 > nmaxCont)
	{
	  alfa            = w1*(nmaxCont+a3)/(w1*nmaxCont+a3);
	  G4double alfa1  = alfa*log(alfa)/(alfa-1.);
	  G4double namean = a3*w1*(alfa-1.)/((w1-1.)*alfa);
	  emean          += namean*e0*alfa1;
	  sig2e          += e0*e0*namean*(alfa-alfa1*alfa1);
	  p3              = a3-namean;
	}

      G4double w2 = alfa*e0;
      G4double w  = (tmax-w2)/tmax;
      //@@@ G4int nb = G4Poisson(p3);
      G4int nb = GPPoisson(fDevStates,fThreadId,p3);
      if(nb > 0) {
	// for (G4int k=0; k<nb; k++) lossc += w2/(1.-w*G4UniformRand());
	for (G4int k=0; k<nb; k++) {
	  lossc += w2/(1.-w*GPUniformRand(fDevStates,fThreadId));
	}
      }
    }

    if(emean > 0.)
      {
	sige   = sqrt(sig2e);
	// lossc += max(0.,G4RandGauss::shoot(emean,sige));
	lossc += fmax(0.,GPRandGauss(fDevStates,fThreadId,emean,sige));
      }

    loss += lossc;

    losstot += loss;
  }
              
  return losstot;

}

FQUALIFIER
G4double GPUniversalFluctuation::Dispersion(GPMaterial* material,
					    G4double kineticEnergy,
					    G4double& tmax,
					    G4double& length)
{
  //  if(!particle) { InitialiseMe(dp->GetDefinition()); }
  InitialiseMe(); 

  electronDensity = GPMaterial_GetElectronDensity(material);

  G4double gam   = (kineticEnergy)/particleMass + 1.0;
  G4double beta2 = 1.0 - 1.0/(gam*gam);

  G4double siga  = (1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                 * electronDensity * chargeSquare;

  return siga;
}

FQUALIFIER
void GPUniversalFluctuation::SetParticleAndCharge(G4double q2)
{
  //  if(part != particle) {
  //    particle       = part;
  //  particleMass   = part->GetPDGMass();
  particleMass   = electron_mass_c2;
  //@@@may extend for photon
  //  }
  chargeSquare = q2;
}
