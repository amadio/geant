#include "GPIoni.h"
//#include "GPRandom.h"
#include "GPVector3.h"
#include <stdio.h>

FQUALIFIER G4double rand_wrapper(curandState* devStates, int id);

FQUALIFIER
void GPIoni::Print() {

	printf(
			"threadId(%d), useLambdaTable(%d), integral(%d), lambdaFactor(%lf), theLambdaTable(%ld), \
mfpKinEnergy(%lf), preStepKinEnergy(%lf), preStepLambda(%lf), fFactor(%lf), \
theNumberOfInteractionLengthLeft(%lf), currentInteractionLength(%lf), isElectron(%d), LPMFlag(%d)\n",
			threadId, int(useLambdaTable), int(integral), lambdaFactor,
			long(theLambdaTable), mfpKinEnergy, preStepKinEnergy, preStepLambda,
			fFactor, theNumberOfInteractionLengthLeft, currentInteractionLength,
			int(isElectron), int(LPMFlag));

}

FQUALIFIER GPIoni::GPIoni() {
	threadId = -1;
	devStates = 0;
	useLambdaTable = false;
	integral = false;

	lambdaFactor = 0.8;
	theLambdaTable = 0;

	mfpKinEnergy = DBL_MAX;
	preStepKinEnergy = 0.0;
	preStepLambda = 0.0;
	fFactor = 1.0;

	theNumberOfInteractionLengthLeft = -1;
	currentInteractionLength = -1;

	isElectron = true;
	LPMFlag = false;
	for (int i = 0; i < nElements; i++)
		xsec[i] = 0.0;
}

FQUALIFIER
void GPIoni::SetCurandStates(curandState* v) {
	devStates = v;
}

FQUALIFIER
void GPIoni::UseIntegral(bool v) {
	integral = v;
}

FQUALIFIER
void GPIoni::InitialiseStep(G4double kineticEnergy) {
	preStepKinEnergy = kineticEnergy;
	if (theNumberOfInteractionLengthLeft < 0.0)
		mfpKinEnergy = DBL_MAX;
}

FQUALIFIER
void GPIoni::SetLambdaTable(GPPhysicsTable* val) {
	theLambdaTable = val;
}

FQUALIFIER
void GPIoni::UseLambdaTable(bool v) {
	useLambdaTable = v;
}

FQUALIFIER
void GPIoni::ResetNumberOfInteractionLengthLeft() {
	theNumberOfInteractionLengthLeft = -log(rand_wrapper(devStates, threadId));
}

FQUALIFIER
void GPIoni::SubtractNumberOfInteractionLengthLeft(
		G4double previousStepSize) {
	if (currentInteractionLength > 0.0) {
		theNumberOfInteractionLengthLeft -= previousStepSize
				/ currentInteractionLength;
		if (theNumberOfInteractionLengthLeft < 0.0) {
			theNumberOfInteractionLengthLeft = perMillion;
		}
	}
}


FQUALIFIER G4double GPIoni::ComputeCrossSectionPerElectron(double kineticEnergy, double cutEnergy, double maxEnergy) {

  G4double cross = 0.0;
  G4double tmax = isElectron ? 0.5*kineticEnergy : kineticEnergy;
  tmax = min(maxEnergy, tmax);

  if(cutEnergy < tmax) {

    G4double xmin  = cutEnergy/kineticEnergy;
    G4double xmax  = tmax/kineticEnergy;
    G4double tau   = kineticEnergy/electron_mass_c2;
    G4double gam   = tau + 1.0;
    G4double gamma2= gam*gam;
    G4double beta2 = tau*(tau + 2)/gamma2;

    //Moller (e-e-) scattering
    if (isElectron) {

      G4double g     = (2.0*gam - 1.0)/gamma2;
      cross = ((xmax - xmin)*(1.0 - g + 1.0/(xmin*xmax)
			      + 1.0/((1.0-xmin)*(1.0 - xmax)))
            - g*log( xmax*(1.0 - xmin)/(xmin*(1.0 - xmax)) ) ) / beta2;

    //Bhabha (e+e-) scattering
    } else {

      G4double y   = 1.0/(1.0 + gam);
      G4double y2  = y*y;
      G4double y12 = 1.0 - 2.0*y;
      G4double b1  = 2.0 - y2;
      G4double b2  = y12*(3.0 + y2);
      G4double y122= y12*y12;
      G4double b4  = y122*y12;
      G4double b3  = b4 + y122;

      cross = (xmax - xmin)*(1.0/(beta2*xmin*xmax) + b2
            - 0.5*b3*(xmin + xmax)
	    + b4*(xmin*xmin + xmin*xmax + xmax*xmax)/3.0)
            - b1*log(xmax/xmin);
    }

    cross *= twopi_mc2_rcl2/kineticEnergy;
  }
  return cross;
}


FQUALIFIER
G4double GPIoni::ComputeCrossSectionPerAtom(double kineticEnergy, double Z,
		G4double cut) {
  return Z*ComputeCrossSectionPerElectron(kineticEnergy,cut, 100.0*TeV);
}

FQUALIFIER
G4double GPIoni::CrossSectionPerVolume(double kineticEnergy,
		G4double cutEnergy, double maxEnergy) {
  G4double electronDensity = 2.06011e+21;
  return electronDensity*ComputeCrossSectionPerElectron(kineticEnergy,cutEnergy,maxEnergy);
}

FQUALIFIER
G4double GPIoni::GetCurrentLambda(double e) {

  G4double x = 0.0;
  if (useLambdaTable && theLambdaTable) {
    //    x = GPIoni::GetLambdaFromTable(e);
    //  return fFactor*((*theLambdaTable)[basedCoupleIndex])->Value(e);
    // basedCoupleIndex 0=gamma, 1=electron, 2=positron
    x = fFactor * theLambdaTable->physicsVectors[1].Value(e);
  } else {
    x = fFactor * CrossSectionPerVolume(e, 0.00099, 100 * TeV);
  }
  return x;
}

FQUALIFIER
void GPIoni::ComputeIntegralLambda(G4double e) {

	// currentCoupleIndex : 1 for PbW)4
	// theEnergyOfCrossSectionMax[0]=1e+07
	// theCrossSectionMax[0]=4.79785e-27
	// theEnergyOfCrossSectionMax[1]=0.00713103
	// theCrossSectionMax[1]=3.4635

	//  mfpKinEnergy  = theEnergyOfCrossSectionMax[currentCoupleIndex];
	mfpKinEnergy = 0.00713103;
	if (e <= mfpKinEnergy) {
		preStepLambda = GetCurrentLambda(e);
	} else {
		G4double e1 = e * lambdaFactor;
		if (e1 > mfpKinEnergy) {
			preStepLambda = GetCurrentLambda(e);
			G4double preStepLambda1 = GetCurrentLambda(e1);
			if (preStepLambda1 > preStepLambda) {
				mfpKinEnergy = e1;
				preStepLambda = preStepLambda1;
			}
		} else {
			//      preStepLambda = fFactor*theCrossSectionMax[currentCoupleIndex];
			preStepLambda = 3.4635;
		}
	}

}

FQUALIFIER
G4double GPIoni::PostStepGetPhysicalInteractionLength(double kineticEnergy,
		G4double previousStepSize, GPForceCondition* condition) {

  // // condition is set to "Not Forced"
  *condition = NotForced;

  G4double x = DBL_MAX;

  // // initialisation of material, mass, charge, model at the beginning of the step
  // const G4ParticleDefinition* currPart = track.GetParticleDefinition();
  // if(isIon) {
  //   if(baseParticle) {
  //     massRatio = baseParticle->GetPDGMass()/currPart->GetPDGMass();
  //   } else {
  //     massRatio = proton_mass_c2/currPart->GetPDGMass();
  //   }
  // }
  // dwjang : massRatio for electron is 1.0
  
  /*
    if(!theDensityFactor || !theDensityIdx) {
    G4cout << "G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength 1: "
    <<  theDensityFactor << "  " << theDensityIdx
    << G4endl;
    G4cout << track.GetDefinition()->GetParticleName() 
    << " e(MeV)= " << track.GetKineticEnergy()
    << " mat " << track.GetMaterialCutsCouple()->GetMaterial()->GetName()
    << G4endl;
    }
  */
  //  DefineMaterial(track.GetMaterialCutsCouple());
  preStepKinEnergy    = kineticEnergy;
  // dwjang : preStepScaledEnergy and preStepKinEnergy are same since massRatio is 1. So just use preStepKinEnergy
  //  preStepScaledEnergy = preStepKinEnergy*massRatio;

  // dwjang : model selection checks whether this model is applicable in the given energy range.
  //  SelectModel(preStepScaledEnergy);
  //  if(!currentModel->IsActive(preStepScaledEnergy)) { return x; }

  // dwjang : for electrons, it doesn't really matter
  // change effective charge of an ion on fly
  // if(isIon) {
  //   G4G4double q2 = currentModel->ChargeSquareRatio(track);
  //   if(q2 != chargeSqRatio) {
  //     chargeSqRatio = q2;
  //     fFactor = q2*biasFactor*(*theDensityFactor)[currentCoupleIndex];
  //     reduceFactor = 1.0/(fFactor*massRatio);
  //   }
  // }
  // q2 = 1, currentCoupleIndex = 1, (*theDensityFactor)[currentCoupleIndex] = 1, fFactor = 1
  // biasFactor = 1, massRatio = 1, reduceFactor = 1, so every value is 1

  // initialisation for sampling of the interaction length 
  if(previousStepSize <= 0.0) { theNumberOfInteractionLengthLeft = -1.0; }
  if(theNumberOfInteractionLengthLeft < 0.0) { mfpKinEnergy = DBL_MAX; }

  // dwjang : not considering bias at this moment
  // // forced biasing only for primary particles
  // if(biasManager) {
  //   if(0 == track.GetParentID()) {
  //     if(0 == track.GetCurrentStepNumber()) {
  //       biasFlag = true; 
  // biasManager->ResetForcedInteraction(); 
  //     }
  //     if(biasFlag && biasManager->ForcedInteractionRegion(currentCoupleIndex)) {
  //       return biasManager->GetStepLimit(currentCoupleIndex, previousStepSize);
  //     }
  //   }
  // }

  // compute mean free path
  //  if(preStepScaledEnergy < mfpKinEnergy) {
  if(preStepKinEnergy < mfpKinEnergy) {
    //    if (integral) { ComputeLambdaForScaledEnergy(preStepScaledEnergy); }
    if (integral) { ComputeIntegralLambda(preStepKinEnergy); }
    //    else  { preStepLambda = GetLambdaForScaledEnergy(preStepScaledEnergy); }
    else  { preStepLambda = GetCurrentLambda(preStepKinEnergy); }
    if(preStepLambda <= 0.0) { mfpKinEnergy = 0.0; }
  }

  // non-zero cross section
  if(preStepLambda > 0.0) { 
    if (theNumberOfInteractionLengthLeft < 0.0) {
      // beggining of tracking (or just after DoIt of this process)
      ResetNumberOfInteractionLengthLeft();
    } else if(currentInteractionLength < DBL_MAX) {
      // subtract NumberOfInteractionLengthLeft
      SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if(theNumberOfInteractionLengthLeft < 0.) {
	theNumberOfInteractionLengthLeft = perMillion;
      }
    }

    // get mean free path and step limit
    currentInteractionLength = 1.0/preStepLambda;
    x = theNumberOfInteractionLengthLeft * currentInteractionLength;

    // zero cross section case
  } else {
    if(theNumberOfInteractionLengthLeft > DBL_MIN && 
       currentInteractionLength < DBL_MAX) {
      // subtract NumberOfInteractionLengthLeft
      SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if(theNumberOfInteractionLengthLeft < 0.) {
	theNumberOfInteractionLengthLeft = perMillion;
      }
    }
    currentInteractionLength = DBL_MAX;
  }
  return x;

}


FQUALIFIER
G4double GPIoni::SampleSecondaries(GPTrack* dp, double tmin,
				 G4double maxEnergy) {

  G4double deltaKinEnergy = 0.0; // to be returned

  G4double kineticEnergy = dp->E;
  //const G4Material* mat = couple->GetMaterial();
  //G4double Zeff = mat->GetElectronDensity()/mat->GetTotNbOfAtomsPerVolume();
  // G4double th   = 0.25*sqrt(Zeff)*keV;
  G4double tmax;
  if(isElectron) { 
    tmax = 0.5*kineticEnergy; 
  } else {
    tmax = kineticEnergy; 
  }
  if(maxEnergy < tmax) { tmax = maxEnergy; }
  if(tmin >= tmax) { return deltaKinEnergy; }

  G4double energy = kineticEnergy + electron_mass_c2;
  G4double totalMomentum = sqrt(kineticEnergy*(energy + electron_mass_c2));
  G4double xmin   = tmin/kineticEnergy;
  G4double xmax   = tmax/kineticEnergy;
  G4double gam    = energy/electron_mass_c2;
  G4double gamma2 = gam*gam;
  G4double beta2  = 1.0 - 1.0/gamma2;
  G4double x, z, q, grej;

  GPVector3 direction(dp->px,dp->py,dp->pz);

  //Moller (e-e-) scattering
  if (isElectron) {

    G4double g = (2.0*gam - 1.0)/gamma2;
    G4double y = 1.0 - xmax;
    grej = 1.0 - g*xmax + xmax*xmax*(1.0 - g + (1.0 - g*y)/(y*y));

    do {
      q = rand_wrapper(devStates,threadId);
      x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      y = 1.0 - x;
      z = 1.0 - g*x + x*x*(1.0 - g + (1.0 - g*y)/(y*y));
    } while(grej * rand_wrapper(devStates,threadId) > z);

  //Bhabha (e+e-) scattering
  } else {

    G4double y   = 1.0/(1.0 + gam);
    G4double y2  = y*y;
    G4double y12 = 1.0 - 2.0*y;
    G4double b1  = 2.0 - y2;
    G4double b2  = y12*(3.0 + y2);
    G4double y122= y12*y12;
    G4double b4  = y122*y12;
    G4double b3  = b4 + y122;

    y    = xmax*xmax;
    grej = 1.0 + (y*y*b4 - xmin*xmin*xmin*b3 + y*b2 - xmin*b1)*beta2; 
    do {
      q = rand_wrapper(devStates,threadId);
      x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      y = x*x;
      z = 1.0 + (y*y*b4 - x*y*b3 + y*b2 - x*b1)*beta2; 
    } while(grej * rand_wrapper(devStates,threadId) > z);
  }

  deltaKinEnergy = x * kineticEnergy;

  // dwjang : wait until deciding how to handle secondaries

  // G4double deltaMomentum = sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  // G4double cost = deltaKinEnergy * (energy + electron_mass_c2) / (deltaMomentum * totalMomentum);
  // G4double sint = (1.0 - cost)*(1. + cost);
  // if(sint > 0.0) { sint = sqrt(sint); }
  // else { sint = 0.0; }

  // G4double phi = twopi * rand_wrapper(devStates,threadId);

  // GPVector3 deltaDirection(sint*cos(phi),sint*sin(phi), cost) ;
  // deltaDirection.RotateUz(direction);

  // // primary change
  // kineticEnergy -= deltaKinEnergy;
  // //  fParticleChange->SetProposedKineticEnergy(kineticEnergy);

  // GPVector3 dir = totalMomentum*direction - deltaMomentum*deltaDirection;
  // direction = dir.Unit();
  // //  fParticleChange->SetProposedMomentumDirection(direction);

  // // create G4DynamicParticle object for delta ray
  // // G4DynamicParticle* delta = new G4DynamicParticle(theElectron,
  // // 						   deltaDirection,deltaKinEnergy);
  // // vdp->push_back(delta);

  return deltaKinEnergy;
}



FQUALIFIER
//G4double GPVEnergyLossProcess_PostStepDoIt(GPTrack* dp)
G4double GPIoni::PostStepDoIt(GPTrack* dp) {

  G4double eloss = 0.0; // to be returned

  // In all cases clear number of interaction lengths
  theNumberOfInteractionLengthLeft = -1.0;

  //  fParticleChange.InitializeForPostStep(track);
  G4double finalT = dp->E;
  G4double lowestKinEnergy = 1.0*MeV;
  if(finalT <= lowestKinEnergy) { return eloss; }

  //  G4double massRatio = 1.0;

  //  G4double postStepScaledEnergy = finalT*massRatio;
  G4double postStepScaledEnergy = finalT;

  //  if(!currentModel->IsActive(postStepScaledEnergy)) { return &fParticleChange; }

  // // forced process - should happen only once per track
  // if(biasFlag) {
  //   if(biasManager->ForcedInteractionRegion(currentCoupleIndex)) {
  //     biasFlag = false;
  //   }
  // }

  // Integral approach
  if (integral) {
    //    G4double lx = GetLambdaForScaledEnergy(postStepScaledEnergy);
    G4double lx = GetCurrentLambda(postStepScaledEnergy);
    if(lx <= 0.0) {
      return eloss;
    } else if(preStepLambda*rand_wrapper(devStates, threadId) > lx) {
      return eloss;
    }
  }

  // SelectModel(postStepScaledEnergy);

  // // define new weight for primary and secondaries
  // G4double weight = fParticleChange.GetParentWeight();
  // if(weightFlag) {
  //   weight /= biasFactor;
  //   fParticleChange.ProposeParentWeight(weight);
  // }

  // const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  // G4G4double tcut = (*theCuts)[currentCoupleIndex];
  // dwjang : weight = 1, tcut = 0.00099 MeV
  G4double tcut = 0.00099;

  // // sample secondaries
  // secParticles.clear();
  // currentModel->SampleSecondaries(&secParticles, currentCouple, dynParticle, tcut);

  // // bremsstrahlung splitting or Russian roulette  
  // if(biasManager) {
  //   if(biasManager->SecondaryBiasingRegion(currentCoupleIndex)) {
  //     weight *= biasManager->ApplySecondaryBiasing(secParticles,currentCoupleIndex);
  //   }
  // }

  eloss = SampleSecondaries(dp, tcut, 100 * TeV);

  // // save secondaries
  // G4int num = secParticles.size();
  // if(num > 0) {
  //   fParticleChange.SetNumberOfSecondaries(num);
  //   for (G4int i=0; i<num; ++i) {

  //     if(secParticles[i]) {
  // 	G4Track* t = new G4Track(secParticles[i], track.GetGlobalTime(), 
  // 				 track.GetPosition());
  // 	t->SetTouchableHandle(track.GetTouchableHandle());
  // 	t->SetWeight(weight); 
  // 	pParticleChange->AddSecondary(t);
  //     }
  //   }
  // }

  // if(0.0 == fParticleChange.GetProposedKineticEnergy() &&
  //    fAlive == fParticleChange.GetTrackStatus()) {
  //   if(particle->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
  //     { fParticleChange.ProposeTrackStatus(fStopButAlive); }
  //   else { fParticleChange.ProposeTrackStatus(fStopAndKill); }
  // }

  // //
  // //   if(-1 < verboseLevel) {
  // //   G4cout << "::PostStepDoIt: Sample secondary; Efin= " 
  // //   << fParticleChange.GetProposedKineticEnergy()/MeV
  // //   << " MeV; model= (" << currentModel->LowEnergyLimit()
  // //   << ", " <<  currentModel->HighEnergyLimit() << ")"
  // //   << "  preStepLambda= " << preStepLambda
  // //   << "  dir= " << track.GetMomentumDirection()
  // //   << "  status= " << track.GetTrackStatus()
  // //   << G4endl;
  // //   }
  // //
  // return &fParticleChange;

  return eloss;
}
