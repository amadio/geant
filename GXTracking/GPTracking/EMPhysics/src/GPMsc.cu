#include "GPMsc.h"
//#include "GPRandom.h"
#include <stdio.h>

FQUALIFIER G4double rand_wrapper(curandState* devStates, int id);
FQUALIFIER G4double rand_normal_wrapper(curandState* devStates, int id, double mean, double sigma);

FQUALIFIER
void GPMsc::Print() {

	printf(
			"threadId(%d), useLambdaTable(%d), integral(%d), lambdaFactor(%lf), theLambdaTable(%ld), \
mfpKinEnergy(%lf), preStepKinEnergy(%lf), preStepLambda(%lf), fFactor(%lf), \
theNumberOfInteractionLengthLeft(%lf), currentInteractionLength(%lf), isElectron(%d), LPMFlag(%d)\n",
			threadId, int(useLambdaTable), int(integral), lambdaFactor,
			long(theLambdaTable), mfpKinEnergy, preStepKinEnergy, preStepLambda,
			fFactor, theNumberOfInteractionLengthLeft, currentInteractionLength,
			int(isElectron), int(LPMFlag));

}

FQUALIFIER GPMsc::GPMsc() :
  threadId(-1),
  devStates(0),
  useLambdaTable(true),
  integral(true),
  samplez(true),
  lambdaFactor(0.8),
  theLambdaTable(0),
  mfpKinEnergy(DBL_MAX),
  preStepKinEnergy(0.0),
  preStepLambda(0.0),
  fFactor(1.0),
  theNumberOfInteractionLengthLeft(-1),
  currentInteractionLength(-1),
  isElectron(true),
  LPMFlag(false),
  steppingAlgorithm(fUseSafety),
  skin(1.0),
  facrange(0.04),
  facgeom(2.5),
  facsafety(0.3),
  dtrl(0.05),
  inside(false),
  insideskin(false),
  currentKinEnergy(0.0),
  tPathLength(0.0),
  currentRange(0.0),
  lambda0(0.0),
  presafety(0.0),
  smallstep(1.0e10),
  tlimitminfix(1.0e-6*mm),
  rangeinit(0.0),
  tlimit(1.0e10*mm),
  tlimitmin(1.0e-5*mm), // 10*tlimitminfix
  stepmin(1.0e-6*mm), // tlimitminfix
  skindepth(1.0e-6*mm), // skin*stepmin
  tgeom(1.0e50*mm),
  geombig(1.0e50*mm),
  geommin(1.0e-3*mm),
  geomlimit(1.0e50*mm),
  fr(0.02),
  masslimite(0.6*MeV),
  mass(proton_mass_c2),
  lambdalimit(1.0*mm),
  par1(0.),par2(0.),par3(0.),
  zPathLength(0.),
  taulim(1.0e-6),
  tausmall(1.0e-16),
  lambdaeff(0.),
  third(1./3.)
{
	for (int i = 0; i < nElements; i++) xsec[i] = 0.0;
}

FQUALIFIER
void GPMsc::SetCurandStates(curandState* v) {
	devStates = v;
}

FQUALIFIER
void GPMsc::UseIntegral(bool v) {
	integral = v;
}

FQUALIFIER
void GPMsc::InitialiseStep(G4double kineticEnergy) {
	preStepKinEnergy = kineticEnergy;
	if (theNumberOfInteractionLengthLeft < 0.0)
		mfpKinEnergy = DBL_MAX;
}

FQUALIFIER
void GPMsc::SetLambdaTable(GPPhysicsTable* val) {
	theLambdaTable = val;
}

FQUALIFIER
void GPMsc::UseLambdaTable(bool v) {
	useLambdaTable = v;
}

FQUALIFIER
void GPMsc::ResetNumberOfInteractionLengthLeft() {
	theNumberOfInteractionLengthLeft = -log(rand_wrapper(devStates, threadId));
}

FQUALIFIER
void GPMsc::SubtractNumberOfInteractionLengthLeft(
		G4double previousStepSize) {
	if (currentInteractionLength > 0.0) {
		theNumberOfInteractionLengthLeft -= previousStepSize
				/ currentInteractionLength;
		if (theNumberOfInteractionLengthLeft < 0.0) {
			theNumberOfInteractionLengthLeft = perMillion;
		}
	}
}

FQUALIFIER
G4double GPMsc::SupressionFunction(double kineticEnergy, double gammaEnergy,
		bool LPMFlag) {
	// supression due to the LPM effect+polarisation of the medium/
	// supression due to the polarisation alone

	G4double totEnergy = kineticEnergy + electron_mass_c2;
	G4double totEnergySquare = totEnergy * totEnergy;

	//  G4double LPMEnergy = LPMconstant*(material->GetRadlen()) ;
	G4double LPMconstant = fine_structure_const * electron_mass_c2
			* electron_mass_c2 / (4. * pi * hbarc);
	G4double LPMEnergy = LPMconstant * 8.92421;
	G4double MigdalConstant = classic_electr_radius * electron_Compton_length
			* electron_Compton_length * 4.0 * pi;

	G4double gammaEnergySquare = gammaEnergy * gammaEnergy;

	//  G4double electronDensity = material->GetElectronDensity();
	G4double electronDensity = 2.06011e+21;

	G4double sp = gammaEnergySquare
			/ (gammaEnergySquare
					+ MigdalConstant * totEnergySquare * electronDensity);
	G4double supr = 1.0;

	if (LPMFlag) {
		G4double s2lpm = LPMEnergy * gammaEnergy / totEnergySquare;

		if (s2lpm < 1.) {

			G4double LPMgEnergyLimit = totEnergySquare / LPMEnergy;
			G4double LPMgEnergyLimit2 = LPMgEnergyLimit * LPMgEnergyLimit;
			G4double splim =
					LPMgEnergyLimit2
							/ (LPMgEnergyLimit2
									+ MigdalConstant * totEnergySquare
											* electronDensity);
			G4double w = 1. + 1. / splim;

			if ((1. - sp) < 1.e-6)
				w = s2lpm * (3. - sp);
			else
				w = s2lpm * (1. + 1. / sp);

			supr = (sqrt(w * w + 4. * s2lpm) - w) / (sqrt(w * w + 4.) - w);
			supr /= sp;
		}

	}
	return supr;
}

FQUALIFIER
G4double GPMsc::PositronCorrFactorSigma(double Z, double kineticEnergy,
		G4double cut) {
	//calculates the correction factor for the energy loss due to bremsstrahlung for positrons
	//the same correction is in the (discrete) bremsstrahlung

	const G4double K = 132.9416 * eV;
	const G4double a1 = 4.15e-1, a3 = 2.10e-3, a5 = 54.0e-5;

	G4double x = log(kineticEnergy / (K * Z * Z)), x2 = x * x, x3 = x2 * x;
	G4double eta = 0.5 + atan(a1 * x + a3 * x3 + a5 * x3 * x2) / pi;
	G4double e0 = cut / kineticEnergy;

	G4double factor = 0.0;
	if (e0 < 1.0) {
		factor = log(1. - e0) / eta;
		factor = exp(factor);
	}
	factor = eta * (1. - factor) / e0;

	return factor;
}

FQUALIFIER
G4double GPMsc::ComputeCrossSectionPerAtom(double kineticEnergy, double Z,
		G4double cut) {
	G4double cross = 0.0;

	if (kineticEnergy < keV || kineticEnergy < cut) {
		return cross;
	}

	const G4double ksi = 2.0;
	//  const G4double alfa=1.00;
	const G4double csigh = 0.127, csiglow = 0.25, asiglow = 0.020 * MeV;
	const G4double Tlim = 10. * MeV;

	const G4double xlim = 1.2;

	const int NZ = 8;
	const int Nsig = 11;
	const float ZZ[NZ] = { 2., 4., 6., 14., 26., 50., 82., 92. };
	const float coefsig[NZ][Nsig] = { { 0.4638, 0.37748, 0.32249, -0.060362,
			-0.065004, -0.033457, -0.004583, 0.011954, 0.0030404, -0.0010077,
			-0.00028131 },

	{ 0.50008, 0.33483, 0.34364, -0.086262, -0.055361, -0.028168, -0.0056172,
			0.011129, 0.0027528, -0.00092265, -0.00024348 },

	{ 0.51587, 0.31095, 0.34996, -0.11623, -0.056167, -0.0087154, 0.00053943,
			0.0054092, 0.00077685, -0.00039635, -6.7818e-05 },

	{ 0.55058, 0.25629, 0.35854, -0.080656, -0.054308, -0.049933, -0.00064246,
			0.016597, 0.0021789, -0.001327, -0.00025983 },

	{ 0.5791, 0.26152, 0.38953, -0.17104, -0.099172, 0.024596, 0.023718,
			-0.0039205, -0.0036658, 0.00041749, 0.00023408 },

	{ 0.62085, 0.27045, 0.39073, -0.37916, -0.18878, 0.23905, 0.095028,
			-0.068744, -0.023809, 0.0062408, 0.0020407 },

	{ 0.66053, 0.24513, 0.35404, -0.47275, -0.22837, 0.35647, 0.13203, -0.1049,
			-0.034851, 0.0095046, 0.0030535 },

	{ 0.67143, 0.23079, 0.32256, -0.46248, -0.20013, 0.3506, 0.11779, -0.1024,
			-0.032013, 0.0092279, 0.0028592 } };

	int iz = 0;
	G4double delz = 1.0e6;
	for (int ii = 0; ii < NZ; ii++) {
		G4double absdelz = fabs(Z - ZZ[ii]);
		if (absdelz < delz) {
			iz = ii;
			delz = absdelz;
		}
	}

	G4double xx = log10(kineticEnergy / MeV);
	G4double fs = 1.0;

	if (xx <= xlim) {
		fs = coefsig[iz][Nsig - 1];
		for (int j = Nsig - 2; j >= 0; j--) {
			fs = fs * xx + coefsig[iz][j];
		}
		if (fs < 0.0) {
			fs = 0.0;
		}
	}

	//  cross = Z*(Z+ksi)*(1.-csigh*exp(log(Z)/4.))*pow(log(kineticEnergy/cut),alfa);
	cross = Z * (Z + ksi) * (1. - csigh * exp(log(Z) / 4.))
			* log(kineticEnergy / cut);

	if (kineticEnergy <= Tlim) {
		cross *= exp(csiglow * log(Tlim / kineticEnergy))
				* (1. + asiglow / (sqrt(Z) * kineticEnergy));
	}

	if (!isElectron)
		cross *= PositronCorrFactorSigma(Z, kineticEnergy, cut);

	cross *= fs / Avogadro;

	if (cross < 0.0) {
		cross = 0.0;
	}

	return cross;
}

FQUALIFIER
G4double GPMsc::CrossSectionPerVolume(double kineticEnergy,
		G4double cutEnergy, double maxEnergy) {

	G4double cross = 0.0;
	G4double tmax = (maxEnergy < kineticEnergy) ? maxEnergy : kineticEnergy;
	G4double cut = (cutEnergy < kineticEnergy) ? cutEnergy : kineticEnergy;
	if (cut >= tmax) {
		return cross;
	}

	// material : PbWO4
	// density : 5.16797e+19
	// electron density : 2.06011e+21
	// radlen : 8.92421
	// effZ, effA : 68.3599, 170.87
	// name: Lead, symbol:Pb, Z: 82, density: 1.0958e+19
	// name: Tungstenm, symbol:W, Z: 74, density: 1.0958e+19
	// name: Oxygen, symbol:O2, Z: 8, density: 4.3832e+19

	//  const int nElements = 3;
	const G4double theAtomNumDensity[nElements] = { 1.0958e+19, 1.0958e+19,
			4.3832e+19 };
	const G4double Z[nElements] = { 82, 74, 8 };

	for (int i = 0; i < nElements; i++) {
		cross += theAtomNumDensity[i]
				* ComputeCrossSectionPerAtom(kineticEnergy, Z[i], cut);
		if (tmax < kineticEnergy) {
			cross -= theAtomNumDensity[i]
					* ComputeCrossSectionPerAtom(kineticEnergy, Z[i], tmax);
		}
	}

	// now compute the correction due to the supression(s)

	G4double kmax = tmax;
	G4double kmin = cut;

	G4double totalEnergy = kineticEnergy + electron_mass_c2;
	//  G4double kp2 = MigdalConstant*totalEnergy*totalEnergy
	//                                             *(material->GetElectronDensity());
	G4double MigdalConstant = classic_electr_radius * electron_Compton_length
			* electron_Compton_length * 4.0 * pi;
	// PbW04 electorn density : 2.06011e+21
	G4double kp2 = MigdalConstant * totalEnergy * totalEnergy * 2.06011e+21;

	G4double fsig = 0.;
	int nmax = 100;
	G4double vmin = log(kmin);
	G4double vmax = log(kmax);
	G4double highKinEnergy = 100.0 * TeV;
	int nn = (int) (nmax * (vmax - vmin) / (log(highKinEnergy) - vmin));
	G4double u, fac, c, v, dv, y;
	if (nn > 0) {

		dv = (vmax - vmin) / nn;
		v = vmin - dv;
		for (int n = 0; n <= nn; n++) {

			v += dv;
			u = exp(v);
			fac = SupressionFunction(kineticEnergy, u, LPMFlag);
			y = u / kmax;
			fac *= (4. - 4. * y + 3. * y * y) / 3.;
			//        fac *= probsup*(u*u/(u*u+kp2))+1.-probsup;
			fac *= (u * u / (u * u + kp2));

			if ((n == 0) || (n == nn))
				c = 0.5;
			else
				c = 1.;

			fac *= c;
			fsig += fac;
		}
		y = kmin / kmax;
		fsig *= dv
				/ (-4. * log(y) / 3. - 4. * (1. - y) / 3. + 0.5 * (1. - y * y));

	} else {

		fsig = 1.;
	}
	if (fsig > 1.)
		fsig = 1.;

	// correct the cross section
	cross *= fsig;

	return cross;

}

FQUALIFIER
G4double GPMsc::GetCurrentLambda(double e) {

	G4double x = 0.0;
	if (useLambdaTable && theLambdaTable) {
		//    x = GPMsc::GetLambdaFromTable(e);
		//  return fFactor*((*theLambdaTable)[basedCoupleIndex])->Value(e);
		// basedCoupleIndex 0=gamma, 1=electron, 2=positron
		x = fFactor * theLambdaTable->physicsVectors[1].Value(e);
	} else {
		x = fFactor * CrossSectionPerVolume(e, 10.0, 100 * TeV);
	}
	return x;
}

FQUALIFIER
void GPMsc::ComputeIntegralLambda(G4double e) {

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
G4double GPMsc::AlongStepGetPhysicalInteractionLength(GPTrack* dp, double currentMinimalStep, GPGPILSelection* selection) {
  // from G4VMultipleScattering.cc

  // get Step limit proposed by the process
  *selection = NotCandidateForSelection;
  G4double x = currentMinimalStep;
  //  DefineMaterial(track.GetMaterialCutsCouple());
  G4double ekin = dp->E;
  //  if(isIon) { ekin *= proton_mass_c2/track.GetParticleDefinition()->GetPDGMass(); }
  //  currentModel = static_cast<G4VMscModel*>(SelectModel(ekin));

  // // define ionisation process
  // if(!currentModel->GetIonisation()) {
  //   currentModel->SetIonisation(G4LossTableManager::Instance()
  // 				->GetEnergyLossProcess(track.GetParticleDefinition()));
  // }  

  //  if(x > 0.0 && ekin > 0.0 && currentModel->IsActive(ekin)) {
  if(x > 0.0 && ekin > 0.0) {
    G4double tPathLength = ComputeTruePathLengthLimit(dp, x);
    if (tPathLength < x) { *selection = CandidateForSelection; }
    //    x = ComputeGeomPathLength(tPathLength);
    x = ComputeGeomPathLength();
  }

  return x;

}


FQUALIFIER G4double GPMsc::ComputeTruePathLengthLimit(GPTrack* dp, double currentMinimalStep) {
  // from G4UrbanMscModel95.cc

  tPathLength = currentMinimalStep;
  //  G4StepPoint* sp = track.GetStep()->GetPreStepPoint();
  //  G4StepStatus stepStatus = sp->GetStepStatus();
  GPStepStatus stepStatus = fAlongStepDoItProc;

  //  const G4DynamicParticle* dp = track.GetDynamicParticle();

  // if(stepStatus == fUndefined) {
  //   inside = false;
  //   insideskin = false;
  //   tlimit = geombig;
  //   SetParticle( dp->GetDefinition() );
  // }

  //  theLambdaTable = theTable;
  //  couple = track.GetMaterialCutsCouple();
  //  currentMaterialIndex = couple->GetIndex();
  currentKinEnergy = dp->E;
  currentRange = GetRange(currentKinEnergy);
  lambda0 = GetCurrentLambda(currentKinEnergy);

  // stop here if small range particle
  if(inside) { return tPathLength; }            
  
  if(tPathLength > currentRange) { tPathLength = currentRange; }

  // @@ dwjang : let me think how to handle this safety
  //  presafety = sp->GetSafety();
  presafety = 0.0;

  // G4cout << "G4Urban2::StepLimit tPathLength= " 
  //   <<tPathLength<<" safety= " << presafety
  //        << " range= " <<currentRange<< " lambda= "<<lambda0
  //   << " Alg: " << steppingAlgorithm <<G4endl;

  // far from geometry boundary
  if(currentRange < presafety)
    {
      inside = true;
      return tPathLength;  
    }

  // standard  version
  //
  if (steppingAlgorithm == fUseDistanceToBoundary)
    {
      //compute geomlimit and presafety 
      // @@ dwjang : need to implement later
      //      G4double geomlimit = ComputeGeomLimit(dp, presafety, currentRange);
      G4double geomlimit = geombig;

      // is it far from boundary ?
      if(currentRange < presafety)
	{
	  inside = true;
	  return tPathLength;   
	}

      smallstep += 1.;
      insideskin = false;

      if((stepStatus == fGeomBoundary) || (stepStatus == fUndefined))
        {
          rangeinit = currentRange;
          if(stepStatus == fUndefined) smallstep = 1.e10;
          else  smallstep = 1.;

          //define stepmin here (it depends on lambda!)
          //rough estimation of lambda_elastic/lambda_transport
          G4double rat = currentKinEnergy/MeV ;
          rat = 1.e-3/(rat*(10.+rat)) ;
          //stepmin ~ lambda_elastic
          stepmin = rat*lambda0;
          skindepth = skin*stepmin;
          //define tlimitmin
          tlimitmin = 10.*stepmin;
          if(tlimitmin < tlimitminfix) tlimitmin = tlimitminfix;
    //G4cout << "rangeinit= " << rangeinit << " stepmin= " << stepmin
    //   << " tlimitmin= " << tlimitmin << " geomlimit= " << geomlimit <<G4endl;
          // constraint from the geometry
          if((geomlimit < geombig) && (geomlimit > geommin))
            {
              // geomlimit is a geometrical step length
              // transform it to true path length (estimation)
              if((1.-geomlimit/lambda0) > 0.)
                geomlimit = -lambda0*log(1.-geomlimit/lambda0)+tlimitmin ;

              if(stepStatus == fGeomBoundary)
                tgeom = geomlimit/facgeom;
              else
                tgeom = 2.*geomlimit/facgeom;
            }
            else
              tgeom = geombig;
        }


      //step limit 
      tlimit = facrange*rangeinit;              

      //lower limit for tlimit
      if(tlimit < tlimitmin) tlimit = tlimitmin;

      if(tlimit > tgeom) tlimit = tgeom;

      //G4cout << "tgeom= " << tgeom << " geomlimit= " << geomlimit  
      //      << " tlimit= " << tlimit << " presafety= " << presafety << G4endl;

      // shortcut
      if((tPathLength < tlimit) && (tPathLength < presafety) &&
         (smallstep >= skin) && (tPathLength < geomlimit-0.999*skindepth))
  return tPathLength;   

      // step reduction near to boundary
      if(smallstep < skin)
  {
    tlimit = stepmin;
    insideskin = true;
  }
      else if(geomlimit < geombig)
  {
    if(geomlimit > skindepth)
      {
        if(tlimit > geomlimit-0.999*skindepth)
    tlimit = geomlimit-0.999*skindepth;
      }
    else
      {
        insideskin = true;
        if(tlimit > stepmin) tlimit = stepmin;
      }
  }

      if(tlimit < stepmin) tlimit = stepmin;

      // randomize 1st step or 1st 'normal' step in volume
      if((stepStatus == fUndefined) ||
         ((smallstep == skin) && !insideskin)) 
        { 
          G4double temptlimit = tlimit;
          if(temptlimit > tlimitmin)
          {
            do {
	      temptlimit = rand_normal_wrapper(devStates,threadId,tlimit,0.3*tlimit);        
               } while ((temptlimit < tlimitmin) || 
                        (temptlimit > 2.*tlimit-tlimitmin));
          }
          else
            temptlimit = tlimitmin;
          if(tPathLength > temptlimit) tPathLength = temptlimit;
        }
      else
        {  
          if(tPathLength > tlimit) tPathLength = tlimit  ; 
        }

    }
    // for 'normal' simulation with or without magnetic field 
    //  there no small step/single scattering at boundaries
  else if(steppingAlgorithm == fUseSafety)
    {
      // compute presafety again if presafety <= 0 and no boundary
      // i.e. when it is needed for optimization purposes

      // @@ dwjang : will implement this part later
      //      if((stepStatus != fGeomBoundary) && (presafety < tlimitminfix)) 
	//  presafety = ComputeSafety(sp->GetPosition(),tPathLength); 

      // is far from boundary
      if(currentRange < presafety)
        {
          inside = true;
          return tPathLength;

        }

      if((stepStatus == fGeomBoundary) || (stepStatus == fUndefined))
      {
        rangeinit = currentRange;
        fr = facrange;
        // 9.1 like stepping for e+/e- only (not for muons,hadrons)
        if(mass < masslimite) 
        {
          if(lambda0 > currentRange)
            rangeinit = lambda0;
          if(lambda0 > lambdalimit)
            fr *= 0.75+0.25*lambda0/lambdalimit;
        }

        //lower limit for tlimit
        G4double rat = currentKinEnergy/MeV ;
        rat = 1.e-3/(rat*(10.+rat)) ;
        tlimitmin = 10.*lambda0*rat;
        if(tlimitmin < tlimitminfix) tlimitmin = tlimitminfix;
      }
      //step limit
      tlimit = fr*rangeinit;               

      if(tlimit < facsafety*presafety)
        tlimit = facsafety*presafety;

      //lower limit for tlimit
      if(tlimit < tlimitmin) tlimit = tlimitmin;
      
      if(tPathLength > tlimit) tPathLength = tlimit;

    }
  
  // version similar to 7.1 (needed for some experiments)
  else
    {
      if (stepStatus == fGeomBoundary)
  {
    if (currentRange > lambda0) tlimit = facrange*currentRange;
    else                        tlimit = facrange*lambda0;

    if(tlimit < tlimitmin) tlimit = tlimitmin;
    if(tPathLength > tlimit) tPathLength = tlimit;
  }
    }
  //G4cout << "tPathLength= " << tPathLength 
  //   << " currentMinimalStep= " << currentMinimalStep << G4endl;
  return tPathLength ;

}


FQUALIFIER
G4double GPMsc::GetRange(double kinEnergy) {
  // from G4VMscModel.hh

  //  G4double localrange = DBL_MAX;
  // G4double q = part->GetPDGCharge()/eplus;
  // localrange = kinEnergy/(dedx*q*q*couple->GetMaterial()->GetDensity()); 

  G4double dedx = 2.0*MeV*cm2/g;
  G4double density = 5.16797e+19; // for electron
  G4double localrange = kinEnergy/(dedx*density);

  return localrange;

}


FQUALIFIER
G4double GPMsc::ComputeGeomPathLength() {

  lambdaeff = lambda0;
  par1 = -1. ;  
  par2 = par3 = 0. ;  

  //  do the true -> geom transformation
  zPathLength = tPathLength;

  // z = t for very small tPathLength
  if(tPathLength < tlimitminfix) return zPathLength;

  // this correction needed to run MSC with eIoni and eBrem inactivated
  // and makes no harm for a normal run
  if(tPathLength > currentRange)
    tPathLength = currentRange ;

  G4double tau   = tPathLength/lambda0 ;

  if ((tau <= tausmall) || insideskin) {
    zPathLength  = tPathLength;
    if(zPathLength > lambda0) zPathLength = lambda0;
    return zPathLength;
  }

  G4double zmean = tPathLength;
  if (tPathLength < currentRange*dtrl) {
    if(tau < taulim) zmean = tPathLength*(1.-0.5*tau) ;
    else             zmean = lambda0*(1.-exp(-tau));
    zPathLength = zmean ;
    return zPathLength;    

  } else if(currentKinEnergy < mass || tPathLength == currentRange)  {
    par1 = 1./currentRange ;
    par2 = 1./(par1*lambda0) ;
    par3 = 1.+par2 ;
    if(tPathLength < currentRange)
      zmean = (1.-exp(par3*log(1.-tPathLength/currentRange)))/(par1*par3) ;
    else {
      zmean = 1./(par1*par3) ;
    }
    zPathLength = zmean ;
    return zPathLength;    

  } else {
    //    G4double T1 = GetEnergy(particle,currentRange-tPathLength,couple);
    G4double T1 = currentKinEnergy;
    G4double lambda1 = GetCurrentLambda(T1);

    par1 = (lambda0-lambda1)/(lambda0*tPathLength) ;
    par2 = 1./(par1*lambda0) ;
    par3 = 1.+par2 ;
    zmean = (1.-exp(par3*log(lambda1/lambda0)))/(par1*par3) ;
  }

  zPathLength = zmean ;

  //  sample z
  if(samplez)
  {
    const G4double  ztmax = 0.999 ;
    G4double zt = zmean/tPathLength ;

    if (tPathLength > stepmin && zt < ztmax)              
    {
      G4double u,cz1;
      if(zt >= third)
      {
        G4double cz = 0.5*(3.*zt-1.)/(1.-zt) ;
        cz1 = 1.+cz ;
        G4double u0 = cz/cz1 ;
        G4double grej ;
        do {
            u = exp(log(rand_wrapper(devStates, threadId))/cz1) ;
            grej = exp(cz*log(u/u0))*(1.-u)/(1.-u0) ;
           } while (grej < rand_wrapper(devStates, threadId)) ;
      }
      else
      {
        u = 2.*zt*rand_wrapper(devStates, threadId);
      }
      zPathLength = tPathLength*u ;
    }
  }

  if(zPathLength > lambda0) { zPathLength = lambda0; }
  //G4cout << "zPathLength= " << zPathLength << " lambda1= " << lambda0 << G4endl;
  return zPathLength;

}


FQUALIFIER
G4double GPMsc::PostStepGetPhysicalInteractionLength(double kineticEnergy,
						   G4double previousStepSize, GPForceCondition* condition) {
  // from G4VMultipleScattering.cc

  *condition = Forced;
  return DBL_MAX;
}


FQUALIFIER
//G4double GPVEnergyLossProcess_PostStepDoIt(GPTrack* dp)
G4double GPMsc::AlongStepDoIt(GPTrack* dp) {
  // from G4VMultipleScattering.cc

  // if(currentModel->IsActive(track.GetKineticEnergy())) {
  //   fParticleChange.ProposeTrueStepLength(currentModel->ComputeTrueStepLength(step.GetStepLength()));
  // } else {
  //   fParticleChange.ProposeTrueStepLength(step.GetStepLength());
  // }
  // return &fParticleChange;

  G4double eloss = 0.0;

  //  ComputeTrueStepLength(stepLength);
  //  ProposeTrueStepLength();

  return eloss;
}


FQUALIFIER
//G4double GPVEnergyLossProcess_PostStepDoIt(GPTrack* dp)
G4double GPMsc::PostStepDoIt(GPTrack* dp) {
  // from G4VMultipleScattering.cc

  // fParticleChange.Initialize(track);
  // if(currentModel->IsActive(track.GetKineticEnergy())) {
  //   currentModel->SampleScattering(track.GetDynamicParticle(),
  // 				   step.GetPostStepPoint()->GetSafety());
  // }
  // return &fParticleChange;

  G4double eloss = 0.0;

  eloss = SampleScattering(dp);

  return eloss;
}


FQUALIFIER G4double GPMsc::SampleScattering(GPTrack* dp)
{
  return 0.0;

  // from G4UrbanMscModel95.cc
  /*
  G4double kineticEnergy = dynParticle->GetKineticEnergy();

  if((kineticEnergy <= 0.0) || (tPathLength <= tlimitminfix) ||
     (tPathLength/tausmall < lambda0)) return;

  G4double cth  = SampleCosineTheta(tPathLength,kineticEnergy);

  // protection against 'bad' cth values
  if(std::fabs(cth) > 1.) return;

  // extra protection agaist high energy particles backscattered 
  if(cth < 1.0 - 1000*tPathLength/lambda0 && kineticEnergy > 20*MeV) { 
    //G4cout << "Warning: large scattering E(MeV)= " << kineticEnergy 
    //     << " s(mm)= " << tPathLength/mm
    //     << " 1-cosTheta= " << 1.0 - cth << G4endl;
    // do Gaussian central scattering
    if(kineticEnergy > GeV && cth < 0.0) {
      G4ExceptionDescription ed;
      ed << dynParticle->GetDefinition()->GetParticleName()
   << " E(MeV)= " << kineticEnergy/MeV
   << " Step(mm)= " << tPathLength/mm
   << " in " << CurrentCouple()->GetMaterial()->GetName()
   << " CosTheta= " << cth 
   << " is too big - the angle is resampled" << G4endl;
      G4Exception("G4UrbanMscModel95::SampleScattering","em0004",
      JustWarning, ed,"");
    }
    do {
      cth = 1.0 + 2*log(G4UniformRand())*tPathLength/lambda0;
    } while(cth < -1.0);
  }

  G4double sth  = sqrt((1.0 - cth)*(1.0 + cth));
  G4double phi  = twopi*G4UniformRand();
  G4double dirx = sth*cos(phi);
  G4double diry = sth*sin(phi);

  G4ThreeVector oldDirection = dynParticle->GetMomentumDirection();
  G4ThreeVector newDirection(dirx,diry,cth);
  newDirection.rotateUz(oldDirection);
  fParticleChange->ProposeMomentumDirection(newDirection);

  if (latDisplasment && safety > tlimitminfix) {

    G4double r = SampleDisplacement();

    // G4cout << "G4UrbanMscModel95::SampleSecondaries: e(MeV)= " << kineticEnergy
    //  << " sinTheta= " << sth << " r(mm)= " << r
    //        << " trueStep(mm)= " << tPathLength
    //        << " geomStep(mm)= " << zPathLength
    //        << G4endl;

    if(r > 0.)
      {
        G4double latcorr = LatCorrelation();
        if(latcorr > r) latcorr = r;

        // sample direction of lateral displacement
        // compute it from the lateral correlation
        G4double Phi = 0.;
        if(std::abs(r*sth) < latcorr)
          Phi  = twopi*G4UniformRand();
        else
        {
          G4double psi = std::acos(latcorr/(r*sth));
          if(G4UniformRand() < 0.5)
            Phi = phi+psi;
          else
            Phi = phi-psi;
        }

        dirx = std::cos(Phi);
        diry = std::sin(Phi);

        G4ThreeVector latDirection(dirx,diry,0.0);
        latDirection.rotateUz(oldDirection);

  ComputeDisplacement(fParticleChange, latDirection, r, safety);
      }
  }
*/
}
