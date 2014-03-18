#include "GPBrem.h"
#include "GPRandom.h"
#include <stdio.h>

FQUALIFIER
void GPBrem::Print() {

	printf(
			"threadId(%d), useLambdaTable(%d), integral(%d), lambdaFactor(%lf), theLambdaTable(%ld), \
mfpKinEnergy(%lf), preStepKinEnergy(%lf), preStepLambda(%lf), fFactor(%lf), \
theNumberOfInteractionLengthLeft(%lf), currentInteractionLength(%lf), isElectron(%d), LPMFlag(%d)\n",
			threadId, int(useLambdaTable), int(integral), lambdaFactor,
			long(theLambdaTable), mfpKinEnergy, preStepKinEnergy, preStepLambda,
			fFactor, theNumberOfInteractionLengthLeft, currentInteractionLength,
			int(isElectron), int(LPMFlag));

}

FQUALIFIER GPBrem::GPBrem() {
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
void GPBrem::SetCurandStates(curandState* v) {
	devStates = v;
}

FQUALIFIER
void GPBrem::UseIntegral(bool v) {
	integral = v;
}

FQUALIFIER
void GPBrem::InitialiseStep(G4double kineticEnergy) {
	preStepKinEnergy = kineticEnergy;
	if (theNumberOfInteractionLengthLeft < 0.0)
		mfpKinEnergy = DBL_MAX;
}

FQUALIFIER
void GPBrem::SetLambdaTable(GPPhysicsTable* val) {
	theLambdaTable = val;
}

FQUALIFIER
void GPBrem::UseLambdaTable(bool v) {
	useLambdaTable = v;
}

FQUALIFIER
void GPBrem::ResetNumberOfInteractionLengthLeft() {
	theNumberOfInteractionLengthLeft = -log(rand_wrapper(devStates, threadId));
}

FQUALIFIER
void GPBrem::SubtractNumberOfInteractionLengthLeft(
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
G4double GPBrem::SupressionFunction(double kineticEnergy, double gammaEnergy,
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
G4double GPBrem::PositronCorrFactorSigma(double Z, double kineticEnergy,
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
G4double GPBrem::ComputeCrossSectionPerAtom(double kineticEnergy, double Z,
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
G4double GPBrem::CrossSectionPerVolume(double kineticEnergy,
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
G4double GPBrem::GetCurrentLambda(double e) {

	G4double x = 0.0;
	if (useLambdaTable && theLambdaTable) {
		//    x = GPBrem::GetLambdaFromTable(e);
		//  return fFactor*((*theLambdaTable)[basedCoupleIndex])->Value(e);
		// basedCoupleIndex 0=gamma, 1=electron, 2=positron
		x = fFactor * theLambdaTable->physicsVectors[1].Value(e);
	} else {
		x = fFactor * CrossSectionPerVolume(e, 10.0, 100 * TeV);
	}
	return x;
}

FQUALIFIER
void GPBrem::ComputeIntegralLambda(G4double e) {

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
G4double GPBrem::PostStepGetPhysicalInteractionLength(double kineticEnergy, G4double previousStepSize,GPForceCondition* condition) {
  // condition is set to "Not Forced"
   *condition = NotForced;

  G4double x = DBL_MAX;

  if (previousStepSize <= 0.0) {
    theNumberOfInteractionLengthLeft = -1.0;
  }

  InitialiseStep(kineticEnergy);

  preStepLambda = 0.0;
  preStepKinEnergy = kineticEnergy;
  mfpKinEnergy = DBL_MAX;

  // compute mean free path

  if (preStepKinEnergy < mfpKinEnergy) {

    if (integral) {
      ComputeIntegralLambda(preStepKinEnergy);
    } else {
      preStepLambda = GetCurrentLambda(preStepKinEnergy);
    }

    if (preStepLambda <= 0.0) {
      mfpKinEnergy = 0.0;
    }
  }

  // non-zero cross section
  if (preStepLambda > 0.0) {

    if (theNumberOfInteractionLengthLeft < 0.0) {
      // beggining of tracking (or just after DoIt of this process)
      ResetNumberOfInteractionLengthLeft();
    } else if (currentInteractionLength < DBL_MAX) {
      // subtract NumberOfInteractionLengthLeft
      SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if (theNumberOfInteractionLengthLeft < 0.)
	theNumberOfInteractionLengthLeft = perMillion;
    }

    // get mean free path and step limit
    currentInteractionLength = 1.0 / preStepLambda;
    x = theNumberOfInteractionLengthLeft * currentInteractionLength;

    // zero cross section case
  } else {

    if (theNumberOfInteractionLengthLeft > DBL_MIN
	&& currentInteractionLength < DBL_MAX) {

      // subtract NumberOfInteractionLengthLeft
      SubtractNumberOfInteractionLengthLeft(previousStepSize);
      if (theNumberOfInteractionLengthLeft < 0.)
	theNumberOfInteractionLengthLeft = perMillion;
    }
    currentInteractionLength = DBL_MAX;
  }

  return x;

}


FQUALIFIER
//G4double GPeBremsstrahlungModel_ScreenFunction1(double ScreenVariable)
G4double GPBrem::ScreenFunction1(double ScreenVariable) {
// compute the value of the screening function 3*PHI1 - PHI2
	G4double screenVal;

	if (ScreenVariable > 1.)
		screenVal = 42.24 - 8.368 * log(ScreenVariable + 0.952);
	else
		screenVal = 42.392 - ScreenVariable * (7.796 - 1.961 * ScreenVariable);

	return screenVal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FQUALIFIER
//G4double GPeBremsstrahlungModel_ScreenFunction2(double ScreenVariable)
G4double GPBrem::ScreenFunction2(double ScreenVariable) {
// compute the value of the screening function 1.5*PHI1 - 0.5*PHI2
	G4double screenVal;

	if (ScreenVariable > 1.)
		screenVal = 42.24 - 8.368 * log(ScreenVariable + 0.952);
	else
		screenVal = 41.734 - ScreenVariable * (6.484 - 1.250 * ScreenVariable);

	return screenVal;
}

FQUALIFIER
//G4double GPModifiedTsai_PolarAngle(double initial_energy, double final_energy, int Z) {
G4double GPBrem::PolarAngle(double initial_energy, double final_energy,
		int Z) {
	// Sample gamma angle (Z - axis along the parent particle).
	// Universal distribution suggested by L. Urban (Geant3 manual (1993) 
	// Phys211) derived from Tsai distribution (Rev Mod Phys 49,421(1977))

	G4double gamma = 1. + initial_energy / electron_mass_c2;
	G4double uMax = gamma * pi;

	const G4double a1 = 0.625;
	const G4double a2 = 1.875;
	const G4double border = 0.25;
	G4double u, theta;

	do {
		u = -log(
				rand_wrapper(devStates, threadId)
						* rand_wrapper(devStates, threadId));

		if (border > rand_wrapper(devStates, threadId))
			u /= a1;
		else
			u /= a2;

	} while (u > uMax);

	theta = u / gamma;

	return theta;

}

FQUALIFIER
G4double GPBrem::GPVEmModel_CrossSectionPerVolume(double ekin, double emin,
		G4double emax) {
	G4double cross = 0.0;
	const G4double theAtomNumDensity[nElements] = { 1.0958e+19, 1.0958e+19,
			4.3832e+19 };

	for (int i = 0; i < nElements; i++) {
		cross += theAtomNumDensity[i]
				* ComputeCrossSectionPerAtom(ekin, emin, emax);
		xsec[i] = cross;
	}
	return cross;
}

FQUALIFIER
//int GPVEmModel_SelectRandomAtom(G4double kinEnergy, double tcut, double tmax) {
int GPBrem::SelectRandomAtom(G4double kinEnergy, double tcut, double tmax) {

	int n = nElements - 1;
	if (n > 0) {
		G4double x = rand_wrapper(devStates, threadId)
				* GPVEmModel_CrossSectionPerVolume(kinEnergy, tcut, tmax);
		for (int i = 0; i < n; ++i) {
			if (x <= xsec[i]) {
				return i;
			}
		}
	}
	return n;
}

//void G4eBremsstrahlungModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp, 
//                                               const G4MaterialCutsCouple* couple,
//                                               const G4DynamicParticle* dp,
//                                               G4double tmin,
//                                               G4double maxEnergy)

FQUALIFIER
G4double GPBrem::SampleSecondaries(GPTrack* dp, double tmin,
		G4double maxEnergy)
// The emitted gamma energy is sampled using a parametrized formula 
// from L. Urban.
// This parametrization is derived from :
//    cross-section values of Seltzer and Berger for electron energies
//    1 keV - 10 GeV,
//    screened Bethe Heilter differential cross section above 10 GeV,
//    Migdal corrections in both case.
//  Seltzer & Berger: Nim B 12:95 (1985)
//  Nelson, Hirayama & Rogers: Technical report 265 SLAC (1985)
//  Migdal: Phys Rev 103:1811 (1956); Messel & Crawford: Pergamon Press (1970)
//
// A modified version of the random number techniques of Butcher&Messel is used
//    (Nuc Phys 20(1960),15).
		{
	G4double gammaEnergy = 0.0;
	G4double kineticEnergy = dp->E;
	G4double tmax = min(maxEnergy, kineticEnergy);
	if (tmin >= tmax) {
		return gammaEnergy;
	}

//
// GEANT4 internal units.
//
	const G4double ah10 = 4.67733E+00, ah11 = -6.19012E-01, ah12 = 2.02225E-02,
			ah20 = -7.34101E+00, ah21 = 1.00462E+00, ah22 = -3.20985E-02, ah30 =
					2.93119E+00, ah31 = -4.03761E-01, ah32 = 1.25153E-02;

	const G4double bh10 = 4.23071E+00, bh11 = -6.10995E-01, bh12 = 1.95531E-02,
			bh20 = -7.12527E+00, bh21 = 9.69160E-01, bh22 = -2.74255E-02, bh30 =
					2.69925E+00, bh31 = -3.63283E-01, bh32 = 9.55316E-03;

	const G4double al00 = -2.05398E+00, al01 = 2.38815E-02, al02 = 5.25483E-04,
			al10 = -7.69748E-02, al11 = -6.91499E-02, al12 = 2.22453E-03, al20 =
					4.06463E-02, al21 = -1.01281E-02, al22 = 3.40919E-04;

	const G4double bl00 = 1.04133E+00, bl01 = -9.43291E-03, bl02 = -4.54758E-04,
			bl10 = 1.19253E-01, bl11 = 4.07467E-02, bl12 = -1.30718E-03, bl20 =
					-1.59391E-02, bl21 = 7.27752E-03, bl22 = -1.94405E-04;

	const G4double tlow = 1. * MeV;

	bool LPMOK = false;

	//  const G4Material* material = couple->GetMaterial();

	//  // select randomly one element constituing the material
	//  const G4Element* anElement = SelectRandomAtom(couple);

	// moved to GPConstants.h
	// const int nElements = 3;
	// const G4double theAtomNumDensity[nElements] = {1.0958e+19, 1.0958e+19, 4.3832e+19};
	// const G4double Z[nElements] = {82, 74, 8};

	// tcut = 0.00099 MeV used above in GPBrem::GetCurrentLambda, taken from eBremsstralung
	int iElement = SelectRandomAtom(kineticEnergy, 0.00099, tmax);

	//  const G4double theAtomNumDensity[nElements] = {1.0958e+19, 1.0958e+19, 4.3832e+19};
	const G4double Z[nElements] = { 82, 74, 8 };

	// from G4IonisParamElm.hh
	//G4double fZ;                 // effective Z
	//G4double fZ3;                // std::pow (Z,1/3)
	//G4double fZZ3;               // std::pow (Z(Z+1),1/3)
	//G4double flogZ3;             // std::log(Z)/3

	// Extract Z factors for this Element
	G4double lnZ = 3. * pow(Z[iElement], 1 / 3.);
	G4double FZ = log(Z[iElement]) * (4. - 0.55 * log(Z[iElement]));
	G4double ZZ = pow(Z[iElement] * (Z[iElement] + 1), 1 / 3.);

	// limits of the energy sampling
	G4double totalEnergy = kineticEnergy + electron_mass_c2;

	G4double xmin = tmin / kineticEnergy;
	G4double xmax = tmax / kineticEnergy;
	G4double kappa = 0.0;
	if (xmax >= 1.) {
		xmax = 1.;
	} else {
		kappa = log(xmax) / log(xmin);
	}
	G4double epsilmin = tmin / totalEnergy;
	G4double epsilmax = tmax / totalEnergy;

	// Migdal factor
	G4double MigdalConstant = classic_electr_radius * electron_Compton_length
			* electron_Compton_length * 4.0 * pi;
	G4double MigdalFactor = (2.06011e+21) * MigdalConstant
			/ (epsilmax * epsilmax);

	G4double x, epsil, greject, migdal, grejmax, q;
	G4double U = log(kineticEnergy / electron_mass_c2);
	G4double U2 = U * U;

	// precalculated parameters
	G4double ah, bh;
	G4double screenfac = 0.0;

	if (kineticEnergy > tlow) {

		G4double ah1 = ah10 + ZZ * (ah11 + ZZ * ah12);
		G4double ah2 = ah20 + ZZ * (ah21 + ZZ * ah22);
		G4double ah3 = ah30 + ZZ * (ah31 + ZZ * ah32);

		G4double bh1 = bh10 + ZZ * (bh11 + ZZ * bh12);
		G4double bh2 = bh20 + ZZ * (bh21 + ZZ * bh22);
		G4double bh3 = bh30 + ZZ * (bh31 + ZZ * bh32);

		ah = 1. + (ah1 * U2 + ah2 * U + ah3) / (U2 * U);
		bh = 0.75 + (bh1 * U2 + bh2 * U + bh3) / (U2 * U);

		// limit of the screening variable
		screenfac = 136. * electron_mass_c2
				/ (pow(Z[iElement], 1 / 3.) * totalEnergy);
		G4double screenmin = screenfac * epsilmin / (1. - epsilmin);

		// Compute the maximum of the rejection function
		G4double F1 = max(ScreenFunction1(screenmin) - FZ, 0.);
		G4double F2 = max(ScreenFunction2(screenmin) - FZ, 0.);
		grejmax = (F1 - epsilmin * (F1 * ah - bh * epsilmin * F2))
				/ (42.392 - FZ);

	} else {

		G4double al0 = al00 + ZZ * (al01 + ZZ * al02);
		G4double al1 = al10 + ZZ * (al11 + ZZ * al12);
		G4double al2 = al20 + ZZ * (al21 + ZZ * al22);

		G4double bl0 = bl00 + ZZ * (bl01 + ZZ * bl02);
		G4double bl1 = bl10 + ZZ * (bl11 + ZZ * bl12);
		G4double bl2 = bl20 + ZZ * (bl21 + ZZ * bl22);

		ah = al0 + al1 * U + al2 * U2;
		bh = bl0 + bl1 * U + bl2 * U2;

		// Compute the maximum of the rejection function
		grejmax = max(1. + xmin * (ah + bh * xmin), 1. + ah + bh);
		G4double xm = -ah / (2. * bh);
		if (xmin < xm && xm < xmax)
			grejmax = max(grejmax, 1. + xm * (ah + bh * xm));
	}

	//
	//  sample the energy rate of the emitted gamma for e- kin energy > 1 MeV
	//

	do {
		if (kineticEnergy > tlow) {
			do {
				q = rand_wrapper(devStates, threadId);
				x = pow(xmin, q + kappa * (1.0 - q));
				epsil = x * kineticEnergy / totalEnergy;
				G4double screenvar = screenfac * epsil / (1.0 - epsil);
				G4double F1 = max(ScreenFunction1(screenvar) - FZ, 0.);
				G4double F2 = max(ScreenFunction2(screenvar) - FZ, 0.);
				migdal = (1. + MigdalFactor) / (1. + MigdalFactor / (x * x));
				greject = migdal * (F1 - epsil * (ah * F1 - bh * epsil * F2))
						/ (42.392 - FZ);
			} while (greject < rand_wrapper(devStates, threadId) * grejmax);

		} else {

			do {
				q = rand_wrapper(devStates, threadId);
				x = pow(xmin, q + kappa * (1.0 - q));
				migdal = (1. + MigdalFactor) / (1. + MigdalFactor / (x * x));
				greject = migdal * (1. + x * (ah + bh * x));
			} while (greject < rand_wrapper(devStates, threadId) * grejmax);
		}
		gammaEnergy = x * kineticEnergy;

		if (LPMFlag) {
			// take into account the supression due to the LPM effect
			if (rand_wrapper(devStates, threadId)
					<= SupressionFunction(kineticEnergy, gammaEnergy, LPMFlag))
				LPMOK = true;
		} else
			LPMOK = true;

	} while (!LPMOK);

	/*

	 GPThreeVector direction;
	 direction.x = dp->px;
	 direction.y = dp->py;
	 direction.z = dp->pz;

	 //
	 // angles of the emitted gamma. ( Z - axis along the parent particle)
	 // use general interface
	 //
	 G4double theta = GPModifiedTsai_PolarAngle( totalEnergy,totalEnergy-gammaEnergy,int(Z[iElement]));

	 G4double sint = sin(theta);

	 G4double phi = twopi * G4UniformRand() ;

	 GPThreeVector gammaDirection;
	 gammaDirection.x = sint*cos(phi);
	 gammaDirection.y = sint*sin(phi);
	 gammaDirection.z = cos(theta);
	 GPThreeVector_rotateUz(&gammaDirection,direction);

	 // create G4DynamicParticle object for the Gamma
	 //  G4DynamicParticle* g = new G4DynamicParticle(theGamma,gammaDirection,gammaEnergy);
	 GPTrack g;
	 g.px = gammaDirection.x;
	 g.py = gammaDirection.y;
	 g.pz = gammaDirection.z;
	 g.E = gammaEnergy;

	 //  vdp->push_back(g);
	 
	 G4double totMomentum = sqrt(kineticEnergy*(totalEnergy + electron_mass_c2));
	 GPThreeVector dir = GPThreeVector_create(totMomentum*direction.x - gammaEnergy*gammaDirection.x,
	 totMomentum*direction.y - gammaEnergy*gammaDirection.y,
	 totMomentum*direction.z - gammaEnergy*gammaDirection.z);
	 dir = GPThreeVector_unit(dir);
	 direction.x = dir.x;
	 direction.y = dir.y;
	 direction.z = dir.z;

	 // energy of primary
	 G4double finalE = kineticEnergy - gammaEnergy;

	 // stop tracking and create new secondary instead of primary
	 if(gammaEnergy > SecondaryThreshold()) { // default bremTh is DBL_MAX
	 fParticleChange->ProposeTrackStatus(fStopAndKill);
	 fParticleChange->SetProposedKineticEnergy(0.0);
	 G4DynamicParticle* el = 
	 new G4DynamicParticle(const_cast<G4ParticleDefinition*>(particle),
	 direction, finalE);
	 vdp->push_back(el);

	 // continue tracking
	 } else {
	 fParticleChange->SetProposedMomentumDirection(direction);
	 fParticleChange->SetProposedKineticEnergy(finalE);
	 }

	 */

	return gammaEnergy;

}

FQUALIFIER
//G4double GPVEnergyLossProcess_PostStepDoIt(GPTrack* dp)
G4double GPBrem::PostStepDoIt(GPTrack* dp) {
	// In all cases clear number of interaction lengths
	theNumberOfInteractionLengthLeft = -1.0;

	G4double eloss = 0.0;

	//  fParticleChange.InitializeForPostStep(track);
	G4double finalT = dp->E;
	G4double lowestKinEnergy = 1.0 * MeV;
	if (finalT <= lowestKinEnergy) {
		return eloss;
	}

	G4double massRatio = 1.0;
	G4double postStepScaledEnergy = finalT * massRatio;

	//  if(!currentModel->IsActive(postStepScaledEnergy)) { return &fParticleChange; }

//  if(-1 < verboseLevel) {
//    G4cout << GetProcessName()
//           << "::PostStepDoIt: E(MeV)= " << finalT/MeV
//           << G4endl;
//  }

	// forced process - should happen only once per track
	//  if(biasFlag) {
	//    if(biasManager->ForcedInteractionRegion(currentCoupleIndex)) {
	//      biasFlag = false;
	//    }
	//  }

	// Integral approach
	if (integral) {
		//    G4double lx = GetLambdaForScaledEnergy(postStepScaledEnergy);
		G4double lx = GetCurrentLambda(postStepScaledEnergy);

		if (lx <= 0.0) {
			return 0;
		} else if (preStepLambda * rand_wrapper(devStates, threadId) > lx) {
			return 0;
		}
	}

	G4double tcut = 0.00099;
	eloss = SampleSecondaries(dp, tcut, 100 * TeV);

	return eloss;
}
