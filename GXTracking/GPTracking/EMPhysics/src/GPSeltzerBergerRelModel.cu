#include "GPSeltzerBergerRelModel.h"
#include "GPRandom.h"
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
  #include <iostream>
#endif
#endif

FQUALIFIER
GPSeltzerBergerRelModel::GPSeltzerBergerRelModel(curandState* devStates,
						 int threadId,
						 GPPhysics2DVector* data)

{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

  //
  fCurrentElement = 0;
  particleMass = electron_mass_c2; 
  kinEnergy=0.0;
  totalEnergy=0.0;
  currentZ=0.0;
  densityFactor=0.0;
  densityCorr=0.0;
  lpmEnergy=0.0;

  z13 =0.0;
  z23 =0.0;
  lnZ =0.0;
  Fel = 0.0;
  Finel = 0.0;
  fCoulomb =0.0;
  fMax = 0.0; 

  lowKinEnergy  = 0.1*eV;

  //model specific

  //cross section data for G4SeltzerBergerModel
  dataSB = data;

  //secondary from either G4SeltzerBergerModel or G4eBremsstrahlungRelModel
  //theGamma = (GXTrack*) malloc(sizeof(GXTrack));

  //G4VEmModel
  pParticleChange = 0;
  energyLimitModels = 1.0*GeV; //model boundary from G4eBremsstrahlung
  secondaryThreshold = DBL_MAX;

  lowLimit = 0.1*keV;
  highLimit = 100.0*TeV;
  eMinActive = 0.0;
  eMaxActive = DBL_MAX;
}

FQUALIFIER GPSeltzerBergerRelModel::~GPSeltzerBergerRelModel()
{
}

//---------------------------------------------------------------------------
//
// G4eBremsstrahlungRelModel
//
//---------------------------------------------------------------------------

FQUALIFIER
GPElement* GPSeltzerBergerRelModel::SelectRandomAtom(GPMaterial* material,
						       G4double kinEnergy,
						       G4double tcut,
						       G4double tmax)
{
  GPElement* theElementVector = GPMaterial_GetElementVector(material);
  G4int n = GPMaterial_GetNumberOfElements(material) - 1;
  fCurrentElement = &(theElementVector[n]);

  if (n > 0) {
    G4double x = GPUniformRand(fDevStates, fThreadId)*
      CrossSectionPerVolume(material,kinEnergy,
			    tcut,tmax);
    for(G4int i=0; i<n; ++i) {
      //@@@CUDA high level of divergence branch
      if (x <= xsec[i]) {
	fCurrentElement = &(theElementVector[i]);
	break;
      }
    }
  }
  return fCurrentElement;
}

G4double GPSeltzerBergerRelModel::CrossSectionPerVolume(GPMaterial* material,
							  G4double ekin,
							  G4double emin,
							  G4double emax)
{
  SetupForMaterial(material, ekin);

  G4double cross = 0.0;
   GPElement* theElementVector = GPMaterial_GetElementVector(material);
  G4double* theAtomNumDensityVector = 
    GPMaterial_GetVecNbOfAtomsPerVolume(material);
  G4int nelm = GPMaterial_GetNumberOfElements(material);

  for (G4int i=0; i<nelm; i++) {
    G4double Z = GPElement_GetZ(&(theElementVector[i]));
    cross += theAtomNumDensityVector[i]*
      ComputeCrossSectionPerAtom(ekin,Z,emin,emax);
    xsec[i] = cross;
  }

  return cross;
}

FQUALIFIER
void GPSeltzerBergerRelModel::SetupForMaterial(GPMaterial* mat, 
						 G4double kineticEnergy)
{
  G4double fMigdalConstant = 
    classic_electr_radius*electron_Compton_length*electron_Compton_length*4.0*pi;
  G4double fLPMconstant = 
    fine_structure_const*electron_mass_c2*electron_mass_c2/(4.*pi*hbarc)*0.5;

  densityFactor = GPMaterial_GetElectronDensity(mat)*fMigdalConstant;
  lpmEnergy = GPMaterial_GetRadlen(mat)*fLPMconstant;

#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
  std::cout<<"GPSeltzer::SetupForMaterial: fMigdalC="<< fMigdalConstant
	<<", fLPM="<< fLPMconstant
	<<", densityFac="<< densityFactor
	<<", lpmEnergy="<< lpmEnergy << std::endl;
#endif
#endif

  // Threshold for LPM effect (i.e. below which LPM hidden by density effect) 
  //@@@@ LPMFlag() is true for RelModel and false for SeltzerBerger, 
  //so use 1 GeV 
  //  if (LPMFlag()) {
  if (kineticEnergy > energyLimitModels) {
    energyThresholdLPM=sqrt(densityFactor)*lpmEnergy;
  } else {
    energyThresholdLPM=1.e39;   // i.e. do not use LPM effect
  }
  // calculate threshold for density effect
  kinEnergy   = kineticEnergy;
  totalEnergy = kineticEnergy + particleMass;
  densityCorr = densityFactor*totalEnergy*totalEnergy;

  // define critical gamma energies (important for integration/dicing)
  klpm=totalEnergy*totalEnergy/lpmEnergy;
  kp=sqrt(densityCorr);

#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
  std::cout<<"GPSeltzer::SetupForMaterial: EneThreshLPM="<< energyThresholdLPM
	<<", kinEne="<< kinEnergy
	<<", totEne="<< totalEnergy
	<<", densityCorr="<< densityCorr
	<<", klpm="<< klpm
	<<", kp="<< kp
	<< std::endl;
#endif
#endif
}

FQUALIFIER
G4double GPSeltzerBergerRelModel::ComputeCrossSectionPerAtom(
                                                       G4double kineticEnergy, 
						       G4double Z, 
						       G4double cutEnergy, 
						       G4double maxEnergy)
{
  if(kineticEnergy < lowKinEnergy) { return 0.0; }

  G4double cut  = fmin(cutEnergy, kineticEnergy);
  G4double tmax = fmin(maxEnergy, kineticEnergy);

  if(cut >= tmax) { return 0.0; }

  SetCurrentElement(Z);

  G4double cross = ComputeXSectionPerAtom(cut);

  // allow partial integration
  if(tmax < kinEnergy) { cross -= ComputeXSectionPerAtom(tmax); }
  
  G4double bremFactor = 
    fine_structure_const*classic_electr_radius*classic_electr_radius*(16./3.);

  cross *= Z*Z*bremFactor;
  
  return cross;
}

FQUALIFIER
void GPSeltzerBergerRelModel::SetCurrentElement(G4double Z)
{
  if(Z != currentZ) {
    currentZ = Z;

    G4int iz = G4int(Z);

    G4double x  = G4double(iz);
    //    z13 = nist->GetZ13(iz);
    //        = g4pow->Z13(iz);
    //        = pz13[Z];
    //        = pow(x,1.0/3.0);

    z13 = pow(x,1.0/3.0);

    z23 = z13*z13;
    //    lnZ = nist->GetLOGZ(iz);
    //        = g4pow->logZ(iz);
    //        = lz[Z];
    //        = log(x);
    lnZ = log(x);

    const G4double Fel_light[5] = {0., 5.31  , 4.79  , 4.74 ,  4.71};
    const G4double Finel_light[5] = {0., 6.144 , 5.621 , 5.805 , 5.924};

    G4double Fel = 0.0;
    G4double Finel = 0.0;

    if (iz <= 4) {
      Fel = Fel_light[iz];  
      Finel = Finel_light[iz] ; 
    }
    else {
      //    const G4double facFel = log(184.15);
      //    const G4double facFinel = log(1194.);
      //    Fel = facFel - lnZ/3. ;
      //    Finel = facFinel - 2.*lnZ/3. ;
      Fel = log(184.15) - lnZ/3. ;
      Finel = log(1194.) - 2.*lnZ/3. ;
    }

    //    fCoulomb = GetCurrentElement()->GetfCoulomb();
    fCoulomb = GPElement_GetfCoulomb(fCurrentElement);
    fMax = Fel-fCoulomb + Finel/currentZ  +  (1.+1./currentZ)/12.;
  }
}

FQUALIFIER
G4double GPSeltzerBergerRelModel::ComputeXSectionPerAtom(G4double cut)
{
  G4double cross = 0.0;

  // number of intervals and integration step 
  G4double vcut = log(cut/totalEnergy);
  G4double vmax = log(kinEnergy/totalEnergy);
  G4int n = (G4int)(0.45*(vmax - vcut)) + 4;
  //  n=1; //  integration test 
  G4double delta = (vmax - vcut)/G4double(n);

  G4double e0 = vcut;
  G4double xs =0.0; 

  // integration
  G4double xgi[8]= {0.0199, 0.1017, 0.2372, 0.4083,
		    0.5917, 0.7628, 0.8983, 0.9801};
  G4double wgi[8]= {0.0506, 0.1112, 0.1569, 0.1813,
		    0.1813, 0.1569, 0.1112, 0.0506};

  for(G4int l=0; l<n; l++) {

    for(G4int i=0; i<8; i++) {

      G4double eg = exp(e0 + xgi[i]*delta)*totalEnergy;

      if(totalEnergy > energyThresholdLPM) {
	xs = ComputeRelDXSectionPerAtom(eg);
      } else {
	//	xs = ComputeDXSectionPerAtom(eg);
        if(kinEnergy > 1.0*GeV) {
          xs = ComputeDXSectionPerAtom_eBremsstrahlungRelModel(eg);
        }
        else {
          xs = ComputeDXSectionPerAtom_SeltzerBergerModel(eg);
        }
      }
      cross += wgi[i]*xs/(1.0 + densityCorr/(eg*eg));
    }
    e0 += delta;
  }

  cross *= delta;

  return cross;
}

FQUALIFIER
G4double GPSeltzerBergerRelModel::ComputeRelDXSectionPerAtom(
                                                       G4double gammaEnergy)
// Ultra relativistic model
//   only valid for very high energies, but includes LPM suppression
//    * complete screening
{
  if(gammaEnergy < 0.0) { return 0.0; }

  G4double y = gammaEnergy/totalEnergy;
  G4double y2 = y*y*.25;
  G4double yone2 = (1.-y+2.*y2);

  // ** form factors complete screening case **      

  // ** calc LPM functions: include ter-mikaelian merging with density effect **
  //  G4double xiLPM, gLPM, phiLPM;  // to be made member variables !!! 

  CalcLPMFunctions(gammaEnergy);

  G4double mainLPM = 
    xiLPM*(y2 * gLPM + yone2*phiLPM) * ( (Fel-fCoulomb) + Finel/currentZ );
  G4double secondTerm = (1.-y)/12.*(1.+1./currentZ);

  G4double cross = mainLPM+secondTerm;
  return cross;
}

FQUALIFIER G4double 
GPSeltzerBergerRelModel::ComputeDXSectionPerAtom_eBremsstrahlungRelModel(
                                                 G4double gammaEnergy)
// Relativistic model
//  only valid for high energies (and if LPM suppression does not play a role)
//  * screening according to thomas-fermi-Model (only valid for Z>5)
//  * no LPM effect
{
  if(gammaEnergy < 0.0) { return 0.0; }

  G4double y = gammaEnergy/totalEnergy;

  G4double main=0.,secondTerm=0.;

  //@@@ use_completescreening(false) in G4eBremsstrahlungRelModel constructor
  //if (use_completescreening|| currentZ<5) {
  if (currentZ<5) {
    // ** form factors complete screening case **      
    main   = (3./4.*y*y - y + 1.) * ( (Fel-fCoulomb) + Finel/currentZ );
    secondTerm = (1.-y)/12.*(1.+1./currentZ);
  }
  else {
    // ** intermediate screening using Thomas-Fermi FF from Tsai 
    // only valid for Z>=5** 
    G4double dd=100.*electron_mass_c2*y/(totalEnergy-gammaEnergy);
    G4double gg=dd/z13;
    G4double eps=dd/z23;
    //    G4double phi1=Phi1(gg,currentZ);
    //    G4double phi1m2=Phi1M2(gg,currentZ);
    //    G4double psi1=Psi1(eps,currentZ);
    //    G4double psi1m2=Psi1M2(eps,currentZ);
    G4double phi1=Phi1(gg);
    G4double phi1m2=Phi1M2(gg);
    G4double psi1=Psi1(eps);
    G4double psi1m2=Psi1M2(eps);
    main   = (3./4.*y*y - y + 1.) * ( (0.25*phi1-1./3.*lnZ-fCoulomb) + 
				      (0.25*psi1-2./3.*lnZ)/currentZ );
    secondTerm = (1.-y)/8.*(phi1m2+psi1m2/currentZ);
  }
  G4double cross = main+secondTerm;
  return cross;
}

FQUALIFIER
void  GPSeltzerBergerRelModel::CalcLPMFunctions(G4double k)
{
  // *** calculate lpm variable s & sprime ***
  // Klein eqs. (78) & (79)

  G4double sprime = sqrt(0.125*k*lpmEnergy/(totalEnergy*(totalEnergy-k)));

  //  G4double s1 = preS1*z23;
  G4double s1 = (1./(184.15*184.15))*z23;
  //  G4double logS1 = 2./3.*lnZ-2.*facFel;
  G4double logS1 = 2./3.*lnZ-2.*log(184.15);
  //  G4double logTS1 = logTwo+logS1;
  G4double logTS1 = log(2.)+logS1;

  xiLPM = 2.;

  if (sprime>1) 
    xiLPM = 1.;
  else if (sprime>sqrt(2.)*s1) {
    G4double h  = log(sprime)/logTS1;
    //    xiLPM = 1+h-0.08*(1-h)*(1-sqr(1-h))/logTS1;
    xiLPM = 1+h-0.08*(1-h)*(1-(1-h)*(1-h))/logTS1;
  }

  G4double s0 = sprime/sqrt(xiLPM); 

  // *** merging with density effect***  should be only necessary in region 
  // "close to" kp, e.g. k<100*kp using Ter-Mikaelian eq. (20.9)
  G4double k2 = k*k;
  s0 *= (1 + (densityCorr/k2) );

  // recalculate Xi using modified s above
  // Klein eq. (75)
  xiLPM = 1.;
  if (s0<=s1) xiLPM = 2.;
  else if ( (s1<s0) && (s0<=1) ) xiLPM = 1. + log(s0)/logS1;
  
  // *** calculate supression functions phi and G ***
  // Klein eqs. (77)
  G4double s2=s0*s0;
  G4double s3=s0*s2;
  G4double s4=s2*s2;

  if (s0<0.1) {
    // high suppression limit
    phiLPM = 6.*s0 - 18.84955592153876*s2 + 39.47841760435743*s3 
      - 57.69873135166053*s4;
    gLPM = 37.69911184307752*s2 - 236.8705056261446*s3 + 807.7822389*s4;
  }
  else if (s0<1.9516) {
    // intermediate suppression
    // using eq.77 approxim. valid s<2.      
    phiLPM = 1.-exp(-6.*s0*(1.+(3.-pi)*s0)
                +s3/(0.623+0.795*s0+0.658*s2));
    if (s0<0.415827397755) {
      // using eq.77 approxim. valid 0.07<s<2
      G4double psiLPM = 1-exp(-4*s0-8*s2/(1+3.936*s0+4.97*s2-0.05*s3+7.50*s4));
      gLPM = 3*psiLPM-2*phiLPM;
    }
    else {
      // using alternative parametrisiation
      G4double pre = 
	-0.16072300849123999 + s0*3.7550300067531581 + s2*-1.7981383069010097 
        + s3*0.67282686077812381 + s4*-0.1207722909879257;
      gLPM = tanh(pre);
    }
  }
  else {
    // low suppression limit valid s>2.
    phiLPM = 1. - 0.0119048/s4;
    gLPM = 1. - 0.0230655/s4;
  }

  // *** make sure suppression is smaller than 1 ***
  // *** caused by Migdal approximation in xi    ***
  if (xiLPM*phiLPM>1. || s0>0.57)  { xiLPM=1./phiLPM; }
}

FQUALIFIER
G4double GPSeltzerBergerRelModel::Phi1(G4double gg)
{
  // Thomas-Fermi FF from Tsai, eq.(3.38) for Z>=5
  //  return 20.863 - 2.*log(1. + sqr(0.55846*gg) )
  return 20.863 - 2.*log(1. + (0.55846*gg)*(0.55846*gg) )
    - 4.*( 1. - 0.6*exp(-0.9*gg) - 0.4*exp(-1.5*gg) );
}

FQUALIFIER
G4double GPSeltzerBergerRelModel::Phi1M2(G4double gg)
{
  // Thomas-Fermi FF from Tsai, eq. (3.39) for Z>=5
  // return Phi1(gg,Z) - 
  return 2./(3.*(1. + 6.5*gg +6.*gg*gg) );
}

FQUALIFIER
G4double GPSeltzerBergerRelModel::Psi1(G4double eps)
{
  // Thomas-Fermi FF from Tsai, eq.(3.40) for Z>=5 
  //  return 28.340 - 2.*log(1. + sqr(3.621*eps) )
  return 28.340 - 2.*log(1. + (3.621*eps)*(3.621*eps) )
    - 4.*( 1. - 0.7*exp(-8*eps) - 0.3*exp(-29.*eps) );
}

FQUALIFIER
G4double GPSeltzerBergerRelModel::Psi1M2(G4double eps)
{
  // Thomas-Fermi FF from Tsai, eq. (3.41) for Z>=5
  return  2./(3.*(1. + 40.*eps +400.*eps*eps) );
}

FQUALIFIER void 
GPSeltzerBergerRelModel::SampleSecondaries_eBremsstrahlungRelModel(
                                             GXTrack* track,
					     GPMaterial* material,
					     G4double cutEnergy,
					     G4double maxEnergy)
{
  G4double kineticEnergy = track->E ; //dp->GetKineticEnergy();
  if(kineticEnergy < lowKinEnergy) { return; }
  G4double cut  = fmin(cutEnergy, kineticEnergy);
  G4double emax = fmin(maxEnergy, kineticEnergy);
  if(cut >= emax) { return; }

  //  SetupForMaterial(particle, couple->GetMaterial(), kineticEnergy);
  SetupForMaterial(material, kineticEnergy);

  GPElement* elm = SelectRandomAtom(material,kineticEnergy,cut,emax);

  SetCurrentElement(GPElement_GetZ(elm));

  kinEnergy   = kineticEnergy;
  totalEnergy = kineticEnergy + particleMass;
  densityCorr = densityFactor*totalEnergy*totalEnergy;

  //G4double fmax= fMax;
  G4bool highe = true;
  if(totalEnergy < energyThresholdLPM) { highe = false; }
 
  G4double xmin = log(cut*cut + densityCorr);
  G4double xmax = log(emax*emax  + densityCorr);
  G4double gammaEnergy, f, x; 

  do {
    //x = exp(xmin + G4UniformRand()*(xmax - xmin)) - densityCorr;
    x = exp(xmin + GPUniformRand(fDevStates,fThreadId)*(xmax-xmin))-densityCorr;

    if(x < 0.0) { x = 0.0; }
    gammaEnergy = sqrt(x);
    if(highe) { f = ComputeRelDXSectionPerAtom(gammaEnergy); }
    else      { //  f = ComputeDXSectionPerAtom(gammaEnergy); 
      f = ComputeDXSectionPerAtom_eBremsstrahlungRelModel(gammaEnergy); 
    }

    if ( f > fMax ) {
      ; //warning counter
    }
    //    ax = fMax*GPUniformRand(fDevStates, fThreadId);
    //  } while (f < ax);
  } while (f < fMax*GPUniformRand(fDevStates, fThreadId));
  //  } while (f < fMax*G4UniformRand());

  //
  // angles of the emitted gamma. ( Z - axis along the parent particle)
  // use general interface
  //

  //  G4ThreeVector gammaDirection = 
  //    GetAngularDistribution()->SampleDirection(dp, totalEnergy-gammaEnergy,
  //                                              G4lrint(currentZ),
  //                                              couple->GetMaterial());
  GPThreeVector refDirection = 
    GPThreeVector_unit(GPThreeVector_create(track->px,track->py,track->pz));
  GPThreeVector gammaDirection = 
    DipBustGenerator_SampleDirection(kineticEnergy,refDirection);

  // create G4DynamicParticle object for the Gamma
  //  G4DynamicParticle* gamma = 
  //    new G4DynamicParticle(theGamma,gammaDirection,gammaEnergy);
  //  vdp->push_back(gamma);

  //@@@G4FWP fill secondary (one photon)
   FillSecondary(track, gammaDirection, gammaEnergy,0);

  //update electron GXTrack

  G4double totMomentum = sqrt(kineticEnergy*(totalEnergy + electron_mass_c2));
  //  G4ThreeVector direction = (totMomentum*dp->GetMomentumDirection()
  //                             - gammaEnergy*gammaDirection).unit();
  GPThreeVector direction = GPThreeVector_unit(GPThreeVector_sub(
			    GPThreeVector_mult(refDirection,totMomentum),
			    GPThreeVector_mult(gammaDirection,gammaEnergy)));

  // energy of primary
  G4double finalE = kineticEnergy - gammaEnergy;

  // stop tracking and create new secondary instead of primary
  if(gammaEnergy > SecondaryThreshold()) {
    pParticleChange->ProposeTrackStatus(fStopAndKill);
    pParticleChange->SetProposedKineticEnergy(0.0);
    //  G4DynamicParticle* el = 
    //    new G4DynamicParticle(const_cast<G4ParticleDefinition*>(particle),
    //                          direction, finalE);
    //  vdp->push_back(el);
    //
    // continue tracking
  } else {
    pParticleChange->SetProposedMomentumDirection(direction);
    pParticleChange->SetProposedKineticEnergy(finalE);
  }
}

//---------------------------------------------------------------------------
//
// G4SeltzerBergerModel
//
//---------------------------------------------------------------------------

FQUALIFIER
G4double GPSeltzerBergerRelModel::ComputeDXSectionPerAtom_SeltzerBergerModel(
                                                 G4double gammaEnergy)
{
  if(gammaEnergy < 0.0 || kinEnergy <= 0.0) { return 0.0; }
  G4double x = gammaEnergy/kinEnergy;
  G4double y = log(kinEnergy/MeV);
  G4int Z = G4int(currentZ);

  //  if(!dataSB[Z]) { ReadData(Z); }
  G4double bremFactor = 
    fine_structure_const*classic_electr_radius*classic_electr_radius*(16./3.);

  G4double invb2 = 
    totalEnergy*totalEnergy/(kinEnergy*(kinEnergy + 2*particleMass));
  G4double cross = dataSB[Z].Value(x,y)*invb2*millibarn/bremFactor;

  //@@@positron is not included yet
  /*
  if(!isElectron) {
    G4double invbeta1 = sqrt(invb2);
    G4double e2 = kinEnergy - gammaEnergy;
    if(e2 > 0.0) {
      G4double invbeta2 = (e2 + particleMass)/sqrt(e2*(e2 + 2*particleMass));
      G4double xxx = twopi*fine_structure_const*currentZ*(invbeta1 - invbeta2);
      if(xxx < expnumlim) { cross = 0.0; }
      else { cross *= exp(xxx); }
    } else {
      cross = 0.0;
    }
    }
  */  
  return cross;
}

FQUALIFIER void 
GPSeltzerBergerRelModel::SampleSecondaries_SeltzerBergerModel(GXTrack* track,
                                                 GPMaterial* material,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double kineticEnergy = track->E ; //dp->GetKineticEnergy();
  G4double cut  = fmin(cutEnergy, kineticEnergy);
  G4double emax = fmin(maxEnergy, kineticEnergy);
  if(cut >= emax) { return; }
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
  std::cout<<"SB model: kinE="<< kineticEnergy <<" cut="<< cut <<" emax="<< emax <<".  Calling SetupForMaterial..."<< std::endl;
#endif
#endif

  //  SetupForMaterial(particle, couple->GetMaterial(), kineticEnergy);
  SetupForMaterial(material, kineticEnergy);

  GPElement* elm = SelectRandomAtom(material,kineticEnergy,cut,emax);
  SetCurrentElement(GPElement_GetZ(elm));
  G4int Z = G4int(currentZ);

  totalEnergy = kineticEnergy + particleMass;
  densityCorr = densityFactor*totalEnergy*totalEnergy;
  G4double totMomentum = sqrt(kineticEnergy*(totalEnergy + electron_mass_c2));

  G4double xmin = log(cut*cut + densityCorr);
  G4double xmax = log(emax*emax  + densityCorr);
  G4double y = log(kineticEnergy/MeV);

#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
  std::cout<<"SB model: Z="<< Z <<", totE="<< totalEnergy <<", densCorr="<< densityCorr
	   <<" xmin="<< xmin <<" xmax="<< xmax <<", y="<< y << std::endl;
#endif
#endif

  G4double gammaEnergy, v; 

  // majoranta
  G4double x0 = cut/kineticEnergy;
  G4double vmax = dataSB[Z].Value(x0, y)*1.02;
  //  G4double invbeta1 = 0;

  const G4double epeaklimit= 300*MeV; 
  const G4double elowlimit = 10*keV; 

  // majoranta corrected for e-
  bool isElectron = true;
  if(isElectron && x0 < 0.97 && 
     ((kineticEnergy > epeaklimit) || (kineticEnergy < elowlimit))) {
    //    G4double ylim = min(ylimit[Z],1.1*dataSB[Z].Value(0.97, y)); 
    //@@@ const G4double emaxlog = 4*log(10.);
    //@@@ ylimit[Z] = v->Value(0.97, emaxlog); 
    G4double ylim = fmin(dataSB[Z].Value(0.97, 4*log(10.)),
                        1.1*dataSB[Z].Value(0.97, y)); 
    if(ylim > vmax) { vmax = ylim; }
  }
  if(x0 < 0.05) { vmax *= 1.2; }

  do {
    //++ncount;
    //    G4double x = exp(xmin + G4UniformRand()*(xmax - xmin)) - densityCorr;
    G4double auxrand = GPUniformRand(fDevStates,fThreadId);
    G4double x = exp(xmin + auxrand*(xmax - xmin))-densityCorr;
    if(x < 0.0) { x = 0.0; }
    gammaEnergy = sqrt(x);
    G4double x1 = gammaEnergy/kineticEnergy;
    v = dataSB[Z].Value(x1, y);

#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
    std::cout<<"SB model: auxrand="<< auxrand
	     <<", gammaEne="<< gammaEnergy
	     <<", x1="<< x1 <<", y="<< y
	     <<", v="<< v
	     << std::endl;
#endif
#endif

    // correction for positrons        
    if(!isElectron) {
      G4double e1 = kineticEnergy - cut;
      G4double invbeta1 = (e1 + particleMass)/sqrt(e1*(e1 + 2*particleMass));
      G4double e2 = kineticEnergy - gammaEnergy;
      G4double invbeta2 = (e2 + particleMass)/sqrt(e2*(e2 + 2*particleMass));
      G4double xxx = twopi*fine_structure_const*currentZ*(invbeta1 - invbeta2); 

      //      if(xxx < expnumlim) { v = 0.0; } //expnumlim = -12.
      if(xxx < -12. ) { v = 0.0; }
      else { v *= exp(xxx); }
    }
   
    //    if (v > 1.05*vmax && nwarn < 20) {
    //      ++nwarn;
    //      if ( 20 == nwarn ) {
    //      }
    //    }
  } while (v < vmax*GPUniformRand(fDevStates, fThreadId));
  //  } while (v < vmax*G4UniformRand());

  //
  // angles of the emitted gamma. ( Z - axis along the parent particle)
  // use general interface
  //

  //  G4ThreeVector gammaDirection = 
  //    GetAngularDistribution()->SampleDirection(dp, totalEnergy-gammaEnergy,
  //                                              Z, couple->GetMaterial());
  GPThreeVector refDirection = 
    GPThreeVector_unit(GPThreeVector_create(track->px,track->py,track->pz));
  GPThreeVector gammaDirection = 
    DipBustGenerator_SampleDirection(kineticEnergy,refDirection);

  // create G4DynamicParticle object for the Gamma
  //  G4DynamicParticle* gamma = 
  //    new G4DynamicParticle(theGamma,gammaDirection,gammaEnergy);
  //  vdp->push_back(gamma);

  //@@@G4FWP fill secondary (one photon)
  FillSecondary(track,gammaDirection,gammaEnergy,0);

  //update electron GXTrack
  //  G4ThreeVector direction = (totMomentum*dp->GetMomentumDirection()
  //                             - gammaEnergy*gammaDirection).unit();
  GPThreeVector direction = GPThreeVector_unit(GPThreeVector_sub(
                            GPThreeVector_mult(refDirection,totMomentum),
                            GPThreeVector_mult(gammaDirection,gammaEnergy)));

  // energy of primary
  G4double finalE = kineticEnergy - gammaEnergy;

  // stop tracking and create new secondary instead of primary
  if(gammaEnergy > SecondaryThreshold()) {
    pParticleChange->ProposeTrackStatus(fStopAndKill);
    pParticleChange->SetProposedKineticEnergy(0.0);
    //@@@never happen: SecondaryThreshold = DBL_MAX
    //    G4DynamicParticle* el = 
    //      new G4DynamicParticle(const_cast<G4ParticleDefinition*>(particle),
    //			    direction, finalE);
    //    vdp->push_back(el);
    //continue tracking
  } else {
    pParticleChange->SetProposedMomentumDirection(direction);
    pParticleChange->SetProposedKineticEnergy(finalE);
  }
}

//---------------------------------------------------------------------------
//
// Combined Methods
//
//---------------------------------------------------------------------------

FQUALIFIER G4double 
GPSeltzerBergerRelModel::ComputeDXSectionPerAtom(G4double gammaEnergy) 
{
  G4double dxsection = 0.0;
  if(kinEnergy < 1.0*GeV) {
    dxsection = ComputeDXSectionPerAtom_SeltzerBergerModel(gammaEnergy);
  }
  else {
    dxsection = ComputeDXSectionPerAtom_eBremsstrahlungRelModel(gammaEnergy);
  }
  return dxsection;
}

FQUALIFIER void 
GPSeltzerBergerRelModel::SampleSecondaries(GXTrack* track,
					   GPMaterial* material,
					   G4double cutEnergy,
					   G4double maxEnergy) 
{
  if(track->E < energyLimitModels ) {
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
    std::cout<<"Calling SeltzerBergerModel::SampleSecondaries"<< std::endl;
#endif
#endif
    SampleSecondaries_SeltzerBergerModel(track,material,cutEnergy,maxEnergy);
  }
  else {
#ifdef GPUDEBUG
#ifndef __CUDA_ARCH__
    std::cout<<"GPSeltzerBerterRelModel: calling eBremsstrahlungRelModel::SampleSecondaries"<< std::endl;
#endif
#endif
     SampleSecondaries_eBremsstrahlungRelModel(track,material,cutEnergy,maxEnergy);
  }
}

//---------------------------------------------------------------------------
//
// G4VEmModel
//
//---------------------------------------------------------------------------

FQUALIFIER 
void GPSeltzerBergerRelModel::SetParticleChange(GPVParticleChange* p)
// 				   GPUniversalFluctuation* f)
{
  if(p && pParticleChange != p) { pParticleChange = p; }
  //  flucModel = f;
}

FQUALIFIER 
void GPSeltzerBergerRelModel::SetEnergyLimitModes(G4double energyLimit)
{
  energyLimitModels = energyLimit;  
}

FQUALIFIER
G4double GPSeltzerBergerRelModel::SecondaryThreshold()
{
  return secondaryThreshold;
}

FQUALIFIER
void GPSeltzerBergerRelModel::SetSecondaryThreshold(G4double val) 
{
  secondaryThreshold = val;
}

FQUALIFIER 
G4double GPSeltzerBergerRelModel::HighEnergyLimit()
{
  return highLimit;
}

FQUALIFIER 
G4double GPSeltzerBergerRelModel::LowEnergyLimit()
{
  return lowLimit;
}

FQUALIFIER 
G4double GPSeltzerBergerRelModel::HighEnergyActivationLimit()
{
  return eMaxActive;
}

FQUALIFIER 
G4double GPSeltzerBergerRelModel::LowEnergyActivationLimit()
{
  return eMinActive;
}

FQUALIFIER
void GPSeltzerBergerRelModel::SetHighEnergyLimit(G4double val)
{
  highLimit = val;
}

FQUALIFIER
void GPSeltzerBergerRelModel::SetLowEnergyLimit(G4double val)
{
  lowLimit = val;
}

FQUALIFIER 
void GPSeltzerBergerRelModel::SetActivationHighEnergyLimit(G4double val)
{
  eMaxActive = val;
}

FQUALIFIER 
void GPSeltzerBergerRelModel::SetActivationLowEnergyLimit(G4double val)
{
  eMinActive = val;
}

FQUALIFIER 
G4bool GPSeltzerBergerRelModel::IsActive(G4double kinEnergy)
{
  return (kinEnergy >= eMinActive && kinEnergy <= eMaxActive);
}

//---------------------------------------------------------------------------
//
// Other Methods
//
//---------------------------------------------------------------------------

//GetAngularDistribution()->SampleDirection 
//see the constructor of G4eBremsstrahlungRelModel where the angle model
//is set: SetAngularDistribution(new G4DipBustGenerator());
//
//G4ThreeVector& 
//G4DipBustGenerator::SampleDirection(const G4DynamicParticle* dp,
//                                    G4double, G4int, const G4Material*)
// change arguments with
// G4double eTkin = dp->GetKineticEnergy();
// refDirection = dp->GetMomentumDirection()
// 
FQUALIFIER GPThreeVector 
GPSeltzerBergerRelModel::DipBustGenerator_SampleDirection(G4double eTkin, 
					     GPThreeVector refDirection)
{
  G4double c, cosTheta, delta, cofA, signc = 1., a, power = 1./3.;


  //  c = 4. - 8.*G4UniformRand();
  c = 4. - 8.*GPUniformRand(fDevStates,fThreadId);
  a = c;
 
  if( c < 0. )
  {
    signc = -1.;
    a     = -c;
  }
  delta  = sqrt(a*a+4.);
  delta += a;
  delta *= 0.5; 

  cofA = -signc*pow(delta, power);

  cosTheta = cofA - 1./cofA;

  G4double tau = eTkin/electron_mass_c2;
  G4double beta = sqrt(tau*(tau + 2.))/(tau + 1.);

  cosTheta = (cosTheta + beta)/(1 + cosTheta*beta);

  G4double sinTheta = sqrt((1 - cosTheta)*(1 + cosTheta));
  //  G4double phi  = twopi*G4UniformRand(); 
  G4double phi  = twopi*GPUniformRand(fDevStates,fThreadId); 

  GPThreeVector fLocalDirection = GPThreeVector_create(sinTheta*cos(phi), 
						       sinTheta*sin(phi),
						       cosTheta);
  //  fLocalDirection.rotateUz(dp->GetMomentumDirection());
  GPThreeVector_rotateUz(&fLocalDirection,refDirection);

  return fLocalDirection;

}

GXTrack& GPSeltzerBergerRelModel::GetSecondary()
{
  return theGamma;
}


FQUALIFIER 
void GPSeltzerBergerRelModel::FillSecondary(GXTrack* track, 
					    GPThreeVector direction,
					    G4double energy, 
					    G4double charge)
{
  theGamma.x  = track->x;
  theGamma.y  = track->y;
  theGamma.z  = track->z;
  theGamma.s  = 0;
  theGamma.px = energy*direction.x;
  theGamma.py = energy*direction.y;
  theGamma.pz = energy*direction.z;
  theGamma.q  = charge;
}
