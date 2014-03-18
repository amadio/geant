#include "GPPhysicsTable.h"
#include <stdio.h>

FQUALIFIER GPPhysicsVector::GPPhysicsVector() {
	type = 2;
	edgeMin = 0;
	edgeMax = 0;
	numberOfNodes = 0;
	useSpline = false;
	dBin = 0;
	baseBin = 0;
	isSecondDerivativeFilled = false;
	for (int i = 0; i < maxPVDataVector; i++) {
		dataVector[i] = 0;
		binVector[i] = 0;
		secDerivative[i] = 0;
	}
}

FQUALIFIER GPPhysicsVector::GPPhysicsVector(GPPhysicsVector& right) {
	useSpline = right.useSpline;
	isSecondDerivativeFilled = right.isSecondDerivativeFilled;
	type = right.type;
	numberOfNodes = right.numberOfNodes;
	edgeMin = right.edgeMin;
	edgeMax = right.edgeMax;
	dBin = right.dBin;
	baseBin = right.baseBin;
	for (int i = 0; i < maxPVDataVector; i++) {
		dataVector[i] = right.dataVector[i];
		binVector[i] = right.binVector[i];
		secDerivative[i] = right.secDerivative[i];
	}
}

FQUALIFIER
bool GPPhysicsVector::SplinePossible() {
	// Initialise second derivative array. If neighbor energy coincide 
	// or not ordered than spline cannot be applied

	for (int i = 0; i < maxPVDataVector; i++)
		secDerivative[i] = 0.0;
	if (!useSpline) {
		return useSpline;
	}
	//  secDerivative.reserve(numberOfNodes);
	for (int j = 0; j < numberOfNodes; ++j) {
		secDerivative[j] = 0.0;
		if (j > 0) {
			if (binVector[j] - binVector[j - 1] <= 0.) {
				useSpline = false;
			}
		}
	}
	return useSpline;
}

FQUALIFIER
void GPPhysicsVector::ComputeSecDerivatives() {
	//  A simplified method of computation of second derivatives 

	if (!SplinePossible()) {
		return;
	}

	if (3 > numberOfNodes) // cannot compute derivatives for less than 4 bins
			{
		useSpline = false;
		return;
	}

	int n = numberOfNodes - 1;

	for (int i = 1; i < n; ++i) {
		secDerivative[i] = 3.0
				* ((dataVector[i + 1] - dataVector[i])
						/ (binVector[i + 1] - binVector[i])
						- (dataVector[i] - dataVector[i - 1])
								/ (binVector[i] - binVector[i - 1]))
				/ (binVector[i + 1] - binVector[i - 1]);
	}
	secDerivative[n] = secDerivative[n - 1];
	secDerivative[0] = secDerivative[1];
}

FQUALIFIER
void GPPhysicsVector::ComputeSecondDerivatives(G4double firstPointDerivative,
		G4double endPointDerivative) {
	//  A standard method of computation of second derivatives 
	//  First derivatives at the first and the last point should be provided
	//  See for example W.H. Press et al. "Numerical reciptes and C"
	//  Cambridge University Press, 1997.

	if (4 > numberOfNodes) // cannot compute derivatives for less than 4 bins
			{
		ComputeSecDerivatives();
		return;
	}

	if (!SplinePossible()) {
		return;
	}

	int n = numberOfNodes - 1;

	//  G4double* u = new double [n];
	G4double u[maxPVDataVector];

	G4double p, sig, un;

	u[0] = (6.0 / (binVector[1] - binVector[0]))
			* ((dataVector[1] - dataVector[0]) / (binVector[1] - binVector[0])
					- firstPointDerivative);

	secDerivative[0] = -0.5;

	// Decomposition loop for tridiagonal algorithm. secDerivative[i]
	// and u[i] are used for temporary storage of the decomposed factors.

	for (int i = 1; i < n; ++i) {
		sig = (binVector[i] - binVector[i - 1])
				/ (binVector[i + 1] - binVector[i - 1]);
		p = sig * (secDerivative[i - 1]) + 2.0;
		secDerivative[i] = (sig - 1.0) / p;
		u[i] = (dataVector[i + 1] - dataVector[i])
				/ (binVector[i + 1] - binVector[i])
				- (dataVector[i] - dataVector[i - 1])
						/ (binVector[i] - binVector[i - 1]);
		u[i] = 6.0 * u[i] / (binVector[i + 1] - binVector[i - 1])
				- sig * u[i - 1] / p;
	}

	sig = (binVector[n - 1] - binVector[n - 2])
			/ (binVector[n] - binVector[n - 2]);
	p = sig * secDerivative[n - 2] + 2.0;
	un = (6.0 / (binVector[n] - binVector[n - 1]))
			* (endPointDerivative
					- (dataVector[n] - dataVector[n - 1])
							/ (binVector[n] - binVector[n - 1])) - u[n - 1] / p;
	secDerivative[n] = un / (secDerivative[n - 1] + 2.0);

	// The back-substitution loop for the triagonal algorithm of solving
	// a linear system of equations.

	for (int k = n - 1; k > 0; --k) {
		secDerivative[k] *= (secDerivative[k + 1]
				- u[k] * (binVector[k + 1] - binVector[k - 1])
						/ (binVector[k + 1] - binVector[k]));
	}
	secDerivative[0] = 0.5 * (u[0] - secDerivative[1]);

}

FQUALIFIER
void GPPhysicsVector::FillSecondDerivatives() {
	// Computation of second derivatives using "Not-a-knot" endpoint conditions
	// B.I. Kvasov "Methods of shape-preserving spline approximation"
	// World Scientific, 2000

	if (5 > numberOfNodes) // cannot compute derivatives for less than 4 points
			{
		ComputeSecDerivatives();
		return;
	}

	if (!SplinePossible()) {
		return;
	}

	int n = numberOfNodes - 1;

	//cout << "GPPhysicsVector_FillSecondDerivatives() n= " << n 
	//     << "   " << this << endl;
	// cout << *this << endl;

	//  G4double* u = new double [n];
	G4double u[maxPVDataVector];

	G4double p, sig;

	u[1] = ((dataVector[2] - dataVector[1]) / (binVector[2] - binVector[1])
			- (dataVector[1] - dataVector[0]) / (binVector[1] - binVector[0]));
	u[1] = 6.0 * u[1] * (binVector[2] - binVector[1])
			/ ((binVector[2] - binVector[0]) * (binVector[2] - binVector[0]));

	// Decomposition loop for tridiagonal algorithm. secDerivative[i]
	// and u[i] are used for temporary storage of the decomposed factors.

	secDerivative[1] = (2.0 * binVector[1] - binVector[0] - binVector[2])
			/ (2.0 * binVector[2] - binVector[0] - binVector[1]);

	for (int i = 2; i < n - 1; ++i) {
		sig = (binVector[i] - binVector[i - 1])
				/ (binVector[i + 1] - binVector[i - 1]);
		p = sig * secDerivative[i - 1] + 2.0;
		secDerivative[i] = (sig - 1.0) / p;
		u[i] = (dataVector[i + 1] - dataVector[i])
				/ (binVector[i + 1] - binVector[i])
				- (dataVector[i] - dataVector[i - 1])
						/ (binVector[i] - binVector[i - 1]);
		u[i] = (6.0 * u[i] / (binVector[i + 1] - binVector[i - 1]))
				- sig * u[i - 1] / p;
	}

	sig = (binVector[n - 1] - binVector[n - 2])
			/ (binVector[n] - binVector[n - 2]);
	p = sig * secDerivative[n - 3] + 2.0;
	u[n - 1] = (dataVector[n] - dataVector[n - 1])
			/ (binVector[n] - binVector[n - 1])
			- (dataVector[n - 1] - dataVector[n - 2])
					/ (binVector[n - 1] - binVector[n - 2]);
	u[n - 1] = 6.0 * sig * u[n - 1] / (binVector[n] - binVector[n - 2])
			- (2.0 * sig - 1.0) * u[n - 2] / p;

	p = (1.0 + sig) + (2.0 * sig - 1.0) * secDerivative[n - 2];
	secDerivative[n - 1] = u[n - 1] / p;

	// The back-substitution loop for the triagonal algorithm of solving
	// a linear system of equations.

	for (int k = n - 2; k > 1; --k) {
		secDerivative[k] *= (secDerivative[k + 1]
				- u[k] * (binVector[k + 1] - binVector[k - 1])
						/ (binVector[k + 1] - binVector[k]));
	}
	secDerivative[n] = (secDerivative[n - 1]
			- (1.0 - sig) * secDerivative[n - 2]) / sig;
	sig = 1.0 - ((binVector[2] - binVector[1]) / (binVector[2] - binVector[0]));
	secDerivative[1] *= (secDerivative[2] - u[1] / (1.0 - sig));
	secDerivative[0] = (secDerivative[1] - sig * secDerivative[2])
			/ (1.0 - sig);

	isSecondDerivativeFilled = true;
}

FQUALIFIER
int GPPhysicsVector::FindBinLocation(G4double energy) {
	return int(log10(energy) / dBin - baseBin);
}

FQUALIFIER
G4double GPPhysicsVector::SplineInterpolation(double energy, int lastBin) {
	// Spline interpolation is used to get the value. If the give energy
	// is in the highest bin, no interpolation will be Done. Because 
	// there is an extra bin hidden from a user at locBin=numberOfBin, 
	// the following interpolation is valid even the current locBin=
	// numberOfBin-1. 

	//  if(0 == secDerivative.size() ) { FillSecondDerivatives(); }
	if (!isSecondDerivativeFilled) {
		FillSecondDerivatives();
	}

	// check bin value
	G4double x1 = binVector[lastBin];
	G4double x2 = binVector[lastBin + 1];
	G4double delta = x2 - x1;

	G4double a = (x2 - energy) / delta;
	G4double b = (energy - x1) / delta;

	// Final evaluation of cubic spline polynomial for return   
	G4double y1 = dataVector[lastBin];
	G4double y2 = dataVector[lastBin + 1];

	G4double res = a * y1 + b * y2
			+ ((a * a * a - a) * secDerivative[lastBin]
					+ (b * b * b - b) * secDerivative[lastBin + 1]) * delta
					* delta / 6.0;

	return res;
}

FQUALIFIER
G4double GPPhysicsVector::LinearInterpolation(double energy, int lastBin) {
	// Linear interpolation is used to get the value. If the give energy
	// is in the highest bin, no interpolation will be Done. Because 
	// there is an extra bin hidden from a user at locBin=numberOfBin, 
	// the following interpolation is valid even the current locBin=
	// numberOfBin-1. 

	G4double intplFactor = (energy - binVector[lastBin])
			/ (binVector[lastBin + 1] - binVector[lastBin]); // Interpol. factor

	return dataVector[lastBin]
			+ (dataVector[lastBin + 1] - dataVector[lastBin]) * intplFactor;
}

FQUALIFIER
G4double GPPhysicsVector::Interpolation(double energy, int lastBin) {
	G4double value = 0;
	if (useSpline) {
		value = SplineInterpolation(energy, lastBin);
	} else {
		value = LinearInterpolation(energy, lastBin);
	}
	return value;
}

FQUALIFIER
void GPPhysicsVector::SetSpline(bool val) {
	useSpline = val;
	if (!isSecondDerivativeFilled) {
		FillSecondDerivatives();
		isSecondDerivativeFilled = true;
	}
}

FQUALIFIER
G4double GPPhysicsVector::Value(double energy) {

	G4double value = 0.0;

	if (energy <= edgeMin) {
		value = dataVector[0];
	} else if (energy >= edgeMax) {
		value = dataVector[numberOfNodes - 1];
	} else {
		int bin = FindBinLocation(energy);
		value = Interpolation(energy, bin);
	}
	return value;
}

FQUALIFIER GPPhysicsTable::GPPhysicsTable() {
	nPhysicsVector = 0;
}

FQUALIFIER
G4double GPPhysicsVector::Energy(size_t binNumber)
{
  return binVector[binNumber];
}

FQUALIFIER
G4double GPPhysicsVector::InvRangeValue(double energy) 
{
  G4double value = 0.0;
  
  if (energy <= edgeMin) {
    value = dataVector[0];
  } else if (energy >= edgeMax) {
    value = dataVector[numberOfNodes - 1];
  } else {
    int bin = InvRangeFindBinLocation(energy);
    value = Interpolation(energy, bin);
  }
  return value;
}

FQUALIFIER
int GPPhysicsVector::InvRangeFindBinLocation(G4double energy) 
{
  G4int n1 = 0;
  G4int n2 = numberOfNodes/2;
  G4int n3 = numberOfNodes - 1;
  while (n1 != n3 - 1) {
    if (energy > binVector[n2]) { 
      n1 = n2; 
    }
    else { 
      n3 = n2; 
    }
    n2 = n1 + (n3 - n1 + 1)/2;
  }
  
  return (int)n1;
}

FQUALIFIER
void GPPhysicsTable::Print() {
	/*
	 printf("%f\n",nPhysicsVector);

	 for(int idx=0; idx<nPhysicsVector; idx++){
	 GPPhysicsVector* pv = &(physicsVectors[idx]);
	 printf("%d\n",pv->type);
	 printf("%f %f %d\n",pv->edgeMin,pv->edgeMax,pv->numberOfNodes);
	 printf("%d\n",pv->numberOfNodes);
	 for(int j=0; j<pv->numberOfNodes; j++) {
	 printf("%f %e\n",pv->binVector[j],pv->dataVector[j]);
	 }
	 }
	 */
}
