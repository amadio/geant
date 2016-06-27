

#include "Spline.h"

// cstdlib is for std::abs(int)
#include <cstdlib>

namespace geant {

// spline does not own fXdata and fYdata
Spline::Spline(double *xdata, double *ydata, int numdata, bool isconstrained) {
  fNPoints = numdata;
  fXdata   = xdata;
  fYdata   = ydata;
  fA.resize(fNPoints);
  fB.resize(fNPoints);
  fC.resize(fNPoints);
  fD.resize(fNPoints);
  fIsConstrained = isconstrained;
  if (fIsConstrained) {
    SetParametersConstrained(); // constrained natural cubic spline
  } else {
    SetParameters();            // natural cubic spline
  }
  // set some members used for the B-search
  if (*(fXdata)>*(fXdata+(fNPoints-1))) {
    fDirection = 1;          // decreasing x grid
    fLowerm    = fNPoints-1;
    fUpperm    = -1;
  } else {
    fDirection = 0;         // increasing x grid
    fLowerm    = -1;
    fUpperm    = fNPoints-1;
  }
}

// spline does not own fXdata and fYdata
void Spline::SetUpSpline(double *xdata, double *ydata, int numdata, bool isconstrained) {
  fNPoints = numdata;
  fXdata   = xdata;
  fYdata   = ydata;
  fA.resize(fNPoints);
  fB.resize(fNPoints);
  fC.resize(fNPoints);
  fD.resize(fNPoints);
  fIsConstrained = isconstrained;
  if (fIsConstrained) {
    SetParametersConstrained(); // constrained natural cubic spline
  } else {
    SetParameters();            // natural cubic spline
  }
  // set some members used for the B-search
  if (*(fXdata)>*(fXdata+(fNPoints-1))) {
    fDirection = 1;          // decreasing x grid
    fLowerm    = fNPoints-1;
    fUpperm    = -1;
  } else {
    fDirection = 0;         // increasing x grid
    fLowerm    = -1;
    fUpperm    = fNPoints-1;
  }
}

// get interpollated Y value at the given X=val point
double Spline::GetValueAt(double val) {
  int    m,ml,mu,mav;
  // check if 'val' is above/below the highes/lowest value
  if (val>=*(fXdata+(fUpperm+fDirection))) {
    m = fUpperm+2*fDirection-1;
  } else if (val<=*(fXdata+(fLowerm+1-fDirection))) {
    m = fLowerm-2*fDirection+1;
  } else {  // Perform a binary search to find the interval val is in
    ml = fLowerm;
    mu = fUpperm;
    while (std::abs(mu-ml)>1) {
      mav = (ml+mu)/2.;
      if (val<*(fXdata+mav))
        mu = mav;
      else
        ml = mav;
      }
    m = mu+fDirection-1;
  }
  if (!fIsConstrained)
    val -= (*(fXdata+m));

 return fA[m]+val*(fB[m]+val*(fC[m]+val*fD[m]));
}


// we need an other GetSplie such that if:
// - we know that x-values are increasing
// - and if we know which is the lower bin index j such that x_j < x <x_j+1
// - NOTE: before calling this, one needs to make sure that x_min < x < x_max
double Spline::GetValueAt(double x, int indx){
  if (!fIsConstrained)
    x -= (*(fXdata+indx));
 return fA[indx] + x*(fB[indx] + x*( fC[indx] + x*fD[indx]));
}


// init the spline parameters
// compute the natural cubic spline paraneters a,b,c,d
void Spline::SetParameters(){
  int    n,m1,m2,m,mr;
  double s,r;
  m1 = 1;
  m2 = fNPoints-1;
  s  = 0.0;
  for (m=0; m < m2; ++m) {
    fD[m] = (*(fXdata+(m+1))) - (*(fXdata+m));
    r     = ((*(fYdata+(m+1))) - (*(fYdata+m)))/fD[m];
    fC[m] = r - s;
    s     = r;
  }
  s     = 0.0;
  r     = 0.0;
  fC[0] = 0.0;
  fC[fNPoints-1] = 0.0;
  for (m=m1; m<m2; ++m) {
    fC[m] = fC[m] + r*fC[m-1];
    fB[m] = 2.0*(  (*(fXdata+(m-1))) - (*(fXdata+(m+1)))  ) - r*s;
    s     = fD[m];
    r     = s/fB[m];
  }
  mr = m2-1;
  for (m=m1; m<m2; ++m) {
    fC[mr] = (fD[mr]*fC[mr+1] - fC[mr])/fB[mr];
    mr     = mr - 1;
  }
  for (m=0; m<m2; ++m) {
   s     = fD[m];
   r     = fC[m+1] - fC[m];
   fD[m] = r/s;
   fC[m] = 3.0*fC[m];
   fB[m] = (  *(fYdata+(m+1))  -  *(fYdata+m)  )/s - (fC[m]+r)*s;
   fA[m] = *(fYdata+m);
  }
}

void Spline::SetParametersConstrained(){
  double f1pxi    = 0.0;
  double f1pxim1  = 0.0;
  double f1ppxi   = 0.0;
  double f1ppxim1 = 0.0;
  for (int i=1;i<fNPoints;++i) {
    f1pxi = 0.0;
    if (i<fNPoints-1) {
      double dum0  = (*(fYdata +(i+1)))-(*(fYdata+i));
      double dum1  = (*(fYdata+i))-(*(fYdata+i-1));
      if (dum0*dum1>0.0) {
        f1pxi = 2.0/( ( (*(fXdata+i+1))-(*(fXdata+i)) )/dum0 + (( (*(fXdata+i))-(*(fXdata+i-1)))/dum1) );
      }
    } else { // last point
      f1pxi = 1.5*( (*(fYdata-fNPoints-1))-(*(fYdata-fNPoints-2)) )/( (*(fXdata+fNPoints-1))-(*(fXdata+fNPoints-2))) - 0.5*f1pxim1;
    }
    if (i==1) { // first point
      f1pxim1 = 1.5*( (*(fYdata+1))-(*(fYdata)) )/( (*(fXdata+1))-(*(fXdata)) ) - 0.5*f1pxi;
    }
    double xdel1 = (*(fXdata+i))-(*(fXdata+i-1));
    double ydel1 = (*(fYdata+i))-(*(fYdata+i-1));
    f1ppxim1     = -2.0*(f1pxi+2.0*f1pxim1)/xdel1 + 6.0*ydel1/(xdel1*xdel1);
    f1ppxi       =  2.0*(2.0*f1pxi+f1pxim1)/xdel1 - 6.0*ydel1/(xdel1*xdel1);
    double xi2   = (*(fXdata+i)) * (*(fXdata+i));
    double xim12 = (*(fXdata+i-1)) * (*(fXdata+i-1));
    fD[i-1] = (f1ppxi-f1ppxim1)/(6.0*xdel1);
    fC[i-1] = 0.5*( (*(fXdata+i))*f1ppxim1-(*(fXdata+i-1))*f1ppxi )/xdel1;
    fB[i-1] = (ydel1-fC[i-1]*(xi2-xim12) - fD[i-1]*( (*(fXdata+i))*xi2 - (*(fXdata+i-1))*xim12 ))/xdel1;
    fA[i-1] = (*(fYdata+i-1)) - fB[i-1]*(*(fXdata+i-1)) - fC[i-1]*xim12 - fD[i-1]*xim12*(*(fXdata+i-1));
    f1pxim1 = f1pxi;
  }
}

} // namespace geant
