#ifndef GXDiscreteSampling_H
#define GXDiscreteSampling_H 1

//----------------------------------------------------------------------------
// Build tables to sample secondary particles of EM Physics processes.  These 
// utility functions are to be called at initialization and not optimized for 
// vectorization. There may be also the more efficient algorithms to build the
// inverse pdf and the alias table -  2014/06 S.Y. Jun
//----------------------------------------------------------------------------

void BuildInversePDF(const G4int n,
		     G4double   *p,
		     G4double   *q)
{
  // Build inverse cumulative distribution 
  //
  // input  :    n  (dimension of discrete outcomes + 1)
  //           p[n] (probability distribution)
  // output :  q[n] (inverse probability distribution) 

  //build cumulative pdf 
  G4double *c = (G4double*) malloc(n*sizeof(G4double)); 

  c[0] = p[0];
  for(int i = 1; i < n ; ++i) {
    c[i] = c[i-1] + p[i]; 
  }

  //temporary array
  G4int *ip = (G4int*) malloc(n*sizeof(G4int)); 

  //likelihood per equal probable event
  const G4double dy = 1.0/(n-1);  

  G4double test = 0.0; 
  for(int i = 0; i < n ; ++i) {
    test = dy*i; 
      for(int j = 0; j < n ; ++j) {
      if(test < c[j]) {
	ip[i] = j;
	break;
      }
    }
  }
  //for the last bin
  ip[n-1] = n-1;
  
  //linear interpolation within the interval with the same index
  unsigned int last, cnt, lcnt;
  last = cnt = lcnt = 0;
  
  for(int i = 0; i < n ; ++i) {
    if(ip[i] > last) {
      for(int j = 0 ; j < (cnt-lcnt) ; ++j ) {
	q[j+lcnt] = dy*(last + (1.0*j/(cnt-lcnt)));
      }
      last = ip[i];
      lcnt = cnt;
    }
    ++cnt;
  }
  //for the last bin
  q[n-1] = 1.0;

  free(c);
  free(ip);
}

void BuildAliasTables(const G4int n,
		      G4double   *p,
		      G4int      *a,
		      G4double   *q)
{
  // Build alias and alias probability
  //    
  // Reference: (1) A.J. Walker, "An Efficient Method for Generating Discrete 
  // Random Variables with General Distributions" ACM Trans. Math. Software, 3,
  // 3, 253-256 (1977) (2) A.L. Edwards, J.A. Rathkopf, and R.K. Smidt, 
  // "Extending the Alias Monte Carlo Sampling Method to General Distributions"
  // UCRL-JC-104791 (1991)
  //
  // input :    n  (dimension of discrete outcomes)
  //          p[n] (probability density function)
  // output:  a[n] (alias)
  //          q[n] (non-alias probability) 

  //temporary array
  G4double *ap = (G4double*)malloc(n*sizeof(G4double)); 

  //likelihood per equal probable event
  const G4double cp = 1.0/(n-1);

  //initialize
  for(int i = 0; i < n ; ++i) {
    a[i] = -1;
    ap[i] = p[i];
  }
  
  //O(n) iterations
  G4int iter = n;
  
  do {
    int donor = 0;
    int recip = 0;
    
    //have the better search algorithm?
    for(int j = donor; j < n ; ++j) {
      if(ap[j] >= cp) {
	donor = j;
	break;
      }
    }
    
    for(int j = recip; j < n ; ++j) {
      if(ap[j] > 0.0 && ap[j] < cp) {
	recip = j;
	break;
      }
    }
    
    //alias and non-alias probability
    a[recip] = donor;
    q[recip] = n*ap[recip];
    
    //update pdf 
    ap[donor] = ap[donor] - (cp-ap[recip]);
    ap[recip] = 0.0;

    --iter;
  }
  while (iter > 0);

  free(ap);
}

G4double KleinNishinaCrossSection(G4double energy0, 
				  G4double energy1) 
{
  // based on Geant4 : G4KleinNishinaCompton
  // input  : energy0 (incomming photon energy)
  //          energy1 (scattered photon energy)
  // output : dsigma  (differential cross section) 

  G4double E0_m = energy0/electron_mass_c2 ;
  G4double epsilon = energy1/energy0;

  G4double onecost = (1.- epsilon)/(epsilon*E0_m);
  G4double sint2   = onecost*(2.-onecost);
  G4double greject = 1. - epsilon*sint2/(1.+ epsilon*epsilon);
  G4double dsigma = (epsilon + 1./epsilon)*greject;

  return dsigma;
}

double MollerBhabhaCrossSection(G4double kineticEnergy, 
				G4double deltaRayEnergy,
				bool isElectron) 
{
  // based on Geant3 : Simulation of the delta-ray production (PHY331-1)
  // input  : energy0 (incomming photon energy)
  //          energy1 (scattered photon energy)
  // output : dsigma  (differential cross section) 

  G4double dcross = 0.0;

  G4double tau = kineticEnergy/electron_mass_c2;
  G4double gam = tau + 1.0;
  G4double gamma2 = gam*gam;
  G4double epsil = deltaRayEnergy/kineticEnergy;

  if(isElectron) {
    //Moller (e-e-) scattering
    G4double fgam = (2*gam-1.0)/gamma2;
    G4double x = 1/epsil;
    G4double y = 1/(1.0-epsil);

    dcross = (1-fgam + x*(x-fgam) + y*y*fgam);
  }
  else {
    //Bhabha (e+e-) scattering
    G4double beta2 = tau*(tau + 2)/gamma2;

    G4double y   = 1.0/(1.0 + gam);
    G4double y2  = y*y;
    G4double y12 = 1.0 - 2.0*y;
    G4double y122= y12*y12;
    
    G4double b1  = 2.0 - y2;
    G4double b2  = y12*(3.0 + y2);
    G4double b4  = y122*y12;
    G4double b3  = b4 + y122;
    
    dcross = 1./(beta2*epsil*epsil) - b1/epsil + b2 - b3*epsil + b4*epsil*epsil;
  }

  return dcross;
}

G4double SeltzerBergerCrossSection(GPPhysics2DVector data,
				   G4double z, 
				   G4double w) 
{
  //based on Geant4 : G4KleinNishinaCompton
  // input  : data (SeltzerBerger data)
  //          z (ratio of photon energy to electron energy)
  //          w (log of the incident electron energy)
  // output : dsigma  (differential cross section) 
  G4double dsigma = data.Value(z,w);

  return dsigma;
}

void BuildPDF_KleinNishina(const G4int n, 
			   G4double    x, 
			   G4double   *p)
{
  // Build PDF based on KleinNishina Formular
  //
  // input  :   n  (dimension)
  //            x  (incident photon energy)
  // output : p[n] (probability distribution) 

  G4double ymin = x/(1+2.0*x/electron_mass_c2);
  G4double dy = (x - ymin)/(n-1);
  G4double yo = ymin + 0.5*dy;
  
  G4double sum = 0.;
  for(int i = 0; i < n ; ++i) {
    G4double y = yo + dy*i;
    G4double xsec = KleinNishinaCrossSection(x,y);
    p[i] = xsec;
    sum += xsec;
  }

  //normalization
  sum = 1.0/sum;
  for(int i = 0; i < n ; ++i) {
    p[i] *= sum;
  }
}

void BuildPDF_MollerBhabha(const G4int n, 
			   G4double    x, 
			   bool        b, 
			   G4double   *p)
{
  // Build PDF based on MollerBhabha Formular
  //
  // input  :   n  (dimension)
  //            x  (incident e-/e+ energy)
  //            b  (electron flag)
  // output : p[n] (probability distribution) 

  //cutoffs for the delta ray production
  G4double ymin = 1.0; 
  G4double ymax = 0.5*x;
  if(!b)   ymax += ymax;

  G4double dy = (ymax - ymin)/(n-1);
  G4double yo = ymin + 0.5*dy;
  
  G4double sum = 0.;
  for(int i = 0; i < n ; ++i) {
    G4double y = yo + dy*i;
    G4double xsec = MollerBhabhaCrossSection(x,y,b);
    p[i] = xsec;
    sum += xsec;
  }

  //normalization
  sum = 1.0/sum;
  for(int i = 0; i < n ; ++i) {
    p[i] *= sum;
  }
}

void BuildPDF_SeltzerBerger(const G4int n, 
                            G4int Z,
			    GPPhysics2DVector* data, 
			    G4double  df, 
			    G4double  xmin, 
			    G4double  xmax, 
			    G4double  x, 
			    G4double *p)
{
  // Build PDF based on SeltzerBerger parameterization (G4SeltzerBergerModel)
  //
  // input  :   n     (dimension)
  //            Z     (atomic number)
  //           *data  (SeltzerBerger cross section data)
  //            df    (electron density * Migdal constant)
  //            xmin  (miminum energy)
  //            xmax  (maxinum energy)
  //            x     (incident electron kinetic energy) in [MeV]
  //
  // output : p[n]    (probability distribution) 

  G4double emin = (xmin < x) ? xmin : x;
  G4double emax = (xmax < x) ? xmax : x;

  //total energy = x + electron_mass_c2 in MeV;
  G4double t = x + 0.510998910;

  //density correction 
  G4double dc = df*t*t;
  
  G4double ymin = log(emin*emin + dc);
  G4double ymax = log(emax*emax + dc);

  G4double dy = (ymax - ymin)/(n-1);
  G4double yo = ymin + 0.5*dy;

  G4double v = 1.0/x;
  G4double w = log(x);
  G4double sum = 0.0;

  for(int i = 0; i < n ; ++i) {
    G4double y = exp(yo+dy*i) - dc;    
    G4double z = (y < 0 ) ? 0 : sqrt(y)*v;

    //cross section based on the Seltzer-Berger Parameterization
    G4double xsec = data[Z].Value(z,w);
    p[i] = xsec;
    sum += xsec;    
  }

  //normalization
  sum = 1.0/sum;
  for(int i = 0; i < n ; ++i) {
    p[i] *= sum;
  }
}

void BuildInversePDF_KleinNishinaModel(G4int Z, 
				       const G4double xmin, 
				       const G4double xmax,
				       const G4int nrow,
				       const G4int ncol,
				       G4double **p,
				       G4double **q)
{
  // Build the 2-dimensional probability density function (p) in [xmin,xmax] 
  // and the inverse pdf (q) with the equal probable bin in [0,1]
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //           q[nrow][ncol] (inverse cumulative distribution) 

  //build pdf  
  G4double dx = (xmax - xmin)/nrow;
  G4double xo =  xmin + 0.5*dx;

  for(int i = 0; i < nrow ; ++i) {
    G4double x = xo + dx*i;
    BuildPDF_KleinNishina(ncol,x,p[i]);
  }

  //invert the cdf to sample the x for the given random draw(0,1] in F^{-1}(y)
  for(int i = 0; i < nrow ; ++i) {
    BuildInversePDF(ncol,p[i],q[i]); 
    
  }
}

void BuildInversePDF_MollerBhabha(G4int Z, 
				  const bool flag,
				  const G4double xmin, 
				  const G4double xmax,
				  const G4int nrow,
				  const G4int ncol,
				  G4double **p,
				  G4double **q)
{
  // Build the 2-dimensional probability density function (p) in [xmin,xmax] 
  // and the inverse pdf (q) with the equal probable bin in [0,1]
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  //           flag (electron flag)
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //           q[nrow][ncol] (inverse cumulative distribution) 

  //build pdf  
  G4double dx = (xmax - xmin)/nrow;
  G4double xo =  xmin + 0.5*dx;

  for(int i = 0; i < nrow ; ++i) {
    G4double x = xo + dx*i;
    BuildPDF_MollerBhabha(ncol,x,flag,p[i]);
  }

  //invert the cdf to sample the x for the given random draw(0,1] in F^{-1}(y)
  for(int i = 0; i < nrow ; ++i) {
    BuildInversePDF(ncol,p[i],q[i]); 
    
  }
}

void BuildInversePDF_SeltzerBerger(G4int Z, 
				   GPPhysics2DVector* data, 
				   G4double df,
				   const G4double xmin, 
				   const G4double xmax,
				   const G4int nrow,
				   const G4int ncol,
				   G4double **p,
				   G4double **q)
{
  // Build the two-dimensional probability density function (pdf) and the 
  // inverse cumulative distribution (inv) in [log(xmin),log(xmax)]
  //
  // input  :  Z    (atomic number)
  //          *data (SeltzerBerger cross section data)
  //           df   (electron density * Migdal constant)
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //           q[nrow][ncol] (inverse cumulative distribution) 

  G4double dx = (log(xmax) - log(xmin))/nrow;
  G4double xo = log(xmin) + 0.5*dx;

  for(int i = 0; i < nrow ; ++i) {
    G4double x = exp(xo+dx*i);
    BuildPDF_SeltzerBerger(ncol,Z,data,df,xmin,xmax,x,p[i]);
  }

  //inverse pdf to sample the x for the given random draw(0,1] in F^{-1}(y)
  for(int i = 0; i < nrow ; ++i) {
    BuildInversePDF(ncol,p[i],q[i]); 
  }
}

void BuildAliasSamplingTable_KleinNishina(G4int Z, 
					  const G4double xmin, 
					  const G4double xmax,
					  const G4int nrow,
					  const G4int ncol,
					  G4double **p,
					  G4int    **a,
					  G4double **q)
{
  // Build the 2-dimensional probability density function (KleinNishina pdf) 
  // in [xmin,xmax] and the alias table (a) and probability (q) with (ncol-1) 
  // equal probable events, each with likelihood 1/(ncol-1)
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //           a[nrow][ncol] (alias table) 
  //           q[nrow][ncol] (non-alias probability) 

  //build pdf  
  G4double dx = (xmax - xmin)/nrow;
  G4double xo =  xmin + 0.5*dx;

  for(int i = 0; i < nrow ; ++i) {
    G4double x = xo + dx*i;
    BuildPDF_KleinNishina(ncol,x,p[i]);
  }
  
  //build alias tables (alias and probability)
  for(int i = 0; i < nrow ; ++i) {
    BuildAliasTables(ncol,p[i],a[i],q[i]);
  }
}

void BuildAliasSamplingTable_MollerBhabha(G4int Z, 
					  const bool flag, 
					  const G4double xmin, 
					  const G4double xmax,
					  const G4int nrow,
					  const G4int ncol,
					  G4double **p,
					  G4int    **a,
					  G4double **q)
{
  // Build the 2-dimensional probability density function (MollerBhabha pdf) 
  // in [xmin,xmax] and the alias table (a) and probability (q) with (ncol-1) 
  // equal probable events, each with likelihood 1/(ncol-1)
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  //           flag (electron flag)
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //           a[nrow][ncol] (alias table) 
  //           q[nrow][ncol] (non-alias probability) 

  //build pdf  
  G4double dx = (xmax - xmin)/nrow;
  G4double xo =  xmin + 0.5*dx;

  for(int i = 0; i < nrow ; ++i) {
    G4double x = xo + dx*i;
    BuildPDF_MollerBhabha(ncol,x,flag,p[i]);
  }
  
  //build alias tables (alias and probability)
  for(int i = 0; i < nrow ; ++i) {
    BuildAliasTables(ncol,p[i],a[i],q[i]);
  }
}

void BuildAliasSamplingTable_SeltzerBerger(G4int Z, 
					   GPPhysics2DVector* data, 
					   G4double df,
					   const G4double xmin, 
					   const G4double xmax,
					   const G4int nrow,
					   const G4int ncol,
					   G4double **p,
					   G4int    **a,
					   G4double **q)
{
  // Build the two-dimensional probability density function (Seltzer-Berger) 
  // in [log(xmin),log(xmax)] and the alias table (a) and probability (q) with 
  // (ncol-1) equal probable events, each with likelihood 1/(ncol-1)
  //
  // input  :  Z    (atomic number)
  //          *data (SeltzerBerger cross section data)
  //           df   (electron density * Migdal constant)
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //           a[nrow][ncol] (alias table) 
  //           q[nrow][ncol] (non-alias probability) 

  G4double dx = (log(xmax) - log(xmin))/nrow;
  G4double xo = log(xmin) + 0.5*dx;

  for(int i = 0; i < nrow ; ++i) {
    G4double x = exp(xo+dx*i);
    BuildPDF_SeltzerBerger(ncol,Z,data,df,xmin,xmax,x,p[i]);
  }

  //build alias tables (alias and probability)
  for(int i = 0; i < nrow ; ++i) {
    BuildAliasTables(ncol,p[i],a[i],q[i]);
  }
}

#endif
