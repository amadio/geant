int sum;
int a[128];
void foo (void)
{
  int i;

  for (i = 0; i < 64; i++)
    {
      sum += a[2*i];
      sum += a[2*i+1];
    }
}


void foo2 (int np)
{
  int i;

  for (i = 0; i < np; i++)
    {
      sum += a[3*i];
      sum += a[3*i+1];
      sum += a[3*i+2];
    }
}


void foo3 (int np, int * result)
{
  int i;

  for (i = 0; i < np; i++)
    {
      result[3*i] = a[3*i];
      result[3*i+1] = a[3*i+1];
      result[3*i+2] = a[3*i+2];
    }
}

void foo4 (int np, double * result)
{
  int i;

  for (i = 0; i < np; i++)
    {
      result[3*i] = a[3*i];
      result[3*i+1] = a[3*i+1];
      result[3*i+2] = a[3*i+2];
    }
}


void foo5 (int np, double * result, double *b )
{
  int i;

  for (i = 0; i < np; i++)
    {
      result[3*i] = b[3*i];
      result[3*i+1] = b[3*i+1];
      result[3*i+2] = b[3*i+2];
    }
}

// gcc DOES NOT vectorize this ( alignment missing )
void foo6 (int np, double *__restrict__ result, const double *__restrict__ b )
{
  int i;

  for (i = 0; i < np; i++)
    {
      result[3*i] = b[3*i];
      result[3*i+1] = b[3*i+1];
      result[3*i+2] = b[3*i+2];
    }
}

// gcc vectorizes this ! (difference is made not by restrict but by assumed aligned
void foo7 (int np, double *__restrict__ result, const double *__restrict__ b )
{
  int i;
#ifndef __INTEL_COMPILER
  const double *x = (const double*) __builtin_assume_aligned (b, 16);
#else
  const double *x=b;
#endif


  for (i = 0; i < np; i++)
    {
      result[3*i] = x[3*i];
      result[3*i+1] = x[3*i+1];
      result[3*i+2] = x[3*i+2];
    }
}

// gcc vectorized this
void foo8 (int np, double *__restrict__ result, const double *__restrict__ b )
{
  int i;
  double f[3];
#ifndef __INTEL_COMPILER
  const double *x = (const double*) __builtin_assume_aligned (b, 16);
#else
  const double *x=b;
#endif

#pragma ivdep
  for (i = 0; i < np; i++)
    {
      result[3*i] = x[3*i]-f[0];
      result[3*i+1] = x[3*i+1]-f[1];
      result[3*i+2] = x[3*i+2]-f[2];
    }
}

// gcc does not vectorize this ; intel vectorizes it
void foo9 (int np, double *__restrict__ result, const double *__restrict__ b )
{
  int i;
  double f[3];

#ifndef __INTEL_COMPILER
  const double *x = (const double*) __builtin_assume_aligned (b, 16);
  double *c1 = (double*) __builtin_assume_aligned (result, 16);
#else
  const double *x = b;
  double *c1=result;
#endif


  for (i = 0; i < np; i++)
    {
      double w,y,z;
      w= x[3*i]-f[0];
      y= x[3*i+1]-f[1];
      z= x[3*i+2]-f[2];
      c1[3*i] = w+y+z;
      c1[3*i+1]=0.;
      c1[3*i+2]=0.;
    }
}


typedef struct
{
  double a,b,c;
} fooT;

// gcc does not vectorize this ; intel vectorizes it
void foo10 (int np, double *__restrict__ result, const double *__restrict__ b )
{
  int i;
  double f[3];

#ifndef __INTEL_COMPILER
  fooT *x1 = (fooT *) __builtin_assume_aligned (b, 16);
#else
  fooT *x1 = (fooT *) b;
  double *c1=result;
#endif


  for (i = 0; i < np; i++)
    {
      double w,y,z;
      w= x1[i].a;
      y= x1[i].b;
      z= x1[i].c;
      result[i] = w+y+z;
    }
}



void foo10 (int np, double *__restrict__ result, const double *__restrict__ b, const double *__restrict__ a, const double *__restrict__ c )
{
  int i;
  double f[3];

#ifndef __INTEL_COMPILER
  //  double *x1 = (double *) __builtin_assume_aligned (a, 16);
  //  double *x2 = (double *) __builtin_assume_aligned (b, 16);
  //  double *x3 = (double *) __builtin_assume_aligned (c, 16);
  double *x1 = (double *) a;
  double *x2 = (double *) b;
  double *x3 = (double *) c;
#else
  double *x1 = (double *) a;
  double *x2 = (double *) b;
  double *x3 = (double *) c;
#endif
  for (i = 0; i < np; i++)
    {
      double w,y,z;
      w= x1[i];
      y= x2[i];
      z= x3[i];
      result[i] = w+y+z;
    }
}

// example with if statements
void foo11 (int np, double *__restrict__ result, const double *__restrict__ b, const double *__restrict__ a, const double *__restrict__ c )
{
  int i;
  double f[3];

#ifndef __INTEL_COMPILER
   double *x1 = (double *) __builtin_assume_aligned (a, 16);
   double *x2 = (double *) __builtin_assume_aligned (b, 16);
   double *x3 = (double *) __builtin_assume_aligned (c, 16);
   //double *x1 = (double *) a;
   //double *x2 = (double *) b;
   //double *x3 = (double *) c;
#else
  double *x1 = (double *) a;
  double *x2 = (double *) b;
  double *x3 = (double *) c;
#endif
  for (i = 0; i < np; i++)
    {
      double w,y,z;
      w= x1[i];
      y= x2[i];
      z= x3[i];

      double a=(z!=0)? z : 1./z; // this vectorizes with compiler option -ffast-math
 
      result[i] = a+y+z;
    }
}




// example with if statements
void foo12 (int np, double *__restrict__ result, const double *__restrict__ b, const double *__restrict__ a, const double *__restrict__ c )
{
  int i;
  double f[3];

#ifndef __INTEL_COMPILER
   double *x1 = (double *) __builtin_assume_aligned (a, 16);
   double *x2 = (double *) __builtin_assume_aligned (b, 16);
   double *x3 = (double *) __builtin_assume_aligned (c, 16);
   //double *x1 = (double *) a;
   //double *x2 = (double *) b;
   //double *x3 = (double *) c;
#else
  double *x1 = (double *) a;
  double *x2 = (double *) b;
  double *x3 = (double *) c;
#endif
  for (i = 0; i < np; i++)
    {
      double w,y,z;
      w= x1[i];
      y= x2[i];
      z= x3[i];

      double a;
      if(y!=0){
	a=(z!=0)? z : 1./z; // this vectorizes with compiler option -ffast-math
      }
      else
	{
	  a=1.;
	}
      result[i] = a+y+z;
    }
}


// example with inner loop
void foo13 (int np, double *__restrict__ result, const double *__restrict__ b, const double *__restrict__ a, const double *__restrict__ c )
{
  int i;
  double f[3];

#ifndef __INTEL_COMPILER
   double *x1 = (double *) __builtin_assume_aligned (a, 16);
   double *x2 = (double *) __builtin_assume_aligned (b, 16);
   double *x3 = (double *) __builtin_assume_aligned (c, 16);
   //double *x1 = (double *) a;
   //double *x2 = (double *) b;
   //double *x3 = (double *) c;
#else
  double *x1 = (double *) a;
  double *x2 = (double *) b;
  double *x3 = (double *) c;
#endif
  double r[3];
  for (i = 0; i < np; i++)
    {
      for(unsigned int j=0;j<3;j++)
	{
	  r[j]=-2*j+x1[i];
	}
      double w,y,z;
      w= x1[i]-r[0];
      y= x2[i]-r[1];
      z= x3[i]-r[2];

      result[i]=w+y+z;
    }
}


// example with double array ( does not get vectorized currently) 
void foo14 (int np, double *__restrict__ result, const double *__restrict__ b, const double *__restrict__ a, const double *__restrict__ c )
{
  int i;
  double f[3];

#ifndef __INTEL_COMPILER
  double *x[3];
  x[0] = (double *) __builtin_assume_aligned (a, 16);
  x[1] = (double *) __builtin_assume_aligned (b, 16);
  x[3] = (double *) __builtin_assume_aligned (c, 16);
#else
  double *x[3];
  x[0] = (double *) a;
  x[1] = (double *) b;
  x[2] = (double *) c;
#endif
  double r[3];
  for (i = 0; i < np; i++)
    {
      //    for(unsigned int j=0;j<3;j++)
      // 	{
      //	  r[j]=-2*j+x[j][i];
      //	}
      // r[0]=x[0][i];
      // r[1]=-2+x[1][i];
      // r[2]=-4+x[2][i];
    
      double w,y,z;
      w= x[0][i]-r[0];
      y= x[1][i]-r[1];
      z= x[2][i]-r[2];

      result[i]=w+y+z;
    }
}

