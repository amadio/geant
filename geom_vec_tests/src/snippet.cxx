#include <iostream>
#include <bitset>

using namespace std;

#include "tbb/tick_count.h"" // timing from Intel TBB 
#include <cassert>

#define NREP 100000000

struct TStopWatch 
{
  tbb::tick_count t1;
  tbb::tick_count t2;
  void Start(){ t1=tbb::tick_count::now(); }
  void Stop(){ t2=tbb::tick_count::now(); }
  void Reset(){ /* */ ;}
  void Print(){  std::cerr << (t2 - t1).seconds() << std::endl; }
  double getDeltaSecs() { return (t2-t1).seconds(); }
};


double foo(double x, double y)
{
  long long zero = (x==0)? 0 : 1;
  
  double s=2.;
  if(x!=0)
    {
      int sign=((long long) x)>0;
      s= y/x * ( sign  - (!sign) );
    }
  return s;
}

typedef union
{
  double d;
  unsigned long long l;
} doublelongunion;


double baz(double x, double y)
{
  doublelongunion u1, u2;
  u2.d=x;
  //  std::bitset<sizeof(double)*8> bits(u2.l);
  //  cerr << "x " << bits << endl;

  unsigned long long bitmask = (x==0.)? 0 : -1; // this sets bitmask to 0000000000 or 1111111111
  unsigned long long flip = -1;


  //std::bitset<sizeof(double)*8> bits(bitmask);
  //cerr << "bitmask " << bits << endl;

  double s=2.;

  double oneoverx=1./x;
  double xtimesx=x*oneoverx;
  cerr << "x " << x << " " << xtimesx << endl;

  u1.d=s;
  u1.l=(bitmask^flip) & u1.l;

  //  std::bitset<sizeof(double)*8> bits2(u2.l);
  //  cerr << "u2.d " << bits2 << endl;

  u2.d=1./x;
  u2.l=(bitmask) & u2.l;

  //  std::bitset<sizeof(double)*8> bits3(u2.l);
  //  cerr << "u2.d " << bits3 << endl;
  int sign=((long long) x)>0;
  double r = u1.d + ( u2.d * y * ( sign  - (!sign) ) );
  return r;
}


double bar(double x, double y)
{
  double s=2.;
  if(x!=0)
    {
      s= (x>0)? y/x : -y/x;
    }
  return s;
}



int main()
{
  volatile double x;
  x=foo(1.2,1.3);
  cerr << x << " " << endl;
  x=bar(1.2,1.3);
  cerr << x << " " << endl;

  x=baz(1.2,1.3);
  cerr << x << " " << endl;


  x=foo(0.,1.3);
  cerr << x << " " << endl;

  x=bar(0.,1.3);
  cerr << x << " " << endl;

  x=baz(0.,1.3);
  cerr << x << " " << endl;


  x=foo(-2.,1.3);
  cerr << x << " " << endl;
  x=bar(-2.,1.3);
  cerr << x << " " << endl;
  x=baz(-2.,1.3);
  cerr << x << " " << endl;

  return 1;

  TStopWatch tt;
  
  tt.Reset();

  double *r= new double[NREP];
  tt.Start();
  for(unsigned int i=0;i<NREP;++i)
    {
      r[i]=foo(0.,1.3*i);
    }
  tt.Stop();
  cerr << tt.getDeltaSecs() << endl;

  tt.Start();
  for(unsigned int i=0;i<NREP;++i)
    {
      r[i]=baz(0.,1.3*i);
    }
  tt.Stop();
  cerr << tt.getDeltaSecs() << endl;


  tt.Start();
  for(unsigned int i=0;i<NREP;++i)
    {
      r[i]=bar(0.,1.3*i);
    }
  tt.Stop();
  cerr << tt.getDeltaSecs() << endl;
}
