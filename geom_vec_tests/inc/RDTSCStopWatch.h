#ifndef RDTSCSTOPWATCH_H
#define RDTSCSTOPWATCH_H

// A "timer" based on CPU cycles which can be read out with the rdtsc instruction
// this can only be interpreted as time if cpu powermanagement is switched off in the bios
// as well as turbo-boost and other features
// to convert to time, please adjust the frequency according to your hardware


// as found on wikipedia for a 64bit system   ( we increase the timer overhead artifially, giving us a controllable handle ):
static inline unsigned long long rdtsc()
{
  unsigned int lo, hi;
  __asm__ __volatile__ (
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "rdtsc\n"
      : "=a" (lo), "=d" (hi)
      :
      : "%ebx", "%ecx" );
  return (unsigned long long)hi << 32 | lo;
}

// does not forbid out-of-order execution -- probably ok after a function call return?
static inline unsigned long long rdtsc_nocpuid()
{
  unsigned int lo, hi;
  __asm__ __volatile__ (
      "rdtsc\n"
      : "=a" (lo), "=d" (hi)
      :
      : "%ebx", "%ecx" );
  return (unsigned long long)hi << 32 | lo;
}



struct RDTSCStopWatch 
{
  static const double inversefreq=1./(3.40017e9); // ( for my 3.4 GHz machine, taken from /proc/cpuinfo )
  unsigned long long t1;
  unsigned long long t2;
  void Start(){ t1= rdtsc(); }
  void Stop(){  t2= rdtsc(); }
  void Stop_nocpuid(){  t2= rdtsc_nocpuid(); }

  void HeatUp(){ for(unsigned int i=0;i<5;++i){ Start();Stop(); } }

  unsigned long long GetDeltaTics(){ return t2-t1; }
  double getDeltaSecs() { return (t2-t1)*inversefreq; }

  double getOverhead(int N)
  {
    HeatUp();
    unsigned long long Taccum=0L;
    for(unsigned int i=0; i < N; i++)
      {
	Start();
	Stop();
	Taccum += (t2-t1);
      }
    return Taccum*(inversefreq)/N;
  }

};

#endif
