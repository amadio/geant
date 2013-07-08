#include "TBBStopWatch.h"
#include <iostream>
#include <iomanip>

StopWatch timer;
StopWatch outertimer;

#define N 100000
double overhead[N];
double total=0.;

int main()
{
  std::setprecision(20);

  timer.HeatUp();
  for(unsigned int i=0;i<N;++i)
    {
      timer.Start();
      timer.Stop();
      std::cerr << timer.getDeltaSecs() << std::endl;
      //      overhead[i]=timer.getDeltaSecs();
    }

  /*
  outertimer.Start();
  for(unsigned int i=0;i<N;++i)
    {
      timer.Start();
      timer.Stop();
      total += timer.getDeltaSecs();
    }
  outertimer.Stop();
  std::cerr << outertimer.getDeltaSecs()/N << " " << total/N << std::endl;
  */
  /* for(unsigned int i=0;i<N;++i)
    {
      std::cerr << overhead[i] << std::endl;
    } */
}
