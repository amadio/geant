#ifndef TBBSTOPWATCH_H
#define TBBSTOPWATCH_H

#include "tbb/tick_count.h" // timing from Intel TBB 

struct StopWatch 
{
  tbb::tick_count t1;
  tbb::tick_count t2;
  void Start(){ t1=tbb::tick_count::now(); }
  void Stop(){ t2=tbb::tick_count::now(); }
  void Reset(){ /* */ ;}
  //  void Print(){  std::cerr << (t2 - t1).seconds() << std::endl; }

  void HeatUp(){ for(unsigned int i=0;i<5;++i) { this->Start();this->Stop(); } }

  double getDeltaSecs() { return (t2-t1).seconds(); }
};

#endif
