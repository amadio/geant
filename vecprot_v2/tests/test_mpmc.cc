#include "../inc/mpmc_bounded_queue.h"
#include "../inc/sync_objects.h"
#include <thread>
#include <iostream>
#include <vector>
#include <sys/time.h>
#include <string.h>
//#include <math.h>
// g++ test_mpmc.cc -o test_mpmc -O3 -pthread -std=c++11 -I`root-config --incdir` `root-config --libs`

int NTHREADS = 4;
long *sumw, *sumr;
long sumwt = 0;
long sumrt = 0;
int MAXCNT = 1000000;  // counts per thread
int multthr = 2;
const std::size_t QSIZE   = 1 << 14;
mpmc_bounded_queue<int> theQ(QSIZE);
dcqueue<int> theDCQ;
int *values;

//std::chrono::milliseconds timespan(1);

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

void f(int n)
{
   int ncnt = MAXCNT/NTHREADS;
   size_t maxsize = 0;
   size_t minsize = 1000000;
   int val;
   sumw[n] = 0;
   sumr[n] = 0;
   for (int cnt = 0; cnt < ncnt; ++cnt) {
      while (!theQ.enqueue(cnt)) {}
      sumw[n] += cnt;
      while (!theQ.dequeue(val)) {}
      sumr[n] += val;
      size_t size = theQ.size();
      if (size>maxsize) maxsize = size;
      if (size<minsize) minsize = size;
   }
   std::cout << "thr " << n << " maxsize = " << maxsize << "  minsize = " << minsize << std::endl;
}

void g(int n)
{
   int ncnt = MAXCNT/NTHREADS;
   int val;
   sumw[n] = 0;
   sumr[n] = 0;
   for (int cnt = 0; cnt < ncnt; ++cnt) {
      theDCQ.push(values+cnt);
      sumw[n] += cnt;
      val = *theDCQ.wait_and_pop();
      sumr[n] += val;
   }
}
 
int main(int argc, char * argv[])
{
   if (argc>1) NTHREADS = atoi(argv[1]);
   std::cout <<  "NTHREADS = " << NTHREADS << std::endl;
   sumw = new long[NTHREADS];
   sumr = new long[NTHREADS];
   int ncnt = MAXCNT/NTHREADS;
   values = new int[ncnt];
   for (auto i=0; i<ncnt; i++) values[i] = i;
   int val = 0;
   for (auto i=0; i<1000; ++i) theQ.enqueue(val);
   std::vector<std::thread> v;
   double wall0 = get_wall_time(); 
   double cpu0  = get_cpu_time();  
   for (auto n = 0; n < NTHREADS; ++n) {
      v.emplace_back(f, n);
   }
   for (auto& t : v) {
      t.join();
   }
   double wall1 = get_wall_time(); 
   double cpu1  = get_cpu_time();  
   for (auto i = 0; i<NTHREADS; ++i) { 
     sumwt += sumw[i];
     sumrt += sumr[i];
   } 
   int extra = 0;
   while (theQ.dequeue(val)) {extra++; sumrt += val;}
   std::cout << "extra = " << extra << std::endl;
   std::cout << "number of reads/writes: " << NTHREADS*ncnt << std::endl;
   std::cout << "sumw = " << sumwt << "  sumr = " << sumrt << std::endl;
   std::cout << "wall time: " << wall1-wall0 << "   cpu time: " << cpu1-cpu0 << std::endl;
   std::cout << "transactions per second: " << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
//   std::cout << "speedup = " << (cpu1-cpu0)/(wall1-wall0) << std::endl;
   
   sumwt = sumrt = 0;
   val = 0;
   for (auto i=0; i<1000; ++i) theDCQ.push(&val);
   std::vector<std::thread> v1;
   wall0 = get_wall_time(); 
   cpu0  = get_cpu_time();  
   for (auto n = 0; n < NTHREADS; ++n) {
      v1.emplace_back(g, n);
   }
   for (auto& t : v1) {
      t.join();
   }
   wall1 = get_wall_time(); 
   cpu1  = get_cpu_time();  
   for (auto i = 0; i<NTHREADS; ++i) { 
     sumwt += sumw[i];
     sumrt += sumr[i];
   } 
   extra = 0;
   int *value1 = 0;
   while ((value1=theDCQ.try_pop())) {extra++; sumrt += *value1;}
   std::cout << "extra = " << extra << std::endl;
   std::cout << "number of reads/writes: " << NTHREADS*ncnt << std::endl;
   std::cout << "sumw = " << sumwt << "  sumr = " << sumrt << std::endl;
   std::cout << "wall time: " << wall1-wall0 << "   cpu time: " << cpu1-cpu0 << std::endl;
   std::cout << "transactions per second: " << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
 //  std::cout << "speedup = " << (cpu1-cpu0)/(wall1-wall0) << std::endl;
   
   exit(0);
}
