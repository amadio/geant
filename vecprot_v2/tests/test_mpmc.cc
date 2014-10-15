#include "../inc/mpmc_bounded_queue.h"
#include "../inc/priority_queue.h"
#include "../inc/sync_objects.h"
#include <boost/lockfree/queue.hpp>
#include "../inc/array_lock_free_queue.h"
#include <thread>
#include <iostream>
#include <vector>
#include <sys/time.h>
#include <string>
#include <fstream>
//#include <math.h>
// g++ test_mpmc.cc -o test_mpmc -O3 -pthread -std=c++11 -I$HOME/boost/include -I`root-config --incdir` `root-config --libs`

int NTHREADS = 4;
long *sumw, *sumr;
long sumwt = 0;
long sumrt = 0;
int MAXCNT = 1000000;  // counts per thread
int ncnt = 0;
int NREP = 1;
const std::size_t QSIZE   = 1 << 16;
mpmc_bounded_queue<int> theQ(QSIZE);
priority_queue<int> thePQ(QSIZE);
boost::lockfree::queue<int> theBQ(1000);
//LockFreeQueue<int> theRBQ(NUM_PROD_THREAD, NUM_CONSUM_THREAD);
//lfringqueue < int, 1000> theRBQ; 
//ArrayLockFreeQueue<int, QSIZE> theRBQ;
ArrayLockFreeQueue<int, 1024> theRBQ;
dcqueue<int> theDCQ;
int *values;
std::ofstream ofs1, ofs2, ofs3, ofs4, ofs5;
//std::chrono::milliseconds timespan(1);

//______________________________________________________________________________
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

//______________________________________________________________________________
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

//______________________________________________________________________________
void f1(int n)
{
// mpmc_atomic   
//   size_t maxsize = 0;
//   size_t minsize = 1000000;
   int val;
   int ifailw=0;
   int ifailr=0;
   sumw[n] = 0;
   sumr[n] = 0;
   for (int cnt = 0; cnt < ncnt; ++cnt) {
      while (!theQ.enqueue(cnt)) {}; //{ifailw++;}
      sumw[n] += cnt;
      while (!theQ.dequeue(val)) {}; //{ifailr++;}
      sumr[n] += val;
//      size_t size = theQ.size();
//      if (size>maxsize) maxsize = size;
//      if (size<minsize) minsize = size;
   }
//   std::cout << "ifailw = " << ifailw << "  ifailr = " << ifailr << std::endl;
//   std::cout << "thr " << n << " maxsize = " << maxsize << "  minsize = " << minsize << std::endl;
}

//______________________________________________________________________________
void f2(int n)
{
// priority_atomic  
   int ncnt = MAXCNT/NTHREADS;
//   size_t maxsize = 0;
//   size_t minsize = 1000000;
   int val;
   sumw[n] = 0;
   sumr[n] = 0;
   for (int cnt = 0; cnt < ncnt; ++cnt) {
      bool priority = ((cnt%3)==0);
      while (!thePQ.push(cnt, priority)) {}
      sumw[n] += cnt;
      thePQ.wait_and_pop(val);
      sumr[n] += val;
//      size_t size = theQ.size();
//      if (size>maxsize) maxsize = size;
//      if (size<minsize) minsize = size;
   }
//   std::cout << "thr " << n << " maxsize = " << maxsize << "  minsize = " << minsize << std::endl;
}

//______________________________________________________________________________
void f3(int n)
{
// boost_lf_queue     
   int ncnt = MAXCNT/NTHREADS;
//   size_t maxsize = 0;
//   size_t minsize = 1000000;
   int val;
   sumw[n] = 0;
   sumr[n] = 0;
   for (int cnt = 0; cnt < ncnt; ++cnt) {
      while (!theBQ.push(cnt)) {}
      sumw[n] += cnt;
      while (!theBQ.pop(val)) {}
      sumr[n] += val;
//      size_t size = theQ.size();
//      if (size>maxsize) maxsize = size;
//      if (size<minsize) minsize = size;
   }
//   std::cout << "thr " << n << " maxsize = " << maxsize << "  minsize = " << minsize << std::endl;
}

//______________________________________________________________________________
void f4(int n)
{
// mutex_dcqueue
   int ncnt = MAXCNT/NTHREADS;
   int val;
   sumw[n] = 0;
   sumr[n] = 0;
   for (int cnt = 0; cnt < ncnt; ++cnt) {
      bool priority = ((cnt%3)==0);
      theDCQ.push(values+cnt, priority);
      sumw[n] += cnt;
      val = *theDCQ.wait_and_pop();
      sumr[n] += val;
   }
}

//______________________________________________________________________________
void f5(int n)
{
// Lock free Circular array Queue   
   int ncnt = MAXCNT/NTHREADS;
//   size_t maxsize = 0;
//   size_t minsize = 1000000;
   int val;
   sumw[n] = 0;
   sumr[n] = 0;
   for (int cnt = 0; cnt < ncnt; ++cnt) {
      while (!theRBQ.push(cnt)) {}
      sumw[n] += cnt;
      while (!theRBQ.pop(val)) {}
      sumr[n] += val;
//      size_t size = theQ.size();
//      if (size>maxsize) maxsize = size;
//      if (size<minsize) minsize = size;
   }
//   std::cout << "thr " << n << " maxsize = " << maxsize << "  minsize = " << minsize << std::endl;

}

//______________________________________________________________________________
void performance()
{
//======== mpmc_atomic
   int val = 0;
   sumwt = sumrt = 0;
   for (auto i=0; i<1000; ++i) theQ.enqueue(val);
   std::vector<std::thread> v;
   double wall0 = get_wall_time(); 
   double cpu0  = get_cpu_time();  
   for (auto n = 0; n < NTHREADS; ++n) {
      v.emplace_back(f1, n);
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
   std::cout << "mpmc atomic queue " << std::endl;
   std::cout << "==================" << std::endl;   
   std::cout << "number of reads/writes: " << NTHREADS*ncnt << std::endl;
   if (sumwt==sumrt) std::cout << "checksum matching: " << sumwt << std::endl;
   else              std::cout << "checksum NOT matching: " << std::endl;
   std::cout << "wall time: " << wall1-wall0 << "   cpu time: " << cpu1-cpu0 << std::endl;
   std::cout << "transactions per second: " << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
   ofs1 << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
//   std::cout << "speedup = " << (cpu1-cpu0)/(wall1-wall0) << std::endl;
   
   sumwt = sumrt = 0;
//======== priority_atomic
   val = 0;
   for (auto i=0; i<1000; ++i) theQ.enqueue(val);
   std::vector<std::thread> v2;
   wall0 = get_wall_time(); 
   cpu0  = get_cpu_time();  
   for (auto n = 0; n < NTHREADS; ++n) {
      v2.emplace_back(f2, n);
   }
   for (auto& t : v2) {
      t.join();
   }
   wall1 = get_wall_time(); 
   cpu1  = get_cpu_time();  
   for (auto i = 0; i<NTHREADS; ++i) { 
     sumwt += sumw[i];
     sumrt += sumr[i];
   } 
   extra = 0;
   while (theQ.dequeue(val)) {extra++; sumrt += val;}
   std::cout << "priority atomic queue " << std::endl;
   std::cout << "==================" << std::endl;   
   std::cout << "number of reads/writes: " << NTHREADS*ncnt << std::endl;
   if (sumwt==sumrt) std::cout << "checksum matching: " << sumwt << std::endl;
   else              std::cout << "checksum NOT matching: " << std::endl;
   std::cout << "wall time: " << wall1-wall0 << "   cpu time: " << cpu1-cpu0 << std::endl;
   std::cout << "transactions per second: " << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
   ofs2 << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
//   std::cout << "speedup = " << (cpu1-cpu0)/(wall1-wall0) << std::endl;
   
   sumwt = sumrt = 0;

//======== boost_lf_queue
   val = 0;
   for (auto i=0; i<1000; ++i) theQ.enqueue(val);
   std::vector<std::thread> v3;
   wall0 = get_wall_time(); 
   cpu0  = get_cpu_time();  
   for (auto n = 0; n < NTHREADS; ++n) {
      v3.emplace_back(f3, n);
   }
   for (auto& t : v3) {
      t.join();
   }
   wall1 = get_wall_time(); 
   cpu1  = get_cpu_time();  
   for (auto i = 0; i<NTHREADS; ++i) { 
     sumwt += sumw[i];
     sumrt += sumr[i];
   } 
   extra = 0;
   while (theBQ.pop(val)) {extra++; sumrt += val;}
   std::cout << "Boost Lock free queue " << std::endl;
   std::cout << "==================" << std::endl;   
   std::cout << "number of reads/writes: " << NTHREADS*ncnt << std::endl;
   if (sumwt==sumrt) std::cout << "checksum matching: " << sumwt << std::endl;
   else              std::cout << "checksum NOT matching: " << std::endl;
   std::cout << "wall time: " << wall1-wall0 << "   cpu time: " << cpu1-cpu0 << std::endl;
   std::cout << "transactions per second: " << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
   ofs3 << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
//   std::cout << "speedup = " << (cpu1-cpu0)/(wall1-wall0) << std::endl;
   
   sumwt = sumrt = 0;

//======== mutex_dcqueue
   val = 0;
   for (auto i=0; i<1000; ++i) theDCQ.push(&val);
   std::vector<std::thread> v1;
   wall0 = get_wall_time(); 
   cpu0  = get_cpu_time();  
   for (auto n = 0; n < NTHREADS; ++n) {
      v1.emplace_back(f4, n);
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
//   extra = 0;
   int *value1 = 0;
   while ((value1=theDCQ.try_pop())) {extra++; sumrt += *value1;}
   std::cout << "mutex_dcqueue" << std::endl;
   std::cout << "===================================" << std::endl;   
   std::cout << "number of reads/writes: " << NTHREADS*ncnt << std::endl;
   if (sumwt==sumrt) std::cout << "checksum matching: " << sumwt << std::endl;
   else              std::cout << "checksum NOT matching: " << std::endl;
   std::cout << "wall time: " << wall1-wall0 << "   cpu time: " << cpu1-cpu0 << std::endl;
   std::cout << "transactions per second: " << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
   ofs4 << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
 //  std::cout << "speedup = " << (cpu1-cpu0)/(wall1-wall0) << std::endl;
//========
   sumwt = sumrt = 0;
   val = 0;
   for (auto i=0; i<1000; ++i) theQ.enqueue(val);
   std::vector<std::thread> v5;
   wall0 = get_wall_time(); 
   cpu0  = get_cpu_time();  
   for (auto n = 0; n < NTHREADS; ++n) {
      v5.emplace_back(f5, n);
   }
   for (auto& t : v5) {
      t.join();
   }
   wall1 = get_wall_time(); 
   cpu1  = get_cpu_time();  
   for (auto i = 0; i<NTHREADS; ++i) { 
     sumwt += sumw[i];
     sumrt += sumr[i];
   } 
   extra = 0;
   while (theRBQ.pop(val)) {extra++; sumrt += val;}
   std::cout << "Lock free Circular array Queue " << std::endl;
   std::cout << "==================" << std::endl;   
   std::cout << "number of reads/writes: " << NTHREADS*ncnt << std::endl;
   if (sumwt==sumrt) std::cout << "checksum matching: " << sumwt << std::endl;
   else              std::cout << "checksum NOT matching: " << std::endl;
   std::cout << "wall time: " << wall1-wall0 << "   cpu time: " << cpu1-cpu0 << std::endl;
   std::cout << "transactions per second: " << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
   ofs5 << 2*NTHREADS*ncnt/(wall1-wall0) << std::endl;
//   std::cout << "speedup = " << (cpu1-cpu0)/(wall1-wall0) << std::endl;
}

//______________________________________________________________________________
int main(int argc, char * argv[])
{
   std::string ARCH;
   if (argc>1) NTHREADS = atoi(argv[1]);
   if (argc>2) ARCH = argv[2];
   if (argc>3) NREP = atoi(argv[3]);
   std::cout <<  "NTHREADS = " << NTHREADS << std::endl;
   sumw = new long[NTHREADS];
   sumr = new long[NTHREADS];
   ncnt = MAXCNT/NTHREADS;
   std::atomic<size_t> atomic_test(0);
   bool lock_free = atomic_test.is_lock_free();
   if (lock_free) std::cout << "std::atomic is lock free" << std::endl;
   else           std::cout << "std::atomic is NOT lock free" << std::endl;
   values = new int[ncnt];
   for (auto i=0; i<ncnt; ++i) values[i] = i;

   std::ifstream if1 ("q_mpmc_atomic.txt");
   bool exist = if1.good();
   if (exist) if1.close();

   ofs1.open ("q_mpmc_atomic.txt", std::ofstream::app);
   ofs2.open ("q_priority_atomic.txt", std::ofstream::app);
   ofs3.open ("q_boost_lockfree.txt", std::ofstream::app);
   ofs4.open ("q_mutex_dcqueue.txt", std::ofstream::app);
   ofs5.open ("q_carray_lockfree.txt", std::ofstream::app);
   if (!exist) {
     ofs1 << "#Transactions per second " << ARCH << std::endl;
     ofs2 << "#Transactions per second " << ARCH << std::endl;
     ofs3 << "#Transactions per second " << ARCH << std::endl;
     ofs4 << "#Transactions per second " << ARCH << std::endl;
     ofs5 << "#Transactions per second " << ARCH << std::endl;
   }

   ofs1 << "NTHREADS " << NTHREADS << std::endl;
   ofs2 << "NTHREADS " << NTHREADS << std::endl;
   ofs3 << "NTHREADS " << NTHREADS << std::endl;
   ofs4 << "NTHREADS " << NTHREADS << std::endl;
   ofs5 << "NTHREADS " << NTHREADS << std::endl;
   for (auto icall=0; icall<NREP; ++icall) performance();
   ofs1.close();
   ofs2.close();
   ofs3.close();
   ofs4.close();
   ofs5.close();
   return 0;
}      
