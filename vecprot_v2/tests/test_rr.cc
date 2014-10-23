#include "../inc/rr_pool.h"
#include <thread>
#include <iostream>
#include <vector>
#include <sys/time.h>
#include <math.h>
// g++ test_rr.cc -o test_rr -O3 -pthread -std=c++11

int NTHREADS = 4;
int multthr = 2;
const int MAXCNT   = 450450;
const int NPERSLOT = 500;
rr_pool<int> pool(multthr*NTHREADS,NPERSLOT, new int(10));
std::atomic_flag lock = ATOMIC_FLAG_INIT;
//std::atomic_int stat[NTHREADS];
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

void do_work(int cycles)
{
   double start = 0.1234;
   for (int i=0; i<cycles; ++i) {
      start = sqrt(start);
   }
}

void f(int n)
{
   std::vector<int*> buffer;
   int ncnt = MAXCNT/NTHREADS;
   for (int cnt = 0; cnt < MAXCNT; ++cnt) {
      do_work(1000);
//      std::this_thread::sleep_for(timespan);
//      int slot = pool.next_slot();
//      stat[slot]++;
      int *obj;
      if (((cnt+1)%100) == 0) {
         // dump the buffer
         while (!buffer.empty()) {
	    obj = buffer.back();
       do_work(1000);
	    pool.release(obj);
	    buffer.pop_back();
	 }
      }
      obj = pool.borrow();
      buffer.push_back(obj);
   }
}
 
int main(int argc, char * argv[])
{
   if (argc>=1) NTHREADS = atoi(argv[1]);
   std::cout <<  "NTHREADS = " << NTHREADS << std::endl;
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
//   for (auto n=0; n<NTHREADS; ++n) {
//      assert(stat[n].load()==MAXCNT);
//       std::cout << n << " : " << stat[n].load() << std::endl;
//   }
   pool.statistics();
   std::cout << "wall time: " << wall1-wall0 << "   cpu time: " << cpu1-cpu0 << std::endl;
   std::cout << "speedup = " << (cpu1-cpu0)/(wall1-wall0) << std::endl;
   exit(0);
}
