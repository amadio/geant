#ifndef GEANT_CONDITIONLOCKER
#define GEANT_CONDITIONLOCKER

#include <mutex>
#include <condition_variable>

// A simple semaphore based on condition. The controlled threads have to
// call the Wait() function and the controller can wake up one waiting
// thread at a time calling StartOne() method, or all by calling StartAll().

struct condition_locker {
  std::mutex              m;
  std::condition_variable cv;
  std::atomic<bool>       start;
  
  condition_locker() : m(), cv(), start(false) {}
  void Wait()     {
                    std::unique_lock<std::mutex> lk(m);
                    while (!start.load()) cv.wait(lk);
                    start.store(false);
                    lk.unlock();
                    cv.notify_one();
                  };
  void StartOne() { 
                    std::unique_lock<std::mutex> lk(m);
                    start.store(true);
                    cv.notify_one();
                  };  
  void StartAll() { 
                    std::unique_lock<std::mutex> lk(m);
                    start.store(true);
                    cv.notify_all();
                  };  
};
#endif
