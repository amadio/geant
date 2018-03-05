#ifndef ET_TIMER_H
#define ET_TIMER_H

#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
/**
 * \file Timer.h  A flexible toolset for timing pieces of code.
 */

#ifdef _MSC_VER
#include <intrin.h>
#pragma intrinsic(__rdtsc)
#endif

#ifndef ET_TIMER_RESULT_VAR
#error undefined ET_TIMER_RESULT_VAR
/**
 * This macro provides the name of a result variable to the Timer tools.
 * The user must define it BEFORE including et/Timer.h. E.g.
 * ~~~
        #define ET_TIMER_RESULT_VAR_TYPE float*
        #define ET_TIMER_RESULT_VAR             a
 * ~~~
 */
#define ET_TIMER_RESULT_VAR
#endif

#ifndef ET_TIMER_RESULT_VAR_TYPE
#error undefined ET_TIMER_RESULT_VAR_TYPE
/**
 * The user must define this macro to provide the type of the result variable
 * (#ET_TIMER_RESULT_VAR) to the Timer tools.
 * ~~~
        #define ET_TIMER_RESULT_VAR_TYPE float*
        #define ET_TIMER_RESULT_VAR             a
 * ~~~
 */
#define ET_TIMER_RESULT_VAR_TYPE
#endif

#ifdef DEBUG
#define NEVER_TRUE k % kPrint == 0
#else
#define NEVER_TRUE k % kPrint == -1
#endif

/**
 * Convenient macro to time a piece of code.
 * Use this macro to time a piece of code and print the results
 * \param TITLE descriptive string for the piece of code to be timed
 * \param CODE  code to be timed
 * \param NREPS request timing for NREPS repetitions of CODE
 * \param PRINT timing results are printed if PRINT is convertible to true
 *
 * E.g. this times a for loop that writes "Hello" ten times, repeating it 100 times, and prints
 * results with title "10 x Hello":
 ~~~
 ET_TIME_THIS( "10 x Hello",
 , for(int i=0; i<10; ++i) {
       std::cout << "Hello" << std::endl;
   }
 , 100
 , true
 )
 ~~~
 \remark If the macro ET_NO_TIMING is defined the code is not timed, only executed.
 */
#ifndef ET_NO_TIMING
#define ET_TIME_THIS(TITLE, CODE, NREPS, PRINT)            \
  {                                                        \
    ET::Timer timer;                                       \
    size_t kPrint = (NREPS < 10 ? 1 : NREPS / 10);         \
    if (NREPS > 1) {                                       \
      timer.start();                                       \
      for (int k = 0; k < NREPS;) {                        \
        CODE++ k;                                          \
        if (NEVER_TRUE) ET::dummy(k, ET_TIMER_RESULT_VAR); \
      }                                                    \
      timer.stop();                                        \
    } else {                                               \
      timer.start();                                       \
      CODE timer.stop();                                   \
    }                                                      \
    timer.push_back(TITLE, NREPS);                         \
    if (PRINT) ET::Timer::print_last();                    \
  }
#else
#define ET_TIME_THIS(TITLE, CODE, NREPS, PRINT) \
  {                                             \
    CODE                                        \
  }
#endif
/**
 * Print all timing results to output stream OUT
 * \param[in] OUT output stream to which the results are sent.
 *
 * E.g.
 ~~~
 ET_PRINT_ALL( std::out )
 ~~~
 \remark If the macro ET_NO_TIMING is defined nothing is printed.
 */
#ifndef ET_NO_TIMING
#define ET_PRINT_ALL(OUT) ET::Timer::print_all(OUT);
#else
#define ET_PRINT_ALL(OUT)
#endif
/**
 * namespace ET for general/generic utilities
 */
namespace ET { //-----------------------------------------------------------------------------
               /**
                * Helper function that prevents that the compiler optimizes away the repetitions.
                */
void dummy(int k, ET_TIMER_RESULT_VAR_TYPE ET_TIMER_RESULT_VAR)
{
  std::cout << "..." << k;
#ifdef DEBUG
  std::cout << std::endl;
#endif
}
//-----------------------------------------------------------------------------
/**
 * Class for timing pieces of code with repetitions, the macro #ET_TIME_THIS
 * provides a very practical way to use it.
 */
class Timer {
public:
  /**
         *  Start the timer.
         */ void start();
  /**
         *  Stop the timer.
         */ void stop();
  /**
         *  Return the number of cycles between the last start() and stop().
         */ unsigned long long cycles() const;
  /**
         *  Return the number of seconds between the last start() and stop().
         */ double seconds() const;
  /**
         *  Compute the average clock frequency for the i-th run.
         */ static double GHz(int i);
  /**
         * Clear all stored results.
         */ static void reset();
  /**
         *  Store results of last run.
         *  \param[in] title descriptive string for the last run
         *  \param[in] nreps number of repetitions for the last run
         */ void push_back(std::string const &title, size_t nreps);
  /**
         *  Index of the run that will be used as a reference for computing the speedup factor.
         */ static size_t speedupReference;
  /**
         *  Print results of i-th run.
         */ static void print_i(int i, std::ostream &o = std::cout, int iSpeedupRef = 0);
  /**
         *  Print results of last run.
         */ static void print_last(std::ostream &o = std::cout);
  /**
         *  Print results of all runs.
         */ static void print_all(std::ostream &o = std::cout);

private:
  union Data {
    unsigned long long a;
    unsigned int b[2];
  } m_start, m_end;
  std::chrono::time_point<std::chrono::system_clock> chrono_start_, chrono_end_;
  static std::vector<std::string> title_list;         // title of each run
  static std::vector<unsigned long long> cycles_list; // #cycles of each run
  static std::vector<double> seconds_list;            // #seconds of each run
  static std::vector<size_t> repeats_list;            // #repetitions of each run
};
//-----------------------------------------------------------------------------
void Timer::print_i(int i, std::ostream &o, int iSpeedupRef)
{
  o << title_list[i] << " : " << cycles_list[i] / repeats_list[i] << " cycles/repetition, " << std::setprecision(4)
    << seconds_list[i] / repeats_list[i] << " seconds/repetition, " << std::setprecision(3)
    << (float)cycles_list[iSpeedupRef] / cycles_list[i] << "x speedup, " << std::setprecision(3) << GHz(i) << " GHz, "
    << repeats_list[i] << " repetitions.";
}
//-----------------------------------------------------------------------------
void Timer::print_last(std::ostream &o)
{
  print_i(title_list.size() - 1);
  o << std::endl;
}
//-----------------------------------------------------------------------------
void Timer::print_all(std::ostream &o)
{
  for (int i = 0; i < title_list.size(); ++i) {
    print_i(i);
    o << '\n';
  }
  o << std::flush;
}
//-----------------------------------------------------------------------------
void Timer::push_back(std::string const &title, size_t nreps)
{
  title_list.push_back(title);
  cycles_list.push_back(this->cycles());
  seconds_list.push_back(this->seconds());
  repeats_list.push_back(nreps);
}
//-----------------------------------------------------------------------------
double Timer::GHz(int i)
{
  return 1.e-9 * cycles_list[i] / seconds_list[i];
}
//-----------------------------------------------------------------------------
void Timer::reset()
{
  title_list.clear();
  cycles_list.clear();
  seconds_list.clear();
  repeats_list.clear();
}
//-----------------------------------------------------------------------------
inline void Timer::start()
{
  this->chrono_start_ = std::chrono::system_clock::now();
#ifdef VC_IMPL_MIC
  asm volatile("xor %%eax,%%eax\n\tcpuid\n\trdtsc" : "=a"(m_start.b[0]), "=d"(m_start.b[1])::"ebx", "ecx");
#elif defined _MSC_VER
  unsigned int tmp;
  m_start.a = __rdtscp(&tmp);
#else
  asm volatile("rdtscp" : "=a"(m_start.b[0]), "=d"(m_start.b[1])::"ecx");
#endif
}

inline void Timer::stop()
{
#ifdef VC_IMPL_MIC
  asm volatile("xor %%eax,%%eax\n\tcpuid\n\trdtsc" : "=a"(m_end.b[0]), "=d"(m_end.b[1])::"ebx", "ecx");
#elif defined _MSC_VER
  unsigned int tmp;
  m_end.a = __rdtscp(&tmp);
#else
  asm volatile("rdtscp" : "=a"(m_end.b[0]), "=d"(m_end.b[1])::"ecx");
#endif
  this->chrono_end_ = std::chrono::system_clock::now();
}

inline unsigned long long Timer::cycles() const
{
  return m_end.a - m_start.a;
}

inline double Timer::seconds() const
{
  std::chrono::duration<double> duration = this->chrono_end_ - this->chrono_start_;
  return duration.count();
}
//-----------------------------------------------------------------------------
// should be in a .cpp file actually, but then it should become a library.
std::vector<std::string> Timer::title_list;
std::vector<unsigned long long> Timer::cycles_list;
std::vector<double> Timer::seconds_list;
std::vector<size_t> Timer::repeats_list;
size_t Timer::speedupReference = 0;
//-----------------------------------------------------------------------------
}
#endif /* ET_TIMER_H */
