#ifndef _STOPWATCH_H_
#define _STOPWATCH_H_

#include "stopwatch_functions.h"

#ifdef _WIN32
typedef StopWatchWin StopWatch;
#else
typedef StopWatchLinux StopWatch;
#endif

#endif // _STOPWATCH_H_

