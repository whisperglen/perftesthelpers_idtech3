
#ifndef TIMING_H_
#define TIMING_H_

#ifdef _MSC_VER
#include "timer_win.h"
#include "cputime_win.h"
#endif

#ifdef __GNUC__
#include "timer_linux.h"
#include "cputime_dummy.h"
#endif

#endif
