
#ifndef LINUXTIMER_H
#define LINUXTIMER_H

#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>        /* Definition of uint64_t */
#include <assert.h>

struct Timer {

	Timer() {
		reset();
    accum_usecs = 0;
	}

	double elapsed_ms() const
	{
    struct timespec now;

    int ercd = clock_gettime(CLOCK_MONOTONIC, &now);
    assert(ercd != -1);

    long long delta = (now.tv_sec * 1000000) + (now.tv_nsec / 1000) - start_usecs;

    return ((double)delta);
	}

  double accum()
  {
    struct timespec now;

    int ercd = clock_gettime(CLOCK_MONOTONIC, &now);
    assert(ercd != -1);

    long long delta = (now.tv_sec * 1000000) + (now.tv_nsec / 1000) - start_usecs;

    accum_usecs += delta;

    return ((double)accum_usecs);
  }

	void reset()
  {
    struct timespec start;

    int ercd = clock_gettime(CLOCK_MONOTONIC, &start);
    assert(ercd != -1);

    start_usecs = start.tv_sec * 1000000 + (start.tv_nsec / 1000);
	}

	long long start_usecs;
  long long accum_usecs;
};

#endif
