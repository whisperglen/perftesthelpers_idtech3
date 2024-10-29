
#ifndef MYTESTIMER_H
#define MYTESTIMER_H

#include <Windows.h>
#include <cassert>

struct Timer {

	double frequency;

	long long start_ticks;
	long long accum_ticks;

	Timer() {
		reset();
		accum_ticks = 0;
	}

	double elapsed_ms() const
	{
		LARGE_INTEGER current_ticks;

		BOOL ercd = QueryPerformanceCounter(&current_ticks);
		assert(ercd != 0);

		LONGLONG delta = current_ticks.QuadPart - start_ticks;

		return (1.0e3 * delta * frequency);
	}

	double accum()
	{
		LARGE_INTEGER current_ticks;

		BOOL ercd = QueryPerformanceCounter(&current_ticks);
		assert(ercd != 0);

		LONGLONG delta = current_ticks.QuadPart - start_ticks;

		accum_ticks += delta;

		return (1.0e3 * accum_ticks * frequency);
	}

	void reset()
	{
		LARGE_INTEGER f;

		BOOL ercd = QueryPerformanceFrequency(&f);
		assert(ercd != 0);

		frequency = 1.0 / f.QuadPart;

		LARGE_INTEGER st;

		ercd = QueryPerformanceCounter(&st);
		assert(ercd != 0);

		start_ticks = st.QuadPart;
	}

};

#endif