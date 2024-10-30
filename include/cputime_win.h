
#ifndef CPUTIME_H
#define CPUTIME_H

#include <Windows.h>
#include <cassert>

struct Cputime {

    HANDLE process;

    long long start_ticks;
    long long accum_ticks;

	Cputime() {
        process = GetCurrentProcess();
		reset();
        accum_ticks = 0;
	}

	double elapsed() const
	{
        ULONG64 cycles = 0;
        BOOL ercd = QueryProcessCycleTime(process, &cycles);
        assert(ercd != 0);

		LONGLONG delta = cycles - start_ticks;

		return ((double)delta / 1000000);
	}

    double accum()
    {
        ULONG64 cycles = 0;
        BOOL ercd = QueryProcessCycleTime(process, &cycles);
        assert(ercd != 0);

        LONGLONG delta = cycles - start_ticks;

        accum_ticks += delta;

        return ((double)accum_ticks / 1000000);
    }

	void reset()
    {
        ULONG64 cycles = 0;
        BOOL ercd = QueryProcessCycleTime(process, &cycles);
        assert(ercd != 0);

        start_ticks = cycles;
    }

};

#endif