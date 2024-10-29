
#ifndef CPUTIME_H
#define CPUTIME_H

#include <Windows.h>
#include <cassert>

struct Cputime {

    double frequency;

    long long start_ticks;
    long long accum_ticks;

	Cputime() {
		reset();
        accum_ticks = 0;
	}

	double elapsed_ms() const
	{
        FILETIME creation, exit, kernel = {0,0}, user = {0, 0};
        BOOL ercd = GetProcessTimes(GetCurrentProcess(), &creation, &exit, &kernel, &user);
        assert(ercd != 0);

        INT64 kval = ((UINT64)kernel.dwHighDateTime << 32) + kernel.dwLowDateTime;
        INT64 uval = ((UINT64)user.dwHighDateTime << 32) + user.dwLowDateTime;

        kval /= 10;
        uval /= 10;


		LONGLONG delta = kval + uval - start_ticks;

		return ((double)delta / 1000);
	}

    double accum()
    {
        FILETIME creation, exit, kernel = {0,0}, user = {0, 0};
        BOOL ercd = GetProcessTimes(GetCurrentProcess(), &creation, &exit, &kernel, &user);
        assert(ercd != 0);

        INT64 kval = ((UINT64)kernel.dwHighDateTime << 32) + kernel.dwLowDateTime;
        INT64 uval = ((UINT64)user.dwHighDateTime << 32) + user.dwLowDateTime;

        kval /= 10;
        uval /= 10;


        LONGLONG delta = kval + uval - start_ticks;

        accum_ticks += delta;

        return ((double)accum_ticks / 1000);
    }

	void reset()
    {
        FILETIME creation, exit, kernel = {0,0}, user = {0, 0};
        BOOL ercd = GetProcessTimes(GetCurrentProcess(), &creation, &exit, &kernel, &user);
        assert(ercd != 0);

        INT64 kval = ((UINT64)kernel.dwHighDateTime << 32) + kernel.dwLowDateTime;
        INT64 uval = ((UINT64)user.dwHighDateTime << 32) + user.dwLowDateTime;

        kval /= 10;
        uval /= 10;

	    start_ticks = kval + uval;
    }

};

#endif