
#ifndef CPUTIME_H
#define CPUTIME_H

struct Cputime {

    Cputime() {	}

    double elapsed_ms() const
    {
	    return ((double)1);
    }

    double accum()
    {
        return ((double)1);
    }

    void reset()
    {
    }

};

#endif
