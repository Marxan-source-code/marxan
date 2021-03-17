#pragma once

namespace marxan {

    typedef struct sanneal
    {
        long int Titns;
        long long iterations;
        long int Tlen;
        double Tinit;    /* Initial Temperature */
        double Tcool;    /* Cooling Factor */
        double temp; /* Current Temperature */
        double tempold;
        int type;    /* Type of annealing. 0 = none, 1 = fixed, 2 = adaptive */
        double sigma; /*Used in adaptive annealing */
        double sum; /* Used in adaptive annealing */
        double sum2; /* used in adaptive annealing */
    } sanneal; /* Annealing Control handler */


} // namespace marxan