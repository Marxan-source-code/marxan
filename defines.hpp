#pragma once

#define DEBUGTRACEFILE
#define PROB2D

// Flags to change to "define" for debugging purposes
#undef EXTRADEBUGTRACE
#undef ANNEALING_TEST
#undef DEBUGCHANGEPEN
#undef DEBUGCALCPENALTIES
#undef DEBUG_PRINTRESVALPROB
#undef DEBUG_COUNTMISSING
#undef DEBUG_HEURISTICS
#undef DEBUG_IIMPHEUR
#undef DEBUG_CLUSTERANALYSIS
#undef DEBUG_CONNECTIONCOST
#undef DEBUG_RESERVECOST
#undef DEBUGCHECKCHANGE
#undef DEBUG_CONNECTIONCOST2
#undef CREATE_R_SCRIPTS
#undef DEBUG_QA
#undef DEBUG_PROB1D
#undef DEBUG_PROB2D


#include <cstdarg>
#include <chrono>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "species.hpp"
#include "connections.hpp"
#include "spu.hpp"
#include "anneal.hpp"
#include "input.hpp"



namespace marxan {
    using namespace std;

    // For printing
    extern string sVersionString;
    extern string sMarxanWebSite;

    // Initialization constants
    extern int verbosity;
    extern FILE* fsavelog;
    extern int savelog;
    extern int asymmetricconnectivity;
    extern int fProb2D, fProb1D, fUserPenalties;
    extern int fOptimiseConnectivityIn, fSpecPROPLoaded;
    extern double delta;

    // some thresholds
    extern double rStartDecThresh, rEndDecThresh, rStartDecMult, rEndDecMult;

    // File names
    extern string sTraceFileName;
    extern string savelogname;

   
    // probability
    extern int iProbFieldPresent;
    extern double rProbabilityWeighting;

    // type definitions for sparse matrix optimisations data structures
    extern map<int, int> PULookup;
    extern map<int, int> SPLookup;

    extern vector<spu> SMGlobal;

    extern vector<spusporder> SMsporder;


    extern vector<spustuff> pu;

    extern scost debugcost_global;


    extern vector<sspecies> specGlobal, bestSpec;

    extern vector<sneighbour> debugnbr;

    extern vector<sconnections> connections;

    extern sanneal anneal_global;

    extern sfname fnames;

    //extern srunoptions runoptions;

}