#pragma once

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
#include "options.hpp"
#include "algorithms.hpp"


namespace marxan {
    using namespace std;

    void initialiseConnollyAnnealing(int puno, int spno, const vector<spustuff>& pu, const vector<sconnections>& connections, vector<sspecies>& spec,
        const vector<spu>& SM, vector<spu_out>& SM_out, double cm, sanneal& anneal, int aggexist,
        vector<int>& R, double prop, int clumptype, int irun, stringstream& logBuffer, rng_engine& rngEngine);
    void initialiseAdaptiveAnnealing(int puno, int spno, double prop, vector<int>& R, const vector<spustuff>& pu, const vector<sconnections>& connections,
        const vector<spu>& SM, vector<spu_out>& SM_out, const double cm, vector<sspecies>& spec, int aggexist, sanneal& anneal, int clumptype, stringstream& logBuffer, rng_engine& rngEngine);
    void thermalAnnealing(int spno, int puno, const vector<sconnections>& connections, vector<int>& R, double cm,
        vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out, scost& reserve,
        long int repeats, int irun, string savename, double misslevel,
        int aggexist, double costthresh, double tpf1, double tpf2, int clumptype, sanneal& anneal, stringstream& logBuffer, rng_engine& rngEngine);
    
} // namespace marxan