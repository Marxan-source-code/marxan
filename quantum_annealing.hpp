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

    void quantumAnnealing(int spno, int puno, const vector<sconnections>& connections, vector<int>& R, double cm,
        vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out, scost& change, scost& reserve,
        long int repeats, int irun, string savename, double misslevel,
        int aggexist, double costthresh, double tpf1, double tpf2, int clumptype, sanneal& anneal, rng_engine& rngEngine);

} // namespace marxan