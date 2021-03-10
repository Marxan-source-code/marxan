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

    void initialiseHillClimbing(int puno, int spno, const vector<spustuff>& pu, const vector<sconnections>& connections, 
        vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, double cm,  int aggexist,
        vector<int>& R, double prop, int clumptype, int irun, stringstream& logBuffer, rng_engine& rngEngine);

    //move if we found a step (flip of pu unit) which generates better solution
    void hill_climbing(int puno, int spno, const vector<spustuff>& pu, const vector<sconnections>& connections,
        vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, vector<int>& R, double cm,
        scost& reserve, double costthresh, double tpf1, double tpf2,
        int clumptype,  int irun, long long iterations, string savename, stringstream& logBuffer, rng_engine& rngEngine);
    
    //consider two steps in a time, move if we found two steps which generates better solution
    void hill_climbing_two_steps(int puno, int spno, const vector<spustuff>& pu, const vector<sconnections>& connections,
        vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, vector<int>& R, double cm,
        scost& reserve, double costthresh, double tpf1, double tpf2,
        int clumptype,  int irun, long long iterations, string savename, stringstream& logBuffer, rng_engine& rngEngine);
} // namespace marxan