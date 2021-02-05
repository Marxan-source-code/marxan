#pragma once

#include <sstream>
#include <string>
#include <vector>
#include "algorithms.hpp"


namespace marxan {
    using namespace std;

    void Heuristics(int spno, int puno, const vector<spustuff>& pu, const vector<sconnections>& connections,
        vector<int>& R, double cm, vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, scost& reserve,
        double costthresh, double tpf1, double tpf2, int imode, int clumptype, stringstream& logBuffer, rng_engine& rngEngine);

} // namespace marxan