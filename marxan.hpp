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

    void setBlockDefinitions(int gspno, int spno, int puno, vector<sgenspec>& gspec, vector<sspecies>& spec, vector<spustuff>& PU, vector<spu>& SM);
    void setDefaultRunOptions(int runopts, srunoptions& runoptions);
    int computePenalties(int puno, int spno, const vector<spustuff>& pu, vector<sspecies>& spec,
        const vector<sconnections>& connections, const vector<spu>& SM, vector<spu_out>& SM_out, vector<int>& PUtemp, int aggexist, double cm, int clumptype, rng_engine& rngEngine);
    int computePenaltiesOptimise(int puno, int spno, vector<spustuff>& pu, vector<sspecies>& spec,
        vector<sconnections>& connections, vector<spu>& SM, vector<spu_out>& SM_out, vector<spusporder>& SMsp,
        vector<int>& PUtemp, int aggexist, double cm, int clumptype, rng_engine& rngEngine);

    //double ConnectionCost2(sconnections &connection, vector<int> &R, int imode, int imode2, double cm);
    void computeReserveValue(int puno, int spno, const vector<int>& R, const vector<spustuff>& pu,
        const vector<sconnections>& connections, const vector<spu>& SM, vector<spu_out>& SM_out,
        double cm, vector<sspecies>& spec, int aggexist, scost& reserve, int clumptype, stringstream& logBuffer);
    void computeChangeScore(int iIteration, int ipu, int spno, int puno, const vector<spustuff>& pu, const vector<sconnections>& connections,
        vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, const vector<int>& R, double cm, int imode,
        scost& change, scost& reserve, double costthresh, double tpf1, double tpf2,
        double timeprop, int clumptype);
    void doChange(int ipu, int puno, vector<int>& R, scost& reserve, scost& change,
        const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out, vector<sspecies>& spec, const vector<sconnections>& connections,
        int imode, int clumptype, stringstream& logBuffer);


    void initialiseConnollyAnnealing(int puno, int spno, const vector<spustuff>& pu, const vector<sconnections>& connections, vector<sspecies>& spec,
        const vector<spu>& SM, vector<spu_out>& SM_out, double cm, sanneal& anneal, int aggexist,
        vector<int>& R, double prop, int clumptype, int irun, stringstream& logBuffer, rng_engine& rngEngine);
    void initialiseAdaptiveAnnealing(int puno, int spno, double prop, vector<int>& R, const vector<spustuff>& pu, const vector<sconnections>& connections,
        const vector<spu>& SM, vector<spu_out>& SM_out, const double cm, vector<sspecies>& spec, int aggexist, sanneal& anneal, int clumptype, stringstream& logBuffer, rng_engine& rngEngine);
    void thermalAnnealing(int spno, int puno, const vector<sconnections>& connections, vector<int>& R, double cm,
        vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out, scost& reserve,
        long int repeats, int irun, string savename, double misslevel,
        int aggexist, double costthresh, double tpf1, double tpf2, int clumptype, sanneal& anneal, stringstream& logBuffer, rng_engine& rngEngine);
    void quantumAnnealing(int spno, int puno, const vector<sconnections>& connections, vector<int>& R, double cm,
        vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out, scost& change, scost& reserve,
        long int repeats, int irun, string savename, double misslevel,
        int aggexist, double costthresh, double tpf1, double tpf2, int clumptype, sanneal& anneal, rng_engine& rngEngine);

    void secondaryExit(void);

    void iterativeImprovement(int puno, int spno, const vector<spustuff>& pu, const vector<sconnections>& connections,
        vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, vector<int>& R, double cm,
        scost& reserve, scost& change, double costthresh, double tpf1, double tpf2,
        int clumptype, int irun, string savename, stringstream& logBuffer, rng_engine& rngEngine);

} // namespace marxan