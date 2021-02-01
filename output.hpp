#pragma once

#include <chrono>
#include <map>
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

    void appendTraceFile(string sMess, ...);

    void displayTimePassed(chrono::high_resolution_clock::time_point start);

    //void SetLogFile(int my_savelog, string my_savelogname);

    void displayErrorMessage(string sMess, ...);
    void displayWarningMessage(string sMess, ...);
    void displayProgress(string sMess, ...);
    void displayProgress1(string sMess, ...);
    void displayProgress2(string sMess, ...);
    void displayProgress3(string sMess, ...);
    void displayShutdownMessage(chrono::high_resolution_clock::time_point start);
    void displayUsage(string programName);
    stringstream displayValueForPUs(int puno, int spno, const vector<int>& R, const scost& reserve,
        const vector<sspecies>& spec, double misslevel);


    void appendSolutionsMatrix(int iRun, int puno, const vector<int>& R, string savename, int iOutputType, int iIncludeHeaders);
    string computeSummary(int puno, int spno, const vector<int>& R, const vector<sspecies>& spec, const scost& reserve,
        int itn, double misslevel, int imode);
    void computeTotalAreas(int puno, int spno, const vector<spustuff>& pu, const vector<sspecies>& spec, const vector<spu>& SM);
    void createLogFile(int my_savelog, string my_savelogname);
    void createSolutionsMatrix(int puno, const vector<spustuff>& pu, string savename_ism, int iOutputType, int iIncludeHeaders);
    void createTraceFile(void);
    void displayStartupMessage(void);
    void writeAsymmetricConnectionFile(int puno, const vector<sconnections>& connections, const vector<spustuff>& pu, sfname fnames);
    void writePenalty(int spno, const vector<sspecies>& spec, string savename, int iOutputType);
    void writePenaltyPlanningUnits(int puno, const vector<spustuff>& pu, const vector<int>& Rtemp, string savename, int iOutputType);
    void writePu(int puno, const vector<spustuff>& pu, string savename);
    void writeR(int iMessage, string sMessage, int puno, const vector<int>& R, const vector<spustuff>& pu, const sfname& fnames);
    void writeRichness(int puno, const vector<spustuff>& pu, string savename, int iOutputType);
    void writeScenario(int puno, int spno, double prop, double cm,
        sanneal& anneal, int seedinit, long int repeats, int clumptype,
        int runopts, int heurotype, double costthresh, double tpf1, double tpf2,
        string savename);
    void writeSecondarySyncFile(void);
    void writeSecondarySyncFileRun(int iSyncRun);
    void writeSolution(int puno, const vector<int>& R, const vector<spustuff>& pu, string savename, int imode, const sfname& fnames);
    void writeSpec(int spno, const vector<sspecies>& spec, string savename);
    void writeSpecies(int spno, vector<sspecies>& spec, string savename, int imode, double misslevel);
    void writeSummary(string savename, const vector<string>& summaries, int imode);
    void writeSumSoln(int puno, const vector<int>& sumsoln, const vector<spustuff>& pu, string savename, int imode);
    void writeTotalAreas(int puno, int spno, const vector<spustuff>& pu, const vector<sspecies>& spec, const vector<spu>& SM, string savename, int iOutputType);
    void writeWeightedConnectivityFile(const sfname& fnames);

} // namespace marxan