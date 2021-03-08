#pragma once

#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "connections.hpp"
#include "species.hpp"
#include "spu.hpp"
#include "anneal.hpp"
#include "options.hpp"

namespace marxan {
    using namespace std;

    typedef struct sfname
    {
        string inputdir;
        string outputdir;
        string specname;
        string puname;
        string puvsprname;
        string matrixspordername;
        string connectionname;
        string blockdefname;
        string bestfieldname;
        string connectionfilesname;
        string rbinarypath;
        string penaltyname;
        int savebest;
        int saverun;
        int savesum;
        int savesen;
        int savespecies;
        int savesumsoln;
        int savelog;
        int savesnapsteps;
        int savesnapchanges;
        int savesnapfrequency;
        int savepenalty;
        int savetotalareas;
        int savedebugtracefile;
        int savesolutionsmatrix;
        int solutionsmatrixheaders;
        int saveannealingtrace;
        int annealingtracerows;
        int saveitimptrace;
        int itimptracerows;
        int saverichness;
        int savespec;
        int savepu;
        int savepuvspr;
        int savematrixsporder;
        int rexecutescript;
        int rimagetype;
        int rimagewidth;
        int rimageheight;
        int rimagefontsize;
        int rclustercount;
    } sfname;


    // input reading
    int readConnections(int puno, vector<sconnections>& connections, const vector<spustuff>& pu,
        const map<int, int>& PULookup, const sfname& fnames);
    void readInputOptions(double& cm, double& prop, sanneal& anneal,
        int& iseed,
        long int& repeats, string& savename, sfname& fnames, string filename,
        srunoptions& runoptions, double& misslevel, int& heurotype, int& clumptype,
        int& itimptype, int& verb,
        double& costthresh, double& tpf1, double& tpf2);
    void readPenalties(vector<sspecies>& spec, int spno, sfname& fnames, map<int, int>& SPLookup);
    int readPlanningUnits(int& puno, vector<spustuff>& pu, const sfname& fnames);
    void readSparseMatrix(vector<spu>& SM, int puno, int spno, vector<spustuff>& pu,
        const map<int, int>& PULookup, const map<int, int>& SPLookup, const sfname& fnames);
    void readSparseMatrixSpOrder(vector<spusporder>& SM, int puno, int spno,
        const map<int, int>& PULookup, const map<int, int>& SPLookup, vector<sspecies>& spec, const sfname& fnames);

    int readSpecies(int& spno, vector<sspecies>& spec, const sfname& fnames);
    int readSpeciesBlockDefinition(int& gspno, vector<sgenspec>& gspec, sfname& fnames);


} // namespace marxan