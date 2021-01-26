#pragma once

#include <csetjmp>
#include <cstdarg>
#include <chrono>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// declare structures, types, and some variables of type
// predefine functions that are called before they're defined
// Contains structures and functions that are used in more than 1 file.

#define DebugFree(x)

namespace marxan {
    using namespace std;

    // For printing
    extern string sVersionString;
    extern string sIanBallEmail;
    extern string sHughPossinghamEmail;
    extern string sMattWattsEmail;
    extern string sMarxanWebSite;

    // Initialization constants
    extern jmp_buf jmpbuf;
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

    // rng
    extern mt19937 rngEngine;

    // probability
    extern int iProbFieldPresent;
    extern double rProbabilityWeighting;

    // type definitions for sparse matrix optimisations data structures
    extern map<int, int> PULookup;
    extern map<int, int> SPLookup;

    typedef struct spu
    {
        double amount; // amount of species in pu
        double prob; // optional field that does ???
        int spindex; // index/id of species
    } spu;

    //Separate output fields for multithreading 
    typedef struct spu_out
    {
        int clump;
        spu_out() { clump = 0; }
    } spu_out;


    extern vector<spu> SMGlobal;

    typedef struct spusporder
    {
        double amount;
        int puindex;
    } spusporder;

    extern vector<spusporder> SMsporder;

    // type definitions for original Ian Ball data structures

    typedef struct spustuff
    {
        int id;
        int status;
        double xloc,yloc;
        double cost;
        double prob;
        int richness,offset,probrichness,proboffset;
    } spustuff;

    extern vector<spustuff> pu;

    typedef struct scost
    {
        double total;
        int pus;
        double connection;
        int missing;
        double penalty;
        double cost;
        double threshpen;
        double shortfall;
        double probability1D;
        double probability2D;
    } scost;

    extern scost debugcost_global;

    typedef struct sclumps
    {
      int clumpid;
      double amount;
      int occs;
      vector<int> head;
    } sclumps; /* Clump nodes for species Clump Structure */

    typedef struct sspecies
    {
        int name;
        int type;
        string sname;
        double target;
        double prop;
        int targetocc;
        double spf;
        double penalty;
        double rUserPenalty;
        double amount;
        double expected1D, expected2D, variance1D, variance2D;
        int occurrence;
        double sepdistance;
        int sepnum;
        int separation;
        int clumps;
        double target2;  // Only clumping species need this
        vector<sclumps> head;  // needed for clumping species
        int richness,offset;
        double Zscore1D, Zscore2D;
        double probability1D, probability2D;
        double ptarget1d, ptarget2d;
    }sspecies;

    extern vector<sspecies> specGlobal, bestSpec;

    /* Connectivity Structure. Fixed connectivity number.*/
    typedef struct sneighbour
    {
        int nbr; // puid of neighbour
        double cost;
        int connectionorigon; // for asymmetric connections, whether connection starts from this node.
        sneighbour(int nbr, double cost, int connectionorigon) { this->nbr = nbr; this->cost = cost; this->connectionorigon = connectionorigon; }
    } sneighbour;
    
    extern vector<sneighbour> debugnbr;

    typedef struct sconnections
    {
        vector<sneighbour> first;
        double fixedcost;
        int nbrno;
    } sconnections;

    extern vector<sconnections> connections;

    typedef struct sanneal
    {
        long int Titns;
        long int iterations;
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

    extern sanneal anneal_global;

    typedef struct sgenspec
    {
        int type;
        int targetocc;
        double target;
        double target2;
        int sepnum;
        double sepdistance;
        double prop;
        double spf;
    } sgenspec;

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

    extern sfname fnames;

    // thermal annealling = simulated annealling
    // Penalties only needed for sa?? - double check
    typedef struct srunoptions
    {
        int CalcPenaltiesOn;
        int HeuristicOn;
        int ThermalAnnealingOn;
        int QuantumAnnealingOn;
        int ItImpOn;
    } srunoptions;

    extern srunoptions runoptions;

void setBlockDefinitions(int gspno,int spno,int puno, vector<sgenspec> &gspec, vector<sspecies> &spec, vector<spustuff> &PU, vector<spu> &SM);
void setDefaultRunOptions(int runopts, srunoptions &runoptions);
int computePenalties(int puno, int spno, const vector<spustuff>& pu, vector<sspecies>& spec,
    const vector<sconnections>& connections, const vector<spu>& SM, vector<spu_out>& SM_out, vector<int>& PUtemp, int aggexist, double cm, int clumptype);
int computePenaltiesOptimise(int puno, int spno, vector<spustuff>& pu, vector<sspecies>& spec,
    vector<sconnections>& connections, vector<spu>& SM, vector<spu_out>& SM_out, vector<spusporder>& SMsp,
    vector<int>& PUtemp, int aggexist, double cm, int clumptype);

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

stringstream displayValueForPUs(int puno, int spno, const vector<int>& R, const scost& reserve,
    const vector<sspecies>& spec, double misslevel);
//void TimePassed(void);

#ifndef annealingheader
#define annealingheader

void initialiseConnollyAnnealing(int puno, int spno, const vector<spustuff>& pu, const vector<sconnections>& connections, vector<sspecies>& spec,
    const vector<spu>& SM, vector<spu_out>& SM_out, double cm, sanneal& anneal, int aggexist,
    vector<int>& R, double prop, int clumptype, int irun, stringstream& logBuffer);
void initialiseAdaptiveAnnealing(int puno, int spno, double prop, vector<int>& R, const vector<spustuff>& pu, const vector<sconnections>& connections,
    const vector<spu>& SM, vector<spu_out>& SM_out, const double cm, vector<sspecies>& spec, int aggexist, sanneal& anneal, int clumptype, stringstream& logBuffer);
void thermalAnnealing(int spno, int puno, const vector<sconnections>& connections, vector<int>& R, double cm,
    vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out, scost& reserve,
    long int repeats, int irun, string savename, double misslevel,
    int aggexist, double costthresh, double tpf1, double tpf2, int clumptype, sanneal& anneal, stringstream& logBuffer);
void quantumAnnealing(int spno, int puno, const vector<sconnections>& connections, vector<int>& R, double cm,
    vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out, scost& change, scost& reserve,
    long int repeats, int irun, string savename, double misslevel,
    int aggexist, double costthresh, double tpf1, double tpf2, int clumptype, sanneal & anneal);
#endif


void ClearClumps(int spno, vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out);

//int NameToPUID(int puno,int name, vector<spustuff> &pu);
//int NameToSPID(int spno,int name,vector<sspecies> &spec);

//void applyUserPenalties(vector<sspecies> &spec,int spno);

//void writeSparseMatrix(int iSMno,int puno, vector<spustuff> &pu, vector<sspecies> &spec, vector<spu> &SM,sfname fnames);
void appendTraceFile(string sMess,...);

void secondaryExit(void);
//void SaveSeed(int iseed);

void iterativeImprovement(int puno, int spno, const vector<spustuff>& pu, const vector<sconnections>& connections,
    vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, vector<int>& R, double cm,
    scost& reserve, scost& change, double costthresh, double tpf1, double tpf2,
    int clumptype, int irun, string savename, stringstream& logBuffer);

// input reading
int readConnections(int& puno, vector<sconnections>& connections, vector<spustuff>& pu,
                    map<int,int>& PULookup, sfname& fnames);
void readInputOptions(double &cm,double &prop,sanneal &anneal,
                      int &iseed,
                      long int &repeats, string& savename, sfname &fname, string filename,
                      int &runopts,double &misslevel,int &heurotype,int &clumptype,
                      int &itimptype, int &verb,
                      double &costthresh,double &tpf1,double &tpf2);
                      void readPenalties(vector<sspecies> &spec,int spno,sfname& fnames,map<int,int> &SPLookup);
int readPlanningUnits(int& puno, vector<spustuff>& pu, sfname& fnames);
void readSparseMatrix(int &iSMSize, vector<spu> &SM, int puno, int spno, vector<spustuff> &pu,
                      map<int,int> &PULookup,map<int,int> &SPLookup,
                      sfname& fnames);
void readSparseMatrixSpOrder(int &iSMSize, vector<spusporder> &SM, int puno, int spno,
                             map<int,int> &PULookup,map<int,int> &SPLookup, vector<sspecies> &spec,
                             sfname &fnames);
int readSpecies(int &spno, vector<sspecies>& spec, sfname& fnames);
int readSpeciesBlockDefinition(int& gspno, vector<sgenspec>& gspec, sfname& fnames);

//void SetVerbosity(int verb);
void displayTimePassed(chrono::high_resolution_clock::time_point start);

//void SetLogFile(int my_savelog, string my_savelogname);

void displayErrorMessage(string sMess,...);
void displayWarningMessage(string sMess,...);
void displayProgress(string sMess,...);
void displayProgress1(string sMess,...);
void displayProgress2(string sMess,...);
void displayProgress3(string sMess,...);
void displayShutdownMessage(chrono::high_resolution_clock::time_point start);
void displayUsage(string programName);

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
void writeR(int iMessage, string sMessage,int puno, const vector<int>& R, const vector<spustuff>& pu, const sfname& fnames);
void writeRichness(int puno, const vector<spustuff>& pu, string savename, int iOutputType);
void writeScenario(int puno,int spno,double prop,double cm,
                    sanneal& anneal,int seedinit,long int repeats,int clumptype,
                    int runopts,int heurotype,double costthresh, double tpf1, double tpf2,
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

typedef struct sseplist{
    int size;
    vector<int> head;
} sseplist;

// clumping
void AddNewPU(int ipu, int isp, const vector<sconnections>& connections, vector<sspecies>& spec, const vector<spustuff>& pu,
    const vector<spu>& SM, vector<spu_out>& SM_out, int clumptype);
void RemPu(int ipu, int isp, const vector<sconnections>& connections, vector<sspecies>& spec, const vector<spustuff>& pu,
    const vector<spu>& SM, vector<spu_out>& SM_out, int clumptype);
int CalcPenaltyType4(int isp, int puno, const vector<spu>& SM, vector<spu_out>& SM_out, const vector<sconnections>& connections,
    vector<sspecies>& spec, const vector<spustuff>& pu, double cm, int clumptype);
int CountSeparation(int isp, const vector<sclumps>& newno,
    const vector<spustuff>& pu, const vector<spu>& SM, const vector<spu_out>& SM_out, const vector<sspecies>& spec, int imode);
int CountSeparation2(int isp, int ipu, const vector<sclumps>& newno, int puno, const vector<int>& R,
    const vector<spustuff>& pu, const vector<spu>& SM, const vector<spu_out>& SM_out, const vector<sspecies>& spec, int imode);
vector<int> makelist(int isp, int ipu, int puno, const vector<int>& R, const vector<sclumps>& newno, const vector<sspecies>& spec,
    const vector<spustuff>& pu, const vector<spu>& SM, const vector<spu_out>& SM_out, int imode);
double NewPenalty4(int ipu, int isp, int puno, const vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, const vector<spu_out>& SM_out,
    const vector<int>& R, const vector<sconnections>& connections, int imode, int clumptype);
double PartialPen4(int isp, double amount, const vector<sspecies>& spec, int clumptype);
int SepDealList(const vector<int>& head, vector<sseplist>& Dist, const vector<spustuff>& pu,
    const vector<sspecies>& spec, int first, int sepnum, double targetdist, int isp);
void SetSpeciesClumps(int puno, const vector<int>& R, vector<sspecies>& spec, const vector<spustuff>& pu,
    const vector<spu>& SM, vector<spu_out>& SM_out, const vector<sconnections>& connections, int clumptype);
void SpeciesAmounts(int spno, int puno, vector<sspecies>& spec, const vector<spustuff>& pu, vector<spu>& SM,
    vector<int>& R, int clumptype);
void SpeciesAmounts4(int isp,vector<sspecies> &spec,int clumptype);

void Heuristics(int spno, int puno, const vector<spustuff>& pu, const vector<sconnections>& connections,
    vector<int>& R, double cm, vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, scost& reserve,
    double costthresh, double tpf1, double tpf2, int imode, int clumptype, stringstream& logBuffer);

double ComputeProbability1D(const vector<double>& ExpectedAmount1D, const vector<double>& VarianceInExpectedAmount1D,
    int spno, vector<sspecies>& spec);
double ComputeProbability2D(vector<double>& ExpectedAmount2D, vector<double>& VarianceInExpectedAmount2D,
    int spno, vector<sspecies>& spec);
double ChangeProbability1D(int iIteration, int ipu, int spno, int puno, vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, int imode);
double ChangeProbability2D(int iIteration, int ipu, int spno, int puno, vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, int imode);
void ComputeP_AllPUsSelected_1D(const string& savename, int puno, int spno, const vector<spustuff>& pu, const vector<spu>& SM, const vector<sspecies>& spec);
void ComputeP_AllPUsSelected_2D(const string& savename, int puno, int spno, const vector<spustuff>& pu, const vector<spu>& SM, const vector<sspecies>& spec);
void ReturnProbabilityAmounts1D(vector<double>& ExpectedAmount1D, vector<double>& VarianceInExpectedAmount1D, int ipu,
    int puno, const vector<spustuff>& pu, const vector<spu>& SM);
void ReturnProbabilityAmounts2D(vector<double>& ExpectedAmount2D, vector<double>& VarianceInExpectedAmount2D, int ipu,
    int puno, const vector<spustuff>& pu, const vector<spu>& SM);

} // namespace marxan