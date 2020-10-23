#pragma once

#include <random>
#include <string>
#include <vector>

// declare structures, types, and some variables of type
// predefine functions that are called before they're defined

#define DebugFree(x)
#ifndef mainheaderfile
#define mainheaderfile

namespace marxan {
    using namespace std;

    int asymmetricconnectivity;
    int fOptimiseConnectivityIn;
    double delta;

    // rng
    mt19937 rngEngine;

    // type definitions for sparse matrix optimisations data structures
    typedef struct binsearch
    {
        int name;
        int index;
    } binsearch;

    vector<binsearch> PULookup;
    vector<binsearch> SPLookup;

    typedef struct spu
    {
        double amount; // amount of species in pu
        double prob; // optional field that does ???
        vector<int> clump; // array of clump terms by thread
        int spindex; // index/id of species
    } spu;

    vector<spu> SMGlobal;

    typedef struct spusporder
    {
        double amount;
        int puindex;
    } spusporder;

    vector<spusporder> SMsporder;

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

    vector<spustuff> pu;

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

    scost debugcost_global;

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

    vector<sspecies> specGlobal, bestSpec;

    /* Connectivity Structure. Fixed connectivity number. Should replace with link list! */
    typedef struct sneighbour
    {
        int nbr;
        double cost;
        int connectionorigon;
    } sneighbour;
    
    vector<sneighbour> debugnbr;

    typedef struct sconnections
    {
        vector<sneighbour> first;
        double fixedcost;
        int nbrno;
    } sconnections;

    vector<sconnections> connections;

    typedef struct sclumps
    {
      int clumpid;
      double amount;
      int occs;
      vector<int> head;
    } sclumps; /* Clump nodes for species Clump Structure */

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

    sanneal anneal_global;

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

    sfname fnames;

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

    srunoptions runoptions;

    /*
    struct slink
    {
        int id;
        struct slink *next;
    };
    */

    typedef struct iimp
    {
        double randomfloat;
        int puindex;
    } iimp;

void addReserve(int puno,vector<spustuff> &pu, vector<int> &R);
void computeSpecProp(int spno, vector<sspecies> &spec, int puno, vector<spustuff> &pu, vector<spu> &SM);
void setBlockDefinitions(int gspno,int spno,int puno, vector<sgenspec> &gspec, vector<sspecies> &spec, vector<spustuff> &PU, vector<spu> &SM);
void setDefaultTargets(int spno, vector<sspecies> &spec);
void setDefaultRunOptions(int runopts, srunoptions &runoptions);
int computePenalties(int puno,int spno, vector<spustuff> &pu, vector<sspecies> &spec,
                     vector<sconnections> &connections, vector<spu> &SM, vector<int> &PUtemp, int aggexist, double cm, int clumptype);
int computePenaltiesOptimise(int puno,int spno, vector<spustuff> &pu, vector<sspecies> &spec,
                             vector<sconnections> &connections, vector<spu> &SM, vector<spusporder> &SMsp,
                             vector<int> &PUtemp, int aggexist, double cm, int clumptype);

double computePlanningUnitValue(int ipu, vector<spustuff> &pu, vector<sconnections> &connections, double cm);
double ConnectionCost2(int ipu, vector<sconnections> &connections, vector<int> &R, int imode, int imode2, double cm);
void computeConnectivityIndices(double &rConnectivityTotal, double &rConnectivityIn,
                                double &rConnectivityEdge, double &rConnectivityOut,
                                int puno, vector<int> &R, vector<sconnections> &connections);
void computeReserveValue(int puno,int spno, vector<int> &R, vector<spustuff> &pu,
                         vector<sconnections> &connections, vector<spu> &SM,
                         double cm, vector<sspecies> &spec, int aggexist, scost &reserve,int clumptype);
void initialiseReserve(int puno,double prop, vector<int> &R);
void computeChangeScore(int iIteration,int ipu,int spno,int puno,vector<spustuff> &pu, vector<sconnections> &connections,
                        vector<sspecies> &spec, vector<spu> &SM, vector<int> &R, double cm, int imode,
                        scost &change, scost &reserve,double costthresh,double tpf1, double tpf2,
                        double timeprop,int clumptype);
double computeSpeciesPlanningUnitPenalty(int ipu,int isp,vector<sspecies> &spec,vector<spustuff> &pu, vector<spu> &SM,int imode);
void doChange(int ipu,int puno,vector<int> &R, scost &reserve, scost &change,
              vector<spustuff> &pu,vector<spu> &SM,vector<sspecies> &spec,vector<sconnections> &connections,
              int imode,int clumptype);

int computeRepresentationMISSLEVEL(int spno,vector<sspecies> &spec,double misslevel,double &shortfall,double &rMinimumProportionMet);
void displayValueForPUs(int puno, int spno,vector<int> &R, scost &reserve,
                        vector<sspecies> &spec,double misslevel);
void TimePassed(void);
void PauseProg(void);
void PauseExit(void);

#endif
#ifndef annealingheader
#define annealingheader


void initialiseConnollyAnnealing(int puno,int spno,vector<spustuff> &pu, vector<sconnections> &connections, vector<sspecies> &spec,
                                 vector<spu> &SM,double cm, sanneal &anneal,int aggexist,
                                 vector<int> &R,double prop,int clumptype,int irun);
void initialiseAdaptiveAnnealing(int puno,int spno,double prop,vector<int> &R,vector<spustuff> &pu,vector<sconnections> &connections,
                                 vector<spu> &SM,double cm,vector<sspecies> &spec,int aggexist,sanneal &anneal,int clumptype);
void thermalAnnealing(int spno, int puno, vector<sconnections> &connections,vector<int> &R, double cm,
                      sspecies &spec, vector<spustuff> &pu, vector<spu> &SM, scost &reserve,
                      long int repeats,int irun,string savename,double misslevel,
                      int aggexist,double costthresh, double tpf1, double tpf2,int clumptype, sanneal &anneal);
void quantumAnnealing(int spno, int puno, vector<sconnections> &connections,vector<int> &R, double cm,
                      sspecies &spec, vector<spustuff> &pu, vector<spu> &SM, scost &change, scost &reserve,
                      long int repeats,int irun,string savename,double misslevel,
                      int aggexist,double costthresh, double tpf1, double tpf2,int clumptype, sanneal &anneal);
#endif

#ifndef fileinheader
#define fileinheader

/*struct snlink{string name; struct snlink *next;};*/

/* Parmtypes */
#define NOTYPE 0
#define INTEGER 1
#define LONGINT 2
#define REAL 3
#define DOUBLE 4
#define STRING 5

vector<int> storeFieldName(string &varlist,int numvars,string sVarName,
    vector<int> &head,string fname);

int NameToPUID(int puno,int name, vector<spustuff> &pu);
int NameToSPID(int spno,int name,vector<sspecies> &spec);

void readInputOptions(double &cm,double &prop,sanneal &anneal,
                      int &iseed,
                      long int &Repeats,char savename[],sfname &fname,char filename[],
                      int &Runopts,double &misslevel,int &heurotype,int &clumptype,
                      int &itimptype, int &verb,
                      double &costthresh,double &tpf1,double &tpf2);

// unused but declared?
//int LoadPuDat(int &puno,spustuff &pu[],sfname fnames);
//int LoadSpecDat(int &spno,sspecies *spec[],sfname fnames);
//int ReadGenSpeciesData(int *gspno,sgenspec *gspec[],sfname fnames);
//int DumpAsymmetricConnectivityFile(int puno,vector<sconnections> &connections,vector<spustuff> &pu,sfname fnames);

int readConnections(int puno,vector<sconnections> &connections,vector<spustuff> &pu,
                    vector<binsearch> &PULookup,sfname &fnames);

//void ReadPUVSPFile22(int puno,int spno,vector<spu> &SM,vector<spustuff> &pu,
//    vector<sspecies> &spec,sfname fnames);
//void ReadPUVSPFileTable(FILE *infile, int puno,int spno,vector<spu> &SM,vector<spustuff> &pu,
//    vector<sspecies> &spec);

void readPenalties(vector<sspecies> &spec,int spno,sfname fnames,vector<binsearch> &SPLookup);
void applyUserPenalties(vector<sspecies> &spec,int spno);

// functions for Matt's Big O notation optimisation
void readSparseMatrix(int &iSMSize, vector<spu> &SM, int puno, int spno, vector<spustuff> &pu,
                      vector<binsearch> &PULookup,vector<binsearch> &SPLookup,
                      sfname fnames);
void writeSparseMatrix(int iSMno,int puno, vector<spustuff> &pu, vector<sspecies> &spec, vector<spu> &SM,sfname fnames);
void readSparseMatrixSpOrder(int &iSMSize, vector<spusporder> &SM, int puno, int spno,
                             vector<binsearch> &PULookup,vector<binsearch> &SPLookup,// vector<sspecies> &spec,
                             sfname &fnames);
void writeBinarySearchArrays(string &sName,sfname &fnames, int puno, int spno, vector<binsearch> &PULookup,
                             vector<binsearch> &SPLookup);

void computeBinarySearch(int puno, int spno, vector<spustuff> &pu, vector<sspecies> &spec,
                         vector<binsearch> &PULookup, vector<binsearch> &SPLookup);
int binarySearchPuIndex(int puno,int name, vector<binsearch> &PULookup);
int binarySearchSpecIndex(int spno,int name, vector<binsearch> &SPLookup);
int returnIndexSpecAtPu(vector<spustuff> &pu, vector<spu> &SM, int iPUIndex, int iSpecIndex);
double returnAmountSpecAtPu(vector<spustuff> &pu, vector<spu> &SM, int iPUIndex, int iSpecIndex);
void appendTraceFile(string sMess,...);

#endif

#ifndef spexioheader
#define spexioheader
void SaveSeed(int iseed);

#endif

#ifndef heuristicheader
#define heuristicheader


double GreedyPen(int ipu, int puno, int spno, vector<sspecies> &spec,vector<int> &R,vector<spustuff> &pu,
                 vector<spu> &SM,int clumptype);
double GreedyScore(int ipu,int puno,int spno, sspecies &spec,vector<spu> &SM,vector<sconnections> &connections,
                   vector<int> &R,vector<spustuff> &pu,double cm,int clumptype);
void SetRareness(int puno, int spno, vector<double> &Rare,vector<int> &R,vector<spustuff> &pu,vector<spu> &SM);
double RareScore(int isp,int ipu,int puno,vector<sspecies> &spec,vector<spu> &SM, vector<int> &R,
    vector<sconnections> &connections, vector<spustuff> &pu, double cm,int clumptype);
double MaxRareScore(int ipu,int puno,vector<sspecies> &spec,vector<spu> &SM,
    vector<int> &R,vector<sconnections> &connections,vector<spustuff> &pu,double cm,vector<double> &Rare,int clumptype);
double BestRareScore(int ipu,int puno,vector<sspecies> &spec,vector<spu> &SM,
    vector<int> &R,vector<sconnections> &connections,vector<spustuff> &pu,double cm,vector<double> &Rare,int clumptype);
double AveRareScore(int ipu,int puno,vector<sspecies> &spec,vector<spu> &SM,
    vector<int> &R,vector<sconnections> &connections,vector<spustuff> &pu,double cm,vector<double> &Rare,int clumptype);
double SumRareScore(int ipu,int puno,vector<sspecies> &spec,vector<spu> &SM,
    vector<int> &R,vector<sconnections> &connections,vector<spustuff> &pu,double cm,vector<double> &Rare,int clumptype);
void SetAbundance(int puno,vector<double> &Rare,vector<spustuff> &pu,vector<spu> &SM);
double Irreplaceability(int ipu,int isp, vector<double> &Rare,vector<spustuff> &pu,vector<spu> &SM,sspecies &spec);
double ProdIrr(int ipu,vector<double> &Rare,vector<spustuff> &pu,vector<spu> &SM,sspecies &spec);
double SumIrr(int ipu,vector<double> &Rare,vector<spustuff> &pu,vector<spu> &SM,sspecies &spec);
void Heuristics(int spno,int puno,vector<spustuff> &pu,vector<sconnections> &connections,
        vector<int> &R, double cm,sspecies &spec,vector<spu> &SM, scost &reserve,
        double costthresh, double tpf1,double tpf2, int imode,int clumptype);

#endif

#ifndef itimpheader
#define itimpheader

/*
slink* ItImpDiscard(int ichoice, slink *list, slink **discard);
slink* ItImpUndiscard(slink *list, slink **discard);
int FindSwap( slink **list,double targetval,int itestchoice,int puuntried,
             int puno,int spno,vector<spustuff> &pu, vector<sconnections> &connections,
             vector<sspecies> &spec,vector<spu> &SM,
             vector<int> &R, double cm, scost &reserve, scost &change,
             double costthresh, double tpf1, double tpf2, int clumptype);
*/

void iterativeImprovement(int puno,int spno,vector<spustuff> &pu, vector<sconnections> &connections,
                          vector<sspecies> &spec,vector<spu> &SM,vector<int> &R, double cm,
                          scost &reserve, scost &change,double costthresh,double tpf1, double tpf2,
                          int clumptype,int irun,string savename);

#endif

#ifndef randomheader
#define randomheader

float returnRandomFloat(void);
void initialiseRandomSeed(int iSeed);
int returnRandom (int num);


#endif

#ifndef outputheader
#define outputheader

void displayStartupMessage(void);
void displayShutdownMessage(void);

void SetVerbosity(int verb);

void displayErrorMessage(string sMess,...);
void displayWarningMessage(string sMess,...);

void displayProgress(string sMess,...);
void displayProgress1(string sMess,...);
void displayProgress2(string sMess,...);
void displayProgress3(string sMess,...);

void displayTimePassed(void);

void SetLogFile(int my_savelog, string my_savelogname);

#endif

typedef struct sseplist{
    int size;
    vector<int> head;
} sseplist;

double probZUT(double z);
double probZLT(double z);

void ComputeP_AllPUsSelected_1D(string &savename,int puno,int spno,vector<spustuff> &pu,vector<spu> &SM,vector<sspecies> &spec);
void ComputeP_AllPUsSelected_2D(string &savename,int puno,int spno,vector<spustuff> &pu,vector<spu> &SM,vector<sspecies> &spec);
double ChangeProbability1D(int iIteration, int ipu, int spno,int puno,vector<sspecies> &spec,vector<spustuff> &pu,vector<spu> &SM,int imode);
double ChangeProbability2D(int iIteration, int ipu, int spno,int puno,vector<sspecies> &spec,vector<spustuff> &pu,vector<spu> &SM,int imode);
double Probability(int ipu, int spno,int puno,vector<sspecies> &spec,vector<spustuff> &pu,vector<spu> &SM);
void ReturnProbabilityAmounts1D(vector<double> &ExpectedAmount1D, vector<double> &VarianceInExpectedAmount1D,int ipu,
                                int puno,vector<spustuff> &pu,vector<spu> &SM);
void ReturnProbabilityAmounts2D(vector<double> &ExpectedAmount2D,vector<double> &VarianceInExpectedAmount2D,int ipu,
                                int puno,vector<spustuff> &pu,vector<spu> &SM);
double ComputeProbability1D(vector<double> &ExpectedAmount1D, vector<double> &VarianceInExpectedAmount1D,
                            int spno,vector<sspecies> &spec);
double ComputeProbability2D(vector<double> &ExpectedAmount2D, vector<double> &VarianceInExpectedAmount2D,
                            int spno,vector<sspecies> &spec);
void createDebugFile(string &sFileName, string &sHeader, sfname &fnames);
void appendDebugFile(string &sFileName, string &sLine, sfname &fnames);
void writePenalty(int spno,vector<sspecies> &spec, string &savename, int iOutputType);
void writePenaltyPlanningUnits(int puno,vector<spustuff> &pu, vector<int> &Rtemp, string &savename,int iOutputType);
void writeSpec(int spno,vector<sspecies> &spec, string &savename);

void createSolutionsMatrix(int puno,vector<spustuff> &pu, string &savename_ism,int iOutputType,int iIncludeHeaders);

void writeSlaveSyncFileRun(int iSyncRun);
void slaveExit(void);

} // namespace marxan