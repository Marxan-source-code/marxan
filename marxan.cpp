// C++ code for Marxan
// version 2.3 introduced multiple connectivity files and their associated weighting file
// version 2.4.3 introduced 1D and 2D probability
// version 3.0.0 is refactoring of code in 2019

#define DEBUGTRACEFILE
#define PROB2D

// Flags to change to "define" for debugging purposes
#undef MEMDEBUG
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

#include <algorithm>
#include <chrono>
#include <ctime>
#include <cfloat>
#include <iostream>
#include <omp.h>
#include <chrono> 

// load the required function definition modules
#include "algorithms.hpp"
#include "computation.hpp"
#include "marxan.hpp"
#include "utils.hpp"

namespace marxan {
using namespace algorithms;
using namespace utils;

// Initialize global constants
jmp_buf jmpbuf;
double delta;
long int RandSeed1;
int iMemoryUsed=0;
int fSpecPROPLoaded = 0;
int iProbFieldPresent = 0;
int iOptimisationCalcPenalties = 1;
int savelog;
int verbosity = 0;
int asymmetricconnectivity = 0;
string sVersionString = "Marxan v 4.0.0 alpha";
string sIanBallEmail = "ian.ball@aad.gov.au";
string sHughPossinghamEmail = "hugh.possingham@tnc.org";
string sMattWattsEmail= "matt.watts@une.edu.au";
string sMarxanWebSite= "http://marxan.net";
string sTraceFileName;
string sApplicationPathName;
string savelogname;
FILE* fsavelog;

// init global structures
sanneal anneal_global;
vector<sconnections> connections;
sfname fnames;
vector<spu> SMGlobal;
vector<spusporder> SMsporder;
vector<spustuff> pu;
map<int, int> PULookup, SPLookup;
vector<sspecies> specGlobal, bestSpec;
mt19937 rngEngine;
srunoptions runoptions;
chrono::high_resolution_clock::time_point startTime;

double rProbabilityWeighting = 1;
double rStartDecThresh = 0.7, rEndDecThresh = 0.3, rStartDecMult = 3, rEndDecMult = 1;
double rQAPROP = 0.5, rQADECAY = 0.0001, rQADECAYB = 0, rQAACCPR = 0;
int iQADECAYTYPE = 0;

int fProb2D = 0, fProb1D = 0, fUserPenalties = 0;
int fOptimiseConnectivityIn = 0;

// number of the best solution
int bestRun = 1;
// score of the best solution
double bestScore;
// R vector of planning unit status for the best solution
vector<int> bestR;

// is marxan being called by another program?
int marxanIsSecondary = 0;

// runs the loop for each "solution" marxan is generating
void executeRunLoop(int iSparseMatrixFileLength, long int repeats,int puno,int spno,double cm,int aggexist,double prop,int clumptype,double misslevel,
                    string savename,double costthresh,double tpf1,double tpf2,int heurotype,int runopts,
                    int itimptype,vector<int> sumsoln)
{
    string bestRunString;
    vector<string> summaries(repeats); // stores individual summaries for each run
    bestR.resize(puno);
    bestScore = DBL_MAX;
    bestSpec.resize(spno); // best species config

    // Locks for bestR and bestScore as it will be read/written by multiple threads
    omp_lock_t bestR_write_lock;
    omp_init_lock(&bestR_write_lock);

    // Locks for solution matrix append
    omp_lock_t solution_matrix_append_lock;
    omp_init_lock(&solution_matrix_append_lock);

    // Locks for sumsoln
    omp_lock_t solution_sum_lock;
    omp_init_lock(&solution_sum_lock);

    // for each repeat run
    int maxThreads = omp_get_max_threads();

    printf("Running multithreaded over number of threads: %d\n", maxThreads);
    displayProgress1("Running multithreaded over number of threads: " + to_string(maxThreads) + "\n");
    #pragma omp parallel for
    for (int i = 1; i <= repeats; i++)
    {
        // Create run specific structures
        int thread = omp_get_thread_num();
        string tempname2;
        stringstream appendLogBuffer; // stores the trace file log
        stringstream runConsoleOutput; // stores the console message for the run. This is needed for more organized printing output due to multithreading.
        sanneal anneal = anneal_global;
        scost reserve = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        scost change = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

        vector<int> R(puno);
        vector<sspecies> spec = specGlobal; // make local copy of original spec

        appendLogBuffer << "\n Start run loop run " << i << endl;
        try {
            if (runoptions.ThermalAnnealingOn)
            {
                // Annealing Setup
                if (anneal.type == 2)
                {
                    appendLogBuffer << "before initialiseConnollyAnnealing run " << i << endl;

                    initialiseConnollyAnnealing(puno, spno, pu, connections, spec, SMGlobal, cm, anneal, aggexist, R, prop, clumptype, i, thread, appendLogBuffer);
                    
                    appendLogBuffer << "after initialiseConnollyAnnealing run " << i << endl;
                }

                if (anneal.type == 3)
                {
                    appendLogBuffer << "before initialiseAdaptiveAnnealing run " << i << endl;

                    initialiseAdaptiveAnnealing(puno, spno, prop, R, pu, connections, SMGlobal, cm, spec, aggexist, anneal, clumptype, thread, appendLogBuffer);

                    appendLogBuffer << "after initialiseAdaptiveAnnealing run " << i << endl;
                }

                if (verbosity > 1) runConsoleOutput << "\n Run " << i << ": Using Calculated Tinit = " << anneal.Tinit << " Tcool = " << anneal.Tcool << "\n";
                anneal.temp = anneal.Tinit;
            }

            appendLogBuffer << "before computeReserveValue run " << i << endl;

            initialiseReserve(prop, pu, R, rngEngine); // Create Initial Reserve
            SpeciesAmounts(spno,puno,spec,pu,SMGlobal,R,clumptype); // Re-added this from v2.4 because spec amounts need to be refreshed when initializing

            if (aggexist)
                ClearClumps(spno, spec, pu, SMGlobal, thread);

            appendLogBuffer << "after computeReserveValue run " << i << endl;

            if (verbosity > 1)
            {
                computeReserveValue(puno, spno, R, pu, connections, SMGlobal, cm, spec, aggexist, reserve, clumptype, thread, appendLogBuffer);
                runConsoleOutput << "Run " << i << " Init: " << displayValueForPUs(puno, spno, R, reserve, spec, misslevel).str();
            }
            if (verbosity > 5)
            {
                displayTimePassed(startTime);
            }

            if (runoptions.ThermalAnnealingOn)
            {
                appendLogBuffer << "before thermalAnnealing run " << i << endl;

                thermalAnnealing(spno, puno, connections, R, cm, spec, pu, SMGlobal, reserve,
                                repeats, i, savename, misslevel,
                                aggexist, costthresh, tpf1, tpf2, clumptype, anneal, thread, appendLogBuffer);

                if (verbosity > 1)
                {
                    computeReserveValue(puno,spno,R,pu,connections,SMGlobal,cm,spec,aggexist,reserve,clumptype, thread, appendLogBuffer);
                    runConsoleOutput << "Run " << i << " ThermalAnnealing: " << displayValueForPUs(puno,spno,R,reserve,spec,misslevel).str();
                }

                appendLogBuffer << "after thermalAnnealing run " << i << endl;
            }

            if (runoptions.QuantumAnnealingOn)
            {
                appendLogBuffer << "before quantumAnnealing run " << i << endl;

                quantumAnnealing(spno, puno, connections, R, cm, spec, pu, SMGlobal, change, reserve,
                                repeats, i, savename, misslevel,
                                aggexist, costthresh, tpf1, tpf2, clumptype, anneal, thread);

                if (verbosity >1)
                {
                    computeReserveValue(puno,spno,R,pu,connections,SMGlobal,cm,spec,aggexist,reserve,clumptype, thread, appendLogBuffer);
                    runConsoleOutput << "Run " << i << "  QuantumAnnealing: " << displayValueForPUs(puno,spno,R,reserve,spec,misslevel).str();

                }

                appendLogBuffer << "after quantumAnnealing run " << i << endl;
            }

            if (runoptions.HeuristicOn)
            {
                appendLogBuffer << "before Heuristics run " << i << endl;

                Heuristics(spno, puno, pu, connections, R, cm, spec, SMGlobal, reserve,
                        costthresh, tpf1, tpf2, heurotype, clumptype, thread, appendLogBuffer);

                if (verbosity > 1 && (runopts == 2 || runopts == 5))
                {
                    computeReserveValue(puno,spno,R,pu,connections,SMGlobal,cm,spec,aggexist,reserve,clumptype, thread, appendLogBuffer);
                    runConsoleOutput << "Run " << i << "  Heuristic: " << displayValueForPUs(puno, spno, R, reserve, spec, misslevel).str();
                }

                appendLogBuffer << "after Heuristics run " << i << endl;
            }

            if (runoptions.ItImpOn)
            {
                appendLogBuffer << "before iterativeImprovement run " << i << endl;

                iterativeImprovement(puno, spno, pu, connections, spec, SMGlobal, R, cm,
                                    reserve, change, costthresh, tpf1, tpf2, clumptype, i, savename, thread, appendLogBuffer);

                if (itimptype == 3)
                    iterativeImprovement(puno, spno, pu, connections, spec, SMGlobal, R, cm,
                                        reserve, change, costthresh, tpf1, tpf2, clumptype, i, savename, thread, appendLogBuffer);

                appendLogBuffer << "after iterativeImprovement run " << i << endl;

                if (aggexist)
                    ClearClumps(spno, spec, pu, SMGlobal, thread);

                if (verbosity > 1)
                {
                    computeReserveValue(puno, spno, R, pu, connections, SMGlobal, cm, spec, aggexist, reserve, clumptype, thread, appendLogBuffer);
                    runConsoleOutput << "Run " << i << " Iterative Improvement: " << displayValueForPUs(puno, spno, R, reserve, spec, misslevel).str();
                }

            } // Activate Iterative Improvement

            appendLogBuffer << "before file output run " << i << endl;
            string fileNumber = to_string(i);
            if (fnames.saverun)
            {
                tempname2 = savename + "_r" + fileNumber + getFileNameSuffix(fnames.saverun);
                writeSolution(puno, R, pu, tempname2, fnames.saverun, fnames);
            }

            if (fnames.savespecies && fnames.saverun)
            {
                tempname2 = savename + "_mv" + fileNumber + getFileNameSuffix(fnames.savespecies);
                writeSpecies(spno, spec, tempname2, fnames.savespecies, misslevel);
            }

            if (fnames.savesum)
            {   // summaries get stored and aggregated to prevent race conditions.
                summaries[i-1] = computeSummary(puno,spno,R,spec,reserve,i,misslevel,fnames.savesum);
            }

            // Print results from run.
            displayProgress1(runConsoleOutput.str());

            // compute and store objective function score for this reserve system
            computeReserveValue(puno, spno, R, pu, connections, SMGlobal, cm, spec, aggexist, change, clumptype, thread, appendLogBuffer);

            // remember the bestScore and bestRun
            if (change.total < bestScore)
            {
                omp_set_lock(&bestR_write_lock);
                // After locking, do another check in case bestScore has changed
                if (change.total < bestScore) {
                    // this is best run so far
                    bestScore = change.total;
                    bestRun = i;
                    // store bestR
                    bestR = R;
                    bestRunString = runConsoleOutput.str();
                    bestSpec = spec; // make a copy of best spec.
                }
                omp_unset_lock(&bestR_write_lock);
            }

            if (fnames.savesolutionsmatrix)
            {
                appendLogBuffer << "before appendSolutionsMatrix savename " << i << endl;
                tempname2 = savename + "_solutionsmatrix"  + getFileNameSuffix(fnames.savesolutionsmatrix);

                omp_set_lock(&solution_matrix_append_lock);
                appendSolutionsMatrix(i, puno, R, tempname2, fnames.savesolutionsmatrix, fnames.solutionsmatrixheaders);
                omp_unset_lock(&solution_matrix_append_lock);

                appendLogBuffer << "after appendSolutionsMatrix savename " << i << endl;
            }

            // Save solution sum
            if (fnames.savesumsoln) {
                omp_set_lock(&solution_sum_lock);
                for (int k = 0; k < puno; k++) {
                    if (R[k] == 1 || R[k] == 2) {
                        sumsoln[k]++;
                    }
                }
                omp_unset_lock(&solution_sum_lock);
            }

            if (aggexist)
                ClearClumps(spno, spec, pu, SMGlobal, thread);

        }
        catch (exception& e) {
            // On exceptions, append exception to log file in addition to existing buffer. 
            displayProgress1(runConsoleOutput.str());
            appendLogBuffer << "Exception occurred on run " << i << ": " << e.what() << endl;
            displayProgress1(appendLogBuffer.str());
            appendTraceFile(appendLogBuffer.str());

            throw(e);
        }

        appendLogBuffer << "after file output run " << i << endl;
        appendLogBuffer << "end run " << i << endl;

        appendTraceFile(appendLogBuffer.str());

        if (marxanIsSecondary == 1)
            writeSecondarySyncFileRun(i);

        if (verbosity > 1)
            displayTimePassed(startTime);
    }

    // Write all summaries for each run.
    if (fnames.savesum)
    {
        string tempname2 = savename + "_sum" + getFileNameSuffix(fnames.savesum);
        writeSummary(tempname2, summaries, fnames.saverun);
    }

    stringstream bestOut;
    bestOut << "\nBest run: " << bestRun << " Best score: " << bestScore << "\n" << bestRunString;
    cout << bestOut.str();
    displayProgress1(bestOut.str());
} // executeRunLoop

int executeMarxan(string sInputFileName)
{
    int iSparseMatrixFileLength = 0, iSparseMatrixFileLength_sporder = 0;
    long int repeats;
    int puno,spno,gspno;
    double cm,prop;
    int runopts,heurotype,clumptype,itimptype;
    string savename,tempname2;
    double misslevel;
    int iseed = time(NULL),seedinit;
    int aggexist=0,sepexist=0;
    vector<int> R_CalcPenalties;
    vector<int> sumsoln;
    double costthresh,tpf1,tpf2;
    long int itemp;
    int isp;
    int maxThreads = omp_get_max_threads();

    // Handle Error driven termination
    if (setjmp(jmpbuf))
        return 1;

    displayStartupMessage();
    startTime = chrono::high_resolution_clock::now(); // set program start time.

    readInputOptions(cm, prop, anneal_global,
                     iseed, repeats, savename, fnames, sInputFileName,
                     runopts, misslevel, heurotype, clumptype, itimptype, verbosity,
                     costthresh, tpf1, tpf2);

    setDefaultRunOptions(runopts, runoptions);

    sTraceFileName = savename + "_TraceFile.txt";
    createTraceFile();
    appendTraceFile("%s begin execution\n\n",sVersionString.c_str());
    appendTraceFile("LoadOptions\n");

    #ifdef DEBUGCHECKCHANGE
    createDebugFile("debug_MarOpt_CheckChange.csv","ipu,puid,R,total,cost,connection,penalty,threshpen,probability\n",fnames);
    #endif

    #ifdef DEBUGCHANGEPEN
    createDebugFile("debug_MarOpt_ChangePen.csv","ipu,puid,isp,spid,penalty,target,targetocc,occurrence,sepnum,amount,newamount,fractionAmount\n",fnames);
    #endif

    #ifdef DEBUGCALCPENALTIES
    createDebugFile("debug_MarZone_CalcPenalties.csv","\n",fnames);
    #endif

    if (fnames.savelog)
    {
        tempname2 = savename + "_log.dat";
        createLogFile(fnames.savelog,tempname2);
    }

    delta = 1e-14;  // This would more elegantly be done as a constant

    // initt global rng engine
    rngEngine = mt19937(iseed);
    RandSeed1 = iseed;
    seedinit = iseed;

    appendTraceFile("RandSeed iseed %i RandSeed1 %li\n",iseed,RandSeed1);

    // read the data files
    displayProgress1("\nEntering in the data files \n");
    displayProgress3("    Reading in the Planning Unit names \n");
    appendTraceFile("before readPlanningUnits\n");

    itemp = readPlanningUnits(puno,pu,fnames);

    appendTraceFile("after readPlanningUnits\n");
    if (iProbFieldPresent == 1)
        appendTraceFile("prob field present\n");
    else
        appendTraceFile("prob field not present\n");

    #ifdef DEBUG_PROB1D
    if (iProbFieldPresent == 1)
        writeProbData(puno,pu,fnames);
    #endif

    displayProgress1("   There are %i Planning units.\n  %i Planning Unit names read in \n",puno,itemp);
    displayProgress3("    Reading in the species file \n");
    appendTraceFile("before readSpecies\n");

    itemp = readSpecies(spno,specGlobal,fnames);

    appendTraceFile("after readSpecies\n");
    displayProgress1("  %i species read in \n",itemp);
    appendTraceFile("before build search arrays\n");

    // create the fast lookup tables for planning units and species names
    populateLookupTables(puno,spno,pu,specGlobal,PULookup,SPLookup);

    appendTraceFile("after build search arrays\n");

    if (fnames.savesumsoln)
        sumsoln.resize(puno);

    connections.resize(puno);

    displayProgress3("    Reading in the Connection file :\n");
    itemp = 0;
    if (!fnames.connectionname.empty())
    {
        appendTraceFile("before readConnections\n");

        if (!fnames.connectionfilesname.empty())
            writeWeightedConnectivityFile(fnames);

        itemp = readConnections(puno,connections,pu,PULookup,fnames);

        appendTraceFile("after readConnections\n");
        if (asymmetricconnectivity)
        {
            appendTraceFile("Asymmetric connectivity is on.\n");
            writeAsymmetricConnectionFile(puno,connections,pu,fnames);
        }
        if (fOptimiseConnectivityIn)
            appendTraceFile("Optimising 'Connectivity In'.\n");
    }
    displayProgress1("  %i connections entered \n",itemp);
    if (asymmetricconnectivity)
        displayProgress1("  Asymmetric connectivity is on.\n");
    if (fOptimiseConnectivityIn)
        displayProgress1("  Optimising 'Connectivity In'.\n");

    displayProgress3("    Reading in the Planning Unit versus Species File \n");
    appendTraceFile("before readSparseMatrix\n");

    readSparseMatrix(iSparseMatrixFileLength,SMGlobal,puno,spno,pu,PULookup,SPLookup,fnames);

    appendTraceFile("after readSparseMatrix\n");
    if (fProb2D == 1)
        appendTraceFile("Prob2D is on\n");
    else
        appendTraceFile("Prob2D is off\n");

    #ifdef DEBUG_PROB2D
    writeSparseMatrix(iSparseMatrixFileLength,puno,pu,spec,SM,fnames);
    #endif

    if (fnames.saverichness)
    {
        tempname2 = savename + "_richness" + getFileNameSuffix(fnames.saverichness);
        writeRichness(puno,pu,tempname2,fnames.saverichness);
    }

    if (!fnames.matrixspordername.empty())
    {
        appendTraceFile("before readSparseMatrixSpOrder\n");

        readSparseMatrixSpOrder(iSparseMatrixFileLength_sporder,SMsporder,puno,spno,PULookup,SPLookup,specGlobal,fnames);

        appendTraceFile("after readSparseMatrixSpOrder\n");

        #ifdef MEMDEBUG
        displayProgress1("after LoadSparseMatrix_sporder\n");
        #endif
    }

    appendTraceFile("before process block definitions\n");

    if  (!fnames.blockdefname.empty())
    {
        displayProgress1("    Reading in the Block Definition File \n");
        vector<sgenspec> gspec; // declare within this scope of usage
        readSpeciesBlockDefinition(gspno,gspec,fnames);
        setBlockDefinitions(gspno,spno,puno,gspec,specGlobal,pu,SMGlobal);
    }

    setDefaultTargets(specGlobal);
    appendTraceFile("after process block definitions\n");

    appendTraceFile("before computeTotalAreas\n");
    computeTotalAreas(puno,spno,pu,specGlobal,SMGlobal);
    appendTraceFile("after computeTotalAreas\n");

    if (fnames.savetotalareas)
    {
        tempname2 = savename + "totalareas" + getFileNameSuffix(fnames.savetotalareas);
        writeTotalAreas(puno,spno,pu,specGlobal,SMGlobal,tempname2,fnames.savepenalty);
    }

    if (fSpecPROPLoaded > 0)
    {
        appendTraceFile("before computeSpecProp\n");

        // species have prop value specified
        computeSpecProp(spno,specGlobal,puno,pu,SMGlobal);

        appendTraceFile("after computeSpecProp\n");
    }

    displayProgress2("Checking to see if there are aggregating or separating species.\n");
    for (isp=0;isp<spno;isp++)
    {
        if (specGlobal[isp].target2>0)
            aggexist = 1;
        if (specGlobal[isp].sepdistance > 0)
            sepexist = 1;
    }

    // if clumping requirements exist, initialize clumping arrays
    if (aggexist) {
        for (spu& term: SMGlobal) {
            term.clump.assign(maxThreads, 0);
        }
    }

    if (fnames.savesen)
    {
        appendTraceFile("before writeScenario\n");

        tempname2 = savename + "_sen.dat";
        writeScenario(puno,spno,prop,cm,anneal_global,seedinit,repeats,clumptype,
                      runopts,heurotype,costthresh,tpf1,tpf2,tempname2);

        appendTraceFile("after writeScenario\n");
    }

    if (verbosity > 1)
        displayTimePassed(startTime);

    // *******  Pre-processing    ************
    displayProgress1("\nPre-processing Section. \n");
    displayProgress2("    Calculating all the penalties \n");

    R_CalcPenalties.resize(puno);

    // load penalties from file if they are present
    if (!fnames.penaltyname.empty())
    {
        fUserPenalties = 1;

        appendTraceFile("before readPenalties\n");
        readPenalties(specGlobal,spno,fnames,SPLookup);
        appendTraceFile("after readPenalties\n");
    }

    if (runoptions.CalcPenaltiesOn == 0)
    {
        // if penalties have not been loaded, then stop error message
        if (fUserPenalties == 0)
        {
            appendTraceFile("Data error: CalcPenalties off but no PENALTYNAME specified, exiting.\n");
            displayProgress1("Data error: CalcPenalties off but no PENALTYNAME specified, exiting.\n");
            exit(1);
        }

        // transfer loaded penalties to correct data structrure
        applyUserPenalties(specGlobal);
    }
    else
    {
        // we are computing penalties
        if (fnames.matrixspordername.empty())
        {
            appendTraceFile("before CalcPenalties\n");

            // we don't have sporder matrix available, so use slow CalcPenalties method
            itemp = computePenalties(puno,spno,pu,specGlobal,connections,SMGlobal,R_CalcPenalties,aggexist,cm,clumptype, 0);

            appendTraceFile("after CalcPenalties\n");
        }
        else
        {
            // we have sporder matrix available, so use optimised CalcPenalties method
            if (iOptimisationCalcPenalties == 1)
            {
                appendTraceFile("before CalcPenaltiesOptimise\n");

                itemp = computePenaltiesOptimise(puno,spno,pu,specGlobal,connections,SMGlobal,SMsporder,R_CalcPenalties,aggexist,cm,clumptype, 0);

                appendTraceFile("after CalcPenaltiesOptimise\n");
            }
            else
            {
                appendTraceFile("before CalcPenalties\n");

                // we have optimise calc penalties switched off, so use slow CalcPenalties method
                itemp = computePenalties(puno,spno,pu,specGlobal,connections,SMGlobal,R_CalcPenalties,aggexist,cm,clumptype,0);

                appendTraceFile("after CalcPenalties\n");
            }
        }
    }

    if (itemp>0)
        displayProgress("%d species cannot meet target%c.\n",itemp,itemp==1? ' ':'s');

    if (runoptions.ThermalAnnealingOn)
    {
        displayProgress2("    Calculating temperatures.\n");
        if (!anneal_global.Titns)
            displayErrorMessage("Initial Temperature is set to zero. Fatal Error \n");

        anneal_global.Tlen = anneal_global.iterations/anneal_global.Titns;
        displayProgress2("  Temperature length %ld \n",anneal_global.Tlen);
        displayProgress2("  iterations %ld, repeats %ld \n",anneal_global.iterations,repeats);
    } // Annealing Preprocessing. Should be moved to SetAnnealingOptions

    if (fnames.savepenalty)
    {
        tempname2 = savename + "_penalty" + getFileNameSuffix(fnames.savepenalty);
        writePenalty(spno,specGlobal,tempname2,fnames.savepenalty);

        tempname2 = savename + "_penalty_planning_units" + getFileNameSuffix(fnames.savepenalty);
        writePenaltyPlanningUnits(puno,pu,R_CalcPenalties,tempname2,fnames.savepenalty);
    }

    //free(R_CalcPenalties);

    if (fnames.savespec)
    {
        tempname2 = savename + "_spec.csv";
        writeSpec(spno,specGlobal,tempname2);
    }

    if (fnames.savepu)
    {
        tempname2 = savename + "_pu.csv";
        writePu(puno,pu,tempname2);
    }

    // If we are in a runmode with only CalcPenalties, we stop/exit here gracefully because we are finished.
    if (runoptions.HeuristicOn == 0 && runoptions.ThermalAnnealingOn == 0 && runoptions.QuantumAnnealingOn == 0 && runoptions.ItImpOn == 0)
    {
        appendTraceFile("end final file output\n");
        appendTraceFile("\nMarxan end execution\n");
        displayShutdownMessage(startTime);

        if (marxanIsSecondary == 1)
            secondaryExit();

        exit(1);
    }

    if (fnames.savesolutionsmatrix)
    {
        tempname2 = savename + "_solutionsmatrix" + getFileNameSuffix(fnames.savesolutionsmatrix);

        #ifdef DEBUG_CLUSTERANALYSIS
        appendTraceFile("before createSolutionsMatrix savename %s\n",savename);
        #endif

        createSolutionsMatrix(puno,pu,tempname2,fnames.savesolutionsmatrix,fnames.solutionsmatrixheaders);

        #ifdef DEBUG_CLUSTERANALYSIS
        appendTraceFile("after createSolutionsMatrix savename %s\n",savename);
        #endif
    }

    if (fProb1D == 1)
    {
        tempname2 = savename + "_ComputeP_AllPUsSelected_1D.csv";
        ComputeP_AllPUsSelected_1D(tempname2,puno,spno,pu,SMGlobal,specGlobal);
    }

    if (fProb2D == 1)
    {
        tempname2 = savename + "_ComputeP_AllPUsSelected_2D.csv";
        ComputeP_AllPUsSelected_2D(tempname2,puno,spno,pu,SMGlobal,specGlobal);
    }
  
    executeRunLoop(iSparseMatrixFileLength, repeats,puno,spno,cm,aggexist,prop,clumptype,misslevel,
                   savename,costthresh,tpf1,tpf2,heurotype,runopts,
                   itimptype,sumsoln);

    appendTraceFile("before final file output\n");

    if (fnames.savebest)
    {
        tempname2 = savename + "best" + getFileNameSuffix(fnames.savebest);
        writeSolution(puno,bestR,pu,tempname2,fnames.savebest,fnames);

        appendTraceFile("Best solution is run %i\n",bestRun);
        displayProgress1("\nBest solution is run %i\n",bestRun);
    }

    if (fnames.savespecies && fnames.savebest)
    {
        tempname2 = savename + "mvbest" + getFileNameSuffix(fnames.savespecies);

         // TODO ADBAI - need to write best spec, NOT spec_global
        writeSpecies(spno,bestSpec,tempname2,fnames.savespecies,misslevel);
    }

    if (fnames.savesumsoln)
    {   
        tempname2 = savename + "_ssoln" + getFileNameSuffix(fnames.savesumsoln);
        writeSumSoln(puno,sumsoln,pu,tempname2,fnames.savesumsoln);
    }

    if (fnames.savesolutionsmatrix)
    {
        if (fnames.rexecutescript)
        {
            if (marxanIsSecondary == 1)
                secondaryExit();
        }
    }

    displayShutdownMessage(startTime);

    appendTraceFile("end final file output\n");
    appendTraceFile("\nMarxan end execution\n");

    return 0;
} // executeMarxan

void SpeciesAmounts(int spno,int puno, vector<sspecies>& spec, vector<spustuff>& pu, vector<spu>& SM,
                    vector<int>& R,int clumptype)
{
    int i, ism, isp, ipu;

    for (isp = 0; isp < spno; isp++)
    {
        spec[isp].amount = 0;
        spec[isp].occurrence = 0;
        if (spec[isp].target2)
            SpeciesAmounts4(isp, spec, clumptype);

        spec[isp].expected1D = 0;
        spec[isp].expected2D = 0;
        spec[isp].variance1D = 0;
        spec[isp].variance2D = 0;
    }

    for (ipu = 0; ipu < puno; ipu++)
        if (pu[ipu].richness)
            if (R[ipu] == 1 || R[ipu] == 2)
                for (i = 0; i < pu[ipu].richness; i++)
                {
                    ism = pu[ipu].offset + i;
                    isp = SM[ism].spindex;
                    if (spec[isp].target2 == 0)
                    {
                        spec[isp].amount += SM[ism].amount;
                        spec[isp].occurrence++;
                    }
                }
} /*** Species Amounts ***/

// apply settings from the block defintion file for species
void setBlockDefinitions(int gspno,int spno,int puno, vector<sgenspec> &gspec, vector<sspecies> &spec, vector<spustuff> &PU, vector<spu>& SM)
{
    int igsp,isp,ipu;
    double totalamount;

    for (igsp=0;igsp<gspno;igsp++)
    {
        if (gspec[igsp].prop > 0) // deal with percentage in a different way
        {
            for (isp=0;isp<spno;isp++)
            {
                if (spec[isp].type == gspec[igsp].type && spec[isp].target < 0)
                {
                    spec[isp].target = computeTotalSpecAmtAllPu(PU, SM, isp) * gspec[igsp].prop;
                } // Setting target with percentage
            }
        }
        if (gspec[igsp].target > 0)
        {
            for (isp=0;isp<spno;isp++)
            {
                if (spec[isp].type == gspec[igsp].type && spec[isp].target < 0)
                    spec[isp].target = gspec[igsp].target;
            }
        }
        if (gspec[igsp].target2 > 0)
        {
            for (isp=0;isp<spno;isp++)
            {
                if (spec[isp].type == gspec[igsp].type && spec[isp].target2 < 0)
                {
                    spec[isp].target2 = gspec[igsp].target2;
                }
            }
        }
        if (gspec[igsp].targetocc > 0)
        {
            for (isp=0;isp<spno;isp++)
            {
                if (spec[isp].type == gspec[igsp].type && spec[isp].targetocc < 0)
                    spec[isp].targetocc = gspec[igsp].targetocc;
            }
        }
        if (gspec[igsp].sepnum > 0)
        {
            for (isp=0;isp<spno;isp++)
            {
                if (spec[isp].type == gspec[igsp].type && spec[isp].sepnum < 0)
                    spec[isp].sepnum = gspec[igsp].sepnum;
            }
        }
        if (gspec[igsp].spf > 0)
        {
            for (isp=0;isp<spno;isp++)
            {
                if (spec[isp].type == gspec[igsp].type && spec[isp].spf < 0)
                    spec[isp].spf = gspec[igsp].spf;
            }
        }
        if (gspec[igsp].sepdistance > 0)
        {
            for (isp=0;isp<spno;isp++)
            {
                if (spec[isp].type == gspec[igsp].type && spec[isp].sepdistance < 0)
                    spec[isp].sepdistance = gspec[igsp].sepdistance;
            }
        }
        // Percentage is not dealt with here yet. To do this I need to identify
        // target species then determine their total abundance then set target
        // according to percentage
    }
} // setBlockDefinitions

// TODO ADBAI - needs refactoring
// set default run options based on the selection algorithm chosen
void setDefaultRunOptions(int runopts, srunoptions &runoptions)
{
    if (runopts < 0)
        return; // runopts < 0 indicates that these are set in some other way

    switch (runopts)
    {
        case 0:
            runoptions.CalcPenaltiesOn = 1;
            runoptions.ThermalAnnealingOn = 1;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 1;
            runoptions.ItImpOn = 0;
            break;
        case 1:
            runoptions.CalcPenaltiesOn = 1;
            runoptions.ThermalAnnealingOn = 1;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 1;
            break;
        case 2:
            runoptions.CalcPenaltiesOn = 1;
            runoptions.ThermalAnnealingOn = 1;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 1;
            runoptions.ItImpOn = 1;
            break;
        case 3:
            runoptions.CalcPenaltiesOn = 0;
            runoptions.ThermalAnnealingOn = 0;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 1;
            runoptions.ItImpOn = 0;
            break;
        case 4:
            runoptions.CalcPenaltiesOn = 0;
            runoptions.ThermalAnnealingOn = 0;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 1;
            break;
        case 5:
            runoptions.CalcPenaltiesOn = 0;
            runoptions.ThermalAnnealingOn = 0;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 1;
            runoptions.ItImpOn = 1;
            break;
        case 6:
            runoptions.CalcPenaltiesOn = 1;
            runoptions.ThermalAnnealingOn = 1;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 0;
            break;
        case 7:
            runoptions.CalcPenaltiesOn = 1;
            runoptions.ThermalAnnealingOn = 0;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 0;
            break;
        case 8:
            runoptions.CalcPenaltiesOn = 0;
            runoptions.ThermalAnnealingOn = 1;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 0;
            break;
        case 9:
            runoptions.CalcPenaltiesOn = 0;
            runoptions.ThermalAnnealingOn = 1;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 1;
            break;
        case 10:
            runoptions.CalcPenaltiesOn = 0;
            runoptions.ThermalAnnealingOn = 0;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 1;
            break;
        case 11:
            runoptions.CalcPenaltiesOn = 1;
            runoptions.ThermalAnnealingOn = 0;
            runoptions.QuantumAnnealingOn = 1;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 0;
            break;
        case 12:
            runoptions.CalcPenaltiesOn = 1;
            runoptions.ThermalAnnealingOn = 0;
            runoptions.QuantumAnnealingOn = 1;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 1;
            break;
        case 13:
            runoptions.CalcPenaltiesOn = 0;
            runoptions.ThermalAnnealingOn = 0;
            runoptions.QuantumAnnealingOn = 1;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 0;
            break;
        case 14:
            runoptions.CalcPenaltiesOn = 0;
            runoptions.ThermalAnnealingOn = 0;
            runoptions.QuantumAnnealingOn = 1;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 1;
            break;
        default:
            runoptions.CalcPenaltiesOn = 0;
            runoptions.ThermalAnnealingOn = 0;
            runoptions.QuantumAnnealingOn = 0;
            runoptions.HeuristicOn = 0;
            runoptions.ItImpOn = 0;
            break;
    }
} // setDefaultRunOptions

// compute initial penalties for species with a greedy algorithm.
// If species has spatial requirements then CalcPenaltyType4 is used instead
int computePenalties(int puno,int spno, vector<spustuff> &pu, vector<sspecies> &spec,
                     vector<sconnections> &connections, vector<spu>& SM, vector<int> &PUtemp, int aggexist, double cm, int clumptype, int thread)
{
    int i,j,ibest,imaxtarget,itargetocc;
    double ftarget,fbest,fbestrat,fcost,ftemp, rAmount, rAmountBest;
    int badspecies = 0,goodspecies = 0;
    initialiseReserve(0,pu,PUtemp, rngEngine); // Initialize reserve to 0 and fixed. 

    for (i=0;i<spno;i++)
    {
        if (spec[i].target2 || spec[i].sepnum)
        {
            j = CalcPenaltyType4(i,puno,SM,connections,spec,pu,cm,clumptype, thread);
            badspecies += (j>0);
            goodspecies += (j<0);

            appendTraceFile("CalcPenalties spname %i penalty %g\n",spec[i].name,spec[i].penalty);

            continue;
        } // Species has aggregation requirements

        computeFixedPenaltyForSpec(PUtemp, pu, SM, connections, i, ftarget, itargetocc, spec[i].penalty, cm, asymmetricconnectivity);

        // Already adequately represented on type 2 planning unit
        if (ftarget >= spec[i].target && itargetocc >= spec[i].targetocc)
        {
            goodspecies++;
            displayProgress2("Species %i (%s) has already met target %.2f\n",
                             spec[i].name,spec[i].sname.c_str(),spec[i].target);

            appendTraceFile("CalcPenalties spname %i penalty %g\n",spec[i].name,spec[i].penalty);

            continue;
        } // Target met in unremovable reserve

        // Reset non fixed pu
        for (int j = 0; j < puno; j++)
            if (PUtemp[j] < 2) PUtemp[j] = 0; 

        do
        {
            fbest =0; imaxtarget = 0; fbestrat = 0, rAmountBest = 0;
            for (j=0;j<puno;j++)
            { // trying to find best pu
                if (PUtemp[j] == 0)
                {
                    rAmount = returnAmountSpecAtPu(pu[j],SM,i).second;
                    if (rAmount>0) {
                        fcost = computePlanningUnitValue(pu[j],connections[j],cm, asymmetricconnectivity);
                        if (fcost == 0)
                            fcost = delta;
                        if (rAmount >= spec[i].target - ftarget && (imaxtarget == 0 || (imaxtarget == 1 && fcost < fbest)))
                        { // can I meet the target cheaply?
                            imaxtarget = 1;
                            ibest = j;
                            fbest = fcost;
                            rAmountBest = rAmount;
                        } else {
                            if (fbestrat < rAmount/fcost)
                            { // finding the cheapest planning unit
                                fbest = fcost;
                                fbestrat = rAmount/fbest;
                                ibest = j;
                                rAmountBest = rAmount;
                            }
                        }

                        #ifdef DEBUGCALCPENALTIES
                        appendTraceFile("CalcPenalties species %i puid %i cost %g\n",spec[i].name,pu[j].id,fcost);
                        #endif
                    }
                }  // Making sure only checking planning units not already used
            }

            if (fbest > 0)
            {
                PUtemp[ibest] = 1;
                ftarget += rAmountBest;
                itargetocc++;
                spec[i].penalty += fbest;

                #ifdef DEBUGCALCPENALTIES
                appendTraceFile("CalcPenalties species %i puid %i ftarget %g fbest %g\n",spec[i].name,pu[ibest].id,ftarget,fbest);
                #endif
            } // Add pu to target
        } while ((ftarget <spec[i].target|| itargetocc < spec[i].targetocc) && fbest > 0); // or no more pu left

        if (fbest == 0) // Could not meet target using all available PUs
        { // If not met target with all available PUs
            displayProgress2("Species %d (%s) cannot reach target %.2f there is only %.2f available.\n",
                           spec[i].name,spec[i].sname.c_str(),spec[i].target,ftarget);
            if (ftarget==0)
                ftarget=delta;  // Protect against divide by zero
            ftemp = 0;
            if (ftarget<spec[i].target)
                ftemp = spec[i].target/ftarget;
            if (itargetocc < spec[i].targetocc && itargetocc)  // If ! itargetocc then also !ftarget
                ftemp += (float) spec[i].targetocc/(float) itargetocc;
            spec[i].penalty = spec[i].penalty * ftemp; // Scale it up
            // This value will be ~ 1/delta when there are no occ's of target species in system
            badspecies++;
        }

        #ifdef DEBUGTRACEFILE
        appendTraceFile("CalcPenalties spname %i penalty %g target %g\n",spec[i].name,spec[i].penalty,spec[i].target);
        #endif
    }  // Penalty for each individual Species
    // Clear clumps in case I needed them for target4 species

    if (aggexist)
        ClearClumps(spno,spec,pu,SM, thread);

    if (goodspecies)
        displayProgress1("%i species are already adequately represented.\n",goodspecies);

    return(badspecies);
}

// compute initial penalties for species with a greedy algorithm.
int computePenaltiesOptimise(int puno,int spno, vector<spustuff> &pu, vector<sspecies> &spec,
                             vector<sconnections> &connections, vector<spu>& SM, vector<spusporder> &SMsp,
                             vector<int> &PUtemp, int aggexist, double cm, int clumptype, int thread)
{
    int i,j,ibest,imaxtarget,itargetocc,ism,ipu;
    double ftarget,fbest,fbestrat,fcost,ftemp, rAmount, r_ibest_amount;
    int badspecies = 0,goodspecies = 0;

    appendTraceFile("CalcPenaltiesOptimise start\n");

    initialiseReserve(puno,pu,PUtemp, rngEngine); // Adds existing reserve to PUtemp

    for (i=0;i<spno;i++)
    {
        appendTraceFile("CalcPenaltiesOptimise spname %i\n",spec[i].name);

        if (spec[i].target2 || spec[i].sepnum)
        {
            j = CalcPenaltyType4(i,puno,SM,connections,spec,pu,cm,clumptype, thread);
            badspecies += (j>0);
            goodspecies += (j<0);
            continue;
        } // Species has aggregation requirements

        ftarget = 0;
        itargetocc = 0;
        spec[i].penalty = 0;

        if (spec[i].richness > 0)
        {
            for (j=0;j<spec[i].richness;j++)  // traverse pu's containing this sp
            { // reset PUtemp and also target
                ism = spec[i].offset + j;
                ipu = SMsp[ism].puindex;

                if (PUtemp[ipu] < 2)
                    PUtemp[ipu] = 0;
                if (PUtemp[ipu] == 2)
                {
                    ftarget += SMsp[ism].amount;
                    itargetocc++;
                    spec[i].penalty += computePlanningUnitValue(pu[ipu],connections[ipu],cm, asymmetricconnectivity);
                }
            }
        }

        // Already adequately represented on type 2 planning unit
        if (ftarget >= spec[i].target && itargetocc >= spec[i].targetocc)
        { // Target met in unremovable reserve
            goodspecies++;
            displayProgress2("Species %i (%s) has already met target %.2f\n",
                             spec[i].name,spec[i].sname.c_str(),spec[i].target);
            continue;
        }

        do
        {
            fbest =0; imaxtarget = 0; fbestrat = 0;
            if (spec[i].richness > 0)
            {
                for (j=0;j<spec[i].richness;j++)  // traverse pu's containing this sp
                { // trying to find best pu
                    ism = spec[i].offset + j;
                    ipu = SMsp[ism].puindex;

                    rAmount = SMsp[ism].amount;
                    if (PUtemp[ipu] == 0)
                    { // Making sure only checking planning units not already used
                        fcost = computePlanningUnitValue(pu[ipu],connections[ipu],cm, asymmetricconnectivity);
                        if (fcost == 0)
                           fcost = delta;
                        if (rAmount >= spec[i].target - ftarget && (imaxtarget == 0
                            || (imaxtarget == 1 && fcost < fbest)))
                        { // can I meet the target cheaply?
                            imaxtarget = 1;
                            ibest = ipu;
                            r_ibest_amount = rAmount;
                            fbest = fcost;
                        } else {
                            if (fbestrat < rAmount/fcost)
                            { // finding the cheapest planning unit
                                fbest = fcost;
                                fbestrat = rAmount/fbest;
                                ibest = ipu;
                                r_ibest_amount = rAmount;
                            }
                        }
                    }
                }
            }

            if (fbest > 0)
            { // Add pu to target
                PUtemp[ibest] = 1;
                ftarget += r_ibest_amount;
                itargetocc++;
                spec[i].penalty += fbest;
            }
        } while ((fbest > 0) && (ftarget <spec[i].target|| itargetocc < spec[i].targetocc));
        // while there is some pu's with this species to test AND a best available pu was found AND targets are not met yet

        if (fbest == 0) // Could not meet target using all available PUs
        { // If not met target with all available PUs
            displayProgress2("Species %d (%s) cannot reach target %.2f there is only %.2f available.\n",
                             spec[i].name,spec[i].sname.c_str(),spec[i].target,ftarget);
            if (ftarget==0)
                ftarget=delta;  // Protect against divide by zero
            ftemp = 0;
            if (ftarget<spec[i].target)
                ftemp = spec[i].target/ftarget;
            if (itargetocc < spec[i].targetocc && itargetocc)  // If ! itargetocc then also !ftarget
                ftemp += (float) spec[i].targetocc/(float) itargetocc;
            spec[i].penalty = spec[i].penalty * ftemp; // Scale it up
            /* This value will be ~ 1/delta when there are no occ's of target species in system*/
            badspecies++;
        }

        appendTraceFile("CalcPenaltiesOptimise spname %i penalty %g\n",spec[i].name,spec[i].penalty);
    }  // Penalty for each individual Species
    
    // Clear clumps in case I needed them for target4 species
    if (aggexist)
        ClearClumps(spno,spec,pu,SM, thread);

    if (goodspecies)
        displayProgress1("%i species are already adequately represented.\n",goodspecies);

    appendTraceFile("CalcPenaltiesOptimise end\n");

    return(badspecies);
} // computePenaltiesOptimise

// compute change in the species representation for adding or removing a single planning unit or set of planning units
double computeChangePenalty(int ipu,int puno, vector<sspecies>& spec, vector<spustuff>& pu, vector<spu>& SM,
                          vector<int>& R, vector<sconnections>& connections, int imode, int clumptype, double& rShortfall, int thread)
{
    int i, ism, isp;
    double fractionAmount, penalty, newamount, tamount;
    double rOldShortfall, rNewAmountHeld, rNewShortfall;
    vector<sclumps> tempSclumps;
#ifdef DEBUGCHANGEPEN
    char debugline[200];
#endif

#ifdef ANNEALING_TEST
    if (ipu == (puno - 1))
        appendTraceFile("computeChangePenalty start\n");
#endif

    rShortfall = 0;
    penalty = 0;

    if (pu[ipu].richness)
    {
        for (i = 0; i < pu[ipu].richness; i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            if (SM[ism].amount)  /** Only worry about PUs where species occurs and target != 0 **/
            {
                fractionAmount = 0;
                newamount = 0; /* Shortfall */

                rOldShortfall = 0;
                rNewShortfall = 0;
                
                if (spec[isp].target > spec[isp].amount && spec[isp].target != 0)
                {
                    fractionAmount = (spec[isp].target - spec[isp].amount) / spec[isp].target;
                    rOldShortfall = spec[isp].target - spec[isp].amount;
                }

                rNewAmountHeld = spec[isp].amount + (SM[ism].amount * imode);
                if (spec[isp].target > rNewAmountHeld)
                    rNewShortfall = spec[isp].target - rNewAmountHeld;
                rShortfall += rNewShortfall - rOldShortfall;

                // does this species have occurrence target?
                if (spec[isp].targetocc > 0)
                {
                    if (spec[isp].targetocc > spec[isp].occurrence)
                        fractionAmount += ((double)spec[isp].targetocc - (double)spec[isp].occurrence) /
                                          (double)spec[isp].targetocc;

                    if (spec[isp].target && spec[isp].targetocc)
                        fractionAmount /= 2;
                }

                if (spec[isp].sepnum)
                    fractionAmount += computeSepPenalty(spec[isp].separation, spec[isp].sepnum);

                if (spec[isp].target2)
                {
                    /* clumping species */
                    /* New Pen 4 includes occurrences, amounts and separation target */
                    newamount = NewPenalty4(ipu, isp, puno, spec, pu, SM, R, connections, imode, clumptype, thread);
                }
                else
                {
                    if (spec[isp].target)
                        newamount = computeSpeciesPlanningUnitPenalty(ipu, isp, spec, pu, SM, imode) / spec[isp].target;
                    if (spec[isp].targetocc)
                    {
                        tamount = (double)(spec[isp].targetocc - spec[isp].occurrence - imode) /
                                  (double)spec[isp].targetocc;
                        newamount += tamount < 0 ? 0 : tamount;
                    }
                    if (spec[isp].target && spec[isp].targetocc)
                        newamount /= 2;
                    if (spec[isp].sepnum)
                        newamount += computeSepPenalty(CountSeparation2(isp, ipu, tempSclumps, puno, R, pu, SM, spec, imode, thread),
                                                spec[isp].sepnum); /* I need a new function here */
#ifdef ANNEALING_TEST
                    if (ipu == (puno - 1))
                    {
                        appendTraceFile("penalty %g spf %g newamount %g fractionAmount %g target %g amount %g\n",
                                        spec[isp].penalty, spec[isp].spf, newamount, fractionAmount, spec[isp].target, spec[isp].amount);
                    }
#endif
                } /* no target2 */

                penalty += spec[isp].penalty * spec[isp].spf * (newamount - fractionAmount);

            }

#ifdef DEBUGCHANGEPEN
            sprintf(debugline, "%i,%i,%i,%i,%g,%g,%i,%i,%i,%g,%g,%g\n",
                    ipu, pu[ipu].id, isp, spec[isp].name, penalty, spec[isp].target,
                    spec[isp].targetocc, spec[isp].occurrence, spec[isp].sepnum,
                    spec[isp].amount, newamount, fractionAmount);
            appendDebugFile("debug_MarOpt_ChangePen.csv", debugline, fnames);
#endif
        }
    }

#ifdef ANNEALING_TEST
    if (ipu == (puno - 1))
    {
        sprintf(debugbuffer, "computeChangePenalty end penalty %g\n", penalty);
        appendTraceFile(debugbuffer);
    }
#endif

    if (isinf(penalty) != 0)
    {
        printf("computeChangePenalty infinite fractionAmount >%g<\n", penalty);
    }

    return (penalty);
} // computeChangePenalty

// compute objective function value of a reserve system
void computeReserveValue(int puno,int spno, vector<int> &R, vector<spustuff> &pu,
                         vector<sconnections> &connections, vector<spu>& SM,
                         double cm, vector<sspecies> &spec, int aggexist, scost &reserve,int clumptype, int thread, stringstream& logBuffer)
{
    vector<sclumps> tempSclumps;
    int i,j;
    double fractionAmount;
    vector<double> ExpectedAmount1D, VarianceInExpectedAmount1D,
           ExpectedAmount2D, VarianceInExpectedAmount2D;
    double rConnectivityValue;
    string sProbDebugFileName; // for debugging only

    // init arrays for prob 1D and 2D
    if (fProb1D == 1)
    {
        ExpectedAmount1D.assign(spno, 0);
        VarianceInExpectedAmount1D.assign(spno, 0);
    }
    if (fProb2D == 1)
    {
        ExpectedAmount2D.assign(spno, 0);
        VarianceInExpectedAmount2D.assign(spno, 0);
    }

    reserve.pus = 0;
    reserve.cost = 0;
    reserve.penalty = 0;
    reserve.connection = 0;
    reserve.shortfall = 0;
    reserve.probability1D = 0;
    reserve.probability2D = 0;

    if (aggexist)
        SetSpeciesClumps(puno,R,spec,pu,SM,connections,clumptype, thread);
        
    // traverse species, computing penalty for each species
    for (i=0;i<spno;i++)
    {
        fractionAmount = 0;
        if (spec[i].target > spec[i].amount)
        {
            if (spec[i].target != 0) {
               fractionAmount = (spec[i].target-spec[i].amount )/ spec[i].target;
            }

            reserve.shortfall += spec[i].target-spec[i].amount;
        }
        
        // does this species have an occurrence target?
        if (spec[i].targetocc > 0)
        {
            if (spec[i].targetocc > spec[i].occurrence)
            {
                fractionAmount += (double) (spec[i].targetocc - spec[i].occurrence)/(double) spec[i].targetocc;
                reserve.shortfall += spec[i].targetocc-spec[i].occurrence;
            }
            if (spec[i].target && spec[i].targetocc)
                fractionAmount /= 2;
        }

        reserve.penalty += fractionAmount * spec[i].penalty * spec[i].spf;

        if (spec[i].sepnum)
        {
            spec[i].separation = CountSeparation2(i,0,tempSclumps,puno,R,pu,SM,spec,0, thread);
            reserve.penalty += computeSepPenalty(spec[i].separation,spec[i].sepnum) *
                                spec[i].spf*spec[i].penalty;
        }
    }

    // traverse planning units, computing planning unit metrics for the reserve
    for (j=0;j<puno;j++)
    {
        // if planning unit is protected
        if (R[j]==1 || R[j] == 2)
        {
            reserve.cost += pu[j].cost;
            reserve.pus += 1;
            rConnectivityValue = ConnectionCost2(j,connections,R,1,0,cm, asymmetricconnectivity, fOptimiseConnectivityIn);
            reserve.connection += rConnectivityValue;

            #ifdef DEBUG_RESERVECOST
            logBuffer << "puid " << pu[j].id << " connectivity " << rConnectivityValue << endl;
            #endif

            if (fProb1D == 1)
                ReturnProbabilityAmounts1D(ExpectedAmount1D,VarianceInExpectedAmount1D,j,puno,pu,SM);
            if (fProb2D == 1)
                ReturnProbabilityAmounts2D(ExpectedAmount2D,VarianceInExpectedAmount2D,j,puno,pu,SM);
        }
    }

    if (fProb1D == 1)
    {
        reserve.probability1D = ComputeProbability1D(ExpectedAmount1D,VarianceInExpectedAmount1D,spno,spec);
    }

    if (fProb2D == 1)
    {
        reserve.probability2D = ComputeProbability2D(ExpectedAmount2D,VarianceInExpectedAmount2D,spno,spec);
    }

    reserve.total = reserve.cost + reserve.connection + reserve.penalty + reserve.probability1D + reserve.probability2D;

    #ifdef DEBUG_PROB1D
    logBuffer << "probability1D " << reserve.probability1D << endl;

    sProbDebugFileName = fnames.outputdir + "output_Prob1DDebug_" + to_string(reserve.cost) + ".csv";
    writeProb1DDebugTable(spno,sProbDebugFileName,
                          ExpectedAmount1D,VarianceInExpectedAmount1D,spec);

    sProbDebugFileName = fnames.outputdir + "output_Prob1DDetailDebug_" + to_string(reserve.cost) + ".csv";
    writeProb1DDetailDebugTable(sProbDebugFileName,puno,spno,pu,SM,R);
    #endif

    #ifdef DEBUG_PROB2D
    logBuffer << "probability2D " << reserve.probability2D << endl;

    sProbDebugFileName = fnames.outputdir + "output_Prob2DDebug_" + to_string(reserve.cost) + ".csv";
    writeProb2DDebugTable(spno,sProbDebugFileName,
                          ExpectedAmount2D,VarianceInExpectedAmount2D,spec);

    sProbDebugFileName = fnames.outputdir + "output_Prob2DDetailDebug_" + to_string(reserve.cost) + ".csv";
    writeProb2DDetailDebugTable(sProbDebugFileName,puno,pu,SM,R);
    #endif
    
    // destroy arrays for prob 1D and 2D
    if (fProb1D == 1)
    {
        for (i=0;i<spno;i++)
        {
            spec[i].expected1D = ExpectedAmount1D[i];
            spec[i].variance1D = VarianceInExpectedAmount1D[i];
        }
    }
    if (fProb2D == 1)
    {
        for (i=0;i<spno;i++)
        {
            spec[i].expected2D = ExpectedAmount2D[i];
            spec[i].variance2D = VarianceInExpectedAmount2D[i];
        }
    }
} // computeReserveValue

// sets cost threshold penalty when "cost threshold" is in use
double thresholdPenalty(double tpf1,double tpf2,double timeprop)
{
    if (tpf2 < 0)
        return(tpf1);
    return(tpf1*exp(tpf2*timeprop));
}

void computeChangeScore(int iIteration,int ipu,int spno,int puno,vector<spustuff> &pu, vector<sconnections> &connections,
                        vector<sspecies> &spec, vector<spu>& SM, vector<int> &R, double cm, int imode,
                        scost &change, scost &reserve,double costthresh,double tpf1, double tpf2,
                        double timeprop,int clumptype, int thread)
// imode = 1 add PU, imode = -1 remove PU
{
    double threshpen = 0;
    int threshtype = 1; /*Debugging line. This should be input parameter not hardwired */
    double tchangeconnection, tresconnection;
#ifdef DEBUGCHECKCHANGE
    char debugline[200];
#endif

    change.cost = pu[ipu].cost * imode; /* Cost of this PU on it's own */
    change.connection = ConnectionCost2(ipu, connections, R, imode, 1, cm, asymmetricconnectivity, fOptimiseConnectivityIn);


    change.penalty = computeChangePenalty(ipu, puno, spec, pu, SM, R, connections, imode, clumptype, change.shortfall, thread);

    if (costthresh)
    {
        // Threshold Penalty for costs
        if (reserve.cost + reserve.connection <= costthresh)
        {
            if (change.cost + change.connection + reserve.cost + reserve.connection <= costthresh)
                threshpen = 0;
            else
                threshpen = (change.cost + change.connection +
                             reserve.cost + reserve.connection - costthresh) *
                            thresholdPenalty(tpf1, tpf2, timeprop);
        }
        else
        {
            if (change.cost + change.connection + reserve.cost + reserve.connection <= costthresh)
                threshpen = (reserve.cost + reserve.connection - costthresh) *
                            thresholdPenalty(tpf1, tpf2, timeprop);
            else
                threshpen = (change.cost + change.connection) *
                            thresholdPenalty(tpf1, tpf2, timeprop);
        }
    }

    change.threshpen = threshpen;

    if (fProb1D == 1)
        change.probability1D = ChangeProbability1D(iIteration, ipu, spno, puno, spec, pu, SM, imode);
    else
        change.probability1D = 0;
    if (fProb2D == 1)
        change.probability2D = ChangeProbability2D(iIteration, ipu, spno, puno, spec, pu, SM, imode);
    else
        change.probability2D = 0;

    change.total = change.cost + change.connection + change.penalty + change.threshpen + change.probability1D + change.probability2D;

#ifdef DEBUGCHECKCHANGE
    sprintf(debugline, "%i,%i,%i,%g,%g,%g,%g,%g,%g,%g\n",
            ipu, pu[ipu].id, R[ipu], change.total, change.cost, change.connection, change.penalty, change.threshpen, change.probability1D, change.probability2D);
    appendDebugFile("debug_MarOpt_CheckChange.csv", debugline, fnames);
#endif
} // computeChangeScore

// compute change in the objective function score for adding or removing a set of planning units
void computeQuantumChangeScore(int spno,int puno, vector<spustuff>& pu, vector<sconnections>& connections,
                               vector<sspecies>& spec, vector<spu>& SM, vector<int>& R,double cm,
                               scost& change, scost& reserve, double costthresh, double tpf1, double tpf2,
                               double timeprop,int clumptype,int iFluctuationCount, vector<int>& PUChosen, int thread)
// imode = 1 add PU, imode = -1 remove PU
{
    // We query a whole bunch of changes in one, passed in by Quantum annealing.
    double threshpen = 0;
    int imode, i, j, threshtype = 1; // Debugging line. This should be input parameter not hardwired
    double tchangeconnection, tresconnection;
#ifdef DEBUG_QA
    char debugline[200];
#endif

#ifdef DEBUG_QA
    appendTraceFile("computeQuantumChangeScore start iFluctuationCount %i\n", iFluctuationCount);
#endif

    change.cost = 0;
    change.connection = 0;
    change.penalty = 0;
    change.shortfall = 0;
    change.probability1D = 0;
    change.probability2D = 0;
    j = -1;

    for (i = 0; i < iFluctuationCount; i++)
    {
        do
            j++;

        while (PUChosen[j] < 1);

        imode = R[j] == 1 ? -1 : 1;

#ifdef DEBUG_QA
        appendTraceFile("computeQuantumChangeScore ipu %i chosen %i imode %i\n", j, PUChosen[j], imode);
#endif

        change.cost += pu[j].cost * imode; /* Cost of this PU on it's own */
        change.connection += ConnectionCost2(j, connections, R, imode, 1, cm, asymmetricconnectivity, fOptimiseConnectivityIn);
        if (threshtype == 1)
        {
            tchangeconnection = change.connection;
            tresconnection = reserve.connection;
            change.connection = 0;
            reserve.connection = 0;
        }

        change.penalty += computeChangePenalty(j, puno, spec, pu, SM, R, connections, imode, clumptype, change.shortfall, thread);

        if (costthresh)
        {
            // Threshold Penalty for costs
            if (reserve.cost + reserve.connection <= costthresh)
            {
                if (change.cost + change.connection + reserve.cost + reserve.connection <= costthresh)
                    threshpen = 0;
                else
                    threshpen = (change.cost + change.connection +
                                 reserve.cost + reserve.connection - costthresh) *
                                thresholdPenalty(tpf1, tpf2, timeprop);
            }
            else
            {
                if (change.cost + change.connection + reserve.cost + reserve.connection <= costthresh)
                    threshpen = (reserve.cost + reserve.connection - costthresh) *
                                thresholdPenalty(tpf1, tpf2, timeprop);
                else
                    threshpen = (change.cost + change.connection) *
                                thresholdPenalty(tpf1, tpf2, timeprop);
            }
        }

        change.threshpen = threshpen;

        if (threshtype == 1)
        {
            change.connection = tchangeconnection;
            reserve.connection = tresconnection;
        }

        if (fProb1D == 1)
            change.probability1D += ChangeProbability1D(-1, j, spno, puno, spec, pu, SM, imode);
        else
            change.probability1D = 0;
        if (fProb2D == 1)
            change.probability2D += ChangeProbability2D(-1, j, spno, puno, spec, pu, SM, imode);
        else
            change.probability2D = 0;
    }

    change.total = change.cost + change.connection + change.penalty + change.threshpen + change.probability1D + change.probability2D;

#ifdef DEBUGCHECKCHANGE
    sprintf(debugline, "%i,%i,%i,%g,%g,%g,%g,%g,%g\n", j, pu[j].id, R[j], change.total, change.cost, change.connection, change.penalty, change.threshpen, change.probability);
    appendDebugFile("debug_MarOpt_CheckChange.csv", debugline, fnames);
#endif

#ifdef DEBUG_QA
    appendTraceFile("computeQuantumChangeScore end\n");
#endif
} // computeQuantumChangeScore

// determines if the change value for changing a single planning unit status is good
// does the change stochastically fall below the current acceptance probability?
int isGoodChange(scost& change,double temp, uniform_real_distribution<double>& float_range)
{
    if (change.total <= 0)
        return 1;
    else
        return (exp(-change.total / temp) > float_range(rngEngine));
}

// determines if the change value for changing status for a set of planning units is good
// does it stochastically fall below the current acceptance probability?
int isGoodQuantumChange(struct scost change,double rProbAcceptance, uniform_real_distribution<double>& float_range)
{
    if (change.total <= 0)
        return 1;
    else
        return (rProbAcceptance > float_range(rngEngine)) ? 1 : 0;
}

// change the status of a single planning unit
void doChange(int ipu,int puno,vector<int> &R, scost &reserve, scost &change,
              vector<spustuff> &pu,vector<spu>& SM,vector<sspecies> &spec,vector<sconnections> &connections,
              int imode,int clumptype, int thread, stringstream& logBuffer)
{
    int i, ism, isp;
    double rAmount;

    R[ipu] = imode == 1 ? 1 : 0;
    reserve.pus += imode;
    reserve.cost += change.cost;
    reserve.connection += change.connection;
    reserve.penalty += change.penalty;
    reserve.probability1D += change.probability1D;
    reserve.probability2D += change.probability2D;
    reserve.shortfall += change.shortfall;

    if (pu[ipu].richness)
    { // Invoke Species Change
        for (i = 0; i < pu[ipu].richness; i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            rAmount = SM[ism].amount;

            if (spec[isp].target2 && rAmount > 0)
            { // Type 4 species and this will impact them
                if (imode == 1)
                {
                    AddNewPU(ipu, isp, connections, spec, pu, SM, clumptype, thread);
                }
                else
                {
                    RemPu(ipu, isp, connections, spec, pu, SM, clumptype, thread);
                }
                if (spec[isp].occurrence < 0)
                {
                    printf("Warning Warning ! isp %i occ %i \n", isp, spec[isp].occurrence);
                }
            }
            else
            { // No clumping species
                spec[isp].occurrence += (rAmount > 0) * imode;
                spec[isp].amount += rAmount * imode;

                if (spec[isp].amount < 0.0001)
                    if (spec[isp].amount > -0.0001)
                        spec[isp].amount = 0;

                if (fProb1D == 1)
                {
                    spec[isp].expected1D += imode * rAmount * (1 - pu[ipu].prob);
                    spec[isp].variance1D += imode * rAmount * rAmount * pu[ipu].prob * (1 - pu[ipu].prob);
                }
                if (fProb2D == 1)
                {
                    spec[isp].expected2D += imode * rAmount * SM[ism].prob;
                    spec[isp].variance2D += imode * rAmount * rAmount * SM[ism].prob * (1 - SM[ism].prob);
                }

#ifdef ANNEALING_TEST
                logBuffer << "doChange ipu " << ipu << " isp " << isp << " spec.amount " << spec[isp].amount << " imode " << imode << endl;
#endif
            }

            if (spec[isp].sepnum > 0) /* Count separation but only if it is possible that it has changed */
                if ((imode == 1 && spec[isp].separation < spec[isp].sepnum) || (imode == -1 && spec[isp].separation > 1))
                {
                    vector<sclumps> tempSclumps;
                    spec[isp].separation = CountSeparation2(isp, 0, tempSclumps, puno, R, pu, SM, spec, 0, thread);
                }
        }
    }

    reserve.total = reserve.cost + reserve.connection + reserve.penalty + reserve.probability1D + reserve.probability2D;
} // doChange

// change the status of a set of planning units
void doQuantumChange(int puno, vector<int>& R, scost& reserve, scost& change,
                     vector<spustuff>& pu,vector<spu>& SM, vector<sspecies> spec, vector<sconnections>& connections,
                     int clumptype, int iFluctuationCount, vector<int>& PUChosen, int thread)
{
    // We accept a whole bunch of changes in one, passed in by Quantum annealing.
    vector<sclumps> tempSclumps;
    int i,j,ipu,ism,isp,imode;
    double rAmount;

    #ifdef DEBUG_QA
    appendTraceFile("doQuantumChange start\n");
    #endif

    reserve.cost += change.cost;
    reserve.connection += change.connection;
    reserve.penalty += change.penalty;
    reserve.probability1D += change.probability1D;
    reserve.probability2D += change.probability2D;

    ipu = -1;
    for (j=0;j<iFluctuationCount;j++)
    {
        do
            ipu++;

        while (PUChosen[ipu] < 1);

        imode = R[ipu] == 1 ? -1 : 1;
        R[ipu] = imode == 1 ? 1 : 0;
        reserve.pus += imode;

        if (pu[ipu].richness)
        { // Invoke Species Change
            for (i=0;i<pu[ipu].richness;i++)
            {
                ism = pu[ipu].offset + i;
                isp = SM[ism].spindex;

                rAmount = SM[ism].amount;

                if (spec[isp].target2 && rAmount > 0)
                { // Type 4 species and this will impact them
                    if (imode == 1)
                    {
                        AddNewPU(ipu,isp,connections,spec,pu,SM,clumptype, thread);
                    } else {
                        RemPu(ipu,isp,connections,spec,pu,SM,clumptype, thread);
                    }
                    if (spec[isp].occurrence < 0)
                    {
                        printf("Warning Warning ! isp %i occ %i \n",isp,spec[isp].occurrence);
                    }
                }
                else
                { // No clumping species
                    spec[isp].occurrence += (rAmount > 0)*imode;
                    spec[isp].amount += rAmount*imode;

                    if (spec[isp].amount < 0.0001)
                    {
                        if (spec[isp].amount > -0.0001)
                            spec[isp].amount = 0;
                    }

                    #ifdef ANNEALING_TEST
                    appendTraceFile("doChange ipu %i isp %i spec.amount %g imode %i\n",
                                         ipu,isp,spec[isp].amount,imode);
                    #endif
                }

                if (spec[isp].sepnum>0) // Count separation but only if it is possible that it has changed
                    if ((imode ==1 && spec[isp].separation < spec[isp].sepnum) || (imode == -1 && spec[isp].separation >1))
                        spec[isp].separation = CountSeparation2(isp,0,tempSclumps,puno,R,pu,SM,spec,0, thread);
            }
        }
    }

    reserve.total = reserve.cost + reserve.connection + reserve.penalty + reserve.probability1D + reserve.probability2D;

    #ifdef DEBUG_QA
    appendTraceFile("doQuantumChange end\n");
    #endif
} // doQuantumChange

// handle command line parameters for the marxan executable
void handleOptions(int argc,char *argv[], string sInputFileName)
{
    if (argc>4)
    {  // if more than one commandline argument then exit
        displayUsage(argv[0]);
        exit(1);
    }

    for (int i=1;i<argc;i++)
    { // Deal with all arguments
        if (argv[i][0] == '/' || argv[i][0] == '-')
        {
            switch(argv[i][1])
            {
                case 'C':
                case 'c':
                case 'S':
                case 's':
                    marxanIsSecondary = 1;
                break;
            default:
                fprintf(stderr,"unknown option %s\n",argv[i]);
                break;
            }
        } else {
            sInputFileName = argv[i]; /* If not a -option then must be input.dat name */
        }
    }
}

void initialiseConnollyAnnealing(int puno,int spno,vector<spustuff> &pu, vector<sconnections> &connections, vector<sspecies> &spec,
                                 vector<spu>& SM,double cm, sanneal &anneal,int aggexist,
                                 vector<int> &R,double prop,int clumptype,int irun, int thread, stringstream& logBuffer)
{
    long int i,ipu,imode, iOldR;
    double deltamin = 0,deltamax = 0;
    double localdelta = 1E-10;
    scost change = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    scost reserve = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    #ifdef DEBUGTRACEFILE
    FILE *fp=nullptr;
    if (verbosity > 4)
    {
        string writename = fnames.outputdir + "debug_maropt_initialiseConnollyAnnealing_" + to_string(irun) + ".csv";
        fp = fopen(writename.c_str(),"w");
        if (fp==NULL)
            displayErrorMessage("cannot create debug_maropt_initialiseConnollyAnnealing file %s\n",writename.c_str());
        fprintf(fp,"i,ipu,puid,old R,imode,R,total,max,min\n");
    }
    #endif

    #ifdef DEBUG_PROB1D
    logBuffer << "initialiseConnollyAnnealing A - before initialise reserve\n";
    #endif

    initialiseReserve(prop,pu,R, rngEngine);

    #ifdef DEBUG_PROB1D
    logBuffer << "initialiseConnollyAnnealing B - after initialise reserve\n";
    #endif

    if (aggexist)
        ClearClumps(spno,spec,pu,SM, thread);

    #ifdef DEBUG_PROB1D
    logBuffer << "initialiseConnollyAnnealing C - before compute reserve\n";
    #endif

    computeReserveValue(puno,spno,R,pu,connections,SM,cm,spec,aggexist, reserve,clumptype, thread, logBuffer);

    #ifdef DEBUG_PROB1D
    logBuffer << "initialiseConnollyAnnealing D - after compute reserve\n";
    #endif

    std::uniform_int_distribution<int> int_range(0, puno-1);
    for (i=1;i<= anneal.iterations/100; i++)
    {
        ipu = int_range(rngEngine);
        iOldR = R[ipu];
        imode = R[ipu]==1?-1:1;

        computeChangeScore(-1,ipu,spno,puno,pu,connections,spec,SM,R,cm,imode,change,reserve,0,0,0,0,clumptype, thread);
        doChange(ipu,puno,R,reserve,change,pu,SM,spec,connections,imode,clumptype, thread, logBuffer);
        if (change.total > deltamax)
            deltamax = change.total;
        if (change.total >localdelta && (deltamin < localdelta || change.total < deltamin))
            deltamin = change.total;

        if (verbosity > 4)
            fprintf(fp,"%li,%li,%i,%li,%li,%i,%g,%g,%g\n",i,ipu,pu[ipu].id,iOldR,imode,R[ipu],change.total,deltamax,deltamin);
                                                   // i,ipu,puid,R,imode,iZone,total,max,min
    }  // Run through this bit for iterations/100 times

    anneal.Tinit = deltamax;
    deltamin *= 0.1;
    anneal.Tcool = exp(log(deltamin/ anneal.Tinit)/(double)anneal.Titns);

    #ifdef DEBUGTRACEFILE
    if (verbosity > 4)
        fclose(fp);
    #endif
} // initialiseConnollyAnnealing

// initialise adaptive annealing (where anneal type = 3)
void initialiseAdaptiveAnnealing(int puno,int spno,double prop,vector<int> &R,vector<spustuff> &pu,vector<sconnections> &connections,
                                 vector<spu>& SM,double cm,vector<sspecies> &spec,int aggexist,sanneal &anneal,int clumptype, int thread, stringstream& logBuffer)
{
    long int i,isamples;
    double sum = 0,sum2 = 0;
    double sigma;
    scost cost;
    double c = 10;  /* An initial temperature acceptance number */
    isamples = 1000; /* Hardwired number of samples to take */

    for (i=0;i<isamples;i++)
    {  /* Generate Random Reserve */
        initialiseReserve(prop, pu, R, rngEngine);
        /* Score Random reserve */
        computeReserveValue(puno,spno,R,pu,connections,SM,cm,spec,aggexist,cost,clumptype, thread, logBuffer);
        /* Add Score to Sum */
        sum += cost.total;
        sum2 += cost.total*cost.total;
    } /* Sample space iterations/100 times */

    sigma = sqrt(sum2 - pow(sum/isamples,2))/(isamples-1);

    anneal.Tinit = c * sigma;
    anneal.sigma = sigma;
    anneal.temp = anneal.Tinit;
    anneal.tempold = anneal.temp;
    anneal.sum = 0;
    anneal.sum2 = 0;

    logBuffer << "Tinit " << anneal.Tinit << " Titns " << anneal.Titns << " Tcool " << anneal.Tcool << endl;
} // initialiseAdaptiveAnnealing

// reduce annealing temperature when anneal type = 3
void reduceTemperature(sanneal& anneal)
{
    double omega = 0.7; /* Control parameter */
    double sigmanew,sigmamod;
    double lambda = 0.7; /* control parameter*/

    sigmanew = (anneal.sum2 - pow((anneal.sum/anneal.Tlen),2))/(anneal.Tlen-1);
    sigmamod = (1-omega)*sigmanew + omega * anneal.sigma *(anneal.temp/anneal.tempold);
    anneal.tempold = anneal.temp;
    anneal.temp = exp(-lambda*anneal.temp/sigmamod);
    anneal.sigma = sigmamod;
    anneal.sum = 0;
    anneal.sum2 = 0;
}

// run simulated thermal annealing selection algorithm
void thermalAnnealing(int spno, int puno, vector<sconnections> &connections,vector<int> &R, double cm,
                      vector<sspecies>& spec, vector<spustuff> &pu, vector<spu>& SM, scost &reserve,
                      long int repeats,int irun,string savename,double misslevel,
                      int aggexist,double costthresh, double tpf1, double tpf2,int clumptype, sanneal &anneal, int thread, stringstream& logBuffer)
{
    scost change = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    long int itime = 0,ipu = -1,i,itemp, snapcount = 0,ichanges = 0, iPreviousR,iGoodChange = 0;
    long int iRowCounter, iRowLimit;
    double rTemperature, rThreshold, rThresholdMultiplier;
    string tempname1,tempname2, sRun = to_string(irun);
    FILE *fp = nullptr,*ttfp = nullptr,*Rfp = nullptr;
    string writename;
    uniform_real_distribution<double> float_range(0.0, 1.0);

    logBuffer << "thermalAnnealing start iterations " << anneal.iterations << "\n";
    if (verbosity > 4)
    {
        writeR(0,"after_Annealing_entered",puno,R,pu,fnames);
        writename = fnames.outputdir + "debug_maropt_annealing_" + sRun + ".csv";
        if ((fp = fopen(writename.c_str(),"w"))==NULL)
            displayErrorMessage("cannot create annealing file %s\n",writename.c_str());
        fprintf(fp,"itime,ipu,puid,R,itemp,newR,iGoodChange,changetotal,changecost,changeconnection,changepen,temp\n");
    }

    if (fnames.saveannealingtrace)
    {
        tempname2 = savename + "_anneal_objective" + sRun + ".csv";
        writename = fnames.outputdir + tempname2;
        if ((ttfp = fopen(writename.c_str(),"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename.c_str());

        fprintf(ttfp,"iteration,threshold,dochange,total,pus,cost,connectivity,penalty,shortfall");
        if (fProb1D == 1)
            fprintf(ttfp,",probability1D");
        if (fProb2D == 1)
            fprintf(ttfp,",probability2D");
        fprintf(ttfp,",puindex\n");

        // write iteration zero
        fprintf(ttfp,"%li,%f,%li,%f,%i,%f,%f,%f,%f",
                itime,costthresh,iGoodChange,reserve.total,
                reserve.pus,reserve.cost,reserve.connection,reserve.penalty,reserve.shortfall);
        if (fProb1D == 1)
            fprintf(ttfp,",%f",reserve.probability1D);
        if (fProb2D == 1)
            fprintf(ttfp,",%f",reserve.probability2D);
        fprintf(ttfp,",%li\n",ipu);
        // iteration,threshold,dochange,total,pus,cost,connectivity,penalty,probability

        tempname2 = savename + "_anneal_zones" + sRun + ".csv";
        writename = fnames.outputdir + tempname2;
        if ((Rfp = fopen(writename.c_str(),"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename.c_str());

        fprintf(Rfp,"configuration");
        for (i = 0;i<puno;i++)
            fprintf(Rfp,",%i",pu[i].id);
        fprintf(Rfp,"\n0");

        for (i = 0;i<puno;i++)
            fprintf(Rfp,",%i",R[i]);
        fprintf(Rfp,"\n");

        iRowCounter = 0;
        if (fnames.annealingtracerows == 0)
            iRowLimit = 0;
        else
            iRowLimit = floor(anneal.iterations / fnames.annealingtracerows);
    }

    displayProgress2("  Main thermalAnnealing Section.\n");

    rThreshold = costthresh;
    costthresh = rThreshold * rStartDecMult;
    rTemperature = 1;
    uniform_int_distribution<int> int_range(0, puno-1);

    for (itime = 1;itime<=anneal.iterations;itime++)
    {
        // Choose random pu. If PU is set > 1 then that pu is fixed and cannot be changed.
        ipu = int_range(rngEngine);
        while (R[ipu] > 1) {
            ipu = int_range(rngEngine);
        }

        itemp = R[ipu] == 1 ? -1 : 1;  /* Add or Remove PU ? */
        computeChangeScore(itime,ipu,spno,puno,pu,connections,spec,SM,R,cm,itemp, change,reserve,
                           costthresh,tpf1,tpf2,(double) itime/ (double) anneal.iterations,clumptype, thread);

        /* Need to calculate Appropriate temperature in isGoodChange or another function */
        /* Upgrade temperature */
        if (itime%anneal.Tlen == 0)
        {
            rTemperature = rTemperature * anneal.Tcool;

            if (rTemperature > rStartDecThresh)
                costthresh = rThreshold * rStartDecMult;
            else
            {
                if (rTemperature < rEndDecThresh)
                    costthresh = rThreshold * rEndDecMult;
                else
                {
                    // map costthresh in the space between (rThreshold * rStartDecMult) and (rThreshold * rEndDecMult)
                    rThresholdMultiplier = (rTemperature - rEndDecThresh) / (rStartDecThresh - rEndDecThresh);
                    costthresh = (rEndDecMult + (rThresholdMultiplier * (rStartDecMult - rEndDecMult))) * rThreshold;
                }
            }
            if (anneal.type == 3)
                reduceTemperature(anneal);
            else
                anneal.temp = anneal.temp*anneal.Tcool;

            displayProgress3("time %ld temp %f Complete %ld%% currval %.4f\n",
                             itime,anneal.temp,(int)itime*100/anneal.iterations,reserve.total);

        } /* reduce temperature */

        if (fnames.savesnapsteps && !(itime % fnames.savesnapfrequency))
        {
            tempname2 = savename + "_snap" + sRun + to_string(++snapcount) + getFileNameSuffix(fnames.savesnapchanges);
            writeSolution(puno,R,pu,tempname2,fnames.savesnapsteps,fnames);
        } /* Save snapshot every savesnapfreq timesteps */

        iPreviousR = R[ipu];
        iGoodChange = isGoodChange(change,anneal.temp, float_range);

        if (iGoodChange)
        {
            ++ichanges;
            
            doChange(ipu,puno,R,reserve,change,pu,SM,spec,connections,itemp,clumptype, thread, logBuffer);

            if (fnames.savesnapchanges && !(ichanges % fnames.savesnapfrequency))
            {
                tempname2 = savename + "_snap" + sRun + to_string(++snapcount) + getFileNameSuffix(fnames.savesnapchanges);
                writeSolution(puno,R,pu,tempname2,fnames.savesnapchanges,fnames);
            } /* Save snapshot every savesnapfreq changes */

        } /* Good change has been made */

        if (anneal.type == 3)
        {
            anneal.sum += reserve.total;
            anneal.sum2 += reserve.total*reserve.total;
        } /* Keep track of scores for averaging stuff */

        if (verbosity > 4)
            fprintf(fp,"%li,%li,%i,%li,%li,%i,%li,%f,%f,%f,%f,%f\n",
                    itime,ipu,pu[ipu].id,iPreviousR,itemp,R[ipu],iGoodChange,change.total,change.cost,change.connection,change.penalty,anneal.temp);

        if (fnames.saveannealingtrace)
        {
            iRowCounter++;
            if (iRowCounter > iRowLimit)
                iRowCounter = 1;

            if (iRowCounter == 1)
            {
                fprintf(Rfp,"%li",itime);

                fprintf(ttfp,"%li,%f,%li,%f,%i,%f,%f,%f,%f",
                        itime,costthresh,iGoodChange,reserve.total,
                        reserve.pus,reserve.cost,reserve.connection,reserve.penalty,reserve.shortfall);
                if (fProb1D == 1)
                    fprintf(ttfp,",%f",reserve.probability1D);
                if (fProb2D == 1)
                    fprintf(ttfp,",%f",reserve.probability2D);
                fprintf(ttfp,",%li\n",ipu);
                // iteration,threshold,dochange,total,pus,cost,connectivity,penalty,probability
                for (i = 0;i<puno;i++)
                    fprintf(Rfp,",%i",R[i]);

                fprintf(Rfp,"\n");
            }
        }

    } /* Run Through Annealing */

    /** Post Processing  **********/
    if (aggexist)
        ClearClumps(spno,spec,pu,SM, thread);

    if (verbosity > 4)
        fclose(fp);

    if (fnames.saveannealingtrace)
    {
        fclose(ttfp);
        fclose(Rfp);
    }
} // thermalAnnealing

void quantumAnnealing(int spno, int puno, vector<sconnections> &connections,vector<int> &R, double cm,
                      vector<sspecies>& spec, vector<spustuff> &pu, vector<spu>& SM, scost &change, scost &reserve,
                      long int repeats,int irun,string savename,double misslevel,
                      int aggexist,double costthresh, double tpf1, double tpf2,int clumptype, sanneal &anneal, int thread)
{
    long int itime,i,j,itemp=0,snapcount,ichanges = 0, iGoodChange;
    long int iRowCounter, iRowLimit, iFluctuationCount;
    double rFluctuationMagnitude, rThreshold, rThresholdMultiplier,
    rAcceptanceProbability;
    string tempname1,tempname2, sRun = to_string(irun);
    FILE *fp = nullptr, *ttfp = nullptr,*Rfp = nullptr;
    string writename, sDecayType;
    vector<int> PUChosen;
    long int iTests = 0, iIterations;
    uniform_real_distribution<double> float_range(0.0, 1.0);

    if (iQADECAYTYPE == 0)
        sDecayType = "EXPONENTIAL";
    else
        sDecayType = "SIGMOIDAL";

    appendTraceFile("quantumAnnealing start iterations %ld decay type %s proportion %f decay A %f decay B %f acceptance probability %f saveannealingtrace %i\n",
                          anneal.iterations,sDecayType.c_str(),rQAPROP,rQADECAY,rQADECAYB,rQAACCPR,fnames.saveannealingtrace);
    if (verbosity > 4)
    {
        writeR(0,"after_Annealing_entered",puno,R,pu,fnames);
        writename = fnames.outputdir + "debug_maropt_annealing_" + sRun + ".csv";
        if ((fp = fopen(writename.c_str(),"w"))==NULL)
            displayErrorMessage("cannot create annealing file %s\n",writename.c_str());
        fprintf(fp,"itime,ipu,puid,R,itemp,newR,iGoodChange,changetotal,changecost,changeconnection,changepen,temp\n");
    }

    if (fnames.saveannealingtrace)
    {
        tempname2 = savename + "_anneal_objective" + sRun + ".csv";
        writename = fnames.outputdir + tempname2;
        if ((ttfp = fopen(writename.c_str(),"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename.c_str());
        fprintf(ttfp,"iteration,threshold,dochange,total,pus,cost,connectivity,penalty");
        if (fProb1D == 1)
            fprintf(ttfp,",probability1D");
        if (fProb2D == 1)
            fprintf(ttfp,",probability2D");
        fprintf(ttfp,",Fmag,Fcount\n");

        tempname2 = savename + "_anneal_zones" + sRun + ".csv";
        writename = fnames.outputdir + tempname2;
        if ((Rfp = fopen(writename.c_str(),"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename.c_str());

        fprintf(Rfp,"configuration");
        for (i = 0;i<puno;i++)
            fprintf(Rfp,",%i",pu[i].id);
        fprintf(Rfp,"\n0");

        for (i = 0;i<puno;i++)
            fprintf(Rfp,",%i",R[i]);
        fprintf(Rfp,"\n");

        iRowCounter = 0;
        if (fnames.annealingtracerows == 0)
            iRowLimit = 0;
        else
            iRowLimit = floor(anneal.iterations / fnames.annealingtracerows);
    }

    displayProgress2("  Main quantumAnnealing Section.\n");

    rThreshold = costthresh;
    costthresh = rThreshold * rStartDecMult;
    rAcceptanceProbability = rQAACCPR; // 1% probability of acceptance of bad moves
    uniform_int_distribution<int> int_range(0, puno-1);

    for (itime = 1;itime<=anneal.iterations;itime++)
    {
        if (iQADECAYTYPE == 0)
        {
            // exponential decay
            rFluctuationMagnitude = exp(-1 * itime / rQADECAY);
        }
        else
        {
            // sigmoidal decay
            rFluctuationMagnitude = 1 / (1 + exp((itime / rQADECAY) - rQADECAYB));
        }
        iFluctuationCount = floor(rFluctuationMagnitude * puno * rQAPROP);

        #ifdef DEBUG_QA
        appendTraceFile("quantumAnnealing rFluctuationMagnitude %f iFluctuationCount %i\n",
                              rFluctuationMagnitude,iFluctuationCount);
        #endif

        if (iFluctuationCount > 0) // we continue if fluctuations are greater than zero
        {
            // we propose to flip the bits on iFluctuationCount PU's
            iTests += iFluctuationCount;
            PUChosen.assign(puno, 0);

            for (i = 0;i<iFluctuationCount;i++)
            {
                do
                {
                    j = int_range(rngEngine);

                    #ifdef DEBUG_QA
                    appendTraceFile("quantumAnnealing j %i PUChosen[j] %i R[j] %i \n",j,PUChosen[j],R[j]);
                    #endif
                }
                while ((PUChosen[j] > 0) || (R[j] > 1));
                // select PU's at random that are not already chosen or locked

                #ifdef DEBUG_QA
                appendTraceFile("quantumAnnealing chose ipu %i\n",j);
                #endif

                PUChosen[j] = 1;
            }

            // compute objective function score with these bits flipped
            computeQuantumChangeScore(spno,puno,pu,connections,spec,SM,R,cm,change,reserve,
                                      costthresh,tpf1,tpf2,(double) itime/ (double) anneal.iterations,
                                      clumptype,iFluctuationCount,PUChosen, thread);

            // we only accept good changes
            if (fnames.savesnapsteps && !(itime % fnames.savesnapfrequency))
            { // Save snapshot every savesnapfreq timesteps
                tempname2 = savename + "_snap" + sRun + to_string(++snapcount) + getFileNameSuffix(fnames.savesnapchanges);
                writeSolution(puno,R,pu,tempname2,fnames.savesnapsteps,fnames);
            }
            if (isGoodQuantumChange(change,rAcceptanceProbability,float_range)==1)
            { // Save snapshot every savesnapfreq changes
                iGoodChange = 1;

                ++ichanges;
                doQuantumChange(puno,R,reserve,change,pu,SM,spec,connections,clumptype,iFluctuationCount,PUChosen, thread);
                if (fnames.savesnapchanges && !(ichanges % fnames.savesnapfrequency))
                {
                    tempname2 = savename + "_snap" + sRun + to_string(++snapcount) + getFileNameSuffix(fnames.savesnapchanges);
                    writeSolution(puno,R,pu,tempname2,fnames.savesnapchanges,fnames);
                }
            } /* Good change has been made */
            else
                iGoodChange = 0;

            if (anneal.type == 3)
            { // Keep track of scores for averaging stuff
                anneal.sum += reserve.total;
                anneal.sum2 += reserve.total*reserve.total;
            }

            if (verbosity > 4)
                fprintf(fp,"%li,%li,%li,%f,%f,%f,%f,%f\n",
                        itime,itemp,iGoodChange,change.total,change.cost,change.connection,change.penalty,anneal.temp);

            if (fnames.saveannealingtrace)
            {
                iRowCounter++;
                if (iRowCounter > iRowLimit)
                    iRowCounter = 1;

                if (iRowCounter == 1)
                {
                    fprintf(Rfp,"%li",itime);

                    fprintf(ttfp,"%li,%f,%li,%f,%i,%f,%f,%f\n",
                            itime,costthresh,iGoodChange,reserve.total,
                            reserve.pus,reserve.cost,reserve.connection,reserve.penalty);
                    if (fProb1D == 1)
                        fprintf(ttfp,",%f",reserve.probability1D);
                    if (fProb2D == 1)
                        fprintf(ttfp,",%f",reserve.probability2D);
                    fprintf(ttfp,",%f,%li\n",rFluctuationMagnitude,iFluctuationCount);
                    // iteration,threshold,dochange,total,pus,cost,connectivity,penalty,probability

                    for (i = 0;i<puno;i++)
                        fprintf(Rfp,",%i",R[i]);

                    fprintf(Rfp,"\n");
                }
            }
        } else {
            // force algorithm to drop out of iterations loop
            iIterations = itime - 1;
            itime = anneal.iterations;
        }
    } /* Run Through Annealing */

    /** Post Processing  **********/
    if (aggexist)
        ClearClumps(spno,spec,pu,SM,thread);

    #ifdef DEBUGTRACEFILE
    if (verbosity > 4)
        fclose(fp);
    #endif

    if (fnames.saveannealingtrace)
    {
        fclose(ttfp);
        fclose(Rfp);
    }
    appendTraceFile("quantumAnnealing end iterations %ld tests %li\n",iIterations,iTests);
} // quantumAnnealing

// marxan is running as a secondary process and has finished. create a sync file so the calling software will know marxan has finished creating the output files
// secondaryExit does not deliver a message prior to exiting, but creates a file so C-Plan/Zonae Cogito/etc knows marxan has exited
void secondaryExit(void)
{
    writeSecondarySyncFile();
}

// iteratively improves a planning unit solutions
// a descent algorithm un-reserves planning units that don't have a negative value when removed
void iterativeImprovement(int puno,int spno,vector<spustuff> &pu, vector<sconnections> &connections,
                          vector<sspecies> &spec,vector<spu>& SM,vector<int> &R, double cm,
                          scost &reserve, scost &change,double costthresh,double tpf1, double tpf2,
                          int clumptype,int irun,string savename, int thread, stringstream& logBuffer)
{
    int puvalid =0,i,j,ipu=0,imode,ichoice, iRowCounter, iRowLimit;
    vector<int> iimparray;
    double debugfloat;
    string tempname2, sRun = to_string(irun);
    FILE *ttfp = nullptr,*Rfp = nullptr;
    string writename;

    logBuffer << "iterativeImprovement start\n";
    // counting pu's we need to test
    for (i=0;i<puno;i++)
    {
        if ((R[i] < 2) && (pu[i].status < 2))
            puvalid++;
    }
    logBuffer << "iterativeImprovement puvalid " << puvalid << "\n";

    if (fnames.saveitimptrace)
    {
        tempname2 = savename + "_itimp_objective" + sRun + ".csv";
        writename = fnames.outputdir + tempname2;
        if ((ttfp = fopen(writename.c_str(),"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename.c_str());
        fprintf(ttfp,"improvement,total,pus,cost,connection,penalty,change total\n");

        tempname2 = savename + "_itimp_zones" + sRun + ".csv";
        writename = fnames.outputdir + tempname2;
        if ((Rfp = fopen(writename.c_str(),"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename.c_str());

        fprintf(Rfp,"configuration");
        for (i = 0;i<puno;i++)
            fprintf(Rfp,",%i",pu[i].id);
        fprintf(Rfp,"\n0");

        for (i = 0;i<puno;i++)
            fprintf(Rfp,",%i",R[i]);
        fprintf(Rfp,"\n");

        iRowCounter = 0;
        if (fnames.itimptracerows == 0)
            iRowLimit = 0;
        else
            iRowLimit = floor(puvalid / fnames.itimptracerows);
    }

    if (puvalid > 0)
    {
        iimparray.resize(puvalid);

        for (i=0;i<puno;i++)
        {
            if ((R[i] < 2) && (pu[i].status < 2))
            {
                iimparray[ipu] = i;
                ipu++;
            }
        }

        logBuffer << "iterativeImprovement after array init\n";

        // shuffle iimp array
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(iimparray.begin(), iimparray.end(), std::default_random_engine(seed));

        /***** Doing the improvements ****/
        for (i=0;i<puvalid;i++)
        {
            ichoice = iimparray[i];

            if ((R[ichoice] < 2) && (pu[ichoice].status < 2))
            {
               imode = R[ichoice] == 1 ? -1 : 1;
               computeChangeScore(-1,ichoice,spno,puno,pu,connections,spec,SM,R,cm,imode,change,reserve,
                                  costthresh,tpf1,tpf2,1,clumptype, thread);
               if (change.total < 0)
               {
                  displayProgress2("It Imp has changed %i with change value %lf \n",ichoice,change.total);
                  doChange(ichoice,puno,R,reserve,change,pu,SM,spec,connections,imode,clumptype, thread, logBuffer);
               }   // I've just made a good change
            }

            if (fnames.saveitimptrace)
            {
               iRowCounter++;
               if (iRowCounter > iRowLimit)
                  iRowCounter = 1;

               if (iRowCounter == 1)
               {
                  fprintf(Rfp,"%i",i);

                  fprintf(ttfp,"%i,%f,%i,%f,%f,%f,%f\n"
                          ,i,reserve.total
                          ,reserve.pus,reserve.cost,reserve.connection,reserve.penalty,
                          change.total); // i,costthresh,pus,cost,connection,penalty

                  for (j = 0;j<puno;j++)
                      fprintf(Rfp,",%i",R[j]);
               }
            }
        } // no untested PUs left
    }

    if (fnames.saveitimptrace)
    {
        fclose(ttfp);
        fclose(Rfp);
    }

    appendTraceFile("iterativeImprovement end\n");
} // iterativeImprovement

} // namespace marxan

// main needs to be defined outside of namespace
int main(int argc, char* argv[]) {
    std::string sInputFileName = "input.dat";

    // store application name
    marxan::sApplicationPathName = argv[0];

    if (argc > 1)
    {
        // handle the program options
        marxan::handleOptions(argc,argv,sInputFileName);
    }

    if (marxan::executeMarxan(sInputFileName)) // Calls the main annealing unit
    {
        if (marxan::marxanIsSecondary == 1)
        {
            marxan::secondaryExit();
        }

        return 1;
    }  // Abnormal Exit

    if (marxan::marxanIsSecondary == 1)
    {
        marxan::secondaryExit();
    }

    return 0;
}