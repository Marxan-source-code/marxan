// C++ code for Marxan
// version 2.3 introduced multiple connectivity files and their associated weighting file
// version 2.4.3 introduced 1D and 2D probability
// version 3.0.0 is refactoring of code in 2019

#include <algorithm>
#include <chrono>
#include <ctime>
#include <cfloat>
#include <iostream>
#include <omp.h>
#include <chrono> 

// load the required function definition modules
#include "utils.hpp"
#include "algorithms.hpp"
#include "computation.hpp"
#include "clumping.hpp"
#include "anneal.hpp"
#include "heuristics.hpp"
#include "probability.hpp"
#include "input.hpp"
#include "output.hpp"

#include "score_change.hpp"
#include "iterative_improvement.hpp"
#include "quantum_annealing.hpp"
#include "thermal_annealing.hpp"
#include "hill_climbing.hpp"

#include "marxan.hpp"

#include "defines.hpp"

namespace marxan {
    using namespace algorithms;
    using namespace utils;

    // Initialize global constants
    double delta;
    long int RandSeed1;
    int iMemoryUsed = 0;
    int fSpecPROPLoaded = 0;
    int iProbFieldPresent = 0;
    int iOptimisationCalcPenalties = 1;
    int savelog;
    int verbosity = 0;
    int asymmetricconnectivity = 0;
    string sVersionString = "Marxan v 4.0.5";
    string sMarxanWebSite = "https://marxansolutions.org/";
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
    chrono::high_resolution_clock::time_point startTime;

    double rProbabilityWeighting = 1;
    double rStartDecThresh = 0.7, rEndDecThresh = 0.3, rStartDecMult = 3, rEndDecMult = 1;

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
    void executeRunLoop(long int repeats, int puno, int spno, double cm, int aggexist, double prop, int clumptype, double misslevel,
        string savename, double costthresh, double tpf1, double tpf2, int heurotype, srunoptions& runoptions,
        int itimptype, vector<int>& sumsoln, rng_engine& rngEngineGlobal)
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
        displayProgress1("Runs will be printed as they complete, and may not be in order due to parallelisation.\n");
        //create seeds for local rng engines
        vector<unsigned int> seeds(repeats);
        for (int run_id = 1; run_id <= repeats; run_id++)
            seeds[run_id - 1] = rngEngineGlobal();

        bool quitting_loop = false;    

        #pragma omp parallel for schedule(dynamic)
        for (int run_id = 1; run_id <= repeats; run_id++)
        {
            if(quitting_loop)
                continue; //skipping iterations. It is not allowed to break or throw out omp for loop.

            // Create run specific structures
            int thread = omp_get_thread_num();
            rng_engine rngEngine(seeds[run_id - 1]);
            string tempname2;
            stringstream appendLogBuffer; // stores the trace file log
            stringstream runConsoleOutput; // stores the console message for the run. This is needed for more organized printing output due to multithreading.
            sanneal anneal = anneal_global;
            scost reserve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            scost change = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            vector<int> R(puno);

            vector<sspecies> spec = specGlobal; // make local copy of original spec

            vector<spu_out> SM_out; // make local copy output part of SMGlobal.
            //if (aggexist)
            SM_out.resize(SMGlobal.size());

            appendLogBuffer << "\n Start run loop run " << run_id << endl;
            try {

                if (runoptions.ThermalAnnealingOn)
                {
                    // Annealing Setup
                    if (anneal.type == 2)
                    {
                        appendLogBuffer << "before initialiseConnollyAnnealing run " << run_id << endl;

                        initialiseConnollyAnnealing(puno, spno, pu, connections, spec, SMGlobal, SM_out, cm, anneal, aggexist, R, prop, clumptype, run_id, appendLogBuffer, rngEngine);

                        appendLogBuffer << "after initialiseConnollyAnnealing run " << run_id << endl;
                    }

                    if (anneal.type == 3)
                    {
                        appendLogBuffer << "before initialiseAdaptiveAnnealing run " << run_id << endl;

                        initialiseAdaptiveAnnealing(puno, spno, prop, R, pu, connections, SMGlobal, SM_out, cm, spec, aggexist, anneal, clumptype, appendLogBuffer, rngEngine);

                        appendLogBuffer << "after initialiseAdaptiveAnnealing run " << run_id << endl;
                    }

                    if (verbosity > 1) runConsoleOutput << "\nRun " << run_id << ": Using Calculated Tinit = " << anneal.Tinit << " Tcool = " << anneal.Tcool << "\n";
                    anneal.temp = anneal.Tinit;
                }

                appendLogBuffer << "before computeReserveValue run " << run_id << endl;

                initialiseReserve(prop, pu, R, rngEngine); // Create Initial Reserve
                SpeciesAmounts(spno, puno, spec, pu, SMGlobal, R, clumptype); // Re-added this from v2.4 because spec amounts need to be refreshed when initializing

                if (aggexist)
                    ClearClumps(spno, spec, pu, SMGlobal, SM_out);

                computeReserveValue(puno, spno, R, pu, connections, SMGlobal, SM_out, cm, spec, aggexist, reserve, clumptype, appendLogBuffer);

                appendLogBuffer << "after computeReserveValue run " << run_id << endl;

                if (verbosity > 1)
                {
                    runConsoleOutput << "Run " << run_id << " Init: " << displayValueForPUs(puno, spno, R, reserve, spec, misslevel).str();
                }
                if (verbosity > 5)
                {
                    displayTimePassed(startTime);
                }

                if (runoptions.ThermalAnnealingOn)
                {
                    appendLogBuffer << "before thermalAnnealing run " << run_id << endl;

                    thermalAnnealing(spno, puno, connections, R, cm, spec, pu, SMGlobal, SM_out, reserve,
                        repeats, run_id, savename, misslevel,
                        aggexist, costthresh, tpf1, tpf2, clumptype, anneal, appendLogBuffer, rngEngine);


                    if (verbosity > 1)
                    {
                        computeReserveValue(puno, spno, R, pu, connections, SMGlobal, SM_out, cm, spec, aggexist, reserve, clumptype, appendLogBuffer);
                        runConsoleOutput << "Run " << run_id << " ThermalAnnealing: " << displayValueForPUs(puno, spno, R, reserve, spec, misslevel).str();
                    }

                    appendLogBuffer << "after thermalAnnealing run " << run_id << endl;
                }

                if (runoptions.HillClimbingOn)
                {
                    appendLogBuffer << "before hill climbing run " << run_id << endl;

                    hill_climbing( puno, spno,  pu, connections, spec,  SMGlobal, SM_out, R,  cm, reserve, costthresh, tpf1, tpf2,
                        clumptype,  run_id, anneal.iterations, savename, appendLogBuffer, rngEngine);

                    if (verbosity > 1)
                    {
                        computeReserveValue(puno, spno, R, pu, connections, SMGlobal, SM_out, cm, spec, aggexist, reserve, clumptype, appendLogBuffer);
                        runConsoleOutput << "Run " << run_id << " Hill climbing: " << displayValueForPUs(puno, spno, R, reserve, spec, misslevel).str();

                    }

                    appendLogBuffer << "after hill climbing run " << run_id << endl;
                }

                if (runoptions.TwoStepHillClimbingOn)
                {
                    appendLogBuffer << "before two step hill climbing run " << run_id << endl;

                    hill_climbing_two_steps( puno, spno,  pu, connections, spec,  SMGlobal, SM_out, R,  cm, reserve, costthresh, tpf1, tpf2,
                        clumptype,  run_id, anneal.iterations, savename, appendLogBuffer, rngEngine);

                    if (verbosity > 1)
                    {
                        computeReserveValue(puno, spno, R, pu, connections, SMGlobal, SM_out, cm, spec, aggexist, reserve, clumptype, appendLogBuffer);
                        runConsoleOutput << "Run " << run_id << " Two Step Hill climbing: " << displayValueForPUs(puno, spno, R, reserve, spec, misslevel).str();

                    }

                    appendLogBuffer << "after hill climbing run " << run_id << endl;
                }

                if (runoptions.HeuristicOn)
                {
                    appendLogBuffer << "before Heuristics run " << run_id << endl;

                    Heuristics(spno, puno, pu, connections, R, cm, spec, SMGlobal, SM_out, reserve,
                        costthresh, tpf1, tpf2, heurotype, clumptype, appendLogBuffer, rngEngine);

                    if (verbosity > 1 && runoptions.ItImpOn)
                    {
                        computeReserveValue(puno, spno, R, pu, connections, SMGlobal, SM_out, cm, spec, aggexist, reserve, clumptype, appendLogBuffer);
                        runConsoleOutput << "Run " << run_id << "  Heuristic: " << displayValueForPUs(puno, spno, R, reserve, spec, misslevel).str();
                    }

                    appendLogBuffer << "after Heuristics run " << run_id << endl;
                }

                if (runoptions.ItImpOn)
                {
                    appendLogBuffer << "before iterativeImprovement run " << run_id << endl;

                    iterativeImprovement(puno, spno, pu, connections, spec, SMGlobal, SM_out, R, cm,
                        reserve, change, costthresh, tpf1, tpf2, clumptype, run_id, savename, appendLogBuffer, rngEngine);

                    if (itimptype == 3)
                        iterativeImprovement(puno, spno, pu, connections, spec, SMGlobal, SM_out, R, cm,
                            reserve, change, costthresh, tpf1, tpf2, clumptype, run_id, savename, appendLogBuffer, rngEngine);

                    appendLogBuffer << "after iterativeImprovement run " << run_id << endl;

                    if (aggexist)
                        ClearClumps(spno, spec, pu, SMGlobal, SM_out);

                    if (verbosity > 1)
                    {
                        computeReserveValue(puno, spno, R, pu, connections, SMGlobal, SM_out, cm, spec, aggexist, reserve, clumptype, appendLogBuffer);
                        runConsoleOutput << "Run " << run_id << " Iterative Improvement: " << displayValueForPUs(puno, spno, R, reserve, spec, misslevel).str();
                    }

                } // Activate Iterative Improvement

                appendLogBuffer << "before file output run " << run_id << endl;
                string fileNumber = utils::intToPaddedString(run_id, 5);
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
                    summaries[run_id - 1] = computeSummary(puno, spno, R, spec, reserve, run_id, misslevel, fnames.savesum);
                }

                // Print results from run.
                displayProgress1(runConsoleOutput.str());

                // compute and store objective function score for this reserve system
                computeReserveValue(puno, spno, R, pu, connections, SMGlobal, SM_out, cm, spec, aggexist, change, clumptype, appendLogBuffer);

                // remember the bestScore and bestRun
                if (change.total < bestScore)
                {
                    omp_set_lock(&bestR_write_lock);
                    // After locking, do another check in case bestScore has changed
                    if (change.total < bestScore) {
                        // this is best run so far
                        bestScore = change.total;
                        bestRun = run_id;
                        // store bestR
                        bestR = R;
                        bestRunString = runConsoleOutput.str();
                        bestSpec = spec; // make a copy of best spec.
                    }
                    omp_unset_lock(&bestR_write_lock);
                }

                if (fnames.savesolutionsmatrix)
                {
                    appendLogBuffer << "before appendSolutionsMatrix savename " << run_id << endl;
                    tempname2 = savename + "_solutionsmatrix" + getFileNameSuffix(fnames.savesolutionsmatrix);

                    omp_set_lock(&solution_matrix_append_lock);
                    appendSolutionsMatrix(run_id, puno, R, tempname2, fnames.savesolutionsmatrix, fnames.solutionsmatrixheaders);
                    omp_unset_lock(&solution_matrix_append_lock);

                    appendLogBuffer << "after appendSolutionsMatrix savename " << run_id << endl;
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
                    ClearClumps(spno, spec, pu, SMGlobal, SM_out);

            }
            catch (const exception& e) {
                // On exceptions, append exception to log file in addition to existing buffer. 
                displayProgress1(runConsoleOutput.str());
                appendLogBuffer << "Exception occurred on run " << run_id << ": " << e.what() << endl;
                displayProgress1(appendLogBuffer.str());
                appendTraceFile(appendLogBuffer.str());
                //cannot throw or break out of omp loop
                quitting_loop = true;
                continue;
            }

            appendLogBuffer << "after file output run " << run_id << endl;
            appendLogBuffer << "end run " << run_id << endl;

            appendTraceFile(appendLogBuffer.str());

            if (marxanIsSecondary == 1)
                writeSecondarySyncFileRun(run_id);

            if (verbosity > 1)
            {
                stringstream done_message;
                done_message << "Run " << run_id << " is finished (out of " << repeats << "). ";
                #pragma omp critical
                {
                    displayProgress1(done_message.str());
                    displayTimePassed(startTime);
                }
            }
        }

        if(quitting_loop)
        {
            displayProgress1("\nRuns were aborted due to error.\n"); 
            throw runtime_error("Runs were aborted due to error.\n");
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
        long int repeats;
        int puno, spno, gspno;
        double cm, prop;
        int heurotype, clumptype, itimptype;
        string savename, tempname2;
        double misslevel;
        int iseed = time(NULL), seedinit;
        int aggexist = 0, sepexist = 0;
        vector<int> R_CalcPenalties;
        vector<int> sumsoln;
        double costthresh, tpf1, tpf2;
        long int itemp;
        int isp;
        int maxThreads = omp_get_max_threads();
        srunoptions runoptions;


        displayStartupMessage();
        startTime = chrono::high_resolution_clock::now(); // set program start time.

        readInputOptions(cm, prop, anneal_global,
            iseed, repeats, savename, fnames, sInputFileName,
            runoptions, misslevel, heurotype, clumptype, itimptype, verbosity,
            costthresh, tpf1, tpf2);

        sTraceFileName = savename + "_TraceFile.txt";
        createTraceFile();
        appendTraceFile("%s begin execution\n\n", sVersionString.c_str());
        appendTraceFile("LoadOptions\n");

#ifdef DEBUGCHECKCHANGE
        createDebugFile("debug_MarOpt_CheckChange.csv", "ipu,puid,R,total,cost,connection,penalty,threshpen,probability\n", fnames);
#endif

#ifdef DEBUGCHANGEPEN
        createDebugFile("debug_MarOpt_ChangePen.csv", "ipu,puid,isp,spid,penalty,target,targetocc,occurrence,sepnum,amount,newamount,fractionAmount\n", fnames);
#endif

#ifdef DEBUGCALCPENALTIES
        createDebugFile("debug_MarZone_CalcPenalties.csv", "\n", fnames);
#endif

        if (fnames.savelog)
        {
            tempname2 = savename + "_log.dat";
            createLogFile(fnames.savelog, tempname2);
        }

        delta = 1e-14;  // This would more elegantly be done as a constant

        // init rng engine
        rng_engine rngEngine(iseed);
        RandSeed1 = iseed;
        seedinit = iseed;

        appendTraceFile("RandSeed iseed %i RandSeed1 %li\n", iseed, RandSeed1);

        // read the data files
        displayProgress1("\nEntering in the data files \n");
        displayProgress3("    Reading in the Planning Unit names \n");
        appendTraceFile("before readPlanningUnits\n");

        itemp = readPlanningUnits(puno, pu, fnames);

        appendTraceFile("after readPlanningUnits\n");
        if (iProbFieldPresent == 1)
            appendTraceFile("prob field present\n");
        else
            appendTraceFile("prob field not present\n");

#ifdef DEBUG_PROB1D
        if (iProbFieldPresent == 1)
            writeProbData(puno, pu, fnames);
#endif

        displayProgress1("   There are %i Planning units.\n  %i Planning Unit names read in \n", puno, itemp);
        displayProgress3("    Reading in the species file \n");
        appendTraceFile("before readSpecies\n");

        itemp = readSpecies(spno, specGlobal, fnames);

        appendTraceFile("after readSpecies\n");
        displayProgress1("  %i species read in \n", itemp);
        appendTraceFile("before build search arrays\n");

        // create the fast lookup tables for planning units and species names
        populateLookupTables(puno, spno, pu, specGlobal, PULookup, SPLookup);

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

            itemp = readConnections(puno, connections, pu, PULookup, fnames);

            appendTraceFile("after readConnections\n");
            if (asymmetricconnectivity)
            {
                appendTraceFile("Asymmetric connectivity is on.\n");
                writeAsymmetricConnectionFile(puno, connections, pu, fnames);
            }
            if (fOptimiseConnectivityIn)
                appendTraceFile("Optimising 'Connectivity In'.\n");
        }
        displayProgress1("  %i connections entered \n", itemp);
        if (asymmetricconnectivity)
            displayProgress1("  Asymmetric connectivity is on.\n");
        if (fOptimiseConnectivityIn)
            displayProgress1("  Optimising 'Connectivity In'.\n");

        displayProgress3("    Reading in the Planning Unit versus Species File \n");
        appendTraceFile("before readSparseMatrix\n");

        readSparseMatrix(SMGlobal, puno, spno, pu, PULookup, SPLookup, fnames);

        appendTraceFile("after readSparseMatrix\n");
        if (fProb2D == 1)
            appendTraceFile("Prob2D is on\n");
        else
            appendTraceFile("Prob2D is off\n");

#ifdef DEBUG_PROB2D
        writeSparseMatrix(iSparseMatrixFileLength, puno, pu, spec, SM, fnames);
#endif

        if (fnames.saverichness)
        {
            tempname2 = savename + "_richness" + getFileNameSuffix(fnames.saverichness);
            writeRichness(puno, pu, tempname2, fnames.saverichness);
        }

        if (!fnames.matrixspordername.empty())
        {
            appendTraceFile("before readSparseMatrixSpOrder\n");

            displayWarningMessage("input.dat option: MATRIXSPORDERNAME is no longer needed, however the supplied file will still be read and processed. Please refer to the version 4 feature changelog.md.");
            readSparseMatrixSpOrder(SMsporder, puno, spno, PULookup, SPLookup, specGlobal, fnames);

            appendTraceFile("after readSparseMatrixSpOrder\n");
        }

        appendTraceFile("before process block definitions\n");

        if (!fnames.blockdefname.empty())
        {
            displayProgress1("    Reading in the Block Definition File \n");
            vector<sgenspec> gspec; // declare within this scope of usage
            readSpeciesBlockDefinition(gspno, gspec, fnames);
            setBlockDefinitions(gspno, spno, puno, gspec, specGlobal, pu, SMGlobal);
        }

        setDefaultTargets(specGlobal);
        appendTraceFile("after process block definitions\n");

        appendTraceFile("before computeTotalAreas\n");
        computeTotalAreas(puno, spno, pu, specGlobal, SMGlobal);
        appendTraceFile("after computeTotalAreas\n");

        if (fnames.savetotalareas)
        {
            tempname2 = savename + "totalareas" + getFileNameSuffix(fnames.savetotalareas);
            writeTotalAreas(puno, spno, pu, specGlobal, SMGlobal, tempname2, fnames.savepenalty);
        }

        if (fSpecPROPLoaded > 0)
        {
            appendTraceFile("before computeSpecProp\n");

            // species have prop value specified
            computeSpecProp(spno, specGlobal, puno, pu, SMGlobal);

            appendTraceFile("after computeSpecProp\n");
        }

        displayProgress2("Checking to see if there are aggregating or separating species.\n");
        for (isp = 0; isp < spno; isp++)
        {
            if (specGlobal[isp].target2 > 0)
                aggexist = 1;
            if (specGlobal[isp].sepdistance > 0)
                sepexist = 1;
        }

        if (fnames.savesen)
        {
            appendTraceFile("before writeScenario\n");

            tempname2 = savename + "_sen.dat";
            writeScenario(puno, spno, prop, cm, anneal_global, seedinit, repeats, clumptype,
                runoptions, heurotype, costthresh, tpf1, tpf2, tempname2);

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
            readPenalties(specGlobal, spno, fnames, SPLookup);
            appendTraceFile("after readPenalties\n");
        }

        if (runoptions.CalcPenaltiesOn == 0)
        {
            // if penalties have not been loaded, then stop error message
            if (fUserPenalties == 0)
            {
                appendTraceFile("Data error: CalcPenalties off but no PENALTYNAME specified, exiting.\n");
                displayProgress1("Data error: CalcPenalties off but no PENALTYNAME specified, exiting.\n");
                exit(EXIT_FAILURE);
            }

            // transfer loaded penalties to correct data structrure
            applyUserPenalties(specGlobal);
        }
        else
        {
            vector<spu_out> SM_out; // make local copy output part of SMGlobal.
            //if (aggexist)
            SM_out.resize(SMGlobal.size());


            // we are computing penalties
            if (fnames.matrixspordername.empty())
            {
                appendTraceFile("before CalcPenalties\n");

                // we don't have sporder matrix available, so use slow CalcPenalties method
                itemp = computePenalties(puno, spno, pu, specGlobal, connections, SMGlobal, SM_out, R_CalcPenalties, aggexist, cm, clumptype, rngEngine);

                appendTraceFile("after CalcPenalties\n");
            }
            else
            {
                // we have sporder matrix available, so use optimised CalcPenalties method
                if (iOptimisationCalcPenalties == 1)
                {
                    appendTraceFile("before CalcPenaltiesOptimise\n");

                    itemp = computePenaltiesOptimise(puno, spno, pu, specGlobal, connections, SMGlobal, SM_out, SMsporder, R_CalcPenalties, aggexist, cm, clumptype, rngEngine);

                    appendTraceFile("after CalcPenaltiesOptimise\n");
                }
                else
                {
                    appendTraceFile("before CalcPenalties\n");

                    // we have optimise calc penalties switched off, so use slow CalcPenalties method
                    itemp = computePenalties(puno, spno, pu, specGlobal, connections, SMGlobal, SM_out, R_CalcPenalties, aggexist, cm, clumptype, rngEngine);

                    appendTraceFile("after CalcPenalties\n");
                }
            }
        }

        if (itemp > 0)
            displayProgress("%d species cannot meet target%c.\n", itemp, itemp == 1 ? ' ' : 's');

        if (runoptions.ThermalAnnealingOn)
        {
            displayProgress2("    Calculating temperatures.\n");
            if (!anneal_global.Titns)
                displayErrorMessage("Initial Temperature is set to zero. Fatal Error \n");

            anneal_global.Tlen = anneal_global.iterations / anneal_global.Titns;
            displayProgress2("  Temperature length %ld \n", anneal_global.Tlen);
            displayProgress2("  iterations %lld, repeats %ld \n", anneal_global.iterations, repeats);
        } // Annealing Preprocessing. Should be moved to SetAnnealingOptions

        if (fnames.savepenalty)
        {
            tempname2 = savename + "_penalty" + getFileNameSuffix(fnames.savepenalty);
            writePenalty(spno, specGlobal, tempname2, fnames.savepenalty);

            tempname2 = savename + "_penalty_planning_units" + getFileNameSuffix(fnames.savepenalty);
            writePenaltyPlanningUnits(puno, pu, R_CalcPenalties, tempname2, fnames.savepenalty);
        }

        //free(R_CalcPenalties);

        if (fnames.savespec)
        {
            tempname2 = savename + "_spec.csv";
            writeSpec(spno, specGlobal, tempname2);
        }

        if (fnames.savepu)
        {
            tempname2 = savename + "_pu.csv";
            writePu(puno, pu, tempname2);
        }

        // If we are in a runmode with only CalcPenalties, we stop/exit here gracefully because we are finished.
        if (runoptions.HeuristicOn == 0 && runoptions.ThermalAnnealingOn == 0 && runoptions.HillClimbingOn == 0 && runoptions.ItImpOn == 0 && runoptions.TwoStepHillClimbingOn == 0)
        {
            appendTraceFile("end final file output\n");
            appendTraceFile("\nMarxan end execution\n");
            displayShutdownMessage(startTime);

            if (marxanIsSecondary == 1)
                secondaryExit();

            exit(EXIT_FAILURE);
        }

        if (fnames.savesolutionsmatrix)
        {
            tempname2 = savename + "_solutionsmatrix" + getFileNameSuffix(fnames.savesolutionsmatrix);

#ifdef DEBUG_CLUSTERANALYSIS
            appendTraceFile("before createSolutionsMatrix savename %s\n", savename);
#endif

            createSolutionsMatrix(puno, pu, tempname2, fnames.savesolutionsmatrix, fnames.solutionsmatrixheaders);

#ifdef DEBUG_CLUSTERANALYSIS
            appendTraceFile("after createSolutionsMatrix savename %s\n", savename);
#endif
        }

        if (fProb1D == 1)
        {
            tempname2 = savename + "_ComputeP_AllPUsSelected_1D.csv";
            ComputeP_AllPUsSelected_1D(tempname2, puno, spno, pu, SMGlobal, specGlobal);
        }

        if (fProb2D == 1)
        {
            tempname2 = savename + "_ComputeP_AllPUsSelected_2D.csv";
            ComputeP_AllPUsSelected_2D(tempname2, puno, spno, pu, SMGlobal, specGlobal);
        }

        try
        {
            executeRunLoop(repeats, puno, spno, cm, aggexist, prop, clumptype, misslevel,
                savename, costthresh, tpf1, tpf2, heurotype, runoptions,
                itimptype, sumsoln, rngEngine);
        }
        catch (const exception& e)
        {
            appendTraceFile("Error during main loop.\n");
            exit(EXIT_FAILURE);
        }

        appendTraceFile("before final file output\n");

        if (fnames.savebest)
        {
            tempname2 = savename + "_best" + getFileNameSuffix(fnames.savebest);
            writeSolution(puno, bestR, pu, tempname2, fnames.savebest, fnames);

            appendTraceFile("Best solution is run %i\n", bestRun);
            displayProgress1("\nBest solution is run %i\n", bestRun);
        }

        if (fnames.savespecies && fnames.savebest)
        {
            tempname2 = savename + "_mvbest" + getFileNameSuffix(fnames.savespecies);

            // TODO ADBAI - need to write best spec, NOT spec_global
            writeSpecies(spno, bestSpec, tempname2, fnames.savespecies, misslevel);
        }

        if (fnames.savesumsoln)
        {
            tempname2 = savename + "_ssoln" + getFileNameSuffix(fnames.savesumsoln);
            writeSumSoln(puno, sumsoln, pu, tempname2, fnames.savesumsoln);
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
        appendTraceFile("\nMarxan end execution. Press any key to continue\n");

        return 0;
    } // executeMarxan


    // handle command line parameters for the marxan executable
    void handleOptions(int argc, char* argv[], string sInputFileName)
    {
        if (argc > 4)
        {  // if more than one commandline argument then exit
            displayUsage(argv[0]);
            exit(1);
        }

        for (int i = 1; i < argc; i++)
        { // Deal with all arguments
            if (argv[i][0] == '/' || argv[i][0] == '-')
            {
                switch (argv[i][1])
                {
                case 'C':
                case 'c':
                case 'S':
                case 's':
                    marxanIsSecondary = 1;
                    break;
                default:
                    fprintf(stderr, "unknown option %s\n", argv[i]);
                    break;
                }
            }
            else {
                sInputFileName = argv[i]; /* If not a -option then must be input.dat name */
            }
        }
    }


    // marxan is running as a secondary process and has finished. create a sync file so the calling software will know marxan has finished creating the output files
    // secondaryExit does not deliver a message prior to exiting, but creates a file so C-Plan/Zonae Cogito/etc knows marxan has exited
    void secondaryExit(void)
    {
        writeSecondarySyncFile();
    }


} // namespace marxan

// main needs to be defined outside of namespace
int main(int argc, char* argv[]) {
    std::string sInputFileName = "input.dat";

    // store application name
    marxan::sApplicationPathName = argv[0];

    if (argc > 1)
    {
        // handle the program options
        marxan::handleOptions(argc, argv, sInputFileName);
    }

    try
    {
        if (marxan::executeMarxan(sInputFileName)) // Calls the main annealing unit
        {
            if (marxan::marxanIsSecondary == 1)
            {
                marxan::secondaryExit();
            }

            return 1;
        }  // Abnormal Exit
    }
    catch(const std::exception& e)
    {
        std::cerr << "Error during Marxan execution."  << '\n';
        exit(EXIT_FAILURE);
    }
    

    if (marxan::marxanIsSecondary == 1)
    {
        marxan::secondaryExit();
    }

    return 0;
}