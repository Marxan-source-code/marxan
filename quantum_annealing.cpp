//Implemented QuantumAnnealing algorithm has the following issue:
//computeQuantumChangeScore has to be redesigned beacause
//Total score !=  Sum of scores of individual pu change
//due mutual effect via species cups, boundaries and clumping.


#define DEBUGTRACEFILE
#define PROB2D

// Flags to change to "define" for debugging purposes
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

namespace marxan {

    using namespace algorithms;
    using namespace utils;

    double rQAPROP = 0.5, rQADECAY = 0.0001, rQADECAYB = 0, rQAACCPR = 0;
    int iQADECAYTYPE = 0;

    // determines if the change value for changing status for a set of planning units is good
    // does it stochastically fall below the current acceptance probability?
    int isGoodQuantumChange(struct scost change, double rProbAcceptance, uniform_real_distribution<double>& float_range, rng_engine& rngEngine)
    {
        if (change.total <= 0)
            return 1;
        else
            return (rProbAcceptance > float_range(rngEngine)) ? 1 : 0;
    }


    // compute change in the objective function score for adding or removing a set of planning units
    void computeQuantumChangeScore(int spno, int puno, const vector<spustuff>& pu, const vector<sconnections>& connections,
        vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, const vector<int>& R, double cm,
        scost& change, scost& reserve, double costthresh, double tpf1, double tpf2,
        double timeprop, int clumptype, int iFluctuationCount, const vector<int>& PUChosen)
        // imode = 1 add PU, imode = -1 remove PU
    {

        //computeQuantumChangeScore has to be redesigned beacause
        //Total score !=  Sum of scores of individual pu change
        //due mutual effect via species cups, boundaries and clumping.
        throw runtime_error("computeQuantumChangeScore has to be redesigned beacause Total score !=  Sum of scores of individual pu change/n"); 

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
            change.connection += ConnectionCost2(connections[j], R, imode, 1, cm, asymmetricconnectivity, fOptimiseConnectivityIn);
            if (threshtype == 1)
            {
                tchangeconnection = change.connection;
                tresconnection = reserve.connection;
                change.connection = 0;
                reserve.connection = 0;
            }

            change.penalty += computeChangePenalty(j, puno, spec, pu, SM, SM_out, R, connections, imode, clumptype, change.shortfall);

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


    // change the status of a set of planning units
    void doQuantumChange(int puno, vector<int>& R, scost& reserve, scost& change,
        const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out, vector<sspecies> spec, const vector<sconnections>& connections,
        int clumptype, int iFluctuationCount, const vector<int>& PUChosen)
    {
        // We accept a whole bunch of changes in one, passed in by Quantum annealing.
        vector<sclumps> tempSclumps;
        int i, j, ipu, ism, isp, imode;
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
        for (j = 0; j < iFluctuationCount; j++)
        {
            do
                ipu++;

            while (PUChosen[ipu] < 1);

            imode = R[ipu] == 1 ? -1 : 1;
            R[ipu] = imode == 1 ? 1 : 0;
            reserve.pus += imode;

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
                            AddNewPU(ipu, isp, connections, spec, pu, SM, SM_out, clumptype);
                        }
                        else {
                            RemPu(ipu, isp, connections, spec, pu, SM, SM_out, clumptype);
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
                        {
                            if (spec[isp].amount > -0.0001)
                                spec[isp].amount = 0;
                        }

#ifdef ANNEALING_TEST
                        appendTraceFile("doChange ipu %i isp %i spec.amount %g imode %i\n",
                            ipu, isp, spec[isp].amount, imode);
#endif
                    }

                    if (spec[isp].sepnum > 0) // Count separation but only if it is possible that it has changed
                        if ((imode == 1 && spec[isp].separation < spec[isp].sepnum) || (imode == -1 && spec[isp].separation > 1))
                            spec[isp].separation = CountSeparation2(isp, 0, tempSclumps, puno, R, pu, SM, SM_out, spec, 0);
                }
            }
        }

        reserve.total = reserve.cost + reserve.connection + reserve.penalty + reserve.probability1D + reserve.probability2D;

#ifdef DEBUG_QA
        appendTraceFile("doQuantumChange end\n");
#endif
    } // doQuantumChange


    void quantumAnnealing(int spno, int puno, const vector<sconnections>& connections, vector<int>& R, double cm,
        vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out, scost& change, scost& reserve,
        long int repeats, int irun, string savename, double misslevel,
        int aggexist, double costthresh, double tpf1, double tpf2, int clumptype, sanneal& anneal, rng_engine& rngEngine)
    {
        long int itime, i, j, itemp = 0, snapcount, ichanges = 0, iGoodChange;
        long int iRowCounter, iRowLimit, iFluctuationCount;
        double rFluctuationMagnitude, rThreshold, rThresholdMultiplier,
            rAcceptanceProbability;
        string tempname1, tempname2, sRun = to_string(irun), paddedRun = utils::intToPaddedString(irun, 5);
        FILE* fp = nullptr, * ttfp = nullptr, * Rfp = nullptr;
        string writename, sDecayType;
        vector<int> PUChosen;
        long int iTests = 0, iIterations;
        uniform_real_distribution<double> float_range(0.0, 1.0);

        if (iQADECAYTYPE == 0)
            sDecayType = "EXPONENTIAL";
        else
            sDecayType = "SIGMOIDAL";

        appendTraceFile("quantumAnnealing start iterations %ld decay type %s proportion %f decay A %f decay B %f acceptance probability %f saveannealingtrace %i\n",
            anneal.iterations, sDecayType.c_str(), rQAPROP, rQADECAY, rQADECAYB, rQAACCPR, fnames.saveannealingtrace);
        if (verbosity > 4)
        {
            writeR(0, "after_Annealing_entered", puno, R, pu, fnames);
            writename = fnames.outputdir + "debug_maropt_annealing_" + sRun + ".csv";
            if ((fp = fopen(writename.c_str(), "w")) == NULL)
                displayErrorMessage("cannot create annealing file %s\n", writename.c_str());
            fprintf(fp, "itime,ipu,puid,R,itemp,newR,iGoodChange,changetotal,changecost,changeconnection,changepen,temp\n");
        }

        if (fnames.saveannealingtrace)
        {
            tempname2 = savename + "_anneal_objective" + paddedRun + ".csv";
            writename = fnames.outputdir + tempname2;
            if ((ttfp = fopen(writename.c_str(), "w")) == NULL)
                displayErrorMessage("cannot create threshold trace file %s\n", writename.c_str());
            fprintf(ttfp, "iteration,threshold,dochange,total,pus,cost,connectivity,penalty");
            if (fProb1D == 1)
                fprintf(ttfp, ",probability1D");
            if (fProb2D == 1)
                fprintf(ttfp, ",probability2D");
            fprintf(ttfp, ",Fmag,Fcount\n");

            tempname2 = savename + "_anneal_zones" + paddedRun + ".csv";
            writename = fnames.outputdir + tempname2;
            if ((Rfp = fopen(writename.c_str(), "w")) == NULL)
                displayErrorMessage("cannot create threshold trace file %s\n", writename.c_str());

            fprintf(Rfp, "configuration");
            for (i = 0; i < puno; i++)
                fprintf(Rfp, ",%i", pu[i].id);
            fprintf(Rfp, "\n0");

            for (i = 0; i < puno; i++)
                fprintf(Rfp, ",%i", R[i]);
            fprintf(Rfp, "\n");

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
        uniform_int_distribution<int> int_range(0, puno - 1);

        for (itime = 1; itime <= anneal.iterations; itime++)
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
                rFluctuationMagnitude, iFluctuationCount);
#endif

            if (iFluctuationCount > 0) // we continue if fluctuations are greater than zero
            {
                // we propose to flip the bits on iFluctuationCount PU's
                iTests += iFluctuationCount;
                PUChosen.assign(puno, 0);

                for (i = 0; i < iFluctuationCount; i++)
                {
                    do
                    {
                        j = int_range(rngEngine);

#ifdef DEBUG_QA
                        appendTraceFile("quantumAnnealing j %i PUChosen[j] %i R[j] %i \n", j, PUChosen[j], R[j]);
#endif
                    } while ((PUChosen[j] > 0) || (R[j] > 1));
                    // select PU's at random that are not already chosen or locked

#ifdef DEBUG_QA
                    appendTraceFile("quantumAnnealing chose ipu %i\n", j);
#endif

                    PUChosen[j] = 1;
                }

                // compute objective function score with these bits flipped
                computeQuantumChangeScore(spno, puno, pu, connections, spec, SM, SM_out, R, cm, change, reserve,
                    costthresh, tpf1, tpf2, (double)itime / (double)anneal.iterations,
                    clumptype, iFluctuationCount, PUChosen);

                // we only accept good changes
                if (fnames.savesnapsteps && !(itime % fnames.savesnapfrequency))
                { // Save snapshot every savesnapfreq timesteps
                    tempname2 = savename + "_snap" + paddedRun + "t" + utils::intToPaddedString(++snapcount, 5) + getFileNameSuffix(fnames.savesnapchanges);
                    writeSolution(puno, R, pu, tempname2, fnames.savesnapsteps, fnames);
                }
                if (isGoodQuantumChange(change, rAcceptanceProbability, float_range, rngEngine) == 1)
                { // Save snapshot every savesnapfreq changes
                    iGoodChange = 1;

                    ++ichanges;
                    doQuantumChange(puno, R, reserve, change, pu, SM, SM_out, spec, connections, clumptype, iFluctuationCount, PUChosen);
                    if (fnames.savesnapchanges && !(ichanges % fnames.savesnapfrequency))
                    {
                        tempname2 = savename + "_snap" + paddedRun + "c" + utils::intToPaddedString(++snapcount, 5) + getFileNameSuffix(fnames.savesnapchanges);
                        writeSolution(puno, R, pu, tempname2, fnames.savesnapchanges, fnames);
                    }
                } /* Good change has been made */
                else
                    iGoodChange = 0;

                if (anneal.type == 3)
                { // Keep track of scores for averaging stuff
                    anneal.sum += reserve.total;
                    anneal.sum2 += reserve.total * reserve.total;
                }

                if (verbosity > 4)
                    fprintf(fp, "%li,%li,%li,%f,%f,%f,%f,%f\n",
                        itime, itemp, iGoodChange, change.total, change.cost, change.connection, change.penalty, anneal.temp);

                if (fnames.saveannealingtrace)
                {
                    iRowCounter++;
                    if (iRowCounter > iRowLimit)
                        iRowCounter = 1;

                    if (iRowCounter == 1)
                    {
                        fprintf(Rfp, "%li", itime);

                        fprintf(ttfp, "%li,%f,%li,%f,%i,%f,%f,%f\n",
                            itime, costthresh, iGoodChange, reserve.total,
                            reserve.pus, reserve.cost, reserve.connection, reserve.penalty);
                        if (fProb1D == 1)
                            fprintf(ttfp, ",%f", reserve.probability1D);
                        if (fProb2D == 1)
                            fprintf(ttfp, ",%f", reserve.probability2D);
                        fprintf(ttfp, ",%f,%li\n", rFluctuationMagnitude, iFluctuationCount);
                        // iteration,threshold,dochange,total,pus,cost,connectivity,penalty,probability

                        for (i = 0; i < puno; i++)
                            fprintf(Rfp, ",%i", R[i]);

                        fprintf(Rfp, "\n");
                    }
                }
            }
            else {
                // force algorithm to drop out of iterations loop
                iIterations = itime - 1;
                itime = anneal.iterations;
            }
        } /* Run Through Annealing */

        /** Post Processing  **********/
        if (aggexist)
            ClearClumps(spno, spec, pu, SM, SM_out);

#ifdef DEBUGTRACEFILE
        if (verbosity > 4)
            fclose(fp);
#endif

        if (fnames.saveannealingtrace)
        {
            fclose(ttfp);
            fclose(Rfp);
        }
        appendTraceFile("quantumAnnealing end iterations %ld tests %li\n", iIterations, iTests);
    } // quantumAnnealing


} // namespace marxan

