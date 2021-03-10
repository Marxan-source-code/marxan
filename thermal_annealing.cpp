
#include <algorithm>
#include <chrono>
#include <ctime>
#include <cfloat>
#include <iostream>
#include <omp.h>
#include <chrono> 

// load the required function definition modules
#include "defines.hpp"
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


    // determines if the change value for changing a single planning unit status is good
    // does the change stochastically fall below the current acceptance probability?
    int isGoodChange(const scost& change, double temp, uniform_real_distribution<double>& float_range, rng_engine& rngEngine)
    {
        if (change.total <= 0)
            return 1;
        else
            return (exp(-change.total / temp) > float_range(rngEngine));
    }

    void initialiseConnollyAnnealing(int puno, int spno, const vector<spustuff>& pu, const vector<sconnections>& connections, vector<sspecies>& spec,
        const vector<spu>& SM, vector<spu_out>& SM_out, double cm, sanneal& anneal, int aggexist,
        vector<int>& R, double prop, int clumptype, int irun, stringstream& logBuffer, rng_engine& rngEngine)
    {
        long long i;
        long int ipu, imode, iOldR;
        double deltamin = 0, deltamax = 0;
        double localdelta = 1E-10;
        scost change = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        scost reserve = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

#ifdef DEBUGTRACEFILE
        FILE* fp = nullptr;
        if (verbosity > 4)
        {
            string writename = fnames.outputdir + "debug_maropt_initialiseConnollyAnnealing_" + to_string(irun) + ".csv";
            fp = fopen(writename.c_str(), "w");
            if (fp == NULL)
                displayErrorMessage("cannot create debug_maropt_initialiseConnollyAnnealing file %s\n", writename.c_str());
            fprintf(fp, "i,ipu,puid,old R,imode,R,total,max,min\n");
        }
#endif

#ifdef DEBUG_PROB1D
        logBuffer << "initialiseConnollyAnnealing A - before initialise reserve\n";
#endif

        initialiseReserve(prop, pu, R, rngEngine);

#ifdef DEBUG_PROB1D
        logBuffer << "initialiseConnollyAnnealing B - after initialise reserve\n";
#endif

        if (aggexist)
            ClearClumps(spno, spec, pu, SM, SM_out);

#ifdef DEBUG_PROB1D
        logBuffer << "initialiseConnollyAnnealing C - before compute reserve\n";
#endif

        computeReserveValue(puno, spno, R, pu, connections, SM, SM_out, cm, spec, aggexist, reserve, clumptype, logBuffer);

#ifdef DEBUG_PROB1D
        logBuffer << "initialiseConnollyAnnealing D - after compute reserve\n";
#endif

        std::uniform_int_distribution<int> int_range(0, puno - 1);
        for (i = 1; i <= anneal.iterations / 100; i++)
        {
            ipu = int_range(rngEngine);
            iOldR = R[ipu];
            imode = R[ipu] == 1 ? -1 : 1;

            computeChangeScore(-1, ipu, spno, puno, pu, connections, spec, SM, SM_out, R, cm, imode, change, reserve, 0, 0, 0, 0, clumptype);
            doChange(ipu, puno, R, reserve, change, pu, SM, SM_out, spec, connections, imode, clumptype, logBuffer);
            if (change.total > deltamax)
                deltamax = change.total;
            if (change.total > localdelta && (deltamin < localdelta || change.total < deltamin))
                deltamin = change.total;

            if (verbosity > 4)
                fprintf(fp, "%li,%li,%i,%li,%li,%i,%g,%g,%g\n", i, ipu, pu[ipu].id, iOldR, imode, R[ipu], change.total, deltamax, deltamin);
            // i,ipu,puid,R,imode,iZone,total,max,min
        }  // Run through this bit for iterations/100 times

        anneal.Tinit = deltamax;
        deltamin *= 0.1;
        anneal.Tcool = exp(log(deltamin / anneal.Tinit) / (double)anneal.Titns);

#ifdef DEBUGTRACEFILE
        if (verbosity > 4)
            fclose(fp);
#endif
    } // initialiseConnollyAnnealing

    // initialise adaptive annealing (where anneal type = 3)
    void initialiseAdaptiveAnnealing(int puno, int spno, double prop, vector<int>& R, const vector<spustuff>& pu, const vector<sconnections>& connections,
        const vector<spu>& SM, vector<spu_out>& SM_out, const double cm, vector<sspecies>& spec, int aggexist, sanneal& anneal, int clumptype, stringstream& logBuffer, rng_engine& rngEngine)
    {
        long int i, isamples;
        double sum = 0, sum2 = 0;
        double sigma;
        scost cost;
        double c = 10;  /* An initial temperature acceptance number */
        isamples = 1000; /* Hardwired number of samples to take */

        for (i = 0; i < isamples; i++)
        {  /* Generate Random Reserve */
            initialiseReserve(prop, pu, R, rngEngine);
            /* Score Random reserve */
            computeReserveValue(puno, spno, R, pu, connections, SM, SM_out, cm, spec, aggexist, cost, clumptype, logBuffer);
            /* Add Score to Sum */
            sum += cost.total;
            sum2 += cost.total * cost.total;
        } /* Sample space iterations/100 times */

        sigma = sqrt(sum2 - pow(sum / isamples, 2)) / (isamples - 1);

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
        double sigmanew, sigmamod;
        double lambda = 0.7; /* control parameter*/

        sigmanew = (anneal.sum2 - pow((anneal.sum / anneal.Tlen), 2)) / (anneal.Tlen - 1);
        sigmamod = (1 - omega) * sigmanew + omega * anneal.sigma * (anneal.temp / anneal.tempold);
        anneal.tempold = anneal.temp;
        anneal.temp = exp(-lambda * anneal.temp / sigmamod);
        anneal.sigma = sigmamod;
        anneal.sum = 0;
        anneal.sum2 = 0;
    }

    // run simulated thermal annealing selection algorithm
    void thermalAnnealing(int spno, int puno, const vector<sconnections>& connections, vector<int>& R, double cm,
        vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out, scost& reserve,
        long int repeats, int irun, string savename, double misslevel,
        int aggexist, double costthresh, double tpf1, double tpf2, int clumptype, sanneal& anneal, stringstream& logBuffer, rng_engine& rngEngine)
    {
        scost change = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        long long itime = 0; 
        long int ipu = -1, i, itemp, snapcount = 0, ichanges = 0, iPreviousR, iGoodChange = 0;
        long int iRowCounter, iRowLimit;
        double rTemperature, rThreshold, rThresholdMultiplier;
        string tempname1, tempname2, sRun = to_string(irun), paddedRun = utils::intToPaddedString(irun, 5);
        FILE* fp = nullptr, * ttfp = nullptr, * Rfp = nullptr;
        string writename;
        uniform_real_distribution<double> float_range(0.0, 1.0);

        logBuffer << "thermalAnnealing start iterations " << anneal.iterations << "\n";
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

            fprintf(ttfp, "iteration,threshold,dochange,total,pus,cost,connectivity,penalty,shortfall");
            if (fProb1D == 1)
                fprintf(ttfp, ",probability1D");
            if (fProb2D == 1)
                fprintf(ttfp, ",probability2D");
            fprintf(ttfp, ",puindex\n");

            // write iteration zero
            fprintf(ttfp, "%li,%f,%li,%f,%i,%f,%f,%f,%f",
                itime, costthresh, iGoodChange, reserve.total,
                reserve.pus, reserve.cost, reserve.connection, reserve.penalty, reserve.shortfall);
            if (fProb1D == 1)
                fprintf(ttfp, ",%f", reserve.probability1D);
            if (fProb2D == 1)
                fprintf(ttfp, ",%f", reserve.probability2D);
            fprintf(ttfp, ",%li\n", ipu);
            // iteration,threshold,dochange,total,pus,cost,connectivity,penalty,probability

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

        displayProgress2("  Main thermalAnnealing Section.\n");

        rThreshold = costthresh;
        costthresh = rThreshold * rStartDecMult;
        rTemperature = 1;
        uniform_int_distribution<int> int_range(0, puno - 1);

        for (itime = 1; itime <= anneal.iterations; itime++)
        {
            // Choose random pu. If PU is set > 1 then that pu is fixed and cannot be changed.
            ipu = int_range(rngEngine);
            while (R[ipu] > 1) {
                ipu = int_range(rngEngine);
            }

            itemp = R[ipu] == 1 ? -1 : 1;  /* Add or Remove PU ? */
            computeChangeScore(itime, ipu, spno, puno, pu, connections, spec, SM, SM_out, R, cm, itemp, change, reserve,
                costthresh, tpf1, tpf2, (double)itime / (double)anneal.iterations, clumptype);

            /* Need to calculate Appropriate temperature in isGoodChange or another function */
            /* Upgrade temperature */
            if (itime % anneal.Tlen == 0)
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
                    anneal.temp = anneal.temp * anneal.Tcool;

                displayProgress3("time %ld temp %f Complete %ld%% currval %.4f\n",
                    itime, anneal.temp, (int)itime * 100 / anneal.iterations, reserve.total);

            } /* reduce temperature */

            if (fnames.savesnapsteps && !(itime % fnames.savesnapfrequency))
            {
                tempname2 = savename + "_snap" + paddedRun + utils::intToPaddedString(++snapcount, 5) + getFileNameSuffix(fnames.savesnapchanges);
                writeSolution(puno, R, pu, tempname2, fnames.savesnapsteps, fnames);
            } /* Save snapshot every savesnapfreq timesteps */

            iPreviousR = R[ipu];
            iGoodChange = isGoodChange(change, anneal.temp, float_range, rngEngine);

            if (iGoodChange)
            {
                ++ichanges;

                doChange(ipu, puno, R, reserve, change, pu, SM, SM_out, spec, connections, itemp, clumptype, logBuffer);

                if (fnames.savesnapchanges && !(ichanges % fnames.savesnapfrequency))
                {
                    tempname2 = savename + "_snap" + paddedRun +  utils::intToPaddedString(++snapcount, 5) + getFileNameSuffix(fnames.savesnapchanges);
                    writeSolution(puno, R, pu, tempname2, fnames.savesnapchanges, fnames);
                } /* Save snapshot every savesnapfreq changes */

            } /* Good change has been made */

            if (anneal.type == 3)
            {
                anneal.sum += reserve.total;
                anneal.sum2 += reserve.total * reserve.total;
            } /* Keep track of scores for averaging stuff */

            if (verbosity > 4)
                fprintf(fp, "%li,%li,%i,%li,%li,%i,%li,%f,%f,%f,%f,%f\n",
                    itime, ipu, pu[ipu].id, iPreviousR, itemp, R[ipu], iGoodChange, change.total, change.cost, change.connection, change.penalty, anneal.temp);

            if (fnames.saveannealingtrace)
            {
                iRowCounter++;
                if (iRowCounter > iRowLimit)
                    iRowCounter = 1;

                if (iRowCounter == 1)
                {
                    fprintf(Rfp, "%li", itime);

                    fprintf(ttfp, "%li,%f,%li,%f,%i,%f,%f,%f,%f",
                        itime, costthresh, iGoodChange, reserve.total,
                        reserve.pus, reserve.cost, reserve.connection, reserve.penalty, reserve.shortfall);
                    if (fProb1D == 1)
                        fprintf(ttfp, ",%f", reserve.probability1D);
                    if (fProb2D == 1)
                        fprintf(ttfp, ",%f", reserve.probability2D);
                    fprintf(ttfp, ",%li\n", ipu);
                    // iteration,threshold,dochange,total,pus,cost,connectivity,penalty,probability
                    for (i = 0; i < puno; i++)
                        fprintf(Rfp, ",%i", R[i]);

                    fprintf(Rfp, "\n");
                }
            }

        } /* Run Through Annealing */

        /** Post Processing  **********/
        if (aggexist)
            ClearClumps(spno, spec, pu, SM, SM_out);

        if (verbosity > 4)
            fclose(fp);

        if (fnames.saveannealingtrace)
        {
            fclose(ttfp);
            fclose(Rfp);
        }
    } // thermalAnnealing


} // namespace marxan
