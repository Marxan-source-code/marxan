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

#include "marxan.hpp"

namespace marxan {
    using namespace algorithms;
    using namespace utils;
    void SpeciesAmounts(int spno, int puno, vector<sspecies>& spec, const vector<spustuff>& pu, vector<spu>& SM,
        vector<int>& R, int clumptype)
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
    void setBlockDefinitions(int gspno, int spno, int puno, vector<sgenspec>& gspec, vector<sspecies>& spec, vector<spustuff>& PU, vector<spu>& SM)
    {
        int igsp, isp, ipu;
        double totalamount;

        for (igsp = 0; igsp < gspno; igsp++)
        {
            if (gspec[igsp].prop > 0) // deal with percentage in a different way
            {
                for (isp = 0; isp < spno; isp++)
                {
                    if (spec[isp].type == gspec[igsp].type && spec[isp].target < 0)
                    {
                        spec[isp].target = computeTotalSpecAmtAllPu(PU, SM, isp) * gspec[igsp].prop;
                    } // Setting target with percentage
                }
            }
            if (gspec[igsp].target > 0)
            {
                for (isp = 0; isp < spno; isp++)
                {
                    if (spec[isp].type == gspec[igsp].type && spec[isp].target < 0)
                        spec[isp].target = gspec[igsp].target;
                }
            }
            if (gspec[igsp].target2 > 0)
            {
                for (isp = 0; isp < spno; isp++)
                {
                    if (spec[isp].type == gspec[igsp].type && spec[isp].target2 < 0)
                    {
                        spec[isp].target2 = gspec[igsp].target2;
                    }
                }
            }
            if (gspec[igsp].targetocc > 0)
            {
                for (isp = 0; isp < spno; isp++)
                {
                    if (spec[isp].type == gspec[igsp].type && spec[isp].targetocc < 0)
                        spec[isp].targetocc = gspec[igsp].targetocc;
                }
            }
            if (gspec[igsp].sepnum > 0)
            {
                for (isp = 0; isp < spno; isp++)
                {
                    if (spec[isp].type == gspec[igsp].type && spec[isp].sepnum < 0)
                        spec[isp].sepnum = gspec[igsp].sepnum;
                }
            }
            if (gspec[igsp].spf > 0)
            {
                for (isp = 0; isp < spno; isp++)
                {
                    if (spec[isp].type == gspec[igsp].type && spec[isp].spf < 0)
                        spec[isp].spf = gspec[igsp].spf;
                }
            }
            if (gspec[igsp].sepdistance > 0)
            {
                for (isp = 0; isp < spno; isp++)
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


    // compute initial penalties for species with a greedy algorithm.
    // If species has spatial requirements then CalcPenaltyType4 is used instead
    int computePenalties(int puno, int spno, const vector<spustuff>& pu, vector<sspecies>& spec,
        const vector<sconnections>& connections, const vector<spu>& SM, vector<spu_out>& SM_out, vector<int>& PUtemp, int aggexist, double cm, int clumptype, rng_engine& rngEngine)
    {
        int i, j, ibest, imaxtarget, itargetocc;
        double ftarget, fbest, fbestrat, fcost, ftemp, rAmount, rAmountBest;
        int badspecies = 0, goodspecies = 0;
        initialiseReserve(0, pu, PUtemp, rngEngine); // Initialize reserve to 0 and fixed. 

        for (i = 0; i < spno; i++)
        {
            if (spec[i].target2 || spec[i].sepnum)
            {
                j = CalcPenaltyType4(i, puno, SM, SM_out, connections, spec, pu, cm, clumptype, rngEngine);
                badspecies += (j > 0);
                goodspecies += (j < 0);

                appendTraceFile("CalcPenalties spname %i penalty %g\n", spec[i].name, spec[i].penalty);

                continue;
            } // Species has aggregation requirements

            computeFixedPenaltyForSpec(PUtemp, pu, SM, connections, i, ftarget, itargetocc, spec[i].penalty, cm, asymmetricconnectivity);

            // Already adequately represented on type 2 planning unit
            if (ftarget >= spec[i].target && itargetocc >= spec[i].targetocc)
            {
                goodspecies++;
                displayProgress2("Species %i (%s) has already met target %.2f\n",
                    spec[i].name, spec[i].sname.c_str(), spec[i].target);

                appendTraceFile("CalcPenalties spname %i penalty %g\n", spec[i].name, spec[i].penalty);

                continue;
            } // Target met in unremovable reserve

            // Reset non fixed pu
            for (int j = 0; j < puno; j++)
                if (PUtemp[j] < 2) PUtemp[j] = 0;

            do
            {
                fbest = 0; imaxtarget = 0; fbestrat = 0, rAmountBest = 0;
                for (j = 0; j < puno; j++)
                { // trying to find best pu
                    if (PUtemp[j] == 0)
                    {
                        rAmount = returnAmountSpecAtPu(pu[j], SM, i).second;
                        if (rAmount > 0) {
                            fcost = computePlanningUnitValue(pu[j], connections[j], cm, asymmetricconnectivity);
                            if (fcost == 0)
                                fcost = delta;
                            if (rAmount >= spec[i].target - ftarget && (imaxtarget == 0 || (imaxtarget == 1 && fcost < fbest)))
                            { // can I meet the target cheaply?
                                imaxtarget = 1;
                                ibest = j;
                                fbest = fcost;
                                rAmountBest = rAmount;
                            }
                            else {
                                if (fbestrat < rAmount / fcost)
                                { // finding the cheapest planning unit
                                    fbest = fcost;
                                    fbestrat = rAmount / fbest;
                                    ibest = j;
                                    rAmountBest = rAmount;
                                }
                            }

#ifdef DEBUGCALCPENALTIES
                            appendTraceFile("CalcPenalties species %i puid %i cost %g\n", spec[i].name, pu[j].id, fcost);
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
                    appendTraceFile("CalcPenalties species %i puid %i ftarget %g fbest %g\n", spec[i].name, pu[ibest].id, ftarget, fbest);
#endif
                } // Add pu to target
            } while ((ftarget < spec[i].target || itargetocc < spec[i].targetocc) && fbest > 0); // or no more pu left

            if (fbest == 0) // Could not meet target using all available PUs
            { // If not met target with all available PUs
                displayProgress2("Species %d (%s) cannot reach target %.2f there is only %.2f available.\n",
                    spec[i].name, spec[i].sname.c_str(), spec[i].target, ftarget);
                if (ftarget == 0)
                    ftarget = delta;  // Protect against divide by zero
                ftemp = 0;
                if (ftarget < spec[i].target)
                    ftemp = spec[i].target / ftarget;
                if (itargetocc < spec[i].targetocc && itargetocc)  // If ! itargetocc then also !ftarget
                    ftemp += (float)spec[i].targetocc / (float)itargetocc;
                spec[i].penalty = spec[i].penalty * ftemp; // Scale it up
                // This value will be ~ 1/delta when there are no occ's of target species in system
                badspecies++;
            }

#ifdef DEBUGTRACEFILE
            appendTraceFile("CalcPenalties spname %i penalty %g target %g\n", spec[i].name, spec[i].penalty, spec[i].target);
#endif
        }  // Penalty for each individual Species
        // Clear clumps in case I needed them for target4 species

        if (aggexist)
            ClearClumps(spno, spec, pu, SM, SM_out);

        if (goodspecies)
            displayProgress1("%i species are already adequately represented.\n", goodspecies);

        return(badspecies);
    }

    // compute initial penalties for species with a greedy algorithm.
    int computePenaltiesOptimise(int puno, int spno, vector<spustuff>& pu, vector<sspecies>& spec,
        vector<sconnections>& connections, vector<spu>& SM, vector<spu_out>& SM_out, vector<spusporder>& SMsp,
        vector<int>& PUtemp, int aggexist, double cm, int clumptype, rng_engine& rngEngine)
    {
        int i, j, ibest, imaxtarget, itargetocc, ism, ipu;
        double ftarget, fbest, fbestrat, fcost, ftemp, rAmount, r_ibest_amount;
        int badspecies = 0, goodspecies = 0;

        appendTraceFile("CalcPenaltiesOptimise start\n");

        initialiseReserve(puno, pu, PUtemp, rngEngine); // Adds existing reserve to PUtemp

        for (i = 0; i < spno; i++)
        {
            appendTraceFile("CalcPenaltiesOptimise spname %i\n", spec[i].name);

            if (spec[i].target2 || spec[i].sepnum)
            {
                j = CalcPenaltyType4(i, puno, SM, SM_out, connections, spec, pu, cm, clumptype, rngEngine);
                badspecies += (j > 0);
                goodspecies += (j < 0);
                continue;
            } // Species has aggregation requirements

            ftarget = 0;
            itargetocc = 0;
            spec[i].penalty = 0;

            if (spec[i].richness > 0)
            {
                for (j = 0; j < spec[i].richness; j++)  // traverse pu's containing this sp
                { // reset PUtemp and also target
                    ism = spec[i].offset + j;
                    ipu = SMsp[ism].puindex;

                    if (PUtemp[ipu] < 2)
                        PUtemp[ipu] = 0;
                    if (PUtemp[ipu] == 2)
                    {
                        ftarget += SMsp[ism].amount;
                        itargetocc++;
                        spec[i].penalty += computePlanningUnitValue(pu[ipu], connections[ipu], cm, asymmetricconnectivity);
                    }
                }
            }

            // Already adequately represented on type 2 planning unit
            if (ftarget >= spec[i].target && itargetocc >= spec[i].targetocc)
            { // Target met in unremovable reserve
                goodspecies++;
                displayProgress2("Species %i (%s) has already met target %.2f\n",
                    spec[i].name, spec[i].sname.c_str(), spec[i].target);
                continue;
            }

            do
            {
                fbest = 0; imaxtarget = 0; fbestrat = 0;
                if (spec[i].richness > 0)
                {
                    for (j = 0; j < spec[i].richness; j++)  // traverse pu's containing this sp
                    { // trying to find best pu
                        ism = spec[i].offset + j;
                        ipu = SMsp[ism].puindex;

                        rAmount = SMsp[ism].amount;
                        if (PUtemp[ipu] == 0)
                        { // Making sure only checking planning units not already used
                            fcost = computePlanningUnitValue(pu[ipu], connections[ipu], cm, asymmetricconnectivity);
                            if (fcost == 0)
                                fcost = delta;
                            if (rAmount >= spec[i].target - ftarget && (imaxtarget == 0
                                || (imaxtarget == 1 && fcost < fbest)))
                            { // can I meet the target cheaply?
                                imaxtarget = 1;
                                ibest = ipu;
                                r_ibest_amount = rAmount;
                                fbest = fcost;
                            }
                            else {
                                if (fbestrat < rAmount / fcost)
                                { // finding the cheapest planning unit
                                    fbest = fcost;
                                    fbestrat = rAmount / fbest;
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
            } while ((fbest > 0) && (ftarget < spec[i].target || itargetocc < spec[i].targetocc));
            // while there is some pu's with this species to test AND a best available pu was found AND targets are not met yet

            if (fbest == 0) // Could not meet target using all available PUs
            { // If not met target with all available PUs
                displayProgress2("Species %d (%s) cannot reach target %.2f there is only %.2f available.\n",
                    spec[i].name, spec[i].sname.c_str(), spec[i].target, ftarget);
                if (ftarget == 0)
                    ftarget = delta;  // Protect against divide by zero
                ftemp = 0;
                if (ftarget < spec[i].target)
                    ftemp = spec[i].target / ftarget;
                if (itargetocc < spec[i].targetocc && itargetocc)  // If ! itargetocc then also !ftarget
                    ftemp += (float)spec[i].targetocc / (float)itargetocc;
                spec[i].penalty = spec[i].penalty * ftemp; // Scale it up
                /* This value will be ~ 1/delta when there are no occ's of target species in system*/
                badspecies++;
            }

            appendTraceFile("CalcPenaltiesOptimise spname %i penalty %g\n", spec[i].name, spec[i].penalty);
        }  // Penalty for each individual Species

        // Clear clumps in case I needed them for target4 species
        if (aggexist)
            ClearClumps(spno, spec, pu, SM, SM_out);

        if (goodspecies)
            displayProgress1("%i species are already adequately represented.\n", goodspecies);

        appendTraceFile("CalcPenaltiesOptimise end\n");

        return(badspecies);
    } // computePenaltiesOptimise

    // compute change in the species representation for adding or removing a single planning unit or set of planning units
    double computeChangePenalty(int ipu, int puno, vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out,
        const vector<int>& R, const vector<sconnections>& connections, int imode, int clumptype, double& rShortfall)
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
                if (SM[ism].amount != 0.0)  /** Only worry about PUs where species occurs and target != 0 **/
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
                        newamount = NewPenalty4(ipu, isp, puno, spec, pu, SM, SM_out, R, connections, imode, clumptype);
                    }
                    else
                    {
                        if (spec[isp].target)
                            newamount = computeSpeciesPlanningUnitPenalty(isp, spec, SM, ism, imode) / spec[isp].target;
                        if (spec[isp].targetocc)
                        {
                            tamount = (double)(spec[isp].targetocc - spec[isp].occurrence - imode) /
                                (double)spec[isp].targetocc;
                            newamount += tamount < 0 ? 0 : tamount;
                        }
                        if (spec[isp].target && spec[isp].targetocc)
                            newamount /= 2;
                        if (spec[isp].sepnum)
                            newamount += computeSepPenalty(CountSeparation2(isp, ipu, tempSclumps, puno, R, pu, SM, SM_out, spec, imode),
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
    void computeReserveValue(int puno, int spno, const vector<int>& R, const vector<spustuff>& pu,
        const vector<sconnections>& connections, const vector<spu>& SM, vector<spu_out>& SM_out,
        double cm, vector<sspecies>& spec, int aggexist, scost& reserve, int clumptype, stringstream& logBuffer)
    {
        vector<sclumps> tempSclumps;
        int i, j;
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
            SetSpeciesClumps(puno, R, spec, pu, SM, SM_out, connections, clumptype);

        // traverse species, computing penalty for each species
        for (i = 0; i < spno; i++)
        {
            fractionAmount = 0;
            if (spec[i].target > spec[i].amount)
            {
                if (spec[i].target != 0) {
                    fractionAmount = (spec[i].target - spec[i].amount) / spec[i].target;
                }

                reserve.shortfall += spec[i].target - spec[i].amount;
            }

            // does this species have an occurrence target?
            if (spec[i].targetocc > 0)
            {
                if (spec[i].targetocc > spec[i].occurrence)
                {
                    fractionAmount += (double)(spec[i].targetocc - spec[i].occurrence) / (double)spec[i].targetocc;
                    reserve.shortfall += spec[i].targetocc - spec[i].occurrence;
                }
                if (spec[i].target && spec[i].targetocc)
                    fractionAmount /= 2;
            }

            reserve.penalty += fractionAmount * spec[i].penalty * spec[i].spf;

            if (spec[i].sepnum)
            {
                spec[i].separation = CountSeparation2(i, 0, tempSclumps, puno, R, pu, SM, SM_out, spec, 0);
                reserve.penalty += computeSepPenalty(spec[i].separation, spec[i].sepnum) *
                    spec[i].spf * spec[i].penalty;
            }
        }

        // traverse planning units, computing planning unit metrics for the reserve
        for (j = 0; j < puno; j++)
        {
            // if planning unit is protected
            if (R[j] == 1 || R[j] == 2)
            {
                reserve.cost += pu[j].cost;
                reserve.pus += 1;
                rConnectivityValue = ConnectionCost2(connections[j], R, 1, 0, cm, asymmetricconnectivity, fOptimiseConnectivityIn);
                reserve.connection += rConnectivityValue;

#ifdef DEBUG_RESERVECOST
                logBuffer << "puid " << pu[j].id << " connectivity " << rConnectivityValue << endl;
#endif

                if (fProb1D == 1)
                    ReturnProbabilityAmounts1D(ExpectedAmount1D, VarianceInExpectedAmount1D, j, puno, pu, SM);
                if (fProb2D == 1)
                    ReturnProbabilityAmounts2D(ExpectedAmount2D, VarianceInExpectedAmount2D, j, puno, pu, SM);
            }
        }

        if (fProb1D == 1)
        {
            reserve.probability1D = ComputeProbability1D(ExpectedAmount1D, VarianceInExpectedAmount1D, spno, spec);
        }

        if (fProb2D == 1)
        {
            reserve.probability2D = ComputeProbability2D(ExpectedAmount2D, VarianceInExpectedAmount2D, spno, spec);
        }

        reserve.total = reserve.cost + reserve.connection + reserve.penalty + reserve.probability1D + reserve.probability2D;

#ifdef DEBUG_PROB1D
        logBuffer << "probability1D " << reserve.probability1D << endl;

        sProbDebugFileName = fnames.outputdir + "output_Prob1DDebug_" + to_string(reserve.cost) + ".csv";
        writeProb1DDebugTable(spno, sProbDebugFileName,
            ExpectedAmount1D, VarianceInExpectedAmount1D, spec);

        sProbDebugFileName = fnames.outputdir + "output_Prob1DDetailDebug_" + to_string(reserve.cost) + ".csv";
        writeProb1DDetailDebugTable(sProbDebugFileName, puno, spno, pu, SM, R);
#endif

#ifdef DEBUG_PROB2D
        logBuffer << "probability2D " << reserve.probability2D << endl;

        sProbDebugFileName = fnames.outputdir + "output_Prob2DDebug_" + to_string(reserve.cost) + ".csv";
        writeProb2DDebugTable(spno, sProbDebugFileName,
            ExpectedAmount2D, VarianceInExpectedAmount2D, spec);

        sProbDebugFileName = fnames.outputdir + "output_Prob2DDetailDebug_" + to_string(reserve.cost) + ".csv";
        writeProb2DDetailDebugTable(sProbDebugFileName, puno, pu, SM, R);
#endif

        // destroy arrays for prob 1D and 2D
        if (fProb1D == 1)
        {
            for (i = 0; i < spno; i++)
            {
                spec[i].expected1D = ExpectedAmount1D[i];
                spec[i].variance1D = VarianceInExpectedAmount1D[i];
            }
        }
        if (fProb2D == 1)
        {
            for (i = 0; i < spno; i++)
            {
                spec[i].expected2D = ExpectedAmount2D[i];
                spec[i].variance2D = VarianceInExpectedAmount2D[i];
            }
        }
    } // computeReserveValue

    // sets cost threshold penalty when "cost threshold" is in use
    double thresholdPenalty(double tpf1, double tpf2, double timeprop)
    {
        if (tpf2 < 0)
            return(tpf1);
        return(tpf1 * exp(tpf2 * timeprop));
    }

    void computeChangeScore(int iIteration, int ipu, int spno, int puno, const vector<spustuff>& pu, const vector<sconnections>& connections,
        vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, const vector<int>& R, double cm, int imode,
        scost& change, scost& reserve, double costthresh, double tpf1, double tpf2,
        double timeprop, int clumptype)
        // imode = 1 add PU, imode = -1 remove PU
    {
        double threshpen = 0;
        int threshtype = 1; /*Debugging line. This should be input parameter not hardwired */
        double tchangeconnection, tresconnection;
#ifdef DEBUGCHECKCHANGE
        char debugline[200];
#endif

        change.cost = pu[ipu].cost * imode; /* Cost of this PU on it's own */
        change.connection = ConnectionCost2(connections[ipu], R, imode, 1, cm, asymmetricconnectivity, fOptimiseConnectivityIn);


        change.penalty = computeChangePenalty(ipu, puno, spec, pu, SM, SM_out, R, connections, imode, clumptype, change.shortfall);

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


    // change the status of a single planning unit
    void doChange(int ipu, int puno, vector<int>& R, scost& reserve, scost& change,
        const vector<spustuff>& pu, const vector<spu>& SM, vector<spu_out>& SM_out, vector<sspecies>& spec, const vector<sconnections>& connections,
        int imode, int clumptype, stringstream& logBuffer)
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
                        AddNewPU(ipu, isp, connections, spec, pu, SM, SM_out, clumptype);
                    }
                    else
                    {
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
                        spec[isp].separation = CountSeparation2(isp, 0, tempSclumps, puno, R, pu, SM, SM_out, spec, 0);
                    }
            }
        }

        reserve.total = reserve.cost + reserve.connection + reserve.penalty + reserve.probability1D + reserve.probability2D;
    } // doChange


} // namespace marxan

