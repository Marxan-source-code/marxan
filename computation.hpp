#pragma once

// All functions here should be unit and should not have dependencies on the core marxan files.

#include <algorithm>

#include "marxan.hpp"
#include "utils.hpp"

namespace marxan {

    // Sep Penalty
    // This returns the penalty for not meeting separation requirments. Feed in sepnum and current
    //  separation and returns a value from 0 to 1 which is an artificial shortfall.
    inline double computeSepPenalty(int ival, int itarget)
    {
        double fval;

        if (!itarget)
            return (0); /* no penalty if no separation requirement*/

        fval = (double)ival / (double)itarget;
        if (!ival)
            fval = 1.0 / (double)itarget;

        return (1 / (7 * fval + 0.2) - (1 / 7.2)); // Gives a nice hyperbole with fval = 1 return 0 and fval = 0 or 0.1 returning almost 1
    } // SepPenalty

    // Supplementary function for computeRepresentationMISSLEVEL that updates shortfall and rMinimumProportionMet based on given target and amount
    inline void checkAndUpdateTargetProportion(double target, double amount, double& shortfall, double& rMinimumProportionMet) {
        if (target > 0)
            if (amount < target)
            {
                shortfall += target - amount;
                double rProportionMet = amount / target;

                if (rProportionMet < rMinimumProportionMet)
                    rMinimumProportionMet = rProportionMet;
            }
    }

    // Sums all connectivity edges for a pu.

    double connectionCost1(const sconnections& connections, double cm, int asymmetricconnectivity);

    /***  Counts the number of species missing from the reserve ****/
    // compute the number of species whose representation fraction is less that the MISSLEVEL parameter
    int computeRepresentationMISSLEVEL(int spno, const vector<sspecies>& spec, double misslevel, double& shortfall, double& rMinimumProportionMet);

    // compute connectivity total, in, edge, out for summary report
    void computeConnectivityIndices(double& rConnectivityTotal, double& rConnectivityIn,
        double& rConnectivityEdge, double& rConnectivityOut,
        int puno, const vector<int>& R, const vector<sconnections>& connections);
        // We record 4 categories for connectivity;
        //  - total, all connections in the region
        //  - in, all connections entirely within the reserve network (ie. both pu's in)
        //  - edge, all connections on the edge of the reserve network (ie. one pu in & one pu out)
        //  - out, all connections not captured in the reserve network (ie. both pu's out)
        //
        // Of these, we previously only recorded "edge", referring to it as boundary length.
        // The proportion of connections captured is given by;
        //  in / total
        //
        // total = in + edge + out


    // Merging repeated area computation code. Recomputes TO, TA etc. from scratch.
    void computeOccurrencesAndAreas(int& puno, const vector<spustuff>& pu, const vector<spu>& SM,
        vector<int>& TotalOccurrences, vector<int>& TO_2, vector<int>& TO_3,
        vector<double>& TotalAreas, vector<double>& TA_2, vector<double>& TA_3);

    // used for writeProb1DDebugTable and writeProb2DDebugTable and other functions needing pr computation
    // used for writeProb1DDebugTable and writeProb2DDebugTable and other functions needing pr computation
    inline void computeProbMeasures(double variance, double target, double prTarget, double expected,
        double& rZ, double& rRawP, int& iHeavisideStepFunction, double& rShortfallPenalty, double& rP) {
        if (variance > 0)
            rZ = (target - expected) / sqrt(variance);
        else
            rZ = 4;

        if (rZ >= 0)
            rRawP = utils::probZUT(rZ);
        else
            rRawP = 1 - utils::probZUT(-1 * rZ);

        if (prTarget > rRawP)
            iHeavisideStepFunction = 1;
        else
            iHeavisideStepFunction = 0;

        if (prTarget > 0)
            rShortfallPenalty = (prTarget - rRawP) / prTarget;
        else
            rShortfallPenalty = 0;

        rP = iHeavisideStepFunction * rShortfallPenalty;
    }

    // Used in ReturnProbabilityAmounts1D and ReturnProbabilityAmounts2D
    void computeExpectedAndVariance(int ipu, const vector<spustuff>& pu, const vector<spu>& SM, vector<double>& variance, vector<double>& expected);

    // compute cost + connectivity for a single planning unit
    inline double computePlanningUnitValue(const spustuff& pu, const sconnections& connections, double cm, int asymmetricconnectivity)
    {
        double theValue;

        theValue = pu.cost;
        theValue += connectionCost1(connections, cm, asymmetricconnectivity);

        return(theValue);
    }

    // Returns both index and species amount for a planning unit
    pair<int, double> returnAmountSpecAtPu(const spustuff& pu, const vector<spu>& SM, int iSpecIndex);

    // compute penalty for a species for changing status of a single planning unit
    inline double computeSpeciesPlanningUnitPenalty(int ipu, int isp, const vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, int imode)
    {
        double newpen = max(0.0, spec[isp].target - spec[isp].amount - returnAmountSpecAtPu(pu[ipu], SM, isp).second * imode);
        return(newpen);
    }

    // compute penalty for a species for changing status of a single planning unit
    inline double computeSpeciesPlanningUnitPenalty(int isp, const vector<sspecies>& spec, const vector<spu>& SM, int ism, int imode)
    {
        double newpen = max(0.0, spec[isp].target - spec[isp].amount - SM[ism].amount * imode);
        return(newpen);
    }

    // Sums the total spec amount across all pu
    double computeTotalSpecAmtAllPu(const vector<spustuff>& PU, const vector<spu>& SM, int speciesInd);

    // compute proportional target for species when prop target is specified
    // use the prop value from the conservation feature file to set a proportion target for species
    void computeSpecProp(int spno, vector<sspecies>& spec, int puno, const vector<spustuff>& pu, const vector<spu>& SM);

    // Computes penalty for a species given a reserve, based on only fixed pu
    void computeFixedPenaltyForSpec(const vector<int>& R, const vector<spustuff>& pu, const vector<spu>& SM, const vector<sconnections>& connections, int spIndex,
        double& ftarget, int& itargetocc, double& penalty, double cm, int asymmetricconnectivity);

    // ********* Connection Cost Type 2 **************
    // **  Requires R[]. imode2 = 0 there is no negative cost for removing connection, we are calling from ReserveCost
    //                         or 1 there is a negative cost for removing connection, we are calling from Annealing
    //                   imode = -1 we are removing the planning unit from a reserve, calling from Annealing
    //                        or 1  we are adding the planning unit to a reserve, or it is already in reserve
    //      It seems that the behaviour of this function is undefined/unsupported if imode2=0 and imode=-1
    double ConnectionCost2(const sconnections& connection, const vector<int>& R, int imode, int imode2, double cm,
        int asymmetricconnectivity, int fOptimiseConnectivityIn);

} // namespace marxan