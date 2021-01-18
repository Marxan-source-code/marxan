#pragma once

// Contains functions relating to computational functions. E.g. computing penalty, shortfall etc.
// Moving to a new file for easier unit testing. 
// All functions here should be unit and should not have dependencies on the core marxan files.

#include <algorithm>

#include "marxan.hpp"
#include "utils.hpp"

namespace marxan {

inline
// Sep Penalty
// This returns the penalty for not meeting separation requirments. Feed in sepnum and current
//  separation and returns a value from 0 to 1 which is an artificial shortfall.
double computeSepPenalty(int ival,int itarget)
{
    double fval;

    if (!itarget)
        return (0); /* no penalty if no separation requirement*/
        
    fval = (double) ival / (double) itarget;
    if (!ival)
        fval = 1.0 /(double) itarget;

    return (1/(7*fval+0.2)-(1/7.2)); // Gives a nice hyperbole with fval = 1 return 0 and fval = 0 or 0.1 returning almost 1
} // SepPenalty

inline
// Supplementary function for computeRepresentationMISSLEVEL that updates shortfall and rMinimumProportionMet based on given target and amount
void checkAndUpdateTargetProportion(double target, double amount, double& shortfall, double& rMinimumProportionMet) {
    if (target > 0)
        if (amount < target)
        {
            shortfall += target - amount;
            double rProportionMet = amount/target;

            if (rProportionMet < rMinimumProportionMet)
                rMinimumProportionMet = rProportionMet;
        }
}

// Sums all connectivity edges for a pu.
inline
double connectionCost1(const sconnections &connections, double cm, int asymmetricconnectivity)
{
    double fcost;

    fcost = connections.fixedcost;
    for (const auto &p : connections.first)
    {
        if (asymmetricconnectivity)
        {
            if (p.connectionorigon)
                fcost += p.cost;
        }
        else
        {
            fcost += p.cost;
        }
    }

    return (fcost * cm);
}

inline
/***  Counts the number of species missing from the reserve ****/
// compute the number of species whose representation fraction is less that the MISSLEVEL parameter
int computeRepresentationMISSLEVEL(int spno,vector<sspecies> &spec,double misslevel,double &shortfall,double &rMinimumProportionMet)
{
    int i,isp = 0;
    double rProportionMet;

    shortfall = 0;
    rMinimumProportionMet = 1;
    for (i=0;i<spno;i++)
    {
        checkAndUpdateTargetProportion(spec[i].target, spec[i].amount, shortfall, rMinimumProportionMet); // check regular target
        checkAndUpdateTargetProportion(spec[i].targetocc, spec[i].occurrence, shortfall, rMinimumProportionMet); // check occurence target

        if (spec[i].target)
        {
            if (spec[i].amount/spec[i].target < misslevel)
            {
                isp++;
                continue;
            }
        }

        if (spec[i].targetocc)
        {
            if ((double)spec[i].occurrence/(double)spec[i].targetocc < misslevel)
            {
                isp++;
                continue;
            }
        }
        if (spec[i].sepdistance && spec[i].separation < 3)
        {
            isp++;  /* count species if not met separation and not already counted */
        }
    }

    return(isp);
}

inline
// compute connectivity total, in, edge, out for summary report
void computeConnectivityIndices(double &rConnectivityTotal, double &rConnectivityIn,
                                double &rConnectivityEdge, double &rConnectivityOut,
                                int puno, vector<int> &R, vector<sconnections> &connections)
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
{
    int i;
    double rFixed;

    for (i=0;i<puno;i++)
    {
        rFixed = connections[i].fixedcost;

        rConnectivityTotal += rFixed;

        if (R[i]==1 || R[i] == 2)
        { // add to 'in' or 'edge'
            rConnectivityEdge += rFixed;

            for(sneighbour p: connections[i].first)
            {
                if (p.nbr > i)
                {
                    if (R[p.nbr] == 1 || R[p.nbr] == 2) // add to 'in'
                        rConnectivityIn += p.cost;
                    else  // add to 'edge'
                        rConnectivityEdge += p.cost;

                    // add to 'total'
                    rConnectivityTotal += p.cost;
                }
            }
        }
        else
        { // add to 'out' or 'edge'
            rConnectivityOut += rFixed;

            for(sneighbour p: connections[i].first)
            {
                if (p.nbr > i)
                {
                    if (R[p.nbr] == 1 || R[p.nbr] == 2) // add to 'edge'
                        rConnectivityEdge += p.cost;
                    else  // add to 'out'
                        rConnectivityOut += p.cost;

                    // add to 'total'
                    rConnectivityTotal += p.cost;
                }
            }
        }
    }
}

inline
// Merging repeated area computation code. Recomputes TO, TA etc. from scratch.
void computeOccurrencesAndAreas (int& puno, vector<spustuff>& pu, vector<spu>& SM, 
    vector<int>& TotalOccurrences, vector<int>& TO_2, vector<int>& TO_3, 
    vector<double>& TotalAreas, vector<double>& TA_2, vector<double>& TA_3) {
    
    int ism, isp;
    for (int ipu=0;ipu<puno;ipu++)
    {
        if (pu[ipu].richness)
        {
            for (int i=0;i<pu[ipu].richness;i++)
            {
                ism = pu[ipu].offset + i;
                isp = SM[ism].spindex;

                TotalOccurrences[isp]++;
                TotalAreas[isp] += SM[ism].amount;

                if (pu[ipu].status == 2)
                {
                    TO_2[isp]++;
                    TA_2[isp] += SM[ism].amount;
                }

                if (pu[ipu].status == 3)
                {
                    TO_3[isp]++;
                    TA_3[isp] += SM[ism].amount;
                }
            }
        }
    }
}

inline 
// used for writeProb1DDebugTable and writeProb2DDebugTable and other functions needing pr computation
void computeProbMeasures(double variance, double target, double prTarget, double expected, 
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

    if (prTarget> 0)
        rShortfallPenalty = (prTarget - rRawP) / prTarget;
    else
        rShortfallPenalty = 0;

    rP = iHeavisideStepFunction * rShortfallPenalty;
}

inline
// Used in ReturnProbabilityAmounts1D and ReturnProbabilityAmounts2D
void computeExpectedAndVariance(int ipu, vector<spustuff> &pu,vector<spu> &SM, vector<double>& variance, vector<double>& expected) {
    int i, ism, isp;

    if (pu[ipu].richness)
    {
        for (i = 0; i < pu[ipu].richness; i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            if (SM[ism].amount)
            {
                expected[isp] += SM[ism].amount * (1 - pu[ipu].prob);
                variance[isp] += SM[ism].amount * SM[ism].amount * pu[ipu].prob * (1 - pu[ipu].prob);
            }
        }
    }
}

// compute cost + connectivity for a single planning unit
inline
double computePlanningUnitValue(const spustuff &pu,const sconnections &connections, double cm, int asymmetricconnectivity)
{
    double theValue;

    theValue = pu.cost;
    theValue += connectionCost1(connections, cm, asymmetricconnectivity);

    return(theValue);
}

// Returns both index and species amount for a planning unit
inline
pair<int,double> returnAmountSpecAtPu(const spustuff &pu, const vector<spu> &SM, int iSpecIndex)
{
    if (pu.richness > 0)
    {
        auto start_it = SM.begin() + pu.offset;
        auto end_it = start_it + pu.richness;
        auto spindex_cmp = [](const spu& lhs, int rhs) -> bool { return lhs.spindex < rhs; };
        auto elem_it = std::lower_bound(start_it, end_it, iSpecIndex, spindex_cmp);
        if (elem_it != end_it)
        {
            size_t index = elem_it - SM.begin();
            return pair<int, double>(index, elem_it->amount);
        }
    }
    return pair<int, double>(-1, 0);
} 

// compute penalty for a species for changing status of a single planning unit
inline
double computeSpeciesPlanningUnitPenalty(int ipu,int isp, const vector<sspecies> &spec, const vector<spustuff> &pu, vector<spu> &SM, int imode)
{
    double newpen = max(0.0, spec[isp].target - spec[isp].amount - returnAmountSpecAtPu(pu[ipu],SM,isp).second*imode);
    return(newpen);
}

// compute penalty for a species for changing status of a single planning unit
inline
double computeSpeciesPlanningUnitPenalty( int isp, const vector<sspecies>& spec, const vector<spu>& SM, int ism, int imode)
{
    double newpen = max(0.0, spec[isp].target - spec[isp].amount - SM[ism].amount * imode);
    return(newpen);
}

// Sums the total spec amount across all pu
inline
double computeTotalSpecAmtAllPu(vector<spustuff> &PU, vector<spu> &SM, int speciesInd) 
{   
    double totalAmount = 0.0;
    for (int ipu=0;ipu<PU.size();ipu++)
        totalAmount += returnAmountSpecAtPu(PU[ipu],SM,speciesInd).second;
    return totalAmount;
}

// compute proportional target for species when prop target is specified
// use the prop value from the conservation feature file to set a proportion target for species
inline
void computeSpecProp(int spno, vector<sspecies> &spec, int puno, vector<spustuff> &pu, vector<spu> &SM)
{
    // compute and set target for species with a prop value
    for (int isp=0;isp<spno;isp++)
    {
        if (spec[isp].prop > 0)
        {
            spec[isp].target = computeTotalSpecAmtAllPu(pu,SM,isp) * spec[isp].prop;
        }
    }
}

// Computes penalty for a species given a reserve, based on only fixed pu
inline
void computeFixedPenaltyForSpec(const vector<int> &R, const vector<spustuff> &pu, const vector<spu> &SM, const vector<sconnections> &connections, int spIndex, 
    double& ftarget, int& itargetocc, double& penalty, double cm, int asymmetricconnectivity) 
{
    ftarget = 0, itargetocc = 0, penalty = 0;

    for (int j=0;j<pu.size();j++)
    {
        if (R[j] == 2)
        {
            ftarget += returnAmountSpecAtPu(pu[j],SM,spIndex).second;
            itargetocc++;
            penalty += computePlanningUnitValue(pu[j],connections[j],cm, asymmetricconnectivity);
        }
    }
}

inline
// ********* Connection Cost Type 2 **************
// **  Requires R[]. imode2 = 0 there is no negative cost for removing connection, we are calling from ReserveCost
//                         or 1 there is a negative cost for removing connection, we are calling from Annealing
//                   imode = -1 we are removing the planning unit from a reserve, calling from Annealing
//                        or 1  we are adding the planning unit to a reserve, or it is already in reserve
//      It seems that the behaviour of this function is undefined/unsupported if imode2=0 and imode=-1
double ConnectionCost2(int ipu, vector<sconnections> &connections, vector<int> &R, int imode, int imode2, double cm, 
    int asymmetricconnectivity, int fOptimiseConnectivityIn)
{
    double fcost  = connections[ipu].fixedcost * imode;
    int R_pu1;

    if (asymmetricconnectivity)
    {
        for (sneighbour &p : connections[ipu].first)
        {
            if (imode2) // calling from Annealing
            {
                // determines if ipu is currently switched on or not
                // if imode==1 then we assume currently switched off, and will switch on.
                if (imode == 1)
                    R_pu1 = 0;
                else
                    R_pu1 = 1;

                if (p.connectionorigon)
                {
                    if (R[p.nbr] == 0)
                    {
                        if (R_pu1 == 1)
                        {
                            fcost -= p.cost;
                        }
                        else
                        {
                            fcost += p.cost;
                        }
                    }
                }
                else
                {
                    if (R[p.nbr] == 1 || R[p.nbr] == 2)
                    {
                        if (R_pu1 == 1)
                        {
                            fcost += p.cost;
                        }
                        else
                        {
                            fcost -= p.cost;
                        }
                    }
                }
            }
            else // calling from ReserveCost
            {
                if (R[p.nbr] == 0)
                    if (p.connectionorigon)
                    {
                        fcost += p.cost;
                    }
            }
        }
    }
    else
    {
        for (sneighbour &p : connections[ipu].first) // treatment for symmetric connectivity
        {
            if (fOptimiseConnectivityIn == 1)
            { // optimise for "Connectivity In"
                if (R[p.nbr] == 1 || R[p.nbr] == 2)
                {
                    fcost += imode * p.cost;
                }
                else
                {
                    fcost += imode * imode2 * p.cost * -1;
                }
            }
            else
            { // optimise for "Connectivity Edge"
                if (R[p.nbr] == 1 || R[p.nbr] == 2)
                {
                    fcost += imode * imode2 * p.cost * -1;
                }
                else
                {
                    fcost += imode * p.cost;
                }
            }
        }
    }

    return (fcost * cm);
}

} // namespace marxan