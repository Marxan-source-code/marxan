
#include "marxan.hpp"
#include "utils.hpp"
#include "computation.hpp"

// functions that read from input files
namespace marxan {

double GreedyPen(int ipu, int puno, int spno, vector<sspecies> &spec,vector<int> &R,vector<spustuff> &pu,
                 vector<spu> &SM,int clumptype, int thread)
{
    double famount = 0.0, fold, newamount;

    for (int i = 0; i < spno; i++)
    {
        fold = (spec[i].target - spec[i].amount);
        if (fold > 0)
        {
            if (spec[i].target2)
                newamount = NewPenalty4(ipu, i, puno, spec, pu, SM, R, connections, 1, clumptype, thread);
            else
                newamount = computeSpeciesPlanningUnitPenalty(ipu, i, spec, pu, SM, 1);

            famount += (newamount - fold) * spec[i].spf;
        } // Add new penalty if species isn't already in the system
    }

    return (famount); // Negative means decrease in amount missing
} // Greedy Species Penalty

double GreedyScore(int ipu,int puno,int spno, vector<sspecies>& spec,vector<spu> &SM,vector<sconnections> &connections,
                   vector<int> &R,vector<spustuff> &pu,double cm,int clumptype, int thread)
{
    double currpen, currcost, currscore;

    currpen = GreedyPen(ipu, puno, spno, spec, R, pu, SM, clumptype, thread);
    currcost = pu[ipu].cost + ConnectionCost2(connections[ipu], R, 1, 1, cm, asymmetricconnectivity, fOptimiseConnectivityIn);
    if (currcost <= 0)
    {
        currscore = -1.0 / delta;
    } // otherwise this 'free pu' will have +score
    else
    {
        currscore = currpen / currcost;
        // multiply by rand (1.000,1.001)
    }
    return (currscore);
} // Score for a planning unit based upon greedy algorithm

/*********** Rarity Settup. Sets up rare score for each species ******/
/**** score is total abundance / smallest species abundance *********/
void SetRareness(int puno, int spno, vector<double> &Rare,vector<int> &R,vector<spustuff> &pu,vector<spu> &SM, stringstream& logBuffer) 
{
    double smallest = 0;
    vector<double> fcount(spno, 0);
    int i, ism, isp, ipu;

#ifdef DEBUG_HEURISTICS
    logBuffer << "SetRareness start\n";
#endif

    for (ipu = 0; ipu < puno; ipu++)
        if (pu[ipu].richness)
            for (i = 0; i < pu[ipu].richness; i++)
            {
                ism = pu[ipu].offset + i;
                isp = SM[ism].spindex;
                if (R[ipu] < 2)
                    fcount[isp] += SM[ism].amount;
            }

    for (isp = 0; isp < spno; isp++)
    {
        if (smallest == 0 || (fcount[isp] < smallest && fcount[isp] > 0))
            smallest = fcount[isp];
        Rare[isp] = fcount[isp];

#ifdef DEBUG_HEURISTICS
        logBuffer << "SetRareness isp " << isp << " Rare " << Rare[isp] << endl;
#endif
    }

    if (smallest == 0)
        displayErrorMessage("Serious Error in calculating Rarenesses. No species detected.\n");

    for (isp = 0; isp < spno; isp++)
        Rare[isp] /= smallest;

#ifdef DEBUG_HEURISTICS
    logBuffer << "SetRareness end\n";
#endif
}  // SetRareness

// RareScore The score for a particular conservation value on a particular PU
double RareScore(int isp,int ipu,int puno,vector<sspecies> &spec,vector<spu> &SM, vector<int> &R,
    vector<sconnections> &connections, vector<spustuff> &pu, double cm,int clumptype, int thread)
{
    double currpen, currcost, currscore;
    double fold, newamount;
    fold = (spec[isp].target - spec[isp].amount);
    if (fold > 0)
    {
        if (spec[isp].target2)
            newamount = NewPenalty4(ipu, isp, puno, spec, pu, SM, R, connections, 1, clumptype, thread);
        else
            newamount = computeSpeciesPlanningUnitPenalty(ipu, isp, spec, pu, SM, 1);
        currpen = newamount - fold;
    } // Add new penalty if species isn't already in the system

    currcost = pu[ipu].cost + ConnectionCost2(connections[ipu], R, 1, 1, cm, asymmetricconnectivity, fOptimiseConnectivityIn);
    if (currcost <= 0)
    {
        currscore = -1.0 / delta;
    } // otherwise this 'free pu' will have +score
    else
    {
        currscore = currpen / currcost;
        // multiply by rand (1.000,1.001)
    }

    return (currscore);
} // RareScore

// Max Rare Score Heuristic. PU scores based on rarest beast on PU
double MaxRareScore(int ipu,int puno,vector<sspecies> &spec,vector<spu> &SM,
    vector<int> &R,vector<sconnections> &connections,vector<spustuff> &pu,double cm,vector<double> &Rare,int clumptype, int thread)
{
    int ism, isp, rareno = -1;
    double rarest = 0.0, rarescore;

    if (pu[ipu].richness)
        for (int i = 0; i < pu[ipu].richness; i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            if (SM[ism].amount && (spec[isp].target > spec[isp].amount || (spec[isp].sepdistance && spec[isp].separation < 3)))
                if (1.0 / Rare[isp] < rarest || rareno < 0)
                {
                    rareno = isp;
                    rarest = Rare[isp];
                } // Determine which is the rarest species
        }

    if (rareno > -1)
        rarescore = RareScore(rareno, ipu, puno, spec, SM, R, connections, pu, cm, clumptype, thread) / rarest;
    else
        rarescore = 1.0 / delta;

    return (rarescore);
} // Max Rare Score

// Best Rarity Score. Determines each species rare score
double BestRareScore(int ipu,int puno,vector<sspecies> &spec,vector<spu> &SM,
    vector<int> &R,vector<sconnections> &connections,vector<spustuff> &pu,double cm,vector<double> &Rare,int clumptype, int thread)
{
    int i, ism, isp, rareno = -1;
    double rarest = 0, rarescore;

    if (pu[ipu].richness)
        for (i = 0; i < pu[ipu].richness; i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            if (SM[ism].amount && (spec[isp].target > spec[isp].amount || (spec[isp].sepdistance && spec[isp].separation < 3)))
            {
                rarescore = RareScore(isp, ipu, puno, spec, SM, R, connections, pu, cm, clumptype, thread) / Rare[isp];
                if (rarescore > rarest || rareno < 0)
                {
                    rarest = rarescore;
                    rareno = isp;
                }
            }
        }

    return (rarescore);
} // Best Rare Score

// Average Rare Score. Rare Score for each scoring species/number scoring species
double AveRareScore(int ipu,int puno,vector<sspecies> &spec,vector<spu> &SM,
    vector<int> &R,vector<sconnections> &connections,vector<spustuff> &pu,double cm,vector<double> &Rare,int clumptype, int thread)
{
    int i, ism, isp, rareno = 0;
    double rarescore = 0;

    if (pu[ipu].richness)
        for (i = 0; i < pu[ipu].richness; i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            if (SM[ism].amount &&
                (spec[isp].target > spec[isp].amount ||
                 (spec[isp].sepdistance && spec[isp].separation < 3)))
            {
                rarescore += RareScore(isp, ipu, puno, spec, SM, R, connections, pu, cm, clumptype, thread) / Rare[isp];
                rareno++;
            }
        }

    return (rarescore / rareno);
} // Ave Rare Score

double SumRareScore(int ipu,int puno,vector<sspecies> &spec,vector<spu> &SM,
    vector<int> &R,vector<sconnections> &connections,vector<spustuff> &pu,double cm,
    vector<double> &Rare,int clumptype, int thread, stringstream& logBuffer)
{
    int i, ism, isp;
    double rarescore = 0;

#ifdef DEBUG_HEURISTICS
    logBuffer << "SumRareScore start\n";
#endif

    if (pu[ipu].richness)
        for (i = 0; i < pu[ipu].richness; i++)
        {
#ifdef DEBUG_HEURISTICS
            logBuffer << "SumRareScore feature " << i << "of " << pu[ipu].richness <<"\n";
#endif

            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;

#ifdef DEBUG_HEURISTICS
            logBuffer << "SumRareScore SMamount "<<SM[ism].amount<<" target " <<spec[isp].target<<" specamount "<<spec[isp].amount<<" rare "<<Rare[isp]<<"\n";
#endif

            if (SM[ism].amount &&
                (spec[isp].target > spec[isp].amount ||
                 (spec[isp].sepdistance && spec[isp].separation < 3)))
                rarescore += RareScore(isp, ipu, puno, spec, SM, R, connections, pu, cm, clumptype, thread) / Rare[isp];
        }

#ifdef DEBUG_HEURISTICS
    logBuffer  << "SumRareScore end\n";
#endif

    return (rarescore);
} // Sum Rare Score

// Set Abundances
void SetAbundance(int puno,vector<double> &Rare,vector<spustuff> &pu,vector<spu> &SM)
{
    int i, j, ism, isp;

    for (i = 0; i < puno; i++)
        if (pu[i].richness)
            for (j = 0; j < pu[i].richness; j++)
            {
                ism = pu[i].offset + j;
                isp = SM[ism].spindex;
                Rare[isp] += SM[ism].amount;
            }
} // Set Abundance

// Irreplaceability For site for species
double Irreplaceability(int ipu,int isp, vector<double> &Rare,vector<spustuff> &pu,vector<spu> &SM, vector<sspecies> &spec) 
{
    double buffer, effamount;

    buffer = Rare[isp] < spec[isp].target ? 0 : Rare[isp] - spec[isp].target;
    if (spec[isp].amount > spec[isp].target)
        return (0);
    effamount = returnAmountSpecAtPu(pu[ipu], SM, isp).second;

    return (buffer < effamount ? 1 : effamount / buffer);
}

// Product Irreplaceability for a single site
double ProdIrr(int ipu,vector<double> &Rare,vector<spustuff> &pu,vector<spu> &SM, vector<sspecies> &spec)
{
    int i, ism, isp;
    double product = 1;

    if (pu[ipu].richness)
        for (i = 0; i < pu[ipu].richness; i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            if (SM[ism].amount && (spec[isp].target - spec[isp].amount) > 0)
                product *= (1 - Irreplaceability(ipu, isp, Rare, pu, SM, spec));
        }

    return (1 - product);
} // Product Irreplaceability

double SumIrr(int ipu,vector<double> &Rare,vector<spustuff> &pu,vector<spu> &SM, vector<sspecies> &spec)
{
    int i, ism, isp;
    double sum = 0;

    if (pu[ipu].richness)
        for (i = 0; i < pu[ipu].richness; i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            if (SM[ism].amount && (spec[isp].target - spec[isp].amount) > 0)
                sum += (Irreplaceability(ipu, isp, Rare, pu, SM, spec));
        }

    return (sum);
} // Sum Irreplaceability

// Main Heuristic Engine
void Heuristics(int spno,int puno,vector<spustuff> &pu,vector<sconnections> &connections,
        vector<int> &R, double cm, vector<sspecies> &spec, vector<spu> &SM, scost &reserve,
        double costthresh, double tpf1,double tpf2, int imode,int clumptype, int thread, stringstream& logBuffer)
// imode = 1: 2: 3: 4:
// imode = 5: 6: Prod Irreplaceability, 7: Sum Irreplaceability
{
    int i, bestpu;
    double bestscore, currscore;
    scost change;
    vector<double> Rare;
    uniform_real_distribution<double> float_range(0.0, 1.0);

#ifdef DEBUG_HEURISTICS
    logBuffer << "Heuristics start mode " << imode << "\n";
#endif

    // Irreplacability
    if (imode >= 6 && imode <= 7)
    {
        Rare.resize(spno);
        SetAbundance(puno, Rare, pu, SM);

#ifdef DEBUG_HEURISTICS
        logBuffer << "Heuristics 6 or 7\n";
#endif
    }

    if (imode >= 2 && imode <= 5) // Rareness Setups
    {
        Rare.resize(spno);
        SetRareness(puno, spno, Rare, R, pu, SM, logBuffer);

#ifdef DEBUG_HEURISTICS
        logBuffer << "Heuristics 2 to 5 after SetRareness\n";
#endif
    }

    do
    {
        bestpu = 0;
        bestscore = 0;
        for (i = 0; i < puno; i++)
        {
            if (!R[i]) // Only look for new PUS
            {
                // Set the score for the given Planning Unit
                currscore = 1; // null if no other mode set
                if (imode == 0)
                    currscore = GreedyScore(i, puno, spno, spec, SM, connections, R, pu, cm, clumptype, thread);
                if (imode == 1)
                {
                    computeChangeScore(-1, i, spno, puno, pu, connections, spec, SM, R, cm, 1, change, reserve,
                                       costthresh, tpf1, tpf2, 1, clumptype, thread);
                    currscore = change.total;
                }
                if (imode == 2)
                {
                    currscore = MaxRareScore(i, puno, spec, SM, R, connections, pu, cm, Rare, clumptype, thread);
                }
                if (imode == 3)
                {
                    currscore = BestRareScore(i, puno, spec, SM, R, connections, pu, cm, Rare, clumptype, thread);
                }
                if (imode == 4)
                {
                    currscore = AveRareScore(i, puno, spec, SM, R, connections, pu, cm, Rare, clumptype, thread);
                }
                if (imode == 5)
                {

#ifdef DEBUG_HEURISTICS
                    logBuffer << "Heuristics pu " << i << " of " << puno << "\n";
#endif

                    currscore = SumRareScore(i, puno, spec, SM, R, connections, pu, cm, Rare, clumptype, thread, logBuffer);

#ifdef DEBUG_HEURISTICS
                    logBuffer << "Heuristics after SumRareScore\n";
#endif
                }
                if (imode == 6)
                {
                    currscore = -ProdIrr(i, Rare, pu, SM, spec);
                }
                if (imode == 7)
                {
                    currscore = -SumIrr(i, Rare, pu, SM, spec);
                }

                currscore *= float_range(rngEngine) * 0.001 + 1.0;
                if (!costthresh || pu[i].cost + reserve.cost <= costthresh)
                    if (currscore < bestscore)
                    {
                        bestpu = i;
                        bestscore = currscore;
                    } // is this better (ie negative) than bestscore?
            }         // I've looked through each pu to find best

            if (bestscore)
            {
                computeChangeScore(-1, bestpu, spno, puno, pu, connections, spec, SM, R, cm, 1, change, reserve,
                                   costthresh, tpf1, tpf2, 1, clumptype, thread);
                doChange(bestpu, puno, R, reserve, change, pu, SM, spec, connections, 1, clumptype, thread, logBuffer);

                // Different Heuristics might have different penalty effects
                // Old fashioned penalty and missing counting
                reserve.missing = 0;
                for (i = 0; i < spno; i++)
                {
                    if (spec[i].amount < spec[i].target)
                    {
                        reserve.missing++;
                    }
                    else if (spec[i].sepdistance && spec[i].separation < 3)
                        reserve.missing++;
                    // Species missing
                } // checking to see who I am missing
            }     // Add Pu as long as I've found one

            if (bestscore)
            {
                displayProgress2("P.U. %i PUs %i score %.6f Cost %.1f Connection %.1f Missing %i",
                                 bestpu, reserve.pus, bestscore, reserve.cost, reserve.connection, reserve.missing);
                if (fProb1D == 1)
                    displayProgress2(" Probability1D %.1f\n", reserve.probability1D);
                if (fProb2D == 1)
                    displayProgress2(" Probability2D %.1f\n", reserve.probability2D);
                displayProgress2(" Penalty %.1f\n", reserve.penalty);
            }
        }
    } while (bestscore); // Repeat until all good PUs have been added

    reserve.total = reserve.cost + reserve.connection + reserve.penalty + reserve.probability1D + reserve.probability2D;

#ifdef DEBUG_HEURISTICS
    logBuffer << "Heuristics end\n";
#endif
} // Heuristics

} // marxan