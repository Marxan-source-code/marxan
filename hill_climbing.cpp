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
#include "defines.hpp"
#include "utils.hpp"
#include "algorithms.hpp"
#include "computation.hpp"
#include "clumping.hpp"
#include "probability.hpp"
#include "input.hpp"
#include "output.hpp"

#include "score_change.hpp"
#include "hill_climbing.hpp"


namespace marxan {
    using namespace algorithms;
    using namespace utils;

    class ImportTraceSaver
    {
        public:
            ImportTraceSaver()
            {
                ttfp = nullptr;
                Rfp = nullptr;
                iRowCounter = 1;
                iRowLimit = 0;
            }

            void init(int irun, int puno, int puvalid, const string& savename, const vector<int>& R)
            {
                string tempname2, paddedRun = utils::intToPaddedString(irun, 5);
                tempname2 = savename + "_itimp_objective" + paddedRun + ".csv";
                writename = fnames.outputdir + tempname2;
                if ((ttfp = fopen(writename.c_str(), "w")) == NULL)
                    displayErrorMessage("cannot create threshold trace file %s\n", writename.c_str());
                fprintf(ttfp, "improvement,total,pus,cost,connection,penalty,change total\n");

                tempname2 = savename + "_itimp_zones" + paddedRun + ".csv";
                writename = fnames.outputdir + tempname2;
                if ((Rfp = fopen(writename.c_str(), "w")) == NULL)
                    displayErrorMessage("cannot create threshold trace file %s\n", writename.c_str());

                fprintf(Rfp, "configuration");
                for (int i = 0; i < puno; i++)
                    fprintf(Rfp, ",%i", pu[i].id);
                fprintf(Rfp, "\n0");

                for (int i = 0; i < puno; i++)
                    fprintf(Rfp, ",%i", R[i]);
                fprintf(Rfp, "\n");

                iRowCounter = 0;
                if (fnames.itimptracerows == 0)
                    iRowLimit = 0;
                else
                    iRowLimit = floor(puvalid / fnames.itimptracerows);

            }

            void append(int i, int puno, const scost& reserve, const scost& change, const vector<int>& R)
            {
                iRowCounter++;
                if (iRowCounter > iRowLimit)
                    iRowCounter = 1;

                if (iRowCounter == 1)
                {
                    fprintf(Rfp, "%i", i);

                    fprintf(ttfp, "%i,%f,%i,%f,%f,%f,%f\n"
                                , i, reserve.total
                                , reserve.pus, reserve.cost, reserve.connection, reserve.penalty,
                                change.total); // i,costthresh,pus,cost,connection,penalty

                    for (int j = 0; j < puno; j++)
                        fprintf(Rfp, ",%i", R[j]);
                }
            }

            ~ImportTraceSaver()
            {
                if(ttfp != nullptr)
                    fclose(ttfp);
                if(Rfp != nullptr)
                    fclose(Rfp);
            }
        private:
            FILE* ttfp, * Rfp;
            string writename;
            int iRowCounter, iRowLimit;
    };

    // iteratively improves a planning unit solutions
    // a descent algorithm un-reserves planning units that don't have a negative value when removed
    void hill_climbing(int puno, int spno, const vector<spustuff>& pu, const vector<sconnections>& connections,
        vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, vector<int>& R, double cm,
        scost& reserve, scost& change, double costthresh, double tpf1, double tpf2,
        int clumptype,  int irun, int iterations, string savename, stringstream& logBuffer, rng_engine& rngEngine)
    {
        int puvalid = 0,  ipu = 0, imode, ichoice;
        vector<int> iimparray;

        logBuffer << "iterativeImprovement start\n";
        // counting pu's we need to test
        for (int i = 0; i < puno; i++)
        {
            if ((R[i] < 2) && (pu[i].status < 2))
                puvalid++;
        }
        logBuffer << "iterativeImprovement puvalid " << puvalid << "\n";

        ImportTraceSaver import_trace_saver;
        if (fnames.saveitimptrace)
            import_trace_saver.init(irun, puno, puvalid, savename, R);

        if (puvalid > 0)
        {
            iimparray.resize(puvalid);

            for (int i = 0; i < puno; i++)
            {
                if ((R[i] < 2) && (pu[i].status < 2))
                {
                    iimparray[ipu] = i;
                    ipu++;
                }
            }

            logBuffer << "iterativeImprovement after array init\n";

             
            for(int itime = 1; itime <= iterations; )
            {
                // shuffle iimp array
                std::shuffle(iimparray.begin(), iimparray.end(), rngEngine);
                bool was_change = false;

                /***** Doing the improvements ****/
                for (int i = 0; i < puvalid && itime <= iterations; i++, itime++)
                {
                    ichoice = iimparray[i];

                    if ((R[ichoice] < 2) && (pu[ichoice].status < 2))
                    {
                        imode = R[ichoice] == 1 ? -1 : 1;
                        computeChangeScore(-1, ichoice, spno, puno, pu, connections, spec, SM, SM_out, R, cm, imode, change, reserve,
                            costthresh, tpf1, tpf2, 1, clumptype);
                        if (change.total < 0)
                        {
                            displayProgress2("It Imp has changed %i with change value %lf \n", ichoice, change.total);
                            doChange(ichoice, puno, R, reserve, change, pu, SM, SM_out, spec, connections, imode, clumptype, logBuffer);
                            was_change = true;
                        }   // I've just made a good change
                    }

                    if (fnames.saveitimptrace)
                        import_trace_saver.append( itime,  puno, reserve, change, R);

                } // no untested PUs left
                if(!was_change)
                {
                    logBuffer << "Converged after "<< itime << " iterations\n";
                    break;
                }
                    
            }
        }
        appendTraceFile("iterativeImprovement end\n");
    } // hill_climbing

    // iteratively improves a planning unit solutions
    // 
    void hill_climbing_two_step(int puno, int spno, const vector<spustuff>& pu, const vector<sconnections>& connections,
        vector<sspecies>& spec, const vector<spu>& SM, vector<spu_out>& SM_out, vector<int>& R, double cm,
        scost& reserve, scost& change, double costthresh, double tpf1, double tpf2,
        int clumptype,  int irun, int iterations, string savename, stringstream& logBuffer, rng_engine& rngEngine)
    {
        int puvalid = 0,  ipu = 0;
        vector<int> iimparray;

        logBuffer << "iterativeImprovement start\n";
        // counting pu's we need to test
        for (int i = 0; i < puno; i++)
        {
            if ((R[i] < 2) && (pu[i].status < 2))
                puvalid++;
        }
        logBuffer << "iterativeImprovement puvalid " << puvalid << "\n";

        ImportTraceSaver import_trace_saver;
        if (fnames.saveitimptrace)
            import_trace_saver.init(irun, puno, puvalid, savename, R);

        if (puvalid > 0)
        {
            iimparray.resize(puvalid);

            for (int i = 0; i < puno; i++)
            {
                if ((R[i] < 2) && (pu[i].status < 2))
                {
                    iimparray[ipu] = i;
                    ipu++;
                }
            }

            logBuffer << "iterativeImprovement after array init\n";

             
            for(int itime = 1; itime <= iterations; )
            {
                // shuffle iimp array
                std::shuffle(iimparray.begin(), iimparray.end(), rngEngine);
                bool was_change = false;

                // ***** Doing the improvements ****  

                for (int i0 = 0; i0 < puvalid && itime <= iterations; i0++)
                {
                    scost change0;
                    int ichoice0 = iimparray[i0];
                    //remenber old score
                    int imode0 = R[ichoice0] == 1 ? -1 : 1;
                    computeChangeScore(-1, ichoice0, spno, puno, pu, connections, spec, SM, SM_out, R, cm, imode0, change0, reserve,
                                costthresh, tpf1, tpf2, 1, clumptype);

                    doChange(ichoice0, puno, R, reserve, change, pu, SM, SM_out, spec, connections, imode0, clumptype, logBuffer);

                    for (int i1 = i0+1; i1 < puvalid && itime <= iterations; i1++, itime++)
                    {
                        scost change1;
                        int ichoice1 = iimparray[i1];

                        int imode1 = R[ichoice1] == 1 ? -1 : 1;
                        computeChangeScore(-1, ichoice1, spno, puno, pu, connections, spec, SM, SM_out, R, cm, imode1, change1, reserve,
                                costthresh, tpf1, tpf2, 1, clumptype);

                        if (change0.total + change1.total < 0)
                        {
                            displayProgress2("It Imp has changed %i with change value %lf \n", ichoice1, change.total);
                            doChange(ichoice1, puno, R, reserve, change1, pu, SM, SM_out, spec, connections, imode1, clumptype, logBuffer);
                            was_change = true;
                        }   // I've just made a good change

                        if (fnames.saveitimptrace)
                            import_trace_saver.append( itime,  puno, reserve, change, R);

                    } // no untested PUs left
                    if(!was_change)
                    {
                        //undo change
                        int imode0_reverted = imode0  == 1 ? -1 : 1;
                        computeChangeScore(-1, ichoice0, spno, puno, pu, connections, spec, SM, SM_out, R, cm, imode0_reverted, change0, reserve,
                                    costthresh, tpf1, tpf2, 1, clumptype);

                        doChange(ichoice0, puno, R, reserve, change0, pu, SM, SM_out, spec, connections, imode0_reverted, clumptype, logBuffer);
                    }

                }
            }
        }

        appendTraceFile("iterativeImprovement end\n");
    } // hill_climbing_two_step
} // namespace marxan

