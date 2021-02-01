
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include "output.hpp"
#include "computation.hpp"

namespace marxan {

    // Constant maps for file writing/headers
    map<int, string> clumptypeMap = {
       {0, "Clumping - default step function"},
       {1, "Clumping - two level step function."},
       {2, "Clumping - rising benefit function"}
    };

    map<int, string> runoptsMap = {
       {0,"Annealing and Heuristic"},
       {1,"Annealing and Iterative Improvement"},
       {2,"Annealing and Both"},
       {3,"Heuristic only"},
       {4,"Iterative Improvement only"},
       {5,"Heuristic and Iterative Improvement"}
    };

    map<int, string> heurotypeMap = {
       {0,"Richness"},
       {1,"Greedy"},
       {2,"Maximum Rarity"},
       {3,"Best Rarity"},
       {4,"Average Rarity"},
       {5,"Summation Rarity"},
       {6,"Product Irreplaceability"},
       {7,"Summation Irreplaceability"}
    };

    // Helper function to return delimiter and file pointer for a file.
    pair<char, FILE*> GetFileAndDelimiter(string filename, int delimiterMode, bool append = false) {
        FILE* fp;
        char sDelimiter;

        fp = append ? fopen(filename.c_str(), "a") : fopen(filename.c_str(), "w");
        if (!fp)
            displayErrorMessage("Cannot save output to %s \n", filename.c_str());

        if (delimiterMode > 1)
            sDelimiter = ',';
        else
            sDelimiter = '\t';

        return pair<char, FILE*>(sDelimiter, fp);
    }

    // debug output for probability 1D
    void writeProb1DDebugTable(int spno, string savename, const vector<double>& ExpectedAmount1D, const vector<double>& VarianceInExpectedAmount1D, const vector<sspecies>& spec)
    {
        FILE* fp;
        int i, iHeavisideStepFunction;
        double rZ, rRawP, rP, rShortfallPenalty;

        fp = fopen(savename.c_str(), "w");

        fprintf(fp, "SPID,EA,VIEA,Z,RawP,ptarget1d,HeavisideSF,ShortfallP,P\n");

        for (i = spno - 1; i >= 0; i--)
        {
            computeProbMeasures(VarianceInExpectedAmount1D[i], spec[i].target, spec[i].ptarget1d, ExpectedAmount1D[i],
                rZ, rRawP, iHeavisideStepFunction, rShortfallPenalty, rP);

            fprintf(fp, "%i,%f,%f,%f,%f,%f,%i,%f,%f\n", spec[i].name, ExpectedAmount1D[i], VarianceInExpectedAmount1D[i], rZ, rRawP, spec[i].ptarget1d, iHeavisideStepFunction, rShortfallPenalty, rP);
        }

        fclose(fp);
    }

    // debug output for probability 2D
    void writeProb2DDebugTable(int spno, string savename, const vector<double>& ExpectedAmount2D, const vector<double>& VarianceInExpectedAmount2D, vector<sspecies> spec)
    {
        FILE* fp;
        int i, iHeavisideStepFunction;
        double rZ, rRawP, rP, rShortfallPenalty;

        fp = fopen(savename.c_str(), "w");

        fprintf(fp, "SPID,CA,EA,VIEA,Z,RawP,ptarget2d,HeavisideSF,ShortfallP,P\n");

        for (i = spno - 1; i >= 0; i--)
        {
            computeProbMeasures(VarianceInExpectedAmount2D[i], spec[i].target, spec[i].ptarget2d, ExpectedAmount2D[i],
                rZ, rRawP, iHeavisideStepFunction, rShortfallPenalty, rP);

            fprintf(fp, "%i,%f,%f,%f,%f,%f,%f,%i,%f,%f\n",
                spec[i].name, spec[i].amount, ExpectedAmount2D[i], VarianceInExpectedAmount2D[i], rZ, rRawP, spec[i].ptarget2d, iHeavisideStepFunction, rShortfallPenalty, rP);
        }

        fclose(fp);
    }

    // debug output for probability 1D
    void writeProb1DDetailDebugTable(string savename, int puno, int spno, const vector<spustuff>& pu, const vector<spu>& SM, const vector<int>& R)
    {
        FILE* fp;
        int i, ipu, ism, isp;
        double rAmount;

        fp = fopen(savename.c_str(), "w");

        fprintf(fp, "PUID,R,richness,PROB");
        for (i = 1; i <= spno; i++)
            fprintf(fp, ",A%i", i);
        for (i = 1; i <= spno; i++)
            fprintf(fp, ",EA%i", i);
        for (i = 1; i <= spno; i++)
            fprintf(fp, ",VIEA%i", i);
        fprintf(fp, "\n");

        for (ipu = puno - 1; ipu >= 0; ipu--)
        {
            vector<double> AMOUNT(spno, 0.0), EA(spno, 0.0), VIEA(spno, 0.0);

            if (pu[ipu].richness)
            {
                for (i = 0; i < pu[ipu].richness; i++)
                {
                    ism = pu[ipu].offset + i;
                    isp = SM[ism].spindex;
                    rAmount = SM[ism].amount;
                    AMOUNT[isp] = rAmount;
                    if (rAmount)
                    {
                        if (R[ipu] == 1 || R[ipu] == 2)
                        {
                            EA[isp] = rAmount * (1 - pu[ipu].prob);
                            VIEA[isp] = rAmount * rAmount * pu[ipu].prob * (1 - pu[ipu].prob);
                        }
                    }
                }
            }

            fprintf(fp, "%i,%i,%i,%f", 100 - ipu, R[ipu], pu[ipu].richness, pu[ipu].prob);
            for (i = spno - 1; i >= 0; i--)
                fprintf(fp, ",%f", AMOUNT[i]);
            for (i = spno - 1; i >= 0; i--)
                fprintf(fp, ",%f", EA[i]);
            for (i = spno - 1; i >= 0; i--)
                fprintf(fp, ",%f", VIEA[i]);
            fprintf(fp, "\n");
        }

        fclose(fp);
    }

    // debug output for probability 2D
    void writeProb2DDetailDebugTable(string savename, int puno, const vector<spustuff>& pu, const vector<spu>& SM, const vector<int>& R)
    {
        FILE* fp;
        int i, ipu, ism, isp;
        double AMOUNT[9], EA[9], VIEA[9], PROB[9], rAmount;

        fp = fopen(savename.c_str(), "w");

        fprintf(fp, "PUID,R,richness,P1,P2,P3,P4,P5,P6,P7,P8,P9,A1,A2,A3,A4,A5,A6,A7,A8,A9,EA1,EA2,EA3,EA4,EA5,EA6,EA7,EA8,EA9,VIEA1,VIEA2,VIEA3,VIEA4,VIEA5,VIEA6,VIEA7,VIEA8,VIEA9\n");

        for (ipu = puno - 1; ipu >= 0; ipu--)
        {
            for (i = 0; i < 9; i++)
            {
                AMOUNT[i] = 0;
                EA[i] = 0;
                VIEA[i] = 0;
                PROB[i] = 0;
            }

            if (pu[ipu].richness)
            {
                for (i = 0; i < pu[ipu].richness; i++)
                {
                    ism = pu[ipu].offset + i;
                    isp = SM[ism].spindex;
                    rAmount = SM[ism].amount;
                    AMOUNT[isp] = rAmount;
                    PROB[isp] = SM[ism].prob;
                    if (rAmount)
                    {
                        if (R[ipu] == 1 || R[ipu] == 2)
                        {
                            EA[isp] = rAmount * PROB[isp];
                            VIEA[isp] = rAmount * rAmount * PROB[isp] * (1 - PROB[isp]);
                        }
                    }
                }
            }

            fprintf(fp, "%i,%i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                107 - ipu, R[ipu], pu[ipu].richness,
                PROB[8], PROB[7], PROB[6], PROB[5], PROB[4], PROB[3], PROB[2], PROB[1], PROB[0],
                AMOUNT[8], AMOUNT[7], AMOUNT[6], AMOUNT[5], AMOUNT[4], AMOUNT[3], AMOUNT[2], AMOUNT[1], AMOUNT[0],
                EA[8], EA[7], EA[6], EA[5], EA[4], EA[3], EA[2], EA[1], EA[0],
                VIEA[8], VIEA[7], VIEA[6], VIEA[5], VIEA[4], VIEA[3], VIEA[2], VIEA[1], VIEA[0]);
        }

        fclose(fp);
    }

    // display startup message: the program title and authors
    void displayStartupMessage(void)
    {
        printf("        %s \n\n   Spatial Prioritization via Annealing\n\n", sVersionString.c_str());
        printf("   Coded by Ian Ball, modified by Matthew Watts\n");
        printf("   Written by Ian Ball and Hugh Possingham\n\n");
        printf("%s\n%s\n%s\n\n", sIanBallEmail.c_str(), sHughPossinghamEmail.c_str(), sMattWattsEmail.c_str());
        printf("   Marxan website\n\n");
        printf("%s\n\n", sMarxanWebSite.c_str());
    }

    // display shutdown message when verbosity > 0
    void displayShutdownMessage(chrono::high_resolution_clock::time_point start)
    {
        if (verbosity > 0)
        {
            printf("\n");
            displayTimePassed(start);
            printf("\n              The End \n");
            if (savelog)
            {
                fprintf(fsavelog, "\n              The End \n");
            }
        }
    }

    // write a sync file: tell calling app we've completed a run
    void writeSecondarySyncFileRun(int iSyncRun)
    {
        FILE* fsync;
        char sSyncFileName[80];

        sprintf(sSyncFileName, "sync%i", iSyncRun);

        fsync = fopen(sSyncFileName, "w");
        fprintf(fsync, sSyncFileName, "%s");
        fclose(fsync);
    }

    // write a sync file: tell calling app we've completed all runs
    void writeSecondarySyncFile(void)
    {
        FILE* fsync;

        fsync = fopen("sync", "w");
        fprintf(fsync, "sync");
        fclose(fsync);
    }

    // displays an error message for any verbosity
    // program is then terminated
    void displayErrorMessage(string sMess, ...)
    {
        extern jmp_buf jmpbuf;
        va_list args;

        va_start(args, sMess);
        vprintf(sMess.c_str(), args);
        if (savelog) vfprintf(fsavelog, sMess.c_str(), args);
        va_end(args);
        longjmp(jmpbuf, 1);
    }

    // displays a warning message when verbosity > 0
    void displayWarningMessage(string sMess, ...)
    {
        va_list args;

        if (verbosity > 0)
        {
            va_start(args, sMess);
            vprintf(sMess.c_str(), args);
            if (savelog) vfprintf(fsavelog, sMess.c_str(), args);
            va_end(args);
        }

    }

    // display a progress message when verbosity > 0
    void displayProgress(string sMess, ...)
    {
        va_list args;

        if (verbosity > 0)
        {
            va_start(args, sMess);
            vprintf(sMess.c_str(), args);
            if (savelog)
            {
                vfprintf(fsavelog, sMess.c_str(), args);
                fflush(fsavelog);
            }
            va_end(args);
        }
    }

    // display a progress message when verbosity > 1
    void displayProgress1(string sMess, ...)
    {
        va_list args;

        if (verbosity > 1)
        {
            va_start(args, sMess);
            vprintf(sMess.c_str(), args);
            if (savelog)
            {
                vfprintf(fsavelog, sMess.c_str(), args);
                fflush(fsavelog);
            }
            va_end(args);
        }
    }

    // display a progress message when verbosity > 2 :(was 2, now 5)
    void displayProgress2(string sMess, ...)
    {
        va_list args;

        if (verbosity > 5)
        {
            va_start(args, sMess);
            vprintf(sMess.c_str(), args);
            if (savelog)
            {
                vfprintf(fsavelog, sMess.c_str(), args);
                fflush(fsavelog);
            }
            va_end(args);
        }
    }

    // display a progress message when verbosity > 3 :(was 3, now 5)
    void displayProgress3(string sMess, ...)
    {
        va_list args;

        if (verbosity > 5)
        {
            va_start(args, sMess);
            vprintf(sMess.c_str(), args);
            if (savelog)
            {
                vfprintf(fsavelog, sMess.c_str(), args);
                fflush(fsavelog);
            }
            va_end(args);
        }
    }

    // create a trace file
    void createTraceFile(void)
    {
        FILE* fdebugtrace;

        if (verbosity > 2)
        {
            if (fnames.savedebugtracefile)
            {
                fdebugtrace = fopen(sTraceFileName.c_str(), "w");
                fflush(fdebugtrace);
                fclose(fdebugtrace);
            }
        }
    }

    // append message to a trace file when verbosity > 2
    void appendTraceFile(string sMess, ...)
    {
        FILE* fdebugtrace;
        va_list args;

        if (verbosity > 2)
        {
            if (fnames.savedebugtracefile)
            {
                va_start(args, sMess);

                fdebugtrace = fopen(sTraceFileName.c_str(), "a");
                vfprintf(fdebugtrace, sMess.c_str(), args);
                fclose(fdebugtrace);

                va_end(args);
            }
        }
    }

    // create a debug file
    void createDebugFile(string sFileName, string sHeader, const sfname& fnames)
    {
        FILE* fdebugtrace;
        string writename = fnames.outputdir + sFileName;

        fdebugtrace = fopen(writename.c_str(), "w");

        fprintf(fdebugtrace, sHeader.c_str(), "%s");
        fflush(fdebugtrace);
        fclose(fdebugtrace);
    }

    // append message to a debug file
    void appendDebugFile(string sFileName, string sLine, const sfname& fnames)
    {
        FILE* fdebugtrace;
        string writename = fnames.outputdir + sFileName;

        fdebugtrace = fopen(writename.c_str(), "a");

        fprintf(fdebugtrace, sLine.c_str(), "%s");
        fclose(fdebugtrace);
    }

    // display how many seconds since program started
    void displayTimePassed(chrono::high_resolution_clock::time_point start)
    {
        auto end = chrono::high_resolution_clock::now();
        int itemp = chrono::duration_cast<std::chrono::seconds>(end - start).count();

        printf("Time passed so far is ");
        if (itemp >= 60 * 60)
        {
            printf(" %i hour%c,%i min%c and %i secs \n",
                itemp / 3600, ((itemp / 3600 == 1) ? ' ' : 's'),
                (itemp / 60) % 60, ((itemp / 60 == 1) ? ' ' : 's'), itemp % 60);
        }
        else {
            if (itemp >= 60)
                printf(" %i min%c and %i secs \n", itemp / 60, ((itemp / 60 == 1) ? ' ' : 's'), itemp % 60);
            else
                printf("%i secs \n", itemp);
        }

        if (savelog)
        {
            fprintf(fsavelog, "Time passed so far is ");
            if (itemp >= 60 * 60)
            {
                fprintf(fsavelog, " %i hour%c,%i min%c and %i secs \n",
                    itemp / 3600, ((itemp / 3600 == 1) ? ' ' : 's'),
                    (itemp / 60) % 60, ((itemp / 60 == 1) ? ' ' : 's'), itemp % 60);
            }
            else {
                if (itemp >= 60)
                    fprintf(fsavelog, " %i min%c and %i secs \n", itemp / 60, ((itemp / 60 == 1) ? ' ' : 's'), itemp % 60);
                else
                    fprintf(fsavelog, "%i secs \n", itemp);
            }
        }
    }

    // create a log file, or reset a log file
    void createLogFile(int my_savelog, string my_savelogname)
    {
        if (savelog)
        { // close and delete old savelog info
            fclose(fsavelog);
        }

        savelog = my_savelog;

        if (savelog)
        {
            savelogname = my_savelogname;
            // Try to open file and complain if it don't work
            fsavelog = fopen(savelogname.c_str(), "w");
            if (fsavelog == NULL)
            {
                savelog = 0;
                displayErrorMessage("Error: Cannot save to log file %s \n", savelogname.c_str());
            }  // open failed

            // Header printing
            fprintf(fsavelog, "        %s \n\n   Spatial Prioritization via Annealing\n\n", sVersionString.c_str());
            fprintf(fsavelog, "   Coded by Ian Ball, modified by Matthew Watts\n");
            fprintf(fsavelog, "   Written by Ian Ball and Hugh Possingham\n\n");
            fprintf(fsavelog, "%s\n%s\n%s\n\n", sIanBallEmail.c_str(), sHughPossinghamEmail.c_str(), sMattWattsEmail.c_str());
            fprintf(fsavelog, "   Marxan website\n\n");
            fprintf(fsavelog, "%s\n\n", sMarxanWebSite.c_str());
        } // save log has just been turned on
    }


    // write an asymmetric connection file
    void writeAsymmetricConnectionFile(int puno, const vector<sconnections>& connections, const vector<spustuff>& pu, sfname fnames)
    {
        int i;
        FILE* fp;
        string writename = fnames.outputdir + "debug_asymmetric_connectivity.csv";

        if ((fp = fopen(writename.c_str(), "w")) == NULL)
        {
            displayProgress1("Warning: Cannot create file %s", writename.c_str());
        }

        fprintf(fp, "idA,idB,connectionorigon\n");
        for (i = 0; i < puno; i++)
        {
            for (const sneighbour& p : connections[i].first)
                fprintf(fp, "%i,%i,%i,%lf\n", pu[i].id, pu[p.nbr].id, p.connectionorigon, p.cost);
        }

        fclose(fp);
    }

    // write an output file from the loaded sparse matrix
    void writeSparseMatrix(int iSMno, int puno, vector<spustuff>& PU, vector<sspecies>& spec, vector<spu>& SM, sfname& fnames)
    {
        FILE* fp;
        string writename = fnames.inputdir + "sm.csv";

        if ((fp = fopen(writename.c_str(), "w")) == NULL)
            displayErrorMessage("cannot create PUvSpecies file %s\n", writename.c_str());

        fputs("species,pu,amount,prob\n", fp);
        for (int i = puno - 1; i >= 0; i--)
        {
            if (PU[i].richness > 0)
            {
                for (int j = 0; j < PU[i].richness; j++)
                    fprintf(fp, "%i,%i,%g,%g\n",
                        spec[SM[PU[i].offset + j].spindex].name,
                        PU[i].id,
                        SM[PU[i].offset + j].amount,
                        SM[PU[i].offset + j].prob);
            }
        }

        fclose(fp);
    }

    double computeTotalConnection2(int puno, const vector<int>& R, const vector<sconnections>& connections) {
        double connectiontemp = 0.0;
        for (int i = 0; i < puno; i++)
        {
            if (R[i] == 1 || R[i] == 2)
            {
                connectiontemp += ConnectionCost2(connections[i], R, 1, 0, 1, asymmetricconnectivity, fOptimiseConnectivityIn);
            }
        }
        return connectiontemp;
    }

    // write a summary file
    void writeSummary(string savename, const vector<string>& summaries, int imode)
    {
        pair<char, FILE*> fileInfo = GetFileAndDelimiter(savename, imode);
        FILE* fp = fileInfo.second; // Imode = 1, REST output, Imode = 2, Arcview output
        char sDelimiter = fileInfo.first;

        fprintf(fp, "\"Run_Number\"%c\"Score\"%c\"Cost\"%c\"Planning_Units\"%c\"Connectivity\"%c\"Connectivity_Total\"%c\"Connectivity_In\"%c\"Connectivity_Edge\"%c\"Connectivity_Out\"%c\"Connectivity_In_Fraction\"%c\"Penalty\"",
            sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter);
        if (fProb1D == 1)
            fprintf(fp, "%c\"Probability1D\"", sDelimiter);
        if (fProb2D == 1)
            fprintf(fp, "%c\"Probability2D\"", sDelimiter);
        fprintf(fp, "%c\"Shortfall\"%c\"Missing_Values\"%c\"MPM\"\n", sDelimiter, sDelimiter, sDelimiter);

        for (const string& line : summaries) {
            fprintf(fp, "%s", line.c_str());
        }

        fclose(fp);
    } // writeSummary

    // caculate and format summary text for a run
    string computeSummary(int puno, int spno, const vector<int>& R, const vector<sspecies>& spec, const scost& reserve,
        int itn, double misslevel, int imode)
    {
        stringstream ss;
        char del = imode > 1 ? ',' : '\t';
        int ino = 0, isp = 0;
        double shortfall, connectiontemp, rMPM, rConnectivityFraction,
            rConnectivityTotal = 0, rConnectivityIn = 0, rConnectivityEdge = 0, rConnectivityOut = 0;

        // Ouput the Summary Statistics
        for (int i = 0; i <= puno - 1; i++)
            if (R[i] == 1 || R[i] == 2)
                ino++;

        isp = computeRepresentationMISSLEVEL(spno, spec, misslevel, shortfall, rMPM);

#ifdef DEBUG_COUNTMISSING
        appendTraceFile("computeSummary shortfall %g\n", shortfall);
#endif

        computeConnectivityIndices(rConnectivityTotal, rConnectivityIn, rConnectivityEdge, rConnectivityOut,
            puno, R, connections);

        if (rConnectivityTotal > 0)
            rConnectivityFraction = rConnectivityIn / rConnectivityTotal;
        else
            rConnectivityFraction = 0;

        connectiontemp = computeTotalConnection2(puno, R, connections);

        ss << itn << del << reserve.total << del << reserve.cost << del << ino << del <<
            connectiontemp << del << rConnectivityTotal << del << rConnectivityIn << del << rConnectivityEdge << del <<
            rConnectivityOut << del << rConnectivityFraction << del << reserve.penalty;

        if (fProb1D == 1)
            ss << del << reserve.probability1D;
        if (fProb2D == 1)
            ss << del << reserve.probability2D;
        ss << del << shortfall << del << isp << del << rMPM << "\n";

        return ss.str();
    } // computeSummary

    // write the contents of the spec data structure.
    // used to validate if the input spec file has been read as intended.
    void writeSpec(int spno, const vector<sspecies>& spec, string savename)
    {
        FILE* fp = GetFileAndDelimiter(savename, 0).second;

        fprintf(fp, "id,name,target,prop,type,spf,target2,targetocc,sepdistance,sepnum,ptarget1d,ptarget2d\n");

        for (int i = 0; i < spno; i++)
            fprintf(fp, "%i,%s,%f,%f,%i,%f,%f,%i,%f,%i,%f,%f\n",
                spec[i].name,
                spec[i].sname.c_str(), spec[i].target,
                spec[i].prop, spec[i].type,
                spec[i].spf, spec[i].target2,
                spec[i].targetocc, spec[i].sepdistance,
                spec[i].sepnum, spec[i].ptarget1d,
                spec[i].ptarget2d);

        fclose(fp);
    }

    // write the contents of the pu data structure. used to validate if the input pu file has been read as intended
    void writePu(int puno, const vector<spustuff>& pu, string savename)
    {
        FILE* fp = GetFileAndDelimiter(savename, 0).second; // 0 is placeholder

        fprintf(fp, "id,status,cost,prob,xloc,yloc\n");

        for (int i = 0; i < puno; i++)
            fprintf(fp, "%i,%i,%f,%f,%f,%f\n",
                pu[i].id, pu[i].status, pu[i].cost, pu[i].prob, pu[i].xloc, pu[i].yloc);

        fclose(fp);
    }

    // write the penalty calculated for each species
    void writePenalty(int spno, const vector<sspecies>& spec, string savename, int iOutputType)
    {
        pair<char, FILE*> fileInfo = GetFileAndDelimiter(savename, iOutputType);
        FILE* fp = fileInfo.second; // Imode = 1, REST output, Imode = 2, Arcview output
        char sDelimiter = fileInfo.first;

        fprintf(fp, "spid%cpenalty\n", sDelimiter);

        // Ouput the Summary Statistics
        for (int i = 0; i < spno; i++)
            fprintf(fp, "%i%c%g\n", spec[i].name, sDelimiter, spec[i].penalty);

        fclose(fp);
    }

    // write the set of planning units used to calculate penalty
    void writePenaltyPlanningUnits(int puno, const vector<spustuff>& pu, const vector<int>& Rtemp, string savename, int iOutputType)
    {
        pair<char, FILE*> fileInfo = GetFileAndDelimiter(savename, iOutputType);
        FILE* fp = fileInfo.second; // Imode = 1, REST output, Imode = 2, Arcview output
        char sDelimiter = fileInfo.first;

        fprintf(fp, "puid%cR\n", sDelimiter);

        // Ouput the Summary Statistics
        for (int i = 0; i < puno; i++)
            fprintf(fp, "%i%c%i\n", pu[i].id, sDelimiter, Rtemp[i]);

        fclose(fp);
    }

    // create a solutions matrix file
    void createSolutionsMatrix(int puno, const vector<spustuff>& pu, string savename_ism, int iOutputType, int iIncludeHeaders)
    {
        pair<char, FILE*> fileInfo = GetFileAndDelimiter(savename_ism, iOutputType);
        FILE* fp = fileInfo.second;
        char sDelimiter = fileInfo.first;

        if (iIncludeHeaders == 1)
        {
            fprintf(fp, "SolutionsMatrix");

            for (int i = (puno - 1); i > (-1); i--)
                fprintf(fp, "%cP%i", sDelimiter, pu[i].id);

            fprintf(fp, "\n");
        }

        fclose(fp);
    }

    // append an entry to a solutions matrix file
    void appendSolutionsMatrix(int iRun, int puno, const vector<int>& R, string savename, int iOutputType, int iIncludeHeaders)
    {
        pair<char, FILE*> fileInfo = GetFileAndDelimiter(savename, iOutputType, true);
        FILE* fp = fileInfo.second;
        char sDelimiter = fileInfo.first;
        int i, iStatus;

        if (iIncludeHeaders == 1)
        {
            fprintf(fp, "S%i%c", iRun, sDelimiter);
        }

        for (i = (puno - 1); i > (-1); i--)
        {
            if (i < (puno - 1))
                fprintf(fp, "%c", sDelimiter);

            iStatus = R[i];
            if (R[i] == 3)
                iStatus = 0;
            if (R[i] == 2)
                iStatus = 1;

            fprintf(fp, "%i", iStatus);
        }

        fprintf(fp, "\n");
        fclose(fp);
    }

    // create a solution file: output_r0001.csv, output_best.csv
    void writeSolution(int puno, const vector<int>& R, const vector<spustuff>& pu, string savename, int imode, const sfname& fnames)
    {
        int i;
        FILE* fp = GetFileAndDelimiter(savename, imode).second; /* Imode = 1, REST output, Imode = 2, Arcview output */

        if (imode == 3)
        {
            fprintf(fp, "PUID,%s\n", fnames.bestfieldname.c_str());
        }
        else {
            if (imode == 2)
                fprintf(fp, "\"planning_unit\",\"solution\"\n");
        }

        for (i = puno - 1; i > -1; i--)
        {
            if (R[i] == 1 || R[i] == 2)
            {
                fprintf(fp, "%i", pu[i].id);
                if (imode > 1)
                    fprintf(fp, ",1");
                fprintf(fp, "\n");
            }
            else {
                fprintf(fp, "%i,0\n", pu[i].id);
            }
        }
        fclose(fp);
    }

    // write scenario file: a text file with input parameters
    void writeScenario(int puno, int spno, double prop, double cm,
        sanneal& anneal, int seedinit, long int repeats, int clumptype,
        int runopts, int heurotype, double costthresh, double tpf1, double tpf2,
        string savename)
    {
        FILE* fp;
        string temp;
        fp = fopen(savename.c_str(), "w");
        if (!fp)
            displayErrorMessage("Cannot save output to %s \n", savename.c_str());

        fprintf(fp, "Number of Planning Units %i\n", puno);
        fprintf(fp, "Number of Conservation Values %i\n", spno);
        fprintf(fp, "Starting proportion %.2f\n", prop);
        fprintf(fp, "Connection modifier %.2f\n\n", cm);

        // print clump type
        fprintf(fp, "%s\n", clumptypeMap[clumptype].c_str());

        fprintf(fp, "Algorithm Used :%s\n", runoptsMap[runopts].c_str());

        if (runopts == 0 || runopts == 3 || runopts == 5)
        {
            if (heurotypeMap.find(heurotype) == heurotypeMap.end()) {
                temp = "Unkown Heuristic Type";
            }
            else {
                temp = heurotypeMap[heurotype];
            }

            fprintf(fp, "Heuristic type : %s\n", temp.c_str());
        }
        else {
            fprintf(fp, "No Heuristic used \n");
        }

        if (runopts <= 2)
        {
            fprintf(fp, "Number of iterations %ld\n", anneal.iterations);
            if (anneal.Tinit >= 0)
            {
                fprintf(fp, "Initial temperature %.2f\n", anneal.Tinit);
                fprintf(fp, "Cooling factor %.6f\n", anneal.Tcool);
            }
            else
            {
                fprintf(fp, "Initial temperature set adaptively\n");
                fprintf(fp, "Cooling factor set adaptively\n");
            }
            fprintf(fp, "Number of temperature decreases %li\n\n", anneal.Titns);
        }
        else {
            fprintf(fp, "Number of iterations N/A\nInitial temperature N/A\nCooling Factor N/A\n");
            fprintf(fp, "Number of temperature decreases N/A\n\n");
        }

        if (costthresh)
        {
            fprintf(fp, "Cost Threshold Enabled: %f\n", costthresh);
            fprintf(fp, "Threshold penalty factor A %.2f\n", tpf1);
            fprintf(fp, "Threshold penalty factor B %.2f\n\n", tpf2);
        }
        else {
            fprintf(fp, "Cost Threshold Disabled\nThreshold penalty factor A N/A\n");
            fprintf(fp, "Threshold penalty factor B N/A\n\n");
        }

        fprintf(fp, "Random Seed %i\n", seedinit);
        fprintf(fp, "Number of runs %ld\n", repeats);
        fclose(fp);
    } // writeScenario

    // write a species file - the missing values file: output_mv1.csv output_mvbest.csv
    void writeSpecies(int spno, vector<sspecies>& spec, string savename, int imode, double misslevel)
    {
        int isp, iHeavisideStepFunction;
        string temp = "";
        double rMPM, rTestMPM, rRawP, rShortfallPenalty, placeholder;

        pair<char, FILE*> fileInfo = GetFileAndDelimiter(savename, imode);
        FILE* fp = fileInfo.second; // Imode = 1, Tab Delimitted Text output, Imode = 2, Arcview output
        char sDelimiter = fileInfo.first;

        fprintf(fp, "\"Conservation Feature\"%c\"Feature Name\"%c\"Target\"%c", sDelimiter, sDelimiter, sDelimiter);
        fprintf(fp, "\"Amount Held\"%c\"Occurrence Target \"%c\"Occurrences Held\"%c", sDelimiter, sDelimiter, sDelimiter);
        fprintf(fp, "\"Separation Target \"%c\"Separation Achieved\"%c\"Target Met\"%c\"MPM\"", sDelimiter, sDelimiter, sDelimiter);

        if (fProb1D == 1)
            fprintf(fp, "%cptarget1d%cEA1D%cVIEA1D%cZ1D%crawP1D%cheavisideSF1D%cshortfallP1D%cP1D", sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter);
        if (fProb2D == 1)
            fprintf(fp, "%cptarget2d%cEA2D%cVIEA2D%cZ2D%crawP2D%cheavisideSF2D%cshortfallP2D%cP2D", sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter);

        fprintf(fp, "\n");

        for (isp = 0; isp < spno; isp++)
        {
            rMPM = 1;

            fprintf(fp, "%i%c%s%c", spec[isp].name, sDelimiter, spec[isp].sname.c_str(), sDelimiter);
            fprintf(fp, "%lf%c%lf%c%i%c%i%c", spec[isp].target, sDelimiter, spec[isp].amount, sDelimiter,
                spec[isp].targetocc, sDelimiter, spec[isp].occurrence, sDelimiter);
            fprintf(fp, "%i%c%i", spec[isp].sepnum, sDelimiter, spec[isp].separation);

            if (spec[isp].target)
            {
                temp = "yes";
                if (spec[isp].amount / spec[isp].target < misslevel)
                    temp = "no";

                rTestMPM = spec[isp].amount / spec[isp].target;
                if (rTestMPM < rMPM)
                    rMPM = rTestMPM;
            }
            if (spec[isp].targetocc)
            {
                temp = "yes";
                if (spec[isp].occurrence / spec[isp].targetocc < misslevel)
                    temp = "no";

                rTestMPM = spec[isp].occurrence / spec[isp].targetocc;
                if (rTestMPM < rMPM)
                    rMPM = rTestMPM;
            }
            if (spec[isp].sepnum)
            {
                temp = "yes";
                if (spec[isp].separation / spec[isp].sepnum < misslevel)
                    temp = "no";
            }
            fprintf(fp, "%c%s", sDelimiter, temp.c_str());
            fprintf(fp, "%c%lf", sDelimiter, rMPM);

            if (fProb1D == 1)
            {
                computeProbMeasures(spec[isp].variance1D, spec[isp].target, spec[isp].ptarget1d, spec[isp].expected1D,
                    spec[isp].Zscore1D, rRawP, iHeavisideStepFunction, rShortfallPenalty, placeholder);

                // "ptarget1d EA1D VIEA1D Z1D rawP1D heavisideSF1D shortfallP1D P1D"
                fprintf(fp, "%c%lf%c%lf%c%lf%c%lf%c%lf%c%i%c%lf%c%lf",
                    sDelimiter, spec[isp].ptarget1d,
                    sDelimiter, spec[isp].expected1D,
                    sDelimiter, spec[isp].variance1D,
                    sDelimiter, spec[isp].Zscore1D,
                    sDelimiter, rRawP,
                    sDelimiter, iHeavisideStepFunction,
                    sDelimiter, rShortfallPenalty,
                    sDelimiter, spec[isp].probability1D);
            }
            if (fProb2D == 1)
            {
                computeProbMeasures(spec[isp].variance2D, spec[isp].target, spec[isp].ptarget2d, spec[isp].expected2D,
                    spec[isp].Zscore2D, rRawP, iHeavisideStepFunction, rShortfallPenalty, placeholder);

                // "ptarget2d EA2D VIEA2D Z2D rawP1D heavisideSF1D shortfallP1D P2D"
                fprintf(fp, "%c%lf%c%lf%c%lf%c%lf%c%lf%c%i%c%lf%c%lf",
                    sDelimiter, spec[isp].ptarget2d,
                    sDelimiter, spec[isp].expected2D,
                    sDelimiter, spec[isp].variance2D,
                    sDelimiter, spec[isp].Zscore2D,
                    sDelimiter, rRawP,
                    sDelimiter, iHeavisideStepFunction,
                    sDelimiter, rShortfallPenalty,
                    sDelimiter, spec[isp].probability2D);
            }

            fprintf(fp, "\n");
        }

        fclose(fp);
    }  // Output missing species information with new information

    // write summed solution file output_ssoln.csv
    void writeSumSoln(int puno, const vector<int>& sumsoln, const vector<spustuff>& pu, string savename, int imode)
    {
        pair<char, FILE*> fileInfo = GetFileAndDelimiter(savename, imode);
        FILE* fp = fileInfo.second; // Imode = 1, REST output, Imode = 2, Arcview output
        char sDelimiter = fileInfo.first;

        if (imode > 1)
        {
            fprintf(fp, "\"planning_unit\",\"number\"\n");
        }

        for (int i = 0; i < puno; i++)
            fprintf(fp, "%i%c%i\n", pu[i].id, sDelimiter, sumsoln[i]);

        fclose(fp);
    }

    // write planning unit richness to a file output_richness.csv
    void writeRichness(int puno, const vector<spustuff>& pu, string savename, int iOutputType)
    {
        pair<char, FILE*> fileInfo = GetFileAndDelimiter(savename, iOutputType);
        FILE* fp = fileInfo.second;
        char sDelimiter = fileInfo.first;

        fprintf(fp, "puid%crichness\n", sDelimiter);

        for (int i = 0; i < puno; i++)
            fprintf(fp, "%i%c%i\n", pu[i].id, sDelimiter, pu[i].richness);

        fclose(fp);
    }

    // compute total area available, reserved, excluded. write it to a file if verbosity > 3.
    void computeTotalAreas(int puno, int spno, const vector<spustuff>& pu, const vector<sspecies>& spec, const vector<spu>& SM)
    {
        if (verbosity > 3)
        {
            vector<int> TotalOccurrences(spno, 0), TO_2(spno, 0), TO_3(spno, 0);
            vector<double> TotalAreas(spno, 0), TA_2(spno, 0), TA_3(spno, 0);
            FILE* TotalAreasFile;

            computeOccurrencesAndAreas(puno, pu, SM,
                TotalOccurrences, TO_2, TO_3,
                TotalAreas, TA_2, TA_3);

            TotalAreasFile = fopen("MarOptTotalAreas.csv", "w");
            fprintf(TotalAreasFile, "spname,spindex,totalarea,reservedarea,excludedarea,targetarea,totalocc,reservedocc,excludedocc,targetocc\n");
            for (int i = 0; i < spno; i++)
                fprintf(TotalAreasFile, "%i,%i,%g,%g,%g,%g,%i,%i,%i,%i\n"
                    , spec[i].name, i, TotalAreas[i], TA_2[i], TA_3[i], spec[i].target
                    , TotalOccurrences[i], TO_2[i], TO_3[i], spec[i].targetocc);
            fclose(TotalAreasFile);
        }
    }

    // compute total area available, reserved, excluded. write it to a file output_totalareas.csv
    void writeTotalAreas(int puno, int spno, const vector<spustuff>& pu, const vector<sspecies>& spec, const vector<spu>& SM, string savename, int iOutputType)
    {
        vector<int> TotalOccurrences(spno, 0), TO_2(spno, 0), TO_3(spno, 0);
        vector<double> TotalAreas(spno, 0), TA_2(spno, 0), TA_3(spno, 0);

        pair<char, FILE*> fileInfo = GetFileAndDelimiter(savename, iOutputType);
        FILE* TotalAreasFile = fileInfo.second;
        char sDelimiter = fileInfo.first;

        computeOccurrencesAndAreas(puno, pu, SM,
            TotalOccurrences, TO_2, TO_3,
            TotalAreas, TA_2, TA_3);

        fprintf(TotalAreasFile, "spname%ctotalarea%creservedarea%cexcludedarea%ctargetarea%ctotalocc%creservedocc%cexcludedocc%ctargetocc\n",
            sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter, sDelimiter);
        for (int i = 0; i < spno; i++)
            fprintf(TotalAreasFile, "%i%c%g%c%g%c%g%c%g%c%i%c%i%c%i%c%i\n",
                spec[i].name, sDelimiter, TotalAreas[i], sDelimiter, TA_2[i], sDelimiter, TA_3[i], sDelimiter,
                spec[i].target, sDelimiter, TotalOccurrences[i], sDelimiter, TO_2[i], sDelimiter, TO_3[i], sDelimiter, spec[i].targetocc);
        fclose(TotalAreasFile);
    }

    // write vector R (status of each planning unit) to file. debug aid for annealing algorithms
    void writeR(int iMessage, string sMessage, int puno, const vector<int>& R, const vector<spustuff>& pu, const sfname& fnames)
    {
        FILE* fp;
        string messagebuffer = sMessage + to_string(iMessage);

        appendTraceFile("writeR %i start\n", iMessage);

        string writename = fnames.inputdir + "debugR_" + messagebuffer + ".csv";
        if ((fp = fopen(writename.c_str(), "w")) == NULL)
            displayErrorMessage("cannot create writeR file %s\n", writename.c_str());

        // write header row
        fprintf(fp, "puid,R\n");

        for (int i = 0; i < puno; i++)
        {
            fprintf(fp, "%i,%i\n", pu[i].id, R[i]);
        }

        fclose(fp);

        appendTraceFile("writeR %i end\n", iMessage);
    }

    // debug output for probability 1D
    void writeChangeProbability1DDebugTable(string savename, int iIteration, int ipu, int spno,
        const vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, int imode)
    {
        FILE* fp;
        int ism, isp;
        vector<double> AMOUNT(spno, 0), DE(spno, 0), DV(spno, 0), NE(spno, 0), NV(spno, 0), NZ(spno, 0), OZ(spno, 0);
        vector<double> OHSF(spno, 0), NHSF(spno, 0), OSFP(spno, 0), NSFP(spno, 0), NRP(spno, 0), ORP(spno, 0);

        if (pu[ipu].richness)
        {
            for (int i = 0; i < pu[ipu].richness; i++)
            {
                ism = pu[ipu].offset + i;
                isp = SM[ism].spindex;
                AMOUNT[isp] = SM[ism].amount;
                if (AMOUNT[isp])
                {
                    DE[isp] = imode * AMOUNT[isp] * (1 - pu[ipu].prob);
                    NE[isp] = spec[isp].expected1D + DE[isp];

                    DV[isp] = imode * AMOUNT[isp] * AMOUNT[isp] * pu[ipu].prob * (1 - pu[ipu].prob);
                    NV[isp] = spec[isp].variance1D + DV[isp];

                    if (NV[isp] > 0)
                        NZ[isp] = (spec[isp].target - NE[isp]) / sqrt(NV[isp]);
                    else
                        NZ[isp] = 4;

                    if (NZ[isp] >= 0)
                        NRP[isp] = utils::probZUT(NZ[isp]);
                    else
                        NRP[isp] = 1 - utils::probZUT(-1 * NZ[isp]);

                    if (spec[isp].variance1D > 0)
                        OZ[isp] = (spec[isp].target - spec[isp].expected1D) / sqrt(spec[isp].variance1D);
                    else
                        OZ[isp] = 4;

                    if (OZ[isp] >= 0)
                        ORP[isp] = utils::probZUT(OZ[isp]);
                    else
                        ORP[isp] = 1 - utils::probZUT(-1 * OZ[isp]);

                    if (spec[isp].ptarget1d > NRP[isp])
                        NHSF[isp] = 1;

                    if (spec[isp].ptarget1d > ORP[isp])
                        OHSF[isp] = 1;

                    if (spec[isp].ptarget1d > 0)
                    {
                        NSFP[isp] = (spec[isp].ptarget1d - NRP[isp]) / spec[isp].ptarget1d;
                        OSFP[isp] = (spec[isp].ptarget1d - ORP[isp]) / spec[isp].ptarget1d;
                    }
                }
            }
        }

        fp = fopen(savename.c_str(), "w");
        fprintf(fp, "PUID,PROB,SPID,AMOUNT,DELTAEXPECTED,DELTAVARIANCE,OLDEXPECTED,NEWEXPECTED,OLDVARIANCE,NEWVARIANCE,NEWZ,OLDZ,NEWRawP,OLDRawP,NEWHEAVISIDE,OLDHEAVISIDE,NEWSHORTFALL,OLDSHORTFALL\n");
        for (int i = spno - 1; i >= 0; i--)
        {
            fprintf(fp, "%i,%f,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                pu[ipu].id, pu[ipu].prob, spec[i].name, AMOUNT[i],
                DE[i], DV[i], spec[i].expected1D, NE[i], NV[i], spec[i].variance1D, NZ[i], OZ[i], NRP[i], ORP[i],
                NHSF[i], OHSF[i], NSFP[i], OSFP[i]);
        }
        fclose(fp);
    }

    // debug output for probability 2D
    void writeChangeProbability2DDebugTable(string savename, int iIteration, int ipu, int spno,
        const vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, int imode)
    {
        FILE* fp;
        int ism, isp;
        vector<double> AMOUNT(spno, 0), DE(spno, 0), DV(spno, 0), NE(spno, 0), NV(spno, 0), NZ(spno, 0), OZ(spno, 0), NP(spno, 0), OP(spno, 0), PROB(spno, 0);
        vector<double> OHSF(spno, 0), NHSF(spno, 0), OSFP(spno, 0), NSFP(spno, 0);

        if (pu[ipu].richness)
        {
            for (int i = 0; i < pu[ipu].richness; i++)
            {
                ism = pu[ipu].offset + i;
                isp = SM[ism].spindex;
                AMOUNT[isp] = SM[ism].amount;
                PROB[isp] = SM[ism].prob;
                if (AMOUNT[isp])
                {
                    DE[isp] = imode * AMOUNT[isp] * PROB[isp];
                    NE[isp] = spec[isp].expected1D + DE[isp];

                    DV[isp] = imode * AMOUNT[isp] * AMOUNT[isp] * PROB[isp] * (1 - PROB[isp]);
                    NV[isp] = spec[isp].variance1D + DV[isp];

                    if (NV[isp] > 0)
                        NZ[isp] = (spec[isp].target - NE[isp]) / sqrt(NV[isp]);
                    else
                        NZ[isp] = 4;

                    if (NZ[isp] >= 0)
                        NP[isp] = utils::probZUT(NZ[isp]);
                    else
                        NP[isp] = 1 - utils::probZUT(-1 * NZ[isp]);

                    if (spec[isp].variance1D > 0)
                        OZ[isp] = (spec[isp].target - spec[isp].expected2D) / sqrt(spec[isp].variance2D);
                    else
                        OZ[isp] = 4;

                    if (OZ[isp] >= 0)
                        OP[isp] = utils::probZUT(OZ[isp]);
                    else
                        OP[isp] = 1 - utils::probZUT(-1 * OZ[isp]);

                    if (spec[isp].ptarget2d > NP[isp])
                        NHSF[isp] = 1;

                    if (spec[isp].ptarget2d > OP[isp])
                        OHSF[isp] = 1;

                    if (spec[isp].ptarget2d > 0)
                    {
                        NSFP[isp] = (spec[isp].ptarget2d - NP[isp]) / spec[isp].ptarget2d;
                        OSFP[isp] = (spec[isp].ptarget2d - OP[isp]) / spec[isp].ptarget2d;
                    }
                }
            }
        }

        fp = fopen(savename.c_str(), "w");
        fprintf(fp, "IPU,PROB,ISP,AMOUNT,DELTAEXPECTED,DELTAVARIANCE,OLDEXPECTED,NEWEXPECTED,OLDVARIANCE,NEWVARIANCE,NEWZ,OLDZ,NEWP,OLDP,NEWHEAVISIDE,OLDHEAVISIDE,NEWSHORTFALL,OLDSHORTFALL\n");
        for (int i = 0; i < spno; i++)
        {
            fprintf(fp, "%i,%f,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                ipu, pu[ipu].prob, i, AMOUNT[i],
                DE[i], DV[i], spec[i].expected1D, NE[i], NV[i], spec[i].variance1D, NZ[i], OZ[i], NP[i], OP[i],
                NHSF[i], OHSF[i], NSFP[i], OSFP[i]);
        }
        fclose(fp);
    }

    // debug output for probability 1D
    void writeProbData(int puno, const vector<spustuff>& pu, const sfname& fnames)
    {
        string writename, sLine, sVarVal;
        FILE* fp;

        writename = fnames.inputdir + "writeProbData.csv";
        if ((fp = fopen(writename.c_str(), "w")) == NULL)
            displayErrorMessage("probability file %s has not been found.\nAborting Program.", writename.c_str());

        fprintf(fp, "puid,prob\n");

        for (int i = 0; i < puno; i++)
        {
            fprintf(fp, "%d,%lf\n", pu[i].id, pu[i].prob);
        }

        fclose(fp);
    }

    // make a backup copy of a connection file
    void copyFile(string sInputFile, string sOutputFile)
    {
        ifstream src(sInputFile.c_str(), ios::binary);
        ofstream dst(sOutputFile.c_str(), ios::binary);

        dst << src.rdbuf();
    }

    // Display statistics for a configuration of planning unit
    // This returns a full formatted string that can be printed
    stringstream displayValueForPUs(int puno, int spno, const vector<int>& R, const scost& reserve,
        const vector<sspecies>& spec, double misslevel)
    {
        stringstream displayValue;
        int i, isp = 0;
        double shortfall, rMPM;

#ifdef DEBUG_PRINTRESVALPROB
        appendTraceFile("displayValueForPUs start\n");
#endif

        isp = computeRepresentationMISSLEVEL(spno, spec, misslevel, shortfall, rMPM);
        double connectiontemp = computeTotalConnection2(puno, R, connections);

#ifdef DEBUG_PRINTRESVALPROB
        appendTraceFile("PrintResVal missing %i connectiontemp %g\n", isp, connectiontemp);
#endif

        displayValue << "Value " << reserve.total << " Cost " << reserve.cost << " PUs " << reserve.pus << " Connection " << connectiontemp << " ";
        displayValue << "Missing " << isp << " Shortfall " << shortfall << " Penalty " << reserve.penalty << " MPM " << rMPM << "\n";

        if (fProb1D == 1)
            displayValue << " Probability1D " << reserve.probability1D << "\n";
        if (fProb2D == 1)
            displayValue << " Probability2D " << reserve.probability2D << "\n";

#ifdef DEBUG_PRINTRESVALPROB
        appendTraceFile("displayValueForPUs end\n");
#endif

        return displayValue;
    }

    // display usage information for the marxan executable
    // displayed when the command line options are not understood
    void displayUsage(string programName)
    {
        fprintf(stderr, "%s usage: %s -[o] -[c] [input file name]\n", programName.c_str(), programName.c_str());
    }

    // if a weighted connectivity file is in use, construct one. written for Maria Beger. Not documented.
    void writeWeightedConnectivityFile(const sfname& fnames)
    {
        string readname1, readname2, writename, sFileName, sWeighting, sId1, sId2, sConnection, sAsymmetric;
        FILE* fpnames, * fpInputConnection, * fpOutputConnection;
        char sLine[500];
        int iRecords, iTotalRecords, iFiles;
        double rWeighting, rConnection;
#ifdef DEBUGTRACEFILE
        char debugbuffer[200];
#endif

        // prepare file names for backing up the connection file
        readname1 = fnames.inputdir + fnames.connectionname;
        writename = fnames.inputdir + fnames.connectionname + "~";

        // back up the connection file
        copyFile(readname1, writename);

        // prepare connection file for output
        writename = readname1;
        if ((fpOutputConnection = fopen(writename.c_str(), "w")) == NULL)
        {
            displayProgress1("Warning: Connection File %s cannot be written ", fnames.connectionname.c_str());
        }

        fprintf(fpOutputConnection, "id1,id2,connection\n");

        // prepare connection file names file for input
        readname1 = fnames.inputdir + fnames.connectionfilesname;
        if ((fpnames = fopen(readname1.c_str(), "r")) == NULL)
        {
            displayProgress1("Warning: Connectivity Files Name File %s not found ", fnames.connectionfilesname.c_str());
        }

        if (fgets(sLine, 500 - 1, fpnames) == NULL)
            displayErrorMessage("Error reading connectivity file.\n");

        // loop through the connectivity files
        iTotalRecords = 0;
        iFiles = 0;
        while (fgets(sLine, 500 - 1, fpnames))
        {
            iFiles++;
            iRecords = 0;

            // read the 3 fields from the line
            sFileName = strtok(sLine, " ,;:^*\"/\t\'\\\n");
            sWeighting = strtok(NULL, " ,;:^*\"/\t\'\\\n");
            sscanf(sWeighting.c_str(), "%lf", &rWeighting);
            sAsymmetric = strtok(NULL, " ,;:^*\"/\t\'\\\n");

            if (rWeighting != 0)
            {
                // prepare current connectivity file for input
                readname2 = fnames.inputdir + sFileName;
                if ((fpInputConnection = fopen(readname2.c_str(), "r")) == NULL)
                {
                    displayProgress1("Warning: Input Connectivity  File %s not found ", sFileName.c_str());
                }

                // read the header row
                if (fgets(sLine, 500 - 1, fpInputConnection) == NULL)
                    displayErrorMessage("Error reading connectivity file.\n");

                // write the input connection file contents to the output connection file with appropriate weightings
                while (fgets(sLine, 500 - 1, fpInputConnection))
                {
                    iRecords++;
                    iTotalRecords++;

                    sId1 = strtok(sLine, " ,;:^*\"/\t\'\\\n");
                    sId2 = strtok(NULL, " ,;:^*\"/\t\'\\\n");
                    sConnection = strtok(NULL, " ,;:^*\"/\t\'\\\n");
                    sscanf(sConnection.c_str(), "%lf", &rConnection);

                    fprintf(fpOutputConnection, "%s,%s,%lf\n", sId1.c_str(), sId2.c_str(), (rConnection * rWeighting));

                    if (sAsymmetric.compare("yes") == 0)
                        fprintf(fpOutputConnection, "%s,%s,%lf\n", sId2.c_str(), sId1.c_str(), (rConnection * rWeighting));
                }

                fclose(fpInputConnection);
            }

#ifdef DEBUGTRACEFILE
            sprintf(debugbuffer, "connectivity file %s weighting %lf asymmetric >%s< records %i\n",
                sFileName, rWeighting, sAsymmetric, iRecords);
            appendTraceFile(debugbuffer);
#endif
        }

#ifdef DEBUGTRACEFILE
        sprintf(debugbuffer, "total files %i records %i\n", iFiles, iTotalRecords);
        appendTraceFile(debugbuffer);
#endif

        fclose(fpOutputConnection);
        fclose(fpnames);

        asymmetricconnectivity = 1;
    }

} // marxan