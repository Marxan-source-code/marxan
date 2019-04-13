// C code for Marxan
#define DEBUGTRACEFILE
#define PROB2D

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

#include <ctype.h>
#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "marxan.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

jmp_buf jmpbuf;
long int RandSeed1;
int iMemoryUsed=0;
int fSpecPROPLoaded = 0;
int iProbFieldPresent = 0;
int iOptimisationCalcPenalties = 1;
int verbosity = 0;
int savelog;
int asymmetricconnectivity = 0;
char sVersionString[80] = "Marxan v 3.0.0";
// version 2.3 introduced multiple connectivity files and their associated weighting file
// version 2.4.3 introduced 1D and 2D probability
// version 3.0.0 is refactoring of code in 2019
char sIanBallEmail[100] = "ian.ball@aad.gov.au";
char sHughPossinghamEmail[100] = "h.possingham@uq.edu.au";
char sMattWattsEmail[100] = "m.watts@uq.edu.au";
char sMarxanWebSite[100] = "http://marxan.net";
char sTraceFileName[1000];
char sApplicationPathName[1000];
char* savelogname;
FILE* fsavelog;

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
// R vector of planning unit status for the current solution
int *R;
// R vector of planning unit status for the best solution
int *bestR;
// is marxan being called by another program?
int marxanIsSlave = 0;

// load the required function definition modules
#include "input.c"
#include "output.c"
#include "clumping.c"
#include "heuristics.c"
#include "probability.c"

// runs the loop for each "solution" marxan is generating
void executeRunLoop(long int repeats,int puno,int spno,double cm,int aggexist,double prop,int clumptype,double misslevel,
                    char savename[],double costthresh,double tpf1,double tpf2,int heurotype,int runopts,
                    int itimptype,int sumsoln[])
{
    long int i, j;
    char tempname2[500];
    int ipu, tid, nthreads;
    
    R = (int *) calloc(puno,sizeof(int));
    bestR = (int *) calloc(puno,sizeof(int));

    // for each repeat run
    for (i = 1;i <= repeats;i++)
    {
        appendTraceFile("start run loop run %ld\n",i);

        displayProgress1("\n");
        displayProgress("Run %ld ",i);

        if (runoptions.ThermalAnnealingOn)
        {
            // Annealing Setup
            if (anneal.type == 2)
            {
                appendTraceFile("before initialiseConnollyAnnealing run %i\n",i);

                initialiseConnollyAnnealing(puno,spno,pu,connections,spec,SM,cm,&anneal,aggexist,R,prop,clumptype,i);

                appendTraceFile("after initialiseConnollyAnnealing run %i\n",i);
            }
            
            if (anneal.type == 3)
            {
                appendTraceFile("before initialiseAdaptiveAnnealing run %i\n",i);

                initialiseAdaptiveAnnealing(puno,spno,prop,R,pu,connections,SM,cm,spec,aggexist,&anneal,clumptype);

                appendTraceFile("after initialiseAdaptiveAnnealing run %i\n",i);
            }
            
            displayProgress1("  Using Calculated Tinit = %.4f Tcool = %.8f \n",
                             anneal.Tinit,anneal.Tcool);
 
            anneal.temp = anneal.Tinit;
        }

        appendTraceFile("before computeReserveValue run %i\n",i);

        displayProgress1("  Creating the initial reserve \n");
        
        initialiseReserve(puno,prop,R);  // Create Initial Reserve
       
        addReserve(puno,pu,R);

        if (aggexist)
            ClearClumps(spno,spec,pu,SM);

        computeReserveValue(puno,spno,R,pu,connections,SM,cm,spec,aggexist,&reserve,clumptype);

        appendTraceFile("after computeReserveValue run %i\n",i);

        if (verbosity > 1)
        {
            displayProgress1("\n  Init:");
            displayValueForPUs(puno,spno,R,reserve,spec,misslevel);
        }
        if (verbosity > 5)
        {
            displayTimePassed();
        }

        if (runoptions.ThermalAnnealingOn)
        {
            appendTraceFile("before thermalAnnealing run %i\n",i);

            thermalAnnealing(spno,puno,connections,R,cm,spec,pu,SM,&change,&reserve,
                             repeats,i,savename,misslevel,
                             aggexist,costthresh,tpf1,tpf2,clumptype);

            appendTraceFile("after thermalAnnealing run %i\n",i);
        }

        if (runoptions.QuantumAnnealingOn)
        {
            appendTraceFile("before quantumAnnealing run %i\n",i);

            quantumAnnealing(spno,puno,connections,R,cm,spec,pu,SM,&change,&reserve,
                             repeats,i,savename,misslevel,
                             aggexist,costthresh,tpf1,tpf2,clumptype);

            appendTraceFile("after quantumAnnealing run %i\n",i);
        }

        if (runoptions.HeuristicOn)
        {
            appendTraceFile("before Heuristics run %i\n",i);

            Heuristics(spno,puno,pu,connections,R,cm,spec,SM,&reserve,
                       costthresh,tpf1,tpf2,heurotype,clumptype);

            if (verbosity > 1 && (runopts == 2 || runopts == 5))
            {
                displayProgress1("  Heuristic:");
                displayValueForPUs(puno,spno,R,reserve,spec,misslevel);
            }

            appendTraceFile("after Heuristics run %i\n",i);
        }

        if (runoptions.ItImpOn)
        {
           appendTraceFile("before iterativeImprovement run %i\n",i);

           iterativeImprovement(puno,spno,pu,connections,spec,SM,R,cm,
                                &reserve,&change,costthresh,tpf1,tpf2,clumptype,i,savename);

           if (itimptype == 3)
               iterativeImprovement(puno,spno,pu,connections,spec,SM,R,cm,
                                   &reserve,&change,costthresh,tpf1,tpf2,clumptype,i,savename);

           appendTraceFile("after iterativeImprovement run %i\n",i);

           if (aggexist)
               ClearClumps(spno,spec,pu,SM);

           if (verbosity > 1)
           {
               computeReserveValue(puno,spno,R,pu,connections,SM,cm,spec,aggexist,&reserve,clumptype);
               displayProgress1("  Iterative Improvement:");
               displayValueForPUs(puno,spno,R,reserve,spec,misslevel);
           }

        } // Activate Iterative Improvement

        appendTraceFile("before file output run %i\n",i);

        if (fnames.saverun)
        {
           if (fnames.saverun == 3)
           {
              sprintf(tempname2,"%s_r%05li.csv",savename,i%10000);
           } else {
               if (fnames.saverun == 2)
                   sprintf(tempname2,"%s_r%05li.txt",savename,i%10000);
               else
                   sprintf(tempname2,"%s_r%05li.dat",savename,i%10000);
           }
           writeSolution(puno,R,pu,tempname2,fnames.saverun,fnames);
        }

        if (fnames.savespecies && fnames.saverun)
        {
           if (fnames.savespecies == 3)
           {
              sprintf(tempname2,"%s_mv%05li.csv",savename,i%10000);
           } else {
               if (fnames.savespecies == 2)
                   sprintf(tempname2,"%s_mv%05li.txt",savename,i%10000);
               else
                   sprintf(tempname2,"%s_mv%05li.dat",savename,i%10000);
           }

           writeSpecies(spno,spec,tempname2,fnames.savespecies,misslevel);
        }

        if (fnames.savesum)
        {
            if (fnames.savesum==3)
            {
                sprintf(tempname2,"%s_sum.csv",savename);
            } else {
            if (fnames.savesum==2)
                sprintf(tempname2,"%s_sum.txt",savename);
            else
                sprintf(tempname2,"%s_sum.dat",savename);
            }
            writeSummary(puno,spno,R,spec,reserve,i,tempname2,misslevel,fnames.savesum);
        }
        
        // compute and store objective function score for this reserve system
        computeReserveValue(puno,spno,R,pu,connections,SM,cm,spec,aggexist,&change,clumptype);
        
        // remember the bestScore and bestRun
        if (i == 1)
        {
            // this is best run so far
            bestScore = change.total;
            bestRun = 1;
            // store bestR
            for (j = 1;j <= puno;j++)
            {
                bestR[j] = R[j];
            }
        } else {
            if (change.total < bestScore)
            {
                // this is best run so far
                bestScore = change.total;
                bestRun = i;
                // store bestR
                for (j = 1;j <= puno;j++)
                {
                    bestR[j] = R[j];
                }
            }
        }
        printf("run: %ld best run: %u score: %f best score: %f\n",i,bestRun,change.total,bestScore);

        if (fnames.savesolutionsmatrix)
        {
            appendTraceFile("before appendSolutionsMatrix savename %s\n",savename);

            if (fnames.savesolutionsmatrix==3)
            {
                sprintf(tempname2,"%s_solutionsmatrix.csv",savename);
            } else {
                if (fnames.savesolutionsmatrix==2)
                {
                    sprintf(tempname2,"%s_solutionsmatrix.txt",savename);
                } else {
                    sprintf(tempname2,"%s_solutionsmatrix.dat",savename);
                }

                appendSolutionsMatrix(i,puno,R,tempname2,fnames.savesolutionsmatrix,fnames.solutionsmatrixheaders);

                appendTraceFile("after appendSolutionsMatrix savename %s\n",savename);
            }
        }
        
        if (aggexist)
            ClearClumps(spno,spec,pu,SM);

        appendTraceFile("after file output run %i\n",i);
        appendTraceFile("end run %i\n",i);

        if (marxanIsSlave == 1)
            writeSlaveSyncFileRun(i);

        if (verbosity > 1)
            displayTimePassed();
    }
} // executeRunLoop

int executeMarxan(char sInputFileName[])
{
    int iSparseMatrixFileLength = 0, iSparseMatrixFileLength_sporder = 0;
    long int repeats;
    int puno,spno,gspno;
    struct sgenspec *gspec;
    double cm,prop;
    int runopts,heurotype,clumptype,itimptype;
    char savename[500],tempname2[500];
    double misslevel;
    int iseed,seedinit;
    int aggexist=0,sepexist=0;
    int *R_CalcPenalties;
    int *sumsoln;
    double costthresh,tpf1,tpf2;
    long int itemp;
    int isp;
    #ifdef DEBUGTRACEFILE
    char debugbuffer[200];
    #endif

    // Handle Error driven termination
    if (setjmp(jmpbuf))
        return 1;

    displayStartupMessage();

    readInputOptions(&cm,&prop,&anneal,
                     &iseed,&repeats,savename,&fnames,sInputFileName,
                     &runopts,&misslevel,&heurotype,&clumptype,&itimptype,&verbosity,
                     &costthresh,&tpf1,&tpf2);

    setDefaultRunOptions(runopts,&runoptions);

    sprintf(sTraceFileName,"%s_TraceFile.txt",savename);
    createTraceFile();
    appendTraceFile("%s begin execution\n\n",sVersionString);
    appendTraceFile("LoadOptions\n");

    #ifdef DEBUGCHECKCHANGE
    createDebugFile("debug_MarOpt_CheckChange.csv","ipu,puid,R,total,cost,connection,penalty,threshpen,probability\n",fnames);
    #endif

    #ifdef DEBUGCHANGEPEN
    createDebugFile("debug_MarOpt_ChangePen.csv","ipu,puid,isp,spid,cost,newamount,famount\n",fnames);
    #endif

    #ifdef DEBUGCALCPENALTIES
    createDebugFile("debug_MarZone_CalcPenalties.csv","\n",fnames);
    #endif

    if (fnames.savelog)
    {
        sprintf(tempname2,"%s_log.dat",savename);
        createLogFile(fnames.savelog,tempname2);
    }

    delta = 1e-14;  // This would more elegantly be done as a constant

    initialiseRandomSeed(iseed);
    seedinit = iseed;

    appendTraceFile("RandSeed iseed %i RandSeed1 %li\n",iseed,RandSeed1);

    // read the data files
    displayProgress1("\nEntering in the data files \n");
    displayProgress3("    Reading in the Planning Unit names \n");
    appendTraceFile("before readPlanningUnits\n");

    itemp = readPlanningUnits(&puno,&pu,fnames);

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

    itemp = readSpecies(&spno,&spec,fnames);

    appendTraceFile("after readSpecies\n");
    displayProgress1("  %i species read in \n",itemp);
    appendTraceFile("before build search arrays\n");

    // create the fast lookup tables for planning units and species names
    computeBinarySearch(puno,spno,pu,spec,&PULookup,&SPLookup);

    appendTraceFile("after build search arrays\n");

    if (fnames.savesumsoln)
        sumsoln = (int *) calloc(puno,sizeof(int));

    connections = (typeconnection *) calloc(puno,sizeof(typeconnection));

    displayProgress3("    Reading in the Connection file :\n");
    itemp = 0;
    if (strcmp("NULL",fnames.connectionname) != 0)
    {
        appendTraceFile("before readConnections\n");

        if (strcmp("NULL",fnames.connectionfilesname) != 0)
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

    readSparseMatrix(&iSparseMatrixFileLength,&SM,puno,spno,pu,PULookup,SPLookup,fnames);

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
        if (fnames.saverichness==3)
        {
            sprintf(tempname2,"%s_richness.csv",savename);
        } else {
            if (fnames.saverichness==2)
                sprintf(tempname2,"%s_richness.txt",savename);
            else
                sprintf(tempname2,"%s_richness.dat",savename);
        }

       writeRichness(puno,pu,tempname2,fnames.saverichness);
    }

    if (strcmp("NULL",fnames.matrixspordername) != 0)
    {
        appendTraceFile("before readSparseMatrixSpOrder\n");

        readSparseMatrixSpOrder(&iSparseMatrixFileLength_sporder,&SMsporder,puno,spno,PULookup,SPLookup,fnames);

        appendTraceFile("after readSparseMatrixSpOrder\n");

        #ifdef MEMDEBUG
        displayProgress1("after LoadSparseMatrix_sporder\n");
        #endif
    }

    appendTraceFile("before process block definitions\n");

    if  (strcmp("NULL",fnames.blockdefname))
    {
        displayProgress1("    Reading in the Block Definition File \n");
        readSpeciesBlockDefinition(&gspno,&gspec,fnames);
        setBlockDefinitions(gspno,spno,puno,gspec,spec,pu,SM);
    }

    setDefaultTargets(spno,spec);

    appendTraceFile("after process block definitions\n");

    appendTraceFile("before computeTotalAreas\n");
    computeTotalAreas(puno,spno,pu,spec,SM);
    appendTraceFile("after computeTotalAreas\n");

    if (fnames.savetotalareas)
    {
        if (fnames.savetotalareas==3)
        {
          sprintf(tempname2,"%s_totalareas.csv",savename);
        } else {
            if (fnames.savetotalareas==2)
                sprintf(tempname2,"%s_totalareas.txt",savename);
            else
                sprintf(tempname2,"%s_totalareas.dat",savename);
        }

        writeTotalAreas(puno,spno,pu,spec,SM,tempname2,fnames.savepenalty);
    }

    if (fSpecPROPLoaded > 0)
    {
        appendTraceFile("before computeSpecProp\n");

        // species have prop value specified
        computeSpecProp(spno,spec,puno,pu,SM);

        appendTraceFile("after computeSpecProp\n");
    }

    displayProgress2("Checking to see if there are aggregating or separating species.\n");
    for (isp=0;isp<spno;isp++)
    {
        if (spec[isp].target2>0)
            aggexist = 1;
        if (spec[isp].sepdistance > 0)
            sepexist = 1;
    }

    if (fnames.savesen)
    {
        appendTraceFile("before writeScenario\n");

        sprintf(tempname2,"%s_sen.dat",savename);
        writeScenario(puno,spno,prop,cm,anneal,seedinit,repeats,clumptype,
                      runopts,heurotype,costthresh,tpf1,tpf2,tempname2);

        appendTraceFile("after writeScenario\n");
    }

    if (verbosity > 1)
        displayTimePassed();

    // *******  Pre-processing    ************
    displayProgress1("\nPre-processing Section. \n");
    displayProgress2("    Calculating all the penalties \n");

    R_CalcPenalties = (int *) calloc(puno,sizeof(int));

    // load penalties from file if they are present
    if (strcmp("NULL",fnames.penaltyname) != 0)
    {
	    fUserPenalties = 1;

        appendTraceFile("before readPenalties\n");

	    readPenalties(spec,spno,fnames,SPLookup);

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
        applyUserPenalties(spec,spno);
	}
	else
	{
        // we are computing penalties
        if (strcmp("NULL",fnames.matrixspordername) == 0)
        {
            appendTraceFile("before CalcPenalties\n");

            // we don't have sporder matrix available, so use slow CalcPenalties method
            itemp = computePenalties(puno,spno,pu,spec,connections,SM,R_CalcPenalties,aggexist,cm,clumptype);

            appendTraceFile("after CalcPenalties\n");
        }
        else
        {
            // we have sporder matrix available, so use optimised CalcPenalties method
            if (iOptimisationCalcPenalties == 1)
            {
                appendTraceFile("before CalcPenaltiesOptimise\n");

                itemp = computePenaltiesOptimise(puno,spno,pu,spec,connections,SM,SMsporder,R_CalcPenalties,aggexist,cm,clumptype);

                appendTraceFile("after CalcPenaltiesOptimise\n");
            }
            else
            {
                appendTraceFile("before CalcPenalties\n");

                // we have optimise calc penalties switched off, so use slow CalcPenalties method
                itemp = computePenalties(puno,spno,pu,spec,connections,SM,R_CalcPenalties,aggexist,cm,clumptype);

                appendTraceFile("after CalcPenalties\n");
            }
        }
    }

    if (itemp>0)
        displayProgress("%d species cannot meet target%c.\n",itemp,itemp==1? ' ':'s');

    if (runoptions.ThermalAnnealingOn)
    {
        displayProgress2("    Calculating temperatures.\n");
        if (!anneal.Titns)
            displayErrorMessage("Initial Temperature is set to zero. Fatal Error \n");

        anneal.Tlen = anneal.iterations/anneal.Titns;
        displayProgress2("  Temperature length %ld \n",anneal.Tlen);
        displayProgress2("  iterations %ld, repeats %ld \n",anneal.iterations,repeats);
    } // Annealing Preprocessing. Should be moved to SetAnnealingOptions

    if (fnames.savepenalty)
    {
        if (fnames.savepenalty==3)
        {
            sprintf(tempname2,"%s_penalty.csv",savename);
        } else {
            if (fnames.savepenalty==2)
                sprintf(tempname2,"%s_penalty.txt",savename);
            else
                sprintf(tempname2,"%s_penalty.dat",savename);
        }

        writePenalty(spno,spec,tempname2,fnames.savepenalty);

        if (fnames.savepenalty==3)
        {
          sprintf(tempname2,"%s_penalty_planning_units.csv",savename);
        } else {
            if (fnames.savepenalty==2)
                sprintf(tempname2,"%s_penalty_planning_units.txt",savename);
            else
                sprintf(tempname2,"%s_penalty_planning_units.dat",savename);
        }

        writePenaltyPlanningUnits(puno,pu,R_CalcPenalties,tempname2,fnames.savepenalty);
    }

    free(R_CalcPenalties);

    if (fnames.savespec)
    {
        sprintf(tempname2,"%s_spec.csv",savename);
        writeSpec(spno,spec,tempname2);
    }

    if (fnames.savepu)
    {
        sprintf(tempname2,"%s_pu.csv",savename);
        writePu(puno,pu,tempname2);
    }

    // If we are in a runmode with only CalcPenalties, we stop/exit here gracefully because we are finished.
    if (runoptions.HeuristicOn == 0)
    {
        if (runoptions.ThermalAnnealingOn == 0)
        {
            if (runoptions.QuantumAnnealingOn == 0)
            {
                if (runoptions.ItImpOn == 0)
                {
                    appendTraceFile("end final file output\n");
                    appendTraceFile("\nMarxan end execution\n");
                    displayShutdownMessage();

                    if (marxanIsSlave == 1)
                        slaveExit();

                    if (aggexist)
                        ClearClumps(spno,spec,pu,SM);  // Remove these pointers for cleanliness sake

                    if (fnames.savelog)
                        createLogFile(0,NULL);  /* tidy up files */

                    exit(1);
			    }
			}
		}
	}

    if (fnames.savesolutionsmatrix)
    {
        #ifdef DEBUG_CLUSTERANALYSIS
        appendTraceFile("before sprintf savename %s\n",savename);
        #endif

        if (fnames.savesolutionsmatrix==3)
        {
          sprintf(tempname2,"%s_solutionsmatrix.csv",savename);
        } else {
            if (fnames.savesolutionsmatrix==2)
                sprintf(tempname2,"%s_solutionsmatrix.txt",savename);
            else
                sprintf(tempname2,"%s_solutionsmatrix.dat",savename);
        }

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
        sprintf(tempname2,"%s_ComputeP_AllPUsSelected_1D.csv",savename);
        ComputeP_AllPUsSelected_1D(tempname2,puno,spno,pu,SM,spec);
    }

    if (fProb2D == 1)
    {
        sprintf(tempname2,"%s_ComputeP_AllPUsSelected_2D.csv",savename);
        ComputeP_AllPUsSelected_2D(tempname2,puno,spno,pu,SM,spec);
    }
  
    executeRunLoop(repeats,puno,spno,cm,aggexist,prop,clumptype,misslevel,
                   savename,costthresh,tpf1,tpf2,heurotype,runopts,
                   itimptype,sumsoln);

    appendTraceFile("before final file output\n");

    if (fnames.savebest)
    {
        if (fnames.savebest == 3)
        {
          sprintf(tempname2,"%s_best.csv",savename);
        } else {
            if (fnames.savebest == 2)
                sprintf(tempname2,"%s_best.txt",savename);
            else
                sprintf(tempname2,"%s_best.dat",savename);
        }

        writeSolution(puno,bestR,pu,tempname2,fnames.savebest,fnames);

        appendTraceFile("Best solution is run %i\n",bestRun);
        displayProgress1("\nBest solution is run %i\n",bestRun);
    }

    if (fnames.savespecies && fnames.savebest)
    {
        if (fnames.savespecies == 3)
        {
          sprintf(tempname2,"%s_mvbest.csv",savename);
        } else {
            if (fnames.savespecies ==2)
                sprintf(tempname2,"%s_mvbest.txt",savename);
            else
                sprintf(tempname2,"%s_mvbest.dat",savename);
        }

        writeSpecies(spno,spec,tempname2,fnames.savespecies,misslevel);
    }

    if (fnames.savesumsoln)
    {
        if (fnames.savesumsoln == 3)
        {
          sprintf(tempname2,"%s_ssoln.csv",savename);
        } else {
            if (fnames.savesumsoln == 2)
                sprintf(tempname2,"%s_ssoln.txt",savename);
            else
                sprintf(tempname2,"%s_ssoln.dat",savename);
        }

        writeSumSoln(puno,sumsoln,pu,tempname2,fnames.savesumsoln);
    }

    if (fnames.savesolutionsmatrix)
    {
        if (fnames.rexecutescript)
        {
            if (fnames.savesolutionsmatrix==3)
            {
                sprintf(tempname2,"%s_solutionsmatrix.csv",savename);
            } else {
                if (fnames.savesolutionsmatrix==2)
                    sprintf(tempname2,"%s_solutionsmatrix.txt",savename);
                else
                    sprintf(tempname2,"%s_solutionsmatrix.dat",savename);
            }

            if (marxanIsSlave == 1)
                slaveExit();
        }
    }

    if (aggexist)
        ClearClumps(spno,spec,pu,SM);  // Remove these pointers for cleanliness sake

    displayShutdownMessage();

    if (fnames.savelog)
        createLogFile(0,NULL);  /* tidy up files */

    appendTraceFile("end final file output\n");
    appendTraceFile("\nMarxan end execution\n");

    return 0;
} // executeMarxan

// returns the 0-base index of a species at a planning unit, if the species doesn't occur here, returns -1
int returnIndexSpecAtPu(struct spustuff PU[], struct spu SM[], int iPUIndex, int iSpecIndex)
{
    int i;

    if (PU[iPUIndex].richness > 0)
    {
        for (i=0;i<PU[iPUIndex].richness;i++)
            if (SM[PU[iPUIndex].offset + i].spindex == iSpecIndex)
                return (PU[iPUIndex].offset + i);
    }

    return -1; // if an index is not otherwise found, we return -1 to indicate no index was found
}

// returns the amount of a species at a planning unit, if the species doesn't occur here, returns 0
double returnAmountSpecAtPu(struct spustuff PU[], struct spu SM[], int iPUIndex, int iSpecIndex)
{
    int i;

    if (PU[iPUIndex].richness > 0)
        for (i=0;i<PU[iPUIndex].richness;i++)
            if (SM[PU[iPUIndex].offset + i].spindex == iSpecIndex)
                return SM[PU[iPUIndex].offset + i].amount;

    return 0;
}

// compute proportional target for species when prop target is specified
// use the prop value from the conservation feature file to set a proportion target for species
void computeSpecProp(int spno,typesp spec[],int puno,struct spustuff pu[],struct spu SM[])
{
    // compute and set target for species with a prop value
    double totalamount;
    int isp, ipu;

    for (isp=0;isp<spno;isp++)
    {
        if (spec[isp].prop > 0)
        {
            for (ipu = 0,totalamount = 0;ipu<puno;ipu++)
                totalamount += returnAmountSpecAtPu(pu,SM,ipu,isp);
            spec[isp].target = totalamount * spec[isp].prop;
        }
    }
}

// apply settings from the block defintion file for species
void setBlockDefinitions(int gspno,int spno,int puno,struct sgenspec gspec[], struct sspecies spec[],
                         struct spustuff PU[], struct spu SM[])
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
                    for (ipu=0,totalamount =0;ipu<puno;ipu++)
                        totalamount += returnAmountSpecAtPu(PU,SM,ipu,isp);
                    spec[isp].target = totalamount * gspec[igsp].prop;
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

// set default targets for species
void setDefaultTargets(int spno, struct sspecies spec[])
{
    int isp;
    for (isp=0;isp<spno;isp++)
    {
        if (spec[isp].target <0)
            spec[isp].target = 0;
        if (spec[isp].target2 < 0)
            spec[isp].target2 = 0;
        if (spec[isp].targetocc < 0)
            spec[isp].targetocc = 0;
        if (spec[isp].sepnum < 0)
            spec[isp].sepnum = 0;
        if (spec[isp].sepdistance < 0)
            spec[isp].sepdistance = 0;
        if (spec[isp].spf < 0)
            spec[isp].spf = 1;
    }
}

// set default run options based on the selection algorithm chosen
void setDefaultRunOptions(int runopts, struct srunoptions *runoptions)
{
    if (runopts < 0)
        return; // runopts < 0 indicates that these are set in some other way

    switch (runopts)
    {
        case 0:
            (*runoptions).CalcPenaltiesOn = 1;
            (*runoptions).ThermalAnnealingOn = 1;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 1;
            (*runoptions).ItImpOn = 0;
            break;
        case 1:
            (*runoptions).CalcPenaltiesOn = 1;
            (*runoptions).ThermalAnnealingOn = 1;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 1;
            break;
        case 2:
            (*runoptions).CalcPenaltiesOn = 1;
            (*runoptions).ThermalAnnealingOn = 1;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 1;
            (*runoptions).ItImpOn = 1;
            break;
        case 3:
            (*runoptions).CalcPenaltiesOn = 0;
            (*runoptions).ThermalAnnealingOn = 0;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 1;
            (*runoptions).ItImpOn = 0;
            break;
        case 4:
            (*runoptions).CalcPenaltiesOn = 0;
            (*runoptions).ThermalAnnealingOn = 0;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 1;
            break;
        case 5:
            (*runoptions).CalcPenaltiesOn = 0;
            (*runoptions).ThermalAnnealingOn = 0;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 1;
            (*runoptions).ItImpOn = 1;
            break;
        case 6:
            (*runoptions).CalcPenaltiesOn = 1;
            (*runoptions).ThermalAnnealingOn = 1;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 0;
            break;
        case 7:
            (*runoptions).CalcPenaltiesOn = 1;
            (*runoptions).ThermalAnnealingOn = 0;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 0;
            break;
        case 8:
            (*runoptions).CalcPenaltiesOn = 0;
            (*runoptions).ThermalAnnealingOn = 1;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 0;
            break;
        case 9:
            (*runoptions).CalcPenaltiesOn = 0;
            (*runoptions).ThermalAnnealingOn = 1;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 1;
            break;
        case 10:
            (*runoptions).CalcPenaltiesOn = 0;
            (*runoptions).ThermalAnnealingOn = 0;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 1;
            break;
        case 11:
            (*runoptions).CalcPenaltiesOn = 1;
            (*runoptions).ThermalAnnealingOn = 0;
            (*runoptions).QuantumAnnealingOn = 1;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 0;
            break;
        case 12:
            (*runoptions).CalcPenaltiesOn = 1;
            (*runoptions).ThermalAnnealingOn = 0;
            (*runoptions).QuantumAnnealingOn = 1;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 1;
            break;
        case 13:
            (*runoptions).CalcPenaltiesOn = 0;
            (*runoptions).ThermalAnnealingOn = 0;
            (*runoptions).QuantumAnnealingOn = 1;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 0;
            break;
        case 14:
            (*runoptions).CalcPenaltiesOn = 0;
            (*runoptions).ThermalAnnealingOn = 0;
            (*runoptions).QuantumAnnealingOn = 1;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 1;
            break;
        default:
            (*runoptions).CalcPenaltiesOn = 0;
            (*runoptions).ThermalAnnealingOn = 0;
            (*runoptions).QuantumAnnealingOn = 0;
            (*runoptions).HeuristicOn = 0;
            (*runoptions).ItImpOn = 0;
            break;
    }
} // setDefaultRunOptions

// for pu's specified as reserved in pu.dat status, make them reserved in a pu status vector
void addReserve(int puno,struct spustuff pu[],int *R)
{
    int i;
    for (i=0;i<puno;i++)
    {
        if (pu[i].status)
            R[i] = pu[i].status; // set planning unit status to pu.dat status
    }
}

// compute initial penalties for species with a greedy algorithm.
// If species has spatial requirements then CalcPenaltyType4 is used instead
int computePenalties(int puno,int spno,struct spustuff pu[],struct sspecies spec[],
                     struct sconnections connections[],struct spu SM[],int PUtemp[],int aggexist,double cm,int clumptype)
{
    int i,j,ibest,imaxtarget,itargetocc;
    double ftarget,fbest,fbestrat,fcost,ftemp, rAmount;
    int badspecies = 0,goodspecies = 0;

    addReserve(puno,pu,PUtemp); // Adds existing reserve to PUtemp

    for (i=0;i<spno;i++)
    {
        if (spec[i].target2 || spec[i].sepnum)
        {
            j = CalcPenaltyType4(i,puno,SM,connections,spec,pu,cm,clumptype);
            badspecies += (j>0);
            goodspecies += (j<0);

            appendTraceFile("CalcPenalties spname %i penalty %g\n",spec[i].name,spec[i].penalty);

            continue;
        } // Species has aggregation requirements

        ftarget = 0;
        itargetocc = 0;
        spec[i].penalty = 0;

        for (j=0;j<puno;j++)
        {
            if (PUtemp[j] < 2)
                PUtemp[j] = 0;
            if (PUtemp[j] == 2)
            {
                ftarget += returnAmountSpecAtPu(pu,SM,j,i);
                itargetocc++;
                spec[i].penalty += computePlanningUnitValue(j,pu,connections,cm);
            }
        } // reset PUtemp and also target

        // Already adequately represented on type 2 planning unit
        if (ftarget >= spec[i].target && itargetocc >= spec[i].targetocc)
        {
            goodspecies++;
            displayProgress2("Species %i (%s) has already met target %.2f\n",
                             spec[i].name,spec[i].sname,spec[i].target);

            appendTraceFile("CalcPenalties spname %i penalty %g\n",spec[i].name,spec[i].penalty);

            continue;
        } // Target met in unremovable reserve

        do
        {
            fbest =0; imaxtarget = 0; fbestrat = 0;
            for (j=0;j<puno;j++)
            { // trying to find best pu
                rAmount = returnAmountSpecAtPu(pu,SM,j,i);
                if (PUtemp[j] == 0 && rAmount>0)
                {
                    fcost = computePlanningUnitValue(j,pu,connections,cm);
                    if (fcost == 0)
                        fcost = delta;
                    if (rAmount >= spec[i].target - ftarget && (imaxtarget == 0 || (imaxtarget == 1 && fcost < fbest)))
                    { // can I meet the target cheaply?
                        imaxtarget = 1;
                        ibest = j;
                        fbest = fcost;
                    } else {
                        if (fbestrat < rAmount/fcost)
                        { // finding the cheapest planning unit
                            fbest = fcost;
                            fbestrat = rAmount/fbest;
                            ibest = j;
                        }
                    }

                    #ifdef DEBUGCALCPENALTIES
                    appendTraceFile("CalcPenalties species %i puid %i cost %g\n",spec[i].name,pu[j].id,fcost);
                    #endif
                }  // Making sure only checking planning units not already used
            }

            if (fbest > 0)
            {
                PUtemp[ibest] = 1;
                ftarget += returnAmountSpecAtPu(pu,SM,ibest,i);
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
                           spec[i].name,spec[i].sname,spec[i].target,ftarget);
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
        ClearClumps(spno,spec,pu,SM);

    if (goodspecies)
        displayProgress1("%i species are already adequately represented.\n",goodspecies);

    return(badspecies);
}

// compute initial penalties for species with a greedy algorithm.
int computePenaltiesOptimise(int puno,int spno,struct spustuff pu[],struct sspecies spec[],
                             struct sconnections connections[],struct spu SM[],struct spusporder SMsp[],
                             int PUtemp[],int aggexist,double cm,int clumptype)
{
    int i,j,ibest,imaxtarget,itargetocc,ism,ipu, iPUsToTest;
    double ftarget,fbest,fbestrat,fcost,ftemp, rAmount, r_ibest_amount;
    int badspecies = 0,goodspecies = 0;

    appendTraceFile("CalcPenaltiesOptimise start\n");

    addReserve(puno,pu,PUtemp); // Adds existing reserve to PUtemp

    for (i=0;i<spno;i++)
    {
        appendTraceFile("CalcPenaltiesOptimise spname %i\n",spec[i].name);

        if (spec[i].target2 || spec[i].sepnum)
        {
            j = CalcPenaltyType4(i,puno,SM,connections,spec,pu,cm,clumptype);
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
                    spec[i].penalty += computePlanningUnitValue(ipu,pu,connections,cm);
                }
            }
        }

        // Already adequately represented on type 2 planning unit
        if (ftarget >= spec[i].target && itargetocc >= spec[i].targetocc)
        { // Target met in unremovable reserve
            goodspecies++;
            displayProgress2("Species %i (%s) has already met target %.2f\n",
                             spec[i].name,spec[i].sname,spec[i].target);
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
                        fcost = computePlanningUnitValue(ipu,pu,connections,cm);
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
                             spec[i].name,spec[i].sname,spec[i].target,ftarget);
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
        ClearClumps(spno,spec,pu,SM);

    if (goodspecies)
        displayProgress1("%i species are already adequately represented.\n",goodspecies);

    appendTraceFile("CalcPenaltiesOptimise end\n");

    return(badspecies);
} // computePenaltiesOptimise

// compute cost + connectivity for a single planning unit
double computePlanningUnitValue(int ipu,struct spustuff pu[],struct sconnections connections[],double cm)
{
    double theValue;

    theValue = pu[ipu].cost;
    theValue += ConnectionCost1(ipu,pu,connections,cm);

    return(theValue);
}

// compute change in the species representation for adding or removing a single planning unit or set of planning units
double computeChangePenalty(int ipu,int puno,struct sspecies spec[],struct spustuff pu[],struct spu SM[],
                          int R[],struct sconnections connections[],int imode,int clumptype,double *rShortfall)
{
    int i, ism, isp;
    double famount, fcost= 0.0,newamount,tamount;
    double rOldShortfall, rNewAmountHeld, rNewShortfall;
    #ifdef DEBUGCHANGEPEN
    char debugline[200];
    #endif

    #ifdef ANNEALING_TEST
    if (ipu == (puno-1))
          appendTraceFile("computeChangePenalty start\n");
    #endif

    *rShortfall = 0;

    if (pu[ipu].richness)
    {
        for (i=0;i<pu[ipu].richness;i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            if (SM[ism].amount)
            {
                famount = 0;
                newamount = 0; /* Shortfall */

                rOldShortfall = 0;
                rNewShortfall = 0;

                if (spec[isp].target > spec[isp].amount)
                {
                    famount = (spec[isp].target-spec[isp].amount)/spec[isp].target;
                    rOldShortfall = spec[isp].target - spec[isp].amount;
                }

                rNewAmountHeld = spec[isp].amount + (SM[ism].amount * imode);
                if (spec[isp].target > rNewAmountHeld)
                    rNewShortfall = spec[isp].target - rNewAmountHeld;
                *rShortfall += rNewShortfall - rOldShortfall;

                if(spec[isp].targetocc > spec[isp].occurrence)
                    famount += ((double) spec[isp].targetocc - (double) spec[isp].occurrence)/
                               (double) spec[isp].targetocc;

                if (spec[isp].target && spec[isp].targetocc)
                    famount /= 2;

                if (spec[isp].sepnum)
                    famount += SepPenalty(spec[isp].separation,spec[isp].sepnum);

                if (spec[isp].target2)
                {
                    /* clumping species */
                    /* New Pen 4 includes occurrences, amounts and separation target */
                    newamount = NewPenalty4(ipu,isp,puno,spec,pu,SM,R,connections,imode,clumptype);
                } else {
                    if (spec[isp].target)
                        newamount = computeSpeciesPlanningUnitPenalty(ipu,isp,spec,pu,SM,imode)/spec[isp].target;
                    if (spec[isp].targetocc)
                    {
                        tamount =  (double) (spec[isp].targetocc - spec[isp].occurrence - imode)    /
                                   (double) spec[isp].targetocc;
                        newamount += tamount<0?0:tamount;
                    }
                    if (spec[isp].target && spec[isp].targetocc)
                        newamount /= 2;
                    if (spec[isp].sepnum)
                        newamount += SepPenalty(CountSeparation2(isp,ipu,NULL,puno,R,pu,SM,spec,imode),
                                                spec[isp].sepnum);  /* I need a new function here */
                    #ifdef ANNEALING_TEST
                    if (ipu == (puno-1))
                    {
                        appendTraceFile("penalty %g spf %g newamount %g famount %g target %g amount %g\n",
                                        spec[isp].penalty,spec[isp].spf,newamount,famount,spec[isp].target,spec[isp].amount);
                    }
                    #endif
                } /* no target2 */
                fcost += spec[isp].penalty*spec[isp].spf*(newamount - famount);
            }  /** Only worry about PUs where species occurs **/

            #ifdef DEBUGCHANGEPEN
            sprintf(debugline,"%i,%i,%i,%i,%g,%g,%g\n",ipu,pu[ipu].id,isp,spec[isp].name,fcost,newamount,famount);
            appendDebugFile("debug_MarOpt_ChangePen.csv",debugline,fnames);
            #endif
        }
    }

    #ifdef ANNEALING_TEST
    if (ipu == (puno-1))
    {
        sprintf(debugbuffer,"computeChangePenalty end fcost %g\n",fcost);
        appendTraceFile(debugbuffer);
    }
    #endif
    return (fcost);
} // computeChangePenalty

// compute objective function value of a reserve system
void computeReserveValue(int puno,int spno,int *R,struct spustuff pu[],
                         struct sconnections connections[],struct spu SM[],
                         double cm, struct sspecies spec[],int aggexist,struct scost *reserve,int clumptype)
{
    int i,j;
    double ftemp;
    double *ExpectedAmount1D, *VarianceInExpectedAmount1D,
           *ExpectedAmount2D, *VarianceInExpectedAmount2D,
           rConnectivityValue;
    #ifdef DEBUG_RESERVECOST
    char debugline[200];
    #endif
    #ifdef DEBUG_PROB2D
    char sProbDebugFileName[500], sProbDebugCost[100];
    #endif
    #ifdef DEBUG_PROB1D
    char sProbDebugFileName[500], sProbDebugCost[100];
    #endif

    // init arrays
    if (fProb1D == 1)
    {
        ExpectedAmount1D = (double *) calloc(spno,sizeof(double));
        VarianceInExpectedAmount1D = (double *) calloc(spno,sizeof(double));
        for (i=0;i<spno;i++)
        {
            ExpectedAmount1D[i] = 0;
            VarianceInExpectedAmount1D[i] = 0;
        }
    }
    if (fProb2D == 1)
    {
        ExpectedAmount2D = (double *) calloc(spno,sizeof(double));
        VarianceInExpectedAmount2D = (double *) calloc(spno,sizeof(double));
        for (i=0;i<spno;i++)
        {
            ExpectedAmount2D[i] = 0;
            VarianceInExpectedAmount2D[i] = 0;
        }
    }

    reserve->pus = 0;
    reserve->cost = 0;
    reserve->penalty = 0;
    reserve->connection = 0;
    reserve->shortfall = 0;
    reserve->probability1D = 0;
    reserve->probability2D = 0;
    if (aggexist)
        SetSpeciesClumps(puno,R,spec,pu,SM,connections,clumptype);
    for (i=0;i<spno;i++)
    {
        ftemp = 0;
        if (spec[i].target > spec[i].amount)
        {
            ftemp = (spec[i].target-spec[i].amount )/ spec[i].target;

            reserve->shortfall += spec[i].target-spec[i].amount;
        }
        if (spec[i].targetocc > spec[i].occurrence)
        {
            ftemp += (double) (spec[i].targetocc - spec[i].occurrence)/(double) spec[i].targetocc;

            reserve->shortfall += spec[i].targetocc-spec[i].occurrence;
        }
        if (spec[i].target && spec[i].targetocc)
            ftemp /= 2;
        reserve->penalty += ftemp * spec[i].penalty * spec[i].spf;

        if (spec[i].sepnum)
        {
            spec[i].separation = CountSeparation2(i,0,NULL,puno,R,pu,SM,spec,0);
            reserve->penalty += SepPenalty(spec[i].separation,spec[i].sepnum) *
                                spec[i].spf*spec[i].penalty;
        }
    }

    for (j=0;j<puno;j++)
        if (R[j]==1 || R[j] == 2)
        {
            reserve->cost += pu[j].cost;
            reserve->pus += 1;
            rConnectivityValue = ConnectionCost2(j,connections,R,1,0,cm);
            reserve->connection += rConnectivityValue;

            #ifdef DEBUG_RESERVECOST
            sprintf(debugline,"puid %i connectivity %lf\n",pu[j].id,rConnectivityValue);
            appendTraceFile(debugline);
            #endif

            if (fProb1D == 1)
                ReturnProbabilityAmounts1D(ExpectedAmount1D,VarianceInExpectedAmount1D,j,puno,pu,SM);
            if (fProb2D == 1)
                ReturnProbabilityAmounts2D(ExpectedAmount2D,VarianceInExpectedAmount2D,j,puno,pu,SM);
        }

    if (fProb1D == 1)
    {
        reserve->probability1D = ComputeProbability1D(ExpectedAmount1D,VarianceInExpectedAmount1D,spno,spec);
    }
    else
        reserve->probability1D = 0;
    if (fProb2D == 1)
        reserve->probability2D = ComputeProbability2D(ExpectedAmount2D,VarianceInExpectedAmount2D,spno,spec);
    else
        reserve->probability2D = 0;

    reserve->total = reserve->cost + reserve->connection + reserve->penalty + reserve->probability1D + reserve->probability2D;

    #ifdef DEBUG_PROB1D
    sprintf(sProbDebugCost,"probability1D %f\n",reserve->probability1D);
    appendTraceFile(sProbDebugCost);

    appendTraceFile("computeReserveValue E\n");

    sprintf(sProbDebugCost,"%f",reserve->cost);
    strcpy(sProbDebugFileName,fnames.outputdir);
    strcat(sProbDebugFileName,"output_Prob1DDebug_");
    strcat(sProbDebugFileName,sProbDebugCost);
    strcat(sProbDebugFileName,".csv");
    writeProb1DDebugTable(spno,sProbDebugFileName,
                          ExpectedAmount1D,VarianceInExpectedAmount1D,spec);


    appendTraceFile("computeReserveValue F\n");

    sprintf(sProbDebugCost,"%f",reserve->cost);
    strcpy(sProbDebugFileName,fnames.outputdir);
    strcat(sProbDebugFileName,"output_Prob1DDetailDebug_");
    strcat(sProbDebugFileName,sProbDebugCost);
    strcat(sProbDebugFileName,".csv");
    writeProb1DDetailDebugTable(sProbDebugFileName,puno,spno,pu,SM,R);

    appendTraceFile("computeReserveValue G\n");
    #endif

    #ifdef DEBUG_PROB2D
    sprintf(sProbDebugCost,"probability2D %f\n",reserve->probability2D);
    appendTraceFile(sProbDebugCost);

    sprintf(sProbDebugCost,"%f",reserve->cost);
    strcpy(sProbDebugFileName,fnames.outputdir);
    strcat(sProbDebugFileName,"output_Prob2DDebug_");
    strcat(sProbDebugFileName,sProbDebugCost);
    strcat(sProbDebugFileName,".csv");

    writeProb2DDebugTable(spno,sProbDebugFileName,
                          ExpectedAmount2D,VarianceInExpectedAmount2D,spec);

    sprintf(sProbDebugCost,"%f",reserve->cost);
    strcpy(sProbDebugFileName,fnames.outputdir);
    strcat(sProbDebugFileName,"output_Prob2DDetailDebug_");
    strcat(sProbDebugFileName,sProbDebugCost);
    strcat(sProbDebugFileName,".csv");

    writeProb2DDetailDebugTable(sProbDebugFileName,puno,pu,SM,R);
    #endif
    // destroy arrays
    if (fProb1D == 1)
    {
        for (i=0;i<spno;i++)
        {
            spec[i].expected1D = ExpectedAmount1D[i];
            spec[i].variance1D = VarianceInExpectedAmount1D[i];
        }
        free(ExpectedAmount1D);
        free(VarianceInExpectedAmount1D);
    }
    if (fProb2D == 1)
    {
        for (i=0;i<spno;i++)
        {
            spec[i].expected2D = ExpectedAmount2D[i];
            spec[i].variance2D = VarianceInExpectedAmount2D[i];
        }
        free(ExpectedAmount2D);
        free(VarianceInExpectedAmount2D);
    }
} // computeReserveValue

// initialise a planning unit vector of status to random 1's and 0's
void initialiseReserve(int puno,double prop, int *R)
{
    int i;
    for (i=0;i<puno;i++)
        R[i] = returnRandomFloat() < prop ? 1:0;
}

// sets cost threshold penalty when "cost threshold" is in use
double thresholdPenalty(double tpf1,double tpf2,double timeprop)
{
    if (tpf2 < 0)
        return(tpf1);
    return(tpf1*exp(tpf2*timeprop));
}

// compute change in the objective function score for adding or removing a single planning unit
void computeChangeScore(int iIteration,int ipu,int spno,int puno,struct spustuff pu[],struct sconnections connections[],
                        struct sspecies spec[],struct spu SM[],int *R,double cm,int imode,
                        struct scost *change, struct scost *reserve,double costthresh,double tpf1, double tpf2,
                        double timeprop,int clumptype)
// imode = 1 add PU, imode = -1 remove PU
{
    double threshpen = 0;
    int threshtype = 1;  /*Debugging line. This should be input parameter not hardwired */
    double tchangeconnection,tresconnection;
    #ifdef DEBUGCHECKCHANGE
    char debugline[200];
    #endif

    change->cost = pu[ipu].cost*imode;    /* Cost of this PU on it's own */
    change->connection = ConnectionCost2(ipu,connections,R,imode,1,cm);
    if (threshtype ==1)
    {
        tchangeconnection = change->connection;
        tresconnection = reserve->connection;
        change->connection = 0;
        reserve->connection = 0;
    }

    change->penalty = computeChangePenalty(ipu,puno,spec,pu,SM,R,connections,imode,clumptype,&change->shortfall);

    if (costthresh)
    {
        // Threshold Penalty for costs
        if (reserve->cost + reserve->connection <= costthresh)
        {
            if (change->cost + change->connection + reserve->cost + reserve->connection <= costthresh)
                threshpen = 0;
            else
                threshpen = (change->cost + change->connection +
                             reserve->cost + reserve->connection - costthresh) *
                             thresholdPenalty(tpf1,tpf2,timeprop);
        }
        else
        {
            if (change->cost + change->connection + reserve->cost + reserve->connection <= costthresh)
                threshpen = (reserve->cost + reserve->connection - costthresh) *
                             thresholdPenalty(tpf1,tpf2,timeprop);
            else
                threshpen = (change->cost + change->connection) *
                             thresholdPenalty(tpf1,tpf2,timeprop);
        }
    }

    change->threshpen = threshpen;

    if (threshtype ==1)
    {
        change->connection = tchangeconnection;
        reserve->connection = tresconnection;
    }

    if (fProb1D == 1)
        change->probability1D = ChangeProbability1D(iIteration,ipu,spno,puno,spec,pu,SM,imode);
    else
        change->probability1D = 0;
    if (fProb2D == 1)
        change->probability2D = ChangeProbability2D(iIteration,ipu,spno,puno,spec,pu,SM,imode);
    else
        change->probability2D = 0;

    change->total = change->cost + change->connection + change->penalty + change->threshpen + change->probability1D + change->probability2D;

    #ifdef DEBUGCHECKCHANGE
    sprintf(debugline,"%i,%i,%i,%g,%g,%g,%g,%g,%g,%g\n",
            ipu,pu[ipu].id,R[ipu],change->total,change->cost,change->connection,change->penalty,change->threshpen,change->probability1D,change->probability2D);
    appendDebugFile("debug_MarOpt_CheckChange.csv",debugline,fnames);
    #endif
} // computeChangeScore

// compute change in the objective function score for adding or removing a set of planning units
void computeQuantumChangeScore(int spno,int puno,struct spustuff pu[],struct sconnections connections[],
                               struct sspecies spec[],struct spu SM[],int *R,double cm,
                               struct scost *change, struct scost *reserve,double costthresh,double tpf1, double tpf2,
                               double timeprop,int clumptype,int iFluctuationCount,int *PUChosen)
// imode = 1 add PU, imode = -1 remove PU
{
    // We query a whole bunch of changes in one, passed in by Quantum annealing.
    double threshpen = 0;
    int imode, i, j, threshtype = 1;  // Debugging line. This should be input parameter not hardwired */
    double tchangeconnection,tresconnection;
    #ifdef DEBUG_QA
    char debugline[200];
    #endif

    #ifdef DEBUG_QA
    appendTraceFile("computeQuantumChangeScore start iFluctuationCount %i\n",iFluctuationCount);
    #endif

    change->cost = 0;
    change->connection = 0;
    change->penalty = 0;
    change->shortfall = 0;
    change->probability1D = 0;
    change->probability2D = 0;
    j=-1;

    for (i=0;i<iFluctuationCount;i++)
    {
        do
            j++;

        while (PUChosen[j] < 1);

        imode = R[j] == 1 ? -1 : 1;

        #ifdef DEBUG_QA
        appendTraceFile("computeQuantumChangeScore ipu %i chosen %i imode %i\n",j,PUChosen[j],imode);
        #endif

        change->cost += pu[j].cost*imode;    /* Cost of this PU on it's own */
        change->connection += ConnectionCost2(j,connections,R,imode,1,cm);
        if (threshtype ==1)
        {
            tchangeconnection = change->connection;
            tresconnection = reserve->connection;
            change->connection = 0;
            reserve->connection = 0;
        }

        change->penalty += computeChangePenalty(j,puno,spec,pu,SM,R,connections,imode,clumptype,&change->shortfall);

        if (costthresh)
        {
            // Threshold Penalty for costs
            if (reserve->cost + reserve->connection <= costthresh)
            {
                if (change->cost + change->connection + reserve->cost + reserve->connection <= costthresh)
                    threshpen = 0;
                else
                    threshpen = (change->cost + change->connection +
                                 reserve->cost + reserve->connection - costthresh) *
                                 thresholdPenalty(tpf1,tpf2,timeprop);
            }
            else
            {
                if (change->cost + change->connection + reserve->cost + reserve->connection <= costthresh)
                    threshpen = (reserve->cost + reserve->connection - costthresh) *
                                 thresholdPenalty(tpf1,tpf2,timeprop);
                else
                    threshpen = (change->cost + change->connection) *
                                 thresholdPenalty(tpf1,tpf2,timeprop);
            }
        }

        change->threshpen = threshpen;

        if (threshtype ==1)
        {
            change->connection = tchangeconnection;
            reserve->connection = tresconnection;
        }

        if (fProb1D == 1)
            change->probability1D += ChangeProbability1D(-1,j,spno,puno,spec,pu,SM,imode);
        else
             change->probability1D = 0;
        if (fProb2D == 1)
            change->probability2D += ChangeProbability2D(-1,j,spno,puno,spec,pu,SM,imode);
        else
            change->probability2D = 0;
    }

    change->total = change->cost + change->connection + change->penalty + change->threshpen + change->probability1D + change->probability2D;

    #ifdef DEBUGCHECKCHANGE
    sprintf(debugline,"%i,%i,%i,%g,%g,%g,%g,%g,%g\n",j,pu[j].id,R[j],change->total,change->cost,change->connection,change->penalty,change->threshpen,change->probability);
    appendDebugFile("debug_MarOpt_CheckChange.csv",debugline,fnames);
    #endif

    #ifdef DEBUG_QA
    appendTraceFile("computeQuantumChangeScore end\n");
    #endif
} // computeQuantumChangeScore

// compute penalty for a species for changing status of a single planning unit
double computeSpeciesPlanningUnitPenalty(int ipu,int isp,struct sspecies spec[],struct spustuff pu[],struct spu SM[],int imode)
{
    double newpen;

    newpen = spec[isp].target - spec[isp].amount - returnAmountSpecAtPu(pu,SM,ipu,isp)*imode;

    if (newpen < 0)
        newpen = 0;

    return(newpen);
}

// determines if the change value for changing a single planning unit status is good
// does the change stochastically fall below the current acceptance probability?
int isGoodChange(struct scost change,double temp)
{
    return (exp(-change.total/temp)> returnRandomFloat()) ? 1 : 0;
}

// determines if the change value for changing status for a set of planning units is good
// does it stochastically fall below the current acceptance probability?
int isGoodQuantumChange(struct scost change,double rProbAcceptance)
{
    if (change.total <= 0)
        return 1;
    else
        return (rProbAcceptance > returnRandomFloat()) ? 1 : 0;
}

// change the status of a single planning unit
void doChange(int ipu,int puno,int *R,struct scost *reserve,struct scost change,
              struct spustuff pu[],struct spu SM[],struct sspecies spec[],struct sconnections connections[],
              int imode,int clumptype)
{
    int i,ism,isp;
    double rAmount;

    R[ipu] = imode == 1 ? 1 : 0;
    reserve->pus += imode;
    reserve->cost += change.cost;
    reserve->connection += change.connection;
    reserve->penalty += change.penalty;
    reserve->probability1D += change.probability1D;
    reserve->probability2D += change.probability2D;
    reserve->shortfall += change.shortfall;

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
                    AddNewPU(ipu,isp,connections,spec,pu,SM,clumptype);
                } else {
                    RemPu(ipu,isp,connections,spec,pu,SM,clumptype);
                }
                if (spec[isp].occurrence < 0)
                {
                    printf("Warning Warning ! isp %i occ %i \n",isp,spec[isp].occurrence);
                }
            } else { // No clumping species
                spec[isp].occurrence += (rAmount > 0)*imode;
                spec[isp].amount += rAmount*imode;

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
                appendTraceFile("doChange ipu %i isp %i spec.amount %g imode %i\n",
                                ipu,isp,spec[isp].amount,imode);
                #endif
            }

            if (spec[isp].sepnum>0) /* Count separation but only if it is possible that it has changed */
               if ((imode ==1 && spec[isp].separation < spec[isp].sepnum) || (imode == -1 && spec[isp].separation >1))
                   spec[isp].separation = CountSeparation2(isp,0,NULL,puno,R,pu,SM,spec,0);
        }
    }

    reserve->total = reserve->cost + reserve->connection + reserve->penalty + reserve->probability1D + reserve->probability2D;
} // doChange

// change the status of a set of planning units
void doQuantumChange(int puno,int *R,struct scost *reserve,struct scost change,
                     struct spustuff pu[],struct spu SM[],struct sspecies spec[],struct sconnections connections[],
                     int clumptype,int iFluctuationCount,int *PUChosen)
{
    // We accept a whole bunch of changes in one, passed in by Quantum annealing.
    int i,j,ipu,ism,isp,imode;
    double rAmount;

    #ifdef DEBUG_QA
    appendTraceFile("doQuantumChange start\n");
    #endif

    reserve->cost += change.cost;
    reserve->connection += change.connection;
    reserve->penalty += change.penalty;
    reserve->probability1D += change.probability1D;
    reserve->probability2D += change.probability2D;

    ipu = -1;
    for (j=0;j<iFluctuationCount;j++)
    {
        do
            ipu++;

        while (PUChosen[ipu] < 1);

        imode = R[ipu] == 1 ? -1 : 1;
        R[ipu] = imode == 1 ? 1 : 0;
        reserve->pus += imode;

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
                        AddNewPU(ipu,isp,connections,spec,pu,SM,clumptype);
                    } else {
                        RemPu(ipu,isp,connections,spec,pu,SM,clumptype);
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
                        spec[isp].separation = CountSeparation2(isp,0,NULL,puno,R,pu,SM,spec,0);
            }
        }
    }

    reserve->total = reserve->cost + reserve->connection + reserve->penalty + reserve->probability1D + reserve->probability2D;

    #ifdef DEBUG_QA
    appendTraceFile("doQuantumChange end\n");
    #endif
} // doQuantumChange

/*****************************************************/
/*********** Post Processing *************************/
/*****************************************************/

/***  Counts the number of species missing from the reserve ****/
// compute the number of species whose representation fraction is less that the MISSLEVEL parameter
int computeRepresentationMISSLEVEL(int spno,struct sspecies spec[],double misslevel,double *shortfall,double *rMinimumProportionMet)
{
    int i,isp = 0;
    double rProportionMet;

    *shortfall = 0;
    *rMinimumProportionMet = 1;
    for (i=0;i<spno;i++)
    {
        rProportionMet = 1;

        if (spec[i].target > 0)
        {
            if (spec[i].amount < spec[i].target)
            {
                *shortfall += spec[i].target - spec[i].amount;
                rProportionMet = spec[i].amount / spec[i].target;

                if (rProportionMet < *rMinimumProportionMet)
                    *rMinimumProportionMet = rProportionMet;

                #ifdef DEBUG_COUNTMISSING
                appendTraceFile("computeRepresentationMISSLEVEL i %i target %g amount %g shortfall %g\n",i,spec[i].target,spec[i].amount,*shortfall);
                #endif
            }
        }
        if (spec[i].targetocc > 0)
        {
            if (spec[i].occurrence < spec[i].targetocc)
            {
                *shortfall += spec[i].targetocc - spec[i].occurrence;
                rProportionMet = spec[i].occurrence / spec[i].targetocc;

                if (rProportionMet < *rMinimumProportionMet)
                    *rMinimumProportionMet = rProportionMet;
            }
        }
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

    #ifdef DEBUG_COUNTMISSING
    appendTraceFile("computeRepresentationMISSLEVEL shortfall %g\n",*shortfall);
    #endif

    return(isp);
}

// compute connectivity total, in, edge, out for summary report
void computeConnectivityIndices(double *rConnectivityTotal,double *rConnectivityIn,
                                double *rConnectivityEdge,double *rConnectivityOut,
                                int puno,int *R,typeconnection connections[])
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
    struct sneighbour *p;

    for (i=0;i<puno;i++)
    {
        rFixed = connections[i].fixedcost;

        *rConnectivityTotal += rFixed;

        if (R[i]==1 || R[i] == 2)
        { // add to 'in' or 'edge'
            *rConnectivityEdge += rFixed;

            p = connections[i].first;
            while (p)
            {
                if (p->nbr > i)
                {
                    if (R[p->nbr] == 1 || R[p->nbr] == 2) // add to 'in'
                        *rConnectivityIn += p->cost;
                    else  // add to 'edge'
                        *rConnectivityEdge += p->cost;

                    // add to 'total'
                    *rConnectivityTotal += p->cost;
                }

                p = p->next;
            }
        }
        else
        { // add to 'out' or 'edge'
            *rConnectivityOut += rFixed;

            p = connections[i].first;
            while (p)
            {
                if (p->nbr > i)
                {
                    if (R[p->nbr] == 1 || R[p->nbr] == 2) // add to 'edge'
                        *rConnectivityEdge += p->cost;
                    else  // add to 'out'
                        *rConnectivityOut += p->cost;

                    // add to 'total'
                    *rConnectivityTotal += p->cost;
                }

                p = p->next;
            }
        }
    }
}

void initialiseConnollyAnnealing(int puno,int spno,struct spustuff pu[],typeconnection connections[],typesp spec[],
                                 struct spu SM[],double cm, struct sanneal *anneal,int aggexist,
                                 int R[],double prop,int clumptype,int irun)
{
    long int i,ipu,imode, iOldR;
    double deltamin = 0,deltamax = 0;
    double localdelta;
    #ifdef DEBUGTRACEFILE
    char sRun[20];
    FILE *fp;
    char *writename;
    #endif

    #ifdef DEBUGTRACEFILE
    appendTraceFile("initialiseConnollyAnnealing start\n");
    if (verbosity > 4)
    {
        sprintf(sRun,"%i",irun);
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debug_maropt_initialiseConnollyAnnealing_.csv") + strlen(sRun) + 2, sizeof(char));
        strcpy(writename,fnames.outputdir);
        strcat(writename,"debug_maropt_initialiseConnollyAnnealing_");
        strcat(writename,sRun);
        strcat(writename,".csv");
        fp = fopen(writename,"w");
        if (fp==NULL)
            displayErrorMessage("cannot create debug_maropt_initialiseConnollyAnnealing file %s\n",writename);
        free(writename);
        fprintf(fp,"i,ipu,puid,old R,imode,R,total,max,min\n");
    }
    #endif

    #ifdef DEBUG_PROB1D
    appendTraceFile("initialiseConnollyAnnealing A\n");
    #endif

    localdelta = 1E-10;

    initialiseReserve(puno,prop,R);

    #ifdef DEBUG_PROB1D
    appendTraceFile("initialiseConnollyAnnealing B\n");
    #endif

    addReserve(puno,pu,R);
    if (aggexist)
        ClearClumps(spno,spec,pu,SM);

    #ifdef DEBUG_PROB1D
    appendTraceFile("initialiseConnollyAnnealing C\n");
    #endif

    computeReserveValue(puno,spno,R,pu,connections,SM,cm,spec,aggexist,&reserve,clumptype);

    #ifdef DEBUG_PROB1D
    appendTraceFile("initialiseConnollyAnnealing D\n");
    #endif

    for (i=1;i<= (*anneal).iterations/100; i++)
    {
        ipu = returnRandom(puno);
        iOldR = R[ipu];
        imode = R[ipu]==1?-1:1;

        computeChangeScore(-1,ipu,spno,puno,pu,connections,spec,SM,R,cm,imode,&change,&reserve,0,0,0,0,clumptype);
        doChange(ipu,puno,R,&reserve,change,pu,SM,spec,connections,imode,clumptype);

        if (change.total > deltamax)
            deltamax = change.total;
        if (change.total >localdelta && (deltamin < localdelta || change.total < deltamin))
            deltamin = change.total;

        if (verbosity > 4)
            fprintf(fp,"%li,%li,%i,%li,%li,%i,%g,%g,%g\n",i,ipu,pu[ipu].id,iOldR,imode,R[ipu],change.total,deltamax,deltamin);
                                                   // i,ipu,puid,R,imode,iZone,total,max,min
    }  // Run through this bit for iterations/100 times

    (*anneal).Tinit = deltamax;
    deltamin *= 0.1;

    (*anneal).Tcool = exp(log(deltamin/ (*anneal).Tinit)/(double)(*anneal).Titns);

    appendTraceFile("initialiseConnollyAnnealing end\n");
    if (verbosity > 4)
        fclose(fp);
} // initialiseConnollyAnnealing

// initialise adaptive annealing (where anneal type = 3)
void initialiseAdaptiveAnnealing(int puno,int spno,double prop,int *R,
    struct spustuff pu[],struct sconnections connections[],
    struct spu SM[],double cm,struct sspecies spec[],int aggexist,struct sanneal *anneal,int clumptype)
{
    long int i,isamples;
    double sum = 0,sum2 = 0;
    double sigma;
    struct scost cost;
    double c = 10;  /* An initial temperature acceptance number */

    isamples = 1000; /* Hardwired number of samples to take */

    for (i=0;i<isamples;i++)
    {  /* Generate Random Reserve */
        initialiseReserve(puno,prop,R);
        addReserve(puno,pu,R);
        /* Score Random reserve */
        computeReserveValue(puno,spno,R,pu,connections,SM,cm,spec,aggexist,&cost,clumptype);
        /* Add Score to Sum */
        sum += cost.total;
        sum2 += cost.total*cost.total;
    } /* Sample space iterations/100 times */

    sigma = sqrt(sum2 - pow(sum/isamples,2))/(isamples-1);

    (*anneal).Tinit = c * sigma;
    (*anneal).sigma = sigma;
    (*anneal).temp = (*anneal).Tinit;
    (*anneal).tempold = (*anneal).temp;
    (*anneal).sum = 0;
    (*anneal).sum2 = 0;

    appendTraceFile("Tinit %g Titns %li Tcool %g\n",(*anneal).Tinit,(*anneal).Titns,(*anneal).Tcool);
} // initialiseAdaptiveAnnealing

// reduce annealing temperature when anneal type = 3
void reduceTemperature(struct sanneal *anneal)
{
    double omega = 0.7; /* Control parameter */
    double sigmanew,sigmamod;
    double lambda = 0.7; /* control parameter*/

    sigmanew = ((*anneal).sum2 - pow(((*anneal).sum/(*anneal).Tlen),2))/((*anneal).Tlen-1);
    sigmamod = (1-omega)*sigmanew + omega * (*anneal).sigma *((*anneal).temp/(*anneal).tempold);
    (*anneal).tempold = (*anneal).temp;
    (*anneal).temp = exp(-lambda*(*anneal).temp/sigmamod);
    (*anneal).sigma = sigmamod;
    (*anneal).sum = 0;
    (*anneal).sum2 = 0;
}

// run simulated thermal annealing selection algorithm
void thermalAnnealing(int spno, int puno, struct sconnections connections[],int R[], double cm,
                      typesp *spec, struct spustuff pu[], struct spu SM[], struct scost *change, struct scost *reserve,
                      long int repeats,int irun,char *savename,double misslevel,
                      int aggexist,
                      double costthresh, double tpf1, double tpf2,int clumptype)
{
    long int itime = 0,ipu = -1,i,itemp,snapcount,ichanges = 0, iPreviousR,iGoodChange = 0;
    long int iRowCounter, iRowLimit;
    double rTemperature, rThreshold, rThresholdMultiplier;
    char tempname1[12],tempname2[100], sRun[20];
    FILE *fp,*ttfp,*Rfp;
    char *writename;

    appendTraceFile("thermalAnnealing start iterations %ld\n",anneal.iterations);
    if (verbosity > 4)
    {
        sprintf(sRun,"%i",irun);
        writeR(0,"after_Annealing_entered",puno,R,pu,fnames);
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debug_maropt_annealing_.csv") + strlen(sRun) + 2, sizeof(char));
        strcpy(writename,fnames.outputdir);
        strcat(writename,"debug_maropt_annealing_");
        strcat(writename,sRun);
        strcat(writename,".csv");
        if ((fp = fopen(writename,"w"))==NULL)
            displayErrorMessage("cannot create annealing file %s\n",writename);
        free(writename);
        fprintf(fp,"itime,ipu,puid,R,itemp,newR,iGoodChange,changetotal,changecost,changeconnection,changepen,temp\n");
    }

    if (fnames.saveannealingtrace)
    {
        sprintf(tempname2,"%s_anneal_objective%05i.csv",savename,irun%10000);
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen(tempname2) + 2, sizeof(char));
        strcat(writename,tempname2);
        if ((ttfp = fopen(writename,"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename);
        free(writename);
        fprintf(ttfp,"iteration,threshold,dochange,total,pus,cost,connectivity,penalty,shortfall");
        if (fProb1D == 1)
            fprintf(ttfp,",probability1D");
        if (fProb2D == 1)
            fprintf(ttfp,",probability2D");
        fprintf(ttfp,",puindex\n");

        // write iteration zero
        fprintf(ttfp,"%li,%f,%li,%f,%i,%f,%f,%f,%f",
                itime,costthresh,iGoodChange,reserve->total,
                reserve->pus,reserve->cost,reserve->connection,reserve->penalty,reserve->shortfall);
        if (fProb1D == 1)
            fprintf(ttfp,",%f",reserve->probability1D);
        if (fProb2D == 1)
            fprintf(ttfp,",%f",reserve->probability2D);
        fprintf(ttfp,",%li\n",ipu);
        // iteration,threshold,dochange,total,pus,cost,connectivity,penalty,probability

        sprintf(tempname2,"%s_anneal_zones%05i.csv",savename,irun%10000);
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen(tempname2) + 2, sizeof(char));
        strcat(writename,tempname2);
        if ((Rfp = fopen(writename,"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename);
        free(writename);
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

    for (itime = 1;itime<=anneal.iterations;itime++)
    {
        do
        { //  Select a PU at random
            ipu = returnRandom(puno);
        } while (R[ipu] > 1);

        itemp = R[ipu] == 1 ? -1 : 1;  /* Add or Remove PU ? */

        computeChangeScore(itime,ipu,spno,puno,pu,connections,spec,SM,R,cm,itemp,change,reserve,
                           costthresh,tpf1,tpf2,(double) itime/ (double) anneal.iterations,clumptype);
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
                reduceTemperature(&anneal);
            else
                anneal.temp = anneal.temp*anneal.Tcool;

            displayProgress3("time %ld temp %f Complete %ld%% currval %.4f\n",
                             itime,anneal.temp,(int)itime*100/anneal.iterations,reserve->total);
        } /* reduce temperature */

        if (fnames.savesnapsteps && !(itime % fnames.savesnapfrequency))
        {
            if (repeats > 1)
                sprintf(tempname1,"_r%05i",irun);
            else
                tempname1[0] = 0;
            if (fnames.savesnapchanges == 3)
                sprintf(tempname2,"%s_snap%st%05li.csv",savename,tempname1,++snapcount%10000);
            else
            if (fnames.savesnapchanges == 2)
                sprintf(tempname2,"%s_snap%st%05li.txt",savename,tempname1,++snapcount%10000);
            else
                sprintf(tempname2,"%s_snap%st%05li.dat",savename,tempname1,++snapcount%10000);
            writeSolution(puno,R,pu,tempname2,fnames.savesnapsteps,fnames);
        } /* Save snapshot every savesnapfreq timesteps */
        iPreviousR = R[ipu];
        if (isGoodChange(*change,anneal.temp)==1)
        {
            iGoodChange = 1;

            ++ichanges;
            doChange(ipu,puno,R,reserve,*change,pu,SM,spec,connections,itemp,clumptype);
            if (fnames.savesnapchanges && !(ichanges % fnames.savesnapfrequency))
            {
                if (repeats > 1)
                {
                    sprintf(tempname1,"_r%05i",irun);
                } else {
                    tempname1[0] = 0;
                }
                if (fnames.savesnapchanges == 3)
                {
                    sprintf(tempname2,"%s_snap%sc%05li.csv",savename,tempname1,++snapcount%10000);
                } else {
                    if (fnames.savesnapchanges == 2)
                        sprintf(tempname2,"%s_snap%sc%05li.txt",savename,tempname1,++snapcount%10000);
                    else
                        sprintf(tempname2,"%s_snap%sc%05li.dat",savename,tempname1,++snapcount%10000);
                }
                writeSolution(puno,R,pu,tempname2,fnames.savesnapchanges,fnames);
            } /* Save snapshot every savesnapfreq changes */

        } /* Good change has been made */
        else
            iGoodChange = 0;

        if (anneal.type == 3)
        {
            anneal.sum += reserve->total;
            anneal.sum2 += reserve->total*reserve->total;
        } /* Keep track of scores for averaging stuff */

        if (verbosity > 4)
            fprintf(fp,"%li,%li,%i,%li,%li,%i,%li,%f,%f,%f,%f,%f\n",
                    itime,ipu,pu[ipu].id,iPreviousR,itemp,R[ipu],iGoodChange,change->total,change->cost,change->connection,change->penalty,anneal.temp);

        if (fnames.saveannealingtrace)
        {
            iRowCounter++;
            if (iRowCounter > iRowLimit)
                iRowCounter = 1;

            if (iRowCounter == 1)
            {
                fprintf(Rfp,"%li",itime);

                fprintf(ttfp,"%li,%f,%li,%f,%i,%f,%f,%f,%f",
                        itime,costthresh,iGoodChange,reserve->total,
                        reserve->pus,reserve->cost,reserve->connection,reserve->penalty,reserve->shortfall);
                if (fProb1D == 1)
                    fprintf(ttfp,",%f",reserve->probability1D);
                if (fProb2D == 1)
                    fprintf(ttfp,",%f",reserve->probability2D);
                fprintf(ttfp,",%li\n",ipu);
                // iteration,threshold,dochange,total,pus,cost,connectivity,penalty,probability
                for (i = 0;i<puno;i++)
                    fprintf(Rfp,",%i",R[i]);

                fprintf(Rfp,"\n");
            }
        }
    } /* Run Through Annealing */

    /** Post Processing  **********/
    if (verbosity > 1)
    {
        computeReserveValue(puno,spno,R,pu,connections,SM,cm,spec,aggexist,reserve,clumptype);
        displayProgress1("  thermalAnnealing:");

        #ifdef DEBUG_PRINTRESVALPROB
        appendTraceFile("before displayValueForPUs thermalAnnealing:\n");
        #endif

        displayValueForPUs(puno,spno,R,*reserve,spec,misslevel);

        #ifdef DEBUG_PRINTRESVALPROB
        appendTraceFile("after displayValueForPUs thermalAnnealing:\n");
        #endif
    }

    if (aggexist)
        ClearClumps(spno,spec,pu,SM);

    if (verbosity > 4)
        fclose(fp);

    if (fnames.saveannealingtrace)
    {
        fclose(ttfp);
        fclose(Rfp);
    }
} // thermalAnnealing

// run simulated quantum annealing selection algorithm
void quantumAnnealing(int spno, int puno, struct sconnections connections[],int R[], double cm,
                      typesp *spec, struct spustuff pu[], struct spu SM[], struct scost *change, struct scost *reserve,
                      long int repeats,int irun,char *savename,double misslevel,
                      int aggexist,
                      double costthresh, double tpf1, double tpf2,int clumptype)
{
    long int itime,i,j,itemp,snapcount,ichanges = 0, iGoodChange;
    long int iRowCounter, iRowLimit, iFluctuationCount;
    double rFluctuationMagnitude, rThreshold, rThresholdMultiplier,
    rAcceptanceProbability;
    char tempname1[12],tempname2[100], sRun[20];
    FILE *fp;
    FILE *ttfp,*Rfp;
    char *writename, sDecayType[20];
    int *PUChosen;
    long int iTests = 0;
    long int iIterations;

    if (iQADECAYTYPE == 0)
        strcpy(sDecayType,"EXPONENTIAL");
    else
        strcpy(sDecayType,"SIGMOIDAL");

    appendTraceFile("quantumAnnealing start iterations %ld decay type %s proportion %f decay A %f decay B %f acceptance probability %f saveannealingtrace %i\n",
                          anneal.iterations,sDecayType,rQAPROP,rQADECAY,rQADECAYB,rQAACCPR,fnames.saveannealingtrace);
    if (verbosity > 4)
    {
        sprintf(sRun,"%i",irun);
        writeR(0,"after_Annealing_entered",puno,R,pu,fnames);
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debug_maropt_annealing_.csv") + strlen(sRun) + 2, sizeof(char));
        strcpy(writename,fnames.outputdir);
        strcat(writename,"debug_maropt_annealing_");
        strcat(writename,sRun);
        strcat(writename,".csv");
        if ((fp = fopen(writename,"w"))==NULL)
            displayErrorMessage("cannot create annealing file %s\n",writename);
        free(writename);
        fprintf(fp,"itime,ipu,puid,R,itemp,newR,iGoodChange,changetotal,changecost,changeconnection,changepen,temp\n");
    }

    if (fnames.saveannealingtrace)
    {
        sprintf(tempname2,"%s_anneal_objective%05i.csv",savename,irun%10000);
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen(tempname2) + 2, sizeof(char));
        strcat(writename,tempname2);
        if ((ttfp = fopen(writename,"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename);
        free(writename);
        fprintf(ttfp,"iteration,threshold,dochange,total,pus,cost,connectivity,penalty");
        if (fProb1D == 1)
            fprintf(ttfp,",probability1D");
        if (fProb2D == 1)
            fprintf(ttfp,",probability2D");
        fprintf(ttfp,",Fmag,Fcount\n");

        sprintf(tempname2,"%s_anneal_zones%05i.csv",savename,irun%10000);
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen(tempname2) + 2, sizeof(char));
        strcat(writename,tempname2);
        if ((Rfp = fopen(writename,"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename);
        free(writename);
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

    PUChosen = (int *) calloc(puno,sizeof(int));

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
            for (i = 0;i<puno;i++)
                PUChosen[i] = 0;

            for (i = 0;i<iFluctuationCount;i++)
            {
                do
                {
                    j = returnRandom(puno);

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
                                      clumptype,iFluctuationCount,PUChosen);

            // we only accept good changes
            if (fnames.savesnapsteps && !(itime % fnames.savesnapfrequency))
            { // Save snapshot every savesnapfreq timesteps
                if (repeats > 1)
                {
                    sprintf(tempname1,"_r%05i",irun);
                } else {
                    tempname1[0] = 0;
                }
                if (fnames.savesnapchanges == 3)
                {
                    sprintf(tempname2,"%s_snap%st%05li.csv",savename,tempname1,++snapcount%10000);
                } else {
                    if (fnames.savesnapchanges == 2)
                        sprintf(tempname2,"%s_snap%st%05li.txt",savename,tempname1,++snapcount%10000);
                    else
                        sprintf(tempname2,"%s_snap%st%05li.dat",savename,tempname1,++snapcount%10000);
                }
                writeSolution(puno,R,pu,tempname2,fnames.savesnapsteps,fnames);
            }
            if (isGoodQuantumChange(*change,rAcceptanceProbability)==1)
            { // Save snapshot every savesnapfreq changes
                iGoodChange = 1;

                ++ichanges;
                doQuantumChange(puno,R,reserve,*change,pu,SM,spec,connections,clumptype,iFluctuationCount,PUChosen);
                if (fnames.savesnapchanges && !(ichanges % fnames.savesnapfrequency))
                {
                    if (repeats > 1)
                    {
                        sprintf(tempname1,"_r%05i",irun);
                    } else {
                        tempname1[0] = 0;
                    }
                    if (fnames.savesnapchanges == 3)
                    {
                        sprintf(tempname2,"%s_snap%sc%05li.csv",savename,tempname1,++snapcount%10000);
                    } else {
                        if (fnames.savesnapchanges == 2)
                            sprintf(tempname2,"%s_snap%sc%05li.txt",savename,tempname1,++snapcount%10000);
                        else
                            sprintf(tempname2,"%s_snap%sc%05li.dat",savename,tempname1,++snapcount%10000);
                    }
                    writeSolution(puno,R,pu,tempname2,fnames.savesnapchanges,fnames);
                }
            } /* Good change has been made */
            else
                iGoodChange = 0;

            if (anneal.type == 3)
            { // Keep track of scores for averaging stuff
                anneal.sum += reserve->total;
                anneal.sum2 += reserve->total*reserve->total;
            }

            if (verbosity > 4)
                fprintf(fp,"%li,%li,%li,%f,%f,%f,%f,%f\n",
                        itime,itemp,iGoodChange,change->total,change->cost,change->connection,change->penalty,anneal.temp);

            if (fnames.saveannealingtrace)
            {
                iRowCounter++;
                if (iRowCounter > iRowLimit)
                    iRowCounter = 1;

                if (iRowCounter == 1)
                {
                    fprintf(Rfp,"%li",itime);

                    fprintf(ttfp,"%li,%f,%li,%f,%i,%f,%f,%f\n",
                            itime,costthresh,iGoodChange,reserve->total,
                            reserve->pus,reserve->cost,reserve->connection,reserve->penalty);
                    if (fProb1D == 1)
                        fprintf(ttfp,",%f",reserve->probability1D);
                    if (fProb2D == 1)
                        fprintf(ttfp,",%f",reserve->probability2D);
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

    free(PUChosen);

    /** Post Processing  **********/
    if (verbosity >1)
    {
        computeReserveValue(puno,spno,R,pu,connections,SM,cm,spec,aggexist,reserve,clumptype);
        displayProgress1("  quantumAnnealing:");

        #ifdef DEBUG_PRINTRESVALPROB
        appendTraceFile("before displayValueForPUs quantumAnnealing:\n");
        #endif

        displayValueForPUs(puno,spno,R,*reserve,spec,misslevel);

        #ifdef DEBUG_PRINTRESVALPROB
        appendTraceFile("after displayValueForPUs quantumAnnealing:\n");
        #endif
    }

    if (aggexist)
        ClearClumps(spno,spec,pu,SM);

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

// Big O notation: optimisation functions by Matt Watts

// helps to do a heap sort
void siftDownBinarySearch(struct binsearch numbers[], int root, int bottom, int array_size)
{
    int done, maxChild;
    typebinsearch temp;

    done = 0;
    while ((root*2 <= bottom) && (!done))
    {
        if (root*2 < array_size)
        {
            if (root*2 == bottom)
            {
                maxChild = root * 2;
            } else {
                if (numbers[root * 2].name > numbers[root * 2 + 1].name)
                {
                    maxChild = root * 2;
                } else {
                    maxChild = root * 2 + 1;
                }
            }

            if (numbers[root].name < numbers[maxChild].name)
            {
                temp = numbers[root];
                numbers[root] = numbers[maxChild];
                numbers[maxChild] = temp;
                root = maxChild;
            } else {
                done = 1;
            }
        } else {
            done = 1;
        }
    }
}

// sort a datastructure with heap sort
void heapSortBinarySearch(struct binsearch numbers[], int array_size)
{
    int i;
    typebinsearch temp;

    for (i = (array_size / 2)-1; i >= 0; i--)
        siftDownBinarySearch(numbers, i, array_size, array_size);

    for (i = array_size-1; i >= 1; i--)
    {
        temp = numbers[0];
        numbers[0] = numbers[i];
        numbers[i] = temp;
        siftDownBinarySearch(numbers, 0, i-1, array_size);
    }
}

// compute binary search arrays for looking up pu's and species fast
void computeBinarySearch(int puno, int spno, struct spustuff PU[], typesp spec[],
                               struct binsearch *PULookup[], struct binsearch *SPLookup[])
{
    int i;

    /* create the lookup arrays for planning unit and species names */
    *PULookup = (struct binsearch *) calloc(puno,sizeof(struct binsearch));
    *SPLookup = (struct binsearch *) calloc(spno,sizeof(struct binsearch));

    /* populate the lookup arrays with planning unit and species names*/
    for (i=0;i<puno;i++)
    {
        (* PULookup)[i].name = PU[i].id;
        (* PULookup)[i].index = i;
    }
    for (i=0;i<spno;i++)
    {
        (* SPLookup)[i].name = spec[i].name;
        (* SPLookup)[i].index = i;
    }

    if (verbosity > 3)
        writeBinarySearchArrays("before",fnames,puno,spno,(* PULookup),(* SPLookup));

    /* sort the lookup arrays by name */
    heapSortBinarySearch((* PULookup),puno);
    heapSortBinarySearch((* SPLookup),spno);

    if (verbosity > 3)
        writeBinarySearchArrays("after",fnames,puno,spno,(* PULookup),(* SPLookup));
}

// use binary search to find a PU index given it's id
int binarySearchPuIndex(int puno,int name, struct binsearch PULookup[])
{
    /* use a binary search to find the index of planning unit "name" */
    int iTop, iBottom, iCentre, iCount;

    iTop = 0;
    iBottom = puno-1;
    iCentre = iTop + floor(puno / 2);

    while ((iTop <= iBottom) && (PULookup[iCentre].name != name))
    {
        if (name < PULookup[iCentre].name)
        {
            iBottom = iCentre - 1;
            iCount = iBottom - iTop + 1;
            iCentre = iTop + floor(iCount / 2);
        } else {
            iTop = iCentre + 1;
            iCount = iBottom - iTop + 1;
            iCentre = iTop + floor(iCount / 2);
        }
    }
    return(PULookup[iCentre].index);
}

// use binary search to find a species index given it's id
int binarySearchSpecIndex(int spno,int name, struct binsearch SPLookup[])
{
    /* use a binary search to find the index of species "name" */
    int iTop, iBottom, iCentre, iCount;

    iTop = 0;
    iBottom = spno-1;
    iCentre = iTop + floor(spno / 2);

    while ((iTop <= iBottom) && (SPLookup[iCentre].name != name))
    {
        if (name < SPLookup[iCentre].name)
        {
            iBottom = iCentre - 1;
            iCount = iBottom - iTop + 1;
            iCentre = iTop + floor(iCount / 2);
        }
        else
        {
            iTop = iCentre + 1;
            iCount = iBottom - iTop + 1;
            iCentre = iTop + floor(iCount / 2);
        }
    }
    return(SPLookup[iCentre].index);
}

// marxan is running as a slave and has finished. create a sync file so the calling software will know marxan has finished creating the output files
// slaveExit does not deliver a message prior to exiting, but creates a file so C-Plan/Zonae Cogito/etc knows marxan has exited
void slaveExit(void)
{
    writeSlaveSyncFile();
}

// add a field name from an input file header to a list
struct snlink *storeFieldName(char **varlist,int numvars,char *sVarName,
                              struct snlink *head,char *fname)
{
    int i,foundit = 0;
    struct snlink *temp,*newlink=NULL;

    for (i=0;(i<numvars && foundit==0);i++)
    {
        if (strcmp(varlist[i],sVarName) == 0)
            foundit++;
    }

    if (head)
    {
        for (temp = head;temp;temp = temp->next)
        {
            if (strcmp(temp->name,sVarName) == 0)
                displayErrorMessage("ERROR: variable %s has been defined twice in data file %s.\n",sVarName,fname);
        }
    }

    newlink = (struct snlink *) malloc(sizeof(struct snlink));
    newlink->next = NULL;
    newlink->name = (char *) calloc(strlen(sVarName)+1,sizeof(char));
    strcpy(newlink->name,sVarName);
    return(newlink);
}

// apply the species penalties nominated in input penalties file for use in the annealing algorithms
void applyUserPenalties(typesp spec[],int spno)
{
    int i;

    for (i=0;i<spno;i++)
        spec[i].penalty = spec[i].rUserPenalty;
}

// helps to do a heap sort
void siftDownIterativeImprovement(struct iimp numbers[], int root, int bottom, int array_size)
{
    int done, maxChild;
    typeiimp temp;

    done = 0;
    while ((root*2 <= bottom) && (!done))
    {
        if (root*2 < array_size)
        {
            if (root*2 == bottom)
            {
                maxChild = root * 2;
            } else {
                if (numbers[root * 2].randomfloat > numbers[root * 2 + 1].randomfloat)
                    maxChild = root * 2;
                else
                    maxChild = root * 2 + 1;
            }

            if (numbers[root].randomfloat < numbers[maxChild].randomfloat)
            {
                 temp = numbers[root];
                 numbers[root] = numbers[maxChild];
                 numbers[maxChild] = temp;
                 root = maxChild;
            } else {
                done = 1;
            }
        } else {
            done = 1;
        }
    }
}

// sort a datastructure with heap sort
void heapSortIterativeImprovement(struct iimp numbers[], int array_size)
{
    int i;
    typeiimp temp;

    for (i = (array_size / 2)-1; i >= 0; i--)
    {
        siftDownIterativeImprovement(numbers, i, array_size, array_size);
    }

    for (i = array_size-1; i >= 1; i--)
    {
        temp = numbers[0];
        numbers[0] = numbers[i];
        numbers[i] = temp;
        siftDownIterativeImprovement(numbers, 0, i-1, array_size);
    }
}

// iteratively improves a planning unit solutions
// a descent algorithm un-reserves planning units that don't have a negative value when removed
void iterativeImprovement(int puno,int spno,struct spustuff pu[], struct sconnections connections[],
                          struct sspecies spec[],struct spu SM[],int R[], double cm,
                          struct scost *reserve,struct scost *change,double costthresh,double tpf1, double tpf2,
                          int clumptype,int irun,char *savename)
{
    int puvalid =0,i,j,ipu=0,imode,ichoice, iRowCounter, iRowLimit;
    struct iimp *iimparray;
    double debugfloat;
    char tempname2[100];
    FILE *ttfp,*Rfp;
    char *writename;

    appendTraceFile("iterativeImprovement start\n");

    // counting pu's we need to test
    for (i=0;i<puno;i++)
    {
        if ((R[i] < 2) && (pu[i].status < 2))
            puvalid++;
    }

    appendTraceFile("iterativeImprovement puvalid %i\n",puvalid);

    if (fnames.saveitimptrace)
    {
        sprintf(tempname2,"%s_itimp_objective%05i.csv",savename,irun%10000);
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen(tempname2) + 2, sizeof(char));
        strcpy(writename,tempname2);
        if ((ttfp = fopen(writename,"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename);
        free(writename);
        fprintf(ttfp,"improvement,total,pus,cost,connection,penalty,change total\n");

        sprintf(tempname2,"%s_itimp_zones%05i.csv",savename,irun%10000);
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen(tempname2) + 2, sizeof(char));
        strcat(writename,tempname2);
        if ((Rfp = fopen(writename,"w"))==NULL)
            displayErrorMessage("cannot create threshold trace file %s\n",writename);
        free(writename);
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
        iimparray = (struct iimp *) calloc(puvalid,sizeof(struct iimp));

        for (i=0;i<puno;i++)
        {
            if ((R[i] < 2) && (pu[i].status < 2))
            {
                iimparray[ipu].puindex = i;
                iimparray[ipu].randomfloat = returnRandomFloat();
                ipu++;
            }
        }

        appendTraceFile("iterativeImprovement after array init\n");

        // sort the iimp array by the randomindex field
        heapSortIterativeImprovement(iimparray,puvalid);

        appendTraceFile("iterativeImprovement after heapSortIterativeImprovement\n");

        /***** Doing the improvements ****/
        for (i=0;i<puvalid;i++)
        {
            ichoice = iimparray[i].puindex;

            if ((R[ichoice] < 2) && (pu[ichoice].status < 2))
            {
               imode = R[ichoice] == 1 ? -1 : 1;
               computeChangeScore(-1,ichoice,spno,puno,pu,connections,spec,SM,R,cm,imode,change,reserve,
                                  costthresh,tpf1,tpf2,1,clumptype);
               if (change->total < 0)
               {
                  displayProgress2("It Imp has changed %i with change value %lf \n",ichoice,change->total);
                  doChange(ichoice,puno,R,reserve,*change,pu,SM,spec,connections,imode,clumptype);
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
                          ,i,reserve->total
                          ,reserve->pus,reserve->cost,reserve->connection,reserve->penalty,
                          change->total); // i,costthresh,pus,cost,connection,penalty

                  for (j = 0;j<puno;j++)
                      fprintf(Rfp,",%i",R[j]);
               }
            }
        }// no untested PUs left

        free(iimparray);
    }

    if (fnames.saveitimptrace)
    {
        fclose(ttfp);
        fclose(Rfp);
    }

    appendTraceFile("iterativeImprovement end\n");
} // iterativeImprovement

// ran1() from numerical recipes: produces a random number between 0 and 1
#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1 + (IM - 1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long RandomIY;
long RandomIV[NTAB];

// random random floating point number
// ran1() from numerical recipes: produces a random number between 0 and 1
float returnRandomFloat(void)
{
    int j;
    long k;
    float temp;

    if(RandSeed1 <= 0 || !RandomIY) // Initialize
    {
        RandSeed1 = -RandSeed1;
        for(j = NTAB+7; j >= 0; j--)
        {
            k = RandSeed1/IQ;
            RandSeed1 = IA * (RandSeed1 - k * IQ) - IR * k;
            if (RandSeed1 < 0)
                RandSeed1 += IM;
            if (j < NTAB)
                RandomIV[j] = RandSeed1;
        }
        RandomIY = RandomIV[0];
    }
    k=RandSeed1/IQ; // The stuff we do on calls after the first
    RandSeed1 = IA * (RandSeed1 - k * IQ) - IR * k;
    if(RandSeed1 < 0)
    {
        RandSeed1 += IM;
    }
    j = RandomIY/NDIV;
    RandomIY=RandomIV[j];
    RandomIV[j] = RandSeed1;
    if ((temp=AM*RandomIY) > RNMX)
    {
        return(RNMX);
    } else {
        return(temp);
    }
}

// initialise random seed
void initialiseRandomSeed(int iSeed) {
    if (iSeed>0)
        RandSeed1 = iSeed;
    else
        RandSeed1 = (long int)time(NULL);
    if(RandSeed1 > 0) RandSeed1 = -RandSeed1;
}

// return random number between 0 and parameter n-1
int returnRandom (int num)
{
    long    temp;

    if (num == 0)
        return(0);

    temp = (int)(returnRandomFloat() * num);
    if (temp == num) return(0);
    else return((int)temp);
}

// handle command line parameters for the marxan executable
void handleOptions(int argc,char *argv[],char sInputFileName[])
{
    int i;

    if (argc>4)
    {  // if more than one commandline argument then exit
        displayUsage(argv[0]);
        exit(1);
    }

    for (i=1;i<argc;i++)
    { // Deal with all arguments
        if (argv[i][0] == '/' || argv[i][0] == '-')
        {
            switch(argv[i][1])
            {
                case 'C':
                case 'c':
                case 'S':
                case 's':
                    marxanIsSlave = 1;
                break;
            default:
                fprintf(stderr,"unknown option %s\n",argv[i]);
                break;
            }
        } else {
            strcpy(sInputFileName,argv[i]); /* If not a -option then must be input.dat name */
        }
    }
}

// c main function for marxan: gets called when the marxan executable is run
int main(int argc,char *argv[])
{
    char sInputFileName[100];

    strcpy(sApplicationPathName,argv[0]);

    // set default input name
    strcpy(sInputFileName,"input.dat");

    if (argc == 1)
    {
        // If no arguments then assume the default file name
        strcpy(sInputFileName,"input.dat");
    }
    else     // handle the program options
        handleOptions(argc,argv,sInputFileName);

    if (executeMarxan(sInputFileName)) // Calls the main annealing unit
    {
        if (marxanIsSlave == 1)
        {
            slaveExit();
        }

        return 1;
    }  // Abnormal Exit
    if (marxanIsSlave == 1)
    {
        slaveExit();
    }

    return 0;
}

