// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

// C file for Marxan

//    Marxan coded by Ian Ball, modified by Matthew Watts
//    Written by Ian Ball and Hugh Possingham

//    ian.ball@aad.gov.au
//    hpossingham@zen.uq.edu.au
//    m.watts@uq.edu.au


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
int marxanisslave = 0;
char sVersionString[80] = "Marxan v 3.0.0";
// version 2.3 introduces multiple connectivity files and their associated weighting file
// version 2.4.3 introduces 1D and 2D probability
char sIanBallEmail[100] = "ian.ball@aad.gov.au";
char sHughPossinghamEmail[100] = "h.possingham@uq.edu.au";
char sMattWattsEmail[100] = "m.watts@uq.edu.au";
char sMarxanWebSite[100] = "http://www.uq.edu.au/marxan";
//char sTraceFileName[80] = "TraceFile_MarOpt.txt";
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

#include "input.c"
#include "output.c"

void ExecuteRunLoop(long int repeats,int puno,int spno,double cm,int aggexist,double prop,int clumptype,double misslevel,
                    char savename[],double costthresh,double tpf1,double tpf2,int heurotype,int runopts,
                    int itimptype,int *iBestRun,int sumsoln[],int marxanisslave)
{
    long int irun;
    int **R2D;
    char tempname2[500];
    int ipu, tid, nthreads;
    double rBestScore;
    
    // NOTE: we have to do something tricky with sumsoln
    // Remember the R for each run in a 2d matrix.
    // Figure out which is best run after threads terminate.
    // Add them to get sumsoln after threads terminate.
    R2D = (int **) calloc(repeats,sizeof(int *));
    for (irun = 1;irun <= repeats;irun++)
        R2D[irun-1] = (int *) calloc(puno,sizeof(int));

    // *******  The larger repetition loop ********
    for (irun = 1;irun <= repeats;irun++)
    {
        appendTraceFile("start run loop run %ld\n",irun);

        displayProgress1("\n");
        displayProgress("Run %ld ",irun);

        if (runoptions.ThermalAnnealingOn)
        {
           if (anneal.type >= 2)
           {
              if (anneal.type == 2 )
              {
                 appendTraceFile("before ConnollyInit run %i\n",irun);

                 ConnollyInit(puno,spno,pu,connections,spec,SM,cm,&anneal,aggexist,R2D[irun-1],prop,clumptype,irun);

                 appendTraceFile("after ConnollyInit run %i\n",irun);
              }
              if (anneal.type == 3)
              {
                 appendTraceFile("before AdaptiveInit run %i\n",irun);

                 AdaptiveInit(puno,spno,prop,R2D[irun-1],pu,connections,SM,cm,spec,aggexist,&anneal,clumptype);

                 appendTraceFile("after AdaptiveInit run %i\n",irun);
              }

              displayProgress1("  Using Calculated Tinit = %.4f Tcool = %.8f \n",
                               anneal.Tinit,anneal.Tcool);
           }  // Using Precalced Temperature Settings
           if (anneal.type == 3)
           {
              // Call annealing init here
           }  // using adaptive annealing type 2 here

           anneal.temp = anneal.Tinit;
        }  // Annealing Setup

        appendTraceFile("before ReserveCost run %i\n",irun);

        displayProgress1("  Creating the initial reserve \n");
        
        InitReserve(puno,prop,R2D[irun-1]);  // Create Initial Reserve
       
        AddReserve(puno,pu,R2D[irun-1]);

        if (aggexist)
           ClearClumps(spno,spec,pu,SM);

        ReserveCost(puno,spno,R2D[irun-1],pu,connections,SM,cm,spec,aggexist,&reserve,clumptype);

        appendTraceFile("after ReserveCost run %i\n",irun);

        if (verbosity > 1)
        {
           displayProgress1("\n  Init:");
           PrintResVal(puno,spno,R2D[irun-1],reserve,spec,misslevel);
        }
        if (verbosity > 5)
        {
           displayTimePassed();
        }

        if (runoptions.ThermalAnnealingOn)
        {
           appendTraceFile("before ThermalAnnealing run %i\n",irun);

           ThermalAnnealing(spno,puno,connections,R2D[irun-1],cm,spec,pu,SM,&change,&reserve,
                            repeats,irun,savename,misslevel,
                            aggexist,costthresh,tpf1,tpf2,clumptype);

           appendTraceFile("after ThermalAnnealing run %i\n",irun);
        }  // End of Thermal Annealing On

        if (runoptions.QuantumAnnealingOn)
        {
           appendTraceFile("before QuantumAnnealing run %i\n",irun);

           QuantumAnnealing(spno,puno,connections,R2D[irun-1],cm,spec,pu,SM,&change,&reserve,
                            repeats,irun,savename,misslevel,
                            aggexist,costthresh,tpf1,tpf2,clumptype);

           appendTraceFile("after QuantumAnnealing run %i\n",irun);
        }  // End of Quantum Annealing On

        if (runoptions.HeuristicOn)
        {
           appendTraceFile("before Heuristics run %i\n",irun);

           Heuristics(spno,puno,pu,connections,R2D[irun-1],cm,spec,SM,&reserve,
                      costthresh,tpf1,tpf2,heurotype,clumptype);

           if (verbosity > 1 && (runopts == 2 || runopts == 5))
           {
              displayProgress1("  Heuristic:");
              PrintResVal(puno,spno,R2D[irun-1],reserve,spec,misslevel);
           }

           appendTraceFile("after Heuristics run %i\n",irun);
        }    // Activate Greedy

        if (runoptions.ItImpOn)
        {
           appendTraceFile("before IterativeImprovement run %i\n",irun);

           IterativeImprovement(puno,spno,pu,connections,spec,SM,R2D[irun-1],cm,
                                &reserve,&change,costthresh,tpf1,tpf2,clumptype,irun,savename);

           if (itimptype == 3)
              IterativeImprovement(puno,spno,pu,connections,spec,SM,R2D[irun-1],cm,
                                   &reserve,&change,costthresh,tpf1,tpf2,clumptype,irun,savename);

           appendTraceFile("after IterativeImprovement run %i\n",irun);

           if (aggexist)
              ClearClumps(spno,spec,pu,SM);

           if (verbosity > 1)
           {
              ReserveCost(puno,spno,R2D[irun-1],pu,connections,SM,cm,spec,aggexist,&reserve,clumptype);
              displayProgress1("  Iterative Improvement:");
              PrintResVal(puno,spno,R2D[irun-1],reserve,spec,misslevel);
           }

        } // Activate Iterative Improvement

        appendTraceFile("before file output run %i\n",irun);

        if (fnames.saverun)
        {
           if (fnames.saverun == 3)
              sprintf(tempname2,"%s_r%05li.csv",savename,irun%10000);
           else
           if (fnames.saverun == 2)
              sprintf(tempname2,"%s_r%05li.txt",savename,irun%10000);
           else
               sprintf(tempname2,"%s_r%05li.dat",savename,irun%10000);
           writeSolution(puno,R2D[irun-1],pu,tempname2,fnames.saverun,fnames);
        }

        if (fnames.savespecies && fnames.saverun)
        {
           if (fnames.savespecies == 3)
              sprintf(tempname2,"%s_mv%05li.csv",savename,irun%10000);
           else
           if (fnames.savespecies == 2)
              sprintf(tempname2,"%s_mv%05li.txt",savename,irun%10000);
           else
               sprintf(tempname2,"%s_mv%05li.dat",savename,irun%10000);

           writeSpecies(spno,spec,tempname2,fnames.savespecies,misslevel);
        }

        if (fnames.savesum)
        {
           if (fnames.savesum==3)
              sprintf(tempname2,"%s_sum.csv",savename);
           else
           if (fnames.savesum==2)
              sprintf(tempname2,"%s_sum.txt",savename);
           else
               sprintf(tempname2,"%s_sum.dat",savename);
           writeSummary(puno,spno,R2D[irun-1],spec,reserve,irun,tempname2,misslevel,fnames.savesum);
        }
        
        // compute and store objective function score for this reserve system
        ReserveCost(puno,spno,R2D[irun-1],pu,connections,SM,cm,spec,aggexist,&change,clumptype);
        rBestScore = change.total;

        if (fnames.savesolutionsmatrix)
        {
           appendTraceFile("before appendSolutionsMatrix savename %s\n",savename);

           if (fnames.savesolutionsmatrix==3)
              sprintf(tempname2,"%s_solutionsmatrix.csv",savename);
           else
           if (fnames.savesolutionsmatrix==2)
              sprintf(tempname2,"%s_solutionsmatrix.txt",savename);
           else
               sprintf(tempname2,"%s_solutionsmatrix.dat",savename);

           appendSolutionsMatrix(irun,puno,R2D[irun-1],tempname2,fnames.savesolutionsmatrix,fnames.solutionsmatrixheaders);

           appendTraceFile("after appendSolutionsMatrix savename %s\n",savename);
        }
        
        if (aggexist)
           ClearClumps(spno,spec,pu,SM);

        appendTraceFile("after file output run %i\n",irun);
        appendTraceFile("end run %i\n",irun);

        if (marxanisslave == 1)
           writeSlaveSyncFileRun(irun);

        if (verbosity > 1)
           displayTimePassed();

    } // *** the repeats  ****
    // *******  The larger repetition loop ********
}

int Marxan(char sInputFileName[])
{
    int iSparseMatrixFileLength = 0, iSparseMatrixFileLength_sporder = 0; //, iProbSparseMatrixFileLength = 0;
    long int repeats;
    int puno,spno,gspno;
    struct sgenspec *gspec;
    double cm,prop;
    int runopts,heurotype,clumptype,itimptype;
    char savename[500],tempname2[500];
    double misslevel;
    int iseed,seedinit;
    int aggexist=0,sepexist=0;
    int *R_CalcPenalties; //,*bestrun,
    int *sumsoln, iBestRun = 1;
    double costthresh,tpf1,tpf2, rBestScore;
    long int itemp;
    int isp;
    #ifdef DEBUGTRACEFILE
    char debugbuffer[200];
    #endif

    // Handle Error driven termination
    if (setjmp(jmpbuf))
       return 1;

    displayStartupMessage();

    SetOptions(&cm,&prop,&anneal,
               &iseed,&repeats,savename,&fnames,sInputFileName,
               &runopts,&misslevel,&heurotype,&clumptype,&itimptype,&verbosity,
               &costthresh,&tpf1,&tpf2);

    SetRunOptions(runopts,&runoptions);

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

    InitRandSeed(iseed);
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
    PrepareBinarySearchArrays(puno,spno,pu,spec,&PULookup,&SPLookup);

    appendTraceFile("after build search arrays\n");

    bestyet = (int *) calloc(puno,sizeof(int));

    if (fnames.savesumsoln)
       sumsoln = (int *) calloc(puno,sizeof(int));

    connections = (typeconnection *) calloc(puno,sizeof(typeconnection));

    displayProgress3("    Reading in the Connection file :\n");
    itemp = 0;
    if (strcmp("NULL",fnames.connectionname) != 0)
    {
       appendTraceFile("before readConnections\n");

       if (strcmp("NULL",fnames.connectionfilesname) != 0)
          PrepareWeightedConnectivityFile(fnames);

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
          sprintf(tempname2,"%s_richness.csv",savename);
       else
       if (fnames.saverichness==2)
          sprintf(tempname2,"%s_richness.txt",savename);
       else
           sprintf(tempname2,"%s_richness.dat",savename);

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
        SetBlockDefs(gspno,spno,puno,gspec,spec,pu,SM);
    }

    SetDefs(spno,spec);

    appendTraceFile("after process block definitions\n");

    appendTraceFile("before computeTotalAreas\n");
    computeTotalAreas(puno,spno,pu,spec,SM);
    appendTraceFile("after computeTotalAreas\n");

    if (fnames.savetotalareas)
    {
       if (fnames.savetotalareas==3)
          sprintf(tempname2,"%s_totalareas.csv",savename);
       else
       if (fnames.savetotalareas==2)
          sprintf(tempname2,"%s_totalareas.txt",savename);
       else
           sprintf(tempname2,"%s_totalareas.dat",savename);

       writeTotalAreas(puno,spno,pu,spec,SM,tempname2,fnames.savepenalty);
    }

    if (fSpecPROPLoaded > 0)
    {
       appendTraceFile("before ApplySpecProp\n");

       // species have prop value specified
       ApplySpecProp(spno,spec,puno,pu,SM);

       appendTraceFile("after ApplySpecProp\n");
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
       MapUserPenalties(spec,spno);
	}
	else
	{
        // we are computing penalties
        if (strcmp("NULL",fnames.matrixspordername) == 0)
        {
           appendTraceFile("before CalcPenalties\n");

           // we don't have sporder matrix available, so use slow CalcPenalties method
           itemp = CalcPenalties(puno,spno,pu,spec,connections,SM,R_CalcPenalties,aggexist,cm,clumptype);

           appendTraceFile("after CalcPenalties\n");
        }
        else
        {
            // we have sporder matrix available, so use optimised CalcPenalties method
            if (iOptimisationCalcPenalties == 1)
            {
               appendTraceFile("before CalcPenaltiesOptimise\n");

               itemp = CalcPenaltiesOptimise(puno,spno,pu,spec,connections,SM,SMsporder,R_CalcPenalties,aggexist,cm,clumptype);

               appendTraceFile("after CalcPenaltiesOptimise\n");
            }
            else
            {
                appendTraceFile("before CalcPenalties\n");

                // we have optimise calc penalties switched off, so use slow CalcPenalties method
                itemp = CalcPenalties(puno,spno,pu,spec,connections,SM,R_CalcPenalties,aggexist,cm,clumptype);

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
          sprintf(tempname2,"%s_penalty.csv",savename);
       else
       if (fnames.savepenalty==2)
          sprintf(tempname2,"%s_penalty.txt",savename);
       else
           sprintf(tempname2,"%s_penalty.dat",savename);

       writePenalty(spno,spec,tempname2,fnames.savepenalty);

       if (fnames.savepenalty==3)
          sprintf(tempname2,"%s_penalty_planning_units.csv",savename);
       else
       if (fnames.savepenalty==2)
          sprintf(tempname2,"%s_penalty_planning_units.txt",savename);
       else
           sprintf(tempname2,"%s_penalty_planning_units.dat",savename);

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
       if (runoptions.ThermalAnnealingOn == 0)
          if (runoptions.QuantumAnnealingOn == 0)
             if (runoptions.ItImpOn == 0)
             {
                appendTraceFile("end final file output\n");
                appendTraceFile("\nMarxan end execution\n");
                displayShutdownMessage();

                if (marxanisslave == 1)
                   SlaveExit();

                if (aggexist)
                   ClearClumps(spno,spec,pu,SM);  // Remove these pointers for cleanliness sake

                if (fnames.savelog)
                   createLogFile(0,NULL);  /* tidy up files */

                exit(1);
			 }

    if (fnames.savesolutionsmatrix)
    {
       #ifdef DEBUG_CLUSTERANALYSIS
       appendTraceFile("before sprintf savename %s\n",savename);
       #endif

       if (fnames.savesolutionsmatrix==3)
          sprintf(tempname2,"%s_solutionsmatrix.csv",savename);
       else
       if (fnames.savesolutionsmatrix==2)
          sprintf(tempname2,"%s_solutionsmatrix.txt",savename);
       else
           sprintf(tempname2,"%s_solutionsmatrix.dat",savename);

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
  
    ExecuteRunLoop(repeats,puno,spno,cm,aggexist,prop,clumptype,misslevel,
                   savename,costthresh,tpf1,tpf2,heurotype,runopts,
                   itimptype,&iBestRun,sumsoln,marxanisslave);

    appendTraceFile("before final file output\n");

    if (fnames.savebest)
    {
       if (fnames.savebest == 3)
          sprintf(tempname2,"%s_best.csv",savename);
       else
       if (fnames.savebest == 2)
          sprintf(tempname2,"%s_best.txt",savename);
       else
           sprintf(tempname2,"%s_best.dat",savename);

       writeSolution(puno,bestyet,pu,tempname2,fnames.savebest,fnames);

       appendTraceFile("Best solution is run %i\n",iBestRun);
       displayProgress1("\nBest solution is run %i\n",iBestRun);
    }

    if (fnames.savespecies && fnames.savebest)
    {
       if (fnames.savespecies ==3)
          sprintf(tempname2,"%s_mvbest.csv",savename);
       else
       if (fnames.savespecies ==2)
          sprintf(tempname2,"%s_mvbest.txt",savename);
       else
           sprintf(tempname2,"%s_mvbest.dat",savename);

       writeSpecies(spno,spec,tempname2,fnames.savespecies,misslevel);
    }

    if (fnames.savesumsoln)
    {
       if (fnames.savesumsoln == 3)
          sprintf(tempname2,"%s_ssoln.csv",savename);
       else
       if (fnames.savesumsoln == 2)
          sprintf(tempname2,"%s_ssoln.txt",savename);
       else
           sprintf(tempname2,"%s_ssoln.dat",savename);

       writeSumSoln(puno,sumsoln,pu,tempname2,fnames.savesumsoln);
    }

    if (fnames.savesolutionsmatrix)
       if (fnames.rexecutescript)
       {
          if (fnames.savesolutionsmatrix==3)
             sprintf(tempname2,"%s_solutionsmatrix.csv",savename);
          else
          if (fnames.savesolutionsmatrix==2)
             sprintf(tempname2,"%s_solutionsmatrix.txt",savename);
          else
              sprintf(tempname2,"%s_solutionsmatrix.dat",savename);

          if (marxanisslave == 1)
             SlaveExit();
       }

    if (aggexist)
       ClearClumps(spno,spec,pu,SM);  // Remove these pointers for cleanliness sake

    displayShutdownMessage();

    if (fnames.savelog)
       createLogFile(0,NULL);  /* tidy up files */

    appendTraceFile("end final file output\n");
    appendTraceFile("\nMarxan end execution\n");

    return 0;
}   // Marxan

// returns the 0-base index of a species at a planning unit, if the species doesn't occur here, returns -1
int rtnIdxSpecAtPu(struct spustuff PU[], struct spu SM[], int iPUIndex, int iSpecIndex)
{
    int i;

    if (PU[iPUIndex].richness > 0)
    {
       for (i=0;i<PU[iPUIndex].richness;i++)
           if (SM[PU[iPUIndex].offset + i].spindex == iSpecIndex)
              return (PU[iPUIndex].offset + i);
    }

    return -1;
}

// returns the clump number of a species at a planning unit, if the species doesn't occur here, returns 0
int rtnClumpSpecAtPu(struct spustuff PU[], struct spu SM[], int iPUIndex, int iSpecIndex)
{
    int i;

    if (PU[iPUIndex].richness > 0)
       for (i=0;i<PU[iPUIndex].richness;i++)
           if (SM[PU[iPUIndex].offset + i].spindex == iSpecIndex)
              return SM[PU[iPUIndex].offset + i].clump;

    return 0;
}

// sets the clump number of a species at a planning unit
void setClumpSpecAtPu(struct spustuff PU[], struct spu SM[], int iPUIndex, int iSpecIndex, int iSetClump)
{
    int i;

    if (PU[iPUIndex].richness > 0)
       for (i=0;i<PU[iPUIndex].richness;i++)
           if (SM[PU[iPUIndex].offset + i].spindex == iSpecIndex)
              SM[PU[iPUIndex].offset + i].clump = iSetClump;
}

// returns the amount of a species at a planning unit, if the species doesn't occur here, returns 0
double rtnAmountSpecAtPu(struct spustuff PU[], struct spu SM[], int iPUIndex, int iSpecIndex)
{
    int i;

    if (PU[iPUIndex].richness > 0)
       for (i=0;i<PU[iPUIndex].richness;i++)
           if (SM[PU[iPUIndex].offset + i].spindex == iSpecIndex)
              return SM[PU[iPUIndex].offset + i].amount;

    return 0;
}

void SetBlockDefs(int gspno,int spno,int puno,struct sgenspec gspec[], struct sspecies spec[],
                  struct spustuff PU[], struct spu SM[])
{
  int igsp,isp,ipu;
  double totalamount;

  for (igsp=0;igsp<gspno;igsp++)
  {
      if (gspec[igsp].prop > 0) // deal with percentage in a different way
         for (isp=0;isp<spno;isp++)
             if (spec[isp].type == gspec[igsp].type && spec[isp].target < 0)
             {
                for (ipu=0,totalamount =0;ipu<puno;ipu++)
                    totalamount += rtnAmountSpecAtPu(PU,SM,ipu,isp);
                spec[isp].target = totalamount * gspec[igsp].prop;
             } // Setting target with percentage
      if (gspec[igsp].target > 0)
         for (isp=0;isp<spno;isp++)
             if (spec[isp].type == gspec[igsp].type && spec[isp].target < 0)
                spec[isp].target = gspec[igsp].target;
      if (gspec[igsp].target2 > 0)
         for (isp=0;isp<spno;isp++)
         {
             if (spec[isp].type == gspec[igsp].type && spec[isp].target2 < 0)
             {
                spec[isp].target2 = gspec[igsp].target2;
             }
         }
      if (gspec[igsp].targetocc > 0)
         for (isp=0;isp<spno;isp++)
             if (spec[isp].type == gspec[igsp].type && spec[isp].targetocc < 0)
                spec[isp].targetocc = gspec[igsp].targetocc;
      if (gspec[igsp].sepnum > 0)
         for (isp=0;isp<spno;isp++)
             if (spec[isp].type == gspec[igsp].type && spec[isp].sepnum < 0)
                spec[isp].sepnum = gspec[igsp].sepnum;
      if (gspec[igsp].spf > 0)
         for (isp=0;isp<spno;isp++)
             if (spec[isp].type == gspec[igsp].type && spec[isp].spf < 0)
                spec[isp].spf = gspec[igsp].spf;
      if (gspec[igsp].sepdistance > 0)
         for (isp=0;isp<spno;isp++)
             if (spec[isp].type == gspec[igsp].type && spec[isp].sepdistance < 0)
                spec[isp].sepdistance = gspec[igsp].sepdistance;
      // Percentage is not dealt with here yet. To do this I need to identify
      // target species then determine their total abundance then set target
      // according to percentage
    } // Loop through each setBlockDef

} // Set Block Defs

// ******** Set Defaults *******
// If '-1' values haven't been set yet then this one will do it
void SetDefs(int spno, struct sspecies spec[])
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
} // Set Defs

// *** Set run options. Takes an integer runopts value and returns flags ***
void SetRunOptions(int runopts, struct srunoptions *runoptions)
{
     if (runopts < 0)
        return; // runopts < 0 indicates that these are set in some other way

     switch (runopts)
     {
            case 0: (*runoptions).CalcPenaltiesOn = 1;
                    (*runoptions).ThermalAnnealingOn = 1;
                    (*runoptions).QuantumAnnealingOn = 0;
                    (*runoptions).HeuristicOn = 1;
                    (*runoptions).ItImpOn = 0;
                    break;
            case 1: (*runoptions).CalcPenaltiesOn = 1;
                    (*runoptions).ThermalAnnealingOn = 1;
                    (*runoptions).QuantumAnnealingOn = 0;
                    (*runoptions).HeuristicOn = 0;
                    (*runoptions).ItImpOn = 1;
                    break;
            case 2: (*runoptions).CalcPenaltiesOn = 1;
                    (*runoptions).ThermalAnnealingOn = 1;
                    (*runoptions).QuantumAnnealingOn = 0;
                    (*runoptions).HeuristicOn = 1;
                    (*runoptions).ItImpOn = 1;
                    break;
            case 3: (*runoptions).CalcPenaltiesOn = 0;
                    (*runoptions).ThermalAnnealingOn = 0;
                    (*runoptions).QuantumAnnealingOn = 0;
                    (*runoptions).HeuristicOn = 1;
                    (*runoptions).ItImpOn = 0;
                    break;
            case 4: (*runoptions).CalcPenaltiesOn = 0;
                    (*runoptions).ThermalAnnealingOn = 0;
                    (*runoptions).QuantumAnnealingOn = 0;
                    (*runoptions).HeuristicOn = 0;
                    (*runoptions).ItImpOn = 1;
                    break;
            case 5: (*runoptions).CalcPenaltiesOn = 0;
                    (*runoptions).ThermalAnnealingOn = 0;
                    (*runoptions).QuantumAnnealingOn = 0;
                    (*runoptions).HeuristicOn = 1;
                    (*runoptions).ItImpOn = 1;
                    break;
            case 6: (*runoptions).CalcPenaltiesOn = 1;
                    (*runoptions).ThermalAnnealingOn = 1;
                    (*runoptions).QuantumAnnealingOn = 0;
                    (*runoptions).HeuristicOn = 0;
                    (*runoptions).ItImpOn = 0;
                    break;
            case 7: (*runoptions).CalcPenaltiesOn = 1;
                    (*runoptions).ThermalAnnealingOn = 0;
                    (*runoptions).QuantumAnnealingOn = 0;
                    (*runoptions).HeuristicOn = 0;
                    (*runoptions).ItImpOn = 0;
                    break;
            case 8: (*runoptions).CalcPenaltiesOn = 0;
                    (*runoptions).ThermalAnnealingOn = 1;
                    (*runoptions).QuantumAnnealingOn = 0;
                    (*runoptions).HeuristicOn = 0;
                    (*runoptions).ItImpOn = 0;
                    break;
            case 9: (*runoptions).CalcPenaltiesOn = 0;
                    (*runoptions).ThermalAnnealingOn = 1;
                    (*runoptions).QuantumAnnealingOn = 0;
                    (*runoptions).HeuristicOn = 0;
                    (*runoptions).ItImpOn = 1;
                    break;
            case 10: (*runoptions).CalcPenaltiesOn = 0;
                     (*runoptions).ThermalAnnealingOn = 0;
                     (*runoptions).QuantumAnnealingOn = 0;
                     (*runoptions).HeuristicOn = 0;
                     (*runoptions).ItImpOn = 1;
                     break;
            case 11: (*runoptions).CalcPenaltiesOn = 1;
                     (*runoptions).ThermalAnnealingOn = 0;
                     (*runoptions).QuantumAnnealingOn = 1;
                     (*runoptions).HeuristicOn = 0;
                     (*runoptions).ItImpOn = 0;
                     break;
            case 12: (*runoptions).CalcPenaltiesOn = 1;
                     (*runoptions).ThermalAnnealingOn = 0;
                     (*runoptions).QuantumAnnealingOn = 1;
                     (*runoptions).HeuristicOn = 0;
                     (*runoptions).ItImpOn = 1;
                     break;
            case 13: (*runoptions).CalcPenaltiesOn = 0;
                     (*runoptions).ThermalAnnealingOn = 0;
                     (*runoptions).QuantumAnnealingOn = 1;
                     (*runoptions).HeuristicOn = 0;
                     (*runoptions).ItImpOn = 0;
                     break;
            case 14: (*runoptions).CalcPenaltiesOn = 0;
                     (*runoptions).ThermalAnnealingOn = 0;
                     (*runoptions).QuantumAnnealingOn = 1;
                     (*runoptions).HeuristicOn = 0;
                     (*runoptions).ItImpOn = 1;
                     break;
            default: (*runoptions).CalcPenaltiesOn = 0;
                     (*runoptions).ThermalAnnealingOn = 0;
                     (*runoptions).QuantumAnnealingOn = 0;
                     (*runoptions).HeuristicOn = 0;
                     (*runoptions).ItImpOn = 0;
                     break;
     }
} // Set Run Options

// ******** Add Reserve to current system *******
void AddReserve(int puno,struct spustuff pu[],int *R)
{
     int i;
     for (i=0;i<puno;i++)
     {
         if (pu[i].status)
            R[i] = pu[i].status; // Change status for status 1, or higher
     }
}    // ******* Add Reserve **********

// ******* Calculate Initial Penalties *************
// This routine calculates the initial penalties or the penalty if you had no representation
// This has been modified in March 2000 to allow a combination of area targets and occurrence
//    targets
// If species has spatial requirements then CalcPenaltyType4 is invoked and will handle that
//    species
int CalcPenalties(int puno,int spno,struct spustuff pu[],struct sspecies spec[],
                  struct sconnections connections[],struct spu SM[],int PUtemp[],int aggexist,double cm,int clumptype)
{
    int i,j,ibest,imaxtarget,itargetocc;
    //int *PUtemp;
    double ftarget,fbest,fbestrat,fcost,ftemp, rAmount;
    int badspecies = 0,goodspecies = 0;

    //PUtemp = (int *) calloc(puno,sizeof(int));

    AddReserve(puno,pu,PUtemp); // Adds existing reserve to PUtemp

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
               ftarget += rtnAmountSpecAtPu(pu,SM,j,i);
               itargetocc++;
               spec[i].penalty += cost(j,pu,connections,cm);
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
          {
              rAmount = rtnAmountSpecAtPu(pu,SM,j,i);
              if (PUtemp[j] == 0 && rAmount>0)
              {
                 fcost = cost(j,pu,connections,cm);
                 if (fcost == 0)
                    fcost = delta;
                 if (rAmount >= spec[i].target - ftarget && (imaxtarget == 0 || (imaxtarget == 1 && fcost < fbest)))
                 {
                    imaxtarget = 1;
                    ibest = j;
                    fbest = fcost;
                 } // can I meet the target cheaply?
                 else
                 if (fbestrat < rAmount/fcost)
                 {
                    fbest = fcost;
                    fbestrat = rAmount/fbest;
                    ibest = j;
                 }  // finding the cheapest planning unit

                 #ifdef DEBUGCALCPENALTIES
                 appendTraceFile("CalcPenalties species %i puid %i cost %g\n",spec[i].name,pu[j].id,fcost);
                 #endif
              }  // Making sure only checking planning units not already used
          }  // trying to find best pu

          if (fbest > 0)
          {
             PUtemp[ibest] = 1;
             ftarget += rtnAmountSpecAtPu(pu,SM,ibest,i);
             itargetocc++;
             spec[i].penalty += fbest;

             #ifdef DEBUGCALCPENALTIES
             appendTraceFile("CalcPenalties species %i puid %i ftarget %g fbest %g\n",spec[i].name,pu[ibest].id,ftarget,fbest);
             #endif
          } // Add pu to target
        } while ((ftarget <spec[i].target|| itargetocc < spec[i].targetocc) && fbest > 0); // or no more pu left

        if (fbest == 0) // Could not meet target using all available PUs
        {
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
        }  // If not met target with all available PUs

        #ifdef DEBUGTRACEFILE
        appendTraceFile("CalcPenalties spname %i penalty %g target %g\n",spec[i].name,spec[i].penalty,spec[i].target);
        #endif

    }  // Penalty for each individual Species
    // Clear clumps in case I needed them for target4 species

    if (aggexist)
        ClearClumps(spno,spec,pu,SM);

    //free(PUtemp);
    //DebugFree(puno*sizeof(int));

    if (goodspecies)
       displayProgress1("%i species are already adequately represented.\n",goodspecies);
    return(badspecies);
}

int CalcPenaltiesOptimise(int puno,int spno,struct spustuff pu[],struct sspecies spec[],
                          struct sconnections connections[],struct spu SM[],struct spusporder SMsp[],
                          int PUtemp[],int aggexist,double cm,int clumptype)
{
    int i,j,ibest,imaxtarget,itargetocc,ism,ipu, iPUsToTest;
    //int *PUtemp;
    double ftarget,fbest,fbestrat,fcost,ftemp, rAmount, r_ibest_amount;
    int badspecies = 0,goodspecies = 0;

    //PUtemp = (int *) calloc(puno,sizeof(int));

    appendTraceFile("CalcPenaltiesOptimise start\n");

    AddReserve(puno,pu,PUtemp); // Adds existing reserve to PUtemp

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
           for (j=0;j<spec[i].richness;j++)  // traverse pu's containing this sp
           {
               ism = spec[i].offset + j;
               ipu = SMsp[ism].puindex;

                 if (PUtemp[ipu] < 2)
                  PUtemp[ipu] = 0;
               if (PUtemp[ipu] == 2)
               {
                  ftarget += SMsp[ism].amount;
                  itargetocc++;
                  spec[i].penalty += cost(ipu,pu,connections,cm);
               }
           } // reset PUtemp and also target

        // Already adequately represented on type 2 planning unit
        if (ftarget >= spec[i].target && itargetocc >= spec[i].targetocc)
        {
           goodspecies++;
           displayProgress2("Species %i (%s) has already met target %.2f\n",
                            spec[i].name,spec[i].sname,spec[i].target);
           continue;
        } // Target met in unremovable reserve

        //iPUsToTest = spec[i].richness;
        do
        {
          fbest =0; imaxtarget = 0; fbestrat = 0;
          if (spec[i].richness > 0)
             for (j=0;j<spec[i].richness;j++)  // traverse pu's containing this sp
             {
                 ism = spec[i].offset + j;
                 ipu = SMsp[ism].puindex;

                 rAmount = SMsp[ism].amount;
                  if (PUtemp[ipu] == 0)
                  {
                    fcost = cost(ipu,pu,connections,cm);
                    if (fcost == 0)
                       fcost = delta;
                    if (rAmount >= spec[i].target - ftarget && (imaxtarget == 0 ||
                             (imaxtarget == 1 && fcost < fbest)))
                    {
                       imaxtarget = 1;
                       ibest = ipu;
                       r_ibest_amount = rAmount;
                       fbest = fcost;
                    } // can I meet the target cheaply?
                    else if (fbestrat < rAmount/fcost)
                         {
                            fbest = fcost;
                            fbestrat = rAmount/fbest;
                            ibest = ipu;
                            r_ibest_amount = rAmount;
                         }  // finding the cheapest planning unit
                 }  // Making sure only checking planning units not already used
             }  // trying to find best pu

          if (fbest > 0)
          {
             PUtemp[ibest] = 1;
             ftarget += r_ibest_amount;
             itargetocc++;
             spec[i].penalty += fbest;
          } // Add pu to target

          //iPUsToTest--;

        } while (/*(iPUsToTest > 0)  ||*/ (fbest > 0) && (ftarget <spec[i].target|| itargetocc < spec[i].targetocc));
        // while there is some pu's with this species to test AND a best available pu was found AND targets are not met yet

        if (fbest == 0) // Could not meet target using all available PUs
        {
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
        }  /* If not met target with all available PUs */

        appendTraceFile("CalcPenaltiesOptimise spname %i penalty %g\n",spec[i].name,spec[i].penalty);
    }  // Penalty for each individual Species
    
    // Clear clumps in case I needed them for target4 species
    if (aggexist)
        ClearClumps(spno,spec,pu,SM);

    if (goodspecies)
         displayProgress1("%i species are already adequately represented.\n",goodspecies);

    appendTraceFile("CalcPenaltiesOptimise end\n");

    return(badspecies);
}// *** Calculate Initial Penalties ***

// ********** Connection Cost Type 1
// ** Total cost of all connections for PU independant of neighbour status
double ConnectionCost1(int ipu,struct spustuff pu[],struct sconnections connections[],double cm)
{
       double fcost;
       struct sneighbour *p;
       #ifdef DEBUG_CONNECTIONCOST
       char debugbuffer[200];
       #endif

       fcost = connections[ipu].fixedcost;
       for (p = connections[ipu].first;p;p=p->next)
           if (asymmetricconnectivity)
           {
              if (p->connectionorigon)
                 fcost += p->cost;
           }
           else
               fcost += p->cost;

       #ifdef DEBUG_CONNECTIONCOST
       sprintf(debugbuffer,"ConnectionCost1 ipu %i connection %g\n",ipu,fcost);
       appendTraceFile(debugbuffer);
       #endif

       return(fcost*cm);
}

// ********** Cost of Planning Unit & all connection costs of planning unit
// **** Used only when calculating penalties
double cost(int ipu,struct spustuff pu[],struct sconnections connections[],double cm)
{
       double fcost;

       fcost = pu[ipu].cost;
       fcost += ConnectionCost1(ipu,pu,connections,cm);

       return(fcost);
}

/************ Change in penalty for adding single PU ******/
double ChangePen(int ipu,int puno,struct sspecies spec[],struct spustuff pu[],struct spu SM[],
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
          appendTraceFile("ChangePen start\n");
       #endif

       *rShortfall = 0;

       if (pu[ipu].richness)
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
                 }
                 else
                 {
                    if (spec[isp].target)
                       newamount = NewPenalty(ipu,isp,spec,pu,SM,imode)/spec[isp].target;
                    if (spec[isp].targetocc)
                    {
                       tamount =  (double) (spec[isp].targetocc - spec[isp].occurrence - imode)    /
                                  (double) spec[isp].targetocc;
                       newamount += tamount<0?0:tamount;
                    }
                    if (spec[isp].target && spec[isp].targetocc)
                       newamount /= 2;
                    if (spec[isp].sepnum)
                       newamount +=SepPenalty(CountSeparation2(isp,ipu,NULL,puno,R,pu,SM,spec,imode),
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

       #ifdef ANNEALING_TEST
       if (ipu == (puno-1))
       {
          sprintf(debugbuffer,"ChangePen end fcost %g\n",fcost);
          appendTraceFile(debugbuffer);
       }
       #endif
       return (fcost);
}  /*** Change in penalty for adding or deleting one PU ****/

/************** Value of a Reserve System ********/
void ComputeP_AllPUsSelected_1D(char savename[],int puno,int spno,struct spustuff pu[],struct spu SM[],struct sspecies spec[])
{
     FILE *fp;
     int i,j,iHeavisideStepFunction,ism,isp;
     double *ExpectedAmount1D, *VarianceInExpectedAmount1D, *TA,
            rProbability, rRawP, rSumProbability = 0, rShortfallPenalty, rZScore;
     char debugbuffer[200];

     // create the output file
     fp = fopen(savename,"w");
     fprintf(fp,"SPID,amount held,ptarget1d,EA1D,VIEA1D,Z1D,rawP1D,heavisideSF1D,shortfallP1D,P1D\n");

     // init arrays
     ExpectedAmount1D = (double *) calloc(spno,sizeof(double));
     VarianceInExpectedAmount1D = (double *) calloc(spno,sizeof(double));
     TA = (double *) calloc(spno,sizeof(double));
     for (i=0;i<spno;i++)
     {
         ExpectedAmount1D[i] = 0;
         VarianceInExpectedAmount1D[i] = 0;
         TA[i] = 0;
     }

     // compute EA, VIEA, TA
     for (j=0;j<puno;j++)
         if (pu[j].richness)
            for (i=0;i<pu[j].richness;i++)
            {
                ism = pu[j].offset + i;
                isp = SM[ism].spindex;

                ExpectedAmount1D[isp] += SM[ism].amount * (1 - pu[j].prob);
                VarianceInExpectedAmount1D[isp] += SM[ism].amount * SM[ism].amount * pu[j].prob * (1 - pu[j].prob);
                TA[isp] += SM[ism].amount;
            }

     // compute probability for each feature
     for (i=0;i<spno;i++)
     {
         if (VarianceInExpectedAmount1D[i] > 0)
            rZScore = (spec[i].target - ExpectedAmount1D[i]) / sqrt(VarianceInExpectedAmount1D[i]);
         else
             rZScore = 4;

         if (rZScore >= 0)
            rRawP = probZUT(rZScore);
         else
             rRawP = 1 - probZUT(-1 * rZScore);

         if (spec[i].ptarget1d > rRawP)
            iHeavisideStepFunction = 1;
         else
             iHeavisideStepFunction = 0;

         if (spec[i].ptarget1d > 0)
            rShortfallPenalty = (spec[i].ptarget1d - rRawP) / spec[i].ptarget1d;
         else
             rShortfallPenalty = 0;

         rProbability = iHeavisideStepFunction * rShortfallPenalty;

         rSumProbability += rProbability;

         fprintf(fp,"%i,%f,%f,%f,%f,%f,%f,%i,%f,%f\n",
                    spec[i].name,TA[i],spec[i].ptarget1d,
                    ExpectedAmount1D[i],VarianceInExpectedAmount1D[i],rZScore,
                    rRawP,iHeavisideStepFunction,rShortfallPenalty,rProbability);
     }

     free(ExpectedAmount1D);
     free(VarianceInExpectedAmount1D);
     free(TA);
     fclose(fp);

     #ifdef DEBUGTRACEFILE
     sprintf(debugbuffer,"ComputeP_AllPUsSelected_1D SumP %f SumP * PW %f\n",rSumProbability,rSumProbability * rProbabilityWeighting);
     appendTraceFile(debugbuffer);
     #endif
}

void ComputeP_AllPUsSelected_2D(char savename[],int puno,int spno,struct spustuff pu[],struct spu SM[],struct sspecies spec[])
{
     FILE *fp;
     int i,j,iHeavisideStepFunction,ism,isp;
     double *ExpectedAmount2D, *VarianceInExpectedAmount2D, *TA,
            rProbability, rRawP, rSumProbability = 0, rShortfallPenalty, rZScore;
     char debugbuffer[200];

     // create the output file
     fp = fopen(savename,"w");
     fprintf(fp,"SPID,amount held,ptarget1d,EA2D,VIEA2D,Z2D,rawP2D,heavisideSF2D,shortfallP2D,P2D\n");

     // init arrays
     ExpectedAmount2D = (double *) calloc(spno,sizeof(double));
     VarianceInExpectedAmount2D = (double *) calloc(spno,sizeof(double));
     TA = (double *) calloc(spno,sizeof(double));
     for (i=0;i<spno;i++)
     {
         ExpectedAmount2D[i] = 0;
         VarianceInExpectedAmount2D[i] = 0;
         TA[i] = 0;
     }

     // compute EA, VIEA, TA
     for (j=0;j<puno;j++)
         if (pu[j].richness)
            for (i=0;i<pu[j].richness;i++)
            {
                ism = pu[j].offset + i;
                isp = SM[ism].spindex;

                ExpectedAmount2D[isp] += SM[ism].amount * SM[ism].prob;
                VarianceInExpectedAmount2D[isp] += SM[ism].amount * SM[ism].amount * SM[ism].prob * (1 - SM[ism].prob);
                TA[isp] += SM[ism].amount;
            }

     // compute probability for each feature
     for (i=0;i<spno;i++)
     {
         if (VarianceInExpectedAmount2D[i] > 0)
            rZScore = (spec[i].target - ExpectedAmount2D[i]) / sqrt(VarianceInExpectedAmount2D[i]);
         else
             rZScore = 4;

         if (rZScore >= 0)
            rRawP = probZUT(rZScore);
         else
             rRawP = 1 - probZUT(-1 * rZScore);

         if (spec[i].ptarget2d > rRawP)
            iHeavisideStepFunction = 1;
         else
             iHeavisideStepFunction = 0;

         if (spec[i].ptarget2d > 0)
            rShortfallPenalty = (spec[i].ptarget2d - rRawP) / spec[i].ptarget2d;
         else
             rShortfallPenalty = 0;

         rProbability = iHeavisideStepFunction * rShortfallPenalty;

         rSumProbability += rProbability;

         fprintf(fp,"%i,%f,%f,%f,%f,%f,%f,%i,%f,%f\n",
                    spec[i].name,TA[i],spec[i].ptarget2d,
                    ExpectedAmount2D[i],VarianceInExpectedAmount2D[i],rZScore,
                    rRawP,iHeavisideStepFunction,rShortfallPenalty,rProbability);
     }

     free(ExpectedAmount2D);
     free(VarianceInExpectedAmount2D);
     free(TA);
     fclose(fp);

     #ifdef DEBUGTRACEFILE
     sprintf(debugbuffer,"ComputeP_AllPUsSelected_2D SumP %f SumP * PW %f\n",rSumProbability,rSumProbability * rProbabilityWeighting);
     appendTraceFile(debugbuffer);
     #endif
}

/************** Value of a Reserve System ********/
void ReserveCost(int puno,int spno,int *R,struct spustuff pu[],
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

     #ifdef DEBUG_RESERVECOST
     //rConnectivityValue = reserve->connection;
     //appendTraceFile("total connectivity %lf\n",rConnectivityValue);
     #endif

     reserve->total = reserve->cost + reserve->connection + reserve->penalty + reserve->probability1D + reserve->probability2D;

     #ifdef DEBUG_PROB1D
     sprintf(sProbDebugCost,"probability1D %f\n",reserve->probability1D);
     appendTraceFile(sProbDebugCost);

     appendTraceFile("ReserveCost E\n");

     sprintf(sProbDebugCost,"%f",reserve->cost);
     strcpy(sProbDebugFileName,fnames.outputdir);
     strcat(sProbDebugFileName,"output_Prob1DDebug_");
     strcat(sProbDebugFileName,sProbDebugCost);
     strcat(sProbDebugFileName,".csv");
     writeProb1DDebugTable(spno,sProbDebugFileName,
                           ExpectedAmount1D,VarianceInExpectedAmount1D,spec);


     appendTraceFile("ReserveCost F\n");

     sprintf(sProbDebugCost,"%f",reserve->cost);
     strcpy(sProbDebugFileName,fnames.outputdir);
     strcat(sProbDebugFileName,"output_Prob1DDetailDebug_");
     strcat(sProbDebugFileName,sProbDebugCost);
     strcat(sProbDebugFileName,".csv");
     writeProb1DDetailDebugTable(sProbDebugFileName,puno,spno,pu,SM,R);

     appendTraceFile("ReserveCost G\n");
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

}  /**** value of a reserve ***/

/*********** Set the Initial Reserve System ************/
void InitReserve(int puno,double prop, int *R)
{
     int i;
     for (i=0;i<puno;i++)
         R[i] = rand1() < prop ? 1:0;
}


/******** Species Amounts ****************/
/******** puts in the effective amount for each species in reserve R */
void SpeciesAmounts(int spno,int puno,struct sspecies spec[],struct spustuff pu[],struct spu SM[],
                    int *R,int clumptype)
{
     int i, ism, isp,ipu;

     for (isp=0;isp<spno;isp++)
     {
         spec[isp].amount = 0;
         spec[isp].occurrence = 0;
         if (spec[isp].target2)
            SpeciesAmounts4(isp,spec,clumptype);

         spec[isp].expected1D = 0;
         spec[isp].expected2D = 0;
         spec[isp].variance1D = 0;
         spec[isp].variance2D = 0;
     }

     for (ipu=0;ipu<puno;ipu++)
         if (pu[ipu].richness)
            if (R[ipu] ==1 || R[ipu] == 2)
               for (i=0;i<pu[ipu].richness;i++)
               {
                   ism = pu[ipu].offset + i;
                   isp = SM[ism].spindex;
                   if (spec[isp].target2 == 0)
                   {
                      spec[isp].amount += SM[ism].amount;
                      spec[isp].occurrence++;
                   }
               }

     #ifdef ANNEALING_TEST
     for (isp=0;isp<spno;isp++)
     {
         appendTraceFile("SpeciesAmounts isp %i spec.amount %g\n",
                              isp,spec[isp].amount);
     }
     #endif

} /*** Species Amounts ***/

/** Threshold penalty used for Check Change ***/
/** only used when there is a cost threshold **/
double ThresholdPenalty(double tpf1,double tpf2,double timeprop)
{
    if (tpf2 < 0)
       return(tpf1);
    return(tpf1*exp(tpf2*timeprop));
} /* Threshold Penalty */

/*** Check Change ******************************/
void CheckChange(int iIteration,int ipu,int spno,int puno,struct spustuff pu[],struct sconnections connections[],
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
     //change->connection = imode * ConnectionCost2(ipu,connections,R,1,0,cm);
     if (threshtype ==1)
     {
        tchangeconnection = change->connection;
        tresconnection = reserve->connection;
        change->connection = 0;
        reserve->connection = 0;
     }

     change->penalty = ChangePen(ipu,puno,spec,pu,SM,R,connections,imode,clumptype,&change->shortfall);

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
                           ThresholdPenalty(tpf1,tpf2,timeprop);
        }
        else
        {
            if (change->cost + change->connection + reserve->cost + reserve->connection <= costthresh)
                threshpen = (reserve->cost + reserve->connection - costthresh) *
                            ThresholdPenalty(tpf1,tpf2,timeprop);
            else
                threshpen = (change->cost + change->connection) *
                            ThresholdPenalty(tpf1,tpf2,timeprop);
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

}  /*** Check Change ***/

void CheckQuantumChange(int spno,int puno,struct spustuff pu[],struct sconnections connections[],
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
     appendTraceFile("CheckQuantumChange start iFluctuationCount %i\n",iFluctuationCount);
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
         appendTraceFile("CheckQuantumChange ipu %i chosen %i imode %i\n",j,PUChosen[j],imode);
         #endif

         change->cost += pu[j].cost*imode;    /* Cost of this PU on it's own */
         change->connection += ConnectionCost2(j,connections,R,imode,1,cm);
         //change->connection = imode * ConnectionCost2(ipu,connections,R,1,0,cm);
         if (threshtype ==1)
         {
            tchangeconnection = change->connection;
            tresconnection = reserve->connection;
            change->connection = 0;
            reserve->connection = 0;
         }

         change->penalty += ChangePen(j,puno,spec,pu,SM,R,connections,imode,clumptype,&change->shortfall);

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
                               ThresholdPenalty(tpf1,tpf2,timeprop);
            }
            else
            {
                if (change->cost + change->connection + reserve->cost + reserve->connection <= costthresh)
                     threshpen = (reserve->cost + reserve->connection - costthresh) *
                                ThresholdPenalty(tpf1,tpf2,timeprop);
                else
                    threshpen = (change->cost + change->connection) *
                                ThresholdPenalty(tpf1,tpf2,timeprop);
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
     appendTraceFile("CheckQuantumChange end\n");
     #endif
}  // Check Quantum Change

// new Animal Penalty
double NewPenalty(int ipu,int isp,struct sspecies spec[],struct spustuff pu[],struct spu SM[],int imode)
{
       double newpen;

       newpen = spec[isp].target - spec[isp].amount - rtnAmountSpecAtPu(pu,SM,ipu,isp)*imode;

       if (newpen < 0)
          newpen = 0;

       return(newpen);
}  // Animal Penalty

int GoodChange(struct scost change,double temp)
{
    return (exp(-change.total/temp)> rand1()) ? 1 : 0;
}

int GoodQuantumChange(struct scost change,double rProbAcceptance)
{
    if (change.total <= 0)
       return 1;
    else
        return (rProbAcceptance > rand1()) ? 1 : 0;
}

void DoChange(int ipu,int puno,int *R,struct scost *reserve,struct scost change,
              struct spustuff pu[],struct spu SM[],struct sspecies spec[],struct sconnections connections[],
              int imode,int clumptype)
{    int i,ism,isp;
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
        for (i=0;i<pu[ipu].richness;i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;

            rAmount = SM[ism].amount;

            if (spec[isp].target2 && rAmount > 0)
            {
               if (imode == 1)
               {
                  AddNewPU(ipu,isp,connections,spec,pu,SM,clumptype);
               }
               else
               {
                   RemPu(ipu,isp,connections,spec,pu,SM,clumptype);
               }
               if (spec[isp].occurrence < 0)
               {
                  printf("Warning Warning ! isp %i occ %i \n",isp,spec[isp].occurrence);
               }
            } /* Type 4 species and this will impact them */
            else
            {
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
                appendTraceFile("DoChange ipu %i isp %i spec.amount %g imode %i\n",
                                     ipu,isp,spec[isp].amount,imode);
                #endif
            }  /* None clumping species */

            if (spec[isp].sepnum>0) /* Count separation but only if it is possible that it has changed */
               if ((imode ==1 && spec[isp].separation < spec[isp].sepnum) || (imode == -1 && spec[isp].separation >1))
                  spec[isp].separation = CountSeparation2(isp,0,NULL,puno,R,pu,SM,spec,0);
        } /* Invoke Species Change */

     reserve->total = reserve->cost + reserve->connection + reserve->penalty + reserve->probability1D + reserve->probability2D;
}

void DoQuantumChange(int puno,int *R,struct scost *reserve,struct scost change,
                     struct spustuff pu[],struct spu SM[],struct sspecies spec[],struct sconnections connections[],
                     int clumptype,int iFluctuationCount,int *PUChosen)
{
     // We accept a whole bunch of changes in one, passed in by Quantum annealing.
     int i,j,ipu,ism,isp,imode;
     double rAmount;

     #ifdef DEBUG_QA
     appendTraceFile("DoQuantumChange start\n");
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
            for (i=0;i<pu[ipu].richness;i++)
            {
                ism = pu[ipu].offset + i;
                isp = SM[ism].spindex;

                rAmount = SM[ism].amount;

                if (spec[isp].target2 && rAmount > 0)
                {
                   if (imode == 1)
                   {
                      AddNewPU(ipu,isp,connections,spec,pu,SM,clumptype);
                   }
                   else
                   {
                       RemPu(ipu,isp,connections,spec,pu,SM,clumptype);
                   }
                   if (spec[isp].occurrence < 0)
                   {
                      printf("Warning Warning ! isp %i occ %i \n",isp,spec[isp].occurrence);
                   }
                } /* Type 4 species and this will impact them */
                else
                {
                    spec[isp].occurrence += (rAmount > 0)*imode;
                    spec[isp].amount += rAmount*imode;

                    if (spec[isp].amount < 0.0001)
                       if (spec[isp].amount > -0.0001)
                          spec[isp].amount = 0;

                    #ifdef ANNEALING_TEST
                    appendTraceFile("DoChange ipu %i isp %i spec.amount %g imode %i\n",
                                         ipu,isp,spec[isp].amount,imode);
                    #endif
                }  // No clumping species

                if (spec[isp].sepnum>0) // Count separation but only if it is possible that it has changed
                   if ((imode ==1 && spec[isp].separation < spec[isp].sepnum) || (imode == -1 && spec[isp].separation >1))
                      spec[isp].separation = CountSeparation2(isp,0,NULL,puno,R,pu,SM,spec,0);
            } // Invoke Species Change
     }

     reserve->total = reserve->cost + reserve->connection + reserve->penalty + reserve->probability1D + reserve->probability2D;

     #ifdef DEBUG_QA
     appendTraceFile("DoQuantumChange end\n");
     #endif
} // Do Quantum Change

/*****************************************************/
/*********** Post Processing *************************/
/*****************************************************/

/***  Counts the number of species missing from the reserve ****/
int CountMissing(int spno,struct sspecies spec[],double misslevel,double *shortfall,double *rMinimumProportionMet)
{
    int i,isp = 0;
    double rProportionMet;

    *shortfall = 0;
    *rMinimumProportionMet = 1;
    for (i=0;i<spno;i++)
    {
        rProportionMet = 1;

        if (spec[i].target > 0)
           if (spec[i].amount < spec[i].target)
           {
              *shortfall += spec[i].target - spec[i].amount;
              rProportionMet = spec[i].amount / spec[i].target;

              if (rProportionMet < *rMinimumProportionMet)
                 *rMinimumProportionMet = rProportionMet;

              #ifdef DEBUG_COUNTMISSING
              appendTraceFile("CountMissing i %i target %g amount %g shortfall %g\n",i,spec[i].target,spec[i].amount,*shortfall);
              #endif
           }
        if (spec[i].targetocc > 0)
           if (spec[i].occurrence < spec[i].targetocc)
           {
              *shortfall += spec[i].targetocc - spec[i].occurrence;
              rProportionMet = spec[i].occurrence / spec[i].targetocc;

              if (rProportionMet < *rMinimumProportionMet)
                 *rMinimumProportionMet = rProportionMet;
           }
        if (spec[i].target)
           if (spec[i].amount/spec[i].target < misslevel)
           {
              isp++;
              continue;
           }

        if (spec[i].targetocc)
           if ((double)spec[i].occurrence/(double)spec[i].targetocc < misslevel)
           {
              isp++;
              continue;
           }
        if (spec[i].sepdistance && spec[i].separation < 3)
        {
           isp++;  /* count species if not met separation and not already counted */
        }
    }

    #ifdef DEBUG_COUNTMISSING
    appendTraceFile("CountMissing shortfall %g\n",*shortfall);
    #endif

    return(isp);
}  /* CountMissing*/

// ********* Connection Cost Type 2 **************
// **  Requires R[]. imode2 = 0 there is no negative cost for removing connection, we are calling from ReserveCost
//                         or 1 there is a negative cost for removing connection, we are calling from Annealing
//                   imode = -1 we are removing the planning unit from a reserve, calling from Annealing
//                        or 1  we are adding the planning unit to a reserve, or it is already in reserve
double ConnectionCost2(int ipu,struct sconnections connections[],int *R,int imode,int imode2,double cm)
{
       double fcost, rDelta;
       struct sneighbour *p;
       int R_pu1,i;

       #ifdef DEBUG_CONNECTIONCOST2
       if (asymmetricconnectivity)
          if (imode2)
          {
             appendTraceFile("ConnectionCost2 start puid %i imode %i imode2 %i\n",pu[ipu].id,imode,imode2);
          }
       #endif

       fcost = connections[ipu].fixedcost*imode;
       p = connections[ipu].first;

       if (asymmetricconnectivity)
       {
          while (p) // treatment for asymmetric connectivity
          {
                if (imode2) // calling from Annealing
                {
                   #ifdef DEBUG_CONNECTIONCOST2
                   rDelta = 0;
                   #endif

                   if (imode == 1)
                      R_pu1 = 0;
                   else
                       R_pu1 = 1;

                   if (p->connectionorigon)
                   {
                      if (R[p->nbr] == 0)
                      {
                         if (R_pu1 == 1)
                         {
                            rDelta = -1*p->cost;
                            fcost += rDelta;
                         }
                         else
                         {
                            rDelta = p->cost;
                            fcost += rDelta;
                         }
                      }
                   }
                   else
                   {
                      if (R[p->nbr] == 1 || R[p->nbr] == 2)
                      {
                         if (R_pu1 == 1)
                         {
                            rDelta = p->cost;
                            fcost += rDelta;
                         }
                         else
                         {
                            rDelta = -1*p->cost;
                            fcost += rDelta;
                         }
                      }
                   }

                   #ifdef DEBUG_CONNECTIONCOST2
                   appendTraceFile("ConnectionCost2 puidnbr %i Rnbr %i connectionorigon %i delta %g\n",
                                        pu[p->nbr].id,R[p->nbr],p->connectionorigon,rDelta);
                   #endif
                }
                else // calling from ReserveCost
                {
                    if (R[p->nbr] == 0)
                       if (p->connectionorigon)
                       {
                          rDelta = p->cost;
                          fcost += rDelta;
                       }
                }

                p = p->next;
          }
       }
       else
       {
           while (p) // treatment for symmetric connectivity
           {
                 if (fOptimiseConnectivityIn == 1)
                 {  // optimise for "Connectivity In"
                     if (R[p->nbr] == 1 || R[p->nbr] == 2)
                     {
                        rDelta = imode*p->cost;
                        fcost += rDelta;
                     }
                     else
                     {
                         rDelta = imode*imode2*p->cost*-1;
                         fcost += rDelta;
                     }
                 }
                 else
                 {   // optimise for "Connectivity Edge"
                     if (R[p->nbr] == 1 || R[p->nbr] == 2)
                     {
                        rDelta = imode*imode2*p->cost*-1;
                        fcost += rDelta;
                     }
                     else
                     {
                         rDelta = imode*p->cost;
                         fcost += rDelta;
                     }
                 }

                 p = p->next;
           }
       }

       #ifdef DEBUG_CONNECTIONCOST2
       if (asymmetricconnectivity)
          if (imode2)
          {
             for (i=puno-1;i>-1;i--)
             {
                 appendTraceFile("puid%i R%i\n",pu[i].id,R[i]);
             }

             appendTraceFile("ConnectionCost2 end puid %i connection %g\n",pu[ipu].id,fcost);
          }
       #endif

       return(fcost*cm);
}/***** Connection Cost Type 2 *****************/

void ComputeConnectivityIndices(double *rConnectivityTotal,double *rConnectivityIn,
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

         if (R[i]==1 || R[i] == 2) // add to 'in' or 'edge'
         {
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
         else // add to 'out' or 'edge'
         {
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

/*************** Reporting Value of a Reserve *******/
void PrintResVal (int puno, int spno,int *R,struct scost reserve,
                  struct sspecies spec[],double misslevel)
{
     int i, isp = 0;
     double connectiontemp = 0, shortfall, rMPM;//, rConnectivityTotal = 0,rConnectivityIn = 0,rConnectivityEdge = 0,rConnectivityOut = 0;

     #ifdef DEBUG_PRINTRESVALPROB
     appendTraceFile("PrintResVal start\n");
     #endif

     isp = CountMissing(spno,spec,misslevel,&shortfall,&rMPM);

     //ComputeConnectivityIndices(rConnectivityTotal,rConnectivityIn,rConnectivityEdge,rConnectivityOut,
     //                           puno,R,connections);

     for (i=0;i<puno;i++)
         if (R[i]==1 || R[i] == 2)
         {
            connectiontemp += ConnectionCost2(i,connections,R,1,0,1);
         }

     #ifdef DEBUG_PRINTRESVALPROB
     appendTraceFile("PrintResVal missing %i connectiontemp %g\n",isp,connectiontemp);
     #endif

     displayProgress("Value %.1f Cost %.1f PUs %i Connection %.1f Missing %i Shortfall %.2f Penalty %.1f MPM %.1f\n",
                     reserve.total,reserve.cost,reserve.pus,connectiontemp,isp,shortfall,reserve.penalty,rMPM);

     if (fProb1D == 1)
        displayProgress(" Probability1D %.1f",reserve.probability1D);
     if (fProb2D == 1)
        displayProgress(" Probability2D %.1f",reserve.probability2D);

     displayProgress("\n");

     #ifdef DEBUG_PRINTRESVALPROB
     appendTraceFile("PrintResVal end\n");
     #endif
} /******* Print Reserve Value *********/

/****** Change the Cost *****/
/* This routine changes the values in the cost structure (such as making them negative ) */
/* Useful for altering the 'change' variable  */
void ChangeCost(struct scost *cost,double changemult)
{
     cost->total *= changemult;
     cost->connection *= changemult;
     cost->penalty *= changemult;
     cost->cost *= changemult;
     cost->probability1D *= changemult;
     cost->probability2D *= changemult;
}

/***************************************************************************/
/************** Annealing Specific Functions *******************************/
/***************************************************************************/

void ConnollyInit(int puno,int spno,struct spustuff pu[],typeconnection connections[],typesp spec[],
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
     appendTraceFile("ConnollyInit start\n");
     if (verbosity > 4)
     {
        sprintf(sRun,"%i",irun);
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen("debug_maropt_ConnollyInit_.csv") + strlen(sRun) + 2, sizeof(char));
        strcpy(writename,fnames.outputdir);
        strcat(writename,"debug_maropt_ConnollyInit_");
        strcat(writename,sRun);
        strcat(writename,".csv");
        fp = fopen(writename,"w");
        if (fp==NULL)
           displayErrorMessage("cannot create debug_maropt_ConnollyInit file %s\n",writename);
        free(writename);
        fprintf(fp,"i,ipu,puid,old R,imode,R,total,max,min\n");
     }
     #endif

     #ifdef DEBUG_PROB1D
     appendTraceFile("ConnollyInit A\n");
     #endif

     localdelta = 1E-10;

     InitReserve(puno,prop,R);

     #ifdef DEBUG_PROB1D
     appendTraceFile("ConnollyInit B\n");
     #endif

     AddReserve(puno,pu,R);
     if (aggexist)
        ClearClumps(spno,spec,pu,SM);

     #ifdef DEBUG_PROB1D
     appendTraceFile("ConnollyInit C\n");
     #endif

     ReserveCost(puno,spno,R,pu,connections,SM,cm,spec,aggexist,&reserve,clumptype);

     #ifdef DEBUG_PROB1D
     appendTraceFile("ConnollyInit D\n");
     #endif

     for (i=1;i<= (*anneal).iterations/100; i++)
     {
         ipu = RandNum(puno);
         iOldR = R[ipu];
         imode = R[ipu]==1?-1:1;

         CheckChange(-1,ipu,spno,puno,pu,connections,spec,SM,R,cm,imode,&change,&reserve,0,0,0,0,clumptype);
         DoChange(ipu,puno,R,&reserve,change,pu,SM,spec,connections,imode,clumptype);

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

     appendTraceFile("ConnollyInit end\n");
     if (verbosity > 4)
        fclose(fp);
} // Init Annealing Schedule According to Connolly Scheme

/******* Adaptive Annealing 2 *********************/
/**** Initial Trial Runs. Run for some small time to establish sigma. ****/
void AdaptiveInit(int puno,int spno,double prop,int *R,
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
        InitReserve(puno,prop,R);
        AddReserve(puno,pu,R);
        /* Score Random reserve */
        ReserveCost(puno,spno,R,pu,connections,SM,cm,spec,aggexist,&cost,clumptype);
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

} /**** Adaptive Annealing Initialisation *****/

/**** Set TInitial from this as well ****/

/**** Function to decrement T and decide if it is time to stop? *****/
void AdaptiveDec(struct sanneal *anneal)
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

} /* Adaptive Decrement. Sets the new temperature based on old values */

void ThermalAnnealing(int spno, int puno, struct sconnections connections[],int R[], double cm,
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

     appendTraceFile("ThermalAnnealing start iterations %ld\n",anneal.iterations);
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
        fprintf(ttfp,"%li,%f,%li,%f,%i,%f,%f,%f,%f"
                ,itime,costthresh,iGoodChange,reserve->total
                ,reserve->pus,reserve->cost,reserve->connection,reserve->penalty,reserve->shortfall);
        if (fProb1D == 1)
           fprintf(ttfp,",%f",reserve->probability1D);
        if (fProb2D == 1)
           fprintf(ttfp,",%f",reserve->probability2D);
        fprintf(ttfp,",%li\n",ipu);
        // iteration,threshold,dochange,total,pus,cost,connectivity,penalty,probability

        sprintf(tempname2,"%s_anneal_zones%05i.csv",savename,irun%10000);
        writename = (char *) calloc(strlen(fnames.outputdir) + strlen(tempname2) + 2, sizeof(char));
        //strcpy(writename,fnames.outputdir);
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

     displayProgress2("  Main ThermalAnnealing Section.\n");

     rThreshold = costthresh;
     costthresh = rThreshold * rStartDecMult;
     rTemperature = 1;

     for (itime = 1;itime<=anneal.iterations;itime++)
     {
         do
         {
           ipu = RandNum(puno);
         } while (R[ipu] > 1); /*  Select a PU at random */

         itemp = R[ipu] == 1 ? -1 : 1;  /* Add or Remove PU ? */

         CheckChange(itime,ipu,spno,puno,pu,connections,spec,SM,R,cm,itemp,change,reserve,
                     costthresh,tpf1,tpf2,(double) itime/ (double) anneal.iterations,clumptype);
         /* Need to calculate Appropriate temperature in GoodChange or another function */
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
               AdaptiveDec(&anneal);
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
         if (GoodChange(*change,anneal.temp)==1)
         {
            iGoodChange = 1;

            ++ichanges;
            DoChange(ipu,puno,R,reserve,*change,pu,SM,spec,connections,itemp,clumptype);
            if (fnames.savesnapchanges && !(ichanges % fnames.savesnapfrequency))
            {
               if (repeats > 1)
                  sprintf(tempname1,"_r%05i",irun);
               else
                   tempname1[0] = 0;
              if (fnames.savesnapchanges == 3)
                 sprintf(tempname2,"%s_snap%sc%05li.csv",savename,tempname1,++snapcount%10000);
              else
              if (fnames.savesnapchanges == 2)
                 sprintf(tempname2,"%s_snap%sc%05li.txt",savename,tempname1,++snapcount%10000);
              else
                  sprintf(tempname2,"%s_snap%sc%05li.dat",savename,tempname1,++snapcount%10000);
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
            fprintf(fp,"%li,%li,%i,%li,%li,%i,%li,%f,%f,%f,%f,%f\n"
                    ,itime,ipu,pu[ipu].id,iPreviousR,itemp,R[ipu],iGoodChange,change->total,change->cost,change->connection,change->penalty,anneal.temp);

         if (fnames.saveannealingtrace)
         {
            iRowCounter++;
            if (iRowCounter > iRowLimit)
               iRowCounter = 1;

            if (iRowCounter == 1)
            {
               fprintf(Rfp,"%li",itime);

               fprintf(ttfp,"%li,%f,%li,%f,%i,%f,%f,%f,%f"
                       ,itime,costthresh,iGoodChange,reserve->total
                       ,reserve->pus,reserve->cost,reserve->connection,reserve->penalty,reserve->shortfall);
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
       ReserveCost(puno,spno,R,pu,connections,SM,cm,spec,aggexist,reserve,clumptype);
       displayProgress1("  ThermalAnnealing:");

       #ifdef DEBUG_PRINTRESVALPROB
       appendTraceFile("before PrintResVal ThermalAnnealing:\n");
       #endif

       PrintResVal(puno,spno,R,*reserve,spec,misslevel);

       #ifdef DEBUG_PRINTRESVALPROB
       appendTraceFile("after PrintResVal ThermalAnnealing:\n");
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
}  // ThermalAnnealing

void QuantumAnnealing(int spno, int puno, struct sconnections connections[],int R[], double cm,
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

     appendTraceFile("QuantumAnnealing start iterations %ld decay type %s proportion %f decay A %f decay B %f acceptance probability %f saveannealingtrace %i\n",
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

     displayProgress2("  Main QuantumAnnealing Section.\n");

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
         //iFluctuationCount = itime;
         //rFluctuationMagnitude = iFluctuationCount / puno;

         #ifdef DEBUG_QA
         appendTraceFile("QuantumAnnealing rFluctuationMagnitude %f iFluctuationCount %i\n",
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
                 j = RandNum(puno);

                 #ifdef DEBUG_QA
                 appendTraceFile("QuantumAnnealing j %i PUChosen[j] %i R[j] %i \n",j,PUChosen[j],R[j]);
                 #endif
               }
               while ((PUChosen[j] > 0) || (R[j] > 1));
               // select PU's at random that are not already chosen or locked

               #ifdef DEBUG_QA
               appendTraceFile("QuantumAnnealing chose ipu %i\n",j);
               #endif


               PUChosen[j] = 1;
           }

           // compute objective function score with these bits flipped
           CheckQuantumChange(spno,puno,pu,connections,spec,SM,R,cm,change,reserve,
                              costthresh,tpf1,tpf2,(double) itime/ (double) anneal.iterations,
                              clumptype,iFluctuationCount,PUChosen);

           // we only accept good changes

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
           if (GoodQuantumChange(*change,rAcceptanceProbability)==1)
           {
              iGoodChange = 1;

              ++ichanges;
              DoQuantumChange(puno,R,reserve,*change,pu,SM,spec,connections,clumptype,iFluctuationCount,PUChosen);
              if (fnames.savesnapchanges && !(ichanges % fnames.savesnapfrequency))
              {
                 if (repeats > 1)
                    sprintf(tempname1,"_r%05i",irun);
                 else
                     tempname1[0] = 0;
                if (fnames.savesnapchanges == 3)
                   sprintf(tempname2,"%s_snap%sc%05li.csv",savename,tempname1,++snapcount%10000);
                else
                if (fnames.savesnapchanges == 2)
                   sprintf(tempname2,"%s_snap%sc%05li.txt",savename,tempname1,++snapcount%10000);
                else
                    sprintf(tempname2,"%s_snap%sc%05li.dat",savename,tempname1,++snapcount%10000);
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
              fprintf(fp,"%li,%li,%li,%f,%f,%f,%f,%f\n"
                      ,itime,itemp,iGoodChange,change->total,change->cost,change->connection,change->penalty,anneal.temp);

           if (fnames.saveannealingtrace)
           {
              iRowCounter++;
              if (iRowCounter > iRowLimit)
                 iRowCounter = 1;

              if (iRowCounter == 1)
              {
                 fprintf(Rfp,"%li",itime);

                 fprintf(ttfp,"%li,%f,%li,%f,%i,%f,%f,%f\n"
                         ,itime,costthresh,iGoodChange,reserve->total
                         ,reserve->pus,reserve->cost,reserve->connection,reserve->penalty);
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
         }
         else
         {
             // force algorithm to drop out of iterations loop
             iIterations = itime - 1;
             itime = anneal.iterations;
         }

     } /* Run Through Annealing */

     free(PUChosen);

     /** Post Processing  **********/
     if (verbosity >1)
     {
       ReserveCost(puno,spno,R,pu,connections,SM,cm,spec,aggexist,reserve,clumptype);
       displayProgress1("  QuantumAnnealing:");

       #ifdef DEBUG_PRINTRESVALPROB
       appendTraceFile("before PrintResVal QuantumAnnealing:\n");
       #endif

       PrintResVal(puno,spno,R,*reserve,spec,misslevel);

       #ifdef DEBUG_PRINTRESVALPROB
       appendTraceFile("after PrintResVal QuantumAnnealing:\n");
       #endif
     }

     if (aggexist) ClearClumps(spno,spec,pu,SM);

     #ifdef DEBUGTRACEFILE
     if (verbosity > 4)
        fclose(fp);
     #endif

     if (fnames.saveannealingtrace)
     {
        fclose(ttfp);
        fclose(Rfp);
     }
     appendTraceFile("QuantumAnnealing end iterations %ld tests %li\n",iIterations,iTests);

}  /* Main Quantum Annealing Function */

/* ANNEALING.C END */
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

/* optimisation functions written by Matt Watts */

void siftDown_bs(struct binsearch numbers[], int root, int bottom, int array_size)
{
     int done, maxChild;
     typebinsearch temp;

     done = 0;
     while ((root*2 <= bottom) && (!done))
     {
           if (root*2 < array_size)
           {
              if (root*2 == bottom)
                 maxChild = root * 2;
              else if (numbers[root * 2].name > numbers[root * 2 + 1].name)
                      maxChild = root * 2;
                   else
                       maxChild = root * 2 + 1;

              if (numbers[root].name < numbers[maxChild].name)
              {
                 temp = numbers[root];
                 numbers[root] = numbers[maxChild];
                 numbers[maxChild] = temp;
                 root = maxChild;
              }
              else
                  done = 1;
           }
           else
               done = 1;
     }
}

void heapSort_bs(struct binsearch numbers[], int array_size)
{
     int i;
     typebinsearch temp;

     for (i = (array_size / 2)-1; i >= 0; i--)
         siftDown_bs(numbers, i, array_size, array_size);

     for (i = array_size-1; i >= 1; i--)
     {
         temp = numbers[0];
         numbers[0] = numbers[i];
         numbers[i] = temp;
         siftDown_bs(numbers, 0, i-1, array_size);
     }
}

void PrepareBinarySearchArrays(int puno, int spno, struct spustuff PU[], typesp spec[],
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
     heapSort_bs((* PULookup),puno);
     heapSort_bs((* SPLookup),spno);

     if (verbosity > 3)
        writeBinarySearchArrays("after",fnames,puno,spno,(* PULookup),(* SPLookup));
}

void TestFastNameToPUID(int puno, struct binsearch PULookup[], struct spustuff PU[], struct sfname fnames)
{
     FILE *fp;
     int i;
     char *writename;

     writename = (char *) calloc(strlen(fnames.inputdir) + strlen("TestFastNameToPUID.csv") + 2, sizeof(char));
     strcpy(writename,fnames.inputdir);
     strcat(writename,"TestFastNameToPUID.csv");
     if ((fp = fopen(writename,"w"))==NULL)
          displayErrorMessage("cannot create TestFastNameToPUID file %s\n",writename);
     free(writename);
     fputs("name,index,bin search index\n",fp);
     for (i=0;i<puno;i++)
          fprintf(fp,"%d,%d,%d\n",PU[i].id,i,FastNameToPUID(puno,PU[i].id,PULookup));
     fclose(fp);
}


int FastNameToPUID(int puno,int name, struct binsearch PULookup[])
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
          }
          else
          {
              iTop = iCentre + 1;
              iCount = iBottom - iTop + 1;
              iCentre = iTop + floor(iCount / 2);
          }
    }
    return(PULookup[iCentre].index);
}

void TestFastNameToSPID(int spno, struct binsearch SPLookup[], typesp spec[], struct sfname fnames)
{
     FILE *fp;
     int i;
     char *writename;

     writename = (char *) calloc(strlen(fnames.inputdir) + strlen("TestFastNameToSPID.csv") + 2, sizeof(char));
     strcpy(writename,fnames.inputdir);
     strcat(writename,"TestFastNameToSPID.csv");
     if ((fp = fopen(writename,"w"))==NULL)
        displayErrorMessage("cannot create TestFastNameToSPID file %s\n",writename);
     free(writename);
     fputs("name,index,bin search index\n",fp);
     for (i=0;i<spno;i++)
         fprintf(fp,"%d,%d,%d\n",spec[i].name,i,FastNameToSPID(spno,spec[i].name,SPLookup));
     fclose(fp);
}


int FastNameToSPID(int spno,int name, struct binsearch SPLookup[])
{
    /* use a binary search to find the index of planning unit "name" */
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

/*********************************************/
/*********  Clump Utilities ******************/
/********************************************/

// Clear a single Clump
void ClearClump(int isp,struct sclumps *target,struct spustuff pu[],
                struct spu SM[])
{
    struct sclumppu *ppu;

     /* Remove all links from this clump */
    while (target->head) {
        ppu = target->head;
        if (rtnClumpSpecAtPu(pu,SM,ppu->puid,isp) == target->clumpid) /* in case pu is in new clump */
            setClumpSpecAtPu(pu,SM,ppu->puid,isp,0);
        target->head = ppu->next;
        free(ppu);
        DebugFree(sizeof(struct sclumppu));
    } /* Remove all links from this clump */
}

// Function does this cut clump?
// Returns the value of the fragmented clumps if the given PU were removed
// If imode = 1 then it will also do a separation count
int ClumpCut(int isp,struct spustuff pu[],
        struct sspecies spec[],struct sclumps *clump,
        struct sclumppu *clumppu,struct sconnections connections[],struct spu SM[],
        double *totalamount,int *totalocc,
        int *iseparation, int imode,int clumptype)
{
    int ineighbour = 0,iclumps = 0;
    struct slink{int id; struct slink *next;} *head = NULL, *newhead,*thead, *clumplist, *clumpcurr;
    struct sneighbour *pnbr;
    struct sclumps *spclump = NULL, *newpclump;
    struct sclumppu *pclumppu;
    double clumpamount, rAmount;
    int iocc;

    *totalamount = 0;
    *totalocc = 0;

    /* Set up spclump for counting Separation */
    if (imode)
    {
        newpclump = (struct sclumps *) malloc(sizeof(struct sclumps));
        newpclump->clumpid = clumppu->puid;
        newpclump->amount = 0;
        newpclump->next = spclump;
        spclump = newpclump;
    }

    /** Generate list of all neighbours and count them **/
    /*First check if there are no neighbours then exit. **/
      /* return null for no clump cut done and need to do separation count */
    if (connections[clumppu->puid].nbrno == 0)
    {
        if (imode)
        {
            *iseparation = CountSeparation(isp,spclump,pu,SM,spec,0);
            free(spclump);
            DebugFree(sizeof(struct sclumps));
        }
        return(0);
    }

    for (pnbr = connections[clumppu->puid].first; pnbr;pnbr = pnbr->next)
    {
        if (rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp) == clump->clumpid)
        {
            ineighbour++;
            newhead = (struct slink *) malloc(sizeof(struct slink));
            newhead->id = pnbr->nbr;
            newhead->next = head;
            head = newhead;
        } // If neighbour is part of the same clump
    } // For cycling through all neighbours

    if (ineighbour <= 1)
    { // One or fewer neighbours
        if (imode)
        { // separation distance called
            for(pclumppu=clump->head;pclumppu;pclumppu=pclumppu->next)
             if (pclumppu != clumppu)
             {
                newpclump = (struct sclumps *) malloc(sizeof(struct sclumps));
                newpclump->clumpid = pclumppu->puid;
                newpclump->amount = clump->amount - rtnAmountSpecAtPu(pu,SM,clumppu->puid,isp);
                newpclump->next = spclump;
                spclump = newpclump;

             } /* found someone in the clump who is not being removed */

             *iseparation = CountSeparation(isp,spclump,pu,SM,spec,0);
        }
        else
            (*iseparation = spec[isp].sepnum);
            
        if (head)
        {
            free(head);
            DebugFree(sizeof(struct slink));
        }
        
        while (spclump)
        {
            newpclump = spclump;
            spclump = spclump->next;
            free(newpclump);
            DebugFree(sizeof(struct sclumps));
        }  // clearing up spclump
        
        rAmount = rtnAmountSpecAtPu(pu,SM,clumppu->puid,isp);
        *totalamount = PartialPen4(isp,clump->amount-rAmount,spec,clumptype);
        *totalocc = (clump->occs - (rAmount > 0))*(*totalamount > 0); // count only if still valid size
        
        return(0);
    }

    // More than one neighbour. Can they form their own clump?
    // Put first neighbour at head of new list
    while (head)
    {
          clumpamount = 0;
          iclumps++;
          clumplist = (struct slink *) malloc(sizeof(struct slink));
          clumplist->next = NULL;
          clumplist->id = head->id;
          clumpcurr = clumplist;
          newhead = head;
          head = head->next;
          free(newhead);  /* move first site from head to clumplist */
          DebugFree(sizeof(struct slink));
          ineighbour--;
          do
          {
            for (pnbr = connections[clumpcurr->id].first;pnbr;pnbr = pnbr->next)
            {
                if (rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp) == clump->clumpid && pnbr->nbr != clumppu->puid) // if neighbour in clump but not cut out one   
                {
                   for (newhead = clumplist;newhead && newhead->id != pnbr->nbr;newhead= newhead->next)
                       ; // Cycle through clumplist looking to see if this fellow is already in it
                       
                   if (!newhead)
                   {
                      newhead = (struct slink *) malloc(sizeof(struct slink));
                      newhead->id = pnbr->nbr;
                      newhead->next = clumpcurr->next;
                      clumpcurr->next = newhead;  /* put this item in my clumplist */
                      
                      // go through neighbour list and see if this one is there
                      for (newhead=head;newhead && newhead->id != pnbr->nbr;newhead = newhead->next)
                          ; // find this item on the neighbour list
                          
                      if (newhead && newhead->id == pnbr->nbr)
                      {
                         ineighbour--;
                         if (newhead == head)
                            head = newhead->next;
                         else
                         {
                             for (thead=head;thead->next != newhead; thead = thead->next)
                                 ; // find link before the one to be removed
                                 
                             thead->next = newhead->next;
                         } // remove link that is not head
                         
                         free(newhead);
                         DebugFree(sizeof(struct slink));
                      } // A new neighbour is taken into account
                   } // Adding a novel neighbour to list
                } // found a neighbour in clump which isn't the one being cut
            } // cycling through every neighbour on this clump

            // point to next one on list but keep clump head where it is
            clumpcurr = clumpcurr->next;

          } while (clumpcurr); // if you've run out of new list then...

          iocc = 0;
          for (newhead=clumplist;newhead;newhead=newhead->next)
          {
              rAmount = rtnAmountSpecAtPu(pu,SM,newhead->id,isp);
              clumpamount += rAmount; // find total amount
              iocc += (rAmount > 0);
          }
          
          *totalamount += PartialPen4(isp,clumpamount,spec,clumptype);
          if (PartialPen4(isp,clumpamount,spec,clumptype))
             *totalocc += iocc;

          if (imode)
             for (newhead=clumplist;newhead;newhead=newhead->next)
             {
                 newpclump = (struct sclumps *)malloc(sizeof(struct sclumps));
                 newpclump->clumpid = newhead->id;
                 newpclump->amount = clumpamount;
                 newpclump->next = spclump;
                 spclump = newpclump;
             } // stick this clump into my clump list for separation purposes

             // clean up all lists
             while (clumplist)
             {
                   clumpcurr = clumplist;
                   clumplist = clumplist->next;
                   free(clumpcurr);
                   DebugFree(sizeof(struct slink));
             } // clean up clumplist
    } // Continue clump formation whilst there are members in the list

    if (imode)
    {
       *iseparation = CountSeparation(isp,spclump,pu,SM,spec,0);
       while (spclump)
       {
             newpclump = spclump;
             spclump = spclump ->next;
             free(newpclump);
             DebugFree(sizeof(struct sclumps));
       } /* clean up separation clump list */
    }
    else
        *iseparation = spec[isp].sepnum;

    while (head)
    {
        newhead = head;
        head = head->next;
        free(newhead);
        DebugFree(sizeof(struct slink));
    } // clean up neighbour list
    
    return(iclumps);
} // Function Clump Cut.. Do I cut ?

// Clear Clumps
// This is for clean up purposes
void ClearClumps(int spno,struct sspecies spec[],struct spustuff pu[],
                 struct spu SM[])
{
     int i;
     struct sclumps *pclump;

     for (i=0;i<spno;i++)
     {
        while (spec[i].head)
        {
              ClearClump(i,spec[i].head,pu,SM);
              pclump = spec[i].head;
              spec[i].head = spec[i].head->next;
              free(pclump);
        }  // Remove each clump
        
        spec[i].clumps = 0;
     } // Clear clump for each species
} // Clear Clumps

// Add New Clump
struct sclumps *AddNewClump(int isp,int ipu,struct sspecies spec[],struct spustuff pu[],struct spu SM[])
{
       int iclumpno = 0;
       struct sclumps *pclump,*pnewclump;
       struct sclumppu *pnewclumppu;
       double rAmount;

       // find good clump number
       pclump = spec[isp].head;
       if (!pclump)
          iclumpno = 1;
    
       while (!iclumpno)
       {
             if (!pclump->next)
             {
                iclumpno = pclump->clumpid+1;
                break;
             } // I've found the end of the list
             
             if (pclump->next->clumpid - pclump->clumpid > 1)
             {
                iclumpno = pclump->clumpid+1;
                continue;
             } // Looking for good number
             
             pclump = pclump->next;
       } // Find first available clump number

       setClumpSpecAtPu(pu,SM,ipu,isp,iclumpno);
       pnewclump = (struct sclumps *) malloc(sizeof(struct sclumps));
       pnewclump->clumpid = iclumpno;
       if (spec[isp].head)
       {
          pnewclump->next = pclump->next;
          pclump->next = pnewclump;
       } // Stick clump into correct location
       else
       {
           spec[isp].head = pnewclump;
           pnewclump->next = NULL;
       } // First clump on the block
    
       // Add first clumppu to this new clump
       pnewclumppu = (struct sclumppu *) malloc(sizeof(struct sclumppu));
       pnewclumppu->puid = ipu;
       pnewclumppu->next = NULL;
       pnewclump->head = pnewclumppu;
       rAmount = rtnAmountSpecAtPu(pu,SM,ipu,isp);
       pnewclump->amount = rAmount;
       pnewclump->occs = (rAmount > 0);

       spec[isp].clumps++;

       return(pnewclump);

}  // Add New Clump

// ADD NEW PU
// Add New Planning Unit for a given Species
void AddNewPU(int ipu,int isp,struct sconnections connections[],struct sspecies spec[],struct spustuff pu[],
              struct spu SM[],int clumptype)
{
     int ineighbours = 0;
     int iclumpno, iClump;
     struct sneighbour *pnbr;
     struct sclumps *pclump, *pnewclump, *ptempclump;
     struct sclumppu *pnewclumppu;
     double ftemp, rAmount;

     pnbr = connections[ipu].first;
     while (pnbr)
     {
           // Check all the neighbours to see if any are already in clumps */
           iClump = rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp);
           if (iClump > 0)
           {
              // Neighbour that is part of clump
              ineighbours++;
              if (ineighbours == 1)
              {
                 // Join to the first clump that is also a neighbour
                 iclumpno = iClump;
                 for (pclump = spec[isp].head; pclump->clumpid != iclumpno;pclump = pclump->next)
                     ;
                     
                 pnewclumppu = (struct sclumppu *) malloc(sizeof(struct sclumppu));
                 pnewclumppu->puid = ipu;
                 pnewclumppu->next = pclump->head;
                 setClumpSpecAtPu(pu,SM,pnewclumppu->puid,isp,iclumpno);
                 pclump->head = pnewclumppu;

                 // Remove old value for this clump
                 ftemp = PartialPen4(isp,pclump->amount,spec,clumptype);
                 spec[isp].amount -= ftemp;
                 spec[isp].occurrence -= pclump->occs *(ftemp > 0);
                 rAmount = rtnAmountSpecAtPu(pu,SM,ipu,isp);
                 pclump->occs += (rAmount > 0);
                 pclump->amount += rAmount;
              } // Adding the pu to the clump
              else
              {
                  // pclump points to the good clump
                  if (pclump->clumpid != iClump)
                  {
                     // Check if this is a different clump
                     // Join this new clump to the old one
                     for (pnewclump= spec[isp].head; pnewclump->clumpid != iClump;pnewclump = pnewclump->next)
                         ;  // point pnewclump to the joining clump
                         
                     // Run through joining clump and tell all pu's their new number
                     for (pnewclumppu = pnewclump->head;pnewclumppu->next;pnewclumppu=pnewclumppu->next)
                         setClumpSpecAtPu(pu,SM,pnewclumppu->puid,isp,pclump->clumpid);
                         
                     setClumpSpecAtPu(pu,SM,pnewclumppu->puid,isp,pclump->clumpid);
                     
                     // cut out this clump and join it to pclump
                     pnewclumppu->next = pclump->head;
                     pclump->head = pnewclump->head;
                     pclump->amount += pnewclump->amount;
                     pclump->occs += pnewclump->occs;
                     ftemp = PartialPen4(isp,pnewclump->amount,spec,clumptype);
                     spec[isp].amount -= ftemp;
                     spec[isp].occurrence -= pnewclump->occs * (ftemp > 0);

                     // Remove clump head and free memory
                     if (pnewclump == spec[isp].head)
                        spec[isp].head = pnewclump->next;
                     else
                     {
                         for (ptempclump = spec[isp].head;ptempclump->next != pnewclump;ptempclump = ptempclump->next)
                             ; // Find clump just before redundant clump
                             
                         ptempclump->next = pnewclump->next;
                     }

                     free(pnewclump);
                     DebugFree(sizeof(struct sclumps));

                  } // Join the two clumps together
              } // Found another neighbour
           }
           pnbr = pnbr->next;
     } // cycling through all the neighbours

     // Adding a New clump
     if (!ineighbours)
     {
        AddNewClump(isp,ipu,spec,pu,SM);
        ftemp = PartialPen4(isp,rAmount,spec,clumptype);
        spec[isp].amount += ftemp;
        spec[isp].occurrence += (ftemp>0);
     } // Adding a new clump

     // Correcting Amount if new clump not added
     if (ineighbours)
     {
        ftemp = PartialPen4(isp,pclump->amount,spec,clumptype);
        spec[isp].amount += ftemp;
        spec[isp].occurrence += pclump->occs * (ftemp > 0);
     }
} // Add New Pu

/************** REM PU ****************************************/
/*********** Remove a planning unit. Note it is similar to CutClump but actually does action **/
/**************************************************************/
void RemPu(int ipu, int isp,struct sconnections connections[], struct sspecies spec[],struct spustuff pu[],
           struct spu SM[],int clumptype)
{
    int ineighbours = 0;
    struct slink{int id;struct slink *next;} *head = NULL, *newhead, *thead;
    struct sclumps *oldclump,*pclump;
    struct sclumppu *cppu,*ppu, *clumpcurr, *tppu;
    struct sneighbour *pnbr;
    double oldamount,newamount = 0.0, rAmount;
    int newoccs;

    for (oldclump = spec[isp].head;oldclump && oldclump->clumpid != rtnClumpSpecAtPu(pu,SM,ipu,isp); oldclump= oldclump->next)
    ; /* Find the correct clump to remove */
    if (!oldclump)
        displayErrorMessage("Serious error in Remove Type 4 species routine\n");

    for(cppu = oldclump->head;cppu->puid != ipu; cppu = cppu->next)
    ; /* Locate the correct clumppu */
    setClumpSpecAtPu(pu,SM,cppu->puid,isp,0);

    oldamount = PartialPen4(isp,oldclump->amount,spec,clumptype);
    spec[isp].amount -= oldamount;
    spec[isp].occurrence -= oldclump->occs * (oldamount > 0);

    for (pnbr = connections[cppu->puid].first;pnbr;pnbr = pnbr->next)
        if(rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp) == oldclump->clumpid) {
            ineighbours++;
            newhead = (struct slink *)malloc(sizeof(struct slink));
            newhead->id = pnbr->nbr;
            newhead->next = head;
            head = newhead;
        } /* Building the neighbour list */

    if (ineighbours <= 1) {
        rAmount = rtnAmountSpecAtPu(pu,SM,ipu,isp);
        oldclump->amount -= rAmount;
        oldclump->occs -= (rAmount > 0);
        newamount = PartialPen4(isp,oldclump->amount,spec,clumptype);
        newoccs = oldclump->occs * (newamount > 0);
        /* remove clumppu */
        if (cppu == oldclump->head) {
            oldclump->head = cppu->next;
        }
        else {
            for (ppu= oldclump->head;ppu->next != cppu; ppu = ppu->next)
            ; /* find preceding clumppu;*/
            ppu->next = cppu->next;
        }
        free(cppu);
        DebugFree(sizeof(struct sclumppu));
        if (ineighbours < 1)
        {
            if (oldclump == spec[isp].head)
                spec[isp].head = oldclump->next;
            else {
                for (pclump = spec[isp].head;pclump->next != oldclump;pclump = pclump->next)
                ;
                pclump->next = oldclump->next;
            } /* find preceeding clump */
            free(oldclump);
            DebugFree(sizeof(struct sclumps));
            spec[isp].clumps--;
        } /* Removing redundant clump */
        spec[isp].amount += newamount;
        spec[isp].occurrence += newoccs;
        if (head) {
            free(head);
            DebugFree(sizeof(struct slink));
        } /* Only need to free head if ineighbours ==1. Then only 1 item in list */
        return;
    } /* I'm not cutting a clump */

    /* Else create new clumps */
    while (head){
        /* take first element as seed for new clump number */
        pclump = AddNewClump(isp,head->id,spec,pu,SM);
        clumpcurr = pclump->head;
        do {
            for(pnbr=connections[clumpcurr->puid].first;pnbr;pnbr=pnbr->next) {
              if(rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp) == oldclump->clumpid)
              {
                  if (oldclump->head->puid == pnbr->nbr) {
                      ppu = oldclump->head;
                      oldclump->head = ppu->next;
                    } /* cut out old clump of head */
                    else{
                        for (tppu= oldclump->head;tppu->next->puid != pnbr->nbr;tppu= tppu->next)
                        ; /* find preceeding pu in clump */
                        ppu = tppu->next;
                        tppu->next = ppu->next;
                    } /* cut from middle of old clump */
                     ppu->next = clumpcurr->next;
                     clumpcurr->next = ppu;
                     setClumpSpecAtPu(pu,SM,ppu->puid,isp,pclump->clumpid);
                     rAmount = rtnAmountSpecAtPu(pu,SM,ppu->puid,isp);
                     pclump->amount += rAmount;
                     pclump->occs += (rAmount>0);
                     /* Check if it is on neighbours list and if so then remove it from that list*/
                    if (head->id == ppu->puid) {
                        newhead = head;
                        head = newhead->next;
                        free(newhead);
                        DebugFree(sizeof(struct slink));
                    }
                    else {
                        for (newhead= head;newhead->next && newhead->next->id != ppu->puid;newhead= newhead->next)
                        ; /* check if next one on list is same person */
                        if (newhead->next && newhead->next->id == ppu->puid) {
                            thead = newhead->next;
                            newhead->next = thead->next;
                            free(thead);
                            DebugFree(sizeof(struct slink));
                        } /* cut out none head element */
                    }
                } /* This one is worth removing */
            } /* Cycling through all neighbours */
            clumpcurr = clumpcurr->next;
        } while (clumpcurr); /* Continue until you've added every conceivable neighbour */
        spec[isp].amount += PartialPen4(isp,pclump->amount,spec,clumptype);
        spec[isp].occurrence += pclump->occs * (PartialPen4(isp,pclump->amount,spec,clumptype)>0);
        newhead = head;
        head = newhead->next;
        free(newhead);
        DebugFree(sizeof(struct slink));
    } /** Account for every neighbour in my list **/

    /* Every neighbour in local list has been used and every clump formed*/
    /* Remove old clump */
    /* Worry about change in amount and hence score */
    if (oldclump == spec[isp].head) {
        spec[isp].head = oldclump->next;
    }
    else {
        for(pclump=spec[isp].head;pclump->next != oldclump;pclump=pclump->next)
        ; /* find neighbouring clump */
        pclump->next = oldclump->next;
    } /* removing old clump */
    ClearClump(isp,oldclump,pu,SM);
    free(oldclump);
    DebugFree(sizeof(struct sclumps));

} /* Remove a Planning Unit ***/

/*************** Setting Clumps for Species Aggregation Rule**********/
void SetSpeciesClumps(int puno,int R[],struct sspecies spec[],struct spustuff pu[],
                      struct spu SM[],struct sconnections connections[],int clumptype)
{
  int i, ipu, isp, ism;

  for (ipu=0;ipu<puno;ipu++)
      if (pu[ipu].richness)
         for (i=0;i<pu[ipu].richness;i++)
         {
             ism = pu[ipu].offset + i;
             isp = SM[ism].spindex;
             if (spec[isp].target2)
             {
                spec[isp].clumps = 0;
                if ((R[ipu]==1 || R[ipu]==2) && SM[ism].amount > 0 && SM[ism].clump  == 0)
                {
                   AddNewPU(ipu,isp,connections,spec,pu,SM,clumptype);
                }// Add a New planning unit
             }// For each type 4 species
         }

} /******** Set Species clumps *************/

/************ Species Amounts Type 4 **************/
/** Assumes Set Species Clumps has been called **/
void SpeciesAmounts4(int isp,struct sspecies spec[],int clumptype)
{
     double ftemp;
     struct sclumps *pclump;

     for (pclump = spec[isp].head;pclump;pclump= pclump->next)
     {
         ftemp = PartialPen4(isp,pclump->amount,spec,clumptype);
         spec[isp].amount += ftemp;
         spec[isp].occurrence += pclump->occs*(ftemp>0);
     }

} /*** Species Amounts 4 **/

/*** Remove Clump Check ***/
/** returns 0 if any member of clump is non-removable, Ie status == 2 **/
int RemClumpCheck(struct sclumps *pclump,struct spustuff pu[])
{  struct sclumppu *pcpu;
    for (pcpu = pclump->head;pcpu;pcpu = pcpu->next)
        if (pu[pcpu->puid].status == 2)
            return(0);
    return(1);
}

/********* Set Penalties for a given Type 4 Species ***/
/* Returns 1 if the species is a 'bad species' and -1 if it is a 'good species' */
/* Also sticks the penalty into spec[isp].penalty */
int CalcPenaltyType4(int isp,int puno, struct spu SM[],struct sconnections connections[],
                     struct sspecies spec[],struct spustuff pu[],double cm,int clumptype)
{   int i,j,ipu,iputotal = 0;
    int ineighbours = 0,iclumpno,badspecies = 0;
    int *R;
    double totalamount,dummy = 0;
    int idummy;
    double cost = 0.0, connection = 0.0, rAmount;
    struct slink {int id; struct slink *next;} *plist,*plisthead = NULL,*pdiscard;
    struct sneighbour *pnbr;
    struct sclumps *pclump, *pnewclump;
    struct sclumppu *pnewclumppu, *pcpu;


    R = (int *) calloc(puno,sizeof(int)); /* needed for separation */
    for (i=0;i<puno;i++)
        R[i] = pu[i].status;
    /*memcpy(R,pustat,sizeof(struct spustuff)*puno);*/

    /*** Step 1. Make a link list of all the possible PUs to be included ****/
    /*** This might change if I change the species v site into link lists ****/
    plisthead = NULL;
    for (i=0;i<puno;i++)
        if (rtnAmountSpecAtPu(pu,SM,i,isp) > 0)
        {
           if (pu[i].status==3) continue; /* not allowed to consider this one */
           if (pu[i].status==2)
           { /* add to clumps and remove from list */
              AddNewPU(i,isp,connections,spec,pu,SM,clumptype);
              continue;
           } /* checking if PU forced into reserve */
           iputotal++;
           plist = (struct slink *) malloc(sizeof(struct slink));
           plist->id = i;
           plist->next = plisthead;  /* Insert on list */
           plisthead = plist;  /* point head to new number */
        } /** Made link list of all sites with this species **/

    /* Check first to see if I've already satisfied targets for this species */
    SpeciesAmounts4(isp,spec,clumptype);
    if (spec[isp].sepnum>0)
       spec[isp].separation = CountSeparation2(isp,0,0,puno,R,pu,SM,spec,0);
    if ((spec[isp].amount >= spec[isp].target) && (spec[isp].occurrence >= spec[isp].targetocc) && (spec[isp].separation >= spec[isp].sepnum))
    {
       spec[isp].amount = 0;
       spec[isp].occurrence = 0;
       spec[isp].separation = 0;
       /** Clean out all the clump numbers for this species.*/
       while (spec[isp].head)
       {
             ClearClump(isp,spec[isp].head,pu,SM);
             pclump = spec[isp].head;
             spec[isp].head = spec[isp].head->next;
             free(pclump);
             DebugFree(sizeof(struct sclumps));
             spec[isp].clumps = 0;
       }  /** Remove each clump ***/
       free(R); /* dummy array for separation */
       DebugFree(puno * sizeof(int));
       return(-1);
    }  /* return when all targets already met. */

    if (iputotal)
    do
    {  /*** take all pu's at random until satisfied or I've run out **/
      /* Pluck a PU out at random */
      ipu = RandNum(iputotal);
      plist = plisthead;
      for (;ipu>0;ipu--)
      {
          plist = plist->next;
      }
      iputotal--;

      /** Add this PU to our system **/
      R[plist->id] = 1;
      AddNewPU(plist->id,isp,connections,spec,pu,SM,clumptype);

      /** Remove the chosen site from my site list **/
      if (plisthead == plist)
      {
         plisthead = plist->next;
      } /* special case for head of list */
      else
      {
          for (pdiscard = plisthead; pdiscard->next != plist; pdiscard = pdiscard->next)
              ; /*** Find link before plist ***/
          pdiscard->next = plist->next;
      } /* remove plist from the list */
      free(plist);
      DebugFree(sizeof(struct slink));

      /*** Check to see if I should continue by calculating current holdings **/
      SpeciesAmounts4(isp,spec,clumptype);
      if (spec[isp].sepnum>0)
         spec[isp].separation = CountSeparation2(isp,0,0,puno,R,pu,SM,spec,0);
    } while ((spec[isp].amount < spec[isp].target || spec[isp].separation < spec[isp].sepnum || spec[isp].occurrence < spec[isp].targetocc)
             && iputotal >= 1  );

    if (spec[isp].amount < spec[isp].target || spec[isp].occurrence < spec[isp].targetocc)
    {
       badspecies = 1;
       displayProgress1("Species %i cannot be fully represented!\n",spec[isp].name);
    } /*** Record the fact that the species is unrepresentable ***/
    if (spec[isp].separation < spec[isp].sepnum && spec[isp].amount >= spec[isp].target && spec[isp].occurrence >= spec[isp].targetocc)
    {
       badspecies = 1;
       displayProgress1("Species %i can only get %i separate valid clumps where %i are wanted!\n",
                        spec[isp].name,spec[isp].separation,spec[isp].sepnum);
    } /*** Record the fact that the species is unrepresentable ***/


    /* Search through the clumps looking for any which can be removed */
    /* But only do this if occurrence target met. Otherwise every single pu is neccessary*/
    if (spec[isp].occurrence >= spec[isp].targetocc)
    {
       pclump = spec[isp].head;
       while (pclump)
       {
             i = 0; /* if i becomes and stays 1 then this clump is removable */
             if (RemClumpCheck(pclump,pu))
                i = 1;
             if (i)
             {
                if (spec[isp].occurrence - pclump->occs >= spec[isp].targetocc)
                   i = 1;  /* if pclump-amount < target2 is caught in next step */
                else
                    i = 0;
             } /* Check is occurrence decrease ok? */
             if (i)
             {
                if ((spec[isp].amount - pclump->amount >= spec[isp].target) || (pclump->amount < spec[isp].target2))
                   i = 1;
                else
                    i = 0;
             } /* Check is amount decrease OK? */
             if (i && spec[isp].sepnum)
             {
                j = CountSeparation2(isp,0,pclump,puno,R,pu,SM,spec,-1);
                if ((j < spec[isp].separation) && (j < spec[isp].sepnum))
                   i = 0;
                else
                    i = 1;
                if (!spec[isp].target2)
                   i = 0; /* cannot elegantly remove clumps if species is listed as non-clumping */
             }
             if (i)   /* This is a clump which can be safely removed */
             {/* cut clump if uneccessary or it is too small */
                if (spec[isp].head == pclump)
                {
                   spec[isp].head = pclump->next;
                }
                else
                {
                    for (pnewclump = spec[isp].head;pnewclump->next != pclump;pnewclump = pnewclump->next)
                        ; /** find clump before pclump **/
                    pnewclump->next = pclump->next;
                }
                while (pclump->head)
                {
                      pnewclumppu = pclump->head;
                      pclump->head = pnewclumppu->next;
                      setClumpSpecAtPu(pu,SM,pnewclumppu->puid,isp,0);
                      free(pnewclumppu);
                      DebugFree(sizeof(struct sclumppu));
                }
                totalamount -= pclump->amount;
                /* cut out clump and progress pclump*/
                pnewclump = pclump;
                pclump = pclump->next;
                free(pnewclump);
                DebugFree(sizeof(struct sclumps));
                spec[isp].clumps--;
             } /** removing unneccessary pclump **/
             else
                 pclump = pclump->next;
       }
    } /*** Remove unneccesary clumps and links****/

    /** Test all PU's to see if any one of them are superfluous **/
    /* But only do this if occurrence target met. Otherwise every single pu is neccessary*/
    if (spec[isp].occurrence >= spec[isp].targetocc)
    {
       pclump = spec[isp].head;
       while (pclump)
       {
             pcpu = pclump->head;
             while (pcpu)
             {     /** Test to see if this pcpu is necessary **/
                   i = 0;
                   if (R[pcpu->puid] != 2)
                      i = 1;
                   if (i)
                   {
                      rAmount = rtnAmountSpecAtPu(pu,SM,pcpu->puid,isp);
                      if ((pclump->amount - rAmount > spec[isp].target2) && (spec[isp].amount - rAmount > spec[isp].target))
                         i = 1;
                      else
                          i = 0;
                   }  /* doesn't drop amount below min clump size or target */
                   if (i)
                   {
                      if (spec[isp].occurrence > spec[isp].targetocc)
                         i = 1;
                      else
                          i = 0;
                   } /* Does it drop occurrences too low? */
                   if (i)
                   {
                      pnewclump = (struct sclumps *)malloc(sizeof(struct sclumps));
                      pnewclump->clumpid = pcpu->puid;  /* sclump used to store clumpPU info */
                      pnewclump->amount = 0;
                      pnewclump->next = NULL;
                      j = CountSeparation2(isp,pcpu->puid,pnewclump,puno,R,pu,SM,spec,-1);
                      free(pnewclump);
                      if ((j < spec[isp].separation) && (j < spec[isp].sepnum))
                         i = 0;
                      else
                          i = 1;
                   } /* How does count sep fare? */
                   if (i)
                   {
                      if (ClumpCut(isp,pu,spec,pclump,pcpu,connections,SM,&dummy,&idummy,&j,0,clumptype))
                         i = 0;
                      else
                          i = 1;
                   } /* Does it cut the clump? these are not allowed to remove */
                   /* Theoretically they could possible be removed */
                   if (i)  /* Is this removable? */
                   {        /* remove pcpu */
                      setClumpSpecAtPu(pu,SM,pcpu->puid,isp,0);
                      totalamount -= rAmount;
                      pclump->amount -= rAmount;
                      if (pcpu == pclump->head)
                      {
                         pclump->head = pcpu->next;
                         free(pcpu);
                         DebugFree(sizeof(struct sclumppu));
                         pcpu = pclump->head;
                      } /* removing first clump */
                      else
                      {
                          for (pnewclumppu = pclump->head;pnewclumppu->next != pcpu;pnewclumppu = pnewclumppu->next)
                              ; /* find previous pcpu */
                          pnewclumppu->next = pcpu->next;
                          free(pcpu);
                          DebugFree(sizeof(struct sclumppu));
                          pcpu = pnewclumppu->next;
                      } /* removing pcpu when it is not the head */
                   }  /** remove unneccessary clumppu **/
                   else
                       pcpu = pcpu->next; /* moving pointer when it is not removable */
             } /* Checking each pcpu in clump */
             pclump = pclump->next;
       }
    } /** Cycle over each pclump **/


    while (plisthead)
    {
          plist = plisthead;
          plisthead = plisthead->next;
          free(plist);
          DebugFree(sizeof(struct slink));
    } /* Cleaing link list */


    /*** Now count the cost of this particular reserve ****/
    /*** For each clump figure out connection cost ***/
    pclump = spec[isp].head;
    while (pclump)
    {
          iclumpno = pclump->clumpid;
          pcpu = pclump->head;
          while (pcpu)
          {
                if (pu[pcpu->puid].status != 2)
                {
                   cost += pu[pcpu->puid].cost;
                   connection += connections[pcpu->puid].fixedcost;
                } /* only count fixed costs if PU not forced into reserve */
                if (connections[pcpu->puid].nbrno)
                {
                   pnbr = connections[pcpu->puid].first;
                   while (pnbr)
                   {
                         if (rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp) != iclumpno)
                            connection += pnbr->cost;
                         pnbr = pnbr->next;
                   } /** Counting each individual connection **/
                } /** Counting connection strength if neccessary **/
                pcpu = pcpu->next;
          } /** Checking each PU in clump **/
          pclump = pclump->next;
    } /*** Count cost for each clump ***/

    /* Finally. Calculate penalty from all of this.*/
    spec[isp].penalty = cost + connection *cm;

    /* Consider case where targets cannot be met */
    totalamount = 0;
    if (spec[isp].amount < spec[isp].target)
       totalamount = spec[isp].target / spec[isp].amount;
    if (spec[isp].occurrence < spec[isp].targetocc)
       totalamount += (float) spec[isp].targetocc/(float) spec[isp].occurrence;
    if (totalamount)
       spec[isp].penalty *= totalamount;  /* Scale it up */

    if (spec[isp].sepdistance)
       spec[isp].separation = 1;
    spec[isp].amount = 0; /* because my routines add it in */
    spec[isp].occurrence = 0;
    /** Clean out all the clump numbers for this species.*/
    while (spec[isp].head)
    {
          ClearClump(isp,spec[isp].head,pu,SM);
          pclump = spec[isp].head;
          spec[isp].head = spec[isp].head->next;
          free(pclump);
          DebugFree(sizeof(struct sclumps));
          spec[isp].clumps = 0;
    }  /** Remove each clump ***/

    free(R); /* dummy array for separation */
    DebugFree(puno * sizeof(int));
    return(badspecies);

} /*** Calculate Penalty for a Type 4 Species ***/

/**** Partial Penalty for type 4 species ***/
double PartialPen4(int isp, double amount,struct sspecies spec[],int clumptype)
{
    if (amount >= spec[isp].target2)
        return (amount);    /* I'm not a partial penalty */
    else
        switch(clumptype)
        {
            case 0:
                return(0.0); /* default step function */
            case 1:
                return(amount/ 2.0); /* nicer step function */
            case 2:
                if (spec[isp].target2)
                   return (amount/spec[isp].target2 * amount);
            default:
                return(0.0);
        }
}  /* Partial Penalty for type 4 species */

/*** Value for Adding a Planning Unit ****/
double ValueAdd(int isp,int ipu,int puno, int R[],
    struct sconnections connections[],struct spustuff pu[],struct spu SM[],struct sspecies spec[],int clumptype)
{  int ineighbours = 0,iclumpid,iseparation;
    struct sneighbour *pnbr;
    struct slink {int clumpid;double amount;
                int occs; struct slink *next;} *head = NULL,*plink;
    struct sclumps *pclump,*sepclump=NULL,*psclump;
    struct sclumppu *ppu;
    double amount,oldamount = 0.0,shortfall;
    int oldoccs = 0,occs, iClump;

        /* Count neighbours */
        if (connections[ipu].nbrno > 0) {
            pnbr = connections[ipu].first;
            while (pnbr) {
                iClump = rtnClumpSpecAtPu(pu,SM,pnbr->nbr,isp);
                if (iClump) {
                    iclumpid = 1;
                    /* Is nbr on my list ?*/
                    for(plink = head;plink;plink=plink->next)
                        if (plink->clumpid == iClump)
                            iclumpid = 0;

                    if (iclumpid){
                        ineighbours++;
                        plink = (struct slink *) malloc(sizeof(struct slink));
                        plink->clumpid = iClump;
                        /* find amount for this clump */
                        for(pclump = spec[isp].head;plink->clumpid != pclump->clumpid;
                                pclump = pclump->next)
                        ; /* find the right clump */
                        plink->amount = pclump->amount;
                        plink->occs = pclump->occs;
                        plink->next = head;
                        head = plink;
                        if (spec[isp].sepnum)
                         for (ppu = pclump->head;ppu;ppu=ppu->next)
                            { psclump = (struct sclumps *) malloc(sizeof(struct sclumps));
                              psclump->clumpid = ppu->puid;
                              psclump->next = sepclump;
                              sepclump = psclump;  /* glue to sep list. Still need amount */
                            } /* stick whole clump onto separation clump for later */
                    } /* new neighbour found */
                } /* neighbour of clump */
                pnbr = pnbr->next;
            } /** count all neighbours if they have a clump **/
        }  /* If There are neighbours */

        if (spec[isp].sepnum) {
          psclump = (struct sclumps *) malloc(sizeof(struct sclumps));
          psclump->clumpid = ipu;
          psclump->next = sepclump;
          sepclump = psclump;
        } /* Add ipu to my sepclump list */

        /* now I know number and names of neighbouring clumps */
        amount = rtnAmountSpecAtPu(pu,SM,ipu,isp);
        occs = (amount > 0);
        for(plink = head;plink;plink = plink->next) {
            amount += plink->amount;
            occs += plink->occs;
            oldamount += PartialPen4(isp,plink->amount,spec,clumptype);
            oldoccs += plink->occs * (PartialPen4(isp,plink->amount,spec,clumptype)>0);
        }

        /* set the sepclump amounts to this new amount */
        if (spec[isp].sepnum)
          for (psclump = sepclump;psclump;psclump = psclump->next)
            psclump->amount = amount;

        amount = PartialPen4(isp,amount,spec,clumptype);
        occs = occs * (amount > 0);

        amount = amount - oldamount; /* amount is change in amount for this species */
        occs = occs - oldoccs;

        if (spec[isp].sepnum) {
           iseparation = CountSeparation2(isp,0,sepclump,puno,R,pu,SM,spec,1);  /* imode = 1 doesn't do anything*/
          while (sepclump) {
            psclump = sepclump;
            sepclump = sepclump->next;
            free(psclump);
            DebugFree(sizeof(struct sclumps));
          }} /* clean up sepcount link list */

        while(head) {
            plink = head;
            head = head->next;
            free(plink);
            DebugFree(sizeof(struct slink));
        }  /* Clean up link list */

        /* Return the effective amount for this species */
        /* Old amount + change in amount + separation penalty for changed structure */

        amount = spec[isp].amount + amount;
        shortfall = 0;
        if (spec[isp].target)
            shortfall = amount >= spec[isp].target ? 0 : (spec[isp].target - amount)/spec[isp].target;
       if (spec[isp].targetocc)  {
            occs = occs + spec[isp].occurrence;
             amount = occs >= spec[isp].targetocc ? 0:
                    ((double)spec[isp].targetocc - (double) occs)/(double)spec[isp].targetocc;
            shortfall += amount;
                    }
        if (spec[isp].target && spec[isp].targetocc)
            shortfall /= 2;
        return(shortfall + SepPenalty(iseparation,spec[isp].sepnum));
} /*** Value for Adding a Planning Unit ****/


/** Value Remove. The amount of species loss for removing a single pu */
double ValueRem(int ipu,int isp,
    struct sspecies spec[],struct sconnections connections[],struct spustuff pu[],struct spu SM[],int clumptype)
{
    double newamount = 0,amount,shortfall=0;
    struct sclumps *pclump;
    struct sclumppu *ppu;
    int iseparation;
    int newocc = 0;

    /* locate the clump and clumppu of the target site ipu */
    for (pclump = spec[isp].head; pclump && pclump->clumpid != rtnClumpSpecAtPu(pu,SM,ipu,isp); pclump = pclump->next)
    ; /* locate correct clump list */

    for (ppu = pclump->head;ppu->puid != ipu; ppu = ppu->next)
    ; /* locate the correct pclump pu */

    /* debugmem2 = debugmem; debugging line */
    if (spec[isp].sepnum)
        ClumpCut(isp,pu,spec,pclump,ppu,connections,SM,&newamount,&newocc,&iseparation,1,clumptype);
    else
        ClumpCut(isp,pu,spec,pclump,ppu,connections,SM,&newamount,&newocc,&iseparation,0,clumptype);

  if (spec[isp].target)
      {
      amount = spec[isp].amount + newamount -PartialPen4(isp,pclump->amount,spec,clumptype) ;
      shortfall = amount > spec[isp].target ? 0 : (spec[isp].target - amount)/spec[isp].target;
      }  /* there is an abundance amount */

  /*if (isp == 16) printf("pclump->occs %i targetocc %i shortfall %.2f\n",
                                  pclump->occs,spec[isp].targetocc,shortfall);*/
  if (spec[isp].targetocc) {  /* Handle the case where there is a targetocc */
     amount = spec[isp].occurrence +newocc - pclump->occs * (PartialPen4(isp,pclump->amount,spec,clumptype)>0);
     if (amount < spec[isp].targetocc)
         shortfall += ((double) spec[isp].targetocc - amount)/(double) spec[isp].targetocc;
     if (spec[isp].target)
         shortfall /= 2;
    }
/*  if (isp ==16) printf("shortfall %.2f occ %i newocc %i pclump->amount %.2f\n",
          shortfall, spec[isp].occurrence,newocc,pclump->amount);*/

  return(shortfall + SepPenalty(iseparation,spec[isp].sepnum));
} /** Value for removing a planning unit ****/


/***************   NewPenalty4   *********************/
/* Calculates the new penalty for adding or removing a PU for species which have
    clumping requirements */


double NewPenalty4(int ipu,int isp,int puno,struct sspecies spec[],struct spustuff pu[],struct spu SM[],
                    int R[],struct sconnections connections[],int imode,int clumptype)
{    double amount;

        if (imode == 1) {
            if (spec[isp].penalty == 0)
                return (0);  /* Targets have all already been met */
            amount = ValueAdd(isp,ipu,puno,R,connections,pu,SM,spec,clumptype);
        }
        else {
              /* determine change in this amount */
            amount = ValueRem(ipu,isp,spec,connections,pu,SM,clumptype);
        } /** removing a planning unit **/
        return(amount);

}  /*** The new penalty for type 4 species ***/

// SlaveExit does not deliver a message prior to exiting, but creates a file so C-Plan knows marxan has exited
void SlaveExit(void)
{
     writeSlaveSyncFile();
}

struct snlink *GetVarName(char **varlist,int numvars,char *sVarName,
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
          for (temp = head;temp;temp = temp->next)
          {
              if (strcmp(temp->name,sVarName) == 0)
                 displayErrorMessage("ERROR: variable %s has been defined twice in data file %s.\n",sVarName,fname);
          }

       newlink = (struct snlink *) malloc(sizeof(struct snlink));
       newlink->next = NULL;
       newlink->name = (char *) calloc(strlen(sVarName)+1,sizeof(char));
       strcpy(newlink->name,sVarName);
       return(newlink);
}

int CheckVarName(char **varlist, int numvars, char *sVarName)
{
    // This routine checks if the variable name occurs in the list. It is similar to GetVarName but does not create list
    int i,foundit = 0;

    for (i=0;i<numvars;++i)
        if(strcmp(varlist[i],sVarName) == 0)
               foundit++;

    return(foundit);
}

void rdsvar(FILE *infile, char varname[], void *address, int parmtype, int crit,int present)
// Reads a variable of parmtype in from infile. Assumes that the next line is the one that has the
// variable in question but will wrap once to find the variable.
//
// I changed this function because it was ignoring the last line of the file (also it was a dogs breakfast).
// it returns the "Final value taken" error message when SAVESOLUTIONSMATRIX parameter is on last line, weird, need to fix this.
{
     int foundit, namelen, check1, check2, gotit;
     char buffer[255] = "\0";    /* for storing the line found in the file */
     namelen = strlen(varname);    /* figure out how long the string is */
     foundit = 0;

     rewind(infile); /* Always search from top of infile */
      /* read first line. I'm in trouble if file is empty*/
     do
     {   /* loop through file looking for varname */
       fgets(buffer,255,infile);
       check1 = 0;
       check2 = 0;

       while (buffer[check1++] == varname[check2++])
             ;

       if (check1 > namelen)
       {  // varname matches upto namelen
          foundit++;
          switch (parmtype)
          {
                 case REAL : gotit = sscanf(&buffer[check1]," %f", (float *) address);
                             break;
                 case DOUBLE : gotit = sscanf(&buffer[check1]," %lf", (double *) address);
                               break;
                 case INTEGER : gotit = sscanf(&buffer[check1]," %d", (int *) address);
                                break;
                 case LONGINT : gotit = sscanf(&buffer[check1]," %ld", (long int *) address);
                                break;
                 case STRING : // trim leading and trailing blanks (allow spaces which are important for directory names
                               check1 += strspn(&buffer[check1]," ,");
                               for (check2 = strlen(&buffer[check1])-1;isspace(buffer[check1+check2]) != 0;check2--)
                                   ; // Find last non space character
                               if (strlen(&buffer[check1]) <2)
                                  buffer[check1] = '\0';
                               buffer[check1 + check2+1] = '\0';

                               strcpy((char *) address,&buffer[check1]);

                               gotit = 1; // So that var check works. This needs further consideration
                               break;
                 default : displayErrorMessage("Invalid parameter type request %d: \n",parmtype);
          }

          if (!gotit)
          {
             displayWarningMessage("WARNING: found bad value for variable %s. Value ignored\n",varname);
             foundit--;
          }
       }
     } while (!(feof(infile)));

     if (!foundit)
        if (crit)
           displayErrorMessage("Unable to find %s in input file.\n",varname);

     if (foundit > 1)
        displayWarningMessage("WARNING variable: %s appears more than once in the input file. Final value taken\n",varname);

     present = foundit;

     return;
}

/********** Set Options ************/
//  PUno is the number of Planning units
//  SPno is the number of species
/* This file uses the freeform format method. Expecting a file with a variable name
 (generally all capital letters) followed by the value for that variable
 */
void SetOptions(double *cm,double *prop,struct sanneal *anneal,
                int *iseed,
                long int *repeats,char savename[],struct sfname *fnames,char filename[],
                int *runopts,double *misslevel,int *heurotype,int *clumptype,
                int *itimptype, int *verb,
                double *costthresh,double *tpf1,double *tpf2)
{
     FILE *fp;
     double version;
     int present;
     char stemp[500];
     #ifdef DEBUGTRACEFILE
     char debugbuffer[200];
     #endif

     verbosity = 1; /* This enables local warning messages */
     /* Setup all of the default parameter variables */
     version = 0.1;
     *cm = 0;
     *prop = 0;
     (*anneal).type = 1;
     (*anneal).iterations = 0;
     (*anneal).Tinit = 1;
     (*anneal).Tcool = 0;
     (*anneal).Titns = 1;
     *iseed = -1;
     *costthresh = 0;
     *tpf1 = 0;
     *tpf2 = 0;
     *repeats = 0;
     (*fnames).saverun = 0;
     (*fnames).savebest = 0;
     (*fnames).savesum = 0;
     (*fnames).savesen = 0;
     (*fnames).savespecies = 0;
     (*fnames).savesumsoln = 0;
     (*fnames).savepenalty = 0;
     (*fnames).savetotalareas = 0;
     (*fnames).savedebugtracefile = 0;
     (*fnames).saverichness = 0;
     (*fnames).savesolutionsmatrix = 0;
     (*fnames).solutionsmatrixheaders = 1;
     (*fnames).savelog = 0;
     (*fnames).saveannealingtrace = 0;
     (*fnames).annealingtracerows = 0;
     (*fnames).saveitimptrace = 0;
     (*fnames).itimptracerows = 0;
     (*fnames).savespec = 0;
     (*fnames).savepu = 0;
     (*fnames).savepuvspr = 0;
     (*fnames).savematrixsporder = 0;
     (*fnames).rimagetype = 0;
     (*fnames).rexecutescript = 0;
     (*fnames).rclustercount = 0;
     (*fnames).rimagewidth = 0;
     (*fnames).rimageheight = 0;
     (*fnames).rimagefontsize = 0;
     asymmetricconnectivity = 0;


     strcpy(savename,"temp");
     *misslevel = 1;
     *heurotype = 1;
     *clumptype = 0;
     verbosity = 1;

     /* Open file and then feed in each variable type */
     if ((fp = fopen(filename,"r"))==NULL)
        displayErrorMessage("input file %s not found\nAborting Program.\n\n",filename);

     rdsvar(fp,"VERSION",&version,DOUBLE,0,present);
     rdsvar(fp,"BLM",cm,DOUBLE,0,present);
     if (present == 0)
        rdsvar(fp,"CM",cm,DOUBLE,0,present);
     rdsvar(fp,"PROP",prop,DOUBLE,0,present);
     rdsvar(fp,"RANDSEED",iseed,INTEGER,0,present);

     /* Annealing Controls */
     rdsvar(fp,"NUMITNS",&(*anneal).iterations,LONGINT,0,present);
     rdsvar(fp,"STARTTEMP",&(*anneal).Tinit,DOUBLE,0,present);
     rdsvar(fp,"COOLFAC",&(*anneal).Tcool,DOUBLE,0,present);
     rdsvar(fp,"NUMTEMP",&(*anneal).Titns,INTEGER,0,present);

     (*anneal).type = 1;
     if ((*anneal).iterations < 1 )
        (*anneal).type = 0;
     if ((*anneal).Tinit < 0)
        (*anneal).type = (int) (-(*anneal).Tinit) + 1;  /* type is negative of Tinit */
     fscanf(fp,"%i",iseed); /* The random seed. -1 to set by clock */

     /* Various controls */
     rdsvar(fp,"NUMREPS",repeats,LONGINT,0,present);
     rdsvar(fp,"COSTTHRESH",costthresh,DOUBLE,0,present);
     rdsvar(fp,"THRESHPEN1",tpf1,DOUBLE,0,present);
     rdsvar(fp,"THRESHPEN2",tpf2,DOUBLE,0,present);

     /* SaveFiles */
     rdsvar(fp,"SCENNAME",savename,STRING,0,present);

     /* SaveFiles New Method */
     rdsvar(fp,"SAVERUN",&(*fnames).saverun,INTEGER,0,present);
     rdsvar(fp,"SAVEBEST",&(*fnames).savebest,INTEGER,0,present);
     rdsvar(fp,"SAVESUMMARY",&(*fnames).savesum,INTEGER,0,present);
     rdsvar(fp,"SAVESCEN",&(*fnames).savesen,INTEGER,0,present);
     rdsvar(fp,"SAVETARGMET",&(*fnames).savespecies,INTEGER,0,present);
     rdsvar(fp,"SAVESUMSOLN",&(*fnames).savesumsoln,INTEGER,0,present);
     rdsvar(fp,"SAVESPEC",&(*fnames).savespec,INTEGER,0,present);
     rdsvar(fp,"SAVEPU",&(*fnames).savepu,INTEGER,0,present);
     rdsvar(fp,"SAVEMATRIXPUORDER",&(*fnames).savepuvspr,INTEGER,0,present);
     rdsvar(fp,"SAVEMATRIXSPORDER",&(*fnames).savematrixsporder,INTEGER,0,present);
     rdsvar(fp,"SAVEPENALTY",&(*fnames).savepenalty,INTEGER,0,present);
     rdsvar(fp,"SAVETOTALAREAS",&(*fnames).savetotalareas,INTEGER,0,present);
     rdsvar(fp,"SAVEDEBUGTRACEFILE",&(*fnames).savedebugtracefile,INTEGER,0,present);
     rdsvar(fp,"SAVERICHNESS",&(*fnames).saverichness,INTEGER,0,present);
     rdsvar(fp,"SAVESOLUTIONSMATRIX",&(*fnames).savesolutionsmatrix,INTEGER,0,present);
     rdsvar(fp,"SOLUTIONSMATRIXHEADERS",&(*fnames).solutionsmatrixheaders,INTEGER,0,present);
     rdsvar(fp,"SAVELOG",&(*fnames).savelog,INTEGER,0,present);
     rdsvar(fp,"SAVESNAPSTEPS",&(*fnames).savesnapsteps,INTEGER,0,present);
     rdsvar(fp,"SAVESNAPCHANGES",&(*fnames).savesnapchanges,INTEGER,0,present);
     rdsvar(fp,"SAVESNAPFREQUENCY",&(*fnames).savesnapfrequency,INTEGER,0,present);
     rdsvar(fp,"SAVEANNEALINGTRACE",&(*fnames).saveannealingtrace,INTEGER,0,present);
     rdsvar(fp,"ANNEALINGTRACEROWS",&(*fnames).annealingtracerows,INTEGER,0,present);
     rdsvar(fp,"SAVEITIMPTRACE",&(*fnames).saveitimptrace,INTEGER,0,present);
     rdsvar(fp,"ITIMPTRACEROWS",&(*fnames).itimptracerows,INTEGER,0,present);
     rdsvar(fp,"RIMAGETYPE",&(*fnames).rimagetype,INTEGER,0,present);
     rdsvar(fp,"REXECUTESCRIPT",&(*fnames).rexecutescript,INTEGER,0,present);
     rdsvar(fp,"RCLUSTERCOUNT",&(*fnames).rclustercount,INTEGER,0,present);
     rdsvar(fp,"RIMAGEWIDTH",&(*fnames).rimagewidth,INTEGER,0,present);
     rdsvar(fp,"RIMAGEHEIGHT",&(*fnames).rimageheight,INTEGER,0,present);
     rdsvar(fp,"RIMAGEFONTSIZE",&(*fnames).rimagefontsize,INTEGER,0,present);
     rdsvar(fp,"ASYMMETRICCONNECTIVITY",&asymmetricconnectivity,INTEGER,0,present);
     rdsvar(fp,"CONNECTIVITYIN",&fOptimiseConnectivityIn,INTEGER,0,present);

     // quantum annealing control parameters
     rdsvar(fp,"QAPROP",&rQAPROP,DOUBLE,0,present);
     rdsvar(fp,"QADECAY",&rQADECAY,DOUBLE,0,present);
     rdsvar(fp,"QADECAYB",&rQADECAYB,DOUBLE,0,present);
     rdsvar(fp,"QADECAYTYPE",&iQADECAYTYPE,INTEGER,0,present);
     rdsvar(fp,"QAACCPR",&rQAACCPR,DOUBLE,0,present);

     if (!(*fnames).savesnapfrequency)
        (*fnames).savesnapfrequency = 1;

     /* Filenames */
     rdsvar(fp,"INPUTDIR",stemp,STRING,1,present);
     if (stemp[strlen(stemp)-1] != '/' && stemp[strlen(stemp)-1] != '\\')
        strcat(stemp,"/");
     (*fnames).inputdir = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).inputdir,stemp);

     rdsvar(fp,"OUTPUTDIR",stemp,STRING,1,present);
     if (stemp[strlen(stemp)-1] != '/' && stemp[strlen(stemp)-1] != '\\')
        strcat(stemp,"/");
     (*fnames).outputdir = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).outputdir,stemp);

     strcpy(stemp,"PU.dat");
     rdsvar(fp,"PUNAME",stemp,STRING,1,present);
     (*fnames).puname = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).puname,stemp);

     strcpy(stemp,"spec.dat");
     rdsvar(fp,"SPECNAME",stemp,STRING,1,present);
     //if (present == 0)
     //   rdsvar(fp,"FEATNAME",stemp,STRING,0,present);
     (*fnames).specname = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).specname,stemp);

     strcpy(stemp,"puvspr2.dat");
     rdsvar(fp,"PUVSPRNAME",stemp,STRING,1,present);
     (*fnames).puvsprname = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).puvsprname,stemp);

     //strcpy(stemp,"NULL");
     //rdsvar(fp,"PUVSPPROBNAME",stemp,STRING,0,present);
     //(*fnames).puvspprobname = (char *) calloc(strlen(stemp)+1,sizeof(char));
     //strcpy((*fnames).puvspprobname,stemp);

     strcpy(stemp,"NULL");
     rdsvar(fp,"MATRIXSPORDERNAME",stemp,STRING,0,present);
     (*fnames).matrixspordername = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).matrixspordername,stemp);

     strcpy(stemp,"NULL");
     rdsvar(fp,"PENALTYNAME",stemp,STRING,0,present);
     (*fnames).penaltyname = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).penaltyname,stemp);

     strcpy(stemp,"NULL");
     rdsvar(fp,"BOUNDNAME",stemp,STRING,0,present);
     if (present == 0)
        rdsvar(fp,"CONNECTIONNAME",stemp,STRING,0,present);
     (*fnames).connectionname = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).connectionname,stemp);

     strcpy(stemp,"NULL");
     rdsvar(fp,"CONNECTIONFILESNAME",stemp,STRING,0,present);
     (*fnames).connectionfilesname = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).connectionfilesname,stemp);

     strcpy(stemp,"NULL");
     rdsvar(fp,"BLOCKDEFNAME",stemp,STRING,0,present);
     (*fnames).blockdefname = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).blockdefname,stemp);

     strcpy(stemp,"SOLUTION");
     rdsvar(fp,"BESTFIELDNAME",stemp,STRING,0,present);
     (*fnames).bestfieldname = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).bestfieldname,stemp);


     strcpy(stemp,"NULL");
     rdsvar(fp,"RBINARYPATHNAME",stemp,STRING,0,present);
     (*fnames).rbinarypath = (char *) calloc(strlen(stemp)+1,sizeof(char));
     strcpy((*fnames).rbinarypath,stemp);

     /* various other controls */
     rdsvar(fp,"RUNMODE",runopts,INTEGER,1,present);
     rdsvar(fp,"MISSLEVEL",misslevel,DOUBLE,0,present);
     rdsvar(fp,"HEURTYPE",heurotype,INTEGER,0,present);
     rdsvar(fp,"CLUMPTYPE",clumptype,INTEGER,0,present);
     rdsvar(fp,"ITIMPTYPE",itimptype,INTEGER,0,present);
     rdsvar(fp,"VERBOSITY",verb,INTEGER,0,present);
     verbosity = *verb;
     rdsvar(fp,"PROBABILITYWEIGHTING",stemp,STRING,0,present);
     sscanf(stemp, "%lf", &rProbabilityWeighting);

     rdsvar(fp,"STARTDECTHRESH",&rStartDecThresh,DOUBLE,0,present);
     rdsvar(fp,"ENDDECTHRESH",&rEndDecThresh,DOUBLE,0,present);
     rdsvar(fp,"STARTDECMULT",&rStartDecMult,DOUBLE,0,present);
     rdsvar(fp,"ENDDECMULT",&rEndDecMult,DOUBLE,0,present);

     #ifdef DEBUGTRACEFILE
     sprintf(debugbuffer,"PROBABILITYWEIGHTING %g\n",rProbabilityWeighting);
     appendTraceFile(debugbuffer);
     sprintf(debugbuffer,"STARTDECTHRESH %g\n",rStartDecThresh);
     appendTraceFile(debugbuffer);
     sprintf(debugbuffer,"ENDDECTHRESH %g\n",rEndDecThresh);
     appendTraceFile(debugbuffer);
     sprintf(debugbuffer,"STARTDECMULT %g\n",rStartDecMult);
     appendTraceFile(debugbuffer);
     sprintf(debugbuffer,"ENDDECMULT %g\n",rEndDecMult);
     appendTraceFile(debugbuffer);
     #endif

     if ((*fnames).outputdir[0] != '0')
     {
        strcpy(stemp,(*fnames).outputdir);
        strcat(stemp,savename);
        strcpy(savename,stemp);
     }
     fclose(fp);

}  /***** Set Options *******/

void PrepareWeightedConnectivityFile(struct sfname fnames)
{
    char *readname1, *readname2, *writename, *sFileName,*sWeighting, *sId1, *sId2, *sConnection, *sAsymmetric;
    FILE *fpnames, *fpInputConnection, *fpOutputConnection;
    char sLine[500];
    int iRecords, iTotalRecords, iFiles;
    double rWeighting, rConnection;
    #ifdef DEBUGTRACEFILE
    char debugbuffer[200];
    #endif

    // prepare file names for backing up the connection file
    readname1 = (char *) calloc(strlen(fnames.connectionname) + strlen(fnames.inputdir)+2, sizeof(char));
    strcpy(readname1,fnames.inputdir);
    strcat(readname1,fnames.connectionname);
    writename = (char *) calloc(strlen(fnames.connectionname) + strlen(fnames.inputdir)+3, sizeof(char));
    strcpy(writename,fnames.inputdir);
    strcat(writename,fnames.connectionname);
    strcat(writename,"~");
    // back up the connection file
    copyFile(readname1,writename);
    free(writename);
    free(readname1);

    // prepare connection file for output
    writename = (char *) calloc(strlen(fnames.connectionname) + strlen(fnames.inputdir)+2, sizeof(char));
    strcpy(writename,fnames.inputdir);
    strcat(writename,fnames.connectionname);
    if ((fpOutputConnection = fopen(writename,"w"))==NULL)
    {
       displayProgress1("Warning: Connection File %s cannot be written ",fnames.connectionname);
       free(writename);
    }
    free(writename);
    fprintf(fpOutputConnection,"id1,id2,connection\n");

    // prepare connection file names file for input
    readname1 = (char *) calloc(strlen(fnames.connectionfilesname) + strlen(fnames.inputdir)+2, sizeof(char));
    strcpy(readname1,fnames.inputdir);
    strcat(readname1,fnames.connectionfilesname);
    if ((fpnames = fopen(readname1,"r"))==NULL)
    {
       displayProgress1("Warning: Connectivity Files Name File %s not found ",fnames.connectionfilesname);
       free(readname1);
    }
    free(readname1);
    fgets(sLine,500-1,fpnames);

    // loop through the connectivity files
    iTotalRecords = 0;
    iFiles = 0;
    while (fgets(sLine,500-1,fpnames))
    {
          iFiles++;
          iRecords = 0;

          // read the 3 fields from the line
          sFileName = strtok(sLine," ,;:^*\"/\t\'\\\n");
          sWeighting = strtok(NULL," ,;:^*\"/\t\'\\\n");
          sscanf(sWeighting,"%lf",&rWeighting);
          sAsymmetric = strtok(NULL," ,;:^*\"/\t\'\\\n");

          if (rWeighting != 0)
          {
             // prepare current connectivity file for input
             readname2 = (char *) calloc(strlen(sFileName) + strlen(fnames.inputdir)+2, sizeof(char));
             strcpy(readname2,fnames.inputdir);
             strcat(readname2,sFileName);
             if ((fpInputConnection = fopen(readname2,"r"))==NULL)
             {
                displayProgress1("Warning: Input Connectivity  File %s not found ",sFileName);
                free(readname2);
             }
             free(readname2);
             // read the header row
             fgets(sLine,500-1,fpInputConnection);

             // write the input connection file contents to the output connection file with appropriate weightings
             while (fgets(sLine,500-1,fpInputConnection))
             {
                   iRecords++;
                   iTotalRecords++;

                   sId1 = strtok(sLine," ,;:^*\"/\t\'\\\n");
                   sId2 = strtok(NULL," ,;:^*\"/\t\'\\\n");
                   sConnection = strtok(NULL," ,;:^*\"/\t\'\\\n");
                   sscanf(sConnection,"%lf",&rConnection);

                   fprintf(fpOutputConnection,"%s,%s,%lf\n",sId1,sId2,(rConnection*rWeighting));

                   if (strcmp("yes",sAsymmetric) == 0)
                      fprintf(fpOutputConnection,"%s,%s,%lf\n",sId2,sId1,(rConnection*rWeighting));
             }

             fclose(fpInputConnection);
          }

          #ifdef DEBUGTRACEFILE
          sprintf(debugbuffer,"connectivity file %s weighting %lf asymmetric >%s< records %i\n",
                              sFileName,rWeighting,sAsymmetric,iRecords);
          appendTraceFile(debugbuffer);
          #endif

          free(sFileName);
          free(sWeighting);
          free(sAsymmetric);
    }

    #ifdef DEBUGTRACEFILE
    sprintf(debugbuffer,"total files %i records %i\n",iFiles,iTotalRecords);
    appendTraceFile(debugbuffer);
    #endif

    fclose(fpOutputConnection);
    fclose(fpnames);

    asymmetricconnectivity = 1;
}


void MapUserPenalties(typesp spec[],int spno)
{
     int i;

     for (i=0;i<spno;i++)
         spec[i].penalty = spec[i].rUserPenalty;
}

/********************************************************/
/************* Greedy Add On ****************************/
/********************************************************/

// Greedy Species Penalty
double GreedyPen(int ipu,int puno, int spno, typesp spec[],int R[],struct spustuff pu[],struct spu SM[],
                 int clumptype)
{
       int i;
       double famount = 0.0, fold,newamount;
       
       for (i = 0;i<spno;i++)
       {
           fold = (spec[i].target - spec[i].amount);
           if (fold > 0)
           {
              if (spec[i].target2)
                 newamount = NewPenalty4(ipu,i,puno,spec,pu,SM,R,connections,1,clumptype);
              else
                  newamount = NewPenalty(ipu,i,spec,pu,SM,1);
                  
              famount += (newamount - fold)*spec[i].spf;
           } // Add new penalty if species isn't already in the system
       }
       
       return(famount);  // Negative means decrease in amount missing
} // Greedy Species Penalty

/********* Greedy Score an alternative to the normal objective function *****/
double GreedyScore(int ipu,int puno, int spno, typesp *spec,struct spu SM[],
                   struct sconnections connections[],int R[],struct spustuff pu[],double cm,int clumptype)
{
       double currpen,currcost,currscore;
       
       currpen = GreedyPen(ipu,puno,spno,spec,R,pu,SM,clumptype);
       currcost = pu[ipu].cost + ConnectionCost2(ipu,connections,R,1,1,cm);
       if (currcost <= 0)
       {
          currscore = -1.0/delta;
       } // otherwise this 'free pu' will have +score
       else
       {
           currscore = currpen/currcost;
           // multiply by rand (1.000,1.001)
       }
       return(currscore);
} // Score for a planning unit based upon greedy algorithm

/*********** Rarity Settup. Sets up rare score for each species ******/
/**** score is total abundance / smallest species abundance *********/
void SetRareness(int puno, int spno, double Rare[],int R[],struct spustuff pu[],struct spu SM[])
{
     double smallest = 0;
     double *fcount;
     int i, ism, isp,ipu;

     #ifdef DEBUG_HEURISTICS
     appendTraceFile("SetRareness start\n");
     #endif

     fcount = (double *) calloc(spno,sizeof(double));

     for (isp=0;isp<spno;isp++)
         fcount[isp] = 0;

     for (ipu=0;ipu<puno;ipu++)
         if (pu[ipu].richness)
            for (i=0;i<pu[ipu].richness;i++)
            {
                ism = pu[ipu].offset + i;
                isp = SM[ism].spindex;
                if (R[ipu] < 2)
                   fcount[isp] += SM[ism].amount;
            }

     for (isp=0;isp<spno;isp++)
     {
         if (smallest == 0 || (fcount[isp] < smallest && fcount[isp] > 0))
            smallest = fcount[isp];
         Rare[isp] = fcount[isp];

         #ifdef DEBUG_HEURISTICS
         appendTraceFile("SetRareness isp %i Rare %g\n",isp,Rare[isp]);
         #endif
     }

     if (smallest == 0)
        displayErrorMessage("Serious Error in calculating Rarenesses. No species detected.\n");

     for (isp=0;isp<spno;isp++)
         Rare[isp] /= smallest;

     free(fcount);

     #ifdef DEBUG_HEURISTICS
     appendTraceFile("SetRareness end\n");
     #endif
}  // SetRareness

// RareScore The score for a particular conservation value on a particular PU
double RareScore(int isp,int ipu,int puno,typesp spec[],struct spu SM[], int R[],
                 struct sconnections connections[],struct spustuff pu[], double cm,int clumptype)
{
       double currpen,currcost,currscore;
       double fold, newamount;
       fold = (spec[isp].target - spec[isp].amount);
       if (fold > 0)
       {
          if (spec[isp].target2)
             newamount = NewPenalty4(ipu,isp,puno,spec,pu,SM,R,connections,1,clumptype);
          else
              newamount = NewPenalty(ipu,isp,spec,pu,SM,1);
          currpen = newamount - fold;
       } // Add new penalty if species isn't already in the system

       currcost = pu[ipu].cost + ConnectionCost2(ipu,connections,R,1,1,cm);
       if (currcost <= 0)
       {
          currscore = -1.0/delta;
       } // otherwise this 'free pu' will have +score
       else
       {
           currscore = currpen/currcost;
           // multiply by rand (1.000,1.001)
       }

       return(currscore);
} // RareScore

// Max Rare Score Heuristic. PU scores based on rarest beast on PU
double MaxRareScore(int ipu,int puno,struct sspecies spec[],struct spu SM[],
                    int R[],struct sconnections connections[],struct spustuff pu[],double cm,double Rare[],int clumptype)
{
       int i, ism, isp,rareno = -1;
       double rarest,rarescore;

       if (pu[ipu].richness)
          for (i=0;i<pu[ipu].richness;i++)
          {
              ism = pu[ipu].offset + i;
              isp = SM[ism].spindex;
              if (SM[ism].amount && (spec[isp].target > spec[isp].amount || (spec[isp].sepdistance && spec[isp].separation < 3)))
                 if (1.0/Rare[isp]< rarest || rareno < 0)
                 {
                    rareno = isp;
                    rarest = Rare[isp];
                 }  // Determine which is the rarest species
          }

       if (rareno > -1)
          rarescore = RareScore(rareno,ipu,puno,spec,SM,R,connections,pu,cm,clumptype)/rarest;
       else
           rarescore = 1.0 / delta;

       return(rarescore);
} // Max Rare Score

// Best Rarity Score. Determines each species rare score
double BestRareScore(int ipu,int puno,struct sspecies spec[],struct spu SM[],
                     int R[],struct sconnections connections[],struct spustuff pu[],double cm,double Rare[],int clumptype)
{
       int i, ism, isp,rareno = -1;
       double rarest = 0,rarescore;

       if (pu[ipu].richness)
          for (i=0;i<pu[ipu].richness;i++)
          {
              ism = pu[ipu].offset + i;
              isp = SM[ism].spindex;
              if (SM[ism].amount && (spec[isp].target > spec[isp].amount || (spec[isp].sepdistance && spec[isp].separation < 3)))
              {
                 rarescore = RareScore(isp,ipu,puno,spec,SM,R,connections,pu,cm,clumptype)/Rare[isp];
                 if (rarescore > rarest || rareno < 0)
                 {
                    rarest = rarescore;
                    rareno = isp;
                 }
              }
          }

       return(rarescore);

} // Best Rare Score

// Average Rare Score. Rare Score for each scoring species/number scoring species
double AveRareScore(int ipu,int puno,struct sspecies spec[],struct spu SM[],
                    int R[],struct sconnections connections[],struct spustuff pu[],double cm,double Rare[],int clumptype)
{
       int i, ism, isp, rareno = 0;
       double rarescore = 0;

       if (pu[ipu].richness)
          for (i=0;i<pu[ipu].richness;i++)
          {
              ism = pu[ipu].offset + i;
              isp = SM[ism].spindex;
              if (SM[ism].amount &&
                  (spec[isp].target > spec[isp].amount ||
                  (spec[isp].sepdistance && spec[isp].separation < 3)))
              {
                 rarescore += RareScore(isp,ipu,puno,spec,SM,R,connections,pu,cm,clumptype)/Rare[isp];
                 rareno++;
              }
          }

       return(rarescore/rareno);
} // Ave Rare Score

// Sum of Rare Score for each scoring species
double SumRareScore(int ipu,int puno,struct sspecies spec[],struct spu SM[],
                    int R[],struct sconnections connections[],struct spustuff pu[],
                    double cm,double Rare[],int clumptype)
{
       int i, ism, isp;
       double rarescore = 0;

       #ifdef DEBUG_HEURISTICS
       appendTraceFile("SumRareScore start\n");
       #endif

       if (pu[ipu].richness)
          for (i=0;i<pu[ipu].richness;i++)
          {
              #ifdef DEBUG_HEURISTICS
              appendTraceFile("SumRareScore feature %i of %i\n",i,pu[ipu].richness);
              #endif

              ism = pu[ipu].offset + i;
              isp = SM[ism].spindex;

              #ifdef DEBUG_HEURISTICS
              appendTraceFile("SumRareScore SMamount %g target %g specamount %g rare %g\n",SM[ism].amount,spec[isp].target,spec[isp].amount,Rare[isp]);
              #endif

              if (SM[ism].amount &&
                  (spec[isp].target > spec[isp].amount ||
                  (spec[isp].sepdistance && spec[isp].separation < 3)))
                 rarescore += RareScore(isp,ipu,puno,spec,SM,R,connections,pu,cm,clumptype)/Rare[isp];
          }

       #ifdef DEBUG_HEURISTICS
       appendTraceFile("SumRareScore end\n");
       #endif

       return(rarescore);
} // Sum Rare Score

// Set Abundances
void SetAbundance(int puno,double Rare[],struct spustuff pu[],struct spu SM[])
{  
     int i,j, ism, isp;

     for (i=0;i<puno;i++)
         if (pu[i].richness)
            for (j=0;j<pu[i].richness;j++)
            {
                ism = pu[i].offset + j;
                isp = SM[ism].spindex;
                Rare[isp] += SM[ism].amount;
            }
} // Set Abundance

// Irreplaceability For site for species
double Irreplaceability(int ipu,int isp, double Rare[],struct spustuff pu[],struct spu SM[],typesp *spec)
{
       double buffer,effamount;
       
       buffer = Rare[isp] < spec[isp].target ? 0 : Rare[isp] - spec[isp].target;
       if (spec[isp].amount > spec[isp].target)
          return(0);
       effamount = rtnAmountSpecAtPu(pu,SM,ipu,isp);
       
       return(buffer<effamount ? 1 : effamount/buffer);
}

// Product Irreplaceability for a single site
double ProdIrr(int ipu,double Rare[],struct spustuff pu[],struct spu SM[],typesp *spec)
{
       int i, ism, isp;
       double product = 1;

       if (pu[ipu].richness)
          for (i=0;i<pu[ipu].richness;i++)
          {
              ism = pu[ipu].offset + i;
              isp = SM[ism].spindex;
              if (SM[ism].amount && (spec[isp].target - spec[isp].amount)> 0)
                 product *= (1-Irreplaceability(ipu,isp,Rare,pu,SM,spec));
          }

       return(1-product);
} // Product Irreplaceability

// Sum Irreplaceability for a single site
double SumIrr(int ipu,double Rare[],struct spustuff pu[],struct spu SM[],typesp *spec)
{
       int i, ism, isp;
       double sum = 0;

       if (pu[ipu].richness)
          for (i=0;i<pu[ipu].richness;i++)
          {
              ism = pu[ipu].offset + i;
              isp = SM[ism].spindex;
              if (SM[ism].amount && (spec[isp].target - spec[isp].amount)> 0)
                 sum += (Irreplaceability(ipu,isp,Rare,pu,SM,spec));
          }

       return(sum);
} // Sum Irreplaceability

// Main Heuristic Engine
void Heuristics(int spno,int puno,struct spustuff pu[],struct sconnections connections[],
                int R[], double cm,typesp *spec,struct spu SM[], struct scost *reserve,
                double costthresh, double tpf1,double tpf2, int imode,int clumptype)
// imode = 1: 2: 3: 4:
// imode = 5: 6: Prod Irreplaceability, 7: Sum Irreplaceability
{
     int i,bestpu;
     double bestscore,currscore;
     struct scost change;
     double *Rare;

     #ifdef DEBUG_HEURISTICS
     appendTraceFile("Heuristics start\n");
     #endif

     // Irreplacability

     if (imode >= 6 && imode <= 7)
     {
        Rare = (double *) calloc(spno,sizeof(double));
        SetAbundance(puno,Rare,pu,SM);

        #ifdef DEBUG_HEURISTICS
        appendTraceFile("Heuristics 6 or 7\n");
        #endif
     }

     if (imode >= 2 && imode <= 5) // Rareness Setups
     {
        Rare = (double *) calloc(spno,sizeof(double));
        SetRareness(puno,spno,Rare,R,pu,SM);

        #ifdef DEBUG_HEURISTICS
        appendTraceFile("Heuristics 2 to 5 after SetRareness\n");
        #endif
     }

     do
     {
       bestpu = 0;
       bestscore = 0;
       for (i=0;i<puno;i++)
           if (!R[i]) // Only look for new PUS
           {
              // Set the score for the given Planning Unit
              currscore = 1; // null if no other mode set
              if (imode == 0)
                 currscore = GreedyScore(i,puno,spno,spec,SM,connections,R,pu,cm,clumptype);
              if (imode == 1)
              {
                 CheckChange(-1,i,spno,puno,pu,connections,spec,SM,R,cm,1,&change,reserve,
                             costthresh,tpf1,tpf2,1, clumptype);
                 currscore = change.total;
              }
              if (imode == 2)
              {
                 currscore = MaxRareScore(i,puno,spec,SM,R,connections,pu,cm,Rare, clumptype);
              }
              if (imode == 3)
              {
                 currscore = BestRareScore(i,puno,spec,SM,R,connections,pu,cm,Rare,clumptype);
              }
              if (imode == 4)
              {
                 currscore = AveRareScore(i,puno,spec,SM,R,connections,pu,cm,Rare,clumptype);
              }
              if (imode == 5)
              {

                 #ifdef DEBUG_HEURISTICS
                 appendTraceFile("Heuristics pu %i of %i\n",i,puno);
                 #endif

                 currscore = SumRareScore(i,puno,spec,SM,R,connections,pu,cm,Rare,clumptype);

                 #ifdef DEBUG_HEURISTICS
                 appendTraceFile("Heuristics after SumRareScore\n");
                 #endif
              }
              if (imode == 6)
              {
                 currscore = -ProdIrr(i,Rare,pu,SM,spec);
              }
              if (imode == 7)
              {
                 currscore = -SumIrr(i,Rare,pu,SM,spec);
              }

              currscore *=(double) rand1()*0.001 + 1.0;
              if (!costthresh || pu[i].cost + reserve->cost <= costthresh)
                 if (currscore < bestscore)
                 {
                    bestpu = i;
                    bestscore = currscore;
                 } // is this better (ie negative) than bestscore?
           } // I've looked through each pu to find best
           if (bestscore)
           {
              CheckChange(-1,bestpu,spno,puno,pu,connections,spec,SM,R,cm,1,&change,reserve,
                          costthresh,tpf1,tpf2,1,clumptype);
              DoChange(bestpu,puno,R,reserve,change,pu,SM,spec,connections,1,clumptype);

              // Different Heuristics might have different penalty effects
              // Old fashioned penalty and missing counting
              reserve->missing = 0;
              for (i=0;i<spno;i++)
              {
                  if (spec[i].amount < spec[i].target)
                  {
                     reserve->missing++;
                  }
                  else
                      if (spec[i].sepdistance && spec[i].separation < 3)
                         reserve->missing++;
                  // Species missing
              } // checking to see who I am missing
           } // Add Pu as long as I've found one
           
           if (bestscore)
           {
              displayProgress2("P.U. %i PUs %i score %.6f Cost %.1f Connection %.1f Missing %i",
                              bestpu,reserve->pus,bestscore,reserve->cost,reserve->connection,reserve->missing);
              if (fProb1D == 1)
                 displayProgress2(" Probability1D %.1f\n",reserve->probability1D);
              if (fProb2D == 1)
                 displayProgress2(" Probability2D %.1f\n",reserve->probability2D);
              displayProgress2(" Penalty %.1f\n",reserve->penalty);
           }

     } while (bestscore); // Repeat until all good PUs have been added
     
     reserve->total = reserve->cost + reserve->connection + reserve->penalty + reserve->probability1D + reserve->probability2D;

     #ifdef DEBUG_HEURISTICS
     appendTraceFile("Heuristics end\n");
     #endif
} // Heuristics

/************ Iterative Improvement *************/

void siftDown_ii(struct iimp numbers[], int root, int bottom, int array_size)
{
     int done, maxChild;
     typeiimp temp;

     done = 0;
     while ((root*2 <= bottom) && (!done))
     {
           if (root*2 < array_size)
           {
              if (root*2 == bottom)
                 maxChild = root * 2;
              else if (numbers[root * 2].randomfloat > numbers[root * 2 + 1].randomfloat)
                      maxChild = root * 2;
                   else
                       maxChild = root * 2 + 1;

              if (numbers[root].randomfloat < numbers[maxChild].randomfloat)
              {
                 temp = numbers[root];
                 numbers[root] = numbers[maxChild];
                 numbers[maxChild] = temp;
                 root = maxChild;
              }
              else
                  done = 1;
           }
           else
               done = 1;
     }
}

void heapSort_ii(struct iimp numbers[], int array_size)
{
     int i;
     typeiimp temp;
     #ifdef DEBUGTRACEFILE
     //char debugbuffer[80];
     #endif

     for (i = (array_size / 2)-1; i >= 0; i--)
     {
         #ifdef DEBUGTRACEFILE
         //sprintf(debugbuffer,"heapSort_ii i %i\n",i);
         //appendTraceFile(debugbuffer);
         #endif

         siftDown_ii(numbers, i, array_size, array_size);
     }

     for (i = array_size-1; i >= 1; i--)
     {
         #ifdef DEBUGTRACEFILE
         //sprintf(debugbuffer,"heapSort_ii i %i\n",i);
         //appendTraceFile(debugbuffer);
         #endif

         temp = numbers[0];
         numbers[0] = numbers[i];
         numbers[i] = temp;
         siftDown_ii(numbers, 0, i-1, array_size);
     }
}

void IterativeImprovement(int puno,int spno,struct spustuff pu[], struct sconnections connections[],
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

     appendTraceFile("IterativeImprovement start\n");

     // counting pu's we need to test
     for (i=0;i<puno;i++)
         if ((R[i] < 2) && (pu[i].status < 2))
            puvalid++;

     appendTraceFile("IterativeImprovement puvalid %i\n",puvalid);

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
            if ((R[i] < 2) && (pu[i].status < 2))
            {
               iimparray[ipu].puindex = i;
               iimparray[ipu].randomfloat = rand1();
               ipu++;
            }

        appendTraceFile("IterativeImprovement after array init\n");

        // sort the iimp array by the randomindex field
        heapSort_ii(iimparray,puvalid);

        appendTraceFile("IterativeImprovement after heapSort_ii\n");

        /***** Doing the improvements ****/
        for (i=0;i<puvalid;i++)
        {
            ichoice = iimparray[i].puindex;

            if ((R[ichoice] < 2) && (pu[ichoice].status < 2))
            {
               imode = R[ichoice] == 1 ? -1 : 1;
               CheckChange(-1,ichoice,spno,puno,pu,connections,spec,SM,R,cm,imode,change,reserve,
                           costthresh,tpf1,tpf2,1,clumptype);
               if (change->total < 0)
               {
                  displayProgress2("It Imp has changed %i with change value %lf \n",ichoice,change->total);
                  DoChange(ichoice,puno,R,reserve,*change,pu,SM,spec,connections,imode,clumptype);
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

     appendTraceFile("IterativeImprovement end\n");
}  // Time Optimised Iterative Improvement

/*********************
 *  ran1() from numerical recipes
     produces a random number between 0 and 1
 */

#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1 + (IM - 1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long    RandomIY;
long    RandomIV[NTAB];
float rand1(void)
{
    int        j;
    long    k;
    float   temp;

    if(RandSeed1 <= 0 || !RandomIY)    /* Initialize */
    {
    RandSeed1 = -RandSeed1;
    for(j = NTAB+7; j >= 0; j--)
    {
        k = RandSeed1/IQ;
        RandSeed1 = IA * (RandSeed1 - k * IQ) - IR * k;
        if(RandSeed1 < 0) RandSeed1 += IM;
        if(j < NTAB) RandomIV[j] = RandSeed1;
    }
    RandomIY=RandomIV[0];
    }
    k=RandSeed1/IQ;        /* The stuff we do on calls after the first */
    RandSeed1 = IA * (RandSeed1 - k * IQ) - IR * k;
    if(RandSeed1 < 0) RandSeed1 += IM;
    j = RandomIY/NDIV;
    RandomIY=RandomIV[j];
    RandomIV[j] = RandSeed1;
    if((temp=AM*RandomIY) > RNMX) return(RNMX);
    else return(temp);
}

void InitRandSeed(int iSeed) {
    if (iSeed>0)
        RandSeed1 = iSeed;
    else
        RandSeed1 = (long int)time(NULL);
    if(RandSeed1 > 0) RandSeed1 = -RandSeed1;

        //return RandSeed1;
}

/* Returns a random number between 0 and num - 1, where num is an int */
int RandNum (int num)
{
    long    temp;

    if(num == 0) return(0);
    temp = (int)(rand1() * num);
    if (temp == num) return(0);
    else return((int)temp);
}

/*****************************************************/
/********* Separation Measure Routines ***************/
/*****************************************************/

// Sep Penalty
// This returns the penalty for not meeting separation requirments. Feed in sepnum and current
//  separation and returns a value from 0 to 1 which is an artificial shortfall.
double SepPenalty(int ival,int itarget)
{
       double fval;

       if (!itarget)
          return (0); /* no penalty if no separation requirement*/
          
       fval = (double) ival / (double) itarget;
       if (!ival)
          fval = 1.0 /(double) itarget;

       return (1/(7*fval+0.2)-(1/7.2)); // Gives a nice hyperbole with fval = 1 return 0 and fval = 0 or 0.1 returning almost 1
} // SepPenalty

int ValidPU(int ipu,int isp,struct sclumps *newno,struct sspecies spec[],struct spustuff pu[],
            struct spu SM[],int imode)
{
    // Returns true if ipu is acceptable as a planning unit
    int i;
    
    i = rtnIdxSpecAtPu(pu,SM,ipu,isp);
    struct sclumps *pclump, *ppu;
    if (newno)
    {
       if (imode == -2)
          if (SM[i].clump == newno->clumpid)
             return(0); // This whole clump is to be removed
             
       for (ppu=newno;ppu;ppu=ppu->next)
       {
           if (ipu == ppu->clumpid)
           {
              if (ppu->amount < spec[isp].target2)
                 return(0);
              else
                  return(1);
           }
       } // ipu is on list of changed pus
    } // Above used only when newno is not NULL
    
    // Find clump
    for (pclump = spec[isp].head;pclump && (SM[i].clump != pclump->clumpid);pclump= pclump->next)
        ; // scan through to find clump
        
    if (pclump)
    {
       if (pclump->amount <spec[isp].target2)
          return(0);
       else
           return(1);
    }
    else
    {
        if (SM[i].amount < spec[isp].target2)
           return(0);
        else
            return(1);
    }
} // Valid PU

int CheckDistance(int i, int j,struct spustuff pu[],double squaretarget)
{
    // compare x1*x2+y1*y2 with squaretarget
    if ((pu[i].xloc-pu[j].xloc)*(pu[i].xloc-pu[j].xloc) + (pu[i].yloc-pu[j].yloc)* (pu[i].yloc-pu[j].yloc) >= squaretarget)
        return(1);
    else
        return(0);
} // Is Distant returns true if PU's are big enough distance apart

int CountSeparation(int isp,struct sclumps *newno,
                struct spustuff pu[],struct spu SM[],typesp spec[],int imode)
{
    // imode 0 = count separation on current
    // imode 1 = count separation if ipu were included
    // imode -1 = count separation if ipu were excluded
    // The following assumes imode = 0 for starters

    struct slink{int id; struct slink *next;} *first = NULL, *second = NULL,*ptemp,*ptest;
    struct sclumps *pclump;
    struct sclumppu *ppu;
    int sepcount = 1,test;
    double targetdist;
  
    targetdist = spec[isp].sepdistance * spec[isp].sepdistance;

    if (targetdist == 0)
       return(3); // Shortcut if sep not apply to this species. This assumes that 3 is highest sep count
       
    // Set up the first list
    if (imode == 1)
       if (ValidPU(newno->clumpid,isp,newno,spec,pu,SM,imode))
       {
          ptemp = (struct slink *) malloc(sizeof(struct slink));
          ptemp->id = newno->clumpid;
          ptemp->next = first;
          first = ptemp;
       }
    for (pclump = spec[isp].head;pclump;pclump = pclump->next)
        for (ppu = pclump->head;ppu;ppu= ppu->next)
            if (ValidPU(ppu->puid,isp,newno,spec,pu,SM,imode))
            {
               ptemp = (struct slink *) malloc(sizeof(struct slink));
               ptemp->id = ppu->puid;
               ptemp->next = first;
               first = ptemp;
            } // Add all valid species bearing PU's to list
    // need to worry about added pu which is not on spec[isp].head list

    // cycle through this list
    while (first)
    {
          test = first->id;
          ptemp = first;
          first = first->next;
          free(ptemp);
          DebugFree(sizeof(struct slink));

          for (ptemp = first;ptemp;ptemp = ptemp->next)
              if (CheckDistance(ptemp->id,test,pu,targetdist))
              {
                 for (ptest=second;ptest;ptest = ptest->next)
                     if (CheckDistance(ptemp->id,ptest->id,pu,targetdist))
                     {
                        // Clear all lists
                        while (first)
                        {
                              ptemp = first;
                              first = ptemp->next;
                              free(ptemp);
                              DebugFree(sizeof(struct slink));
                        }
                        while (second)
                        {
                              ptemp = second;
                              second = ptemp->next;
                              free(ptemp);
                              DebugFree(sizeof(struct slink));
                        }
                        return(3);
                     } // I have succeeded in finding what I'm looking for

                 ptest = (struct slink *) malloc(sizeof(struct slink));
                 ptest->id = ptemp->id;
                 ptest->next = second;
                 second = ptest;
                 sepcount = 2; // there is a separation of at least 2. This should be used to cut down calls to this function

              } // I am sufficient distance from test location

              while (second)
              {
                    ptemp = second;
                    second = ptemp->next;
                    free(ptemp);
                    DebugFree(sizeof(struct slink));
              } // clear second between tests
    } // finished scanning through list. first is neccessarily empty now
    
    while (second)
    {
          ptemp = second;
          second = ptemp->next;
          free(ptemp);
          DebugFree(sizeof(struct slink));
    }

    return(sepcount);
} // CountSeparation

// Make List
// This makes a list of all the valid PUs which occur on the reserve and on which species
// isp is present (or NULL), in the form of a slink link list

struct slink *makelist(int isp,int ipu,int puno,int R[],
                       struct sclumps *newno,struct sspecies spec[],
                       struct spustuff pu[],struct spu SM[],int imode)
       // imode: 0 : as is. -1 ipu being removed, +1 ipu being added
{
       struct sclumps *pclump;
       struct sclumppu *ppu;
       struct slink *ptemp,*head=NULL;
       int i;
       double rAmount = rtnAmountSpecAtPu(pu,SM,ipu,isp);

       if (spec[isp].target2)
       {
          // deal with clumping species differently from non-clumping
          if (imode == 1)
             if (ValidPU(newno->clumpid,isp,newno,spec,pu,SM,imode))
             {
                ptemp = (struct slink *) malloc(sizeof(struct slink));
                ptemp->id = newno->clumpid;
                ptemp->next = head;
                head = ptemp;
             }
             
          for (pclump = spec[isp].head;pclump;pclump = pclump->next)
              for (ppu = pclump->head;ppu;ppu= ppu->next)
                  if (ValidPU(ppu->puid,isp,newno,spec,pu,SM,imode))
                  {
                     ptemp = (struct slink *) malloc(sizeof(struct slink));
                     ptemp->id = ppu->puid;
                     ptemp->next = head;
                     head = ptemp;
                  }  // Add all valid species bearing PU's to list
       } // if target2
       else
       {
           // non clumping species
           if ((imode ==1) && rAmount)
           {
              ptemp = (struct slink *)malloc(sizeof(struct slink));
              ptemp->id = ipu;
              ptemp->next = head;
              head = ptemp;
           } // deal with imode == 1 case
           
           for (i=0;i<puno;i++)
               if (((R[i] == 1 || R[i] == 2) && rAmount) && !(imode == -1 && ipu == i))
               {
                  ptemp = (struct slink *)malloc(sizeof(struct slink));
                  ptemp->id = i;
                  ptemp->next = head;
                  head = ptemp;
               }
       } // non clumping species

       return(head);
} // Makelist

// Sep Deal List
// This funciton is called by count separation2. It takes a link list of sites and 'deals' them
//  out on to the seplist
int SepDealList(struct slink *head, typeseplist *Dist,struct spustuff *pu,
                struct sspecies spec[],int first,int sepnum,double targetdist,
                int isp)
// Currsep is the current separation maximum it is 0 up to sepnum
// first is only needed if maximum is at 0, sepnum is the target separation
{  
    int placefound,currtarget,bestsep=0;
    int currsep;
    struct slink *temp;

    while (head)
    {
          placefound = 0;
          currtarget = first;
          currsep = sepnum;
          do
          {
            if (CheckDistance(head->id,currtarget,pu,targetdist))
            {
               currsep++;
               if (currsep == spec[isp].sepnum-1)
               {
                  while (head)
                  {
                        temp = head->next;
                        head->next = Dist[currsep].head;
                        Dist[currsep].head = head;
                        head = temp;
                  } // glue remaining list on to bottom of Dist. ignoring size and tail as useless
                
                  return(currsep);
               } // Just found valid separation
               
               if (Dist[currsep].head)
                  currtarget = Dist[currsep].head->id;
               else
               {
                   placefound = 1;
                   Dist[currsep].head = head;
                   Dist[currsep].tail = head;
                   Dist[currsep].size++;
                   head = head->next;
                   Dist[currsep].tail->next = NULL;
               } // I'm at the end of the line
            } // Good distance
            else
            {
                placefound = 1;
                Dist[currsep].tail->next = head;
                Dist[currsep].tail = head;
                Dist[currsep].size++;
                head = head->next;
                Dist[currsep].tail->next = NULL;
            } // bad distance
          } while (!placefound); // Doing each individual
          if (currsep > bestsep)
             bestsep = currsep;
    }
    
    return(bestsep);
}

/* This is a modified form of count separation where the user can specify any
    maximum separation distance rather than just assuming a sep distance of three */
/* ipu and newno used when imode <> 0. When counting as if ipu were added or removed
    ipu used for non-clumping and newno for clumping species */

int CountSeparation2(int isp,int ipu,struct sclumps *newno,int puno,int R[],
                     struct spustuff pu[],struct spu SM[],typesp spec[],int imode)
{   
    typeseplist *Dist;
    struct slink *head = NULL,*temp;
    int sepcount,bestsep = 0,i,currcol;
    double targetdist;
    
    targetdist = spec[isp].sepdistance * spec[isp].sepdistance;

    if (targetdist == 0)
       return(spec[isp].sepnum); // Shortcut if sep not apply to this species

    // Set up array for counting separation
    Dist = (typeseplist *) calloc(spec[isp].sepnum,sizeof(typeseplist));
    // First scan through sites. Grab first valid and place rest in lists
    head = makelist(isp,ipu,puno,R,newno,spec,pu,SM,imode);

    if (!head)
    {
       free(Dist);
       return(0);
    } // There was nothing to put in the list


    Dist[0].head = head;
    Dist[0].size = 1;
    Dist[0].tail = head;
    head = head->next;
    Dist[0].tail->next = NULL;
    if (!head)
    {
       free(Dist[0].head);
       free(Dist);
       return(1);
    }  // There was only one item in the list

    // Deal out link list
    sepcount = SepDealList(head,Dist,pu,spec,Dist[0].head->id,0,targetdist,isp);
    if (sepcount >= spec[isp].sepnum-1)
    {
       // clean up arrays
       for (i=0;i<spec[isp].sepnum;i++)
           while (Dist[i].head)
           {
                 temp = Dist[i].head;
                 Dist[i].head = Dist[i].head->next;
                 free(temp);
           }
           
       free(Dist);
       return(spec[isp].sepnum);
    } // I'm at maximum separation
    
    bestsep = sepcount;

    do
    {
      // The main Loop
      for (currcol=0;Dist[currcol+1].head && currcol < spec[isp].sepnum-2;currcol++)
          ;
          
      if (currcol == 0)
      {
         if (Dist[0].size < spec[isp].sepnum)
         {
            while (Dist[0].head)
            {
                  temp = Dist[0].head;
                  Dist[0].head = Dist[0].head->next;
                  free(temp);
            }
            free(Dist);
            return(bestsep + 1);
         } // cannot increase separation terminate function
         else
         {
             temp = Dist[0].head;
             Dist[0].head = Dist[0].head->next;
             head = Dist[0].head->next;
             Dist[0].head->next = NULL;
             Dist[0].size = 1;
             Dist[0].tail = Dist[0].head;
             free(temp);
             sepcount = SepDealList(head,Dist,pu,spec,Dist[0].head->id,0,targetdist,isp);
         }
      } // Deal with first column
      else
      {
          if (Dist[currcol].size + currcol  < spec[isp].sepnum)
          {
             Dist[currcol-1].tail->next = Dist[currcol].head;
             Dist[currcol-1].tail = Dist[currcol].tail;
             Dist[currcol-1].tail->next = NULL;
             Dist[currcol-1].size += Dist[currcol].size;
             Dist[currcol].head = NULL;
             Dist[currcol].size = 0;
             Dist[currcol].tail = NULL;
             sepcount = 0;
          } // list is not long enough to increase sepcount
          else
          {
              Dist[currcol-1].tail->next = Dist[currcol].head;
              Dist[currcol-1].tail = Dist[currcol].head;
              Dist[currcol-1].size++;
              Dist[currcol].head = Dist[currcol].head->next;
              head = Dist[currcol].head->next;
              Dist[currcol].head->next = NULL;
              Dist[currcol-1].tail->next = NULL;
              Dist[currcol].tail = Dist[currcol].head;
              Dist[currcol].size = 1;
              sepcount = SepDealList(head,Dist,pu,spec,
              Dist[currcol].head->id,currcol,targetdist,isp);
          } // else this column might be long enough
      } // Deal with columns other than the first
      
      if (sepcount > bestsep)
         bestsep = sepcount;
    } while(bestsep < spec[isp].sepnum-1); // Main loop.

    // clean up arrays
    for (i=0;i<spec[isp].sepnum;i++)
        while (Dist[i].head)
        {
              temp = Dist[i].head;
              Dist[i].head = Dist[i].head->next;
              free(temp);
        }
        
    free(Dist);
    return(bestsep+1);
} // CountSeparation 2

/* This function produces a usage statement for when the command line is
not correct. Usually when there is the wrong number (ie > 1) arguments
passed to the program */

void Usage(char *programName)
{
    fprintf(stderr,"%s usage: %s -[o] -[c] [input file name]\n",programName,programName);
}  // Usage


/* returns the index of the first argument that is not an option; i.e.
does not start with a dash or a slash
*/
void HandleOptions(int argc,char *argv[],char sInputFileName[])
{
     int i;

     if (argc>4)
     {  // if more than one commandline argument then exit
        Usage(argv[0]);
        exit(1);
     }

     for (i=1;i<argc;i++)
     {   // Deal with all arguments
         if (argv[i][0] == '/' || argv[i][0] == '-')
         {
            switch(argv[i][1])
            {
                        case 'C':
                        case 'c':
                        case 'S':
                        case 's':
                                marxanisslave = 1;
                                break;
            default:
                fprintf(stderr,"unknown option %s\n",argv[i]);
                break;
            }
         }
         else
             strcpy(sInputFileName,argv[i]); /* If not a -option then must be input.dat name */
     }
}

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
        HandleOptions(argc,argv,sInputFileName);

    if (Marxan(sInputFileName))        // Calls the main annealing unit
    {
        if (marxanisslave == 1)
        {
           SlaveExit();
        }

        return 1;
    }  // Abnormal Exit
    if (marxanisslave == 1)
    {
       SlaveExit();
    }

    return 0;
}

// use the prop value from the conservation feature file to set a proportion target for species
void ApplySpecProp(int spno,typesp spec[],int puno,struct spustuff pu[],struct spu SM[])
{
     // compute and set target for species with a prop value
     double totalamount;
     int isp, ipu;

     for (isp=0;isp<spno;isp++)
         if (spec[isp].prop > 0)
         {
            for (ipu = 0,totalamount = 0;ipu<puno;ipu++)
                totalamount += rtnAmountSpecAtPu(pu,SM,ipu,isp);
            spec[isp].target = totalamount * spec[isp].prop;
         }
}

// Change in probability 1D for adding or deleting one PU
double ChangeProbability1D(int iIteration, int ipu, int spno,int puno,struct sspecies spec[],struct spustuff pu[],struct spu SM[],int imode)
{
       int i, ism, isp, iNewHeavisideStepFunction, iOrigHeavisideStepFunction;
       double rSumAmount, rOrigExpected, rOrigVariance, rNewExpected, rNewVariance, rOriginalZScore;
       double rNewZScore, rProbOriginal, rProbNew, rSumProbability = 0;
       double rDeltaExpected, rDeltaVariance, rNewShortfallPenalty, rOrigShortfallPenalty;
       char debugbuffer[200];

       #ifdef DEBUG_PROB1D
       char sDebugFileName[300], sIteration[20];
       if (iIteration < 10)
       {
          sprintf(sIteration,"%i",iIteration);
          strcpy(sDebugFileName,fnames.outputdir);
          strcat(sDebugFileName,"output_ChangeProbability1DDebug_");
          strcat(sDebugFileName,sIteration);
          strcat(sDebugFileName,".csv");
          writeChangeProbability1DDebugTable(sDebugFileName,iIteration,ipu,spno,spec,pu,SM,imode);
       }
       #endif

       if (pu[ipu].richness)
          for (i=0;i<pu[ipu].richness;i++)
          {
              ism = pu[ipu].offset + i;
              isp = SM[ism].spindex;
              if (SM[ism].amount)
              {
                 // compute new probability
                 rDeltaExpected = imode * SM[ism].amount * (1 - pu[ipu].prob);
                 rNewExpected = spec[isp].expected1D + rDeltaExpected;

                 rDeltaVariance = imode * SM[ism].amount * SM[ism].amount * pu[ipu].prob * (1 - pu[ipu].prob);
                 rNewVariance = spec[isp].variance1D + rDeltaVariance;

                 if (rNewVariance > 0)
                    rNewZScore = (spec[isp].target - rNewExpected) / sqrt(rNewVariance);
                 else
                     rNewZScore = 4;

                 spec[isp].Zscore1D = rNewZScore;

                 if (rNewZScore >= 0)
                    rProbNew = probZUT(rNewZScore);
                 else
                     rProbNew = 1 - probZUT(-1 * rNewZScore);

                 spec[isp].probability1D = rProbNew;

                 // compute original probability
                 rSumAmount = spec[isp].amount;

                 if (spec[isp].variance1D > 0)
                    rOriginalZScore = (spec[isp].target - spec[isp].expected1D) / sqrt(spec[isp].variance1D);
                 else
                     rOriginalZScore = 4;

                 if (rOriginalZScore >= 0)
                    rProbOriginal = probZUT(rOriginalZScore);
                 else
                     rProbOriginal = 1 - probZUT(-1 * rOriginalZScore);

                 if (spec[i].ptarget1d > rProbNew)
                    iNewHeavisideStepFunction = 1;
                 else
                     iNewHeavisideStepFunction = 0;

                 if (spec[i].ptarget1d > rProbOriginal)
                    iOrigHeavisideStepFunction = 1;
                 else
                     iOrigHeavisideStepFunction = 0;

                 if (spec[i].ptarget1d > 0)
                 {
                    rNewShortfallPenalty = (spec[i].ptarget1d - rProbNew) / spec[i].ptarget1d;
                    rOrigShortfallPenalty = (spec[i].ptarget1d - rProbOriginal) / spec[i].ptarget1d;
                 }
                 else
                 {
                     rNewShortfallPenalty = 0;
                     rOrigShortfallPenalty = 0;
                 }

                 // change in probability
                 rSumProbability += (iNewHeavisideStepFunction * rNewShortfallPenalty) - (iOrigHeavisideStepFunction * rOrigShortfallPenalty);
              }
          }

       return (rSumProbability * rProbabilityWeighting);
}  // Change in probability for adding or deleting one PU

// Change in probability 2D for adding or deleting one PU
double ChangeProbability2D(int iIteration, int ipu, int spno,int puno,struct sspecies spec[],struct spustuff pu[],struct spu SM[],int imode)
{
    int i, ism, isp, iNewHeavisideStepFunction, iOrigHeavisideStepFunction;
    double rSumAmount, rOrigExpected, rOrigVariance, rNewExpected, rNewVariance, rOriginalZScore;
    double rProb, rNewZScore, rProbOriginal, rProbNew, rSumProbability = 0;
    double rDeltaExpected, rDeltaVariance, rNewShortfallPenalty, rOrigShortfallPenalty;

    #ifdef DEBUG_PROB2D
    char sDebugFileName[300], sIteration[20];
    if (iIteration < 10)
    {
       sprintf(sIteration,"%i",iIteration);
       strcpy(sDebugFileName,fnames.outputdir);
       strcat(sDebugFileName,"output_ChangeProbability2DDebug_");
       strcat(sDebugFileName,sIteration);
       strcat(sDebugFileName,"_.csv");
       writeChangeProbability2DDebugTable(sDebugFileName,iIteration,ipu,spno,spec,pu,SM,imode);
    }
    #endif

    if (pu[ipu].richness)
       for (i=0;i<pu[ipu].richness;i++)
       {
           ism = pu[ipu].offset + i;
           isp = SM[ism].spindex;
           if (SM[ism].amount)
           {
              rProb = SM[ism].prob;

              // compute new probability
              rDeltaExpected = imode * SM[ism].amount * rProb;
              rNewExpected = spec[isp].expected2D + rDeltaExpected;

              rDeltaVariance = imode * SM[ism].amount * SM[ism].amount * rProb * (1 - rProb);
              rNewVariance = spec[isp].variance2D + rDeltaVariance;

              if (rNewVariance > 0)
                 rNewZScore = (spec[isp].target - rNewExpected) / sqrt(rNewVariance);
              else
                  rNewZScore = 4;

              spec[isp].Zscore2D = rNewZScore;

              if (rNewZScore >= 0)
                 rProbNew = probZUT(rNewZScore);
              else
                  rProbNew = 1 - probZUT(-1 * rNewZScore);

              spec[isp].probability2D = rProbNew;

              // compute original probability
              rSumAmount = spec[isp].amount;

              if (spec[isp].variance2D > 0)
                 rOriginalZScore = (spec[isp].target - spec[isp].expected2D) / sqrt(spec[isp].variance2D);
              else
                  rOriginalZScore = 4;

              if (rOriginalZScore >= 0)
                 rProbOriginal = probZUT(rOriginalZScore);
              else
                  rProbOriginal = 1 - probZUT(-1 * rOriginalZScore);

              if (spec[i].ptarget2d > rProbNew)
                 iNewHeavisideStepFunction = 1;
              else
                  iNewHeavisideStepFunction = 0;

              if (spec[i].ptarget2d > rProbOriginal)
                 iOrigHeavisideStepFunction = 1;
              else
                  iOrigHeavisideStepFunction = 0;

              if (spec[i].ptarget2d > 0)
              {
                 rNewShortfallPenalty = (spec[i].ptarget2d - rProbNew) / spec[i].ptarget2d;
                 rOrigShortfallPenalty = (spec[i].ptarget2d - rProbOriginal) / spec[i].ptarget2d;
              }
              else
              {
                  rNewShortfallPenalty = 0;
                  rOrigShortfallPenalty = 0;
              }

              // change in probability
              rSumProbability += (iNewHeavisideStepFunction * rNewShortfallPenalty) - (iOrigHeavisideStepFunction * rOrigShortfallPenalty);
           }
        }

    return (rSumProbability * rProbabilityWeighting);
}  // Change in probability for adding or deleting one PU

// accumulate ExpectedAmount and VarianceInExpectedAmount for each species at this planning unit
void ReturnProbabilityAmounts1D(double ExpectedAmount1D[],double VarianceInExpectedAmount1D[],int ipu,
                                int puno,struct spustuff pu[],struct spu SM[])
{
     int i, ism, isp;

     if (pu[ipu].richness)
        for (i=0;i<pu[ipu].richness;i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            if (SM[ism].amount)
            {
               ExpectedAmount1D[isp] += SM[ism].amount * (1 - pu[ipu].prob);

               VarianceInExpectedAmount1D[isp] += SM[ism].amount * SM[ism].amount * pu[ipu].prob * (1 - pu[ipu].prob);
            }
        }
}
void ReturnProbabilityAmounts2D(double ExpectedAmount2D[],double VarianceInExpectedAmount2D[],int ipu,
                                int puno,struct spustuff pu[],struct spu SM[])
{
     int i, ism, isp;

     if (pu[ipu].richness)
        for (i=0;i<pu[ipu].richness;i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            if (SM[ism].amount)
            {
               ExpectedAmount2D[isp] += SM[ism].amount * SM[ism].prob;

               VarianceInExpectedAmount2D[isp] += SM[ism].amount * SM[ism].amount * SM[ism].prob * (1 - SM[ism].prob);
            }
        }
}

double ComputeProbability1D(double ExpectedAmount1D[],double VarianceInExpectedAmount1D[],
                            int spno,struct sspecies spec[])
{
       // compute Probability for all reserved planning units
       int i, iHeavisideStepFunction;
       double rProbability, rSumProbability = 0, rShortfallPenalty;

       appendTraceFile("ComputeProbability1D start\n");

       for (i=0;i<spno;i++)
       {
           if (VarianceInExpectedAmount1D[i] > 0)
              spec[i].Zscore1D = (spec[i].target - ExpectedAmount1D[i]) / sqrt(VarianceInExpectedAmount1D[i]);
           else
               spec[i].Zscore1D = 4;

           if (spec[i].Zscore1D >= 0)
              rProbability = probZUT(spec[i].Zscore1D);
           else
               rProbability = 1 - probZUT(-1 * spec[i].Zscore1D);

           if (spec[i].ptarget1d > rProbability)
              iHeavisideStepFunction = 1;
           else
               iHeavisideStepFunction = 0;

           if (spec[i].ptarget1d > 0)
              rShortfallPenalty = (spec[i].ptarget1d - rProbability) / spec[i].ptarget1d;
           else
               rShortfallPenalty = 0;

           spec[i].probability1D = iHeavisideStepFunction * rShortfallPenalty;

           rSumProbability += spec[i].probability1D;

           #ifdef DEBUG_PROB1D
           appendTraceFile("ComputeProbability1D i %i p %lf one minus p %lf pst %lf pst2 %lf\n",
                                i,rProbability,(1-rProbability),spec[i].probability1D,(1-spec[i].probability1D));
           #endif
       }

       appendTraceFile("ComputeProbability1D sump %lf\n",rSumProbability);

       return rSumProbability * rProbabilityWeighting;
}

double ComputeProbability2D(double ExpectedAmount2D[],double VarianceInExpectedAmount2D[],
                            int spno,struct sspecies spec[])
{
       // compute Probability for all reserved planning units
       int i, iHeavisideStepFunction;
       double rProbability, rSumProbability = 0, rShortfallPenalty;

       appendTraceFile("ComputeProbability2D start\n");

       for (i=0;i<spno;i++)
       {
           if (VarianceInExpectedAmount2D[i] > 0)
              spec[i].Zscore2D = (spec[i].target - ExpectedAmount2D[i]) / sqrt(VarianceInExpectedAmount2D[i]);
           else
               spec[i].Zscore2D = 4;

           if (spec[i].Zscore2D >= 0)
              rProbability = probZUT(spec[i].Zscore2D);
           else
               rProbability = 1 - probZUT(-1 * spec[i].Zscore2D);

           if (spec[i].ptarget2d > rProbability)
              iHeavisideStepFunction = 1;
           else
               iHeavisideStepFunction = 0;

           if (spec[i].ptarget2d > 0)
              rShortfallPenalty = (spec[i].ptarget2d - rProbability) / spec[i].ptarget2d;
           else
               rShortfallPenalty = 0;

           spec[i].probability2D = iHeavisideStepFunction * rShortfallPenalty;

           rSumProbability += spec[i].probability2D;

           #ifdef DEBUG_PROB2D
           appendTraceFile("ComputeProbability2D i %i p %lf one minus p %lf psf %lf psf2 %lf\n",
                                i,rProbability,(1-rProbability),spec[i].probability2D,(1-spec[i].probability2D));
           #endif
       }

       appendTraceFile("ComputeProbability2D sump %lf\n",rSumProbability);

       return rSumProbability * rProbabilityWeighting;
}

double probZUT(double z)
/*
Probability that a standard normal random variable has value >= z
(i.e. the area under the standard normal curve for Z in [z,+inf]

Originally adapted by Gary Perlman from a polynomial approximation in:
Ibbetson D, Algorithm 209
Collected Algorithms of the CACM 1963 p. 616
Adapted (returns upper tail instead of lower tail)

This function is not copyrighted
*/
{
       double Z_MAX = 5;
       double y, x, w;
       
       if (z == 0.0)
          x = 0.0;
       else
       {
           y = 0.5 * fabs (z);
           if (y >= (Z_MAX * 0.5))
              x = 1.0;
           else
               if (y < 1.0)
               {
                  w = y*y;
                  x = ((((((((0.000124818987 * w
                      -0.001075204047) * w +0.005198775019) * w
                      -0.019198292004) * w +0.059054035642) * w
                      -0.151968751364) * w +0.319152932694) * w
                      -0.531923007300) * w +0.797884560593) * y * 2.0;
               }
               else
               {
                   y -= 2.0;
                   x = (((((((((((((-0.000045255659 * y
                       +0.000152529290) * y -0.000019538132) * y
                       -0.000676904986) * y +0.001390604284) * y
                       -0.000794620820) * y -0.002034254874) * y
                       +0.006549791214) * y -0.010557625006) * y
                       +0.011630447319) * y -0.009279453341) * y
                       +0.005353579108) * y -0.002141268741) * y
                       +0.000535310849) * y +0.999936657524;
               }
       }
       
       return (z < 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}

double probZLT(double z)
{
       return 1.0 - probZUT(z);
}
