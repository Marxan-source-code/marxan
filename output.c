// functions that write to output files or display console output

// debug output for probability 1D
void writeProb1DDebugTable(int spno,char savename[],double ExpectedAmount1D[],double VarianceInExpectedAmount1D[], struct sspecies spec[])
{
     FILE *fp;
     int i, iHeavisideStepFunction;
     double rZ, rRawP, rP, rShortfallPenalty;

     fp = fopen(savename,"w");

     fprintf(fp,"SPID,EA,VIEA,Z,RawP,ptarget1d,HeavisideSF,ShortfallP,P\n");

     for (i=spno-1;i>=0;i--)
     {
         if (VarianceInExpectedAmount1D[i] > 0)
            rZ = (spec[i].target - ExpectedAmount1D[i]) / sqrt(VarianceInExpectedAmount1D[i]);
         else
             rZ = 4;

         if (rZ >= 0)
            rRawP = probZUT(rZ);
         else
             rRawP = 1 - probZUT(-1 * rZ);

         if (spec[i].ptarget1d > rRawP)
            iHeavisideStepFunction = 1;
         else
             iHeavisideStepFunction = 0;

         if (spec[i].ptarget1d > 0)
            rShortfallPenalty = (spec[i].ptarget1d - rRawP) / spec[i].ptarget1d;
         else
             rShortfallPenalty = 0;

         rP = iHeavisideStepFunction * rShortfallPenalty;

         fprintf(fp,"%i,%f,%f,%f,%f,%f,%i,%f,%f\n",spec[i].name,ExpectedAmount1D[i],VarianceInExpectedAmount1D[i],rZ,rRawP,spec[i].ptarget1d,iHeavisideStepFunction,rShortfallPenalty,rP);
     }

     fclose(fp);
}

// debug output for probability 2D
void writeProb2DDebugTable(int spno,char savename[],double ExpectedAmount2D[],double VarianceInExpectedAmount2D[], struct sspecies spec[])
{
     FILE *fp;
     int i, iHeavisideStepFunction;
     double rZ, rRawP, rP, rShortfallPenalty;

     fp = fopen(savename,"w");

     fprintf(fp,"SPID,CA,EA,VIEA,Z,RawP,ptarget2d,HeavisideSF,ShortfallP,P\n");

     for (i=spno-1;i>=0;i--)
     {
         if (VarianceInExpectedAmount2D[i] > 0)
            rZ = (spec[i].target - ExpectedAmount2D[i]) / sqrt(VarianceInExpectedAmount2D[i]);
         else
             rZ = 4;

         if (rZ >= 0)
            rRawP = probZUT(rZ);
         else
             rRawP = 1 - probZUT(-1 * rZ);

         if (spec[i].ptarget2d > rRawP)
            iHeavisideStepFunction = 1;
         else
             iHeavisideStepFunction = 0;

         if (spec[i].ptarget2d > 0)
            rShortfallPenalty = (spec[i].ptarget2d - rRawP) / spec[i].ptarget2d;
         else
             rShortfallPenalty = 0;

         rP = iHeavisideStepFunction * rShortfallPenalty;

         fprintf(fp,"%i,%f,%f,%f,%f,%f,%f,%i,%f,%f\n"
                   ,spec[i].name,spec[i].amount,ExpectedAmount2D[i],VarianceInExpectedAmount2D[i],rZ,rRawP,spec[i].ptarget2d,iHeavisideStepFunction,rShortfallPenalty,rP);
     }

     fclose(fp);
}

// debug output for probability 1D
void writeProb1DDetailDebugTable(char savename[],int puno,int spno,struct spustuff pu[],struct spu SM[],int R[])
{
     FILE *fp;
     int i,ipu,ism,isp;
     double *AMOUNT, *EA, *VIEA, rAmount;

     AMOUNT = (double *) calloc(spno,sizeof(double));
     EA = (double *) calloc(spno,sizeof(double));
     VIEA = (double *) calloc(spno,sizeof(double));

     fp = fopen(savename,"w");

     fprintf(fp,"PUID,R,richness,PROB");
     for (i=1;i<=spno;i++)
         fprintf(fp,",A%i",i);
     for (i=1;i<=spno;i++)
         fprintf(fp,",EA%i",i);
     for (i=1;i<=spno;i++)
         fprintf(fp,",VIEA%i",i);
     fprintf(fp,"\n");

     for (ipu=puno-1;ipu>=0;ipu--)
     {

         for (i=0;i<spno;i++)
         {
             AMOUNT[i] = 0;
             EA[i] = 0;
             VIEA[i] = 0;
         }

         if (pu[ipu].richness)
            for (i=0;i<pu[ipu].richness;i++)
            {
                ism = pu[ipu].offset + i;
                isp = SM[ism].spindex;
                rAmount = SM[ism].amount;
                AMOUNT[isp] = rAmount;
                if (rAmount)
                   if (R[ipu]==1 || R[ipu] == 2)
                   {
                      EA[isp] = rAmount * (1 - pu[ipu].prob);
                      VIEA[isp] = rAmount * rAmount * pu[ipu].prob * (1 - pu[ipu].prob);
                   }
            }

         fprintf(fp,"%i,%i,%i,%f",100-ipu,R[ipu],pu[ipu].richness,pu[ipu].prob);
         for (i=spno-1;i>=0;i--)
             fprintf(fp,",%f",AMOUNT[i]);
         for (i=spno-1;i>=0;i--)
             fprintf(fp,",%f",EA[i]);
         for (i=spno-1;i>=0;i--)
             fprintf(fp,",%f",VIEA[i]);
         fprintf(fp,"\n");
     }

     fclose(fp);

     free(AMOUNT);
     free(EA);
     free(VIEA);
}

// debug output for probability 2D
void writeProb2DDetailDebugTable(char savename[],int puno,struct spustuff pu[],struct spu SM[],int R[])
{
     FILE *fp;
     int i,ipu,ism,isp;
     double AMOUNT[9],EA[9], VIEA[9], PROB[9], rAmount;

     fp = fopen(savename,"w");

     fprintf(fp,"PUID,R,richness,P1,P2,P3,P4,P5,P6,P7,P8,P9,A1,A2,A3,A4,A5,A6,A7,A8,A9,EA1,EA2,EA3,EA4,EA5,EA6,EA7,EA8,EA9,VIEA1,VIEA2,VIEA3,VIEA4,VIEA5,VIEA6,VIEA7,VIEA8,VIEA9\n");

     for (ipu=puno-1;ipu>=0;ipu--)
     //for (ipu=0;ipu<puno;ipu++)
     {

         for (i=0;i<9;i++)
         {
             AMOUNT[i] = 0;
             EA[i] = 0;
             VIEA[i] = 0;
             PROB[i] = 0;
         }

         if (pu[ipu].richness)
            for (i=0;i<pu[ipu].richness;i++)
            {
                ism = pu[ipu].offset + i;
                isp = SM[ism].spindex;
                rAmount = SM[ism].amount;
                AMOUNT[isp] = rAmount;
                PROB[isp] = SM[ism].prob;
                if (rAmount)
                   if (R[ipu]==1 || R[ipu] == 2)
                   {
                      EA[isp] = rAmount * PROB[isp];
                      VIEA[isp] = rAmount * rAmount * PROB[isp] * (1 - PROB[isp]);
                   }
            }

         fprintf(fp,"%i,%i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"
                   ,107-ipu,R[ipu],pu[ipu].richness
                   ,PROB[8],PROB[7],PROB[6],PROB[5],PROB[4],PROB[3],PROB[2],PROB[1],PROB[0]
                   ,AMOUNT[8],AMOUNT[7],AMOUNT[6],AMOUNT[5],AMOUNT[4],AMOUNT[3],AMOUNT[2],AMOUNT[1],AMOUNT[0]
                   ,EA[8],EA[7],EA[6],EA[5],EA[4],EA[3],EA[2],EA[1],EA[0]
                   ,VIEA[8],VIEA[7],VIEA[6],VIEA[5],VIEA[4],VIEA[3],VIEA[2],VIEA[1],VIEA[0]);
     }

     fclose(fp);
}

// debug output for binary search arrays
void writeBinarySearchArrays(char *sName,struct sfname fnames, int puno, int spno, struct binsearch PULookup[],
                             struct binsearch SPLookup[])
{
    FILE *pufp, *specfp;
    int i;
    char *writename;

    writename = (char *) calloc(strlen(fnames.inputdir) + strlen(sName) + strlen("pu.csv") + 2, sizeof(char));
    strcpy(writename,fnames.inputdir);
    strcat(writename,sName);
    strcat(writename,"pu.csv");
    if ((pufp = fopen(writename,"w"))==NULL)
         displayErrorMessage("cannot create BinarySearchArrays pu file %s\n",writename);
    free(writename);
    fputs("name,index\n",pufp);
    for (i=0;i<puno;i++){
        fprintf(pufp,"%d,%d\n",PULookup[i].name,PULookup[i].index);
    }
    fclose(pufp);

    writename = (char *) calloc(strlen(fnames.inputdir) + strlen(sName) + strlen("spec.csv") + 2, sizeof(char));
    strcpy(writename,fnames.inputdir);
    strcat(writename,sName);
    strcat(writename,"spec.csv");
    if ((specfp = fopen(writename,"w"))==NULL)
         displayErrorMessage("cannot create BinarySearchArrays spec file %s\n",writename);
    free(writename);
    fputs("name,index\n",specfp);
    for (i=0;i<spno;i++){
        fprintf(specfp,"%d,%d\n",SPLookup[i].name,SPLookup[i].index);
    }
    fclose(specfp);
}

// display startup message: the program title and authors
void displayStartupMessage(void)
{
     printf("        %s \n\n   Marine Reserve Design via Annealing\n\n",sVersionString);
     printf("   Coded by Ian Ball, modified by Matthew Watts\n");
     printf("   Written by Ian Ball and Hugh Possingham\n\n");
     printf("%s\n%s\n%s\n\n",sIanBallEmail,sHughPossinghamEmail,sMattWattsEmail);
     printf("   Marxan website\n\n");
     printf("%s\n\n",sMarxanWebSite);

}

// display shutdown message when verbosity > 0
void displayShutdownMessage(void)
{
     if (verbosity > 0)
     {
        printf("\n");
        displayTimePassed();
        printf("\n              The End \n");
        if (savelog)
        {
           fprintf(fsavelog,"\n              The End \n");
        }
     }
}

// write a sync file: tell calling app we've completed a run
void writeSlaveSyncFileRun(int iSyncRun)
{
     FILE* fsync;
     char sSyncFileName[80];

     sprintf(sSyncFileName,"sync%i",iSyncRun);

     fsync = fopen(sSyncFileName,"w");
     fprintf(fsync,sSyncFileName,"%s");
     fclose(fsync);
}

// write a sync file: tell calling app we've completed all runs
void writeSlaveSyncFile(void)
{
     FILE* fsync;

     fsync = fopen("sync","w");
     fprintf(fsync,"sync");
     fclose(fsync);
}

// displays an error message for any verbosity
// program is then terminated
void displayErrorMessage(char sMess[],...)
{
     extern jmp_buf jmpbuf;
     va_list args;

     va_start(args,sMess);
     vprintf(sMess,args);
     if (savelog) vfprintf(fsavelog,sMess,args);
        va_end(args);
     longjmp(jmpbuf,1);
}

// displays a warning message when verbosity > 0
void displayWarningMessage(char sMess[],...)
{
    va_list args;

    if (verbosity > 0)
    {
        va_start(args,sMess);
        vprintf(sMess,args);
        if (savelog) vfprintf(fsavelog,sMess,args);
        va_end(args);
    }

}

// display a progress message when verbosity > 0
void displayProgress(char sMess[],...)
{
    va_list args;

    if (verbosity > 0)
    {
        va_start(args,sMess);
        vprintf(sMess,args);
        if (savelog)
        {
                   vfprintf(fsavelog,sMess,args);
                   fflush(fsavelog);
        }
        va_end(args);
    }
}

// display a progress message when verbosity > 1
void displayProgress1(char sMess[],...)
{
    va_list args;

    if (verbosity > 1)
    {
        va_start(args,sMess);
        vprintf(sMess,args);
        if (savelog)
        {
            vfprintf(fsavelog,sMess,args);
            fflush(fsavelog);
        }
        va_end(args);
    }
}

// display a progress message when verbosity > 2 :(was 2, now 5)
void displayProgress2(char sMess[],...)
{
    va_list args;

    if (verbosity > 5)
    {
        va_start(args,sMess);
        vprintf(sMess,args);
        if (savelog)
        {
            vfprintf(fsavelog,sMess,args);
            fflush(fsavelog);
        }
        va_end(args);
    }
}

// display a progress message when verbosity > 3 :(was 3, now 5)
void displayProgress3(char sMess[],...)
{
    va_list args;

    if (verbosity > 5)
    {
        va_start(args,sMess);
        vprintf(sMess,args);
        if (savelog)
        {
            vfprintf(fsavelog,sMess,args);
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
           fdebugtrace = fopen(sTraceFileName,"w");
           fflush(fdebugtrace);
           fclose(fdebugtrace);
        }
     }
}

// append message to a trace file when verbosity > 2
void appendTraceFile(char sMess[],...)
{
     FILE* fdebugtrace;
     va_list args;

     if (verbosity > 2)
     {
        if (fnames.savedebugtracefile)
        {
           va_start(args,sMess);
        
           {
                   fdebugtrace = fopen(sTraceFileName,"a");
                   vfprintf(fdebugtrace,sMess,args);
                   fclose(fdebugtrace);
           }
           va_end(args);
        }
     }
}

// create a debug file
void createDebugFile(char sFileName[],char sHeader[],struct sfname fnames)
{
     FILE* fdebugtrace;
     char *writename;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen(sFileName) + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,sFileName);
     fdebugtrace = fopen(writename,"w");
     free(writename);

     fprintf(fdebugtrace,sHeader,"%s");
     fflush(fdebugtrace);
     fclose(fdebugtrace);
}

// append message to a debug file
void appendDebugFile(char sFileName[],char sLine[],struct sfname fnames)
{
     FILE* fdebugtrace;
     char *writename;

     writename = (char *) calloc(strlen(fnames.outputdir) + strlen(sFileName) + 2, sizeof(char));
     strcpy(writename,fnames.outputdir);
     strcat(writename,sFileName);
     fdebugtrace = fopen(writename,"a");
     free(writename);

     fprintf(fdebugtrace,sLine,"%s");
     fclose(fdebugtrace);
}

// display how many seconds since program started
void displayTimePassed(void)
{
     int itemp;
     itemp = (int) clock()/CLOCKS_PER_SEC;
     printf("Time passed so far is ");
     if (itemp >= 60*60)
        printf(" %i hour%c,%i min%c and %i secs \n",
               itemp/3600,((itemp/3600==1)?' ':'s'),
               (itemp/60)%60,((itemp/60==1)?' ':'s'),itemp%60);
     else
     {
         if (itemp >=60 )
            printf(" %i min%c and %i secs \n",itemp/60,((itemp/60==1)?' ':'s'),itemp%60);
         else
             printf("%i secs \n",itemp);
     }

     if (savelog)
     {
        {
                fprintf(fsavelog,"Time passed so far is ");
                if (itemp >= 60*60)
                   fprintf(fsavelog," %i hour%c,%i min%c and %i secs \n",
                           itemp/3600,((itemp/3600==1)?' ':'s'),
                           (itemp/60)%60,((itemp/60==1)?' ':'s'),itemp%60);
                else
                {
                    if (itemp >=60 )
                       fprintf(fsavelog," %i min%c and %i secs \n",itemp/60,((itemp/60==1)?' ':'s'),itemp%60);
                    else
                        fprintf(fsavelog,"%i secs \n",itemp);
                }
        }
     }
}

// create a log file, or reset a log file
void createLogFile(int my_savelog, char* my_savelogname)
{
     if (savelog)
     {
        fclose(fsavelog);
        free(savelogname);
     } // close and delete old savelog info

     savelog = my_savelog;

     if (savelog)
     {
        savelogname = calloc(strlen(my_savelogname)+1,sizeof(char));
        strcpy(savelogname,my_savelogname);
        // Try to open file and complain if it don't work
        fsavelog = fopen(savelogname,"w");
        if (fsavelog==NULL)
        {
           free(savelogname);
           savelog = 0;
           displayErrorMessage("Error: Cannot save to log file %s \n",savelogname);
        }  // open failed

        // Header printing
        fprintf(fsavelog,"        %s \n\n   Marine Reserve Design via Annealing\n\n",sVersionString);
        fprintf(fsavelog,"   Coded by Ian Ball, modified by Matthew Watts\n");
        fprintf(fsavelog,"   Written by Ian Ball and Hugh Possingham\n\n");
        fprintf(fsavelog,"%s\n%s\n%s\n\n",sIanBallEmail,sHughPossinghamEmail,sMattWattsEmail);
        fprintf(fsavelog,"   Marxan website\n\n");
        fprintf(fsavelog,"%s\n\n",sMarxanWebSite);
     } // save log has just been turned on
}

// write an asymmetric connection file
void writeAsymmetricConnectionFile(int puno,struct sconnections connections[],struct spustuff pu[],struct sfname fnames)
{
    int i;
    FILE *fp;
    char *writename;
    struct sneighbour *p;

    writename = (char *) calloc(22 + strlen(fnames.outputdir)+2, sizeof(char));
    strcpy(writename,fnames.outputdir);
    strcat(writename,"debug_asymmetric_connectivity.csv");
    if ((fp = fopen(writename,"w"))==NULL)
    {
       displayProgress1("Warning: Cannot create file %s",writename);
       free(writename);
    }
    free(writename);

    fprintf(fp,"idA,idB,connectionorigon\n");
    for (i=0;i<puno;i++)
    {
        for (p = connections[i].first;p;p=p->next)
            fprintf(fp,"%i,%i,%i,%lf\n",pu[i].id,pu[p->nbr].id,p->connectionorigon,p->cost);
    }

    fclose(fp);
}

// write an output file from the loaded sparse matrix
void writeSparseMatrix(int iSMno,int puno, struct spustuff PU[], struct sspecies spec[], struct spu SM[],struct sfname fnames)
{
     FILE *fp;
     int i,j;
     char *writename;

     writename = (char *) calloc(strlen(fnames.inputdir) + strlen("sm.csv") + 2, sizeof(char));
     strcpy(writename,fnames.inputdir);
     strcat(writename,"sm.csv");
     if ((fp = fopen(writename,"w"))==NULL)
        displayErrorMessage("cannot create PUvSpecies file %s\n",writename);
     free(writename);

     fputs("species,pu,amount,prob\n",fp);
     for (i=puno-1;i>=0;i--)
     {
         if (PU[i].richness > 0)
         {
            for (j=0;j<PU[i].richness;j++)
                fprintf(fp,"%i,%i,%g,%g\n"
                          ,spec[SM[PU[i].offset + j].spindex].name
                          ,PU[i].id
                          ,SM[PU[i].offset + j].amount
                          ,SM[PU[i].offset + j].prob);
         }
     }

     fclose(fp);
}

// write a summary file
//  imode = 1   Output Summary Stats only
//  imode = 2   Output Everything
void writeSummary(int puno,int spno,int R[],struct sspecies spec[],struct scost reserve,
                   int itn,char savename[],double misslevel,int imode)
{
     FILE *fp;  // Imode = 1, REST output, Imode = 2, Arcview output
     int i,ino=0,isp=0;
     double shortfall,connectiontemp,rMPM, rConnectivityFraction,
            rConnectivityTotal = 0,rConnectivityIn = 0,rConnectivityEdge = 0,rConnectivityOut = 0;
     char sDelimiter[20];

     {
             if (itn==1)
                fp = fopen(savename,"w");
             else
                 fp = fopen(savename,"a");

             if (imode > 1)
                strcpy(sDelimiter,",");
             else
                 strcpy(sDelimiter,"\t");

             if (!fp)
                displayErrorMessage("Cannot save output to %s \n",savename);

             // Ouput the Summary Statistics
             for (i=0;i<=puno-1;i++)
                if (R[i]==1 || R[i]==2)
                   ino ++;
             isp = CountMissing(spno,spec,misslevel,&shortfall,&rMPM);

             #ifdef DEBUG_COUNTMISSING
             appendTraceFile("OutputSummary shortfall %g\n",shortfall);
             #endif

             ComputeConnectivityIndices(&rConnectivityTotal,&rConnectivityIn,&rConnectivityEdge,&rConnectivityOut,
                                        puno,R,connections);

             if (rConnectivityTotal > 0)
                rConnectivityFraction = rConnectivityIn / rConnectivityTotal;
             else
                 rConnectivityFraction = 0;

             for (i=0,connectiontemp = 0;i<puno;i++)
                 if (R[i]==1 || R[i] == 2)
                 {
                    connectiontemp += ConnectionCost2(i,connections,R,1,0,1);
                 } // Find True (non modified) connection
             if (itn==1)
             {
                fprintf(fp,"\"Run_Number\"%s\"Score\"%s\"Cost\"%s\"Planning_Units\"%s\"Connectivity\"%s\"Connectivity_Total\"%s\"Connectivity_In\"%s\"Connectivity_Edge\"%s\"Connectivity_Out\"%s\"Connectivity_In_Fraction\"%s\"Penalty\"",
                           sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter);
                if (fProb1D == 1)
                   fprintf(fp,"%s\"Probability1D\"",sDelimiter);
                if (fProb2D == 1)
                   fprintf(fp,"%s\"Probability2D\"",sDelimiter);
                fprintf(fp,"%s\"Shortfall\"%s\"Missing_Values\"%s\"MPM\"\n",sDelimiter,sDelimiter,sDelimiter);
             }
             fprintf(fp,"%i%s%f%s%f%s%i%s%f%s%f%s%f%s%f%s%f%s%f%s%f"
                     ,itn,sDelimiter,reserve.total,sDelimiter,reserve.cost,sDelimiter,ino,sDelimiter
                     ,connectiontemp,sDelimiter,rConnectivityTotal,sDelimiter,rConnectivityIn,sDelimiter,rConnectivityEdge,sDelimiter
                     ,rConnectivityOut,sDelimiter,rConnectivityFraction,sDelimiter,reserve.penalty);
             if (fProb1D == 1)
                fprintf(fp,"%s%f",sDelimiter,reserve.probability1D);
             if (fProb2D == 1)
                fprintf(fp,"%s%f",sDelimiter,reserve.probability2D);
             fprintf(fp,"%s%f%s%i%s%f\n",sDelimiter,shortfall,sDelimiter,isp,sDelimiter,rMPM);
             fclose(fp);
     }
     return;

} /** Output Summary ***/

// write the contents of the spec data structure.
// used to validate if the input spec file has been read as intended.
void writeSpec(int spno,struct sspecies spec[],char savename[])
{
     FILE *fp;
     int i;

     fp = fopen(savename,"w");
     if (!fp)
        displayErrorMessage("Cannot save output to %s \n",savename);

     // char *sname

     fprintf(fp,"id,name,target,prop,type,spf,target2,targetocc,sepdistance,sepnum,ptarget1d,ptarget2d\n");

     for (i=0;i<spno;i++)
         fprintf(fp,"%i,%s,%f,%f,%i,%f,%f,%i,%f,%i,%f,%f\n",
                    spec[i].name,
                    spec[i].sname,spec[i].target,
                    spec[i].prop,spec[i].type,
                    spec[i].spf,spec[i].target2,
                    spec[i].targetocc,spec[i].sepdistance,
                    spec[i].sepnum,spec[i].ptarget1d,
                    spec[i].ptarget2d);

     fclose(fp);
}

// write the contents of the pu data structure. used to validate if the input pu file has been read as intended
void writePu(int puno,struct spustuff pu[],char savename[])
{
     FILE *fp;
     int i;

     fp = fopen(savename,"w");
     if (!fp)
        displayErrorMessage("Cannot save output to %s \n",savename);

     fprintf(fp,"id,status,cost,prob,xloc,yloc\n");

     for (i=0;i<puno;i++)
         fprintf(fp,"%i,%i,%f,%f,%f,%f\n",
                    pu[i].id,pu[i].status,pu[i].cost,pu[i].prob,pu[i].xloc,pu[i].yloc);

     fclose(fp);
}

// write the penalty calculated for each species
void writePenalty(int spno,struct sspecies spec[],char savename[],int iOutputType)
{
     FILE *fp;  // Imode = 1, REST output, Imode = 2, Arcview output
     int i;
     char sDelimiter[20];

     fp = fopen(savename,"w");
     if (!fp)
        displayErrorMessage("Cannot save output to %s \n",savename);

     if (iOutputType > 1)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"\t");

     fprintf(fp,"spid%spenalty\n",sDelimiter);

     // Ouput the Summary Statistics
     for (i=0;i<spno;i++)
         fprintf(fp,"%i%s%g\n",spec[i].name,sDelimiter,spec[i].penalty);

     fclose(fp);
} // Output Penalty

// write the set of planning units used to calculate penalty
void writePenaltyPlanningUnits(int puno,struct spustuff pu[],int Rtemp[],char savename[],int iOutputType)
{
     FILE *fp;  // Imode = 1, REST output, Imode = 2, Arcview output
     int i;
     char sDelimiter[20];

     fp = fopen(savename,"w");
     if (!fp)
        displayErrorMessage("Cannot save output to %s \n",savename);

     if (iOutputType > 1)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"\t");

     fprintf(fp,"puid%sR\n",sDelimiter);

     // Ouput the Summary Statistics
     for (i=0;i<puno;i++)
         fprintf(fp,"%i%s%i\n",pu[i].id,sDelimiter,Rtemp[i]);

     fclose(fp);
}

// create a solutions matrix file
void createSolutionsMatrix(int puno,struct spustuff pu[],char savename_ism[],int iOutputType,int iIncludeHeaders)
{
     FILE *fp;
     int i;
     char sDelimiter[20];

     fp = fopen(savename_ism,"w");
     if (!fp)
        displayErrorMessage("Cannot save output to %s \n",savename_ism);

     if (iIncludeHeaders == 1)
     {
        if (iOutputType > 1)
           strcpy(sDelimiter,",");
        else
            strcpy(sDelimiter,"    ");

        fprintf(fp,"SolutionsMatrix");

        for (i=(puno-1);i>(-1);i--)
            fprintf(fp,"%sP%i",sDelimiter,pu[i].id);

        fprintf(fp,"\n");
     }

     fclose(fp);
}

// append an entry to a solutions matrix file
void appendSolutionsMatrix(int iRun,int puno,int R[],char savename[],int iOutputType,int iIncludeHeaders)
{
     FILE *fp;
     int i, iStatus;
     char sDelimiter[20];

     if (iOutputType > 1)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"\t");

     {
             fp = fopen(savename,"a");
             if (!fp)
                displayErrorMessage("Cannot save output to %s \n",savename);

             if (iIncludeHeaders == 1)
             {
                fprintf(fp,"S%i%s",iRun,sDelimiter);
             }

             for (i=(puno-1);i>(-1);i--)
             {
                 if (i < (puno-1))
                    fprintf(fp,"%s",sDelimiter);

                 iStatus = R[i];
                 if (R[i] == 3)
                    iStatus = 0;
                 if (R[i] == 2)
                    iStatus = 1;

                 fprintf(fp,"%i",iStatus);
             }

             fprintf(fp,"\n");
             fclose(fp);
     }
}

// create a solution file: output_r0001.csv, output_best.csv
void writeSolution(int puno,int R[],struct spustuff pu[],char savename[],int imode,struct sfname fnames)
{
     FILE *fp;  /* Imode = 1, REST output, Imode = 2, Arcview output */
     int i;

     fp = fopen(savename,"w");
     if (!fp)
        displayErrorMessage("Cannot save output to %s \n",savename);

     if (imode == 3)
        fprintf(fp,"PUID,%s\n",fnames.bestfieldname);
     else
         if (imode == 2)
            fprintf(fp,"\"planning_unit\",\"solution\"\n");

     for (i=puno-1;i>-1;i--)
         if (R[i]==1 || R[i] == 2)
         {
            fprintf(fp,"%i",pu[i].id);
            if (imode > 1)
               fprintf(fp,",1");
            fprintf(fp,"\n");
         }
         else
             fprintf(fp,"%i,0\n",pu[i].id);
     fclose(fp);
}

// write scenario file: a text file with input parameters
void writeScenario(int puno,int spno,double prop,double cm,
                    struct sanneal anneal,int seedinit,long int repeats,int clumptype,
                    int runopts,int heurotype,double costthresh, double tpf1, double tpf2,
                    char savename[])
{
     FILE *fp;
     char temp[40];
     fp = fopen(savename,"w");
     if (!fp)
        displayErrorMessage("Cannot save output to %s \n",savename);

     fprintf(fp,"Number of Planning Units %i\n",puno);
     fprintf(fp,"Number of Conservation Values %i\n",spno);
     fprintf(fp,"Starting proportion %.2f\n",prop);
     fprintf(fp,"Connection modifier %.2f\n\n",cm);
     switch (clumptype)
     {
            case 0:strcpy(temp,"Clumping - default step function");break;
            case 1: strcpy(temp,"Clumping - two level step function.");break;
            case 2: strcpy(temp,"Clumping - rising benefit function");break;
     }
     fprintf(fp,"%s\n",temp);

     /* Use character array here and set up the name of the algorithm used */
     switch (runopts)
     {
            case 0: strcpy(temp,"Annealing and Heuristic");break;
            case 1: strcpy(temp,"Annealing and Iterative Improvement");break;
            case 2: strcpy(temp,"Annealing and Both");break;
            case 3: strcpy(temp,"Heuristic only");break;
            case 4: strcpy(temp,"Iterative Improvement only");break;
            case 5: strcpy(temp,"Heuristic and Iterative Improvement");
     }
     fprintf(fp,"Algorithm Used :%s\n",temp);
     if (runopts == 0 || runopts == 3 || runopts == 5)
     {
        switch (heurotype)
        {
               case 0: strcpy(temp,"Richness");break;
               case 1: strcpy(temp,"Greedy");break;
               case 2: strcpy(temp,"Maximum Rarity");break;
               case 3: strcpy(temp,"Best Rarity");break;
               case 4: strcpy(temp,"Average Rarity");break;
               case 5: strcpy(temp,"Summation Rarity");break;
               case 6: strcpy(temp,"Product Irreplaceability");break;
               case 7: strcpy(temp,"Summation Irreplaceability");break;
               default: strcpy(temp,"Unkown Heuristic Type");
        }
        fprintf(fp,"Heuristic type : %s\n",temp);
     }
     else
         fprintf(fp,"No Heuristic used \n");

     if (runopts <=2)
     {
        fprintf(fp,"Number of iterations %ld\n",anneal.iterations);
        if (anneal.Tinit >= 0)
        {
           fprintf(fp,"Initial temperature %.2f\n",anneal.Tinit);
           fprintf(fp,"Cooling factor %.6f\n",anneal.Tcool);
        }
        else
        {
            fprintf(fp,"Initial temperature set adaptively\n");
            fprintf(fp,"Cooling factor set adaptively\n");
        }
        fprintf(fp,"Number of temperature decreases %li\n\n",anneal.Titns);
     }
     else
     {
         fprintf(fp,"Number of iterations N/A\nInitial temperature N/A\nCooling Factor N/A\n");
         fprintf(fp,"Number of temperature decreases N/A\n\n");
     }


     if (costthresh)
     {
        fprintf(fp,"Cost Threshold Enabled: %f\n",costthresh);
        fprintf(fp,"Threshold penalty factor A %.2f\n",tpf1);
        fprintf(fp,"Threshold penalty factor B %.2f\n\n",tpf2);
     }
     else
     {
         fprintf(fp,"Cost Threshold Disabled\nThreshold penalty factor A N/A\n");
         fprintf(fp,"Threshold penalty factor B N/A\n\n");
     }

     fprintf(fp,"Random Seed %i\n",seedinit);
     fprintf(fp,"Number of runs %ld\n",repeats);
     fclose(fp);
}  // Output Scenario

// write a species file - the missing values file: output_mv1.csv output_mvbest.csv
void writeSpecies(int spno,struct sspecies spec[],char savename[],int imode,double misslevel)
{
     FILE *fp; // Imode = 1, Tab Delimitted Text output, Imode = 2, Arcview output
     int isp, iHeavisideStepFunction;
     char temp[4];
     double rMPM, rTestMPM, rRawP, rShortfallPenalty;
     char sDelimiter[20];

     fp = fopen(savename,"w");
     if (!fp)
        displayErrorMessage("Cannot save output to %s \n",savename);

     if (imode > 1)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"\t");

     fprintf(fp,"\"Conservation Feature\"%s\"Feature Name\"%s\"Target\"%s",sDelimiter,sDelimiter,sDelimiter);
     fprintf(fp,"\"Amount Held\"%s\"Occurrence Target \"%s\"Occurrences Held\"%s",sDelimiter,sDelimiter,sDelimiter);
     fprintf(fp,"\"Separation Target \"%s\"Separation Achieved\"%s\"Target Met\"%s\"MPM\"",sDelimiter,sDelimiter,sDelimiter);

     if (fProb1D == 1)
        fprintf(fp,"%sptarget1d%sEA1D%sVIEA1D%sZ1D%srawP1D%sheavisideSF1D%sshortfallP1D%sP1D",sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter);
     if (fProb2D == 1)
        fprintf(fp,"%sptarget2d%sEA2D%sVIEA2D%sZ2D%srawP2D%sheavisideSF2D%sshortfallP2D%sP2D",sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter);

     fprintf(fp,"\n");

     for (isp=0;isp<spno;isp++)
     {
         rMPM = 1;

         fprintf(fp,"%i%s%s%s",spec[isp].name,sDelimiter,spec[isp].sname,sDelimiter);
         fprintf(fp,"%lf%s%lf%s%i%s%i%s",spec[isp].target,sDelimiter,spec[isp].amount,sDelimiter,
                    spec[isp].targetocc,sDelimiter,spec[isp].occurrence,sDelimiter);
         fprintf(fp,"%i%s%i",spec[isp].sepnum,sDelimiter,spec[isp].separation);
         strcpy(temp,"");  // use MISSLEVEL when computing target met
         if (spec[isp].target)
         {
            strcpy(temp,"yes");
            if (spec[isp].amount/spec[isp].target < misslevel)
               strcpy(temp,"no");

            rTestMPM = spec[isp].amount/spec[isp].target;
            if (rTestMPM < rMPM)
               rMPM = rTestMPM;
         }
         if (spec[isp].targetocc)
         {
            strcpy(temp,"yes");
            if (spec[isp].occurrence/spec[isp].targetocc < misslevel)
               strcpy(temp,"no");

            rTestMPM = spec[isp].occurrence/spec[isp].targetocc;
            if (rTestMPM < rMPM)
               rMPM = rTestMPM;
         }
         if (spec[isp].sepnum)
         {
            strcpy(temp,"yes");
            if (spec[isp].separation/spec[isp].sepnum < misslevel)
               strcpy(temp,"no");
         }
         fprintf(fp,"%s%s",sDelimiter,temp);
         fprintf(fp,"%s%lf",sDelimiter,rMPM);

         if (fProb1D == 1)
         {
           if (spec[isp].variance1D > 0)
              spec[isp].Zscore1D = (spec[isp].target - spec[isp].expected1D) / sqrt(spec[isp].variance1D);
           else
               spec[isp].Zscore1D = 4;

           if (spec[isp].Zscore1D >= 0)
              rRawP = probZUT(spec[isp].Zscore1D);
           else
               rRawP = 1 - probZUT(-1 * spec[isp].Zscore1D);

           if (spec[isp].ptarget1d > rRawP)
              iHeavisideStepFunction = 1;
           else
               iHeavisideStepFunction = 0;

           if (spec[isp].ptarget1d > 0)
              rShortfallPenalty = (spec[isp].ptarget1d - rRawP) / spec[isp].ptarget1d;
           else
               rShortfallPenalty = 0;

            // "ptarget1d EA1D VIEA1D Z1D rawP1D heavisideSF1D shortfallP1D P1D"
            fprintf(fp,"%s%lf%s%lf%s%lf%s%lf%s%lf%s%i%s%lf%s%lf"
                      ,sDelimiter,spec[isp].ptarget1d
                      ,sDelimiter,spec[isp].expected1D
                      ,sDelimiter,spec[isp].variance1D
                      ,sDelimiter,spec[isp].Zscore1D
                      ,sDelimiter,rRawP
                      ,sDelimiter,iHeavisideStepFunction
                      ,sDelimiter,rShortfallPenalty
                      ,sDelimiter,spec[isp].probability1D);
         }
         if (fProb2D == 1)
         {
           if (spec[isp].variance2D > 0)
              spec[isp].Zscore2D = (spec[isp].target - spec[isp].expected2D) / sqrt(spec[isp].variance2D);
           else
               spec[isp].Zscore2D = 4;

           if (spec[isp].Zscore2D >= 0)
              rRawP = probZUT(spec[isp].Zscore2D);
           else
               rRawP = 1 - probZUT(-1 * spec[isp].Zscore2D);

           if (spec[isp].ptarget2d > rRawP)
              iHeavisideStepFunction = 1;
           else
               iHeavisideStepFunction = 0;

           if (spec[isp].ptarget2d > 0)
              rShortfallPenalty = (spec[isp].ptarget2d - rRawP) / spec[isp].ptarget2d;
           else
               rShortfallPenalty = 0;

            // "ptarget2d EA2D VIEA2D Z2D rawP1D heavisideSF1D shortfallP1D P2D"
            fprintf(fp,"%s%lf%s%lf%s%lf%s%lf%s%lf%s%i%s%lf%s%lf"
                      ,sDelimiter,spec[isp].ptarget2d
                      ,sDelimiter,spec[isp].expected2D
                      ,sDelimiter,spec[isp].variance2D
                      ,sDelimiter,spec[isp].Zscore2D
                      ,sDelimiter,rRawP
                      ,sDelimiter,iHeavisideStepFunction
                      ,sDelimiter,rShortfallPenalty
                      ,sDelimiter,spec[isp].probability2D);
         }

         fprintf(fp,"\n");
     }

     fclose(fp);
}  // Output missing species information with new information

// write summed solution file output_ssoln.csv
void writeSumSoln(int puno,int sumsoln[],struct spustuff pu[],char savename[],int imode)
{
     FILE *fp;  // Imode = 1, REST output, Imode = 2, Arcview output
     int i;
     char sDelimiter[20];

     fp = fopen(savename,"w");
     if (!fp)
        displayErrorMessage("Cannot save output to %s \n",savename);

     if (imode > 1)
     {
        fprintf(fp,"\"planning_unit\",\"number\"\n");
        strcpy(sDelimiter,",");
     }
     else
         strcpy(sDelimiter,"\t");

     for (i=0;i<puno;i++)
         fprintf(fp,"%i%s%i\n",pu[i].id,sDelimiter,sumsoln[i]);

     fclose(fp);
}

// write planning unit richness to a file output_richness.csv
void writeRichness(int puno,struct spustuff pu[],char savename[],int iOutputType)
{
     FILE *fp;
     int i;
     char sDelimiter[20];

     fp = fopen(savename,"w");
     if (!fp)
        displayErrorMessage("Cannot save output to %s \n",savename);

     if (iOutputType > 1)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"    ");

     fprintf(fp,"puid%srichness\n",sDelimiter);


     for (i=0;i<puno;i++)
         fprintf(fp,"%i%s%i\n",pu[i].id,sDelimiter,pu[i].richness);

     fclose(fp);
}

// compute total area available, reserved, excluded. write it to a file if verbosity > 3.
void computeTotalAreas(int puno,int spno,struct spustuff pu[],struct sspecies spec[],struct spu SM[])
{
     int ipu, i, ism, isp, *TotalOccurrences, *TO_2, *TO_3;
     double *TotalAreas, *TA_2, *TA_3;
     FILE* TotalAreasFile;

     if (verbosity > 3)
     {
        TotalOccurrences = (int *) calloc(spno,sizeof(int));
        TO_2 = (int *) calloc(spno,sizeof(int));
        TO_3 = (int *) calloc(spno,sizeof(int));
        TotalAreas = (double *) calloc(spno,sizeof(double));
        TA_2 = (double *) calloc(spno,sizeof(double));
        TA_3 = (double *) calloc(spno,sizeof(double));

        for (i=0;i<spno;i++)
        {
            TotalAreas[i] = 0;
            TA_2[i] = 0;
            TA_3[i] = 0;
        }

        for (ipu=0;ipu<puno;ipu++)
            if (pu[ipu].richness)
               for (i=0;i<pu[ipu].richness;i++)
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

        TotalAreasFile = fopen("MarOptTotalAreas.csv","w");
        fprintf(TotalAreasFile,"spname,spindex,totalarea,reservedarea,excludedarea,targetarea,totalocc,reservedocc,excludedocc,targetocc\n");
        for (i=0;i<spno;i++)
            fprintf(TotalAreasFile,"%i,%i,%g,%g,%g,%g,%i,%i,%i,%i\n"
                                  ,spec[i].name,i,TotalAreas[i],TA_2[i],TA_3[i],spec[i].target
                                  ,TotalOccurrences[i],TO_2[i],TO_3[i],spec[i].targetocc);
        fclose(TotalAreasFile);

        free(TotalOccurrences);
        free(TO_2);
        free(TO_3);
        free(TotalAreas);
        free(TA_2);
        free(TA_3);
     }
}

// compute total area available, reserved, excluded. write it to a file output_totalareas.csv
void writeTotalAreas(int puno,int spno,struct spustuff pu[],struct sspecies spec[],struct spu SM[],char savename[],int iOutputType)
{
     int ipu, i, ism, isp, *TotalOccurrences, *TO_2, *TO_3;
     double *TotalAreas, *TA_2, *TA_3;
     FILE* TotalAreasFile;
     char sDelimiter[20];

     TotalOccurrences = (int *) calloc(spno,sizeof(int));
     TO_2 = (int *) calloc(spno,sizeof(int));
     TO_3 = (int *) calloc(spno,sizeof(int));
     TotalAreas = (double *) calloc(spno,sizeof(double));
     TA_2 = (double *) calloc(spno,sizeof(double));
     TA_3 = (double *) calloc(spno,sizeof(double));

     for (i=0;i<spno;i++)
     {
         TotalAreas[i] = 0;
         TA_2[i] = 0;
         TA_3[i] = 0;
     }

     for (ipu=0;ipu<puno;ipu++)
         if (pu[ipu].richness)
            for (i=0;i<pu[ipu].richness;i++)
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

     if (iOutputType > 1)
        strcpy(sDelimiter,",");
     else
         strcpy(sDelimiter,"\t");

     TotalAreasFile = fopen(savename,"w");
     fprintf(TotalAreasFile,"spname%stotalarea%sreservedarea%sexcludedarea%stargetarea%stotalocc%sreservedocc%sexcludedocc%stargetocc\n"
                           ,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter,sDelimiter);
     for (i=0;i<spno;i++)
         fprintf(TotalAreasFile,"%i%s%g%s%g%s%g%s%g%s%i%s%i%s%i%s%i\n"
                               ,spec[i].name,sDelimiter,TotalAreas[i],sDelimiter,TA_2[i],sDelimiter,TA_3[i],sDelimiter
                               ,spec[i].target,sDelimiter,TotalOccurrences[i],sDelimiter,TO_2[i],sDelimiter,TO_3[i],sDelimiter,spec[i].targetocc);
     fclose(TotalAreasFile);

     free(TotalOccurrences);
     free(TO_2);
     free(TO_3);
     free(TotalAreas);
     free(TA_2);
     free(TA_3);
}

// write vector R (status of each planning unit) to file. debug aid for annealing algorithms
void writeR(int iMessage,char sMessage[],int puno,int R[],struct spustuff pu[],struct sfname fnames)
{
     FILE *fp;
     char *writename;
     int i;
     char messagebuffer[80];

     sprintf(messagebuffer,"%s%i",sMessage,iMessage);

     appendTraceFile("writeR %i start\n",iMessage);

     writename = (char *) calloc(strlen(fnames.inputdir) + strlen("debugR_.csv") + strlen(messagebuffer) + 2, sizeof(char));
     strcpy(writename,fnames.inputdir);
     strcat(writename,"debugR_");
     strcat(writename,messagebuffer);
     strcat(writename,".csv");
     if ((fp = fopen(writename,"w"))==NULL)
        displayErrorMessage("cannot create writeR file %s\n",writename);
     free(writename);

     // write header row
     fprintf(fp,"puid,R\n");

     for (i=0;i<puno;i++)
     {
         fprintf(fp,"%i,%i\n",pu[i].id,R[i]);
     }

     fclose(fp);

     appendTraceFile("writeR %i end\n",iMessage);
}

// debug output for probability 1D
void writeChangeProbability1DDebugTable(char savename[],int iIteration,int ipu,int spno,struct sspecies spec[],struct spustuff pu[],struct spu SM[],int imode)
{
     FILE *fp;
     int i,ism,isp;
     double *AMOUNT,*DE,*DV,*NE,*NV,*NZ,*OZ;
     double *OHSF, *NHSF, *OSFP, *NSFP, *NRP, *ORP;

     AMOUNT = (double *) calloc(spno,sizeof(double));
     DE = (double *) calloc(spno,sizeof(double));
     DV = (double *) calloc(spno,sizeof(double));
     NE = (double *) calloc(spno,sizeof(double));
     NV = (double *) calloc(spno,sizeof(double));
     NZ = (double *) calloc(spno,sizeof(double));
     OZ = (double *) calloc(spno,sizeof(double));
     OHSF = (double *) calloc(spno,sizeof(double));
     NHSF = (double *) calloc(spno,sizeof(double));
     OSFP = (double *) calloc(spno,sizeof(double));
     NSFP = (double *) calloc(spno,sizeof(double));
     NRP = (double *) calloc(spno,sizeof(double));
     ORP = (double *) calloc(spno,sizeof(double));

     for (i=0;i<spno;i++)
     {
         AMOUNT[i] = 0;
         DE[i] = 0;
         DV[i] = 0;
         NE[i] = 0;
         NV[i] = 0;
         NZ[i] = 0;
         OZ[i] = 0;
         NRP[i] = 0;
         ORP[i] = 0;
         OHSF[i] = 0;
         NHSF[i] = 0;
         OSFP[i] = 0;
         NSFP[i] = 0;
     }

     if (pu[ipu].richness)
        for (i=0;i<pu[ipu].richness;i++)
        {
            ism = pu[ipu].offset + i;
            isp = SM[ism].spindex;
            AMOUNT[isp] = SM[ism].amount;
            if (AMOUNT[isp])
            {
               DE[isp] = imode * AMOUNT[isp] * (1-pu[ipu].prob);
               NE[isp] = spec[isp].expected1D + DE[isp];

               DV[isp] = imode * AMOUNT[isp] * AMOUNT[isp] * pu[ipu].prob * (1 - pu[ipu].prob);
               NV[isp] = spec[isp].variance1D + DV[isp];

               if (NV[isp] > 0)
                  NZ[isp] = (spec[isp].target - NE[isp]) / sqrt(NV[isp]);
               else
                   NZ[isp] = 4;

               if (NZ[isp] >= 0)
                  NRP[isp] = probZUT(NZ[isp]);
               else
                   NRP[isp] = 1 - probZUT(-1 * NZ[isp]);

               if (spec[isp].variance1D > 0)
                  OZ[isp] = (spec[isp].target - spec[isp].expected1D) / sqrt(spec[isp].variance1D);
               else
                   OZ[isp] = 4;

               if (OZ[isp] >= 0)
                  ORP[isp] = probZUT(OZ[isp]);
               else
                   ORP[isp] = 1 - probZUT(-1 * OZ[isp]);

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

     fp = fopen(savename,"w");
     fprintf(fp,"PUID,PROB,SPID,AMOUNT,DELTAEXPECTED,DELTAVARIANCE,OLDEXPECTED,NEWEXPECTED,OLDVARIANCE,NEWVARIANCE,NEWZ,OLDZ,NEWRawP,OLDRawP,NEWHEAVISIDE,OLDHEAVISIDE,NEWSHORTFALL,OLDSHORTFALL\n");
     for (i=spno-1;i>=0;i--)
     {
         fprintf(fp,"%i,%f,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"
                   ,pu[ipu].id,pu[ipu].prob,spec[i].name,AMOUNT[i]
                   ,DE[i],DV[i],spec[i].expected1D,NE[i],NV[i],spec[i].variance1D,NZ[i],OZ[i],NRP[i],ORP[i]
                   ,NHSF[i],OHSF[i],NSFP[i],OSFP[i]);
     }
     fclose(fp);

     free(AMOUNT);
     free(DE);
     free(DV);
     free(NE);
     free(NV);
     free(NZ);
     free(OZ);
     free(OHSF);
     free(NHSF);
     free(OSFP);
     free(NSFP);
     free(NRP);
     free(ORP);
}

// debug output for probability 2D
void writeChangeProbability2DDebugTable(char savename[],int iIteration,int ipu,int spno,struct sspecies spec[],struct spustuff pu[],struct spu SM[],int imode)
{
     FILE *fp;
     int i,ism,isp;
     double *AMOUNT, *DE, *DV, *NE, *NV, *NZ, *OZ, *NP, *OP, *PROB;
     double *OHSF, *NHSF, *OSFP, *NSFP;

     AMOUNT = (double *) calloc(spno,sizeof(double));
     DE = (double *) calloc(spno,sizeof(double));
     DV = (double *) calloc(spno,sizeof(double));
     NE = (double *) calloc(spno,sizeof(double));
     NV = (double *) calloc(spno,sizeof(double));
     NZ = (double *) calloc(spno,sizeof(double));
     OZ = (double *) calloc(spno,sizeof(double));
     NP = (double *) calloc(spno,sizeof(double));
     OP = (double *) calloc(spno,sizeof(double));
     PROB = (double *) calloc(spno,sizeof(double));
     OHSF = (double *) calloc(spno,sizeof(double));
     NHSF = (double *) calloc(spno,sizeof(double));
     OSFP = (double *) calloc(spno,sizeof(double));
     NSFP = (double *) calloc(spno,sizeof(double));

     for (i=0;i<spno;i++)
     {
         AMOUNT[i] = 0;
         DE[i] = 0;
         DV[i] = 0;
         NE[i] = 0;
         NV[i] = 0;
         NZ[i] = 0;
         OZ[i] = 0;
         NP[i] = 0;
         OP[i] = 0;
         PROB[i] = 0;
         OHSF[i] = 0;
         NHSF[i] = 0;
         OSFP[i] = 0;
         NSFP[i] = 0;
     }

     if (pu[ipu].richness)
        for (i=0;i<pu[ipu].richness;i++)
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
                  NP[isp] = probZUT(NZ[isp]);
               else
                   NP[isp] = 1 - probZUT(-1 * NZ[isp]);

               if (spec[isp].variance1D > 0)
                  OZ[isp] = (spec[isp].target - spec[isp].expected2D) / sqrt(spec[isp].variance2D);
               else
                   OZ[isp] = 4;

               if (OZ[isp] >= 0)
                  OP[isp] = probZUT(OZ[isp]);
               else
                   OP[isp] = 1 - probZUT(-1 * OZ[isp]);

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

     fp = fopen(savename,"w");
     fprintf(fp,"IPU,PROB,ISP,AMOUNT,DELTAEXPECTED,DELTAVARIANCE,OLDEXPECTED,NEWEXPECTED,OLDVARIANCE,NEWVARIANCE,NEWZ,OLDZ,NEWP,OLDP,NEWHEAVISIDE,OLDHEAVISIDE,NEWSHORTFALL,OLDSHORTFALL\n");
     for (i=0;i<spno;i++)
     {
         fprintf(fp,"%i,%f,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"
                   ,ipu,pu[ipu].prob,i,AMOUNT[i]
                   ,DE[i],DV[i],spec[i].expected1D,NE[i],NV[i],spec[i].variance1D,NZ[i],OZ[i],NP[i],OP[i]
                   ,NHSF[i],OHSF[i],NSFP[i],OSFP[i]);
     }
     fclose(fp);

     free(AMOUNT);
     free(DE);
     free(DV);
     free(NE);
     free(NV);
     free(NZ);
     free(OZ);
     free(NP);
     free(OP);
     free(PROB);
     free(OHSF);
     free(NHSF);
     free(OSFP);
     free(NSFP);
}

// debug output for probability 1D
void writeProbData(int puno,struct spustuff pu[],struct sfname fnames)
{
    int i;
    char *writename;
    char sLine[600];
    char *sVarVal;
    FILE *fp;

    writename = (char *) calloc(strlen("writeProbData.csv") + strlen(fnames.inputdir)+2, sizeof(char));
    strcpy(writename,fnames.inputdir);
    strcat(writename,"writeProbData.csv");
    if((fp = fopen(writename,"w"))==NULL)
        displayErrorMessage("probability file %s has not been found.\nAborting Program.",writename);
    free(writename);

    fprintf(fp,"puid,prob\n");

    for (i=0;i<puno;i++)
    {
        fprintf(fp,"%d,%lf\n",pu[i].id,pu[i].prob);
    }

    fclose(fp);
}

// make a backup copy of a connection file
void copyFile(char sInputFile[],char sOutputFile[])
{
     FILE *fpInputFile, *fpOutputFile;
     char ch;

     if ((fpInputFile = fopen(sInputFile, "rb"))!=NULL)
     // if the input file does not exist, then exit gracefully
     {
        fpOutputFile = fopen(sOutputFile, "wb");

        while (!feof(fpInputFile))
        {
              ch = fgetc(fpInputFile);
              if (!feof(fpInputFile))
                 fputc(ch, fpOutputFile);
        }

        fclose(fpInputFile);
        fclose(fpOutputFile);
     }
}

// display statistics for a configuration of planning unit
void displayValueForPUs(int puno, int spno,int *R,struct scost reserve,
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
}

// display usage information for the marxan executable
// displayed when the command line options are not understood
void displayUsage(char *programName)
{
    fprintf(stderr,"%s usage: %s -[o] -[c] [input file name]\n",programName,programName);
}

