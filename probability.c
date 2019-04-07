// functions relating to probability 1D and 2D

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
