
#include <cmath> 

#include "computation.hpp"
#include "marxan.hpp"

namespace marxan {

// functions relating to probability 1D and 2D

/************** Value of a Reserve System ********/
void ComputeP_AllPUsSelected_1D(const string &savename, int puno ,int spno, const vector<spustuff> &pu, const vector<spu> &SM, const vector<sspecies> &spec) {
   FILE *fp;
   int i,j,iHeavisideStepFunction,ism,isp;
   double rProbability, rRawP, rSumProbability = 0, rShortfallPenalty, rZScore;

   // create the output file
   fp = fopen(savename.c_str(),"w");
   fprintf(fp,"SPID,amount held,ptarget1d,EA1D,VIEA1D,Z1D,rawP1D,heavisideSF1D,shortfallP1D,P1D\n");

   // init arrays
   vector<double> ExpectedAmount1D(spno, 0), VarianceInExpectedAmount1D(spno, 0), TA(spno, 0);

   // compute EA, VIEA, TA
   for (j=0;j<puno;j++)
   {
      if (pu[j].richness)
      {
         for (i=0;i<pu[j].richness;i++)
         {
               ism = pu[j].offset + i;
               isp = SM[ism].spindex;

               ExpectedAmount1D[isp] += SM[ism].amount * (1 - pu[j].prob);
               VarianceInExpectedAmount1D[isp] += SM[ism].amount * SM[ism].amount * pu[j].prob * (1 - pu[j].prob);
               TA[isp] += SM[ism].amount;
         }
      }
   }

   // compute probability for each feature
   for (i=0;i<spno;i++)
   {
      computeProbMeasures(VarianceInExpectedAmount1D[i], spec[i].target, spec[i].ptarget1d, ExpectedAmount1D[i],
         rZScore, rRawP, iHeavisideStepFunction, rShortfallPenalty, rProbability);

      rSumProbability += rProbability;

      fprintf(fp,"%i,%f,%f,%f,%f,%f,%f,%i,%f,%f\n",
                  spec[i].name,TA[i],spec[i].ptarget1d,
                  ExpectedAmount1D[i],VarianceInExpectedAmount1D[i],rZScore,
                  rRawP,iHeavisideStepFunction,rShortfallPenalty,rProbability);
   }

   fclose(fp);

   #ifdef DEBUGTRACEFILE
   sprintf(debugbuffer,"ComputeP_AllPUsSelected_1D SumP %f SumP * PW %f\n",rSumProbability,rSumProbability * rProbabilityWeighting);
   appendTraceFile(debugbuffer);
   #endif
}

void ComputeP_AllPUsSelected_2D(const string &savename, int puno, int spno, const vector<spustuff> &pu, const vector<spu> &SM, const vector<sspecies> &spec) {
   FILE *fp;
   int i,j,iHeavisideStepFunction,ism,isp;
   double rProbability, rRawP, rSumProbability = 0, rShortfallPenalty, rZScore;

   // create the output file
   fp = fopen(savename.c_str(),"w");
   fprintf(fp,"SPID,amount held,ptarget1d,EA2D,VIEA2D,Z2D,rawP2D,heavisideSF2D,shortfallP2D,P2D\n");

   // init arrays
   vector<double> ExpectedAmount2D(spno, 0), VarianceInExpectedAmount2D(spno, 0), TA(spno, 0);

   // compute EA, VIEA, TA
   for (j=0;j<puno;j++)
   {
      if (pu[j].richness)
      {
         for (i=0;i<pu[j].richness;i++)
         {
               ism = pu[j].offset + i;
               isp = SM[ism].spindex;

               ExpectedAmount2D[isp] += SM[ism].amount * SM[ism].prob;
               VarianceInExpectedAmount2D[isp] += SM[ism].amount * SM[ism].amount * SM[ism].prob * (1 - SM[ism].prob);
               TA[isp] += SM[ism].amount;
         }
      }
   }

   // compute probability for each feature
   for (i=0;i<spno;i++)
   {
      computeProbMeasures(VarianceInExpectedAmount2D[i], spec[i].target, spec[i].ptarget2d, ExpectedAmount2D[i],
         rZScore, rRawP, iHeavisideStepFunction, rShortfallPenalty, rProbability);

      rSumProbability += rProbability;

      fprintf(fp,"%i,%f,%f,%f,%f,%f,%f,%i,%f,%f\n",
                  spec[i].name,TA[i],spec[i].ptarget2d,
                  ExpectedAmount2D[i],VarianceInExpectedAmount2D[i],rZScore,
                  rRawP,iHeavisideStepFunction,rShortfallPenalty,rProbability);
   }

   fclose(fp);

   #ifdef DEBUGTRACEFILE
   char debugbuffer[200];
   sprintf(debugbuffer,"ComputeP_AllPUsSelected_2D SumP %f SumP * PW %f\n",rSumProbability,rSumProbability * rProbabilityWeighting);
   appendTraceFile(debugbuffer);
   #endif
}

// Change in probability 1D for adding or deleting one PU
double ChangeProbability1D(int iIteration, int ipu, int spno,int puno, vector<sspecies> &spec, const vector<spustuff> &pu, const vector<spu> &SM, int imode)
{
   int i, ism, isp, iNewHeavisideStepFunction, iOrigHeavisideStepFunction;
   double rNewExpected, rNewVariance, rOriginalZScore;
   double rNewZScore, rProbOriginal, rProbNew, rSumProbability = 0;
   double rDeltaExpected, rDeltaVariance, rNewShortfallPenalty, rOrigShortfallPenalty, rProbOld;

   #ifdef DEBUG_PROB1D
   string sDebugFileName;
   if (iIteration < 10)
   {
      sprintf(sIteration,"%i",iIteration);
      sDebugFileName = fnames.outputdir + "output_ChangeProbability1DDebug_" + to_string(iIteration) + ".csv";
      writeChangeProbability1DDebugTable(sDebugFileName.c_str(),iIteration,ipu,spno,spec,pu,SM,imode);
   }
   #endif

   if (pu[ipu].richness)
   {
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
            
            // Get pr measures for new state
            computeProbMeasures(rNewVariance, spec[isp].target, spec[isp].ptarget1d, rNewExpected,
               rNewZScore, rProbNew, iNewHeavisideStepFunction, rNewShortfallPenalty, rProbNew);

            spec[isp].Zscore1D = rNewZScore;
            spec[isp].probability1D = rProbNew;

            // Get pr measures for old state
            computeProbMeasures(spec[isp].variance1D, spec[isp].target, spec[isp].ptarget1d, spec[isp].expected1D,
               rOriginalZScore, rProbOriginal, iOrigHeavisideStepFunction, rOrigShortfallPenalty, rProbOld);

            // change in probability
            rSumProbability += rProbNew - rProbOld;
         }
      }
   }

   return (rSumProbability * rProbabilityWeighting);
}  // Change in probability for adding or deleting one PU

// Change in probability 2D for adding or deleting one PU
double ChangeProbability2D(int iIteration, int ipu, int spno, int puno, vector<sspecies> &spec, const vector<spustuff> &pu, const vector<spu> &SM, int imode)
{
   int i, ism, isp, iNewHeavisideStepFunction, iOrigHeavisideStepFunction;
   double rNewExpected, rNewVariance, rOriginalZScore;
   double rProb, rNewZScore, rProbOriginal, rProbNew, rSumProbability = 0;
   double rDeltaExpected, rDeltaVariance, rNewShortfallPenalty, rOrigShortfallPenalty, rProbOld;

   #ifdef DEBUG_PROB2D
   string sDebugFileName;
   if (iIteration < 10)
   {
      sDebugFileName = fnames.outputdir + "output_ChangeProbability2DDebug_" + to_string(iIteration) + ".csv";
      writeChangeProbability2DDebugTable(sDebugFileName.c_str(),iIteration,ipu,spno,spec,pu,SM,imode);
   }
   #endif

   if (pu[ipu].richness)
   {
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

            // Get pr measures for new state
            computeProbMeasures(rNewVariance, spec[isp].target, spec[isp].ptarget2d, rNewExpected,
               rNewZScore, rProbNew, iNewHeavisideStepFunction, rNewShortfallPenalty, rProbNew);

            spec[isp].Zscore2D = rNewZScore;
            spec[isp].probability2D = rProbNew;

            // Get pr measures for old state
            computeProbMeasures(spec[isp].variance2D, spec[isp].target, spec[isp].ptarget2d, spec[isp].expected2D,
               rOriginalZScore, rProbOriginal, iOrigHeavisideStepFunction, rOrigShortfallPenalty, rProbOld);

            // change in probability
            rSumProbability += rProbNew - rProbOld;
         }
      }
   }

   return (rSumProbability * rProbabilityWeighting);
}

// accumulate ExpectedAmount and VarianceInExpectedAmount for each species at this planning unit
void ReturnProbabilityAmounts1D(vector<double> &ExpectedAmount1D, vector<double> &VarianceInExpectedAmount1D, int ipu,
                                int puno, const vector<spustuff> &pu, const vector<spu> &SM)
{
   computeExpectedAndVariance(ipu, pu, SM, VarianceInExpectedAmount1D, ExpectedAmount1D);
}

void ReturnProbabilityAmounts2D(vector<double> &ExpectedAmount2D, vector<double> &VarianceInExpectedAmount2D,int ipu,
                                int puno, const vector<spustuff> &pu, const vector<spu> &SM)
{
   computeExpectedAndVariance(ipu, pu, SM, VarianceInExpectedAmount2D, ExpectedAmount2D);
}

double ComputeProbability1D(const vector<double> &ExpectedAmount1D, const vector<double> &VarianceInExpectedAmount1D,
                            int spno,vector<sspecies> &spec) 
{
   // compute Probability for all reserved planning units
   int i, iHeavisideStepFunction;
   double rProbability, rSumProbability = 0, rShortfallPenalty;

   appendTraceFile("ComputeProbability1D start\n");

   for (i=0;i<spno;i++)
   {
      computeProbMeasures(VarianceInExpectedAmount1D[i], spec[i].target, spec[i].ptarget1d, ExpectedAmount1D[i],
               spec[i].Zscore1D, rProbability, iHeavisideStepFunction, rShortfallPenalty, spec[i].probability1D);

      rSumProbability += spec[i].probability1D;

      #ifdef DEBUG_PROB1D
      appendTraceFile("ComputeProbability1D i %i p %lf one minus p %lf pst %lf pst2 %lf\n",
                     i,rProbability,(1-rProbability),spec[i].probability1D,(1-spec[i].probability1D));
      #endif
   }

   appendTraceFile("ComputeProbability1D sump %lf\n",rSumProbability);
   return rSumProbability * rProbabilityWeighting;
}

double ComputeProbability2D(vector<double> &ExpectedAmount2D, vector<double> &VarianceInExpectedAmount2D,
                            int spno, vector<sspecies> &spec) 
{
   // compute Probability for all reserved planning units
   int i, iHeavisideStepFunction;
   double rProbability, rSumProbability = 0, rShortfallPenalty;

   appendTraceFile("ComputeProbability2D start\n");

   for (i=0;i<spno;i++)
   {
      computeProbMeasures(VarianceInExpectedAmount2D[i], spec[i].target, spec[i].ptarget2d, ExpectedAmount2D[i],
               spec[i].Zscore2D, rProbability, iHeavisideStepFunction, rShortfallPenalty, spec[i].probability2D);

      rSumProbability += spec[i].probability2D;

      #ifdef DEBUG_PROB2D
      appendTraceFile("ComputeProbability2D i %i p %lf one minus p %lf psf %lf psf2 %lf\n",
                     i,rProbability,(1-rProbability),spec[i].probability2D,(1-spec[i].probability2D));
      #endif
   }

   appendTraceFile("ComputeProbability2D sump %lf\n",rSumProbability);

   return rSumProbability * rProbabilityWeighting;
}

} // marxan