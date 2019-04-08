// functions relating to Heuristics

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
       effamount = returnAmountSpecAtPu(pu,SM,ipu,isp);
       
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
                 computeChangeScore(-1,i,spno,puno,pu,connections,spec,SM,R,cm,1,&change,reserve,
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
              computeChangeScore(-1,bestpu,spno,puno,pu,connections,spec,SM,R,cm,1,&change,reserve,
                                 costthresh,tpf1,tpf2,1,clumptype);
              doChange(bestpu,puno,R,reserve,change,pu,SM,spec,connections,1,clumptype);

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

