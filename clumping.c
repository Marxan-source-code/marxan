// functions relating to species clumping

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

// *********  Clump Utilities ******************

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
                newpclump->amount = clump->amount - returnAmountSpecAtPu(pu,SM,clumppu->puid,isp);
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
        
        rAmount = returnAmountSpecAtPu(pu,SM,clumppu->puid,isp);
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
              rAmount = returnAmountSpecAtPu(pu,SM,newhead->id,isp);
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
       rAmount = returnAmountSpecAtPu(pu,SM,ipu,isp);
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
                 rAmount = returnAmountSpecAtPu(pu,SM,ipu,isp);
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
        rAmount = returnAmountSpecAtPu(pu,SM,ipu,isp);
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
                     rAmount = returnAmountSpecAtPu(pu,SM,ppu->puid,isp);
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
        if (returnAmountSpecAtPu(pu,SM,i,isp) > 0)
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
      ipu = returnRandom(iputotal);
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
                      rAmount = returnAmountSpecAtPu(pu,SM,pcpu->puid,isp);
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
        amount = returnAmountSpecAtPu(pu,SM,ipu,isp);
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

int ValidPU(int ipu,int isp,struct sclumps *newno,struct sspecies spec[],struct spustuff pu[],
            struct spu SM[],int imode)
{
    // Returns true if ipu is acceptable as a planning unit
    int i;
    
    i = returnIndexSpecAtPu(pu,SM,ipu,isp);
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
       double rAmount = returnAmountSpecAtPu(pu,SM,ipu,isp);

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

