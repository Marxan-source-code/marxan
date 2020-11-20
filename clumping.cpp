
#include <algorithm>
#include <stdexcept>

#include "computation.hpp"

// functions relating to species clumping

namespace marxan {

   // returns the clump number of a species at a planning unit, if the species doesn't occur here, returns 0
   int rtnClumpSpecAtPu(vector<spustuff> &pu, vector<spu> &SM, int iPUIndex, int iSpecIndex, int thread)
   {
      if (pu[iPUIndex].richness > 0)
         for (int i=0;i<pu[iPUIndex].richness;i++)
            if (SM[pu[iPUIndex].offset + i].spindex == iSpecIndex)
               return SM[pu[iPUIndex].offset + i].clump[thread]; //TODO verify thread change

      return 0;
   }

   // sets the clump number of a species at a planning unit
   void setClumpSpecAtPu(vector<spustuff> &pu, vector<spu> &SM, int iPUIndex, int iSpecIndex, int iSetClump, int thread)
   {
      if (pu[iPUIndex].richness > 0)
         for (int i=0;i<pu[iPUIndex].richness;i++)
            if (SM[pu[iPUIndex].offset + i].spindex == iSpecIndex)
               SM[pu[iPUIndex].offset + i].clump[thread] = iSetClump; //TODO verify thread change
   }

   void ClearClump(int isp, sclumps &target, vector<spustuff> &pu, vector<spu> &SM, int thread) {

      /* Remove all links from this clump */
      for (int i : target.head) {
         if (rtnClumpSpecAtPu(pu,SM,i,isp, thread) == target.clumpid) {/* in case pu is in new clump */
            setClumpSpecAtPu(pu,SM,i,isp,0, thread);
         }
      }

      target.head.clear();
   }

   int ClumpCut(int isp,vector<spustuff> &pu,
         vector<sspecies> &spec, sclumps &clump,
         int &clumppu, vector<sconnections> &connections, vector<spu> &SM,
         double &totalamount,int &totalocc,
         int &iseparation, int imode, int clumptype, int thread) {

      int ineighbour = 0,iclumps = 0;
      vector<int> newhead, head;
      vector<sclumps> spclump, newpclump;
      double clumpamount, rAmount;
      int iocc;

      totalamount = 0;
      totalocc = 0;

      /* Set up spclump for counting Separation */
      if (imode)
      {  
         sclumps temp;
         temp.amount = 0;
         temp.clumpid = clumppu;
         spclump.push_back(temp);
         newpclump = spclump;
      }

      /** Generate list of all neighbours and count them **/
      /*First check if there are no neighbours then exit. **/
      /* return null for no clump cut done and need to do separation count */
      if (connections[clumppu].nbrno == 0)
      {
         if (imode)
         {
            iseparation = CountSeparation(isp,spclump,pu,SM,spec,0, thread);
         }
         return 0;
      }

      for (sneighbour pnbr: connections[clumppu].first)
      {
         if (rtnClumpSpecAtPu(pu,SM,pnbr.nbr,isp,thread) == clump.clumpid)
         {
            ineighbour++;
            newhead.push_back(pnbr.nbr);
         } // If neighbour is part of the same clump
      } // For cycling through all neighbours
      head = newhead;

      if (ineighbour <= 1)
      { // One or fewer neighbours
         if (imode)
         { // separation distance called
            for (int pclumpu: clump.head) {
               if (pclumpu != clumppu) {
                  sclumps temp;
                  temp.clumpid = pclumpu;
                  temp.amount = clump.amount - returnAmountSpecAtPu(pu[clumppu],SM,isp).second;
                  spclump.push_back(temp);
               }
            }

            iseparation = CountSeparation(isp,spclump,pu,SM,spec,0, thread);
         } else {
            iseparation = spec[isp].sepnum;
         }
         
         rAmount = returnAmountSpecAtPu(pu[clumppu],SM,isp).second;
         totalamount = PartialPen4(isp,clump.amount-rAmount,spec,clumptype);
         totalocc = (clump.occs - (rAmount > 0))*(totalamount > 0); // count only if still valid size
         
         return 0;
      }

      // More than one neighbour. Can they form their own clump?
      // Put first neighbour at head of new list
      // Since we are using vectors instead of LL, this new logic might be slightly slower.
      // Reserve size for clumplist to allow iteration while modifying (without new allocation)
      // TODO - revisit. 
      while(!head.empty()) {
         vector<int> clumplist;
         clumplist.reserve(head.size());

         clumpamount = 0;
         iclumps++;
         clumplist.push_back(head.front());
         head.erase(head.begin()); //remove first element from head
         ineighbour--;

         for (int id2: clumplist) {
            for (sneighbour pnbr :connections[id2].first)
            {
               // if neighbour in clump but not cut out one  
               if (rtnClumpSpecAtPu(pu,SM,pnbr.nbr,isp, thread) == clump.clumpid && pnbr.nbr != clumppu)  
               {
                  if (find(clumplist.begin(), clumplist.end(), pnbr.nbr) == clumplist.end()) {
                     // add item to clumplist if not already in 
                     clumplist.push_back(pnbr.nbr);

                     auto nbrFound = find(head.begin(), head.end(), pnbr.nbr);
                     if (nbrFound != head.end()) {
                        ineighbour--;
                        head.erase(nbrFound);
                     }
                  }
               }

            } // cycling through every neighbour on this clump
         } 

         iocc = 0;
         for (int id2: clumplist) {
            rAmount = returnAmountSpecAtPu(pu[id2],SM,isp).second;
            clumpamount += rAmount;
            iocc += (rAmount > 0);
         }

         totalamount += PartialPen4(isp,clumpamount,spec,clumptype);

         if (PartialPen4(isp,clumpamount,spec,clumptype)) 
         {
            totalocc += iocc;
         }

         if (imode) {
            for (int id2: clumplist) {
               sclumps temp;
               temp.clumpid = id2;
               temp.amount = clumpamount;
               spclump.push_back(temp);
            }
         } 
      }  // Continue clump formation whilst there are members in the list

      if (imode) {
         iseparation = CountSeparation(isp,spclump,pu,SM,spec,0, thread);
      }
      else
      {
         iseparation = spec[isp].sepnum;
      }
      
      return iclumps;
   }
   
   /*************** Setting Clumps for Species Aggregation Rule**********/
   void SetSpeciesClumps(int puno,vector<int> &R,vector<sspecies> &spec,vector<spustuff> &pu,
                     vector<spu> &SM,vector<sconnections> &connections,int clumptype, int thread) {
      int isp, ism;
      for (int ipu=0;ipu<puno;ipu++)
      {
         if (pu[ipu].richness)
         {
            for (int i=0;i<pu[ipu].richness;i++)
            {
               ism = pu[ipu].offset + i;
               isp = SM[ism].spindex;
               if (spec[isp].target2)
               {
                  spec[isp].clumps = 0;
                  if ((R[ipu]==1 || R[ipu]==2) && SM[ism].amount > 0 && SM[ism].clump[thread] == 0)
                  {
                     AddNewPU(ipu,isp,connections,spec,pu,SM,clumptype, thread);
                  } // Add a New planning unit
               } // For each type 4 species
            }
         }
      }
   }

   /************ Species Amounts Type 4 **************/
   /** Assumes Set Species Clumps has been called **/
   void SpeciesAmounts4(int isp,vector<sspecies> &spec,int clumptype) {
      double ftemp;
      struct sclumps *pclump;

      for (sclumps pclump: spec[isp].head) {
         ftemp = PartialPen4(isp,pclump.amount,spec,clumptype);
         spec[isp].amount += ftemp;
         spec[isp].occurrence += pclump.occs*(ftemp>0);
      }
   }

   // Clear Clumps
   // This is for clean up purposes
   void ClearClumps(int spno,vector<sspecies> &spec,vector<spustuff> &pu, vector<spu> &SM, int thread) {
      for (int i=0;i<spno;i++)
      {
         for (sclumps& clump: spec[i].head) {
            ClearClump(i,clump,pu,SM, thread);
         }

         spec[i].head.clear();
         spec[i].clumps = 0;
      } // Clear clump for each species
   }

   // Add New Clump
   sclumps AddNewClump(int isp,int ipu,vector<sspecies> &spec,vector<spustuff> &pu, vector<spu> &SM, int thread) {
      int iclumpno = 0;
      double rAmount;

      if (spec[isp].head.empty()) {
         iclumpno = 1;
      }

      int ind = 0;
      bool found = false;
      sclumps prev;
      prev.clumpid = -1; // placeholder initial value
      for (sclumps& pclump: spec[isp].head) {
         if (prev.clumpid == -1) {
            prev = pclump;
            continue;
         }
         else {
            if (pclump.clumpid - prev.clumpid > 1)
            {
                  iclumpno = prev.clumpid+1;
                  found = true;
                  break;
            } // Looking for good number
         }
         ind ++;
         prev = pclump;
      }

      if (!found){
         iclumpno = spec[isp].head.back().clumpid+1;
      }

      setClumpSpecAtPu(pu,SM,ipu,isp,iclumpno, thread);

      // Set up new clump to insert
      sclumps temp;
      temp.clumpid = iclumpno;
      temp.head.push_back(ipu);
      rAmount = returnAmountSpecAtPu(pu[ipu],SM,isp).second;
      temp.amount = rAmount;
      temp.occs = (rAmount > 0);

      if (!spec[isp].head.empty())
      {
         spec[isp].head.insert(spec[isp].head.begin() + ind, temp);
      } // Stick clump into correct location
      else
      {
         spec[isp].head.push_back(temp);
      } // First clump on the block

      spec[isp].clumps++;
      return temp; // TODO ADBAI - flesh out return value
   }

   // Add New Planning Unit for a given Species
   void AddNewPU(int ipu,int isp,vector<sconnections> &connections,vector<sspecies> &spec,vector<spustuff> &pu,
               vector<spu> &SM, int clumptype, int thread) {
      int ineighbours = 0;
      int iclumpno, iClump;
      sclumps temp;
      struct sclumppu *pnewclumppu;
      double ftemp, rAmount;

      for (sneighbour pnbr: connections[ipu].first) {
         // Check all the neighbours to see if any are already in clumps */
         iClump = rtnClumpSpecAtPu(pu,SM,pnbr.nbr,isp, thread);
         if (iClump > 0) {
            // Neighbour that is part of clump
            int ind;
            ineighbours++;
            if (ineighbours == 1) {
               ind = 0; // reset
               for (sclumps tempClump : spec[isp].head) {
                  if (tempClump.clumpid == iClump) {
                     temp = tempClump;
                     break;
                  }
                  ind++;
               }

               setClumpSpecAtPu(pu,SM,ipu,isp,iClump,thread);
               spec[isp].head[ind].head.push_back(ipu);

               // Remove old value for this clump
               ftemp = PartialPen4(isp,spec[isp].head[ind].amount,spec,clumptype);
               spec[isp].amount -= ftemp;
               spec[isp].occurrence -= spec[isp].head[ind].occs *(ftemp > 0);
               rAmount = returnAmountSpecAtPu(pu[ipu],SM,isp).second;
               spec[isp].head[ind].occs += (rAmount > 0);
               spec[isp].head[ind].amount += rAmount;
            }
            else {
               // TODO - revisit this branch logic. Removing LL implementation may have slowed things.
               // look at previous found index clump
               sclumps &pclump = spec[isp].head[ind];
               if (pclump.clumpid != iClump)
               {
                  // Check if this is a different clump
                  // Join this new clump to the old one
                  int ind2 = 0;

                  for (sclumps tempClump: spec[isp].head) {
                     if (tempClump.clumpid == iClump) {
                        break;
                     }
                     ind2++;
                  }

                  sclumps& pnewclump = spec[isp].head[ind2];
                  for (int puid: spec[isp].head[ind2].head) {
                     setClumpSpecAtPu(pu,SM,puid,isp,pclump.clumpid, thread);
                  }
                     
                  // cut out this clump and join it to pclump
                  // this operation copies over the pu
                  pnewclump.head.insert(pnewclump.head.end(), pclump.head.begin(), pclump.head.end());

                  // merge clumps
                  pclump.head = pnewclump.head; // copy assignment
                  pclump.amount += pnewclump.amount;
                  pclump.occs += pnewclump.occs;

                  ftemp = PartialPen4(isp,pnewclump.amount,spec,clumptype);
                  spec[isp].amount -= ftemp;
                  spec[isp].occurrence -= pnewclump.occs * (ftemp > 0);

                  // Remove pnewclump from spec[isp]
                  spec[isp].head.erase(spec[isp].head.begin() + ind2);

               } // Join the two clumps together
            }
         }
      }

      // Adding a New clump
      if (!ineighbours)
      {
         AddNewClump(isp,ipu,spec,pu,SM, thread);
         ftemp = PartialPen4(isp,rAmount,spec,clumptype);
         spec[isp].amount += ftemp;
         spec[isp].occurrence += (ftemp>0);
      } // Adding a new clump

      // Correcting Amount if new clump not added
      if (ineighbours)
      {
         ftemp = PartialPen4(isp,temp.amount,spec,clumptype);
         spec[isp].amount += ftemp;
         spec[isp].occurrence += temp.occs * (ftemp > 0);
      }
   }

   /************** REM PU ****************************************/
   /*********** Remove a planning unit. Note it is similar to CutClump but actually does action **/
   /**************************************************************/
   void RemPu(int ipu, int isp,vector<sconnections> &connections, vector<sspecies> &spec,vector<spustuff> &pu,
         vector<spu> &SM,int clumptype, int thread) {
      
      int ineighbours = 0;
      //struct slink{int id;struct slink *next;} *head = NULL, *newhead, *thead;
      //struct sclumps *oldclump,*pclump;
      sclumps pclump;
      //struct sclumppu *cppu,*ppu, *clumpcurr, *tppu;
      int cppu;
      struct sneighbour *pnbr;
      double oldamount,newamount = 0.0, rAmount;
      int newoccs;

       /* Find the correct clump to remove */
      bool found = false;
      int foundInd = 0;
      int idToFind =  rtnClumpSpecAtPu(pu,SM,ipu,isp, thread);
      for (sclumps tempClump : spec[isp].head) {
         if (tempClump.clumpid == idToFind) {
            found = true;
            break;
         }
         foundInd++;
      }

      sclumps& oldclump = spec[isp].head[foundInd];

      if (!found)
         displayErrorMessage("Serious error in Remove Type 4 species routine\n");

      /* Locate the correct clumppu */
      auto cppuIt = find(oldclump.head.begin(), oldclump.head.end(), ipu);
      cppu = *cppuIt;
      setClumpSpecAtPu(pu,SM,cppu,isp,0, thread);

      oldamount = PartialPen4(isp, oldclump.amount, spec, clumptype);
      spec[isp].amount -= oldamount;
      spec[isp].occurrence -= oldclump.occs * (oldamount > 0);

      /* Building the neighbour list */
      vector<int> head;
      for (sneighbour pnbr: connections[cppu].first) {
         ineighbours++;
         head.push_back(pnbr.nbr);
      }

      if (ineighbours <= 1) {
         rAmount = returnAmountSpecAtPu(pu[ipu],SM,isp).second;
         oldclump.amount -= rAmount;
         oldclump.occs -= (rAmount > 0);
         newamount = PartialPen4(isp,oldclump.amount,spec,clumptype);
         newoccs = oldclump.occs * (newamount > 0);

         // remove clummpu
         oldclump.head.erase(cppuIt);
         // remove old clump
         spec[isp].head.erase(spec[isp].head.begin() + foundInd);

         if (ineighbours < 1)
         {
            // remove old clump
            spec[isp].head.erase(spec[isp].head.begin() + foundInd);
            spec[isp].clumps--;
         } /* Removing redundant clump */

         spec[isp].amount += newamount;
         spec[isp].occurrence += newoccs;

         return;
      } /* I'm not cutting a clump */

      /* Else create new clumps */
      vector<int> modifiedHead = head; // copy assignment
      for (int id: head) {
         /* take first element as seed for new clump number */
         pclump = AddNewClump(isp,id,spec,pu,SM, thread);

         // TODO - revisit this logic
         /* Continue until you've added every conceivable neighbour */
         for (int id2: pclump.head) {
            for (sneighbour& pnbr: connections[id2].first) {
               if (rtnClumpSpecAtPu(pu,SM,pnbr.nbr,isp, thread) == oldclump.clumpid) {
                  auto eraseIt = find(oldclump.head.begin(), oldclump.head.end(), pnbr.nbr);
                  oldclump.head.erase(eraseIt);

                  setClumpSpecAtPu(pu,SM,id2,isp,pclump.clumpid, thread);
                  rAmount = returnAmountSpecAtPu(pu[id2],SM,isp).second;
                  pclump.amount += rAmount;
                  pclump.occs += (rAmount>0);
                  
                  /* Check if it is on neighbours list and if so then remove it from that list*/
                  auto idIt = find(modifiedHead.begin(), modifiedHead.end(), id2);
                  if (idIt != modifiedHead.end()) {
                     modifiedHead.erase(idIt);
                  }
               }
            }
         }

         spec[isp].amount += PartialPen4(isp,pclump.amount,spec,clumptype);
         spec[isp].occurrence += pclump.occs * (PartialPen4(isp,pclump.amount,spec,clumptype)>0);
      }

      /* Every neighbour in local list has been used and every clump formed*/
      /* Remove old clump */
      /* Worry about change in amount and hence score */
      spec[isp].head.erase(spec[isp].head.begin() + foundInd);
      ClearClump(isp,oldclump,pu,SM, thread);
   }

   /*** Remove Clump Check ***/
   /** returns 0 if any member of clump is non-removable, Ie status == 2 **/
   int RemClumpCheck(sclumps &pclump,vector<spustuff> &pu) {
      for (int pcpu: pclump.head) {
         if (pu[pcpu].status == 2)
            return 0;
      }

      return 1;
   }

   /********* Set Penalties for a given Type 4 Species ***/
   /* Returns 1 if the species is a 'bad species' and -1 if it is a 'good species' */
   /* Also sticks the penalty into spec[isp].penalty */
   int CalcPenaltyType4(int isp,int puno, vector<spu> &SM,vector<sconnections> &connections,
      vector<sspecies> &spec,vector<spustuff> &pu,double cm,int clumptype, int thread) {
      
      vector<sclumps> newno;
      int j,ipu,iputotal = 0;
      int ineighbours = 0,iclumpno,badspecies = 0;
      double totalamount,dummy = 0;
      int idummy;
      double cost = 0.0, connection = 0.0, rAmount = 0.0;

      vector<int> R;
      R.resize(puno); /* needed for separation */
      for (int i=0;i<puno;i++)
         R[i] = pu[i].status;


      /*** Step 1. Make a list of all the possible PUs to be included ****/
      vector<int> plisthead;
      for (int i=0;i<puno;i++) {
         if (returnAmountSpecAtPu(pu[i],SM,isp).second > 0)
         {
            if (pu[i].status==3) continue; /* not allowed to consider this one */
            if (pu[i].status==2)
            { /* add to clumps and remove from list */
               AddNewPU(i,isp,connections,spec,pu,SM,clumptype, thread);
               continue;
            } /* checking if PU forced into reserve */
            iputotal++;
            plisthead.push_back(i);
         } /** Made list of all sites with this species **/
      }

      /* Check first to see if I've already satisfied targets for this species */
      SpeciesAmounts4(isp,spec,clumptype);
      if (spec[isp].sepnum>0)
         spec[isp].separation = CountSeparation2(isp,0,newno,puno,R,pu,SM,spec,0,thread);
      if ((spec[isp].amount >= spec[isp].target) && (spec[isp].occurrence >= spec[isp].targetocc) && (spec[isp].separation >= spec[isp].sepnum))
      {
         spec[isp].amount = 0;
         spec[isp].occurrence = 0;
         spec[isp].separation = 0;

         /** Clean out all the clump numbers for this species.*/
         for (sclumps& c: spec[isp].head) {
            ClearClump(isp,c,pu,SM, thread);
         }
         spec[isp].clumps = 0;

         return(-1);
      }  /* return when all targets already met. */

      if (iputotal) {
         // shuffle vector to emulate randomly choosing
         shuffle(plisthead.begin(), plisthead.end(), rngEngine);

         do {
            int chosenPu = plisthead.back();
            plisthead.pop_back();

            // Add pu to system
            R[chosenPu] = 1;
            AddNewPU(chosenPu,isp,connections,spec,pu,SM,clumptype, thread);

            /*** Check to see if I should continue by calculating current holdings **/
            SpeciesAmounts4(isp,spec,clumptype);
            if (spec[isp].sepnum>0)
               spec[isp].separation = CountSeparation2(isp,0,newno,puno,R,pu,SM,spec,0, thread);

         } while ((spec[isp].amount < spec[isp].target || spec[isp].separation < spec[isp].sepnum || spec[isp].occurrence < spec[isp].targetocc)
                  && iputotal >= 1  );
      }

      if (spec[isp].amount < spec[isp].target || spec[isp].occurrence < spec[isp].targetocc)
      {
         badspecies = 1;
         displayProgress1("Species %i cannot be fully represented!\n", spec[isp].name);
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
         int count = 0;
         vector<sclumps> newClumps;
         while (count < spec[isp].head.size()) {
            sclumps& pclump = spec[isp].head[count];
            int i = 0;
            if (RemClumpCheck(pclump,pu)) {
               i = 1;
            }
            if (i)
            {
               if (spec[isp].occurrence - pclump.occs < spec[isp].targetocc)
                  i = 0;  
            } /* Check is occurrence decrease ok? */
            if (i)
            {
               if (!((spec[isp].amount - pclump.amount >= spec[isp].target) || (pclump.amount < spec[isp].target2)))
                  i = 0;
            } /* Check is amount decrease OK? */
            if (i && spec[isp].sepnum)
            {
               vector<sclumps> tempSepVector(spec[isp].head.begin() + count, spec[isp].head.end()); // In LL implementation, this would have been current point onwards.
               j = CountSeparation2(isp,0, tempSepVector,puno,R,pu,SM,spec,-1, thread); //TODO - revise vector that gets passed in here
               if ((j < spec[isp].separation) && (j < spec[isp].sepnum))
                  i = 0;

               if (!spec[isp].target2)
                  i = 0; /* cannot elegantly remove clumps if species is listed as non-clumping */
            }

            if (i)   /* This is a clump which can be safely removed */
            {  /* cut clump if uneccessary or it is too small */
               for (int puid: pclump.head) {
                  setClumpSpecAtPu(pu,SM,puid,isp,0, thread);
               }

               totalamount -= pclump.amount;
               spec[isp].clumps--;

               // physically remove clump
               spec[isp].head.erase(spec[isp].head.begin() + count);
            }
            else {
               count++; //increment only if we didn't remove
            }
         }
      } /*** Remove unneccesary clumps and links****/

      /** Test all PU's to see if any one of them are superfluous **/
      /* But only do this if occurrence target met. Otherwise every single pu is neccessary*/
      if (spec[isp].occurrence >= spec[isp].targetocc)
      {
         for (sclumps& pclump: spec[isp].head) {
            vector<int> newPu;
            for (int oldPu: pclump.head) {
               /** Test to see if this oldPu is necessary **/
               int i = 0;
               if (R[oldPu] != 2)
                  i = 1;
               if (i)
               {
                  rAmount = returnAmountSpecAtPu(pu[oldPu],SM,isp).second;
                  if ((pclump.amount - rAmount > spec[isp].target2) && (spec[isp].amount - rAmount > spec[isp].target))
                     i = 1;
                  else
                     i = 0;
               }  /* doesn't drop amount below min clump size or target */
               if (i)
               {
                  if (spec[isp].occurrence <= spec[isp].targetocc)
                     i = 0;
               } /* Does it drop occurrences too low? */
               if (i)
               {
                  vector<sclumps> temp;
                  sclumps pnewclump = {oldPu, 0, 0};
                  temp.push_back(pnewclump);
                  j = CountSeparation2(isp,oldPu,temp,puno,R,pu,SM,spec,-1, thread);
                  if ((j < spec[isp].separation) && (j < spec[isp].sepnum))
                     i = 0;
               } /* How does count sep fare? */
               if (i)
               {
                  if (ClumpCut(isp,pu,spec,pclump,oldPu,connections,SM, dummy, idummy, j,0,clumptype, thread))
                     i = 0;
               } /* Does it cut the clump? these are not allowed to remove */
               if (i)  /* Is this removable? */
               {        /* remove pcpu */
                  setClumpSpecAtPu(pu,SM,oldPu,isp,0, thread);
                  totalamount -= rAmount;
                  pclump.amount -= rAmount;
                  // pu effectively removed by not including in newPu
               }  /** remove unneccessary clumppu **/
               else {
                  newPu.push_back(oldPu); // add to new list so we include it.
               }
            }
            pclump.head = newPu;
         }
      } /** Cycle over each pclump **/

      /*** Now count the cost of this particular reserve ****/
      /*** For each clump figure out connection cost ***/
      for (sclumps& pclump: spec[isp].head) { 
         iclumpno = pclump.clumpid;
         for (int oldPu: pclump.head) {
            if (pu[oldPu].status != 2)
               {
                  cost += pu[oldPu].cost;
                  connection += connections[oldPu].fixedcost;
               } /* only count fixed costs if PU not forced into reserve */
               if (connections[oldPu].nbrno)
               {
                  for (sneighbour& pnbr: connections[oldPu].first)
                  {
                     if (rtnClumpSpecAtPu(pu,SM,pnbr.nbr,isp,thread) != iclumpno)
                        connection += pnbr.cost;
                  } /** Counting each individual connection **/
            } /** Counting connection strength if neccessary **/
         }
      }

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
      for (sclumps& clump: spec[isp].head)
      {
         ClearClump(isp,clump,pu,SM,thread);
      }  /** Remove each clump ***/
      spec[isp].clumps = 0;

      return(badspecies);
   }

   /**** Partial Penalty for type 4 species ***/
   double PartialPen4(int isp,double amount, vector<sspecies> &spec,int clumptype) {
      if (amount >= spec[isp].target2)
      {
         return (amount);    /* I'm not a partial penalty */
      }
      else
      {
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
      }
   }

   /*** Value for Adding a Planning Unit ****/
   double ValueAdd(int isp,int ipu,int puno, vector<int> &R,vector<sconnections> &connections,vector<spustuff> &pu,
                  vector<spu> &SM,vector<sspecies> &spec,int clumptype, int thread) {

      int ineighbours = 0,iclumpid,iseparation;
      vector<sneighbour> pnbr;
      vector<sclumps> head;
      vector<sclumps> sepclump;
      double amount,oldamount = 0.0,shortfall;
      int oldoccs = 0,occs, iClump;

      /* Count neighbours */
      if (connections[ipu].nbrno > 0) {
         pnbr = connections[ipu].first;
         for (sneighbour& neighbour: pnbr) {
            iClump = rtnClumpSpecAtPu(pu,SM,neighbour.nbr,isp, thread);
            if (iClump) {
               iclumpid = 1;
               
               /* Is nbr on my list ?*/
               if (!head.empty()) {
                  for (sclumps tempClump: head) {
                     if (tempClump.clumpid == iClump)
                        iclumpid = 0;
                  }
               }

               if (iclumpid) {
                  ineighbours++;

                  sclumps correctClump;
                  /* find the right clump */
                  for (sclumps tempClump2: spec[isp].head) {
                     if (tempClump2.clumpid == iClump) {
                        correctClump = tempClump2;
                     }
                  }

                  head.push_back(correctClump);

                  if (spec[isp].sepnum) {
                     for (int ppu: correctClump.head) {
                        sclumps newClump;
                        newClump.clumpid = ppu;
                        sepclump.push_back(newClump); /* glue to sep list. Still need amount */
                     } /* stick whole clump onto separation clump for later */
                  } 
               }
            }
         }
      }  /* If There are neighbours */

      if (spec[isp].sepnum) {
         sclumps psclump;
         psclump.clumpid = ipu;
         sepclump.push_back(psclump);
      } /* Add ipu to my sepclump list */

      /* now I know number and names of neighbouring clumps */
      amount = returnAmountSpecAtPu(pu[ipu],SM,isp).second;
      occs = (amount > 0);
      for(sclumps plink: head) {
            amount += plink.amount;
            occs += plink.occs;
            oldamount += PartialPen4(isp,plink.amount,spec,clumptype);
            oldoccs += plink.occs * (PartialPen4(isp,plink.amount,spec,clumptype)>0);
      }

      /* set the sepclump amounts to this new amount */
      if (spec[isp].sepnum)
      {
         for (auto& psclump: sepclump) {
            psclump.amount = amount;
         }
      }

      amount = PartialPen4(isp,amount,spec,clumptype);
      occs = occs * (amount > 0);
      amount = amount - oldamount; /* amount is change in amount for this species */
      occs = occs - oldoccs;

      if (spec[isp].sepnum) {
         iseparation = CountSeparation2(isp,0,sepclump,puno,R,pu,SM,spec,1, thread);  /* imode = 1 doesn't do anything*/
      }

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
      return(shortfall + computeSepPenalty(iseparation,spec[isp].sepnum));
   }

   /** Value Remove. The amount of species loss for removing a single pu */
   double ValueRem(int ipu,int isp,vector<sspecies> &spec,vector<sconnections> &connections,
                  vector<spustuff> &pu,vector<spu> &SM,int clumptype, int thread) {
      double newamount = 0,amount,shortfall=0;
      sclumps pclump;
      int ppu;
      int iseparation;
      int newocc = 0;

      /* locate the clump and clumppu of the target site ipu */
      for (sclumps clumpTemp : spec[isp].head) {
         if (clumpTemp.clumpid == rtnClumpSpecAtPu(pu,SM,ipu,isp, thread)) {
            pclump = clumpTemp;
            break;
         }
      }

      /* locate the correct pclump pu */
      for (int id : pclump.head) {
         if (id == ipu) {
            ppu = id;
            break;
         }
      }

      /* debugmem2 = debugmem; debugging line */
      if (spec[isp].sepnum)
         ClumpCut(isp,pu,spec,pclump,ppu,connections,SM, newamount, newocc, iseparation,1, clumptype, thread);
      else
         ClumpCut(isp,pu,spec,pclump,ppu,connections,SM, newamount, newocc, iseparation,0, clumptype, thread);

      if (spec[isp].target)
      {
         amount = spec[isp].amount + newamount - PartialPen4(isp, pclump.amount,spec,clumptype) ;
         shortfall = amount > spec[isp].target ? 0 : (spec[isp].target - amount)/spec[isp].target;
      }  /* there is an abundance amount */

      /*if (isp == 16) printf("pclump->occs %i targetocc %i shortfall %.2f\n",
                                       pclump->occs,spec[isp].targetocc,shortfall);*/
      if (spec[isp].targetocc) {  /* Handle the case where there is a targetocc */
         amount = spec[isp].occurrence +newocc - pclump.occs * (PartialPen4(isp,pclump.amount,spec,clumptype)>0);
         if (amount < spec[isp].targetocc)
               shortfall += ((double) spec[isp].targetocc - amount)/(double) spec[isp].targetocc;
         if (spec[isp].target)
               shortfall /= 2;
         }
      /*  if (isp ==16) printf("shortfall %.2f occ %i newocc %i pclump->amount %.2f\n",
               shortfall, spec[isp].occurrence,newocc,pclump->amount);*/

      return(shortfall + computeSepPenalty(iseparation,spec[isp].sepnum));
   }

   /***************   NewPenalty4   *********************/
   /* Calculates the new penalty for adding or removing a PU for species which have
      clumping requirements */  
   double NewPenalty4(int ipu,int isp,int puno,vector<sspecies> &spec,vector<spustuff> &pu,vector<spu> &SM,
               vector<int> &R,vector<sconnections> &connections,int imode,int clumptype, int thread) {
      double amount;

      if (imode == 1) {
         if (spec[isp].penalty == 0)
               return (0);  /* Targets have all already been met */
         amount = ValueAdd(isp,ipu,puno,R,connections,pu,SM,spec,clumptype, thread);
      }
      else {
            /* determine change in this amount */
         amount = ValueRem(ipu,isp,spec,connections,pu,SM,clumptype, thread);
      } /** removing a planning unit **/
      return amount;
   }

   int ValidPU(int ipu,int isp, vector<sclumps> &newno,vector<sspecies> &spec,vector<spustuff> &pu,
            vector<spu> &SM,int imode, int thread) {

      // Returns true if ipu is acceptable as a planning unit
      int i =  returnAmountSpecAtPu(pu[ipu],SM,isp).first;
      sclumps pclump;
      bool found = false;

      if (!newno.empty())
      {
         if (imode == -2)
         {
            if (SM[i].clump[thread] == newno[0].clumpid)
               return(0); // This whole clump is to be removed
         }
         
         // iterate through newno to find if there's any clumpId equal to given ipu.
         // if there is check if amount meets target. 
         for (sclumps clump: newno) {
            if (ipu == clump.clumpid) {
               if (clump.amount < spec[isp].target2) {
                  return 0;
               }
               else {
                  return 1;
               }
            }
         } // ipu is on list of changed pus
      }
      
      // Find clump
      for (sclumps clump: spec[isp].head) {
         if (SM[i].clump[thread] == clump.clumpid) {
            pclump = clump;
            found = true;
         }
      } // scan through to find clump

      if (found)
      {
         if (pclump.amount <spec[isp].target2)
            return 0;
         else
            return 1;
      }
      else
      {
         if (SM[i].amount < spec[isp].target2)
            return 0;
         else
            return 1;
      }
   }

   bool CheckDistance(int i, int j, vector<spustuff> &pu, double squaretarget) {
      // compare x1*x2+y1*y2 with squaretarget
      if ((pu[i].xloc-pu[j].xloc)*(pu[i].xloc-pu[j].xloc) + (pu[i].yloc-pu[j].yloc)* (pu[i].yloc-pu[j].yloc) >= squaretarget)
         return true;
      else
         return false;
   } // Is Distant returns true if PU's are big enough distance apart

   int CountSeparation(int isp, vector<sclumps> &newno,
   vector<spustuff> &pu,vector<spu> &SM,vector<sspecies> &spec,int imode, int thread)
   {
      // imode 0 = count separation on current
      // imode 1 = count separation if ipu were included
      // imode -1 = count separation if ipu were excluded
      // The following assumes imode = 0 for starters

      vector<int> ptemp;
      int sepcount = 1;
      double targetdist = spec[isp].sepdistance * spec[isp].sepdistance;

      if (targetdist == 0)
         return 3; // Shortcut if sep not apply to this species. This assumes that 3 is highest sep count
         
      // Set up the first list
      if (imode == 1) 
      {
         if (newno.empty()) {
            throw runtime_error("CountSeparation Error: newno is empty and imode==1");
         }

         if (ValidPU(newno[0].clumpid,isp,newno,spec,pu,SM,imode, thread))
         {
            ptemp.push_back(newno[0].clumpid);
         }
      }

      for (sclumps pclump : spec[isp].head) {
         for (int ppu : pclump.head) {
            if (ValidPU(ppu,isp,newno,spec,pu,SM,imode, thread)) {
               ptemp.push_back(ppu);
            }
         } // Add all valid species bearing PU's to list
      } // need to worry about added pu which is not on spec[isp].head list

      // TODO - revisit this logic to make sure
      // Check distance between each combination of ids until target dist met.
      for (int i = 0; i < ptemp.size(); i++) {
         for (int j = i; j < ptemp.size(); j++) {
            if (CheckDistance(ptemp[i],ptemp[j],pu,targetdist))
            {
               for (int k = j; k < ptemp.size(); k++) {
                  if (CheckDistance(ptemp[j], ptemp[k], pu, targetdist))
                  {
                     return 3;
                  } // I have succeeded in finding what I'm looking for

                  sepcount = 2; // there is a separation of at least 2. This should be used to cut down calls to this function
               }
            }
         }
      }

      return(sepcount);
   } // CountSeparation

   /* This is a modified form of count separation where the user can specify any
    maximum separation distance rather than just assuming a sep distance of three */
   /* ipu and newno used when imode <> 0. When counting as if ipu were added or removed
      ipu used for non-clumping and newno for clumping species */
   int CountSeparation2(int isp,int ipu, vector<sclumps> &newno,int puno,vector<int> &R,
      vector<spustuff> &pu,vector<spu> &SM,vector<sspecies> &spec,int imode, int thread) {

      vector<sseplist> Dist;
      vector<int> head;
      int sepcount,bestsep = 0,i,currcol;
      double targetdist;
      
      targetdist = spec[isp].sepdistance * spec[isp].sepdistance;

      if (targetdist == 0)
         return(spec[isp].sepnum); // Shortcut if sep not apply to this species

      // Set up array for counting separation
      Dist.resize(spec[isp].sepnum);

      // First scan through sites. Grab first valid and place rest in lists
      head = makelist(isp,ipu,puno,R,newno,spec,pu,SM,imode, thread);

      if (head.empty())
      {
         return(0);
      } // There was nothing to put in the list


      Dist[0].head = head;
      Dist[0].size = 1;

      //head = head->next;
      head.erase(head.begin());

      if (head.empty())
      {
         return(1);
      }  // There was only one item in the list

      // Deal out link list
      sepcount = SepDealList(head,Dist,pu,spec, Dist[0].head[0], 0, targetdist, isp);
      if (sepcount >= spec[isp].sepnum-1)
      {
         return(spec[isp].sepnum);
      } // I'm at maximum separation
      
      bestsep = sepcount;

      do
      {
         // Note - I have no clue what this is intended to do.
         // The main Loop
         for (currcol=0;!Dist[currcol+1].head.empty() && currcol < spec[isp].sepnum-2;currcol++)
            ;
            
         if (currcol == 0)
         {
            if (Dist[0].size < spec[isp].sepnum)
            {
               return(bestsep + 1);
            } // cannot increase separation terminate function
            else
            {
               Dist[0].head.erase(Dist[0].head.begin());
               head = Dist[0].head;
               Dist[0].size = 1;
               sepcount = SepDealList(head,Dist,pu,spec,Dist[0].head[0],0,targetdist,isp);
            }
         } // Deal with first column
         else
         {
            if (Dist[currcol].size + currcol  < spec[isp].sepnum)
            {
               Dist[currcol-1].size += Dist[currcol].size;
               Dist[currcol].head.clear();
               Dist[currcol].size = 0;
               sepcount = 0;
            } // list is not long enough to increase sepcount
            else
            {
               Dist[currcol-1].size++;
               Dist[currcol].head.erase(Dist[currcol].head.begin());
               head = Dist[currcol].head;
               Dist[currcol].size = 1;
               sepcount = SepDealList(head,Dist,pu,spec,Dist[currcol].head[0],currcol,targetdist,isp);
            } // else this column might be long enough
         } // Deal with columns other than the first
         
         if (sepcount > bestsep)
         {
            bestsep = sepcount;
         }

      } while(bestsep < spec[isp].sepnum-1); // Main loop.

      return(bestsep+1);
   }

   // Make List
   // This makes a list of all the valid PUs which occur on the reserve and on which species
   // isp is present (or NULL), in the form of a slink link list
   vector<int> makelist(int isp,int ipu, int puno,vector<int> &R, vector<sclumps> &newno,vector<sspecies> &spec,
                       vector<spustuff> &pu,vector<spu> &SM,int imode, int thread) {

      // Note! Changes made so that we append to end instead of inserting to beginning of linked list. 
      vector<int> head;
      if (spec[isp].target2)
      {
         // deal with clumping species differently from non-clumping
         int tempClumpId;
         if (imode == 1)
         {
            if (newno.empty()) {
               throw runtime_error("Error makelist - newno is null when imode==1.");
            }

            if (ValidPU(newno[0].clumpid,isp,newno,spec,pu,SM,imode, thread))
            {
               head.push_back(newno[0].clumpid);
            }
         }

         for (sclumps& pclump: spec[isp].head) {
            for (int ppu: pclump.head) {
               if (ValidPU(ppu,isp,newno,spec,pu,SM,imode, thread))
               {
                  head.push_back(ppu);
               }  // Add all valid species bearing PU's to list
            }
         }
      } // if target2
      else
      {
         double rAmount = returnAmountSpecAtPu(pu[ipu],SM,isp).second;
         // non clumping species
         if ((imode ==1) && rAmount)
         {
            head.push_back(ipu);
         } // deal with imode == 1 case
         
         for (int i=0;i<puno;i++)
            if (((R[i] == 1 || R[i] == 2) && rAmount) && !(imode == -1 && ipu == i))
            {
               head.push_back(i);
            }
      } // non clumping species

      return(head);
   }

   // Sep Deal List
   // This funciton is called by count separation2. It takes a link list of sites and 'deals' them
   //  out on to the seplist
   int SepDealList(vector<int> &head, vector<sseplist> &Dist, vector<spustuff> &pu,
      vector<sspecies> &spec,int first,int sepnum,double targetdist,int isp) {
      
      // Currsep is the current separation maximum it is 0 up to sepnum
      // first is only needed if maximum is at 0, sepnum is the target separation
      int placefound,currtarget,bestsep=0;
      int currsep;
      struct slink *temp;

      int i = 0;
      for (int id: head) {
         placefound = 0;
         currtarget = first;
         currsep = sepnum;

         do
         {
            if (CheckDistance(id,currtarget,pu,targetdist))
            {
               currsep++;
               if (currsep == spec[isp].sepnum-1)
               {
                  for (int j = i; j < head.size(); j++) {
                     Dist[currsep].head.push_back(head[j]);
                  } // glue remaining list on to bottom of Dist. ignoring size and tail as useless
               
                  return(currsep);
               } // Just found valid separation
               
               if (!Dist[currsep].head.empty())
                  currtarget = Dist[currsep].head[0];
               else
               {
                  placefound = 1;
                  Dist[currsep].head.push_back(head[i]);
                  Dist[currsep].size++;
               } // I'm at the end of the line
            } // Good distance
            else
            {
               placefound = 1;
               Dist[currsep].head.push_back(head[i]);
               Dist[currsep].size++;
            } // bad distance
         } while (!placefound); // Doing each individual

         if (currsep > bestsep)
         {
            bestsep = currsep;
         }

         i+=1;
      }

      return(bestsep);
   }

   // ********* Connection Cost Type 2 **************
   // **  Requires R[]. imode2 = 0 there is no negative cost for removing connection, we are calling from ReserveCost
   //                         or 1 there is a negative cost for removing connection, we are calling from Annealing
   //                   imode = -1 we are removing the planning unit from a reserve, calling from Annealing
   //                        or 1  we are adding the planning unit to a reserve, or it is already in reserve
   double ConnectionCost2(int ipu, vector<sconnections> &connections, vector<int> &R, int imode, int imode2, double cm) {
      double fcost, rDelta;
      int R_pu1,i;

      #ifdef DEBUG_CONNECTIONCOST2
      if (asymmetricconnectivity)
         if (imode2)
         {
            appendTraceFile("ConnectionCost2 start puid %i imode %i imode2 %i\n",pu[ipu].id,imode,imode2);
         }
      #endif

      fcost = connections[ipu].fixedcost*imode;

      if (asymmetricconnectivity)
      {
         for(sneighbour& p: connections[ipu].first) {
            if (imode2) // calling from Annealing
            {
               #ifdef DEBUG_CONNECTIONCOST2
               rDelta = 0;
               #endif

               if (imode == 1)
                  R_pu1 = 0;
               else
                  R_pu1 = 1;

               if (p.connectionorigon)
               {
                  if (R[p.nbr] == 0)
                  {
                     if (R_pu1 == 1)
                     {
                        rDelta = -1*p.cost;
                        fcost += rDelta;
                     }
                     else
                     {
                        rDelta = p.cost;
                        fcost += rDelta;
                     }
                  }
               }
               else
               {
                  if (R[p.nbr] == 1 || R[p.nbr] == 2)
                  {
                     if (R_pu1 == 1)
                     {
                        rDelta = p.cost;
                        fcost += rDelta;
                     }
                     else
                     {
                        rDelta = -1*p.cost;
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
               if (R[p.nbr] == 0)
                  if (p.connectionorigon)
                  {
                     rDelta = p.cost;
                     fcost += rDelta;
                  }
            }
         }
      }
      else
      {
         for (sneighbour &p : connections[ipu].first) // treatment for symmetric connectivity
         {
            if (fOptimiseConnectivityIn == 1)
            { // optimise for "Connectivity In"
               if (R[p.nbr] == 1 || R[p.nbr] == 2)
               {
                  rDelta = imode * p.cost;
                  fcost += rDelta;
               }
               else
               {
                  rDelta = imode * imode2 * p.cost * -1;
                  fcost += rDelta;
               }
            }
            else
            { // optimise for "Connectivity Edge"
               if (R[p.nbr] == 1 || R[p.nbr] == 2)
               {
                  rDelta = imode * imode2 * p.cost * -1;
                  fcost += rDelta;
               }
               else
               {
                  rDelta = imode * p.cost;
                  fcost += rDelta;
               }
            }
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
   }

} // marxan