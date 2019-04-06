// read the planning units file pu.dat
int readPlanningUnits(int *puno,struct spustuff *pu[],struct sfname fnames)
{
    FILE *fp;
    struct pustruct{int id;double cost; int status; double xloc; double yloc; double prob;};
    char *readname;
    char sLine[600];
    char *varlist[6] = {"id","cost","status","xloc","yloc","prob"};
    int numvars = 6,ivars,i=0,j;
    char *sVarName,*sVarVal;
    struct snlink *head = NULL, *temp = NULL;
    struct spustuff putemp;
    struct spulink{struct spustuff stemp;struct spulink *next;} *spuhead = NULL, *newspulink;

    readname = (char *) calloc(strlen(fnames.puname) + strlen(fnames.inputdir)+2, sizeof(char));

    #ifdef MEMDEBUG
    iMemoryUsed += (strlen(fnames.puname) + strlen(fnames.inputdir)+2) * sizeof(char);
    displayProgress1("memory used %i\n",iMemoryUsed);
    #endif

    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.puname);
    if((fp = fopen(readname,"r"))==NULL)
        displayErrorMessage("Planning Unit file %s has not been found.\nAborting Program.",readname);
    free(readname);

    // Scan header
    fgets(sLine,500-1,fp);

    sVarName = strtok(sLine," ,;:^*\"/\t\'\\\n");
    head = GetVarName(varlist,numvars,sVarName,head,fnames.puname);

    ivars = 1;
    temp = head;
    while ((sVarName = strtok(NULL," ,;:^*\"/\t\'\\\n")) != NULL)
    {
          ivars++;
          temp->next = GetVarName(varlist,numvars,sVarName,head,fnames.puname);
          temp = temp->next;
    }  // tokking out all the variable names from the header line. There are numVars of them.

    /* While there are still lines left feed information into temporary link list */
    while (fgets(sLine,500-1,fp))
    {
          i++;
          putemp.id = -1; /* Set defaults for any missing values */
          putemp.cost = 1;
          putemp.status = 0;
          putemp.xloc = -1;
          putemp.yloc = -1;
          putemp.prob = 0;
          for (temp = head;temp;temp = temp->next)
          {
              if (temp == head)
                 sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
              else
                  sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
              if (strcmp("id",temp->name)==0)
              {
                  sscanf(sVarVal,"%d",&putemp.id);
              }
              else
                  if (strcmp("status",temp->name)==0)
                  {
                     sscanf(sVarVal,"%d",&putemp.status);
                  }
                  else
                      if (strcmp("xloc",temp->name)==0)
                      {
                         sscanf(sVarVal,"%lf",&putemp.xloc);
                      }
                      else
                          if (strcmp("yloc",temp->name)==0)
                          {
                             sscanf(sVarVal,"%lf",&putemp.yloc);
                          }
                          else
                              if (strcmp("cost",temp->name)==0)
                              {
                                 sscanf(sVarVal,"%lf",&putemp.cost);
                              }
                              else
                                  if (strcmp("prob",temp->name)==0)
                                  {
                                     sscanf(sVarVal,"%lf",&putemp.prob);
                                     iProbFieldPresent = 1;
                                  }

          } /* looking for ivar different input variables */

          if (putemp.id == -1)
             displayErrorMessage("ERROR: Missing planning unit id for line %d. \n",i);
          /* Stick everything from putemp into link list */
          newspulink = (struct spulink *) malloc(sizeof(struct spulink));
          newspulink->stemp.id = putemp.id;
          newspulink->stemp.status = putemp.status;
          newspulink->stemp.cost = putemp.cost;
          newspulink->stemp.xloc = putemp.xloc;
          newspulink->stemp.yloc = putemp.yloc;
          newspulink->stemp.prob = putemp.prob;
          newspulink->next = spuhead;
          spuhead = newspulink;

    } /* while still lines in data file */

    fclose(fp);

    /* Create array to store the information */
    *puno = i;
    *pu = (struct spustuff *) calloc(*puno,sizeof(struct spustuff));

    #ifdef MEMDEBUG
    iMemoryUsed += (*puno) * sizeof(struct spustuff);
    displayProgress1("memory used %i\n",iMemoryUsed);
    #endif

    for (i=*puno-1;i>=0;i--)
    {
        (* pu)[i].id = spuhead->stemp.id;
        (* pu)[i].cost = spuhead->stemp.cost;
        (* pu)[i].status = spuhead->stemp.status;
        (* pu)[i].xloc = spuhead->stemp.xloc;
        (* pu)[i].yloc = spuhead->stemp.yloc;
        (* pu)[i].prob = spuhead->stemp.prob;
        (* pu)[i].richness = 0;
        (* pu)[i].offset = 0;
        (* pu)[i].probrichness = 0;
        (* pu)[i].proboffset = 0;

        newspulink = spuhead;
        spuhead = spuhead->next;
        free(newspulink);
    }

    if (iProbFieldPresent == 1)
    {
       fProb1D = 1;
    }

    return(i);
} // readPlanningUnits

// read species file: spec.dat
int readSpecies(int *spno,struct sspecies *spec[],struct sfname fnames)
{
    FILE *fp;
    int n=0;
    struct snlink *snhead= NULL,*temp;
    struct slink{struct sspecies stemp;struct slink *next;} *head = NULL,*newlink;
    struct sspecies spectemp;
    char *readname;
    char speciesname[255];
    char sLine[500];
    char *varlist[12] = {"id","type","target","spf",
                         "target2","sepdistance","sepnum","name",
                         "targetocc","prop","ptarget1d","ptarget2d"};
    int numvars = 12,ivars;
    char *sVarName,*sVarVal;

    readname = (char *) calloc(strlen(fnames.specname) + strlen(fnames.inputdir)+2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.specname);
    if((fp = fopen(readname,"r"))==NULL)
        displayErrorMessage("Species file %s has not been found.\nAborting Program.",readname);
    free(readname);


    // Scan header
    fgets(sLine,500-1,fp);

    sVarName = strtok(sLine," ,;:^*\"/\t\'\\\n");
    snhead = GetVarName(varlist,numvars,sVarName,snhead,fnames.specname);
    ivars = 1;
    temp = snhead;
    while ((sVarName = strtok(NULL," ,;:^*\"/\t\'\\\n")) != NULL) {
        ivars++;
        temp->next = GetVarName(varlist,numvars,sVarName,snhead,fnames.specname);
        temp = temp->next;
    }  // tokking out all the variable names from the header line. There are numVars of them

    // While there are still lines left feed information into temporary link list
    while (fgets(sLine,500-1,fp))
    {
          n++;
          // Clear important species stats
          spectemp.name = -1;
          spectemp.target = -1;
          spectemp.type = -1;
          spectemp.spf = -1;
          spectemp.target2 = -1;
          spectemp.targetocc = -1;
          spectemp.sepdistance = -1;
          spectemp.sepnum = -1;
          spectemp.prop = -1;
          spectemp.ptarget1d = -1;
          spectemp.ptarget2d = -1;
          spectemp.rUserPenalty = -1;
          speciesname[0] = '\0';
          for (temp = snhead;temp;temp = temp->next)
          {
              if (temp == snhead)
                 sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
              else
                  sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");

              if (strcmp("id",temp->name)==0)
              {
                 sscanf(sVarVal,"%d",&spectemp.name);
              }
              else
                  if (strcmp("type",temp->name)==0)
                  {
                     sscanf(sVarVal,"%d",&spectemp.type);
                  }
                  else
                      if (strcmp("target",temp->name)==0)
                      {
                         sscanf(sVarVal,"%lf",&spectemp.target);
                      }
                      else
                          if (strcmp("prop",temp->name)==0)
                          {
                             sscanf(sVarVal,"%lf",&spectemp.prop);
                             if (spectemp.prop > 0)
                                fSpecPROPLoaded = 1;
                          }
                          else
                              if (strcmp("spf",temp->name)==0)
                              {
                                 sscanf(sVarVal,"%lf",&spectemp.spf);
                              }
                              else
                                  if (strcmp("sepdistance",temp->name)==0)
                                  {
                                     sscanf(sVarVal,"%lf",&spectemp.sepdistance);
                                  }
                                  else
                                      if (strcmp("sepnum",temp->name)==0)
                                      {
                                         sscanf(sVarVal,"%d",&spectemp.sepnum);
                                      }
                                      else
                                          if (strcmp("target2",temp->name)==0)
                                          {
                                             sscanf(sVarVal,"%lf",&spectemp.target2);
                                          }
                                          else
                                              if (strcmp("targetocc",temp->name)==0)
                                              {
                                                 sscanf(sVarVal,"%d",&spectemp.targetocc);
                                              }
                                              else
                                                  if (strcmp("ptarget1d",temp->name)==0)
                                                  {
                                                     sscanf(sVarVal,"%lf",&spectemp.ptarget1d);
                                                  }
                                                  else
                                                      if (strcmp("ptarget2d",temp->name)==0)
                                                      {
                                                         sscanf(sVarVal,"%lf",&spectemp.ptarget2d);
                                                      }
        } // looking for ivar different input variables
        newlink = (struct slink *) malloc(sizeof(struct slink));
        newlink->stemp.name = spectemp.name;
        newlink->stemp.target = spectemp.target;
        newlink->stemp.prop = spectemp.prop;
        newlink->stemp.spf = spectemp.spf;
        newlink->stemp.type = spectemp.type;
        newlink->stemp.targetocc = spectemp.targetocc;
        newlink->stemp.rUserPenalty = spectemp.rUserPenalty;
        newlink->stemp.target2 = spectemp.target2;
        newlink->stemp.sepdistance = spectemp.sepdistance;
        newlink->stemp.sepnum = spectemp.sepnum;
        newlink->stemp.sname = (char *) calloc(strlen(speciesname)+1,sizeof(char));
        strcpy(newlink->stemp.sname,speciesname);
        newlink->stemp.ptarget1d = spectemp.ptarget1d;
        newlink->stemp.ptarget2d = spectemp.ptarget2d;
        newlink->next = head;
        head = newlink;
    } // Scanning through each line of file

    fclose(fp);

    // Now do as Name.dat in forming species array
    *spno = n;
    *spec = (struct sspecies *) calloc(*spno,sizeof(struct sspecies));
    // put each link into namelist and free it
    n = 0;
    //n = *spno-1;
    while (head)
    {
          (* spec)[n].name = head->stemp.name;
          (* spec)[n].type = head->stemp.type;
          (* spec)[n].target = head->stemp.target;
          (* spec)[n].prop = head->stemp.prop;
          (* spec)[n].spf = head->stemp.spf;
          (* spec)[n].target2 = head->stemp.target2;
          (* spec)[n].targetocc = head->stemp.targetocc;
          (* spec)[n].sepdistance = head->stemp.sepdistance;
          (* spec)[n].rUserPenalty = head->stemp.rUserPenalty;
          (* spec)[n].sepnum = head->stemp.sepnum;
          (* spec)[n].richness = 0;
          (* spec)[n].probability1D = 0;
          (* spec)[n].probability2D = 0;
          (* spec)[n].Zscore1D = 0;
          (* spec)[n].Zscore2D = 0;
          (* spec)[n].ptarget1d = head->stemp.ptarget1d;
          (* spec)[n].ptarget2d = head->stemp.ptarget2d;
          (* spec)[n].sname = (char *) calloc(strlen(head->stemp.sname)+1,sizeof(char));
          strcpy((* spec)[n].sname,head->stemp.sname);
          n++;
          newlink = head;
          head = head->next;
          free(newlink->stemp.sname);
          free(newlink);
    }

    return(n);
}  // readSpecies

// read species block definition file
// allows attribute to be applied to species based upon a type
int readSpeciesBlockDefinition(int *gspno,struct sgenspec *gspec[],struct sfname fnames)
{
    FILE *fp;
    char *readname;
    char sLine[500];
    char *varlist[8] = {"type","target","target2","targetocc",
                        "sepnum","sepdistance","prop","spf"};
    int numvars = 8,ivars,i=0;
    char *sVarName,*sVarVal;
    struct snlink *head = NULL, *temp = NULL;
    struct sgenspec gstemp;
    struct sgslink{struct sgenspec stemp;struct sgslink *next;} *gshead = NULL, *newgslink;

    /* Find and Open File */
    readname = (char *) calloc(strlen(fnames.blockdefname) + strlen(fnames.inputdir)+2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.blockdefname);
    if((fp = fopen(readname,"r"))==NULL) {
        displayWarningMessage("Warning: Block Definition File %s not found.\n",fnames.blockdefname);
        free(readname);
        return(0);
    }
    free(readname);

    /* Scan header */
    fgets(sLine,500-1,fp);
    sVarName = strtok(sLine," ,;:^*\"/\t\'\\\n");
    head = GetVarName(varlist,numvars,sVarName,head,fnames.blockdefname);
    ivars = 1;
    temp = head;
    while ((sVarName = strtok(NULL," ,;:^*\"/\t\'\\\n")) != NULL) {
        ivars++;
        temp->next = GetVarName(varlist,numvars,sVarName,head,fnames.blockdefname);
        temp = temp->next;
    }  /* tokking out all the variable names from the header line. There are numVars of them*/
    /* While there are still lines left feed information into temporary link list */
    while (fgets(sLine,500-1,fp)) {
        i++;
        gstemp.type = -1; /* Set defaults for any missing values */
        gstemp.targetocc = -1;
        gstemp.target = -1;
        gstemp.target2 = -1;
        gstemp.sepnum = -1;
        gstemp.sepdistance = -1;
        gstemp.prop = -1;
        for (temp = head;temp->next;temp = temp->next)
         {
             if (temp == head)
                   sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
            else
                sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");

            if (strcmp("type",temp->name)==0) {
                sscanf(sVarVal,"%d",&gstemp.type);
            }
            else if (strcmp("targetocc",temp->name)==0){
                sscanf(sVarVal,"%d",&gstemp.targetocc);
            }
            else if (strcmp("target",temp->name)==0){
                sscanf(sVarVal,"%lf",&gstemp.target);
            }
            else if (strcmp("target2",temp->name)==0){
                sscanf(sVarVal,"%lf",&gstemp.target2);
            }
            else if (strcmp("sepnum",temp->name)==0){
                sscanf(sVarVal,"%d",&gstemp.sepnum);
            }
            else if (strcmp("sepdistance",temp->name)==0){
                sscanf(sVarVal,"%lf",&gstemp.sepdistance);
            }
            else if (strcmp("prop",temp->name)==0){
                sscanf(sVarVal,"%lf",&gstemp.prop);
            }
            else if (strcmp("spf",temp->name)==0){
                sscanf(sVarVal,"%lf",&gstemp.spf);
            }
            else 
            {
                displayWarningMessage("Cannot find >%s< \n",temp->name);
                displayErrorMessage("Serious error in GenSpecies data reading function.\n");
            }

        } /* looking for ivar different input variables */

        if (gstemp.type== -1)
            displayErrorMessage("ERROR: Missing Gen Species type for line %d. \n",i);
        /* Stick everything from gstemp into link list */
        newgslink = (struct sgslink *) malloc(sizeof(struct sgslink));
        newgslink->stemp.type = gstemp.type;
        newgslink->stemp.targetocc = gstemp.targetocc;
        newgslink->stemp.target = gstemp.target;
        newgslink->stemp.target2 = gstemp.target2;
        newgslink->stemp.sepnum = gstemp.sepnum;
        newgslink->stemp.sepdistance = gstemp.sepdistance;
        newgslink->stemp.prop = gstemp.prop;
        newgslink->next = gshead;
        gshead = newgslink;
    } /* while still lines in data file */

    fclose(fp);

    /* Now do as Name.dat in forming species array */
    *gspno = i;
    *gspec = (struct sgenspec *) calloc(*gspno,sizeof(struct sgenspec));
    /* put each link into namelist and free it */
    i = 0;
    while (gshead) {
        (* gspec)[i].type = gshead->stemp.type;
        (* gspec)[i].targetocc = gshead->stemp.targetocc;
        (* gspec)[i].target = gshead->stemp.target;
        (* gspec)[i].target2 = gshead->stemp.target2;
        (* gspec)[i].sepnum = gshead->stemp.sepnum;
        (* gspec)[i].sepdistance = gshead->stemp.sepdistance;
        (* gspec)[i].prop = gshead->stemp.prop;
        i++;
        newgslink = gshead;
        gshead = gshead->next;
        free(newgslink);
    }

    return(i);
} // readSpeciesBlockDefinition

// read connections file: read bound.dat (boundaries) or connections.dat (connections)
// boundaries are a subset of asymmetric connections
int readConnections(int puno,struct sconnections connections[],struct spustuff pu[],
                    struct binsearch PULookup[],struct sfname fnames)
{
    FILE *fp;
    int id1,id2;
    double fcost;
    int icount = 0,idup = 0;
    int bad;
    struct sneighbour *p;
    char *readname;
    int numvars = 3,ivars;
    char *sVarName,*sVarVal;
    char sLine[500];

    readname = (char *) calloc(strlen(fnames.connectionname) + strlen(fnames.inputdir)+2, sizeof(char));
    strcpy(readname,fnames.inputdir);
    strcat(readname,fnames.connectionname);
    fp = fopen(readname,"r");
    if (fp==NULL)
    {
       displayProgress1("Warning: Connection File %s not found ",fnames.connectionname);
       free(readname);
       return(0);
    }
    free(readname);
    fgets(sLine,500-1,fp);

    while (fgets(sLine,500-1,fp))
    {
          icount++;

          sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
          sscanf(sVarVal,"%d",&id1);
          sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
          sscanf(sVarVal,"%d",&id2);
          sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
          sscanf(sVarVal,"%lf",&fcost);

          id1 = FastNameToPUID(puno,id1,PULookup);
          id2 = FastNameToPUID(puno,id2,PULookup);

          if (id1==id2)
          {  // irremovable connection
             if (asymmetricconnectivity)
                connections[id1].fixedcost = 0;
             else
                 connections[id1].fixedcost += fcost;

             continue;
          }
          if (id1>=0 && id1<puno)
          {  // Is I a sensible number ?
             p = connections[id1].first;
             bad = 0;
             while (p)
             {
                   if (p->nbr == id2)
                      bad = 1;

                   p = p->next;
             }

             if (asymmetricconnectivity)
                bad = 0;

             if (bad)
                displayProgress3("Double connection definition %i %i \n",id1,id2);
             else
             {
                 connections[id1].nbrno++;
                 p = (struct sneighbour *) malloc (sizeof(struct sneighbour));
                 p->cost = fcost;
                 p->nbr = id2;
                 p->next = connections[id1].first;
                 if (asymmetricconnectivity)
                 {
                    p->connectionorigon = 1;
                 }
                 else
                     p->connectionorigon = 1;
                 connections[id1].first = p;
             }
          }
          else
              displayErrorMessage("A connection is out of range %f %i %i \n",fcost,id1,id2);

          if (id2>=0 && id2<puno)
          {  /* Is I a sensible number ?*/
             p = connections[id2].first;
             bad = 0;
             while (p)
             {
                   if (p->nbr == id1)
                      bad = 1;

                   p = p->next;
             }

             if (asymmetricconnectivity)
                bad = 0;

             if (bad && verbosity > 4)
                displayProgress3("Double connection definition %i %i \n",id1,id2);

             if (bad)
                idup++;
             else
             {
                connections[id2].nbrno++;
                p = (struct sneighbour *) malloc (sizeof(struct sneighbour));

                #ifdef MEMDEBUG
                iMemoryUsed += sizeof(struct sneighbour);
                #endif

                p->cost = fcost;
                p->nbr = id1;
                p->next = connections[id2].first;
                p->connectionorigon = 1;

                if (asymmetricconnectivity)
                   p->connectionorigon = 0;
                else
                    p->connectionorigon = 1;

                connections[id2].first = p;
             }
          }
          else displayErrorMessage("A connection is out of range %f %i %i \n",fcost,id1,id2);
    }

    fclose(fp);

    if (idup)
       displayProgress1("There were %i duplicate connection definitions.\n",idup);

    return(icount);

} // readConnections


// read the planning unit versus species sparse matrix: ordered by planning units
void readSparseMatrix(int *iSMSize, struct spu *SM[], int puno, int spno, struct spustuff PU[],
                      struct binsearch PULookup[],struct binsearch SPLookup[],
                      struct sfname fnames)
{
     FILE *fp;
     char *readname,sLine[500],*sVarName,*sVarVal;
     int i, _spid, _puid, iInternalSMSize = 0, iBigMatrixSize, iLength;
     double amount, rDensity, rInternalSMSize, rBigMatrixSize, rProbability = 1;
     int iP, iR, iO, iB, iLastPUID;
     char cP, cR, cO, cB;

     readname = (char *) calloc(strlen(fnames.puvsprname) + strlen(fnames.inputdir)+2, sizeof(char));

     strcpy(readname,fnames.inputdir);
     strcat(readname,fnames.puvsprname);
     if((fp = fopen(readname,"r"))==NULL)
         displayErrorMessage("PU v Species file %s not found\nAborting Program.",readname);
     free(readname);

     // read through the file first to see how many lines
     fgets(sLine,500-1,fp);
     while (fgets(sLine,500-1,fp))
           iInternalSMSize++;

     rewind(fp);
     fgets(sLine,500-1,fp);

     // scan the first line to see if the prob field is tagged on the end
     // 3 = regular marxan matrix
     // 4 = prob2d marxan matrix
     cP = 'p';
     cR = 'r';
     cO = 'o';
     cB = 'b';
     iP = cP;
     iR = cR;
     iO = cO;
     iB = cB;
     iLength = strlen(sLine);
     if ((sLine[iLength-5] == iP) && (sLine[iLength-4] == iR) && (sLine[iLength-3] == iO) && (sLine[iLength-2] == iB))
        fProb2D = 1;

     *iSMSize = iInternalSMSize;

     // create the sparse matrix
     *SM = (struct spu *) calloc(iInternalSMSize,sizeof(struct spu));

     iLastPUID = -1;

     // init with zero values
     for (i=0;i<iInternalSMSize;i++)
     {
         fgets(sLine,500-1,fp);

         sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
         sscanf(sVarVal,"%d",&_spid);
         sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
         sscanf(sVarVal,"%d",&_puid);
         sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
         sscanf(sVarVal,"%lf",&amount);

         if (fProb2D == 1)
         {
            sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
            sscanf(sVarVal,"%lf",&rProbability);
         }

         if (_puid < iLastPUID)
         {
	        // error condition exists, file is not in ascending order for PUID
            AppendDebugTraceFile("Error: PU v Species file %s is not in ascending order for PUID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
	        displayErrorMessage("Error: PU v Species file %s is not in ascending order for PUID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
		 }

         iLastPUID = _puid;

         _puid = FastNameToPUID(puno,_puid,PULookup);
         _spid = FastNameToSPID(spno,_spid,SPLookup);

         /* increment richness for planning unit containing this feature */
         PU[_puid].richness += 1;
         /* if planning units richness is one, set its offset */
         if (PU[_puid].richness == 1)
            PU[_puid].offset = i;

         (* SM)[i].prob = rProbability;
         (* SM)[i].amount = amount;
         (* SM)[i].clump = 0;
         (* SM)[i].spindex = _spid;
     }

     fclose(fp);

     iBigMatrixSize = puno * spno;
     rInternalSMSize = iInternalSMSize;
     rBigMatrixSize = iBigMatrixSize;
     rDensity = rInternalSMSize / rBigMatrixSize * 100;

     displayProgress1("%i conservation values counted, %i big matrix size, %g%% density of matrix \n",
                      iInternalSMSize,iBigMatrixSize,rDensity);
} // readSparseMatrix

void readPenalties(typesp spec[],int spno,struct sfname fnames,struct binsearch SPLookup[])
{
     FILE *fp;
     char *readname,*sVarVal,sLine[500];
     int i,iSPID;
     double rPenalty;

     readname = (char *) calloc(strlen(fnames.penaltyname) + strlen(fnames.inputdir)+2, sizeof(char));

     strcpy(readname,fnames.inputdir);
     strcat(readname,fnames.penaltyname);
     if((fp = fopen(readname,"r"))==NULL)
         displayErrorMessage("Penalty file %s not found\nAborting Program.",readname);
     free(readname);

     for (i=0;i<spno;i++)
         spec[i].rUserPenalty = 0;

     fgets(sLine,500-1,fp);
     while (fgets(sLine,500-1,fp))
     {
           sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
           sscanf(sVarVal,"%d",&iSPID);
           sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
           sscanf(sVarVal,"%lf",&rPenalty);

           i = FastNameToSPID(spno,iSPID,SPLookup);

           spec[i].rUserPenalty = rPenalty;

           AppendDebugTraceFile("readPenalties spname %i user penalty %g\n",spec[i].name,rPenalty);
	 }

     fclose(fp);
}

// read the planning unit versus species sparse matrix: ordered by species
void readSparseMatrixSpOrder(int *iSMSize, struct spusporder *SM[], int puno, int spno,
                             struct binsearch PULookup[],struct binsearch SPLookup[],
                             struct sfname fnames)
{
     FILE *fp;
     char *readname,sLine[500],*sVarName,*sVarVal;
     int i, _spid,spid, _puid, iInternalSMSize = 0, iBigMatrixSize, iLastSPID;
     double amount, rDensity, rInternalSMSize, rBigMatrixSize;

     readname = (char *) calloc(strlen(fnames.matrixspordername) + strlen(fnames.inputdir)+2, sizeof(char));

     strcpy(readname,fnames.inputdir);
     strcat(readname,fnames.matrixspordername);
     if((fp = fopen(readname,"r"))==NULL)
         displayErrorMessage("PU v Species file %s not found\nAborting Program.",readname);
     free(readname);

     // read through the file first to see how many lines
     fgets(sLine,500-1,fp);
     while (fgets(sLine,500-1,fp))
        iInternalSMSize++;

     rewind(fp);
     fgets(sLine,500-1,fp);

     *iSMSize = iInternalSMSize;

     // create the sparse matrix
     *SM = (struct spusporder *) calloc(iInternalSMSize,sizeof(struct spusporder));

     iLastSPID = -1;
     // planning unit richness and offset are already set to zero

     // init with zero values
     for (i=0;i<iInternalSMSize;i++)
     {

         fgets(sLine,500-1,fp);

         sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
         sscanf(sVarVal,"%d",&_spid);
         sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
         sscanf(sVarVal,"%d",&_puid);
         sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
         sscanf(sVarVal,"%lf",&amount);

         if (_spid < iLastSPID)
         {
	        // error condition exists, file is not in ascending order for SPID
            AppendDebugTraceFile("Error: PU v Species file %s is not in ascending order for SPID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
	        displayErrorMessage("Error: PU v Species file %s is not in ascending order for SPID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
		 }

         iLastSPID = _spid;

         _puid = FastNameToPUID(puno,_puid,PULookup);
         spid = FastNameToSPID(spno,_spid,SPLookup);

         // increment richness for planning unit containing this feature
         spec[spid].richness += 1;
         // if planning units richness is one, set its offset
         if (spec[spid].richness == 1)
            spec[spid].offset = i;

         (* SM)[i].amount = amount;
         (* SM)[i].puindex = _puid;
     }

     fclose(fp);

     iBigMatrixSize = puno * spno;
     rInternalSMSize = iInternalSMSize;
     rBigMatrixSize = iBigMatrixSize;
     rDensity = rInternalSMSize / rBigMatrixSize * 100;

     displayProgress1("%i conservation values counted, %i big matrix size, %g%% density of matrix \n",
                      iInternalSMSize,iBigMatrixSize,rDensity);
} // readSparseMatrixSpOrder

