// functions that read from input files

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
    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading planning units.\n");

    sVarName = strtok(sLine," ,;:^*\"/\t\'\\\n");
    head = storeFieldName(varlist,numvars,sVarName,head,fnames.puname);

    ivars = 1;
    temp = head;
    while ((sVarName = strtok(NULL," ,;:^*\"/\t\'\\\n")) != NULL)
    {
          ivars++;
          temp->next = storeFieldName(varlist,numvars,sVarName,head,fnames.puname);
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

    return(*puno);
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
    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading species.\n");

    sVarName = strtok(sLine," ,;:^*\"/\t\'\\\n");
    snhead = storeFieldName(varlist,numvars,sVarName,snhead,fnames.specname);
    ivars = 1;
    temp = snhead;
    while ((sVarName = strtok(NULL," ,;:^*\"/\t\'\\\n")) != NULL) {
        ivars++;
        temp->next = storeFieldName(varlist,numvars,sVarName,snhead,fnames.specname);
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
    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading block definition.\n");

    sVarName = strtok(sLine," ,;:^*\"/\t\'\\\n");
    head = storeFieldName(varlist,numvars,sVarName,head,fnames.blockdefname);
    ivars = 1;
    temp = head;
    while ((sVarName = strtok(NULL," ,;:^*\"/\t\'\\\n")) != NULL) {
        ivars++;
        temp->next = storeFieldName(varlist,numvars,sVarName,head,fnames.blockdefname);
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

    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading connections.\n");

    while (fgets(sLine,500-1,fp))
    {
          icount++;

          sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
          sscanf(sVarVal,"%d",&id1);
          sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
          sscanf(sVarVal,"%d",&id2);
          sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
          sscanf(sVarVal,"%lf",&fcost);

          id1 = binarySearchPuIndex(puno,id1,PULookup);
          id2 = binarySearchPuIndex(puno,id2,PULookup);

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
    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading sparse matrix.\n");

    while (fgets(sLine,500-1,fp))
        iInternalSMSize++;

    rewind(fp);

    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading sparse matrix.\n");

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
        if (fgets(sLine,500-1,fp) == NULL)
            displayErrorMessage("Error reading sparse matrix.\n");

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
            appendTraceFile("Error: PU v Species file %s is not in ascending order for PUID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
            displayErrorMessage("Error: PU v Species file %s is not in ascending order for PUID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
        }

        iLastPUID = _puid;

        _puid = binarySearchPuIndex(puno,_puid,PULookup);
        _spid = binarySearchSpecIndex(spno,_spid,SPLookup);

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

    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading penalties.\n");

    while (fgets(sLine,500-1,fp))
    {
        sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
        sscanf(sVarVal,"%d",&iSPID);
        sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
        sscanf(sVarVal,"%lf",&rPenalty);

        i = binarySearchSpecIndex(spno,iSPID,SPLookup);

        spec[i].rUserPenalty = rPenalty;

        appendTraceFile("readPenalties spname %i user penalty %g\n",spec[i].name,rPenalty);
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
    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading sparse matrix sporder.\n");

    while (fgets(sLine,500-1,fp))
        iInternalSMSize++;

    rewind(fp);
    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading sparse matrix sporder.\n");

    *iSMSize = iInternalSMSize;

    // create the sparse matrix
    *SM = (struct spusporder *) calloc(iInternalSMSize,sizeof(struct spusporder));

    iLastSPID = -1;
    // planning unit richness and offset are already set to zero

    // init with zero values
    for (i=0;i<iInternalSMSize;i++)
    {
        if (fgets(sLine,500-1,fp) == NULL)
            displayErrorMessage("Error reading sparse matrix sporder.\n");

        sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
        sscanf(sVarVal,"%d",&_spid);
        sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
        sscanf(sVarVal,"%d",&_puid);
        sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
        sscanf(sVarVal,"%lf",&amount);

        if (_spid < iLastSPID)
        {
            // error condition exists, file is not in ascending order for SPID
            appendTraceFile("Error: PU v Species file %s is not in ascending order for SPID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
            displayErrorMessage("Error: PU v Species file %s is not in ascending order for SPID at record %i.\nAborting Program.",fnames.puvsprname,i+1);
        }

        iLastSPID = _spid;

        _puid = binarySearchPuIndex(puno,_puid,PULookup);
        spid = binarySearchSpecIndex(spno,_spid,SPLookup);

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

// read value for a single parameter specified in input.dat parameter file
void readInputOption(FILE *infile, char varname[], void *address, int parmtype, int crit,int present)
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

    rewind(infile); // search from top of infile

    do
    {   /* loop through file looking for varname */
        if (fgets(buffer,255,infile) == NULL)
            //displayErrorMessage("Error reading input parameter %s.\n",varname)
            ;
        check1 = 0;
        check2 = 0;

        while (buffer[check1++] == varname[check2++])
            ;

        if (check1 > namelen)
        {  // varname matches upto namelen
            foundit++;
            switch (parmtype)
            {
                case REAL :
                    gotit = sscanf(&buffer[check1]," %f", (float *) address);
                    break;
                case DOUBLE :
                    gotit = sscanf(&buffer[check1]," %lf", (double *) address);
                    break;
                case INTEGER :
                    gotit = sscanf(&buffer[check1]," %d", (int *) address);
                    break;
                case LONGINT :
                    gotit = sscanf(&buffer[check1]," %ld", (long int *) address);
                    break;
                case STRING :
                    // trim leading and trailing blanks (allow spaces which are important for directory names
                    check1 += strspn(&buffer[check1]," ,");
                    for (check2 = strlen(&buffer[check1])-1;isspace(buffer[check1+check2]) != 0;check2--)
                        ; // Find last non space character
                    if (strlen(&buffer[check1]) <2)
                        buffer[check1] = '\0';
                    buffer[check1 + check2+1] = '\0';

                    strcpy((char *) address,&buffer[check1]);

                    gotit = 1; // So that var check works. This needs further consideration
                    break;
                default :
                    displayErrorMessage("Invalid parameter type request %d: \n",parmtype);
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
} // readInputOption

// read values for the parameters specified in input.dat parameter file
void readInputOptions(double *cm,double *prop,struct sanneal *anneal,
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

    readInputOption(fp,"VERSION",&version,DOUBLE,0,present);
    readInputOption(fp,"BLM",cm,DOUBLE,0,present);
    if (present == 0)
        readInputOption(fp,"CM",cm,DOUBLE,0,present);
    readInputOption(fp,"PROP",prop,DOUBLE,0,present);
    readInputOption(fp,"RANDSEED",iseed,INTEGER,0,present);

    /* Annealing Controls */
    readInputOption(fp,"NUMITNS",&(*anneal).iterations,LONGINT,0,present);
    readInputOption(fp,"STARTTEMP",&(*anneal).Tinit,DOUBLE,0,present);
    readInputOption(fp,"COOLFAC",&(*anneal).Tcool,DOUBLE,0,present);
    readInputOption(fp,"NUMTEMP",&(*anneal).Titns,INTEGER,0,present);

    (*anneal).type = 1;
    if ((*anneal).iterations < 1 )
        (*anneal).type = 0;
    if ((*anneal).Tinit < 0)
        (*anneal).type = (int) (-(*anneal).Tinit) + 1;  /* type is negative of Tinit */
    if (fscanf(fp,"%i",iseed) == 0)  /* The random seed. -1 to set by clock */
        displayErrorMessage("Error reading input parameter random seed.\n");

    /* Various controls */
    readInputOption(fp,"NUMREPS",repeats,LONGINT,0,present);
    readInputOption(fp,"COSTTHRESH",costthresh,DOUBLE,0,present);
    readInputOption(fp,"THRESHPEN1",tpf1,DOUBLE,0,present);
    readInputOption(fp,"THRESHPEN2",tpf2,DOUBLE,0,present);

    /* SaveFiles */
    readInputOption(fp,"SCENNAME",savename,STRING,0,present);

    /* SaveFiles New Method */
    readInputOption(fp,"SAVERUN",&(*fnames).saverun,INTEGER,0,present);
    readInputOption(fp,"SAVEBEST",&(*fnames).savebest,INTEGER,0,present);
    readInputOption(fp,"SAVESUMMARY",&(*fnames).savesum,INTEGER,0,present);
    readInputOption(fp,"SAVESCEN",&(*fnames).savesen,INTEGER,0,present);
    readInputOption(fp,"SAVETARGMET",&(*fnames).savespecies,INTEGER,0,present);
    readInputOption(fp,"SAVESUMSOLN",&(*fnames).savesumsoln,INTEGER,0,present);
    readInputOption(fp,"SAVESPEC",&(*fnames).savespec,INTEGER,0,present);
    readInputOption(fp,"SAVEPU",&(*fnames).savepu,INTEGER,0,present);
    readInputOption(fp,"SAVEMATRIXPUORDER",&(*fnames).savepuvspr,INTEGER,0,present);
    readInputOption(fp,"SAVEMATRIXSPORDER",&(*fnames).savematrixsporder,INTEGER,0,present);
    readInputOption(fp,"SAVEPENALTY",&(*fnames).savepenalty,INTEGER,0,present);
    readInputOption(fp,"SAVETOTALAREAS",&(*fnames).savetotalareas,INTEGER,0,present);
    readInputOption(fp,"SAVEDEBUGTRACEFILE",&(*fnames).savedebugtracefile,INTEGER,0,present);
    readInputOption(fp,"SAVERICHNESS",&(*fnames).saverichness,INTEGER,0,present);
    readInputOption(fp,"SAVESOLUTIONSMATRIX",&(*fnames).savesolutionsmatrix,INTEGER,0,present);
    readInputOption(fp,"SOLUTIONSMATRIXHEADERS",&(*fnames).solutionsmatrixheaders,INTEGER,0,present);
    readInputOption(fp,"SAVELOG",&(*fnames).savelog,INTEGER,0,present);
    readInputOption(fp,"SAVESNAPSTEPS",&(*fnames).savesnapsteps,INTEGER,0,present);
    readInputOption(fp,"SAVESNAPCHANGES",&(*fnames).savesnapchanges,INTEGER,0,present);
    readInputOption(fp,"SAVESNAPFREQUENCY",&(*fnames).savesnapfrequency,INTEGER,0,present);
    readInputOption(fp,"SAVEANNEALINGTRACE",&(*fnames).saveannealingtrace,INTEGER,0,present);
    readInputOption(fp,"ANNEALINGTRACEROWS",&(*fnames).annealingtracerows,INTEGER,0,present);
    readInputOption(fp,"SAVEITIMPTRACE",&(*fnames).saveitimptrace,INTEGER,0,present);
    readInputOption(fp,"ITIMPTRACEROWS",&(*fnames).itimptracerows,INTEGER,0,present);
    readInputOption(fp,"RIMAGETYPE",&(*fnames).rimagetype,INTEGER,0,present);
    readInputOption(fp,"REXECUTESCRIPT",&(*fnames).rexecutescript,INTEGER,0,present);
    readInputOption(fp,"RCLUSTERCOUNT",&(*fnames).rclustercount,INTEGER,0,present);
    readInputOption(fp,"RIMAGEWIDTH",&(*fnames).rimagewidth,INTEGER,0,present);
    readInputOption(fp,"RIMAGEHEIGHT",&(*fnames).rimageheight,INTEGER,0,present);
    readInputOption(fp,"RIMAGEFONTSIZE",&(*fnames).rimagefontsize,INTEGER,0,present);
    readInputOption(fp,"ASYMMETRICCONNECTIVITY",&asymmetricconnectivity,INTEGER,0,present);
    readInputOption(fp,"CONNECTIVITYIN",&fOptimiseConnectivityIn,INTEGER,0,present);

    // quantum annealing control parameters
    readInputOption(fp,"QAPROP",&rQAPROP,DOUBLE,0,present);
    readInputOption(fp,"QADECAY",&rQADECAY,DOUBLE,0,present);
    readInputOption(fp,"QADECAYB",&rQADECAYB,DOUBLE,0,present);
    readInputOption(fp,"QADECAYTYPE",&iQADECAYTYPE,INTEGER,0,present);
    readInputOption(fp,"QAACCPR",&rQAACCPR,DOUBLE,0,present);

    if (!(*fnames).savesnapfrequency)
        (*fnames).savesnapfrequency = 1;

    /* Filenames */
    readInputOption(fp,"INPUTDIR",stemp,STRING,1,present);
    if (stemp[strlen(stemp)-1] != '/' && stemp[strlen(stemp)-1] != '\\')
        strcat(stemp,"/");
    (*fnames).inputdir = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).inputdir,stemp);

    readInputOption(fp,"OUTPUTDIR",stemp,STRING,1,present);
    if (stemp[strlen(stemp)-1] != '/' && stemp[strlen(stemp)-1] != '\\')
        strcat(stemp,"/");
    (*fnames).outputdir = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).outputdir,stemp);

    strcpy(stemp,"PU.dat");
    readInputOption(fp,"PUNAME",stemp,STRING,1,present);
    (*fnames).puname = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).puname,stemp);

    strcpy(stemp,"spec.dat");
    readInputOption(fp,"SPECNAME",stemp,STRING,1,present);
    (*fnames).specname = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).specname,stemp);

    strcpy(stemp,"puvspr2.dat");
    readInputOption(fp,"PUVSPRNAME",stemp,STRING,1,present);
    (*fnames).puvsprname = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).puvsprname,stemp);

    strcpy(stemp,"NULL");
    readInputOption(fp,"MATRIXSPORDERNAME",stemp,STRING,0,present);
    (*fnames).matrixspordername = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).matrixspordername,stemp);

    strcpy(stemp,"NULL");
    readInputOption(fp,"PENALTYNAME",stemp,STRING,0,present);
    (*fnames).penaltyname = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).penaltyname,stemp);

    strcpy(stemp,"NULL");
    readInputOption(fp,"BOUNDNAME",stemp,STRING,0,present);
    if (present == 0)
        readInputOption(fp,"CONNECTIONNAME",stemp,STRING,0,present);
    (*fnames).connectionname = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).connectionname,stemp);

    strcpy(stemp,"NULL");
    readInputOption(fp,"CONNECTIONFILESNAME",stemp,STRING,0,present);
    (*fnames).connectionfilesname = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).connectionfilesname,stemp);

    strcpy(stemp,"NULL");
    readInputOption(fp,"BLOCKDEFNAME",stemp,STRING,0,present);
    (*fnames).blockdefname = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).blockdefname,stemp);

    strcpy(stemp,"SOLUTION");
    readInputOption(fp,"BESTFIELDNAME",stemp,STRING,0,present);
    (*fnames).bestfieldname = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).bestfieldname,stemp);


    strcpy(stemp,"NULL");
    readInputOption(fp,"RBINARYPATHNAME",stemp,STRING,0,present);
    (*fnames).rbinarypath = (char *) calloc(strlen(stemp)+1,sizeof(char));
    strcpy((*fnames).rbinarypath,stemp);

    /* various other controls */
    readInputOption(fp,"RUNMODE",runopts,INTEGER,1,present);
    readInputOption(fp,"MISSLEVEL",misslevel,DOUBLE,0,present);
    readInputOption(fp,"HEURTYPE",heurotype,INTEGER,0,present);
    readInputOption(fp,"CLUMPTYPE",clumptype,INTEGER,0,present);
    readInputOption(fp,"ITIMPTYPE",itimptype,INTEGER,0,present);
    readInputOption(fp,"VERBOSITY",verb,INTEGER,0,present);
    verbosity = *verb;
    readInputOption(fp,"PROBABILITYWEIGHTING",stemp,STRING,0,present);
    sscanf(stemp, "%lf", &rProbabilityWeighting);

    readInputOption(fp,"STARTDECTHRESH",&rStartDecThresh,DOUBLE,0,present);
    readInputOption(fp,"ENDDECTHRESH",&rEndDecThresh,DOUBLE,0,present);
    readInputOption(fp,"STARTDECMULT",&rStartDecMult,DOUBLE,0,present);
    readInputOption(fp,"ENDDECMULT",&rEndDecMult,DOUBLE,0,present);

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
} // readInputOptions
