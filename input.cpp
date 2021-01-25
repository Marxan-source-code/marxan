
#include <algorithm>
#include <cstring>
#include <ctime>
#include <fstream>
#include <locale>
#include <map>
#include <sstream>

#include "marxan.hpp"
#include "utils.hpp"

// functions that read from input files
namespace marxan {

string storeFieldName(vector<string>& varlist,int numvars,string sVarName,
    vector<string> &head,string fname) {

    // add a field name from an input file header to a list
    int foundit = 0;
    int temp, newlink;

    for (int i=0;(i<numvars && foundit==0);i++)
    {
        if (varlist[i].compare(sVarName) == 0)
            foundit++;
    }

    if (!head.empty() && find(head.begin(), head.end(), sVarName) != head.end())
    {
        displayErrorMessage("ERROR: variable %s has been defined twice in data file %s.\n",sVarName.c_str(),fname.c_str());
    }

    return(sVarName);
}

vector<string> GetFieldNames(string readname, string fname, FILE*& fp, int& ivars, vector<string>& varList) {
    char sLine[500];
    vector<string> fieldNames;
    int numvars = varList.size();

    /* Find and Open File */
    if((fp =  fopen(readname.c_str(),"r"))==NULL) {
        displayWarningMessage("Warning: File %s not found.\n", fname.c_str());
        return fieldNames; // return empty list
    }

    /* Scan header */
    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading %s.\n", fname.c_str());

    char* temp = strtok(sLine," ,;:^*\"/\t\'\\\n");
    string sVarName(temp);
    fieldNames.push_back(storeFieldName(varList,numvars,sVarName,fieldNames,fname));
    ivars = 1;

    while ((temp = strtok(NULL," ,;:^*\"/\t\'\\\n")) != NULL) {
        sVarName = temp;
        ivars++;
        fieldNames.push_back(storeFieldName(varList,numvars,sVarName,fieldNames,fname));
    }  /* tokking out all the variable names from the header line. There are numVars of them*/

    return fieldNames;
}

// read the planning units file pu.dat
int readPlanningUnits(int& puno, vector<spustuff>& pu, sfname& fnames)
{
    FILE *fp;
    string readname;
    char sLine[600];
    vector<string> varlist = {"id","cost","status","xloc","yloc","prob"};
    int ivars,i = 0;
    char* sVarVal;

    readname = fnames.inputdir + fnames.puname;

    #ifdef MEMDEBUG
    iMemoryUsed += (strlen(fnames.puname) + strlen(fnames.inputdir)+2) * sizeof(char);
    displayProgress1("memory used %i\n",iMemoryUsed);
    #endif

    vector<string> head = GetFieldNames(readname, fnames.puname, fp, ivars, varlist);
    if(head.empty())
        displayErrorMessage("Planning Unit file %s has not been found.\nAborting Program.",readname.c_str());

    /* While there are still lines left feed information into temporary list */
    while (fgets(sLine,500-1,fp))
    {
        i++;
        spustuff putemp;
        putemp.id = -1; /* Set defaults for any missing values */
        putemp.cost = 1;
        putemp.status = 0;
        putemp.xloc = -1;
        putemp.yloc = -1;
        putemp.prob = 0;
        putemp.richness = 0;
        putemp.offset = 0;
        putemp.probrichness = 0;
        putemp.proboffset = 0;

        for (string temp : head)
        {
            if (temp == head.front())
                sVarVal = strtok(sLine, " ,;:^*\"/\t\'\\\n");
            else
                sVarVal = strtok(NULL, " ,;:^*\"/\t\'\\\n");

            if (temp.compare("id") == 0)
            {
                sscanf(sVarVal, "%d", &putemp.id);
            }
            else if (temp.compare("status") == 0)
            {
                sscanf(sVarVal, "%d", &putemp.status);
            }
            else if (temp.compare("xloc") == 0)
            {
                sscanf(sVarVal, "%lf", &putemp.xloc);
            }
            else if (temp.compare("yloc") == 0)
            {
                sscanf(sVarVal, "%lf", &putemp.yloc);
            }
            else if (temp.compare("cost") == 0)
            {
                sscanf(sVarVal, "%lf", &putemp.cost);
            }
            else if (temp.compare("prob") == 0)
            {
                sscanf(sVarVal, "%lf", &putemp.prob);
                iProbFieldPresent = 1;
            }
        } /* looking for ivar different input variables */

        if (putemp.id == -1)
            displayErrorMessage("ERROR: Missing planning unit id for line %d. \n",i);

        pu.push_back(putemp);

    } /* while still lines in data file */

    fclose(fp);

    /* Create array to store the information */
    puno = i;

    #ifdef MEMDEBUG
    iMemoryUsed += (*puno) * sizeof(struct spustuff);
    displayProgress1("memory used %i\n",iMemoryUsed);
    #endif

    if (iProbFieldPresent == 1)
    {
        fProb1D = 1;
    }

    return puno;
} // readPlanningUnits

// read species file: spec.dat
int readSpecies(int &spno, vector<sspecies>& spec, sfname& fnames)
{
    FILE *fp;
    int n=0;
    string readname;
    char sLine[500];
    vector<string> varlist = {"id","type","target","spf",
                         "target2","sepdistance","sepnum","name",
                         "targetocc","prop","ptarget1d","ptarget2d"};
    int ivars;
    char* sVarVal;

    readname = fnames.inputdir + fnames.specname;
    vector<string> snhead = GetFieldNames(readname, fnames.specname, fp, ivars, varlist);
    if(snhead.empty())
        displayErrorMessage("Species file %s has not been found.\nAborting Program.",readname.c_str());

    // While there are still lines left feed information into temporary link list
    while (fgets(sLine,500-1,fp))
    {
        n++;
        // Clear important species stats
        sspecies spectemp;
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
        spectemp.richness = 0;
        spectemp.probability1D = 0;
        spectemp.probability2D = 0;
        spectemp.Zscore1D = 0;
        spectemp.Zscore2D = 0;

        for (string temp : snhead)
        {
            if (temp == snhead.front())
                sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
            else
                sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");

            if (temp.compare("id")==0)
            {
                sscanf(sVarVal,"%d",&spectemp.name);
            }
            else if (temp.compare("type")==0)
            {
                sscanf(sVarVal,"%d",&spectemp.type);
            }
            else if (temp.compare("target")==0)
            {
                sscanf(sVarVal,"%lf",&spectemp.target);
            }
            else if (temp.compare("prop") == 0)
            {
                sscanf(sVarVal, "%lf", &spectemp.prop);
                if (spectemp.prop > 0)
                    fSpecPROPLoaded = 1;
            }
            else if (temp.compare("spf") == 0)
            {
                sscanf(sVarVal, "%lf", &spectemp.spf);
            }
            else if (temp.compare("sepdistance") == 0)
            {
                sscanf(sVarVal, "%lf", &spectemp.sepdistance);
            }
            else if (temp.compare("sepnum") == 0)
            {
                sscanf(sVarVal, "%d", &spectemp.sepnum);
            }
            else if (temp.compare("target2") == 0)
            {
                sscanf(sVarVal, "%lf", &spectemp.target2);
            }
            else if (temp.compare("targetocc") == 0)
            {
                sscanf(sVarVal, "%d", &spectemp.targetocc);
            }
            else if (temp.compare("ptarget1d") == 0)
            {
                sscanf(sVarVal, "%lf", &spectemp.ptarget1d);
            }
            else if (temp.compare("ptarget2d") == 0)
            {
                sscanf(sVarVal, "%lf", &spectemp.ptarget2d);
            }
        } // looking for ivar different input variables
        
        spec.push_back(spectemp);
    } // Scanning through each line of file

    fclose(fp);
    spno = n;
    return(n);
}  // readSpecies

// read species block definition file
// allows attribute to be applied to species based upon a type
int readSpeciesBlockDefinition(int& gspno, vector<sgenspec>& gspec, sfname& fnames)
{
    FILE *fp;
    string readname;
    char sLine[500];
    vector<string> varlist = {"type","target","target2","targetocc",
                        "sepnum","sepdistance","prop","spf"};
    int ivars,i=0;
    char* sVarVal;

    /* Find and Open File */
    readname = fnames.inputdir + fnames.blockdefname;
    vector<string> head = GetFieldNames(readname, fnames.blockdefname, fp, ivars, varlist);

    /* While there are still lines left feed information into temporary link list */
    while (fgets(sLine,500-1,fp)) {
        i++;
        sgenspec gstemp;
        gstemp.type = -1; /* Set defaults for any missing values */
        gstemp.targetocc = -1;
        gstemp.target = -1;
        gstemp.target2 = -1;
        gstemp.sepnum = -1;
        gstemp.sepdistance = -1;
        gstemp.prop = -1;

        for (string temp : head)
        {
            if (temp == head.front())
                sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
            else
                sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");

            if (temp.compare("type")==0) {
                sscanf(sVarVal,"%d",&gstemp.type);
            }
            else if (temp.compare("targetocc")==0){
                sscanf(sVarVal,"%d",&gstemp.targetocc);
            }
            else if (temp.compare("target")==0){
                sscanf(sVarVal,"%lf",&gstemp.target);
            }
            else if (temp.compare("target2")==0){
                sscanf(sVarVal,"%lf",&gstemp.target2);
            }
            else if (temp.compare("sepnum")==0){
                sscanf(sVarVal,"%d",&gstemp.sepnum);
            }
            else if (temp.compare("sepdistance")==0){
                sscanf(sVarVal,"%lf",&gstemp.sepdistance);
            }
            else if (temp.compare("prop")==0){
                sscanf(sVarVal,"%lf",&gstemp.prop);
            }
            else if (temp.compare("spf")==0){
                sscanf(sVarVal,"%lf",&gstemp.spf);
            }
            else 
            {
                displayWarningMessage("Cannot find >%s< \n",temp.c_str());
                displayErrorMessage("Serious error in GenSpecies data reading function.\n");
            }
        } /* looking for ivar different input variables */

        if (gstemp.type== -1)
            displayErrorMessage("ERROR: Missing Gen Species type for line %d. \n",i);

        gspec.push_back(gstemp);

    } /* while still lines in data file */

    fclose(fp);

    gspno = i;
    return(i);
} // readSpeciesBlockDefinition

// read connections file: read bound.dat (boundaries) or connections.dat (connections)
// boundaries are a subset of asymmetric connections
int readConnections(int& puno, vector<sconnections>& connections, vector<spustuff>& pu,
                    map<int,int>& PULookup, sfname& fnames)
{
    FILE *fp;
    int id1,id2;
    double fcost;
    int icount = 0,idup = 0, bad;
    string readname;
    int ivars;
    char *sVarVal;
    char sLine[500];

    readname = fnames.inputdir + fnames.connectionname;

    fp = fopen(readname.c_str(),"r");
    if (fp==NULL)
    {
       displayProgress1("Warning: Connection File %s not found ", fnames.connectionname.c_str());
       return(0);
    }

    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading connections.\n");

    while (fgets(sLine, 500 - 1, fp))
    {
        icount++;

        sVarVal = strtok(sLine, " ,;:^*\"/\t\'\\\n");
        sscanf(sVarVal, "%d", &id1);
        sVarVal = strtok(NULL, " ,;:^*\"/\t\'\\\n");
        sscanf(sVarVal, "%d", &id2);
        sVarVal = strtok(NULL, " ,;:^*\"/\t\'\\\n");
        sscanf(sVarVal, "%lf", &fcost);

        id1 = PULookup[id1];
        id2 = PULookup[id2];

        if (id1 == id2)
        { // irremovable connection
            if (asymmetricconnectivity)
                connections[id1].fixedcost = 0;
            else
                connections[id1].fixedcost += fcost;

            continue;
        }

        if (id1 >= 0 && id1 < puno)
        { // Is I a sensible number ?
            bad = 0;
            for (sneighbour p: connections[id1].first)
            {
                if (p.nbr == id2)
                    bad = 1;
            }

            if (asymmetricconnectivity)
                bad = 0;

            if (bad)
                displayProgress3("Double connection definition %i %i \n", id1, id2);
            else
            {
                connections[id1].nbrno++;
                sneighbour p(id2, fcost, 1);

                connections[id1].first.push_back(p);
            }
        }
        else
            displayErrorMessage("A connection is out of range %f %i %i \n", fcost, id1, id2);

        if (id2 >= 0 && id2 < puno)
        { /* Is I a sensible number ?*/
            bad = 0;
            for (sneighbour p: connections[id2].first)
            {
                if (p.nbr == id1)
                    bad = 1;
            }

            if (asymmetricconnectivity)
                bad = 0;

            if (bad && verbosity > 4)
                displayProgress3("Double connection definition %i %i \n", id1, id2);

            if (bad)
                idup++;
            else
            {
                connections[id2].nbrno++;
                sneighbour p(id1, fcost, 1);

#ifdef MEMDEBUG
                iMemoryUsed += sizeof(struct sneighbour);
#endif

                if (asymmetricconnectivity)
                    p.connectionorigon = 0;

                connections[id2].first.push_back(p);
            }
        }
        else
            displayErrorMessage("A connection is out of range %f %i %i \n", fcost, id1, id2);
    }

    fclose(fp);

    if (idup)
       displayProgress1("There were %i duplicate connection definitions.\n",idup);

    return(icount);

} // readConnections

// functions for Matt's Big O notation optimisation
void readSparseMatrix(int &iSMSize, vector<spu> &SM, int puno, int spno, vector<spustuff> &pu,
                      map<int,int> &PULookup,map<int,int> &SPLookup,
                      sfname& fnames) {
    vector<map<int,spu>> SMTemp; // temporarily storing in this structure prevents the need for ordering.
    FILE *fp;
    string readname;
    char sLine[500],*sVarVal;
    int i, _spid, _puid, iInternalSMSize = 0, iBigMatrixSize, iLength;
    double amount, rDensity, rInternalSMSize, rBigMatrixSize, rProbability = 1;
    int iLastPUID;

    readname = fnames.inputdir + fnames.puvsprname;
    if((fp = fopen(readname.c_str(),"r"))==NULL)
        displayErrorMessage("PU v Species file %s not found\nAborting Program.",readname.c_str());

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
    iLength = strlen(sLine);
    if ((sLine[iLength-5] == 'p') && (sLine[iLength-4] == 'r') && (sLine[iLength-3] == 'o') && (sLine[iLength-2] == 'b'))
        fProb2D = 1;

    iSMSize = iInternalSMSize;

    // create the sparse matrix
    SMTemp.resize(puno);
    SM.resize(iInternalSMSize);

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

        int old_puid;
        try {
            old_puid = _puid;
            _puid = PULookup.at(_puid);
        }
        catch (out_of_range ex) {
            displayWarningMessage("Puid %d found in puvspr file but not found in pu file. Ignoring.\n", old_puid);
        }

        int old_spid;
        try {
            old_spid = _spid;
            _spid = SPLookup.at(_spid);
        }
        catch (out_of_range ex) {
            displayWarningMessage("Spid %d found in puvspr file but not found in spec file. Ignoring.\n", old_spid);
        }

        /* increment richness for planning unit containing this feature */
        pu[_puid].richness += 1;
        SMTemp[_puid][_spid].prob = rProbability;
        SMTemp[_puid][_spid].amount = amount;
    }

    fclose(fp);

    int j = 0;
    // Populate real SM using SMTemp
    for (int i = 0; i < puno; i++) {
        pu[i].offset = j;
        for (auto &[spindex, val] : SMTemp[i]) {
            SM[j].amount = val.amount;
            SM[j].prob = val.prob;
            SM[j].spindex = spindex;
            j++;
        }
    }

    iBigMatrixSize = puno * spno;
    rInternalSMSize = iInternalSMSize;
    rBigMatrixSize = iBigMatrixSize;
    rDensity = rInternalSMSize / rBigMatrixSize * 100;

    displayProgress1("%i conservation values counted, %i big matrix size, %g%% density of matrix \n",
                     iInternalSMSize,iBigMatrixSize,rDensity);
} // readSparseMatrix

void readPenalties(vector<sspecies> &spec,int spno,sfname& fnames,map<int,int> &SPLookup) {
    FILE *fp;
    string readname;
    char *sVarVal,sLine[500];
    int i,iSPID;
    double rPenalty;

    readname = fnames.inputdir + fnames.penaltyname;
    if((fp = fopen(readname.c_str(),"r"))==NULL)
        displayErrorMessage("Penalty file %s not found\nAborting Program.",readname.c_str());

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

        i = SPLookup[iSPID];
        spec[i].rUserPenalty = rPenalty;

        appendTraceFile("readPenalties spname %i user penalty %g\n",spec[i].name,rPenalty);
    }

    fclose(fp);
}

// read the planning unit versus species sparse matrix: ordered by species
void readSparseMatrixSpOrder(int &iSMSize, vector<spusporder> &SM, int puno, int spno,
                             map<int,int> &PULookup,map<int,int> &SPLookup, vector<sspecies> &spec,
                             sfname &fnames) {
    vector<map<int,spusporder>> SMTemp;
    FILE *fp;
    string readname;
    char sLine[500],*sVarVal;
    int _spid, _puid, iInternalSMSize = 0, iBigMatrixSize, iLastSPID;
    double amount, rDensity, rInternalSMSize, rBigMatrixSize;

    readname = fnames.inputdir + fnames.matrixspordername;
    if((fp = fopen(readname.c_str(),"r"))==NULL)
        displayErrorMessage("PU v Species file %s not found\nAborting Program.",readname.c_str());

    // read through the file first to see how many lines
    if (fgets(sLine,500-1,fp) == NULL)
        displayErrorMessage("Error reading sparse matrix sporder.\n");

    // create the sparse matrix
    SMTemp.resize(spno);

    // planning unit richness and offset are already set to zero
    // init with zero values
    while (fgets(sLine,500-1,fp))
    {
        sVarVal = strtok(sLine," ,;:^*\"/\t\'\\\n");
        sscanf(sVarVal,"%d",&_spid);
        sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
        sscanf(sVarVal,"%d",&_puid);
        sVarVal = strtok(NULL," ,;:^*\"/\t\'\\\n");
        sscanf(sVarVal,"%lf",&amount);

        int old_puid;
        try {
            old_puid = _puid;
            _puid = PULookup.at(_puid);
        }
        catch (out_of_range ex) {
            displayWarningMessage("Puid %d found in puvspr file but not found in pu file. Ignoring.\n", old_puid);
        }

        int old_spid;
        try {
            old_spid = _spid;
            _spid = SPLookup.at(_spid);
        }
        catch (out_of_range ex) {
            displayWarningMessage("Spid %d found in puvspr file but not found in spec file. Ignoring.\n", old_spid);
        }

        // increment richness for planning unit containing this feature
        spec[_spid].richness += 1;
        SMTemp[_spid][_puid].amount = amount;
    }

    fclose(fp);

    // Fill the SM vector
    int j = 0;
    for (int i = 0; i < spno; i++) {
        spec[i].offset = j;
        for (auto& [puindex, value]: SMTemp[i]) {
            spusporder temp;
            temp.amount = value.amount;
            temp.puindex = puindex;
            SM.push_back(temp);
            j++;
        }
    }

    iBigMatrixSize = puno * spno;
    rInternalSMSize = iInternalSMSize;
    rBigMatrixSize = iBigMatrixSize;
    rDensity = rInternalSMSize / rBigMatrixSize * 100;

    displayProgress1("%i conservation values counted, %i big matrix size, %g%% density of matrix \n",
                     iInternalSMSize,iBigMatrixSize,rDensity);
}

// read value for a single parameter specified in input.dat parameter file
template<class T>
void readInputOption(vector<string>& infile, string varname, T& value, int crit, int& present)
// Given a varname, scans the infile to see if this field was configured with primitive type T. 
// Here the lines of the file are supplied as a vector of strings.
// If configured, returns the value in "value". 
// crit = whether variable is mandatory
// present = stores whether variable was found.
{
    int foundit = 0;
    string modifiedVar = varname + " "; // pad with space to prevent fields that are  substrings of other fields

    for (string line : infile)
    {   /* loop through file looking for varname */
        // if line contains varname
        if (line.find(modifiedVar)!=string::npos) {
            // try to parse the value next to the variable by erasing first modifiedVar characters of line
            string varValue = line.erase(0, modifiedVar.size());
            utils::trim(varValue);

            if (varValue.empty()) { // variable defined but no value - keep looping.
                displayWarningMessage("WARNING: found bad/empty value for variable %s. Value ignored\n",varname.c_str());
                continue;
            }

            // read string representation to specified type T
            stringstream linestream(varValue);
            linestream >> value; 

            if (linestream.fail()) {
                // some error occurred when reading in value 
                displayErrorMessage("Invalid parameter type for given param %s. Expecting type: %s\n", varname.c_str(), typeid(value).name());
                continue;
            }

            foundit++;
        }
    }

    if (!foundit)
        if (crit)
           displayErrorMessage("Unable to find a valid value for %s in input file.\n",varname.c_str());

    if (foundit > 1)
        displayWarningMessage("WARNING variable: %s appears more than once in the input file. Final value taken\n",varname.c_str());

    present = foundit;
} // readInputOption

// read values for the parameters specified in input.dat parameter file
void readInputOptions(double &cm,double &prop,sanneal &anneal,
                      int &iseed,
                      long int &repeats, string& savename, sfname &fname, string filename,
                      int &runopts,double &misslevel,int &heurotype,int &clumptype,
                      int &itimptype, int &verb,
                      double &costthresh,double &tpf1,double &tpf2)
{
    vector<string> fileLines;
    double version;
    int present;
    #ifdef DEBUGTRACEFILE
    char debugbuffer[200];
    #endif

    verbosity = 1; /* This enables local warning messages */
    /* Setup all of the default parameter variables */
    version = 0.1;
    cm = 0;
    prop = 0;
    anneal = {1, 0, 1, 0 ,1};
    iseed = -1;
    costthresh = 0;
    tpf1 = 0;
    tpf2 = 0;
    repeats = 0;
    asymmetricconnectivity = 0;

    savename = "temp";
    misslevel = 1;
    heurotype = 1;
    clumptype = 0;
    verbosity = 1;

    /* Open file and store all lines in vector */
    ifstream fp(filename.c_str());
    // Check if object is valid
    if(!fp)
    {
        displayErrorMessage("input file %s not found\nAborting Program.\n\n",filename.c_str());
    }

    string line;
    while (getline(fp, line))
    {
        // Store non empty lines
        if(line.size() > 0)
            fileLines.push_back(line);
    }
    fp.close();

    readInputOption(fileLines,"VERSION",version, 0,present);
    readInputOption(fileLines,"BLM",cm,0,present);
    if (present == 0)
        readInputOption(fileLines,"CM",cm,0,present);
    readInputOption(fileLines,"PROP",prop,0,present);

    readInputOption(fileLines,"RANDSEED",iseed,0,present);
    if (present==0)  /* The random seed. -1 to set by clock */
        displayErrorMessage("Error reading input parameter random seed.\n");

    /* Annealing Controls */
    readInputOption(fileLines,"NUMITNS",anneal.iterations,0,present);
    readInputOption(fileLines,"STARTTEMP",anneal.Tinit,0,present);
    readInputOption(fileLines,"COOLFAC",anneal.Tcool,0,present);
    readInputOption(fileLines,"NUMTEMP",anneal.Titns,0,present);

    anneal.type = 1;
    if (anneal.iterations < 1 )
        anneal.type = 0;
    if (anneal.Tinit < 0)
        anneal.type = (int) (-anneal.Tinit) + 1;  /* type is negative of Tinit */

    /* Various controls */
    readInputOption(fileLines,"NUMREPS",repeats, 0,present);
    readInputOption(fileLines,"COSTTHRESH",costthresh,0,present);
    readInputOption(fileLines,"THRESHPEN1",tpf1,0,present);
    readInputOption(fileLines,"THRESHPEN2",tpf2,0,present);

    /* SaveFiles */
    readInputOption(fileLines,"SCENNAME",savename,0,present);

    /* SaveFiles New Method */
    readInputOption(fileLines,"SAVERUN",fnames.saverun,0,present);
    readInputOption(fileLines,"SAVEBEST",fnames.savebest,0,present);
    readInputOption(fileLines,"SAVESUMMARY",fnames.savesum,0,present);
    readInputOption(fileLines,"SAVESCEN",fnames.savesen,0,present);
    readInputOption(fileLines,"SAVETARGMET",fnames.savespecies,0,present);
    readInputOption(fileLines,"SAVESUMSOLN",fnames.savesumsoln,0,present);
    readInputOption(fileLines,"SAVESPEC",fnames.savespec,0,present);
    readInputOption(fileLines,"SAVEPU",fnames.savepu,0,present);
    readInputOption(fileLines,"SAVEMATRIXPUORDER",fnames.savepuvspr,0,present);
    readInputOption(fileLines,"SAVEMATRIXSPORDER",fnames.savematrixsporder,0,present);
    readInputOption(fileLines,"SAVEPENALTY",fnames.savepenalty,0,present);
    readInputOption(fileLines,"SAVETOTALAREAS",fnames.savetotalareas,0,present);
    readInputOption(fileLines,"SAVEDEBUGTRACEFILE",fnames.savedebugtracefile,0,present);
    readInputOption(fileLines,"SAVERICHNESS",fnames.saverichness,0,present);
    readInputOption(fileLines,"SAVESOLUTIONSMATRIX",fnames.savesolutionsmatrix,0,present);
    readInputOption(fileLines,"SOLUTIONSMATRIXHEADERS",fnames.solutionsmatrixheaders,0,present);
    readInputOption(fileLines,"SAVELOG",fnames.savelog,0,present);
    readInputOption(fileLines,"SAVESNAPSTEPS",fnames.savesnapsteps,0,present);
    readInputOption(fileLines,"SAVESNAPCHANGES",fnames.savesnapchanges,0,present);
    readInputOption(fileLines,"SAVESNAPFREQUENCY",fnames.savesnapfrequency,0,present);
    readInputOption(fileLines,"SAVEANNEALINGTRACE",fnames.saveannealingtrace,0,present);
    readInputOption(fileLines,"ANNEALINGTRACEROWS",fnames.annealingtracerows,0,present);
    readInputOption(fileLines,"SAVEITIMPTRACE",fnames.saveitimptrace,0,present);
    readInputOption(fileLines,"ITIMPTRACEROWS",fnames.itimptracerows,0,present);
    readInputOption(fileLines,"RIMAGETYPE",fnames.rimagetype,0,present);
    readInputOption(fileLines,"REXECUTESCRIPT",fnames.rexecutescript,0,present);
    readInputOption(fileLines,"RCLUSTERCOUNT",fnames.rclustercount,0,present);
    readInputOption(fileLines,"RIMAGEWIDTH",fnames.rimagewidth,0,present);
    readInputOption(fileLines,"RIMAGEHEIGHT",fnames.rimageheight,0,present);
    readInputOption(fileLines,"RIMAGEFONTSIZE",fnames.rimagefontsize,0,present);
    readInputOption(fileLines,"ASYMMETRICCONNECTIVITY",asymmetricconnectivity,0,present);
    readInputOption(fileLines,"CONNECTIVITYIN",fOptimiseConnectivityIn,0,present);

    // quantum annealing control parameters
    // TODO ADBAI - enable when ready
    /*
    readInputOption(fileLines,"QAPROP",&rQAPROP,0,present);
    readInputOption(fileLines,"QADECAY",&rQADECAY,0,present);
    readInputOption(fileLines,"QADECAYB",&rQADECAYB,0,present);
    readInputOption(fileLines,"QADECAYTYPE",&iQADECAYTYPE,0,present);
    readInputOption(fileLines,"QAACCPR",&rQAACCPR,0,present);
    */

    if (!fnames.savesnapfrequency)
        fnames.savesnapfrequency = 1;

    /* Filenames */
    readInputOption(fileLines,"INPUTDIR",fnames.inputdir,1,present);
    fnames.inputdir = utils::cleanDirectoryString(fnames.inputdir);
    readInputOption(fileLines,"OUTPUTDIR",fnames.outputdir,1,present);
    fnames.outputdir = utils::cleanDirectoryString(fnames.outputdir);
    readInputOption(fileLines,"PUNAME",fnames.puname,1,present);
    if (!present)
        fnames.puname = "PU.dat";
    readInputOption(fileLines,"SPECNAME",fnames.specname,1,present);
    if (!present)
        fnames.puname = "spec.dat";
    readInputOption(fileLines,"PUVSPRNAME",fnames.puvsprname,1,present);
    if (!present)
        fnames.puname = "puvspr2.dat";

    readInputOption(fileLines,"MATRIXSPORDERNAME",fnames.matrixspordername,0,present);
    readInputOption(fileLines,"PENALTYNAME",fnames.penaltyname,0,present);
    readInputOption(fileLines,"BOUNDNAME",fnames.connectionname,0,present);
    if (present == 0)
        readInputOption(fileLines,"CONNECTIONNAME",fnames.connectionname,0,present);

    readInputOption(fileLines,"CONNECTIONFILESNAME",fnames.connectionfilesname,0,present);
    readInputOption(fileLines,"BLOCKDEFNAME",fnames.blockdefname,0,present);

    readInputOption(fileLines,"BESTFIELDNAME",fnames.bestfieldname,0,present);
    if (!present)
        fnames.bestfieldname = "SOLUTION";

    readInputOption(fileLines,"RBINARYPATHNAME",fnames.rbinarypath,0,present);

    /* various other controls */
    readInputOption(fileLines,"RUNMODE",runopts,1,present);
    readInputOption(fileLines,"MISSLEVEL",misslevel,0,present);
    readInputOption(fileLines,"HEURTYPE",heurotype,0,present);
    readInputOption(fileLines,"CLUMPTYPE",clumptype,0,present);
    readInputOption(fileLines,"ITIMPTYPE",itimptype,0,present);
    readInputOption(fileLines,"VERBOSITY",verbosity,0,present);
    readInputOption(fileLines,"PROBABILITYWEIGHTING",rProbabilityWeighting,0,present);
    readInputOption(fileLines,"STARTDECTHRESH",rStartDecThresh,0,present);
    readInputOption(fileLines,"ENDDECTHRESH",rEndDecThresh,0,present);
    readInputOption(fileLines,"STARTDECMULT",rStartDecMult,0,present);
    readInputOption(fileLines,"ENDDECMULT",rEndDecMult,0,present);

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

    if (!fnames.outputdir.empty())
    {
        savename = fnames.outputdir + savename;
    }
} // readInputOptions

} // marxan