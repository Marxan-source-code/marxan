#pragma once

#include <string>
#include <vector>


namespace marxan {
    using namespace std;

    typedef struct sclumps
    {
        int clumpid;
        double amount;
        int occs;
        vector<int> head;
    } sclumps; /* Clump nodes for species Clump Structure */

    typedef struct sspecies
    {
        int name;
        int type;
        string sname;
        double target;
        double prop;
        int targetocc;
        double spf;
        double penalty;
        double rUserPenalty;
        double amount;
        double expected1D, expected2D, variance1D, variance2D;
        int occurrence;
        double sepdistance;
        int sepnum;
        int separation;
        int clumps;
        double target2;  // Only clumping species need this
        vector<sclumps> head;  // needed for clumping species
        int richness, offset;
        double Zscore1D, Zscore2D;
        double probability1D, probability2D;
        double ptarget1d, ptarget2d;
    }sspecies;


    typedef struct sgenspec
    {
        int type;
        int targetocc;
        double target;
        double target2;
        int sepnum;
        double sepdistance;
        double prop;
        double spf;
    } sgenspec;


} // namespace marxan