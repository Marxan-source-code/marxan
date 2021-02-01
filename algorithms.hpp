#pragma once

// Contains functions relating to marxan algorithmic execution, like sa and sqa
// All functions here should be unit and should not have dependencies on the core marxan files.

#include <utility>
#include <vector>
#include <map>
#include <random>
#include "species.hpp"
#include "spu.hpp"


namespace marxan {
    namespace algorithms {
        // sets a given vector to 0,1 depending on proportion. prop = proportion of 1s.
        // if pu.status exists then it takes precedence. 
        void initialiseReserve(double prop, const vector<spustuff>& pu, vector<int>& R, mt19937& rngEngine);

        // apply the species penalties nominated in input penalties file for use in the annealing algorithms
        void applyUserPenalties(vector<sspecies>& spec);

        // compute tables for looking up pu's and species fast
        void populateLookupTables(int puno, int spno, vector<spustuff>& pu, vector<sspecies>& spec,
            map<int, int>& PULookup, map<int, int>& SPLookup);

        // set default targets for species
        void setDefaultTargets(vector<sspecies>& spec);
    } // namespace algorithms
} // namespace marxan