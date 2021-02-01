#pragma once

// Contains functions relating to marxan algorithmic execution, like sa and sqa
// Moving to a new file for easier unit testing. 
// All functions here should be unit and should not have dependencies on the core marxan files.

#include <algorithm>
#include <utility>

#include "algorithms.hpp"

namespace marxan {
    namespace algorithms {
        // sets a given vector to 0,1 depending on proportion. prop = proportion of 1s.
        // if pu.status exists then it takes precedence. 
        void initialiseReserve(double prop, const vector<spustuff>& pu, vector<int>& R, mt19937& rngEngine)
        {
            uniform_real_distribution<double> float_range(0.0, 1.0);
            for (int i = 0; i < R.size(); i++)
                if (pu[i].status)
                    R[i] = pu[i].status; // set planning unit status to pu.dat status
                else
                    R[i] = float_range(rngEngine) < prop ? 1 : 0;
        }

        // apply the species penalties nominated in input penalties file for use in the annealing algorithms
        void applyUserPenalties(vector<sspecies>& spec)
        {
            for (int i = 0; i < spec.size(); i++)
                spec[i].penalty = spec[i].rUserPenalty;
        }

        // compute tables for looking up pu's and species fast
        void populateLookupTables(int puno, int spno, vector<spustuff>& pu, vector<sspecies>& spec,
            map<int, int>& PULookup, map<int, int>& SPLookup)
        {
            /* create the lookup map for planning unit and species names */
            /* populate the lookup map with planning unit and species names*/
            for (int i = 0; i < puno; i++)
            {
                PULookup[pu[i].id] = i;
            }
            for (int i = 0; i < spno; i++)
            {
                SPLookup[spec[i].name] = i;
            }
        }

        // set default targets for species
        void setDefaultTargets(vector<sspecies>& spec)
        {
            for (int isp = 0; isp < spec.size(); isp++)
            {
                spec[isp].target = max(0.0, spec[isp].target);
                spec[isp].target2 = max(0.0, spec[isp].target2);
                spec[isp].targetocc = max(0, spec[isp].targetocc);
                spec[isp].sepnum = max(0, spec[isp].sepnum);
                spec[isp].sepdistance = max(0.0, spec[isp].sepdistance);

                if (spec[isp].spf < 0)
                    spec[isp].spf = 1;
            }
        }
    } // namespace algorithms
} // namespace marxan