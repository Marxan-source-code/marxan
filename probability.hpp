#pragma once

#include <map>
#include <string>
#include <vector>

#include "species.hpp"
#include "spu.hpp"

namespace marxan {
    using namespace std;

    double ComputeProbability1D(const vector<double>& ExpectedAmount1D, const vector<double>& VarianceInExpectedAmount1D,
        int spno, vector<sspecies>& spec);
    double ComputeProbability2D(vector<double>& ExpectedAmount2D, vector<double>& VarianceInExpectedAmount2D,
        int spno, vector<sspecies>& spec);
    double ChangeProbability1D(int iIteration, int ipu, int spno, int puno, vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, int imode);
    double ChangeProbability2D(int iIteration, int ipu, int spno, int puno, vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, int imode);
    void ComputeP_AllPUsSelected_1D(const string& savename, int puno, int spno, const vector<spustuff>& pu, const vector<spu>& SM, const vector<sspecies>& spec);
    void ComputeP_AllPUsSelected_2D(const string& savename, int puno, int spno, const vector<spustuff>& pu, const vector<spu>& SM, const vector<sspecies>& spec);
    void ReturnProbabilityAmounts1D(vector<double>& ExpectedAmount1D, vector<double>& VarianceInExpectedAmount1D, int ipu,
        int puno, const vector<spustuff>& pu, const vector<spu>& SM);
    void ReturnProbabilityAmounts2D(vector<double>& ExpectedAmount2D, vector<double>& VarianceInExpectedAmount2D, int ipu,
        int puno, const vector<spustuff>& pu, const vector<spu>& SM);

} // namespace marxan