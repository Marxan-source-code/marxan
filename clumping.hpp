#pragma once

#include <vector>
#include "species.hpp"
#include "spu.hpp"
#include "connections.hpp"


namespace marxan {
    using namespace std;

    typedef struct sseplist {
        int size;
        vector<int> head;
    } sseplist;

    // clumping
    void AddNewPU(int ipu, int isp, const vector<sconnections>& connections, vector<sspecies>& spec, const vector<spustuff>& pu,
        const vector<spu>& SM, vector<spu_out>& SM_out, int clumptype);
    void RemPu(int ipu, int isp, const vector<sconnections>& connections, vector<sspecies>& spec, const vector<spustuff>& pu,
        const vector<spu>& SM, vector<spu_out>& SM_out, int clumptype);
    int CalcPenaltyType4(int isp, int puno, const vector<spu>& SM, vector<spu_out>& SM_out, const vector<sconnections>& connections,
        vector<sspecies>& spec, const vector<spustuff>& pu, double cm, int clumptype);
    int CountSeparation(int isp, const vector<sclumps>& newno,
        const vector<spustuff>& pu, const vector<spu>& SM, const vector<spu_out>& SM_out, const vector<sspecies>& spec, int imode);
    int CountSeparation2(int isp, int ipu, const vector<sclumps>& newno, int puno, const vector<int>& R,
        const vector<spustuff>& pu, const vector<spu>& SM, const vector<spu_out>& SM_out, const vector<sspecies>& spec, int imode);
    vector<int> makelist(int isp, int ipu, int puno, const vector<int>& R, const vector<sclumps>& newno, const vector<sspecies>& spec,
        const vector<spustuff>& pu, const vector<spu>& SM, const vector<spu_out>& SM_out, int imode);
    double NewPenalty4(int ipu, int isp, int puno, const vector<sspecies>& spec, const vector<spustuff>& pu, const vector<spu>& SM, const vector<spu_out>& SM_out,
        const vector<int>& R, const vector<sconnections>& connections, int imode, int clumptype);
    double PartialPen4(int isp, double amount, const vector<sspecies>& spec, int clumptype);
    int SepDealList(const vector<int>& head, vector<sseplist>& Dist, const vector<spustuff>& pu,
        const vector<sspecies>& spec, int first, int sepnum, double targetdist, int isp);
    void SetSpeciesClumps(int puno, const vector<int>& R, vector<sspecies>& spec, const vector<spustuff>& pu,
        const vector<spu>& SM, vector<spu_out>& SM_out, const vector<sconnections>& connections, int clumptype);
    void SpeciesAmounts(int spno, int puno, vector<sspecies>& spec, const vector<spustuff>& pu, vector<spu>& SM,
        vector<int>& R, int clumptype);
    void SpeciesAmounts4(int isp, vector<sspecies>& spec, int clumptype);
} // namespace marxan