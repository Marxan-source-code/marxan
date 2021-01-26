
#include "../computation.hpp"
#include "CppUTest/TestHarness.h"

using namespace marxan;

TEST_GROUP(ComputationTestsGroup)
{
    // Ctrl+F the function name to find tests for the function.
};

// computeSepPenalty basic test
TEST(ComputationTestsGroup, computeSepPenalty_test)
{
    // Ensure 0 sep requirement returns 0 penalty
    CHECK_EQUAL(0, computeSepPenalty(111, 0));
}

// Test 1 specie
TEST(ComputationTestsGroup, computeRepresentationMISSLEVEL_test_simple_target)
{
    int spno = 1;
    vector<sspecies> spec(1);
    double shortfall = 0.0, rMinimumProportionMet = 0.0;

    // Define the species targets and amounts
    spec[0].target = 10;
    spec[0].amount = 5;

    // Simple case where proportion is actually 0.5, so modifying the misslevel should yield different counts
    CHECK_EQUAL(0, computeRepresentationMISSLEVEL(spno, spec, 0.25, shortfall, rMinimumProportionMet));
    CHECK_EQUAL(1, computeRepresentationMISSLEVEL(spno, spec, 0.75, shortfall, rMinimumProportionMet));
    CHECK_EQUAL(0.5, rMinimumProportionMet);
    CHECK_EQUAL(5, shortfall);
}

// Test 1 specie with both targetocc and target
TEST(ComputationTestsGroup, computeRepresentationMISSLEVEL_test_both_targets1)
{
    int spno = 1;
    vector<sspecies> spec(1);
    double shortfall = 0.0, rMinimumProportionMet = 0.0;

    // Define the species targets and amounts
    spec[0].target = 10;
    spec[0].amount = 5;
    spec[0].targetocc = 4;
    spec[0].occurrence = 1;

    CHECK_EQUAL(0, computeRepresentationMISSLEVEL(spno, spec, 0.25, shortfall, rMinimumProportionMet));
    CHECK_EQUAL(1, computeRepresentationMISSLEVEL(spno, spec, 0.75, shortfall, rMinimumProportionMet));
    CHECK_EQUAL(0.25, rMinimumProportionMet); // occ proportion is lower
    CHECK_EQUAL(8, shortfall);
}

// Test 1 specie with both targetocc and target but targetocc is not met
TEST(ComputationTestsGroup, computeRepresentationMISSLEVEL_test_both_targets2)
{
    int spno = 1;
    vector<sspecies> spec(1);
    double shortfall = 0.0, rMinimumProportionMet = 0.0;

    // Define the species targets and amounts
    spec[0].target = 10;
    spec[0].amount = 15;
    spec[0].targetocc = 4;
    spec[0].occurrence = 1;

    CHECK_EQUAL(1, computeRepresentationMISSLEVEL(spno, spec, 0.3, shortfall, rMinimumProportionMet));
    CHECK_EQUAL(0.25, rMinimumProportionMet); // occ proportion is lower
    CHECK_EQUAL(3, shortfall);
}

TEST(ComputationTestsGroup, computeConnectivityIndices_test)
{
    // TODO
}

TEST(ComputationTestsGroup, computePlanningUnitValue_test)
{
    sconnections connection;
    connection.fixedcost = 10.0;
    connection.nbrno = 2;

    CHECK_EQUAL(10.0, connectionCost1(connection, 1.0, 0));

    // construct neighbours
    sneighbour n1 = { 0, 2.0, 1 };
    sneighbour n2 = { 0, 4.4, 0 };
    connection.first.push_back(n1);
    connection.first.push_back(n2);

    // Test connection cost first
    CHECK_EQUAL(16.4, connectionCost1(connection, 1.0, 0));
    CHECK_EQUAL(12.0, connectionCost1(connection, 1.0, 1));

    // Test computePlanningUnitValue 
    spustuff pu;
    pu.cost = 7.0;
    CHECK_EQUAL(23.4, computePlanningUnitValue(pu, connection, 1.0, 0));
    CHECK_EQUAL(19.0, computePlanningUnitValue(pu, connection, 1.0, 1));
}

TEST(ComputationTestsGroup, returnAmountSpecAtPu_test)
{
    // Test case with 2 pu, 1 specie
    // pu0 has 50, pu1 has 80 of the specie. 
    vector<spu> SM;
    spu sm1;
    sm1.amount = 50.0;
    sm1.spindex = 0;
    spu sm2 = sm1;
    sm2.amount = 80.0;
    SM.push_back(sm1);
    SM.push_back(sm2);

    spustuff pu1;
    pu1.id = 0;
    pu1.richness = 1;
    pu1.offset = 0;
    pair<int, double> result1 = returnAmountSpecAtPu(pu1, SM, 0);
    pair<int, double> result2 = returnAmountSpecAtPu(pu1, SM, 1); // non existent specie

    CHECK_EQUAL(0, result1.first);
    CHECK_EQUAL(50.0, result1.second);
    CHECK_EQUAL(-1, result2.first);
    CHECK_EQUAL(0, result2.second);

    spustuff pu2 = pu1;
    pu2.offset = 1;
    pair<int, double> result3 = returnAmountSpecAtPu(pu2, SM, 0);

    CHECK_EQUAL(1, result3.first);
    CHECK_EQUAL(80.0, result3.second);
}

TEST(ComputationTestsGroup, ConnectionCost2_test_simple)
{
    int asymmetricconnectivity = 0, fOptimiseConnectivityIn = 0, ipu = 0;
    vector<sconnections> connections;
    vector<int> R;

    // Set up scenario with only 1 pu
    sconnections c1;
    c1.fixedcost = 50.0;
    c1.nbrno = 0;
    connections.push_back(c1);

    R.push_back(1);

    // Cost should only be fixed cost
    CHECK_EQUAL(50.0, ConnectionCost2(connections[ipu], R, 1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
    CHECK_EQUAL(-50.0, ConnectionCost2(connections[ipu], R, -1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));

    // with cm multiplier
    CHECK_EQUAL(750.0, ConnectionCost2(connections[ipu], R, 1, 1, 15.0, asymmetricconnectivity, fOptimiseConnectivityIn));

    // asymmetric
    asymmetricconnectivity = 1;
    CHECK_EQUAL(50.0, ConnectionCost2(connections[ipu], R, 1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
}

TEST(ComputationTestsGroup, ConnectionCost2_test_2_neighbours)
{
    int asymmetricconnectivity = 0, fOptimiseConnectivityIn = 0, ipu = 0;
    vector<sconnections> connections;
    vector<int> R(2, 1); // Both pu are switched on.

    // Set up scenario with 2 pu that are neighbours
    sconnections c1;
    c1.fixedcost = 50.0;
    c1.nbrno = 1;
    c1.first.push_back(sneighbour(1, 30.0, 1));
    connections.push_back(c1);

    sconnections c2;
    c2.fixedcost = 2.0;
    c2.nbrno = 1;
    c2.first.push_back(sneighbour(0, 30.0, 1)); // symmetric connectivity
    connections.push_back(c2);

    // The first cost should be fixed cost, but minus the neighbouring cost because imode2=1 (so we treat as if removing connection)
    // The second cost is similar but imode2=0 so we exclude connection cost. Therefore no change in cost
    // The third cost is removing a pu and including connection cost change
    // The behaviour of this function is undefined if imode2=0 and imode=-1
    CHECK_EQUAL(50.0 - 30.0, ConnectionCost2(connections[ipu], R, 1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
    CHECK_EQUAL(50.0, ConnectionCost2(connections[ipu], R, 1, 0, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
    CHECK_EQUAL(-50.0 + 30.0, ConnectionCost2(connections[ipu], R, -1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));

    // perform another set of checks but this time for ipu=1
    ipu = 1;
    CHECK_EQUAL(2.0 - 30.0, ConnectionCost2(connections[ipu], R, 1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
    CHECK_EQUAL(2.0, ConnectionCost2(connections[ipu], R, 1, 0, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
    CHECK_EQUAL(-2 + 30.0, ConnectionCost2(connections[ipu], R, -1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));

    // repeat but now pu 0 is removed
    R[0] = 0;
    CHECK_EQUAL(2.0 + 30.0, ConnectionCost2(connections[ipu], R, 1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
    CHECK_EQUAL(2.0 + 30.0, ConnectionCost2(connections[ipu], R, 1, 0, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
    CHECK_EQUAL(-2 - 30.0, ConnectionCost2(connections[ipu], R, -1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
}

// Asymmetric connectivity means we incur the cost if the destination pu is active, but origin isn't.
// This function tests asymmetricconnectivity cost for imode2=1
TEST(ComputationTestsGroup, ConnectionCost2_test_asymmetricconnectivity1) {
    int asymmetricconnectivity = 1, fOptimiseConnectivityIn = 0, ipu = 0;
    vector<sconnections> connections;
    vector<int> R(2, 1); // Both pu are switched on.

    // Set up scenario with 2 pu that are neighbours.
    // pu 0 is the origin.
    sconnections c1;
    c1.fixedcost = 50.0;
    c1.nbrno = 1;
    c1.first.push_back(sneighbour(1, 30.0, 0));
    connections.push_back(c1);

    sconnections c2;
    c2.fixedcost = 2.0;
    c2.nbrno = 1;
    c2.first.push_back(sneighbour(0, 30.0, 1));
    connections.push_back(c2);

    // In the case where both PU are turned on, the connection cost is not incurred
    // This is apparent for imode=1 as we are turning on the origin node, thereby removing the cost.
    CHECK_EQUAL(50.0 - 30.0, ConnectionCost2(connections[ipu], R, 1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
    CHECK_EQUAL(-50.0 + 30.0, ConnectionCost2(connections[ipu], R, -1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));

    // perform another set of checks but this time for ipu=1
    ipu = 1;
    // imode==1 should result in no connection change since origin (pu=0) was on already. 
    CHECK_EQUAL(2.0, ConnectionCost2(connections[ipu], R, 1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
    CHECK_EQUAL(-2.0, ConnectionCost2(connections[ipu], R, -1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));

    // repeat but now pu 0 is removed
    R[0] = 0;

    // Now that origin is removed, switching on pu1 will incur the connection cost
    CHECK_EQUAL(2.0 + 30.0, ConnectionCost2(connections[ipu], R, 1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
    CHECK_EQUAL(-2.0 - 30.0, ConnectionCost2(connections[ipu], R, -1, 1, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
}


// Asymmetric connectivity means we incur the cost if the destination pu is active, but origin isn't.
// This function tests asymmetricconnectivity cost for imode2=0
// Note that the behaviour of imode2=0 and imode1=-1 is NOT DEFINED! Therefore we do not test it.
TEST(ComputationTestsGroup, ConnectionCost2_test_asymmetricconnectivity2) {
    int asymmetricconnectivity = 1, fOptimiseConnectivityIn = 0, ipu = 0;
    vector<sconnections> connections;
    vector<int> R(2, 1); // Both pu are switched on.

    // Set up scenario with 2 pu that are neighbours.
    // pu 0 is the origin.
    sconnections c1;
    c1.fixedcost = 50.0;
    c1.nbrno = 1;
    c1.first.push_back(sneighbour(1, 30.0, 0));
    connections.push_back(c1);

    sconnections c2;
    c2.fixedcost = 2.0;
    c2.nbrno = 1;
    c2.first.push_back(sneighbour(0, 30.0, 1));
    connections.push_back(c2);

    // For imode=1 as we are turning on the origin node, thereby removing the cost. 
    // However imode2 is on so negative connections are not included.
    CHECK_EQUAL(50.0, ConnectionCost2(connections[ipu], R, 1, 0, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));

    // perform another set of checks but this time for ipu=1
    ipu = 1;
    // imode==1 should result in no connection change since origin (pu=0) was on already. 
    CHECK_EQUAL(2.0, ConnectionCost2(connections[ipu], R, 1, 0, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));

    // repeat but now pu 0 is removed
    R[0] = 0;

    // Now that origin is removed, switching on pu1 will incur the connection cost
    CHECK_EQUAL(2.0 + 30.0, ConnectionCost2(connections[ipu], R, 1, 0, 1.0, asymmetricconnectivity, fOptimiseConnectivityIn));
}

TEST(ComputationTestsGroup, computeSpecProp_test) {
    // Construct 2 specie
    vector<sspecies> spec(3);
    spec[0].prop = 0.1;
    spec[1].prop = 0.5;
    spec[2].prop = 0;

    // 3 pu with 2 species each
    vector<spustuff> pu(3);
    pu[0].richness = 2;
    pu[0].offset = 0;
    pu[1].richness = 2;
    pu[1].offset = 1;
    pu[2].richness = 2;
    pu[2].offset = 2;

    // Construct SM
    /*
    vector<map<int,spu>> SM(3);
    SM[0][0].amount = 100.0;
    SM[0][1].amount = 10.0;
    SM[1][0].amount = 50.0;
    SM[1][1].amount = 20.0;
    SM[2][0].amount = 5.0;
    SM[2][2].amount = 200.0;
    */
    //Temporal fix test compilaton
    vector<spu> SM(6);
    SM[0].amount = 100.0;
    SM[1].amount = 10.0;
    SM[2].amount = 50.0;
    SM[3].amount = 20.0;
    SM[4].amount = 5.0;
    SM[5].amount = 200.0;

    computeSpecProp(3, spec, 3, pu, SM);

    CHECK_EQUAL(15.5, spec[0].target);
    CHECK_EQUAL(15, spec[1].target);
    CHECK_EQUAL(0, spec[2].target);
}

// This function is a fixed value, we just need to check the 0 cases.
TEST(ComputationTestsGroup, computeSepPenalty_test_zero) {
    CHECK_EQUAL(0.0, computeSepPenalty(1, 0));
    CHECK(computeSepPenalty(0, 1) != 0);
}
