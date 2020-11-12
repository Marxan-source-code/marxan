
#include "../computation.hpp"
#include "CppUTest/TestHarness.h"

using namespace marxan;

TEST_GROUP(ComputationTestsGroup)
{
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
   sneighbour n1 = {0, 2.0, 1};
   sneighbour n2 = {0, 4.4, 0};
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
   pair<int,double> result1 = returnAmountSpecAtPu(pu1, SM, 0);
   pair<int,double> result2 = returnAmountSpecAtPu(pu1, SM, 1); // non existent specie

   CHECK_EQUAL(0, result1.first);
   CHECK_EQUAL(50.0, result1.second);
   CHECK_EQUAL(-1, result2.first);
   CHECK_EQUAL(0, result2.second);

   spustuff pu2 = pu1;
   pu2.offset = 1;
   pair<int,double> result3 = returnAmountSpecAtPu(pu2, SM, 0);

   CHECK_EQUAL(1, result3.first);
   CHECK_EQUAL(80.0, result3.second);
}

