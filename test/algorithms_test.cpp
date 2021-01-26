
#include "../algorithms.hpp"
#include "CppUTest/TestHarness.h"

using namespace marxan;
using namespace marxan::algorithms;

TEST_GROUP(AlgorithmsTestsGroup)
{
};

// test init reserve sets vector correctly
TEST(AlgorithmsTestsGroup, initialiseReserve_test)
{
    mt19937 rngEngine(1); //arbitrary seed
    std::vector<int> v;
    std::vector<spustuff> pu;
    initialiseReserve(0.5, pu, v, rngEngine);
    CHECK(v.empty());

    v.resize(5);
    pu.resize(5);
    for (spustuff& p : pu) {
        p.status = 0;
    }
    initialiseReserve(0.5, pu, v, rngEngine);
    for (int term : v) {
        CHECK(term == 1 || term == 0);
    }

    // check proportion working
    initialiseReserve(1.01, pu, v, rngEngine);
    for (int term : v) {
        CHECK(term == 1);
    }

    initialiseReserve(0, pu, v, rngEngine);
    for (int term : v) {
        CHECK(term == 0);
    }

    // check pu override working
    for (spustuff& p : pu) {
        p.status = 2;
    }
    initialiseReserve(0.5, pu, v, rngEngine);
    for (int term : v) {
        CHECK(term == 2);
    }
}
