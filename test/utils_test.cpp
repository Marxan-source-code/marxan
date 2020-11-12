
#include "../utils.hpp"
#include "CppUTest/TestHarness.h"

using namespace marxan::utils;

TEST_GROUP(UtilsTestsGroup)
{
};

// Test mode and suffix align
TEST(UtilsTestsGroup, getFileNameSuffix_test)
{
    CHECK_EQUAL(".csv", getFileNameSuffix(3));
    CHECK_EQUAL(".txt", getFileNameSuffix(2));
    CHECK_EQUAL(".dat", getFileNameSuffix(0)); // all else
    CHECK_EQUAL(".dat", getFileNameSuffix(1));
    CHECK_EQUAL(".dat", getFileNameSuffix(4));
}

TEST(UtilsTestsGroup, trim_test)
{
    std::string original = " original s ";
    trim(original);
    CHECK_EQUAL("original s", original);
}