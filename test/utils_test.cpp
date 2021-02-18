
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

TEST(UtilsTestsGroup, is_like_numerical_data_test)
{
    CHECK_TRUE(is_like_numerical_data("1.0000,   12.0000,    0.0003"))
    CHECK_TRUE(is_like_numerical_data("1,3,0.202500000596"))
    CHECK_FALSE(is_like_numerical_data("id,prop,target,spf,name"))
    CHECK_FALSE(is_like_numerical_data("2,0.2,1,SG"))
}

TEST(UtilsTestsGroup, intToPaddedString_test)
{
    CHECK_EQUAL("01", intToPaddedString(1,2));
    CHECK_EQUAL("0100", intToPaddedString(100,4));
    CHECK_EQUAL("100", intToPaddedString(100,2));
    CHECK_EQUAL("1", intToPaddedString(1,1));
    CHECK_EQUAL("0000056", intToPaddedString(56,7));
    CHECK_EQUAL("56", intToPaddedString(56,0));
    CHECK_EQUAL("1", intToPaddedString(1,0));
    CHECK_EQUAL("0", intToPaddedString(0,0));
}