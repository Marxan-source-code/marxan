
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

TEST(UtilsTestsGroup, guess_delimeter_test)
{
    CHECK_EQUAL(guess_delimeter("1.0000\t12.0000\t0.0003"), '\t')
    CHECK_EQUAL(guess_delimeter("1,3,0.202500000596"), ',')
    CHECK_EQUAL(guess_delimeter("id,prop,target,spf,name"), ',')
    CHECK_EQUAL(guess_delimeter("2;0,2;1;SG"), ';')
    CHECK_EQUAL(guess_delimeter("2,0.2,1,SG\tUSA"), ',')
}

TEST(UtilsTestsGroup, formatted_string_stream_test)
{
    formatted_string_stream ss("2,0.2,1,SG", ',');
    int i1, i2;
    double d1;
    std::string s1;
    ss>>i1>>d1>>i2>>s1; 
    CHECK_EQUAL(i1,2);
    CHECK_TRUE(std::abs(d1-0.2)<0.00001);
    CHECK_EQUAL(i2,1);
    CHECK_TRUE(s1 == "SG");
}

TEST(UtilsTestsGroup, get_tokens_test)
{
    std::vector<std::string> v1 = { "id","prop","target","spf","name" };
    CHECK_TRUE( get_tokens("id,prop,target,spf,name", ',') ==  v1 );
    CHECK_TRUE( get_tokens("id;prop;target;spf;name", ';') == v1 );
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