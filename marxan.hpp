#pragma once

#include <cstdarg>
#include <chrono>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "species.hpp"
#include "connections.hpp"
#include "spu.hpp"
#include "anneal.hpp"
#include "input.hpp"
#include "options.hpp"
#include "algorithms.hpp"


namespace marxan {
    using namespace std;
    void setBlockDefinitions(int gspno, int spno, int puno, vector<sgenspec>& gspec, vector<sspecies>& spec, vector<spustuff>& PU, vector<spu>& SM);
    void setDefaultRunOptions(int runopts, srunoptions& runoptions);
    void secondaryExit(void);

} // namespace marxan