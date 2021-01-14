# Marxan
In development 2020

# How to Build (windows)
Marxan has now been refactored to use c++17 for more modernised code and more extensive standard library. MinGW64 is needed to build on windows. Ideally we would compile with latest c++ (c++20 at the time of writing) but I could not find a windows build toolchain that was up to date with gcc 10.2 that actually worked. 

Command to build is 
```
g++ -O3 -std=c++17 -static -fopenmp marxan.cpp clumping.cpp heuristics.cpp input.cpp output.cpp probability.cpp -o bin/marxan
```

All libraries are statically linked so we ship 1 executable. 

# How to Build (Mac with Intel processor)

OpenMP library is needed to build on Mac. Set up steps using Homebrew:
```
g++
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install libomp
```
Command to build is: 
```
g++ marxan.cpp clumping.cpp heuristics.cpp input.cpp output.cpp probability.cpp -lomp -o bin/marxan -Xclang -fopenmp -std=c++17  -O3
```
Libraries are not statically linked. 

# How to run unit tests (Linux or WSL only)
CppuTest is used as the testing framework. Only compilable on linux or WSL. 

Command to build tests is:
```

```

