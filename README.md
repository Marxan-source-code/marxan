# Marxan

## Latest Release
### **Marxan 4 is here for you to try**
 
Marxan has been rebuilt in C++ and we are ready for you to test the new version 4.  
 
**Upgrades include:**  
- support for multithreading and parallelization (full utilisation of CPU)  
- improved computational speed and efficiency  
- additional error reporting  
 
Download for Windows / MacOS / Linux: see releases section. 

You can access the code for the new version in a separate branch of this repo: /marxan4  
We will officially announce the release of the new version and merge changes to main soon, and are looking for your feedback in the meantime.  
Please feel free to use our [Google Group](https://groups.google.com/g/marxan). You can also create GitHub issues on this repo or e-mail marxancloud@gmail.com  

# Releases
## v4.0.3
- Windows [x86-64](https://github.com/Marxan-source-code/marxan/releases/download/v4.0.3/build4.0.3Windows.zip)
- Linux [x86-64](https://github.com/Marxan-source-code/marxan/releases/download/v4.0.3/build4.0.3Linux.zip)
- MacOS 10.15 (Catalina) [x86-64](https://github.com/Marxan-source-code/marxan/releases/download/v4.0.3/Marxan-4.0.3-macOS.zip)
- MacOS 10.13 (High Sierra) [x86-64](https://github.com/Marxan-source-code/marxan/releases/download/v4.0.3/Marxan-4.0.3-MacOS-10.13.zip)
# Test Data
- From Google Drive: [MarxanData.zip](https://drive.google.com/file/d/1VGN4S5L_F80Ds2JlSy0MclxS5-KpkwlH/view?usp=sharing)
- From Releases: [MarxanData.zip](https://github.com/Marxan-source-code/marxan/releases/download/v4.0.3/MarxanData.zip)

# How to Build (Windows)
Marxan has now been refactored to use c++17 for more modernised code and more extensive standard library. MinGW64 is needed to build on windows. Ideally we would compile with latest c++ (c++20 at the time of writing) but I could not find a windows build toolchain that was up to date with gcc 10.2 that actually worked. 

Command to build is 
```
g++ -O3 -std=c++17 -static -fopenmp marxan.cpp clumping.cpp heuristics.cpp input.cpp output.cpp probability.cpp algorithms.cpp computation.cpp utils.cpp -o bin/marxan
```

All libraries are statically linked so we ship 1 executable. 

# How to Build (Mac with Intel processor)

OpenMP library is needed to build on Mac. Set up steps using Homebrew:
```
g++
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install libomp
```
Command to build without statically linked libraries is: 
```
g++ marxan.cpp clumping.cpp heuristics.cpp input.cpp output.cpp probability.cpp  algorithms.cpp computation.cpp utils.cpp -lomp -o bin/marxan -Xclang -fopenmp -std=c++17  -O3
```
Command to build with statically linked openMP library is: 
```
g++ marxan.cpp clumping.cpp heuristics.cpp input.cpp output.cpp probability.cpp  algorithms.cpp computation.cpp utils.cpp /usr/local/opt/libomp/lib/libomp.a -o bin/marxan -Xclang -fopenmp -std=c++17  -O3
```
# How to run (Mac with M1 processor)
Install translator from Intel architecture using terminal command:
```
softwareupdate --install-rosetta
```
Follow the steps for Mac with x86-64 processor.
# How to run (Mac with x86-64 processor)
From directory with marxan excutable, assign executable flag to marxan
```
chmod +x ./marxan
```
Try to run marxan:
```
./marxan
```
System will prevent it from running. 
Go to System preferences -> Security & Privacy -> Allow apps downloaded from:  
Allow to run marxan.

# How to run unit tests (Linux or WSL only)
CppuTest is used as the testing framework. Only compilable on linux or WSL. 

Command to build tests is:

- cd to the /test directory
- call "make" to compile the tests
- ./test to run


