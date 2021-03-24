# Marxan

## Latest Release
### **Marxan 4**
 
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
## v4.0.5
- Windows [x86-64](https://github.com/Marxan-source-code/marxan/releases/download/v4.0.5/marxan4.0.5-windows.zip)
- Linux [x86-64](https://github.com/Marxan-source-code/marxan/releases/download/v4.0.5/marxan4.0.5-linux.zip)
- MacOS [x86-64](https://github.com/Marxan-source-code/marxan/releases/download/v4.0.5/marxan-4.0.5-MacOS-10.15-x86-64.zip) (Created in MacOS 10.15 Catalina)
- MacOS [M1](https://github.com/Marxan-source-code/marxan/releases/download/v4.0.5/marxan-4.0.5-MacOS-11-M1.zip) (Created in MacOS 11 Big Sur)
# Test Data
- From Google Drive: [MarxanData.zip](https://drive.google.com/file/d/1VGN4S5L_F80Ds2JlSy0MclxS5-KpkwlH/view?usp=sharing)
- From Releases: [MarxanData.zip](https://github.com/Marxan-source-code/marxan/releases/download/v4.0.3/MarxanData.zip)

# How to Run (Windows)
Instructions on how to run are similar in Linux (minus the .bat file). Ensure your input files are in the structure:

```
root/
- /input (contains pu.dat, spec.dat ..etc.)
- /output (empty directory)
```

1. Copy the Marxan_x64 executable and bat file (if windows) into the root of your data folder, outside of /input and /output
2. Run the .bat file (if windows) or double click on the executable. If you want the terminal to stay open upon finishing, calling from the .bat file is needed.

# How to Build (Windows)
Marxan has now been refactored to use c++17 for more modernised code and more extensive standard library. MinGW64 is needed to build on windows. Ideally we would compile with latest c++ (c++20 at the time of writing) but I could not find a windows build toolchain that was up to date with gcc 10.2 that actually worked. 

Command to build is 
```
g++ -O3 -std=c++17 -static -fopenmp *.cpp -o bin/Marxan_x64
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
g++ *.cpp -lomp -o bin/Marxan_x64 -Xclang -fopenmp -std=c++17  -O3
```
Command to build with statically linked openMP library is: 
```
g++ *.cpp  /usr/local/opt/libomp/lib/libomp.a -o bin/Marxan_x64 -Xclang -fopenmp -std=c++17  -O3
```
# How to run (Mac with M1 processor)
M1 build runs noticable faster.
To be able to run x-86 builds, install translator from Intel architecture using terminal command:
```
softwareupdate --install-rosetta
```
Follow the steps for Mac with x86-64 processor.
# How to run (Mac with x86-64 processor)
From directory with marxan excutable, assign executable flag to marxan
```
chmod +x ./Marxan_x64
```
Try to run marxan:
```
./Marxan_x64
```
If system prevents it from running due to Unidentified Developer. 
Go to System preferences -> Security & Privacy -> Allow apps downloaded from:  
Allow to run marxan.

# How to run unit tests (Linux or WSL only)
CppuTest is used as the testing framework. Only compilable on linux or WSL. 

Command to build tests is:

- cd to the /test directory
- call "make" to compile the tests
- ./test to run


