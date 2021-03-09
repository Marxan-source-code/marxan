#COMPILER=g++
#FLAGS = -O3 -std=c++17 -fopenmp -Wall
BIN_DIR = bin/

all: directories marxan

directories:
	mkdir -p $(BIN_DIR)

marxan: marxan.cpp
	g++ -O3 -std=c++17 -static -fopenmp *.cpp -o bin/marxan

