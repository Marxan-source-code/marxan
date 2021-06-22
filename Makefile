#COMPILER=g++
#FLAGS = -O3 -std=c++17 -fopenmp -Wall
TARGET_EXEC := marxan
BUILD_DIR := ./bin

all: $(BUILD_DIR)/$(TARGET_EXEC)

$(BUILD_DIR)/marxan: marxan.cpp
	mkdir -p $(dir $@)
	g++ -O3 -std=c++17 -static -fopenmp *.cpp -o $(BUILD_DIR)/marxan
