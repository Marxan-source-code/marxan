CXXFLAGS ?= -O3 -std=c++17 -fopenmp
TARGET_EXEC := marxan
BUILD_DIR := ./bin

all: $(BUILD_DIR)/$(TARGET_EXEC)

$(BUILD_DIR)/marxan: marxan.cpp
	mkdir -p $(dir $@)
	$(CXX) -static $(CXXFLAGS) $(CFLAGS) *.cpp -o $(BUILD_DIR)/$(TARGET_EXEC)

.PHONY: clean
clean:
	rm -r $(BUILD_DIR)
