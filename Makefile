prefix		?= /usr
bindir		?= $(prefix)/bin

CXXFLAGS ?= -O3 -std=c++17 -fopenmp
TARGET_EXEC := marxan
BUILD_DIR := ./bin

all: $(BUILD_DIR)/$(TARGET_EXEC)

$(BUILD_DIR)/$(TARGET_EXEC): marxan.cpp
	mkdir -p $(dir $@)
	$(CXX) -static $(CXXFLAGS) $(CFLAGS) *.cpp -o $(BUILD_DIR)/$(TARGET_EXEC)

install: $(BUILD_DIR)/$(TARGET_EXEC)
	install -d $(DESTDIR)/$(bindir)
	install -m 755 $(BUILD_DIR)/$(TARGET_EXEC) $(DESTDIR)/$(bindir)

clean:
	rm -r $(BUILD_DIR)

.PHONY: clean install
