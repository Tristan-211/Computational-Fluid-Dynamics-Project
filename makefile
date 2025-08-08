# Detect platform
ifeq ($(OS),Windows_NT)
    # Use MSVC on Windows
    CXX = cl
    CXXFLAGS = /std:c++17 /nologo /EHsc /W4 /Iinclude /fp:strict
    RM = rm -rf
    MKDIR = mkdir -p $(OBJ_DIR)
    OUT_EXT = .exe
    OBJ_EXT = .obj
    DEP_EXT = .d
    COMPILE = $(CXX) $(CXXFLAGS) /c $< /Fo$@
    LINK = $(CXX) $(CXXFLAGS) /Fe:$(TARGET)$(OUT_EXT) $(OBJECTS)
else
    # Use g++ on Linux/macOS
    CXX = g++
    CXXFLAGS = -std=c++17 -Wall -MMD -MP -g -Iinclude -O2 -ffloat-store -fexcess-precision=standard
    RM = rm -rf
    MKDIR = mkdir -p $(OBJ_DIR)
    OUT_EXT =
    OBJ_EXT = .o
    DEP_EXT = .d
    COMPILE = $(CXX) $(CXXFLAGS) -c $< -o $@
    LINK = $(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)
endif

SRC_DIR = src
OBJ_DIR = build

SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%$(OBJ_EXT),$(SOURCES))
DEPS := $(OBJECTS:$(OBJ_EXT)=.$(DEP_EXT))

TARGET = cfd_project

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINK)

$(OBJ_DIR)/%$(OBJ_EXT): $(SRC_DIR)/%.cpp
	$(MKDIR)
	$(COMPILE)

# Only include dependency files when using GCC (not MSVC)
ifeq ($(OS),Windows_NT)
else
-include $(DEPS)
endif

clean:
	$(RM) $(OBJ_DIR) $(TARGET)$(OUT_EXT)

.PHONY: all clean
