# Detect platform
ifeq ($(OS),Windows_NT)
SHELL := cmd
.SHELLFLAGS := /C
endif

ifeq ($(OS),Windows_NT)
    CXX = cl
    #CXXFLAGS = /std:c++17 /nologo /EHsc /W4 /Iinclude /fp:strict /O2 /Zi /Zo /Oy- /openmp:llvm
    CXXFLAGS = /std:c++17 /nologo /EHsc /W4 /Iinclude /O2 /arch:AVX2 /openmp:llvm /Zi /Zo /DNDEBUG
    MKDIR = if not exist "$(OBJ_DIR)" mkdir "$(OBJ_DIR)"
    OUT_EXT = .exe
    OBJ_EXT = .obj
    DEP_EXT = .d
    COMPILE = $(CXX) $(CXXFLAGS) /c $< /Fo$@
    #LINK = $(CXX) /nologo $(OBJECTS) /Fe:$(TARGET)$(OUT_EXT) /link /DEBUG:FULL /INCREMENTAL:NO
    LINK = $(CXX) /nologo $(OBJECTS) /Fe:$(TARGET)$(OUT_EXT) /link /LTCG /OPT:REF /OPT:ICF /INCREMENTAL:NO
    

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

# --- cross-platform clean ---
ifeq ($(OS),Windows_NT)
CLEAN_CMD = \
    if exist "$(OBJ_DIR)" rmdir /S /Q "$(OBJ_DIR)" && \
    if exist "$(TARGET)$(OUT_EXT)" del /Q "$(TARGET)$(OUT_EXT)" && \
    if exist "$(TARGET).pdb" del /Q "$(TARGET).pdb" && \
    if exist "$(TARGET).ilk" del /Q "$(TARGET).ilk"
else
CLEAN_CMD = rm -rf "$(OBJ_DIR)" "$(TARGET)$(OUT_EXT)"
endif

.PHONY: clean
clean:
	-$(CLEAN_CMD)