# Apache Ant.
ANT := ant

# Fortran compiler.
CC := gfortran

# Source code.
SOURCE_DIR := src/fortran

# Binary directory.
BUILD_DIR := bin

# Source files.
FREDMETHOD_SOURCE := binningdanny.f90 correctpcr.f90 divergencematrix.f90 fredMethod.f90 readsynec.f90 removegaps.f90
HILLCLIMBING_SOURCE := hillclimb.f90 driftCI.f90 npopCI.f90 omegaCI.f90 sigmaCI.f90 demarcationsCI.f90
NELDER_MEAD_SOURCE := NelderMead.f90
METHODS_SOURCE := Methods.f90


# Determine the operating system.
OS := ${shell uname -s}

# Set operating system specific variables.
ifeq (${OS}, Linux)
  ifeq (${shell getconf LONG_BIT}, 64)
    BINARY_EXT := .amd64
  else
    BINARY_EXT := .i386
  endif
else ifeq (${OS}, Darwin)
  BINARY_EXT := .app
else ifneq (,${findstring windows, ${OS}})
  BINARY_EXT := .exe
else ifneq (,${findstring CYGWIN, ${OS}})
  BINARY_EXT := .exe
endif

# Create a list of binary files to build.
FREDMETHOD_BINARIES := $(patsubst %.f90, $(BUILD_DIR)/%$(BINARY_EXT), $(FREDMETHOD_SOURCE))
HILLCLIMBING_BINARIES := $(patsubst %.f90, $(BUILD_DIR)/%$(BINARY_EXT), $(HILLCLIMBING_SOURCE))

# Create a list of object files to build.
NELDER_MEAD_OBJECT := $(patsubst %.f90, $(BUILD_DIR)/%.o, $(NELDER_MEAD_SOURCE))
METHODS_OBJECT := $(patsubst %.f90, $(BUILD_DIR)/%.o, $(METHODS_SOURCE))

# List of phony build targets.
.PHONY: build pre-build post-build clean java fortran

# Build the Java and Fortran code.
build: pre-build java fortran post-build

# Pre build.
pre-build: $(BUILD_DIR)

# Post build.
post-build:
	rm -f $(METHODS_OBJECT) $(NELDER_MEAD_OBJECT)

# Clean the build directory.
clean:
	${ANT} clean
	rm -f $(FREDMETHOD_BINARIES) $(HILLCLIMBING_BINARIES) $(METHODS_OBJECT) $(NELDER_MEAD_OBJECT)

# Build the Java portion of the program.
java:
	${ANT}

# Build the Fortran portion of the program.
fortran: $(FREDMETHOD_BINARIES) $(HILLCLIMBING_BINARIES)

# Build the fredmethod binaries.
$(FREDMETHOD_BINARIES) : $(BUILD_DIR)/%$(BINARY_EXT) : $(SOURCE_DIR)/%.f90 $(METHODS_OBJECT)
	$(CC) $^ $(CFLAGS) -o $@

# Build the hillclimbing binaries.
$(HILLCLIMBING_BINARIES) : $(BUILD_DIR)/%$(BINARY_EXT) : $(SOURCE_DIR)/%.f90 $(METHODS_OBJECT) $(NELDER_MEAD_OBJECT)
	$(CC) $^ $(CFLAGS) -o $@

# Build the methods object.
$(METHODS_OBJECT) : $(BUILD_DIR)/%.o : $(SOURCE_DIR)/%.f90
	$(CC) $^ -c $(CFLAGS) -o $@

# Build the nelder mead object.
$(NELDER_MEAD_OBJECT) : $(BUILD_DIR)/%.o : $(SOURCE_DIR)/%.f90
	$(CC) $^ -c $(CFLAGS) -o $@

# Make sure the build directory is there.
$(BUILD_DIR):
	mkdir $(BUILD_DIR)
