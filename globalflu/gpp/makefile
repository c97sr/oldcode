# Set main options (can be at command line)
# options for SYSTEM = max | windows
# Assumes that source files are in ../src/.

ifndef $(SYSTEM) 
SYSTEM=mac
endif
ifndef $(MODE)
MODE=debug
endif

# Set main compile option
CXX = g++

# Set local paths
SRCLIB = ../src

# Set header paths and library paths for gsl if required
ifeq ($(SYSTEM),mac)
GSLHEADER = /opt/local/include
GSLLIB = /opt/local/lib
else ifeq ($(SYSTEM),windows)
GSLHEADER = /usr/include
GSLLIB = /lib
else ifeq ($(SYSTEM),linux)
GSLHEADER = /usr/include
GSLLIB = /lib
endif

INCLUDEDIRS = -I$(SRCLIB) -I$(GSLHEADER)
LDFLAGS =   -L$(GSLLIB) -lgsl -lgslcblas

# System variable to tell the compiler where to look for .cpps and .hs
VPATH = $(SRCLIB)

# Set standatd compile variables
ifeq ($(MODE),run)
CXXFLAGS = -O3 $(INCLUDEDIRS)
else ifeq ($(MODE),debug)
 
# Possible flags -Wall for all warnings
CXXFLAGS = -Wall -ggdb $(INCLUDEDIRS) 
endif

# Define the object files
OBJFILES = sim_aux.o

all : sim.exe test.exe

sim.exe : $(OBJFILES) sim.o sim.h sim_aux.cpp
	$(CXX) $(OBJFILES) sim.o -o $@ $(LDFLAGS) 

test.exe : $(OBJFILES) test.o sim.h sim_aux.cpp
	$(CXX) $(OBJFILES) sim.o -o $@ $(LDFLAGS) 

clean:
	rm *.o *.exe
