# Set main options (can be at command line)
# options for SYSTEM = mac | windows | linux
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
SRCLIB = ../cpp

# Set header paths and library paths for gsl if required
ifeq ($(SYSTEM),mac)
GSLHEADER = /usr/include
RHEADER = /Library/Frameworks/R.framework/Resources/include
GSLLIB = /usr/lib
RCPPGSLHEADER = /Library/Frameworks/R.framework/Versions/3.0/Resources/library/RcppGSL/include
RCPPHEADER = /Library/Frameworks/R.framework/Versions/3.0/Resources/library/Rcpp/include
else ifeq ($(SYSTEM),windows)
GSLHEADER = /usr/include
GSLLIB = /lib
RHEADER = /Library/Frameworks/R.framework/Resources/include
else ifeq ($(SYSTEM),linux)
GSLHEADER = /usr/include
GSLLIB = /lib
RHEADER = /Library/Frameworks/R.framework/Resources/include
endif

INCLUDEDIRS = -I$(SRCLIB) -I$(GSLHEADER) -I$(RHEADER) -I$(RCPPGSLHEADER) -I$(RCPPHEADER)
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

all : parameters_test.exe sim.exe test.exe 

# You need to tell it about the .o s and the .hs. It figures out the .cpp s from the .o s
parameters_test.exe : parameters_test.o parameters.o utility.o parameters.h utility.h
	$(CXX) utility.o parameters_test.o parameters.o -o $@ $(LDFLAGS) 

# You need to tell it about the .o s and the .hs. It figures out the .cpp s from the .o s
asciipop_test.exe : asciipop_test.o asciipop.o utility.o asciipop.h parameters.o parameters.h utility.h
	$(CXX) utility.o asciipop_test.o asciipop.o parameters.o -o $@ $(LDFLAGS) 

# You need to tell it about the .o s and the .hs. It figures out the .cpp s from the .o s
farms_test.exe : farms_test.o farms.h farms.o
	$(CXX) farms_test.o farms.o -o $@ $(LDFLAGS) 

sim.exe : $(OBJFILES) sim.o sim.h sim_aux.cpp
	$(CXX) $(OBJFILES) sim.o -o $@ $(LDFLAGS) 

test.exe : $(OBJFILES) test.o sim.h sim_aux.cpp
	$(CXX) $(OBJFILES) sim.o -o $@ $(LDFLAGS) 

clean:
	rm *.o *.exe
