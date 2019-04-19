SHELL := bash
CXX := g++
CPPFLAGS := -std=c++14 -Iinclude
CXXFLAGS := -Wall -O3 -flto -fmax-errors=3 $(CPPFLAGS)
# CXXFLAGS := -Wall -g -fmax-errors=3 $(CPPFLAGS)
LDFLAGS :=
LDLIBS :=

BLD := .build
EXT := .cc

VPATH = $(BLD)

.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean)))

ROOT_CXXFLAGS := $(shell root-config --cflags | sed 's/ -std=c++[^ ]\+ / /')
ROOT_LDFLAGS  := $(shell root-config --ldflags)
ROOT_LDLIBS   := $(shell root-config --libs)

SRCS := $(shell find -L src -type f -name '*$(EXT)')
DEPS := $(patsubst src/%$(EXT),$(BLD)/%.d,$(SRCS))

GREP_EXES := grep -rl '^ *int \+main *(' src --include='*$(EXT)'
EXES := $(patsubst src%$(EXT),bin%, $(shell $(GREP_EXES)))

all: $(EXES)

-include $(DEPS)

# -------------------------------------------------------------------
bin/poly_coords: linalg.o
bin/test_fit: linalg.o wls.o
bin/unbinned: linalg.o
bin/toy_conventional: linalg.o wls.o

C_test_fit := -fopenmp $(ROOT_CXXFLAGS)
L_test_fit := -fopenmp $(ROOT_LDLIBS) -lMinuit
C_unbinned := -fopenmp $(ROOT_CXXFLAGS)
L_unbinned := -fopenmp $(ROOT_LDLIBS) -lgsl -lgslcblas
C_toy_conventional := -fopenmp $(ROOT_CXXFLAGS)
L_toy_conventional := -fopenmp $(ROOT_LDLIBS) -lMinuit
# -------------------------------------------------------------------

$(DEPS): $(BLD)/%.d: src/%$(EXT)
	@mkdir -pv $(dir $@)
	$(CXX) $(CPPFLAGS) $(C_$*) -MM -MT '$(@:.d=.o)' $< -MF $@

$(BLD)/%.o:
	@mkdir -pv $(dir $@)
	$(CXX) $(CXXFLAGS) $(C_$*) -c $(filter %$(EXT),$^) -o $@

bin/%: $(BLD)/%.o
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

endif

clean:
	@rm -rfv $(BLD) bin

