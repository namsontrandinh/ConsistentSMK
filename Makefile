# =====================================================
# Compiler configuration
# =====================================================

CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Isrc

# OpenMP: required for REVENUE and IC
OMPFLAGS = -fopenmp

# Libraries
LDLIBS   = -lpthread


# =====================================================
# Targets
# =====================================================

.PHONY: all clean debug

all: maxcut revenue ic preproc preproc_ic


# =====================================================
# MAX-CUT objective
# =====================================================

maxcut: src/main.cpp src/sfunctions.h src/sfunctions_impl.h
	$(CXX) src/main.cpp -o maxcut \
	    $(CXXFLAGS) \
	    -DSFUNC_MAXCUT \
	    $(LDLIBS)


# =====================================================
# REVENUE objective
# OpenMP enabled
# =====================================================

revenue: src/main.cpp src/sfunctions.h src/sfunctions_impl.h
	$(CXX) src/main.cpp -o revenue \
	    $(CXXFLAGS) \
	    $(OMPFLAGS) \
	    -DSFUNC_REVENUE \
	    $(LDLIBS)


# =====================================================
# INFLUENCE MAXIMIZATION - IC objective
# OpenMP enabled
# =====================================================

ic: src/main.cpp src/sfunctions.h src/sfunctions_impl.h
	$(CXX) src/main.cpp -o ic \
	    $(CXXFLAGS) \
	    $(OMPFLAGS) \
	    -DSFUNC_IC \
	    $(LDLIBS)


# =====================================================
# Preprocess general graph
# edges.txt -> graph.bin
# =====================================================

preproc: src/data/preprocess.cpp
	$(CXX) src/data/preprocess.cpp -o preproc \
	    $(CXXFLAGS) \
	    $(LDLIBS)


# =====================================================
# Preprocess IC graph
# Always directed + normalize incoming probabilities
# =====================================================

preproc_ic: src/data/preprocess_ic.cpp
	$(CXX) src/data/preprocess_ic.cpp -o preproc_ic \
	    $(CXXFLAGS) \
	    $(LDLIBS)


# =====================================================
# Debug MAX-CUT
# =====================================================

debug: src/main.cpp src/sfunctions.h src/sfunctions_impl.h
	$(CXX) src/main.cpp -o maxcut_debug \
	    -std=c++17 \
	    -O0 \
	    -g \
	    -ggdb3 \
	    -Wall \
	    -Isrc \
	    -DSFUNC_MAXCUT \
	    $(LDLIBS)


# =====================================================
# Clean
# =====================================================

clean:
	rm -f maxcut revenue ic preproc preproc_ic maxcut_debug