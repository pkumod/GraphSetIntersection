CXX = g++ -std=c++11
CXXFLAGS = -O3 -Wall -Wno-unused-result -mavx2 -lmetis -pthread

# SOURCES = $(shell find . | grep -e ".cpp" | grep -v "tc.cpp" | grep -v "reorder.cpp"| grep -v "mc.cpp")
# OBJECTS = $(SOURCES:%.cpp=%.o)
SOURCES = $(shell find . | grep -e ".hpp")
OBJECTS = $(SOURCES:%.hpp=%.o)
OBJECTS := $(OBJECTS) roaring.o

all: tc mc sm reorder test_intersection_algos

tc: $(OBJECTS) tc.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

mc: $(OBJECTS) mc.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

sm: $(OBJECTS) sm.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

reorder: $(OBJECTS) reorder.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

test_intersection_algos: $(OBJECTS) exp/test_intersection_algos.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

%.o: %.cpp %.hpp util.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

roaring.o: roaring/roaring.c roaring/roaring.h roaring/roaring.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf tc mc sm reorder test_intersection_algos $(OBJECTS)