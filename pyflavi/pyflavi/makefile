# determine the platform
UNAME_S := $(shell uname -s)

# CXX and CXXFLAGS                                                                                                                                                         
ifeq ($(UNAME_S),Darwin)
  CXX := clang++
  CXXFLAGS:=-Ofast -Wall -shared -std=c++11 -undefined dynamic_lookup -fPIC
else
  CXX := g++
  CXXFLAGS:=-Ofast -Wall -shared -std=c++11 -undefined -fPIC
endif

# pybind11 include paths
INCPYBIND := $(shell python3 -m pybind11 --includes)

# all include paths
INC:=$(INCPYBIND) -I/src -I${CURDIR}/src -I${CURDIR}/src/model -I${CURDIR}/../src/utils

# use output of python3-config --extension-suffix
EXT:=.cpython-38-darwin.so

DEPS:=src/pyflavi.cpp ../src/utils/utils.cpp ../src/utils/WeibullAgeSampler.cpp ../src/utils/event_vec.cpp ../src/model/individual.cpp src/bind_pyflavi.cpp

lib/pyflavi.$(EXT): $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) $(DEPS) -o lib/pyflavi$(EXT)

#bin/main: ../src/utils/utils.cpp ../src/utils/WeibullAgeSampler.cpp ../src/utils/event_vec.cpp ../#src/model/individual.cpp src/pyflavi.cpp src/main.cpp
#		$(CXX) -O3 -Wall -std=c++11 $(INC) ../src/utils/utils.cpp ../src/utils/#WeibullAgeSampler.cpp ../src/utils/event_vec.cpp ../src/model/individual.cpp src/pyflavi.cpp src/#main.cpp -o bin/main

clean:
	@rm lib/*$(EXT)

.PHONY=clean
