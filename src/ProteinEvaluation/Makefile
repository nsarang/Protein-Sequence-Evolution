TARGETS  := gscores gprofs driver
CC       := gcc-8
CXX      := g++-8
# INCLUDE  := .
# SRC      := .
OBJ      := .
# SRCS    := $(wildcard $(SRC)/*.c)
SRCS     := protein.cpp protein_profile.cpp evaluator.cpp utility.cpp mapping.cpp dDFIRE2.cpp statistics.cpp deep.cpp
OBJS     := $(subst .cpp,.o,$(SRCS))
# EXE      := $(BIN)/main
CXXFLAGS += -std=c++17 -g -Wall -pedantic
# CFLAGS   := -I$(INCLUDE)
LDFLAGS  += -lstdc++fs -pthread
# LDLIBS   := -lm
MKDIR    := md
RMDIR    := rd /S /Q
RM       := rm -f


ifeq (, $(shell which $(CXX)))
$(error "$(CXX) is required")
endif


.PHONY: all clean

all: $(TARGETS)

$(TARGETS): % : %.o $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o "$(BINDIR)/$@" $(LDLIBS)

$(OBJ)/%.o: $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJ)/*.o