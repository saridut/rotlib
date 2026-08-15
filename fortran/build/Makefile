#Makefile for GCC compilers

SHELL = /bin/bash

DEBUG = true#{true, false}  (NOTE: no whitespace at end)

FC = gfortran

#*************************************************

ifeq ($(DEBUG), true)
FCFLAGS = -c -O0 -Wall -fcheck=all -fbacktrace -std=f2008ts \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
else
FCFLAGS = -c -O3 -std=f2008ts \
		  -funroll-loops \
		  -flto          \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
endif

SRC_DIR = ../src
BUILD_DIR = .
LIB_DIR = ../lib
TEST_DIR = ../test
TEST_EXECUTABLE = $(TEST_DIR)/test_rotlib
LIB_NAME = $(LIB_DIR)/librotlib.a

vpath %.f90 $(SRC_DIR) $(TEST_DIR)

INCS = -I$(SRC_DIR) -I$(TEST_DIR)

SRC_OBJ_NAMES = \
		  m_precision        \
		  m_constants_math   \
		  m_utils_math       \
		  m_rotlib

TEST_OBJ_NAMES = test_main

SRC_OBJECTS = $(addsuffix .o,$(SRC_OBJ_NAMES))
TEST_OBJECTS = $(addsuffix .o,$(TEST_OBJ_NAMES))

.PHONY : test
test : $(TEST_EXECUTABLE)

.PHONY : lib
lib : $(LIB_NAME)

.PHONY : doc
doc :
	ford ../rotlib_f.md 

$(TEST_EXECUTABLE) : $(SRC_OBJECTS) $(TEST_OBJECTS)
	$(FC) $(INCS) $^ -o $@

$(SRC_OBJECTS) : %.o : %.f90
	$(FC) $(INCS) $(FCFLAGS) $^ -o $@

$(TEST_OBJECTS) : %.o : %.f90
	$(FC) $(INCS) $(FCFLAGS) $^ -o $@

$(LIB_NAME) : $(SRC_OBJECTS)
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	ar -rcs $@ $^

.PHONY : clean clobber
clean:
	-rm -rf *.o *.mod
clobber:
	-rm -rf $(TEST_EXECUTABLE) $(LIB_NAME) *.o *.mod
