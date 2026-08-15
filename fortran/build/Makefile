#Makefile for Taylor with GCC compilers

SHELL = /bin/bash

DEBUG = true#{true, false}  (NOTE: no whitespace at end)

FC = gfortran

#*************************************************

ifeq ($(DEBUG), true)
FCFLAGS = -c -O0 -Wall -fcheck=all -fbacktrace -std=f2008ts \
		  -fdefault-double-8 \
		  -fdefault-real-8 \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
else
FCFLAGS = -c -O3 -std=f2008ts \
		  -fdefault-double-8 \
		  -fdefault-real-8 \
		  -funroll-loops \
		  -flto          \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
endif

SRC_DIR = ../src
BUILD_DIR = .
TEST_DIR = ../test
EXECUTABLE = $(TEST_DIR)/test_rotlib

vpath %.f90 $(SRC_DIR) $(TEST_DIR)
#vpath $(EXECUTABLE) $(INSTALL_DIR)

INCS = -I$(SRC_DIR) -I$(TEST_DIR)

OBJ_NAMES = \
		  m_precision        \
		  m_constants_math   \
		  m_utils_math       \
		  m_rotlib           \
		  test_main

OBJECTS = $(addsuffix .o,$(OBJ_NAMES))

#move : $(EXECUTABLE)
#	if [ -e "$(EXECUTABLE)" ]; then mv $(EXECUTABLE) $(INSTALL_DIR);fi 
#
$(EXECUTABLE) : $(OBJECTS)
	$(FC) $(INCS) $^ -o $@

$(OBJECTS) : %.o : %.f90
	$(FC) $(INCS) $(FCFLAGS) $^ -o $@

docs :
	ford ../rotlib_f.md 

.PHONY : clean clobber
clean:
	-rm -rf *.o *.mod
clobber:
	-rm -rf $(EXECUTABLE) *.o *.mod
