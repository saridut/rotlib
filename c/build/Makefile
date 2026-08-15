#Makefile for Taylor with GCC compilers

SHELL = /bin/bash

DEBUG = true#{true, false}  (NOTE: no whitespace at end)

CC = gcc

#*************************************************

ifeq ($(DEBUG), true)
CFLAGS = -c -O0 -Wall -std=c23 \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
else
CFLAGS = -c -O3 -std=c23 \
		  -funroll-loops \
		  -flto          \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
endif

SRC_DIR = ../src
BUILD_DIR = .
TEST_DIR = ../test
EXECUTABLE = $(TEST_DIR)/test_rotlib

vpath %.c $(SRC_DIR) $(TEST_DIR)

INCS = -I$(SRC_DIR) -I$(TEST_DIR)

OBJ_NAMES = \
		  utils_math \
		  rotlib     \
		  test_main

OBJECTS = $(addsuffix .o,$(OBJ_NAMES))

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(INCS) $^ -o $@

$(OBJECTS) : %.o : %.c
	$(CC) $(INCS) $(CFLAGS) $^ -o $@

docs :
	doxygen ../Doxyfile 

.PHONY : clean clobber
clean:
	-rm -rf *.o *.mod
clobber:
	-rm -rf $(EXECUTABLE) *.o *.mod
