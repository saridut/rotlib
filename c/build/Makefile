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
LIB_DIR = ../lib
TEST_DIR = ../test
DOC_DIR = ../../docs/c_api
TEST_EXECUTABLE = $(TEST_DIR)/test_rotlib
LIB_NAME = $(LIB_DIR)/librotlib.a

vpath %.c $(SRC_DIR) $(TEST_DIR)

INCS = -I$(SRC_DIR) -I$(TEST_DIR)

SRC_OBJ_NAMES = utils_math rotlib
TEST_OBJ_NAMES = test_main

SRC_OBJECTS = $(addsuffix .o,$(SRC_OBJ_NAMES))
TEST_OBJECTS = $(addsuffix .o,$(TEST_OBJ_NAMES))

.PHONY : test
test : $(TEST_EXECUTABLE)

.PHONY : lib
lib : $(LIB_NAME)

.PHONY : doc
doc :
	if [ ! -d $(DOC_DIR) ]; then mkdir -p $(DOC_DIR); fi
	doxygen ../Doxyfile 

$(TEST_EXECUTABLE) : $(SRC_OBJECTS) $(TEST_OBJECTS)
	$(CC) $(INCS) $^ -o $@

$(SRC_OBJECTS) : %.o : %.c
	$(CC) $(INCS) $(CFLAGS) $^ -o $@

$(TEST_OBJECTS) : %.o : %.c
	$(CC) $(INCS) $(CFLAGS) $^ -o $@

$(LIB_NAME) : $(SRC_OBJECTS)
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	ar -rcs $@ $^

.PHONY : clean clobber
clean:
	-rm -rf *.o
clobber:
	-rm -rf $(TEST_EXECUTABLE) $(LIB_NAME) *.o
