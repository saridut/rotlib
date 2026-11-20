# Minimal makefile for Sphinx documentation

SHELL = /bin/bash

ALL_DOCS = ../../docs
PACKAGE_NAME = rotlib
MODULE_PATH = ../$(PACKAGE_NAME)
SOURCE_DIR = source
BUILD_DIR  = build
OUT_DIR = $(ALL_DOCS)/python_api

.PHONY: doc
doc: $(ALL_DOCS)/rotconv.md
	if [ ! -d $(OUT_DIR) ]; then mkdir -p $(OUT_DIR); fi
	cp $^ $(SOURCE_DIR)
	sphinx-build -b html "$(SOURCE_DIR)" "$(OUT_DIR)"

.PHONY: rst
rst:
	sphinx-apidoc -M -e -f -d 2 -o "$(SOURCE_DIR)" "$(MODULE_PATH)"

.PHONY: clean
clean :
	rm -rf "$(OUT_DIR)"/*

.PHONY: clobber
clobber :
	rm -rf "$(SOURCE_DIR)"/$(PACKAGE_NAME)*.rst \
		"$(SOURCE_DIR)"/modules.rst \
		"$(BUILD_DIR)"/* \
		"$(OUT_DIR)"/*
