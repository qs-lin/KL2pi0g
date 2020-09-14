# -*- Makefile -*-

###
# Package Name
#  in usual case, the following description works:
#    PACKAGE = $(shell basename `pwd`)

PACKAGE      = $(shell basename `pwd`)

###
# Name of Executables
#  list up names of executables which will be compiled

MAINS        = g5anaKL2pi0g 
MAINS        += read

###
# INSTALL_LIBS
#  set 'no' if you do NOT want to copy produced shared object ".so"
#  under $MY_TOP_DIR/lib of your environment
#  in usual case, set 'yes' if the package is located under $MY_TOP_DIR/AnalysisLibrary

INSTALL_LIBS =



## optional flags ##

###
# Extra Libraries
#  if the package depends on external libraries, you need to add those libraries explicitly
#    ex) EXTRALIBS += -L/use/lib/foo -lFoo

#EXTRALIBS     += -L/home/had/jaylin/workspace/share/CLHEP/2.1.1.0/lib -lCLHEP-2.1.1.0

###
# Compilation Flags for ROOT cint
#  you can give extra options for the compilation of ROOT dictionaries with ROOT cint

ROOTCINTFLAGS +=

###
# LD_LINK_LIBS
#  you can associate other libraries to the shared object (.so) which will be generated in this Makefile
#  so that one can resolve all dependency to the associated libraries by the reference only to the shared object

LD_LINK_LIBS  += $(EXTRALIBS)

###
# INSTALL_BINS
#  DO NOT change for E14 Analysis Suite

INSTALL_BINS = no

include $(MY_TOP_DIR)/config/Makefile.pre

#include $(E14_CONFIG_DIR)/Makefile.g4
include $(E14_CONFIG_DIR)/Makefile.bin
