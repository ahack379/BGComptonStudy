#!/bin/bash

temp=make_rootmap.$$.temp

rootlibmap() {
 ROOTMAP=$1
 SOFILE=$2
 LINKDEF=$3
 shift 3
 DEPS=$*
 if [[ -e $SOFILE && -e $LINKDEF ]]; then
     rlibmap -f -o $ROOTMAP -l $SOFILE -d $DEPS -c $LINKDEF 2>> $temp
     rm -f $temp
 fi
}

######################################################
# Package_Name
rootlibmap libBGComptonStudy_Package_Name.rootmap libBGComptonStudy_Package_Name.so $LARLITE_USERDEVDIR/BGComptonStudy/Package_Name/LinkDef.h \
    libLArLite_Analysis.so













