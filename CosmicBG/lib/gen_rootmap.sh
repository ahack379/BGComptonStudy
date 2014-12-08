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
# CosmicBG
rootlibmap libBGComptonStudy_CosmicBG.rootmap libBGComptonStudy_CosmicBG.so $LARLITE_USERDEVDIR/BGComptonStudy/CosmicBG/LinkDef.h \
 	libMCPartInfo_MCPartGetter.so libLArLite_LArUtil.so libBasicTool_GeoAlgo.so libLArLite_Analysis.so 
	
	
	
	














