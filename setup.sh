#!/bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
setupATLAS

mkdir build
cd build
rm -rf * 
#acmSetup --sourcedir=../source AnalysisBase,21.2.8
acmSetup --sourcedir=../source AnalysisBase,21.2.10
## will auto cmake 


### acm sparse_clone_project athena

#acm compile

#cd ..
