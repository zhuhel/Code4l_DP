#!/bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
setupATLAS

mkdir build
cd build
rm -rf * 
acmSetup --sourcedir=../source AnalysisBase,21.2.64
cp ../setup_bld.sh .
#asetup AnalysisBase,21.2.26
#source /afs/cern.ch/work/h/hezhu/public/workplace/DarkPh/code4l_dp/build/x86_64-slc7-gcc62-opt/setup.sh

