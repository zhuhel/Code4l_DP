Original code from Bing Li (bing.li@cern.ch)
NOTE: Code is under development. Do not change releaser, etc.
Otherwise could have quite some technical WARNING/ERROR.

Rel. 21.2.8, CMake based
General info, refer to
https://twiki.cern.ch/twiki/bin/view/AtlasComputing/RootCoreToCMake

###############################################
Structure:
==> ./source: analysis code
==> ./build: compilation
==> ./run: run folder

###############################################
First time setup:
==> source setup.sh

Compilation
==> cd build
==> make

###############################################
How to run
==> cd run
==> ./test_run.sh input.list
==> Please refer to input_mc16_364250.list to see the input format
==> Check your results under /output
==> Histograms under data-hist_output
==> Trees under data-tree_output

###############################################
#Grid submission
#==> source ~/setups/rucio.sh
#==> lsetup panda
