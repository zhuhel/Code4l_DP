The is a framework used for (at least) 4 leptons selection for Dark Photon analysis. 
For more details of selections, one can find in the paper of https://arxiv.org/pdf/1710.07635.pdf

Original code from Bing Li (bing.li@cern.ch)

NOTE: Code is under development. Do not change releaser, etc.
Otherwise could have quite some technical WARNING/ERROR.

Rel. 21.2.8, CMake based
General info, refer to
https://twiki.cern.ch/twiki/bin/view/AtlasComputing/RootCoreToCMake

### Structure:
* ./source: analysis code
* ./build: compilation
* ./run: run folder

### First time setup:
```bash
source setup.sh
```

### Compilation
```bash
cd build
make
# Sometimes the 'make' may be failed, just *repeat make* it again (a known issue for cmake packages?).
```

### How to run
```bash
cd run
./test_run.sh input.list
```
* Check the cutflow under /cutflow
* Check your results under /output
* Histograms under /data-hist_output
* Trees under /data-tree_output

### Grid submission
```bash
localSetupRucioClients
lsetup panda
grid_submission datasetName
```
* The codes of grid submission are under ./source/MyAnalysis/util/
* Remember to change the *User* inside the code to your own
