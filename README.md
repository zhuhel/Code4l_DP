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
mkdir build
cd build
make
```

### How to run
```bash
cd run
./test_run.sh input.list
```
* Check your results under /output
* Histograms under /data-hist_output
* Trees under /data-tree_output

