#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "SampleHandler/DiskListLocal.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoop/CondorDriver.h"
#include "EventLoop/LocalDriver.h"
#include "EventLoop/TorqueDriver.h"
#include "EventLoopGrid/PrunDriver.h"

#include <PATInterfaces/SystematicCode.h>
#include "PATInterfaces/CorrectionCode.h"
#include "PATInterfaces/SystematicsUtil.h"

#include "TChain.h"

#include "MyAnalysis/MyxAODAnalysis.h"

#include <iostream>     // std::cout
#include <fstream>

using namespace std;

bool cmdline(int argc, char** argv, map<string,string> &opts);
void usage();

int main (int argc, char **argv) {
  map<string, string> opts;
  if (!cmdline(argc, argv, opts)) return 0;

  string outDir_in = opts["outdir"];
  if (outDir_in == "") {
    cout << "Name of output directory required!" << endl;
    return 1;
  }
  //string outDir = string(getenv("WORK")) + "/workarea/outData/" + outDir_in; // Laser was here
  string outDir = outDir_in;

  xAOD::Init().ignore();

  SH::SampleHandler sh;

  if (opts["in"] == "") {
    cout << "Name of input file(s) required!" << endl;
    return 2;
  }

  vector<string> ins;

  //<Get input list>
  ifstream filelist;
  filelist.open(opts["in"]);
  if(!filelist.good()) {
    cout<<"ERROR: Cannot find the input filelist, now quit!"<<endl;
    return 2;
  }
  string file;
  while(!filelist.eof()) {
    getline(filelist,file);        
    if(file.size()==0) continue; //remove the blank lines
    cout<<"Add file \""<<file<<"\""<<endl;
    //if dcache, use dcap://dcap.aglt2.org as prefix
    if(file.find("/pnfs")!=string::npos)
      file =  "dcache:" + file;
    if(file.find("/eos")!=string::npos)
      file =  "root://eosatlas/" + file;
    ins.push_back(file);
  } 

  TChain chain("CollectionTree");
  for (auto &f : ins) {
    chain.Add(f.c_str());
  }
  sh.add(SH::makeFromTChain("xAOD", chain));

  sh.setMetaString("nc_tree", "CollectionTree");

  sh.print();

  EL::Job job;
  job.sampleHandler(sh);

  MyxAODAnalysis* alg = new MyxAODAnalysis("physics");
  job.algsAdd( alg );

  //job.options()->setDouble(EL::Job::optMaxEvents, 5000);

  EL::DirectDriver driver;
  driver.submit(job, outDir);

  return 0;
}

bool cmdline(int argc, char** argv, map<string,string> &opts) {
  opts.clear();

  // defaults
  opts["outdir"] = "";
  opts["in"] = "";
  opts["debug"] = "0";

  for (int i=1;i<argc;i++) {

    string opt=argv[i];

    if (opt=="--help") {usage(); return false;}

    if (0 != opt.find("-")) {
      cout<<"ERROR: options start with '-'!"<<endl;
      return false;
    }
    opt.erase(0,1);
    if (opts.find(opt) == opts.end()) {
      cout<<"ERROR: invalid option '"<<opt<<"'!"<<endl;
      return false;
    }
    string nxtopt=argv[i+1];
    if (0 == nxtopt.find("-") || i+1 >= argc) {
      cout<<"ERROR: option '"<<opt<<"' requires value!"<<endl;
      return false;
    }

    opts[opt] = nxtopt;
    i++;
  }

  return true;
}

void usage()
{
  cout<<"USAGE: run [-option value]\n\n"
    <<"options [default]:\n\n"
    <<"-outdir (required!)\n"
    <<"-in (required!)\n"
    <<"-debug <0/1> [0]\n"
    <<endl;

  return;
}
