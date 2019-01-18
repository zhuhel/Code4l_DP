#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "SampleHandler/DiskListLocal.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/PrunDriver.h"

#include <PATInterfaces/SystematicCode.h>
#include "PATInterfaces/CorrectionCode.h"
#include "PATInterfaces/SystematicsUtil.h"

#include "MyAnalysis/MyxAODAnalysis.h"

#include <iostream>     // std::cout
#include <fstream> 

using namespace std;

int main( int argc, char* argv[] ) {

   // Take the submit directory from the input if provided:
   std::string submitDir = "outputI";
   if( argc > 1 ) submitDir = argv[ 1 ];

   // Set up the job for xAOD access:
   xAOD::Init().ignore();
   
   // Construct the samples to run on:
   SH::SampleHandler sh;
   //SH::scanRucio (sh, "data17_13TeV.periodI.physics_Main.PhysCont.DAOD_HIGG2D1.grp17_v01_p3213/");
   SH::scanRucio (sh, "data17_13TeV.periodI.physics_Main.PhysCont.DAOD_STDM3.grp15_v11_p3372_p3388_p3402/");

   // Set the name of the input TTree. It's always "CollectionTree"
   // for xAOD files.
   sh.setMetaString( "nc_tree", "CollectionTree" );
   
   // Print what we found:
   sh.print();
   
   // Create an EventLoop job:
   EL::Job job;
   job.sampleHandler( sh );
   
   // Add our analysis to the job:
   MyxAODAnalysis* alg = new MyxAODAnalysis("physics");
 
   job.algsAdd( alg );
   
   // Run the job using the local/direct driver:
   //EL::DirectDriver driver;
   EL::PrunDriver driver;
   driver.options()->setString("nc_outputSampleName", "user.bili.Data17_periodI.v1p3");
   //driver.options()->setString(EL::Job::optGridNFilesPerJob, "MAX"); //By default, split in as few jobs as possible
   driver.options()->setDouble("nc_nFilesPerJob", 10); 
   driver.options()->setDouble("nc_mergeOutput", 0);  
   //sh.get("mc14_13TeV.110401.PowhegPythia_P2012_ttbar_nonallhad.merge.DAOD_STDM4.e2928_s1982_s2008_r5787_r5853_p1807/")->SetMetaDouble(EL::Job::optGridNFilesPerJob, 1); //For this particular sample, split into one job per input file
   //driver.options()->setDouble(EL::Job::optGridMergeOutput, 1); //run merging jobs for all samples before downloading (recommended) 
   //driver.submit( job, submitDir );
   driver.submitOnly (job, submitDir);
   
   
   return 0;
}
