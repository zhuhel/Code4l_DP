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
   //std::string submitDir = "output_364250";
   std::string User = "hezhu";
   std::string dataset = "mc16_13TeV.364250.Sherpa_222_NNPDF30NNLO_llll.deriv.DAOD_STDM3.e5894_s3126_r9364_r9315_p3371/";
   if( argc > 1 ) dataset = argv[ 1 ];
   std::string submitDir = "output_"+dataset;

   // set the output name with user and date
   time_t now = time(0);
   tm *ltm = localtime(&now);
   char* datestr = Form("%d%d%d", 1900+ltm->tm_year, 1+ltm->tm_mon, ltm->tm_mday);
   //std::string outstr = Form("user.hezhu.%s.%s", datestr, dataset.c_str());
   // too long, need to split
   std::stringstream ss(dataset);
   std::string item;
   char delim = '.';
   vector<string> tokens;
   while (std::getline(ss, item, delim)) {
     //cout << item << endl;
     tokens.push_back(item);
   }
   string mdata = "";
   mdata = tokens.at(0)+"."+tokens.at(1);
   std::string outstr = Form("user.%s.%s.%s", User.c_str(), datestr, mdata.c_str());


   // Set up the job for xAOD access:
   xAOD::Init().ignore();
   
   // Construct the samples to run on:
   SH::SampleHandler sh;
   SH::scanRucio (sh, dataset);

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
   driver.options()->setString("nc_outputSampleName", outstr);
   cout << "Submitting job with name: " << outstr << endl;
   //driver.options()->setString(EL::Job::optGridNFilesPerJob, "MAX"); //By default, split in as few jobs as possible
   driver.options()->setDouble("nc_nFilesPerJob", 10); 
   driver.options()->setDouble("nc_mergeOutput", 0);  
   //sh.get("mc14_13TeV.110401.PowhegPythia_P2012_ttbar_nonallhad.merge.DAOD_STDM4.e2928_s1982_s2008_r5787_r5853_p1807/")->SetMetaDouble(EL::Job::optGridNFilesPerJob, 1); //For this particular sample, split into one job per input file
   //driver.options()->setDouble(EL::Job::optGridMergeOutput, 1); //run merging jobs for all samples before downloading (recommended) 
   //driver.submit( job, submitDir );
   driver.submitOnly (job, submitDir);
   
   
   return 0;
}
