/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack



#ifdef ROOTCORE

//
// includes
//

#include <QuickAna/ValidationEL.h>

#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include <QuickAna/MessageCheck.h>
#include <QuickAna/QuickAna.h>
#include <QuickAna/ValidationHists.h>
#include <fstream>

//
// method implementations
//

// this is needed to distribute the algorithm to the workers
ClassImp(ana::ValidationEL)

namespace ana
{
  ValidationEL :: ValidationEL ()
  {
    // Here you put any code for the base initialization of variables,
    // e.g. initialize all pointers to 0.  Note that you should only put
    // the most basic initialization here, since this method will be
    // called on both the submission and the worker node.  Most of your
    // initialization code will go into histInitialize() and
    // initialize().
  }



  ValidationEL :: ~ValidationEL ()
  {
    delete hists;
  }



  EL::StatusCode ValidationEL :: setupJob (EL::Job& job)
  {
    // Here you put code that sets up the job on the submission object
    // so that it is ready to work with your algorithm, e.g. you can
    // request the D3PDReader service or add output files.  Any code you
    // put here could instead also go into the submission script.  The
    // sole advantage of putting it here is that it gets automatically
    // activated/deactivated when you add/remove the algorithm from your
    // job, which may or may not be of value to you.
    job.useXAOD();
    return EL::StatusCode::SUCCESS;
  }



  EL::StatusCode ValidationEL :: histInitialize ()
  {
    // Here you do everything that needs to be done at the very
    // beginning on each worker node, e.g. create histograms and output
    // trees.  This method gets called before any input files are
    // connected.

    return EL::StatusCode::SUCCESS;
  }



  EL::StatusCode ValidationEL :: fileExecute ()
  {
    // Here you do everything that needs to be done exactly once for every
    // single file, e.g. collect a list of all lumi-blocks processed
    ++ m_filesRead;
    return EL::StatusCode::SUCCESS;
  }



  EL::StatusCode ValidationEL :: changeInput (bool /*firstFile*/)
  {
    // Here you do everything you need to do when we change input files,
    // e.g. resetting branch addresses on trees.  If you are using
    // D3PDReader or a similar service this method is not needed.
    return EL::StatusCode::SUCCESS;
  }



  EL::StatusCode ValidationEL :: initialize ()
  {
    // Here you do everything that you need to do after the first input
    // file has been connected and before the first event is processed,
    // e.g. create additional histograms based on which variables are
    // available in the input files.  You can also create all of your
    // histograms and trees in here, but be aware that this method
    // doesn't get called if no events are processed.  So any objects
    // you create here won't be available in the output if you have no
    // input events.
    ANA_CHECK_SET_TYPE (EL::StatusCode);

    std::unique_ptr<QuickAna> myQuickAna (new QuickAna ("quickana"));
    myQuickAna->setConfig (*this);
    quickAna = std::move (myQuickAna);
    ANA_CHECK (quickAna->initialize());
    hists = new ValidationHists (wk(), quickAna->internalConfiguration());

    for (auto sys : quickAna->affectingSystematics())
      ATH_MSG_INFO ("found systematic: \"" << sys.name() << "\"");
    return EL::StatusCode::SUCCESS;
  }



  EL::StatusCode ValidationEL :: execute ()
  {
    // Here you do everything that needs to be done on every single
    // events, e.g. read input variables, apply cuts, and fill
    // histograms and trees.  This is where most of your actual analysis
    // code will go.

    ANA_CHECK_SET_TYPE (EL::StatusCode);

    ++ m_eventsRead;

    ANA_CHECK (quickAna->process (*wk()->xaodEvent()));

    hists->fillHists ("val/", quickAna->objects(), quickAna->weight());
    return EL::StatusCode::SUCCESS;
  }



  EL::StatusCode ValidationEL :: postExecute ()
  {
    // Here you do everything that needs to be done after the main event
    // processing.  This is typically very rare, particularly in user
    // code.  It is mainly used in implementing the NTupleSvc.
    return EL::StatusCode::SUCCESS;
  }



  EL::StatusCode ValidationEL :: finalize ()
  {
    // This method is the mirror image of initialize(), meaning it gets
    // called after the last event has been processed on the worker node
    // and allows you to finish up any objects you created in
    // initialize() before they are written to disk.  This is actually
    // fairly rare, since this happens separately for each worker node.
    // Most of the time you want to do your post-processing on the
    // submission node after all your histogram outputs have been
    // merged.  This is different from histFinalize() in that it only
    // gets called on worker nodes that processed input events.
    return EL::StatusCode::SUCCESS;
  }



  EL::StatusCode ValidationEL :: histFinalize ()
  {
    // This method is the mirror image of histInitialize(), meaning it
    // gets called after the last event has been processed on the worker
    // node and allows you to finish up any objects you created in
    // histInitialize() before they are written to disk.  This is
    // actually fairly rare, since this happens separately for each
    // worker node.  Most of the time you want to do your
    // post-processing on the submission node after all your histogram
    // outputs have been merged.  This is different from finalize() in
    // that it gets called on all worker nodes regardless of whether
    // they processed input events.

    if (m_writeSummary)
    {
      std::ofstream stream ("AthSummary.txt");
      stream << "Files read: " << m_filesRead << std::endl;
      stream << "Events Read: " << m_eventsRead << std::endl;
      stream << "{ \"events\":{\"read\":" << m_eventsRead << "} " << "}" <<
	std::endl;
    }

    return EL::StatusCode::SUCCESS;
  }
}
#endif
