/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack



#ifndef QUICK_ANA__VALIDATION_EL_H
#define QUICK_ANA__VALIDATION_EL_H


#include <QuickAna/Global.h>

#ifdef ROOTCORE

#include <EventLoop/Algorithm.h>
#include <QuickAna/Configuration.h>
#include <memory>

class TH1;

namespace ana
{
  class ValidationEL : public EL::Algorithm, public ana::Configuration
  {
    /// description: the quickana tool handle
  public:
    std::unique_ptr<IQuickAna> quickAna; //!

    /// description: the histograms we create and fill
  public:
    ValidationHists *hists = nullptr; //!

    /// \brief whether to write AthSummary.txt
  public:
    bool m_writeSummary = false;

    /// \brief the number of files read
  public:
    unsigned m_filesRead = 0; //!

    /// \brief the number of events read
  public:
    unsigned m_eventsRead = 0; //!



    // this is a standard constructor
    ValidationEL ();

    // this is a standard destructor
    ~ValidationEL ();

    // these are the functions inherited from Algorithm
    virtual EL::StatusCode setupJob (EL::Job& job);
    virtual EL::StatusCode fileExecute ();
    virtual EL::StatusCode histInitialize ();
    virtual EL::StatusCode changeInput (bool firstFile);
    virtual EL::StatusCode initialize ();
    virtual EL::StatusCode execute ();
    virtual EL::StatusCode postExecute ();
    virtual EL::StatusCode finalize ();
    virtual EL::StatusCode histFinalize ();

    // this is needed to distribute the algorithm to the workers
    ClassDef(ValidationEL, 1);
  };
}

#endif

#endif
