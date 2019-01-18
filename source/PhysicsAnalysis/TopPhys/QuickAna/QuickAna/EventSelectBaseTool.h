/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef QUICK_ANA__EVENT_SELECT_BASE_TOOL_H
#define QUICK_ANA__EVENT_SELECT_BASE_TOOL_H

//        Copyright Iowa State University 2016.
//                  Author: Nils Krumnack
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// Please feel free to contact me (nils.erik.krumnack@cern.ch) for bug
// reports, feature suggestions, praise and complaints.


#include <QuickAna/Global.h>

#include <QuickAna/AnaTool.h>
#include <QuickAna/SelectionData.h>

namespace ana
{

  /// @brief An AnaTool used to initialize the EventSelect object
  ///
  class EventSelectBaseTool : virtual public AnaTool
  {
    //
    // public interface
    //

    ASG_TOOL_CLASS (EventSelectBaseTool, ana::IAnaTool)

  public:

    /// effects: standard constructor
    /// guarantee: strong
    /// failures: out of memory II
    EventSelectBaseTool (const std::string& name);

    /// \copydoc IAnaTool::setObjectType
    virtual StatusCode
    setObjectType (ObjectType type, const std::string& workingPoint) override;

    /// returns: the processing step we are working on
    /// guarantee: no-fail
    virtual AnalysisStep step () const override;

    /// \copydoc IAnaTool::inputTypes
    virtual unsigned inputTypes () const override;

    /// \copydoc IAnaTool::outputTypes
    virtual unsigned outputTypes () const override;

  public:
    virtual StatusCode
    useConfiguration (const InternalConfiguration& configuration)
      override;

    /// effects: apply this analysis tool to the objects
    /// guarantee: basic
    /// failures: tool dependent
    virtual StatusCode execute (IEventObjects& objects) override;

    /// \brief the selection fields
  private:
    std::vector<SG::AuxElement::Accessor<SelectType>> m_select;
  };
}

#endif
