/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack



//
// includes
//

#include <QuickAna/AnaToolRetrieve.h>

#include <QuickAna/EventData.h>
#include <QuickAna/IEventObjects.h>
#include <QuickAna/InternalConfiguration.h>

//
// method implementations
//

namespace ana
{
  AnaToolRetrieve ::
  AnaToolRetrieve (const std::string& name)
    : AsgTool (name), AnaTool (name),
      m_type (OBJECT_JET)
  {}



  AnaToolRetrieve ::
  AnaToolRetrieve (const std::string& name, const std::string& containerName)
    : AsgTool (name), AnaTool (name), m_containerName (containerName)
  {}



  StatusCode AnaToolRetrieve ::
  setObjectType (ObjectType type, const std::string& workingPoint)
  {
    m_type = type;
    if (!workingPoint.empty())
    {
      ATH_MSG_ERROR ("can't have working point in retrieve step");
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }



  AnalysisStep AnaToolRetrieve ::
  step () const
  {
    return STEP_RETRIEVE;
  }



  unsigned AnaToolRetrieve ::
  inputTypes () const
  {
    return 0;
  }



  unsigned AnaToolRetrieve ::
  outputTypes () const
  {
    return 1 << m_type;
  }



  void AnaToolRetrieve ::
  fillEventDataSource (EventData& event) const
  {
    for (auto& sys : event.getSystematics())
    {
      event.setSourceName (m_type, sys, m_containerName);
    }
  }



  StatusCode AnaToolRetrieve ::
  execute (IEventObjects& objects)
  {
    ATH_CHECK (objects.addCopy (m_type, m_containerName));
    return StatusCode::SUCCESS;
  }



  StatusCode AnaToolRetrieve ::
  getInitialConfiguration (InternalConfiguration& conf)
  {
    ATH_CHECK (AnaTool::getInitialConfiguration (conf));
    conf.setInputName (m_type, m_containerName);
    return StatusCode::SUCCESS;
  }
}
