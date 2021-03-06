/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack



//
// includes
//

#include <QuickAna/EventInfoIsData.h>

#include <RootCoreUtils/Assert.h>
#include <xAODEventInfo/EventInfo.h>

//
// method implementations
//

namespace ana
{
  void EventInfoIsData ::
  testInvariant () const
  {
    RCU_INVARIANT (this != nullptr);
  }



  EventInfoIsData ::
  EventInfoIsData (const std::string& name, bool isData)
    : AsgTool (name), AnaTool (name), m_isData (isData)
  {
    RCU_NEW_INVARIANT (this);
  }



  StatusCode EventInfoIsData ::
  execute (IEventObjects& /*objects*/)
  {
    RCU_CHANGE_INVARIANT (this);

    if (!m_checkedInput)
    {
      const xAOD::EventInfo *eventInfo = nullptr;
      ATH_CHECK (evtStore()->retrieve (eventInfo, "EventInfo"));
      bool isSim = eventInfo->eventType (xAOD::EventInfo::IS_SIMULATION);
      if (isSim == m_isData)
      {
	ATH_MSG_FATAL ("isDataFlag = " << m_isData << " is inconsistent with EventInfo isSim = " << isSim);
	return StatusCode::FAILURE;
      }
      m_checkedInput = true;
    }
    return StatusCode::SUCCESS;
  }



  StatusCode EventInfoIsData ::
  setObjectType (ObjectType /*type*/, const std::string& /*workingPoint*/)
  {
    return StatusCode::SUCCESS;
  }



  AnalysisStep EventInfoIsData ::
  step () const
  {
    return STEP_RETRIEVE;
  }



  unsigned EventInfoIsData ::
  inputTypes () const
  {
    return 0;
  }



  unsigned EventInfoIsData ::
  outputTypes () const
  {
    return 0;
  }
}
