/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack



//
// includes
//

#include <QuickAna/EventObjects.h>

#include <QuickAna/MessageCheck.h>
#include <QuickAna/xAODHelpers.h>
#include <RootCoreUtils/Assert.h>

//
// method implementations
//

namespace ana
{
  EventObjects ::
  EventObjects ()
  {
    clear (nullptr, "");
  }


  void EventObjects ::
  clear (StoreType *store, const std::string& postfix)
  {
    m_postfix = postfix;
    m_store = store;
    for (unsigned iter = 0; iter < ObjectTypeInfo::numTypes; ++ iter)
    {
      m_objects[iter] = nullptr;
    }
    m_eventWeight = 1;
  }



  StatusCode EventObjects ::
  addCopy (ObjectType type, const std::string& name)
  {
    using namespace msgStore;

    RCU_REQUIRE (type < ObjectTypeInfo::numTypes);
    if (m_objects[type])
    {
      ANA_MSG_ERROR ("trying to initialize " + ObjectTypeInfo::name[type] + " twice");
      return StatusCode::FAILURE;
    }
    std::string anaName = "ana_" + ObjectTypeInfo::name[type] + "_" + m_postfix;
    ANA_CHECK (xAODHelpers::makeShallow (*m_store, anaName, type, m_objects[type], name));
    return StatusCode::SUCCESS;
  }



  void *EventObjects ::
  getVoid (ObjectType type) const
  {
    RCU_REQUIRE (type < ObjectTypeInfo::numTypes);
    return m_objects[type];
  }



  StatusCode EventObjects ::
  addNew (ObjectType type)
  {
    using namespace msgStore;

    RCU_REQUIRE (type < ObjectTypeInfo::numTypes);
    if (m_objects[type])
    {
      ANA_MSG_ERROR ("trying to initialize " + ObjectTypeInfo::name[type] + " twice");
      return StatusCode::FAILURE;
    }
    std::string anaName = "ana_" + ObjectTypeInfo::name[type] + "_" + m_postfix;
    ANA_CHECK (xAODHelpers::makeNew (*m_store, anaName, type, m_objects[type]));
    return StatusCode::SUCCESS;
  }



  //
  // legacy code
  //



  float EventObjects ::
  eventWeight () const
  {
    return m_eventWeight;
  }
}

