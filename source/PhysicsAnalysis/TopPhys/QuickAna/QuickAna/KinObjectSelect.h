/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef QUICK_ANA__KIN_OBJECT_SELECT_H
#define QUICK_ANA__KIN_OBJECT_SELECT_H

//        Copyright Iowa State University 2014.
//                  Author: Nils Krumnack
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// Please feel free to contact me (nils.erik.krumnack@cern.ch) for bug
// reports, feature suggestions, praise and complaints.


// This module still needs to be documented.  The interface provided
// in this module is intended for experts only.  The module is
// considered to be in the pre-alpha stage.



#include <QuickAna/Global.h>

#include <QuickAna/xAODInclude.h>

namespace ana
{
  /// TODO: needs documentation
  class KinObjectSelect
  {
    //
    // public interface
    //

    /// effects: construct a kinematic selector based on the given formula
    /// guarantee: strong
    /// failures: out of memory II
    /// failures: formula errors
  public:
    KinObjectSelect (const std::string& formula);


    /// effects: apply this selection to the given particle
    /// guarantee: strong
    /// failures: xAOD errors
  public:
    bool select (const xAOD::IParticle& particle) const;

    //
    // private interface
    //

    /// description: the function we use internally
  private:
    std::function<bool(const xAOD::IParticle&)> m_function;
  };
}

#endif
