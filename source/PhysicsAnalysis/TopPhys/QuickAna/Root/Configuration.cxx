/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//        Copyright Iowa State University 2014.
//                  Author: Nils Krumnack
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// Please feel free to contact me (nils.erik.krumnack@cern.ch) for bug
// reports, feature suggestions, praise and complaints.


//
// includes
//

#include <QuickAna/Configuration.h>

//
// method implementations
//

namespace ana
{
  Configuration ::
  Configuration ()
    : eventinfoDef ("none"), eventSelectDef ("none"),
      electronDef ("none"), photonDef ("none"), muonDef ("none"),
      tauDef ("none"), jetDef ("none"),truthjetDef("none"), pfJetDef ("none"), fatJetDef ("none"),
      metDef ("none"), met2Def ("none"), cleanDef ("none"), orDef ("none"), triggerDef(""),
      schedulerDef ("basic")
  {}



  void Configuration ::
  setConfig (const Configuration& conf)
  {
    *this = conf;
  }
}
