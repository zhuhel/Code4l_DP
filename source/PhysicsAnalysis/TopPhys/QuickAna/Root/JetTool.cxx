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

#include <QuickAna/JetTool.h>

#include <JetCalibTools/JetCalibrationTool.h>
#include <JetResolution/JERSmearingTool.h>
#include <JetResolution/JERTool.h>
#include <JetSelectorTools/JetCleaningTool.h>
#include <JetUncertainties/JetUncertaintiesTool.h>
#include <QuickAna/AnaToolRetrieve.h>
#include <QuickAna/DefinitionArgs.h>
#include <QuickAna/DefinitionMaker.h>
#include <QuickAna/InternalConfiguration.h>
#include <QuickAna/MessageCheck.h>
#include <xAODBTaggingEfficiency/BTaggingEfficiencyTool.h>
#include <xAODBTaggingEfficiency/BTaggingSelectionTool.h>
#include <JetMomentTools/JetVertexTaggerTool.h>
#include <JetJvtEfficiency/JetJvtEfficiency.h>
#include <AthContainers/ConstDataVector.h>

static const float GeV = 1000.;
static const float TeV = 1e6;

// BTag config options currently depend on release
namespace
{
#if ROOTCORE_RELEASE_SERIES <= 23
  const char* btagAlgDefault = "MV2c20";
  const std::string bTagCalibFile =
    "xAODBTaggingEfficiency/13TeV/2016-Winter-13TeV-MC15-CDI-March10_v1.root";
  const char *jesFile = "JES_data2016_data2015_Recommendation_Dec2016.config";
  const std::string uncertConfigFile = "JES_2015/Moriond2016/JES2015_SR_Scenario1.config";
  const char *mcType = "MC15";
#elif ROOTCORE_RELEASE_SERIES == 24
  const char* btagAlgDefault = "MV2c10";
  const std::string bTagCalibFile =
    "xAODBTaggingEfficiency/13TeV/2016-20_7-13TeV-MC15-CDI-2017-04-24_v1.root";
  const char *jesFile = "JES_data2016_data2015_Recommendation_Dec2016.config";
  const std::string uncertConfigFile = "JES_2016/Moriond2017/JES2016_SR_Scenario1.config";
  const char *mcType = "MC15";
#else
  const char* btagAlgDefault = "MV2c10";
  const std::string bTagCalibFile =
    "xAODBTaggingEfficiency/13TeV/2017-21-13TeV-MC16-CDI-2018-06-29_v1.root";
  const char *jesFile = "JES_data2017_2016_2015_Recommendation_Aug2018_rel21.config";
  const char *jesFile_AFII = "JES_MC16Recommendation_AFII_EMTopo_April2018_rel21.config";
  const std::string uncertConfigFile = "rel21/Moriond2018/R4_StrongReduction_Scenario1.config";
  const char *mcType = "MC16";
#endif
}

//
// method implementations
//

namespace ana
{
  JetToolCorrect ::
  JetToolCorrect (const std::string& name)
    : AsgTool (name), AnaToolCorrect<xAOD::JetContainer> (name),
      m_calibration_tool ("calibration", this),
      m_uncertainties_tool ("uncertainties", this),
      m_resolution_tool ("resolution", this),
      m_smearing_tool ("smearing", this),
      m_jvt_tool ("jvt", this),
      m_jvtEffTool("jvt_eff", this),
    m_fjvtEffTool("fjvt_eff", this),
      m_bsel_tool ("btag", this),
      m_bsel_OR_tool ("btag_OR", this),
      m_cleaning_tool ("cleaning", this)
  {
    declareProperty("EnableBTagging", m_enableBTagging = true);
    declareProperty("BTagger", m_btagger = btagAlgDefault);
    declareProperty("BTagWP", m_btagWP = "FixedCutBEff_77");
    // Set to a large negative number to disable
    // TODO: add a better way to disable
    declareProperty("BTagWP_OR", m_btagWP_OR = "FixedCutBEff_85");
  }



  StatusCode JetToolCorrect ::
  useInitialConfiguration (const InternalConfiguration& conf)
  {
    ATH_CHECK (AnaTool::useInitialConfiguration (conf));

    m_isData = conf.isData();
    m_isAFII = conf.isAFII();

    m_jetContainer = conf.inputName (OBJECT_JET);
    if (m_jetContainer.empty())
    {
      ATH_MSG_ERROR ("can't apply correction without jet collection name");
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }



  StatusCode JetToolCorrect ::
  initialize()
  {
    ATH_MSG_DEBUG("initialize");

    // JES MC15 calibration recommendations from:
    //  https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtmissRecommendationsMC15

    // Strip off the "Jets" suffix in the jet collection name
    // @TODO: update AnaToolHandle tool creation mechanism
    ATH_CHECK (ASG_MAKE_ANA_TOOL (m_calibration_tool, JetCalibrationTool));
    const auto jetCollection = m_jetContainer.substr(0, m_jetContainer.size()-4);
#if ROOTCORE_RELEASE_SERIES == 24
    const std::string configFile = m_isAFII ? "JES_MC15Prerecommendation_AFII_June2015.config"
                                            : jesFile;
    const std::string calibSeq = m_isData ? "JetArea_Residual_Origin_EtaJES_GSC_Insitu"
                                          : "JetArea_Residual_Origin_EtaJES_GSC";
#else
    const std::string configFile = m_isAFII ? jesFile_AFII : jesFile;

    const std::string calibSeq = m_isData ? "JetArea_Residual_EtaJES_GSC_Insitu"
                                          : "JetArea_Residual_EtaJES_GSC_Smear";

    ATH_CHECK( m_calibration_tool.setProperty("CalibArea", "00-04-81") );
#endif
    ATH_CHECK( m_calibration_tool.setProperty("JetCollection", jetCollection) );
    ATH_CHECK( m_calibration_tool.setProperty("ConfigFile", configFile) );
    ATH_CHECK( m_calibration_tool.setProperty("CalibSequence", calibSeq) );
    ATH_CHECK( m_calibration_tool.setProperty("IsData", m_isData) );
    ATH_CHECK( m_calibration_tool.initialize() );

    // JES MC15 uncertainty recommendations from:
    //  https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetUncertainties2015Moriond2016

    // NOTE: using the strongly reduced uncertainty breakdown. JetEtMiss says
    // that analyses MUST evaluate all four strongly reduced scenarios before
    // adopting this one!
    // const std::string uncertConfigFile = "JES_2015/ICHEP2016/JES2015_SR_Scenario1.config";
    //const std::string uncertConfigFile = "JES_2015/Moriond2016/JES2015_19NP.config";

    // @TODO: update AnaToolHandle tool creation mechanism
    ATH_CHECK( ASG_MAKE_ANA_TOOL(m_uncertainties_tool, JetUncertaintiesTool) );
    ATH_CHECK( m_uncertainties_tool.setProperty("JetDefinition", jetCollection) );
    ATH_CHECK( m_uncertainties_tool.setProperty("MCType", m_isAFII ? "AFII" : mcType) );
    ATH_CHECK( m_uncertainties_tool.setProperty("ConfigFile", uncertConfigFile) );
#if ROOTCORE_RELEASE_SERIES == 25
    ATH_CHECK( m_uncertainties_tool.setProperty("CalibArea", "CalibArea-03") );
#endif
    ATH_CHECK( m_uncertainties_tool.initialize() );
    registerTool( &*m_uncertainties_tool);

    // JER MC15 pre-recommendations

    // Specify the JER file - Don't want behavior to change with the same tag of
    // QuickAna just because the package underlying was updated (JER in this case)
    const std::string jerFile =
      "JetResolution/Prerec2015_xCalib_2012JER_ReducedTo9NP_Plots_v2.root";
    ATH_CHECK (ASG_MAKE_ANA_TOOL (m_resolution_tool, JERTool));
    ATH_CHECK( m_resolution_tool.setProperty("PlotFileName", jerFile) );
    ATH_CHECK( m_resolution_tool.setProperty("CollectionName", m_jetContainer) );
    ATH_CHECK( m_resolution_tool.initialize() );

    // @TODO: update AnaToolHandle tool creation mechanism
    ATH_CHECK( ASG_MAKE_ANA_TOOL(m_smearing_tool, JERSmearingTool) );
    ATH_CHECK( m_smearing_tool.setProperty("isMC", !m_isData) );
    ATH_CHECK( m_smearing_tool.setProperty("ApplyNominalSmearing", false) );
    ATH_CHECK( m_smearing_tool.setProperty("SystematicMode", "Simple") );
    ATH_CHECK( m_smearing_tool.setProperty
                 ("JERTool", m_resolution_tool) );
    ATH_CHECK( m_smearing_tool.initialize() );
    registerTool( &*m_smearing_tool );

    // JVT tool
    const std::string jvtFile = "JetMomentTools/JVTlikelihood_20140805.root";
    // @TODO: update AnaToolHandle tool creation mechanism
    ATH_CHECK( ASG_MAKE_ANA_TOOL(m_jvt_tool, JetVertexTaggerTool) );
    ATH_CHECK( m_jvt_tool.setProperty("JVTFileName", jvtFile) );
    ATH_CHECK( m_jvt_tool.initialize() );

    // JVT efficiency SF
    //  From https://twiki.cern.ch/twiki/bin/view/AtlasProtected/JVTCalibration
      const std::string jvtEffile = "JetJvtEfficiency/Moriond2018/JvtSFFile_EMTopoJets.root";// changed by shuzhou 09/17/2018
      const std::string fjvtEFFile = "JetJvtEfficiency/Moriond2018/fJvtSFFile.root";
    ATH_CHECK( ASG_MAKE_ANA_TOOL(m_jvtEffTool, CP::JetJvtEfficiency) );
    // The default working point is the recommended one
    ATH_CHECK( m_jvtEffTool.setProperty("SFFile", jvtEffile) );
    ATH_CHECK( m_jvtEffTool.initialize() );
    registerTool (&*m_jvtEffTool);
      
      // FJVT Tool
      ATH_CHECK( ASG_MAKE_ANA_TOOL(m_fjvtEffTool, CP::JetJvtEfficiency) );
      // The default working point is the recommended one
      //  ATH_CHECK( m_jvtEffTool.setProperty("WorkingPoint","Default") );
      ATH_CHECK( m_fjvtEffTool.setProperty("SFFile", fjvtEFFile) );
      ATH_CHECK( m_fjvtEffTool.initialize() );
      registerTool (&*m_fjvtEffTool);

    // b-tagging tools
    if(m_enableBTagging) {
      //const std::string bTagCalibFile =
      //   "xAODBTaggingEfficiency/13TeV/2016-20_7-13TeV-MC15-CDI-May31_v1.root";
      // @TODO: update AnaToolHandle tool creation mechanism
      ATH_CHECK( ASG_MAKE_ANA_TOOL(m_bsel_tool, BTaggingSelectionTool) );
      ATH_CHECK( m_bsel_tool.setProperty("TaggerName", m_btagger) );
      ATH_CHECK( m_bsel_tool.setProperty("OperatingPoint", m_btagWP) );
      ATH_CHECK( m_bsel_tool.setProperty("JetAuthor", m_jetContainer) );
      ATH_CHECK( m_bsel_tool.setProperty("FlvTagCutDefinitionsFileName", bTagCalibFile) );
      ATH_CHECK( m_bsel_tool.initialize() );

      ATH_CHECK( ASG_MAKE_ANA_TOOL(m_bsel_OR_tool, BTaggingSelectionTool) );
      ATH_CHECK( m_bsel_OR_tool.setProperty("TaggerName", m_btagger) );
      ATH_CHECK( m_bsel_OR_tool.setProperty("OperatingPoint", m_btagWP_OR) );
      ATH_CHECK( m_bsel_OR_tool.setProperty("JetAuthor", m_jetContainer) );
      ATH_CHECK( m_bsel_OR_tool.setProperty("FlvTagCutDefinitionsFileName", bTagCalibFile) );
      ATH_CHECK( m_bsel_OR_tool.initialize() );
    }

    // Jet cleaning tool for decoration
    ATH_CHECK (ASG_MAKE_ANA_TOOL (m_cleaning_tool, JetCleaningTool));
    ATH_CHECK (m_cleaning_tool.setProperty("CutLevel", "LooseBad"));
    ATH_CHECK (m_cleaning_tool.setProperty("DoUgly", false));
    ATH_CHECK (m_cleaning_tool.initialize());

    registerCut (SelectionStep::MET, "calibration_tool", cut_calibration_tool);
    registerCut (SelectionStep::MET, "uncertainties_tool", cut_uncertainties_tool);
    registerCut (SelectionStep::MET, "smearing_tool", cut_smearing_tool);

    // Only decorate jets with the information, so that event-level
    // cleaning can be performed later
    registerCut (SelectionStep::MANUAL, "cleaning_tool", cut_cleaning_tool);

    return StatusCode::SUCCESS;
  }



  StatusCode JetToolCorrect ::
  correctObject (xAOD::Jet& jet)
  {
    ATH_MSG_DEBUG("correctObject");
    // jet calibration
    QA_CHECK_CUT( cut_calibration_tool, m_calibration_tool->applyCorrection(jet) );

    // jet energy scale uncertainties
    QA_CHECK_CUT( cut_uncertainties_tool, m_uncertainties_tool->applyCorrection(jet) );

    // jet resolution smearing
    QA_CHECK_CUT( cut_smearing_tool, m_smearing_tool->applyCorrection(jet) );

    // JVT Recalculation after calibration
    float jvt = m_jvt_tool->updateJvt(jet);
    // Update "Jvt" of the jet - required by the MET tool
    jet.auxdata<float>("Jvt") = jvt;
      
      bool jvt_pass =true;
      //bool fjvt_pass = true;
      if(jet.pt()>20.*GeV&&jet.pt()<60.*GeV&&std::abs(jet.eta())<2.4){
          jvt_pass=m_jvtEffTool->passesJvtCut(jet);
          
      }
      
      //else if(jet.pt()>20.*GeV&&jet.pt()<60.*GeV&&std::abs(jet.eta())>2.5){
      //  jvt_pass=m_fjvtEffTool->passesJvtCut(jet);
      
      //}
      else{
          jvt_pass=true;
      }
      //jvt_pass=true;
      //fjvt_pass=true;
      jet.auxdata<char>("Jvt_pass") = jvt_pass;

    // B-tagging section
    if (m_enableBTagging)
    {
      // BTagging valid in pt > 20, |eta| < 2.5
      bool inBTagKinRange = jet.pt() > 20.*GeV && std::abs(jet.eta()) < 2.5;

      // B-Jet criteria
      bool isbjet = ( inBTagKinRange && jvt_pass && m_bsel_tool->accept(jet) );
      jet.auxdecor<char>("bjet") = isbjet;

      // Apply the dedicated bjet decoration for overlap removal as well.
      // Working point can be different from standard one.
      bool isbjet_OR = isbjet;
      if (m_btagWP_OR != m_btagWP) {
        isbjet_OR = ( inBTagKinRange && jvt_pass && m_bsel_OR_tool->accept(jet) );
      }
      jet.auxdecor<char>("bjet_OR") = isbjet_OR;

      // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/JetUncertainties2015Prerec#Tool_requirements_and_assumption
      // Explains that this decoration *must* be called IsBjet for the
      // nuisance parameter breakdown to work.
      // - Why can't we just adopt this as our convention, then??
      jet.auxdecor<char>("IsBjet") = isbjet;
    }

    // We only clean, by default, jets that might've passed our JVT selection.
    // This is too hard-coded, ugly.
    bool is_clean = jet.isAvailable<char>("DFCommonJets_jetClean_LooseBad") ?
                    static_cast<bool>(jet.auxdata<char>("DFCommonJets_jetClean_LooseBad")) :
                    ( jet.pt() < 20.*GeV || (jet.pt()<60.*GeV && !jvt_pass) ||
                   m_cleaning_tool->keep(jet) );
    cut_cleaning_tool.setPassedIf ( is_clean );

    // Also decorate the jet with the information, so that
    // event-level cleaning can be performed later.
    jet.auxdecor<char>("clean_jet") = is_clean;

    return StatusCode::SUCCESS;
  }



  JetToolSelect ::
  JetToolSelect (const std::string& name)
    : AsgTool (name), AnaToolSelect<xAOD::JetContainer> (name),
      m_jvt_cut_step (SelectionStep::MET)
  {
  }



  StatusCode JetToolSelect ::
  initialize()
  {
    ATH_MSG_DEBUG("initialize");

    // For SUSY: Make sure to *not* cut jets before they hit the overlap removal!
    // Just decorate.
    registerCut (m_jvt_cut_step, "jvt", cut_jvt);

    return StatusCode::SUCCESS;
  }



  StatusCode JetToolSelect ::
  selectObject (xAOD::Jet& jet)
  {
    ATH_MSG_DEBUG("selectObject");

    // Apply the recommended selection
    cut_jvt.setPassedIf ( jet.getAttribute<char>("Jvt_pass") );

    return StatusCode::SUCCESS;
  }



  JetToolWeight ::
  JetToolWeight (const std::string& name)
    : AsgTool (name), AnaToolWeight<xAOD::JetContainer> (name),
      m_btagging_eff_tool ("btagging_eff", this),
      m_jvtEffTool("jvt_eff", this),
      m_anaSelect ("ana_select"),
      m_anaWeight ("ana_weight")
  {
    declareProperty("BTagger", m_btagger = btagAlgDefault);
    declareProperty("BTagWP", m_btagWP = "-0_4434");
  }


  unsigned JetToolWeight ::
  inputTypes () const
  {
    return (1 << OBJECT_EVENTINFO) | (1 << OBJECT_JET);
  }


  unsigned JetToolWeight ::
  outputTypes () const
  {
    return (1 << OBJECT_EVENTINFO) | (1 << OBJECT_JET);
  }


  StatusCode JetToolWeight ::
  initialize()
  {
    ATH_MSG_DEBUG("initialize");
      // JVT tool
      const std::string jvtFile = "JetMomentTools/JVTlikelihood_20140805.root";
      const std::string jvtEFFile = "JetJvtEfficiency/Moriond2018/JvtSFFile_EMTopoJets.root";// changed by shuzhou 09/17/2018
      // @TODO: update AnaToolHandle tool creation mechanism
      //const std::string fjvtEFFile = "JetJvtEfficiency/Moriond2018/fJvtSFFile.root";
     // ATH_CHECK( ASG_MAKE_ANA_TOOL(m_jvt_tool, JetVertexTaggerTool) );
     // ATH_CHECK( m_jvt_tool.setProperty("JVTFileName", jvtFile) );
      
     // ATH_CHECK( m_jvt_tool.initialize() );
      
      // JVT efficiency SF
      //  From https://twiki.cern.ch/twiki/bin/view/AtlasProtected/JVTCalibration
      // @TODO update tool creation mechanism
      ATH_CHECK( ASG_MAKE_ANA_TOOL(m_jvtEffTool, CP::JetJvtEfficiency) );
      // The default working point is the recommended one
      //  ATH_CHECK( m_jvtEffTool.setProperty("WorkingPoint","Default") );
      ATH_CHECK( m_jvtEffTool.setProperty("SFFile", jvtEFFile) );
      ATH_CHECK( m_jvtEffTool.initialize() );
      registerTool (&*m_jvtEffTool);
      
      // FJVT Tool
      //ATH_CHECK( ASG_MAKE_ANA_TOOL(m_fjvtEffTool, CP::JetJvtEfficiency) );
      // The default working point is the recommended one
      //  ATH_CHECK( m_jvtEffTool.setProperty("WorkingPoint","Default") );
      //ATH_CHECK( m_fjvtEffTool.setProperty("SFFile", fjvtEFFile) );
      //ATH_CHECK( m_fjvtEffTool.initialize() );
      //registerTool (&*m_fjvtEffTool);
      
      if(m_btagWP.empty()){
          
          return StatusCode::SUCCESS;
          
      }
      
      // BTag efficiency SF
      // Recommendations come from
      //  https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTagCalib2015
      
      //const std::string btagFile =
      //  "xAODBTaggingEfficiency/13TeV/2016-20_7-13TeV-MC15-CDI-May31_v1.root";
      
      // @TODO update tool creation mechanism
      ATH_CHECK( ASG_MAKE_ANA_TOOL(m_btagging_eff_tool, BTaggingEfficiencyTool) );
      ATH_CHECK( m_btagging_eff_tool.setProperty("TaggerName", m_btagger) );
      // really stupid that these have different formats
      ATH_CHECK( m_btagging_eff_tool.setProperty("OperatingPoint", m_btagWP) );
      ATH_CHECK( m_btagging_eff_tool.setProperty("JetAuthor", "AntiKt4EMTopoJets") );
      ATH_CHECK( m_btagging_eff_tool.setProperty("ScaleFactorFileName", bTagCalibFile) );
      ATH_CHECK( m_btagging_eff_tool.setProperty("SystematicsStrategy", "Envelope") );
      ATH_CHECK( m_btagging_eff_tool.initialize() );
      // register for systematics
      registerTool (&*m_btagging_eff_tool);
      
      
      return StatusCode::SUCCESS;
  }


  StatusCode JetToolWeight ::
  execute (IEventObjects& objects)
  {
      

    ConstDataVector<xAOD::JetContainer> jvtjets(SG::VIEW_ELEMENTS);
      const xAOD::JetContainer* TruthJets = 0;
      const xAOD::JetContainer* MJets = 0;
      TruthJets=objects.truthjets();
      MJets=objects.jets();
      
      m_jvtEffTool->tagTruth(MJets,TruthJets);
    for (auto object : *objects.jets())
    {
      float weight = 1;
      if (m_anaSelect (*object))
      {
        ATH_CHECK (this->objectWeight (*object, weight));
        if (!(weight > 0))
        {
          ATH_MSG_WARNING ("object weight of " << weight <<
                             " is not allowed: pt=" << object->pt() <<
                             " eta=" << object->eta() <<
                             " phi=" << object->phi());
          //return StatusCode::FAILURE;
        }
        jvtjets.push_back( object );
      } else
      {
        weight = 1;
      }
      m_anaWeight (*object) = weight;
    }

    float totalSF=1.;

    CP::CorrectionCode ret = m_jvtEffTool->applyAllEfficiencyScaleFactor( jvtjets.asDataVector() , totalSF );

    switch (ret) {
    case CP::CorrectionCode::Error:
      ATH_MSG_ERROR( "Failed to retrieve SF for jet in SUSYTools_xAOD::JVT_SF" );
    case CP::CorrectionCode::OutOfValidityRange:
      ATH_MSG_VERBOSE( "No valid SF for jet in SUSYTools_xAOD::JVT_SF" );
    default:
      ATH_MSG_VERBOSE( " Retrieve SF for jet container in SUSYTools_xAOD::JVT_SF with value " << totalSF );
    }
    objects.eventinfo()->auxdata<float>("JVT_SF") = totalSF;


    return StatusCode::SUCCESS;
  }


  StatusCode JetToolWeight ::
  objectWeight (const xAOD::Jet& jet, float& weight)
  {
    ATH_MSG_DEBUG("objectWeight");
      float jvt_sSF=1.0;
      if(jet.pt() > 20.*GeV && jet.pt() < 60.*GeV&&std::abs(jet.eta()) < 2.4){
          if(jet.auxdata<char>("Jvt_pass")){
              QA_CHECK_WEIGHT( float, jvt_sSF,
                              m_jvtEffTool->getEfficiencyScaleFactor(jet, jvt_sSF) );
          }
          else{
              QA_CHECK_WEIGHT( float, jvt_sSF,
                              m_jvtEffTool->getInefficiencyScaleFactor(jet, jvt_sSF) );
              
          }
      }
      else{
          jvt_sSF=1.0;
      }
      
      
      if(m_btagWP.empty()){
          weight=weight*jvt_sSF;
          return StatusCode::SUCCESS;
          
      }

    // Apply btag efficiency SF.
    // The btag tool only supports pt > 20 and |eta| < 2.5
    if (jet.pt() > 20.*GeV && jet.pt() <= 1.*TeV && std::abs(jet.eta()) < 2.5)
    {
      // If the jet is b-tagged, use the efficiency weight
      if (jet.auxdecor<char>("IsBjet")) {
        QA_CHECK_WEIGHT( float, weight,
                         m_btagging_eff_tool->getScaleFactor(jet, weight) );
      }
      // If it is not, then use the *inefficiency* weight
      else {
        QA_CHECK_WEIGHT( float, weight,
                         m_btagging_eff_tool->getInefficiencyScaleFactor(jet, weight) );
      }
    }

    return StatusCode::SUCCESS;
  }



  // Maker function for jet tools
  StatusCode makeJetTool (DefinitionArgs& args,
                          const std::string& collection,
                          const SelectionStep& jvt_step,
                          const std::string& btagWP)
  {
    using namespace msgObjectDefinition;

    const bool useBTagging = !btagWP.empty();
    const auto config = args.configuration();

    args.add (std::unique_ptr<IAnaTool>
      (new AnaToolRetrieve (args.prefix() + "_retrieve", collection)));

    std::unique_ptr<JetToolCorrect> correctTool
      (new JetToolCorrect (args.prefix() + "_correct"));
    ANA_CHECK( correctTool->setProperty("EnableBTagging", useBTagging) );
    ANA_CHECK( correctTool->setProperty("BTagWP", btagWP) );
    args.add (std::move (correctTool));

    std::unique_ptr<JetToolSelect> selectTool
      (new JetToolSelect (args.prefix() + "_select"));
    selectTool->m_jvt_cut_step = jvt_step;
    args.add (std::move (selectTool));

    // Only apply jet weights to MC.
    // Also, at the moment we're only applying btag-related weights, so
    // disable the weight tool if we're not using btagging
    if (config->isData() == false)
    {
      std::unique_ptr<JetToolWeight> weightTool
        (new JetToolWeight (args.prefix() + "_weight"));
      ANA_CHECK( weightTool->setProperty("BTagWP", btagWP));
      args.add (std::move (weightTool));
    }
    return StatusCode::SUCCESS;
  }

  QUICK_ANA_JET_DEFINITION_MAKER ("antikt04_noBtag",
    makeJetTool (args, "AntiKt4EMTopoJets"))

  QUICK_ANA_JET_DEFINITION_MAKER ("AntiKt4EMTopoJets antikt04",
    makeJetTool (args, "AntiKt4EMTopoJets", SelectionStep::MET,
                 "FixedCutBEff_77"))

  QUICK_ANA_JET_DEFINITION_MAKER ("antikt04_HZZ",
    makeJetTool (args, "AntiKt4EMTopoJets", SelectionStep::ANALYSIS,
                 "FixedCutBEff_85"))

  QUICK_ANA_JET_DEFINITION_MAKER( "AntiKt4EMTopo_SUSY",
    makeJetTool (args, "AntiKt4EMTopoJets", SelectionStep::ANALYSIS,
                 "FixedCutBEff_77"))

  QUICK_ANA_JET_DEFINITION_MAKER( "antikt04_DiMu",
    makeJetTool (args, "AntiKt4EMTopoJets", SelectionStep::MET,
                 "FixedCutBEff_50"))
}
