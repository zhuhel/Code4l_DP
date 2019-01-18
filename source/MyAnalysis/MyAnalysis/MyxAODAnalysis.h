#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include <EventLoop/Algorithm.h>
#include "AsgTools/ToolHandle.h"

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "AnalysisVar.h"


#include "PATInterfaces/SystematicRegistry.h"
#include "PATInterfaces/SystematicCode.h"
#include "PATInterfaces/CorrectionCode.h"
#include "PATInterfaces/SystematicsUtil.h"
#include "PATInterfaces/SystematicVariation.h" 


#include <TTree.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <time.h>
#include <TMacro.h>
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
using namespace std;

namespace ORUtils{
  class EleMuSharedTrkOverlapTool;
}

namespace CP{
  class IPileupReweightingTool;
}

namespace ana{
  class QuickAna;
}

namespace TrigConf{
  class xAODConfigTool;
}

namespace Trig{
  class TrigDecisionTool;
}


class MyxAODAnalysis : public EL::Algorithm, public AnalysisVar
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  public:

    time_t start, end; //!

    xAOD::TEvent *m_event;  //!
    //xAOD::EventInfo* eventInfo; //!
    // count number of events
    bool isMC;
    string set;

    int m_tag1; //!
    int m_tag2; //!
    int m_tag3; //!

    int m_eventCounter; //!
    int m_filter; //!
    int m_nlep;//!
    float m_4mfilter; //!
    float m_4efilter; //!
    float m_2e2mfilter; //!

    float m_2e2mfidu0; //!
    float m_2e2mfidu1; //!
    float m_2e2mfidu2; //!
    float m_2e2mfidu3; //!
    float m_2e2mfidu4; //!

    float m_4mfidu0; //!
    float m_4mfidu1; //!
    float m_4mfidu2; //!
    float m_4mfidu3; //!
    float m_4mfidu4; //!

    float m_4efidu0; //!
    float m_4efidu1; //!
    float m_4efidu2; //!
    float m_4efidu3; //!
    float m_4efidu4; //!

    float sumOfWeights; //!
    float sumOfWeightsSquared; //!

    //// tree varaibles
    TTree* tree; //!


    ana::QuickAna *quickAna; //!
    std::vector<CP::SystematicSet> m_sysList; //!

    //  VOmuon goodm;
    //  VOelectron goode;
    //  VOjet goodj;
    //  VOmet goodmet;

#ifndef __CINT__
    //  CP::MuonCalibrationAndSmearingTool *m_muonCalibrationAndSmearingTool; //! 
    //  CP::MuonSelectionTool *m_muonSelectionTool; //!
    //  AsgElectronLikelihoodTool* myLikelihood; //!
    TrigConf::xAODConfigTool *configTool;//!
    Trig::TrigDecisionTool *trigDecTool; //!
    //OverlapRemovalTool *orTool;//!
    ORUtils::EleMuSharedTrkOverlapTool *orTool;//!
#endif // not __CINT__

    ToolHandle<CP::IPileupReweightingTool>  m_pileReweight;//!
    GoodRunsListSelectionTool *m_grl; //!

    // variables that don't get filled at submission time should be
    // protected from being send from the submission node to the worker
    // node (done by the //!)
  public:
    // Tree *myTree; //!
    // TH1 *myHist; //!
    //std::vector<string>* m_setSysList; //!

    // this is a standard constructor
    MyxAODAnalysis ();
    MyxAODAnalysis (string );

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
    ClassDef(MyxAODAnalysis, 1);

    void ClearFlags(MapType2_Int& map);
    void ClearWeight(MapType2_Double& map);
    void ClearVariables(MapType2_Int& map);
    void ClearVariables(MapType2_Long& map);
    void ClearVariables(MapType2_Float& map);
    void ClearVariables(MapType2_Double& map);
    void ClearVariables(MapType2_TLorentzVector& map);
    void ClearVariables(MapType2_VTLorentzVector& map);
    void ClearVariables(MapType2_Double2D& map);
    void ClearVariables(MapType2_VInt& map);
    void ClearVariables(MapType2_VDouble& map);
    void ClearVariables(MapType2_VFloat& map);
    void InitStrVec(vector<string>& out, string in, string de=",");
    void InitObjSTEP(MapType_VString& STEP, string obj, string steps);
    void InitSetting(MapType2_Int& setmap, string setname, string settings);
    void CreateCountingMap();

    void InitHistVar(string varlist, int nbin, double xmin, double xmax, string cutstep="");
    void InitHist2DVar(string varlist, int nxbin, double xmin, double xmax,
        int nybin, double ymin, double ymax, string cutstep="");
    void InitHist2DVVar(string varlist, int nxbin, double xmin, double xmax,
        int nybin, double ymin, double ymax, string cutstep=""); 
    void InitTreeVar(string varlist, string type);
    void InitVVar(string varlist, int nbin=0, double xmin=0, double xmax=0, string cutstep="xAOD");
    void AddVarIntoTree(TTree *tree, string SYS="NOMINAL", bool isMC=true);
    void ReadTriggers(string& TRIG_list);
    void CreateTreeVarHistoMap(TFile *file);
    void FillHistograms(string sysname);


};

#endif
