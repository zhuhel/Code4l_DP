#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "TSystem.h"
using namespace std;

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

#include "xAODMuon/Muon.h"
#include "xAODEventInfo/EventInfo.h"

#include <QuickAna/QuickAna.h>

#include "MyAnalysis/MyxAODAnalysis.h"
#include "MyAnalysis/AnalysisVar.h"

#include "PathResolver/PathResolver.h"
#include "PileupReweighting/PileupReweightingTool.h"
#include "AssociationUtils/EleMuSharedTrkOverlapTool.h"

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

ort::inputAccessor_t selectAcc("selected");
ort::inputDecorator_t selectDec("selected");
ort::outputAccessor_t overlapAcc("overlaps");
ort::objLinkAccessor_t objLinkAcc("overlapObject");

EL::StatusCode MyxAODAnalysis :: initialize ()
{


  time(&start);

  SETNAME.push_back("physics");
  InitSetting(SETTING,"physics","docorr=1,dosys=0,doweight=1");

  InitStrVec(CHN,"eeee,eemm,mmmm,incl");
  //InitStrVec(STEP_cut,"xAOD,Trigger,FourMore,Quad,Kine,TrigMatch,dR4l,JPsiVeto,EleLoose,Iso,IP,OnShell,Filter,TwoJet,TwoSigJet,dEtaJJ,MJJ");
  //InitStrVec(STEP_cut,"xAOD,FourMore,Quad,Kine,nCalo,dR4l,JPsiVeto,EleLoose,IP,Iso,OnShell,Filter,TwoJet,TwoSigJet,dEtaJJ,MJJ,Cleaning");
  InitStrVec(STEP_cut,"xAOD,Cleaning,FourMore,VetoPair,VetoQuad,Quad,dR4l,JPsiVeto,EleLoose,IP,Iso,Filter");
  //InitStrVec(STEP_cut,"xAOD,FourMore,Quad,Kine,JPsiVeto,EleLoose,Iso,IP,OnShell,Filter,TwoJet,TwoSigJet,dEtaJJ,MJJ");


  // --> original simple step cut for obj mu
  InitObjSTEP(STEP_obj,"mu","All,Pt,Eta,Tool,Z0,D0,OverLap");
  // --> original simple step cut for obj ele
  InitObjSTEP(STEP_obj,"ele","All,Pt,Eta,ObjQ,Z0,ID,OverLap");
  InitObjSTEP(STEP_obj,"jet","All,Eta,PtEta,OverLap,JVT");

  //InitHistVar("mu", 100, 0, 50, "xAOD,TwoJet");
  //InitHistVar("MZZ,PtZZ,MZ1,MZ2,PtZ1,PtZ2", 1000, 0, 1000, "Filter,TwoJet,TwoSigJet,dEtaJJ,MJJ");
  //InitHistVar("PtL1,PtL2,PtL3,PtL4", 1000, 0, 1000, "Filter,TwoJet,TwoSigJet,dEtaJJ,MJJ");
  //InitHistVar("EtaL1,EtaL2,EtaL3,EtaL4", 200, -10, 10, "Filter,TwoJet,TwoSigJet,dEtaJJ,MJJ");
  //InitHistVar("PhiL1,PhiL2,PhiL3,PhiL4", 200, -10, 10, "Filter,TwoJet,TwoSigJet,dEtaJJ,MJJ");
  //InitHistVar("MJJ", 5000, 0, 5000, "TwoJet,TwoSigJet,dEtaJJ,MJJ");
  //InitHistVar("MJ1,MJ2,PtJ1,PtJ2", 1000, 0, 1000, "TwoJet,TwoSigJet,dEtaJJ,MJJ");
  //InitHistVar("EtaJ1,EtaJ2", 200, -10, 10, "TwoJet,TwoSigJet,dEtaJJ,MJJ");
  //InitHistVar("EtaJ1xJ2", 1000, -50, 50, "TwoJet,TwoSigJet,dEtaJJ,MJJ");
  //InitHistVar("PhiJ1,PhiJ2", 200, -10, 10, "TwoJet,TwoSigJet,dEtaJJ,MJJ");
  //InitHistVar("dEtaJJ,dPhiJJ", 200, -10, 10, "TwoJet,TwoSigJet,dEtaJJ,MJJ");

  //InitHistVar("YZ1Star,YZ2Star,Z1Cen,Z2Cen", 1000, -50, 50, "TwoJet,TwoSigJet,dEtaJJ,MJJ");
  //InitHistVar("PtJJOHtJJ,PtZZJJOHtZZJJ", 1000, 0, 10, "TwoJet,TwoSigJet,dEtaJJ,MJJ");

  //InitHistVar("gen_weight", 200,-10,10,"xAOD");
  //InitHistVar("mqq,mjj",2000,0,2000,"xAOD");
  //InitHistVar("Ptj1,Ptj2,Ptq1,Ptq2",1000,0,1000,"xAOD");
  //InitHistVar("Etaj1,Etaj2,Etaq1,Etaq2",2000,-10,10,"xAOD");
  //InitHistVar("Phij1,Phij2,Phiq1,Phiq2",2000,-10,10,"xAOD");

  m_event = wk()->xaodEvent();
  Info("initialize()", "Number of events = %lli", m_event->getEntries() );


  const xAOD::EventInfo* eventInfo = 0;
  if( !m_event->retrieve( eventInfo, "EventInfo").isSuccess() )
  {
    Error("execute ()", "Failed to retrieve EventInfo. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  isMC = eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION );

  sumOfWeights = 0.;
  sumOfWeightsSquared = 0.;

  if(isMC) {

    TTree *MetaData = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));

    if (!MetaData) {
      Error("fileExecute()", "MetaData not found! Exiting.");
      return EL::StatusCode::FAILURE;
    }

    MetaData->LoadTree(0);
    bool m_isDerivation = !MetaData->GetBranch("StreamAOD");
    if(m_isDerivation ){

      const xAOD::CutBookkeeperContainer* completeCBC = 0;
      if(!m_event->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess()){
        Error("initializeEvent()","Failed to retrieve CutBookkeepers from MetaData! Exiting.");
        return EL::StatusCode::FAILURE;
      }

      TString xStream="";
      const xAOD::CutBookkeeper* allEventsCBK=0;
      int maxcycle = -1;

      for ( auto cbk :  *completeCBC ) {
        if ( cbk->name() == "AllExecutedEvents" && TString(cbk->inputStream()).Contains("StreamDAOD")){
          xStream = TString(cbk->inputStream()).ReplaceAll("Stream","");
          std::cout << "xStream = " << xStream << "  (i.e. indentified DxAOD flavour)" << std::endl;
        }

        if ( cbk->name() == "AllExecutedEvents" && cbk->cycle() > maxcycle && cbk->inputStream() == "StreamAOD"){
          maxcycle = cbk->cycle();
          allEventsCBK = cbk;
        }
      }

      if (allEventsCBK) {
        sumOfWeights        = allEventsCBK->sumOfEventWeights();
        sumOfWeightsSquared = allEventsCBK->sumOfEventWeightsSquared();
      }
    }
  }



  m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
  //const char* GRLFilePath_2015 = "GoodRunsLists/data15_13TeV/20160720/data15_13TeV.periodAllYear_DetStatus-v79-repro20-02_DQDefects-00-02-02_PHYS_StandardGRL_All_Good_25ns.xml";
  const char* GRLFilePath_2015 = "GoodRunsLists/data15_13TeV/20170619/data15_13TeV.periodAllYear_DetStatus-v89-pro21-02_Unknown_PHYS_StandardGRL_All_Good_25ns.xml";
  const char* fullGRLFilePath_2015 = gSystem->ExpandPathName (GRLFilePath_2015);
  //const char* GRLFilePath_2016 = "GoodRunsLists/data16_13TeV/20170720/data16_13TeV.periodAllYear_DetStatus-v88-pro20-21_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml";
  const char* GRLFilePath_2016 = "GoodRunsLists/data16_13TeV/20180129/data16_13TeV.periodAllYear_DetStatus-v89-pro21-01_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml";
  const char* fullGRLFilePath_2016 = gSystem->ExpandPathName (GRLFilePath_2016);
  const char* GRLFilePath_2017 = "GoodRunsLists/data17_13TeV/20171130/data17_13TeV.periodAllYear_DetStatus-v97-pro21-13_Unknown_PHYS_StandardGRL_All_Good_25ns_Triggerno17e33prim.xml";
  const char* fullGRLFilePath_2017 = gSystem->ExpandPathName (GRLFilePath_2017);
  std::vector<std::string> vecStringGRL;
  vecStringGRL.push_back(fullGRLFilePath_2015);
  vecStringGRL.push_back(fullGRLFilePath_2016);
  vecStringGRL.push_back(fullGRLFilePath_2017);
  m_grl->setProperty( "GoodRunsListVec", vecStringGRL);
  m_grl->setProperty("PassThrough", false); // if true (default) will ignore result of GRL and will just pass all events
  m_grl->initialize();

  std::vector<std::string> prwFiles;
  string PileConfigFile = "dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root";
  prwFiles.push_back(PileConfigFile);
  std::vector<std::string> lumicalcFiles;
  string LumiCalcFile1 = "GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root";
  string LumiCalcFile2 = "GoodRunsLists/data16_13TeV/20170720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-009.root";
  string LumiCalcFile3 = "GoodRunsLists/data17_13TeV/20171130/physics_25ns_Triggerno17e33prim.lumicalc.OflLumi-13TeV-001.root";
  lumicalcFiles.push_back(LumiCalcFile1);
  lumicalcFiles.push_back(LumiCalcFile2);
  lumicalcFiles.push_back(LumiCalcFile3);

  string TRIG_list = "";
  ReadTriggers( TRIG_list );

  quickAna =new ana::QuickAna("quickana");
  if(!isMC) quickAna->isDataFlag="true";
  quickAna->eventinfoDef = "default";
  quickAna->triggerDef = TRIG_list;
  quickAna->muonDef="darkph";
  quickAna->electronDef="darkph hzhinv_medium";
  //quickAna->muonDef="smzz4l";
  //quickAna->electronDef="smzz4l_veryloose smzz4l hzhinv_medium";
  quickAna->jetKine = "pt > 30e3 && eta < 4.5 && eta >-4.5";
  quickAna->jetDef="antikt04_noBtag"; // no btag
  //quickAna->jetDef="antikt04_HZZ"; // btag setup in HZZ analysis
  quickAna->photonDef="none";
  quickAna->tauDef="none";
  quickAna->metDef="none";
  //quickAna->orDef="default";
  quickAna->orDef="zzllll";
  // disable Pileup reweighting, no prw files yet
  //quickAna->muDataFiles=lumicalcFiles;
  //quickAna->muMcFiles=prwFiles;
  quickAna->initialize().ignore();

  m_sysList= CP::make_systematics_vector (quickAna->recommendedSystematics());
  CP::SystematicCode::enableFailure();


  SYSNAME.push_back("NOCORR");
  int Nsyst = 0;
  for (auto sysListItr : m_sysList){
    cout << "name sys is : " << sysListItr.name() << endl;
    //if(Nsyst > 1) break;
    if(sysListItr.name()=="") SYSNAME.push_back("NOMINAL");
    else SYSNAME.push_back(sysListItr.name());
    Nsyst++;
  }


  CreateCountingMap();

  TFile *file1 = wk()->getOutputFile ("tree_output");
  for(int i=0; i<(int)SYSNAME.size(); i++) {

    string sys=SYSNAME[i];

    if(!isMC && sys!="NOMINAL") continue;

    if(sys!="NOCORR" && (!SETTING["physics"]["docorr"])) continue;

    if(sys!="NOMINAL" && SETTING["physics"]["docorr"] && (!SETTING["physics"]["dosys"])) continue;

    if(sys=="NOCORR" && SETTING["physics"]["docorr"] && SETTING["physics"]["dosys"]) continue;

    //  if(sys!="NOMINAL" && SYS_leaf[sys]=="ON") continue;

    string tree_name = "tree_" + sys;

    Tree[sys] = new TTree(tree_name.c_str(), "output tree");
    Tree[sys]->SetDirectory (file1);
    AddVarIntoTree(Tree[sys], sys, isMC);

  }


  TFile *file2 = wk()->getOutputFile ("hist_output");
  CreateTreeVarHistoMap(file2);

  m_eventCounter = 0;
  m_filter = 0;
  m_nlep = 0;
  m_4mfilter = 0;
  m_4efilter = 0;
  m_2e2mfilter = 0;

  m_2e2mfidu0 = 0;
  m_2e2mfidu1 = 0;
  m_2e2mfidu2 = 0;
  m_2e2mfidu3 = 0;
  m_2e2mfidu4 = 0;

  m_4mfidu0 = 0;
  m_4mfidu1 = 0;
  m_4mfidu2 = 0;
  m_4mfidu3 = 0;
  m_4mfidu4 = 0;

  m_4efidu0 = 0;
  m_4efidu1 = 0;
  m_4efidu2 = 0;
  m_4efidu3 = 0;
  m_4efidu4 = 0;

  m_tag1 = 0;
  m_tag2 = 0;
  m_tag3 = 0;

  return EL::StatusCode::SUCCESS;

}


void MyxAODAnalysis :: ClearFlags(MapType2_Int& map) {

  MapType2_Int::iterator it;
  for(it=map.begin(); it!=map.end(); it++) {
    string chn=(*it).first;
    MapType_Int::iterator it2;
    for(it2=(*it).second.begin(); it2!=(*it).second.end(); it2++) {
      string cut=(*it2).first;
      map[chn][cut] = 0;
    }
  }
}

void MyxAODAnalysis :: ClearWeight(MapType2_Double& map) {
  MapType2_Double::iterator it;

  MapType_Double::iterator iit;

  for(it=map.begin(); it!=map.end(); it++)
    for(iit=(*it).second.begin(); iit!=(*it).second.end(); iit ++)
      map[(*it).first][(*iit).first]=1.0;
}



void MyxAODAnalysis :: CreateCountingMap() {

  MapType_VString::iterator it;
  for(it=STEP_obj.begin(); it!=STEP_obj.end(); it++) {
    string obj=(*it).first;
    for(int i=0; i<(int)(*it).second.size(); i++) {
      string cut=STEP_obj[obj][i];
      COUNT ini={0.,0.};

      for(int j=0; j<(int)SYSNAME.size(); j++) {
        CNT_obj[SYSNAME[j]][obj][cut]=ini;
      }

    }
  }

  for(int i=0; i<(int)CHN.size(); i++) {
    string chn=CHN[i];
    for(int j=0; j<(int)STEP_cut.size(); j++) { 
      string cut=STEP_cut[j];
      FLAG_cut_temp[chn][cut]=0;
      FLAG_cut[chn][cut]=0;
      COUNT ini={0.,0.};
      Evt_Weight[chn][cut]=1.0;

      for(int k=0; k<(int)SYSNAME.size(); k++) { 
        CNT_cut[SYSNAME[k]][chn][cut]=ini;
      }
    }
  }
}


void MyxAODAnalysis :: InitStrVec(vector<string>& out, string in, string de) {
  int pos=0, pos_pre=0;
  while(true) {
    pos=in.find(de,pos_pre);
    if(pos==-1) {out.push_back(in.substr(pos_pre,in.size()-pos_pre)); break;}
    else  out.push_back(in.substr(pos_pre,pos-pos_pre));
    pos_pre=pos+1;
  }
}



void MyxAODAnalysis :: InitObjSTEP(MapType_VString& STEP, string obj, string steps) {
  vector<string> str, objstr;
  InitStrVec(str, steps, ",");

  for(int i=0; i<(int)str.size(); i++) {
    string objs = str[i];
    objstr.push_back(objs);
  }

  STEP[obj]=objstr;
}

void MyxAODAnalysis :: ClearVariables(MapType2_Int& map) {

  MapType2_Int::iterator it;
  for(it=map.begin(); it!=map.end(); it++) {
    string varname = (*it).first;
    map[varname]["Value"]=-9999;
  }
}

void MyxAODAnalysis :: ClearVariables(MapType2_Long& map) {

  MapType2_Long::iterator it;
  for(it=map.begin(); it!=map.end(); it++) {
    string varname = (*it).first;
    map[varname]["Value"]=-9999;
  }
}


void MyxAODAnalysis :: ClearVariables(MapType2_Float& map) {

  MapType2_Float::iterator it;
  for(it=map.begin(); it!=map.end(); it++) {
    string varname = (*it).first;
    map[varname]["Value"]=-9999.;
  }
}

void MyxAODAnalysis :: ClearVariables(MapType2_TLorentzVector& map) {

  MapType2_TLorentzVector::iterator it;
  for(it=map.begin(); it!=map.end(); it++) {
    string varname = (*it).first;
    TLorentzVector *tlv = new TLorentzVector();
    tlv->SetPtEtaPhiM(0.,0.,0.,0.);
    map[varname]["Value"]=*tlv;
  }
}

void MyxAODAnalysis :: ClearVariables(MapType2_Double& map) {

  MapType2_Double::iterator it;
  for(it=map.begin(); it!=map.end(); it++) {
    string varname = (*it).first;
    map[varname]["Value"]=-9999.;
  }
}

void MyxAODAnalysis :: ClearVariables(MapType2_Double2D& map) {

  MapType2_Double2D::iterator it;
  pair<double, double> init (-9999.0, -9999.0);
  for(it=map.begin(); it!=map.end(); it++) {
    string varname = (*it).first;
    map[varname]["Value"] = init;
  }
}

void MyxAODAnalysis :: ClearVariables(MapType2_VInt& map) {

  MapType2_VInt::iterator it;
  for(it=map.begin(); it!=map.end(); it++) {
    string varname = (*it).first;
    map[varname]["Value"].clear();
  }
}


void MyxAODAnalysis :: ClearVariables(MapType2_VDouble& map) {

  MapType2_VDouble::iterator it;
  for(it=map.begin(); it!=map.end(); it++) {
    string varname = (*it).first;
    map[varname]["Value"].clear();
    map[varname]["Weight"].clear();
  }
}

void MyxAODAnalysis :: ClearVariables(MapType2_VFloat& map) {

  MapType2_VFloat::iterator it;
  for(it=map.begin(); it!=map.end(); it++) {
    string varname = (*it).first;
    map[varname]["Value"].clear();
  }
}

void MyxAODAnalysis :: ClearVariables(MapType2_VTLorentzVector& map) {

  MapType2_VTLorentzVector::iterator it;
  for(it=map.begin(); it!=map.end(); it++) {
    string varname = (*it).first;
    map[varname]["Value"].clear();
  }
}


void MyxAODAnalysis :: InitTreeVar(string varlist, string type) {
  vector<string> variables;
  InitStrVec(variables, varlist, ",");

  if(type=="I") {
    for(int i=0; i<(int)variables.size(); i++) {
      TreeIntVar[variables[i]]["Value"]=-9999;
    }
  }

  if(type=="F") {
    for(int i=0; i<(int)variables.size(); i++) {
      TreeFltVar[variables[i]]["Value"]=-9999.0;
    }
  }

  if(type=="D") {
    for(int i=0; i<(int)variables.size(); i++) {
      TreeDouVar[variables[i]]["Value"]=-9999.0;
    }
  }

  if(type=="O") {
    for(int i=0; i<(int)variables.size(); i++) {
      TreeBoolVar[variables[i]]["Value"]=false;
    }
  }

  if(type=="TLV") {
    for(int i=0; i<(int)variables.size(); i++) {
      TLorentzVector *tlv = new TLorentzVector();
      tlv->SetPtEtaPhiM(0.,0.,0.,0.);
      TreeTLVVar[variables[i]]["Value"]=*tlv;
    }
  }

  if(type=="vector<TLV>") {
    for(int i=0; i<(int)variables.size(); i++) {
      TreeTLVVVar[variables[i]]["Value"].clear();
    }
  }

  if(type=="l" || type=="L") {
    for(int i=0; i<(int)variables.size(); i++) {
      TreeLngVar[variables[i]]["Value"]=-9999;
    }
  }

  if(type=="vector<F>") {
    for(int i=0; i<(int)variables.size(); i++) {
      TreeFltVVar[variables[i]]["Value"].clear();
    }
  }

  if(type=="vector<I>") {
    for(int i=0; i<(int)variables.size(); i++) {
      TreeIntVVar[variables[i]]["Value"].clear();
    }
  }



}

void MyxAODAnalysis :: InitHistVar(string varlist, int nbin, double xmin, double xmax, string cutstep) {
  vector<string> variables;
  InitStrVec(variables, varlist, ",");

  vector<string> cuts;
  if(cutstep=="") cuts.clear();
  else if(cutstep=="All") cuts = STEP_cut;
  else if(cutstep.compare(0,1,"-")==0) {
    cutstep = cutstep.erase(0,1);
    for(int i=0; i<(int)STEP_cut.size(); i++) {
      if(cutstep != STEP_cut[i]) cuts.push_back(STEP_cut[i]);
      else if(cutstep == STEP_cut[i]) {cuts.push_back(STEP_cut[i]); break;}
    }
  }else {
    InitStrVec(cuts, cutstep, ",");
  }

  for(int i=0; i<(int)variables.size(); i++) {

    HistVar[variables[i]]["Value"]=-9999.0;
    HistVar[variables[i]]["NBins"]=nbin;
    HistVar[variables[i]]["Xmin"]=xmin;
    HistVar[variables[i]]["Xmax"]=xmax;
    //Var[variables[i]]["Vector"]=0;

    for(int j=0; j<(int)cuts.size(); j++) {
      HistVar[variables[i]][cuts[j]]=1;
    }
  }

}

void MyxAODAnalysis :: InitHist2DVar(string varlist, int nxbin, double xmin, double xmax,
    int nybin, double ymin, double ymax, string cutstep) {
  vector<string> variables;
  InitStrVec(variables, varlist, ",");

  vector<string> cuts;
  if(cutstep=="") cuts.clear();
  else if(cutstep=="All") cuts = STEP_cut;
  else if(cutstep.compare(0,1,"-")==0) {
    cutstep = cutstep.erase(0,1);
    for(int i=0; i<(int)STEP_cut.size(); i++) {
      if(cutstep != STEP_cut[i]) cuts.push_back(STEP_cut[i]);
      else if(cutstep == STEP_cut[i]) {cuts.push_back(STEP_cut[i]); break;}
    }
  }else {
    InitStrVec(cuts, cutstep, ",");
  }

  pair<double, double> init (-9999.0, -9999.0);
  for(int i=0; i<(int)variables.size(); i++) {

    Hist2DVar[variables[i]]["Value"]=init;
    Hist2DVar[variables[i]]["NBins"]=make_pair(nxbin, nybin);
    Hist2DVar[variables[i]]["Min"]=make_pair(xmin, ymin);
    Hist2DVar[variables[i]]["Max"]=make_pair(xmax, ymax);

    for(int j=0; j<(int)cuts.size(); j++) {
      Hist2DVar[variables[i]][cuts[j]]=make_pair(1,1);
    }
  }
}


void MyxAODAnalysis :: InitVVar(string varlist, int nbin, double xmin, double xmax, string cutstep) {
  vector<string> variables;
  InitStrVec(variables, varlist, ",");

  vector<string> cuts;
  if(cutstep=="") cuts.clear();
  else if(cutstep=="All") cuts = STEP_cut;
  else if(cutstep.compare(0,1,"-")==0) {
    cutstep = cutstep.erase(0,1);
    for(int i=0; i<(int)STEP_cut.size(); i++) {
      if(cutstep != STEP_cut[i]) cuts.push_back(STEP_cut[i]);
      else if(cutstep == STEP_cut[i]) {cuts.push_back(STEP_cut[i]); break;}
    }
  }else {
    InitStrVec(cuts, cutstep, ",");
  }

  for(int i=0; i<(int)variables.size(); i++) {

    VVar[variables[i]]["Value"].clear();
    VVar[variables[i]]["NBins"].push_back(nbin);
    VVar[variables[i]]["Xmin"].push_back(xmin);
    VVar[variables[i]]["Xmax"].push_back(xmax);
    VVar[variables[i]]["Weight"].clear();

    for(int j=0; j<(int)cuts.size(); j++) {
      VVar[variables[i]][cuts[j]].push_back(1);
    }
  }

}

void MyxAODAnalysis :: InitHist2DVVar(string varlist, int nxbin, double xmin, double xmax,int nybin, double ymin, double ymax, string cutstep) {
  vector<string> variables;
  InitStrVec(variables, varlist, ",");

  vector<string> cuts;
  if(cutstep=="") cuts.clear();
  else if(cutstep=="All") cuts = STEP_cut;
  else if(cutstep.compare(0,1,"-")==0) {
    cutstep = cutstep.erase(0,1);
    for(int i=0; i<(int)STEP_cut.size(); i++) {
      if(cutstep != STEP_cut[i]) cuts.push_back(STEP_cut[i]);
      else if(cutstep == STEP_cut[i]) {cuts.push_back(STEP_cut[i]); break;}
    }
  }else {
    InitStrVec(cuts, cutstep, ",");
  }

  for(int i=0; i<(int)variables.size(); i++) {

    V2DVar[variables[i]]["Value"].clear();
    V2DVar[variables[i]]["NBins"].push_back(make_pair(nxbin, nybin));
    V2DVar[variables[i]]["Min"].push_back(make_pair(xmin, ymin));
    V2DVar[variables[i]]["Max"].push_back(make_pair(xmax, ymax));

    for(int j=0; j<(int)cuts.size(); j++) {
      V2DVar[variables[i]][cuts[j]].push_back(make_pair(1,1));
    }
  }
}

void MyxAODAnalysis :: ReadTriggers(string& TRIG_list) {

  ifstream file;
  string varfile = PathResolverFindCalibFile("MyAnalysis/Triggers.txt");
  file.open(varfile.c_str(), ios::out);

  //FILE* fp;
  //char result_buf[1000];
  //fp = popen("find $ROOTCOREBIN/data/MyAnalysis/ -name Triggers.txt", "r");
  //fgets(result_buf, sizeof(result_buf), fp);

  //string varfile(result_buf);
  //size_t len = varfile.size();
  //varfile.erase(len-1);
  //ifstream file;
  //file.open(varfile.c_str(), ios::out);

  if (file.is_open())  {
    char line[256];
    while (!file.eof() )  {
      string varname, type;

      file.getline (line,200);
      string sline(line);

      if(sline.find("HLT")!=string::npos) {
        size_t spa_pos = sline.find(" ");
        string hltname = sline.substr(0, spa_pos);
        string status = sline.substr(spa_pos);
        boost::algorithm::trim(hltname);
        boost::algorithm::trim(status);
        TRIG_list = TRIG_list+" "+hltname;
        TRIG_vec[hltname] = status;
      }
    }
  }

  boost::algorithm::trim(TRIG_list);

}



void MyxAODAnalysis :: AddVarIntoTree(TTree *tree, string SYS, bool isMC) {

  string varfile = PathResolverFindCalibFile("MyAnalysis/MiniTree.txt");
  ifstream file;
  file.open(varfile.c_str(), ios::out);

  //FILE* fp;
  //char result_buf[1000];
  //fp = popen("find $ROOTCOREBIN/data/MyAnalysis/ -name MiniTree.txt", "r");
  //fgets(result_buf, sizeof(result_buf), fp);

  //string varfile(result_buf);
  //size_t len = varfile.size();
  //varfile.erase(len-1);
  //ifstream file;
  //file.open(varfile.c_str(), ios::out);

  if (file.is_open())  {
    char line[256];
    while (!file.eof() )  {
      string varname, type;

      file.getline (line,100);
      string sline(line);

      if(sline.find("#")!=string::npos) continue; // use # for comment line in MiniTree.txt
      std::size_t pos = sline.find_last_of(" ");
      varname = sline.substr(pos+1);

      if(sline.find("Int_t")!=string::npos) {
        type = "I";
        //varname = sline.substr(9);
        type = varname+"/"+type;
        InitTreeVar(varname,"I");
        tree->Branch(varname.c_str(),&TreeIntVar[varname]["Value"], type.c_str());
      }
      else if(sline.find("Float_t")!=string::npos) {
        type = "F";
        //varname = sline.substr(9);
        type = varname+"/"+type;
        InitTreeVar(varname, "F");
        tree->Branch(varname.c_str(),&TreeFltVar[varname]["Value"], type.c_str());
      }
      else if(sline.find("Double_t")!=string::npos) {
        type = "D";
        //varname = sline.substr(9);
        type = varname+"/"+type;
        InitTreeVar(varname, "D");
        tree->Branch(varname.c_str(),&TreeDouVar[varname]["Value"], type.c_str());
      }
      else if(sline.find("Bool_t")!=string::npos) {
        type = "O";
        //varname = sline.substr(9);
        type = varname+"/"+type;
        InitTreeVar(varname, "O");
        tree->Branch(varname.c_str(),&TreeBoolVar[varname]["Value"], type.c_str());
      }
      else if(sline.find("TLorentzVector_t")!=string::npos) {
        type = "TLV";
        //varname = sline.substr(9);
        type = varname+"/"+type;
        InitTreeVar(varname, "TLV");
        //tree->Branch(varname.c_str(),&TreeTLVVar[varname]["Value"], type.c_str());
        tree->Branch(varname.c_str(),&TreeTLVVar[varname]["Value"]);
      }
      else if(sline.find("ULong")!=string::npos) {
        type = "l";
        //varname = sline.substr(9);
        type = varname+"/"+type;
        InitTreeVar(varname, "l");
        tree->Branch(varname.c_str(),&TreeLngVar[varname]["Value"], type.c_str());
      }
      else if(sline.find("vector<F>")!=string::npos) {
        //varname = sline.substr(11);
        InitTreeVar(varname, "vector<F>");
        tree->Branch(varname.c_str(), &TreeFltVVar[varname]["Value"]);
      }
      else if(sline.find("vector<I>")!=string::npos) {
        //varname = sline.substr(11);
        InitTreeVar(varname, "vector<I>");
        tree->Branch(varname.c_str(), &TreeIntVVar[varname]["Value"]);
      }
      else if(sline.find("vector<TLV>")!=string::npos) {
        //varname = sline.substr(11);
        InitTreeVar(varname, "vector<TLV>");
        tree->Branch(varname.c_str(), &TreeTLVVVar[varname]["Value"]);
      }

    }
  }


  if(isMC && SYS=="NOMINAL") {
    MapType_String::iterator it;
    for(it=SYS_leaf.begin(); it!=SYS_leaf.end(); it++) {
      string varname=(*it).first;
      if(SYS_leaf[varname]=="OFF") continue;
      string leafname = "weight_"+varname;
      string type = leafname+"/D";
      tree->Branch(leafname.c_str(), &SYS_weight[varname], type.c_str());
    }
  }

}

void MyxAODAnalysis :: CreateTreeVarHistoMap(TFile* file) {

  for(int k=0; k<(int)SYSNAME.size(); k++) {
    string sysname = SYSNAME[k];

    if(sysname!="NOCORR" && (!SETTING["physics"]["docorr"])) continue;

    if(sysname!="NOMINAL" && SETTING["physics"]["docorr"] && (!SETTING["physics"]["dosys"])) continue;

    if(sysname=="NOCORR" && SETTING["physics"]["docorr"] && SETTING["physics"]["dosys"]) continue;

    MapType2_Double::iterator it;
    for(it=HistVar.begin(); it!=HistVar.end(); it++) {
      string varname=(*it).first;
      int nbin = int((*it).second["NBins"]);
      double xlow = double((*it).second["Xmin"]);
      double xhigh = double((*it).second["Xmax"]);

      if(nbin == 0) continue;

      for(int i=0; i<(int)CHN.size(); i++) {
        string chn=CHN[i];

        for(int j=0; j<(int)STEP_cut.size(); j++) {
          string cut=STEP_cut[j];
          if(int((*it).second[cut])!=1) continue;

          string histo_name = sysname + "_" + chn + "_" + cut + "_" + varname;
          TH1F *histo_pointer = new TH1F(histo_name.c_str(),histo_name.c_str(),nbin,xlow,xhigh);
          histo[sysname][chn][cut][varname]=histo_pointer;
          histo[sysname][chn][cut][varname]->SetDirectory(file);
        }
      }
    }

    MapType2_Double2D::iterator it2D;
    for(it2D=Hist2DVar.begin(); it2D!=Hist2DVar.end(); it2D++) {
      string varname=(*it2D).first;

      int nxbin = int((*it2D).second["NBins"].first);
      double xlow = double((*it2D).second["Min"].first);
      double xhigh = double((*it2D).second["Max"].first);
      int nybin = int((*it2D).second["NBins"].second);
      double ylow = double((*it2D).second["Min"].second);
      double yhigh = double((*it2D).second["Max"].second);

      if(nxbin == 0) continue;

      for(int i=0; i<(int)CHN.size(); i++) {
        string chn=CHN[i];

        for(int j=0; j<(int)STEP_cut.size(); j++) {
          string cut=STEP_cut[j];

          pair<double, double> ncuts = (*it2D).second[cut];
          if(int(ncuts.first)!=1) continue;

          string histo_name = sysname + "_" + chn + "_" + cut + "_" + varname;
          TH2F *histo_pointer = new TH2F(histo_name.c_str(),histo_name.c_str(),nxbin,xlow,xhigh,nybin,ylow,yhigh);
          histo_2D[sysname][chn][cut][varname]=histo_pointer;
          histo_2D[sysname][chn][cut][varname]->SetDirectory(file);
        }
      }
    }


    MapType2_VDouble::iterator itv;
    for(itv=VVar.begin(); itv!=VVar.end(); itv++) {
      string varname=(*itv).first;
      int nbin = int((*itv).second["NBins"][0]);
      double xlow = double((*itv).second["Xmin"][0]);
      double xhigh = double((*itv).second["Xmax"][0]);

      if(nbin == 0) continue;

      for(int i=0; i<(int)CHN.size(); i++) {
        string chn=CHN[i];

        for(int j=0; j<(int)STEP_cut.size(); j++) {
          string cut=STEP_cut[j];
          if(int((*itv).second[cut].size())!=1) continue;
          if(int((*itv).second[cut][0])!=1) continue;

          string histo_name = sysname + "_" + chn + "_" + cut + "_" + varname;
          TH1F *histo_pointer = new TH1F(histo_name.c_str(),histo_name.c_str(),nbin,xlow,xhigh);
          histo[sysname][chn][cut][varname]=histo_pointer;
          histo[sysname][chn][cut][varname]->SetDirectory(file);
        }
      }
    }
  }

}

void MyxAODAnalysis :: FillHistograms(string sysname) {

  MapType2_Double::iterator it;
  for(it=HistVar.begin(); it!=HistVar.end(); it++) {
    string varname=(*it).first;

    if((*it).second["NBins"] == 0) continue;

    for(int i=0; i<(int)CHN.size(); i++) {
      string chn=CHN[i];

      for(int j=0; j<(int)STEP_cut.size(); j++) {
        string cut=STEP_cut[j];
        if(int((*it).second[cut])==1 && FLAG_cut[chn][cut]) 
          histo[sysname][chn][cut][varname]->Fill(HistVar[varname]["Value"]);
      } 
    }
  } 

  MapType2_Double2D::iterator it2D;
  for(it2D=Hist2DVar.begin(); it2D!=Hist2DVar.end(); it2D++) {
    string varname=(*it2D).first;

    pair<double, double> nbins = (*it2D).second["NBins"];
    if(nbins.first == 0) continue;

    for(int i=0; i<(int)CHN.size(); i++) {
      string chn=CHN[i];

      for(int j=0; j<(int)STEP_cut.size(); j++) {
        string cut=STEP_cut[j];
        pair<double, double> ncuts = (*it2D).second[cut];
        if(int(ncuts.first)==1 && FLAG_cut[chn][cut]) {
          pair<double, double> values = (*it2D).second["Value"];
          histo_2D[sysname][chn][cut][varname]->Fill(values.first, values.second);
        }
      }
    }
  }


  MapType2_VDouble::iterator itv;
  for(itv=VVar.begin(); itv!=VVar.end(); itv++) {
    string varname=(*itv).first;

    if((*itv).second["NBins"][0] == 0) continue;

    for(int i=0; i<(int)CHN.size(); i++) {
      string chn=CHN[i];

      for(int j=0; j<(int)STEP_cut.size(); j++) {
        string cut=STEP_cut[j];
        if(int((*itv).second[cut].size())!=1) continue;
        if(int((*itv).second[cut][0])==1 && FLAG_cut[chn][cut]) {
          for(int k=0; k<(int)VVar[varname]["Value"].size(); k++)
            if(VVar[varname]["Weight"].size()>0)
              histo[sysname][chn][cut][varname]->Fill(VVar[varname]["Value"][k], VVar[varname]["Weight"][k]);
            else histo[sysname][chn][cut][varname]->Fill(VVar[varname]["Value"][k]);
        }
      }
    }
  }

}


void MyxAODAnalysis :: InitSetting(MapType2_Int& setmap, string setname, string settings) {
  vector<string> vsettings;
  InitStrVec(vsettings, settings, ",");
  for(int i=0; i<(int)vsettings.size(); i++) {
    vector<string> pairs;
    InitStrVec(pairs,vsettings[i],"=");
    if(pairs.size()<2) {
      cout<<"Error in setting: can not parse "<<vsettings[i]<<endl;
      exit(-1);
    }
    setmap[setname][pairs[0]]=atoi(pairs[1].c_str());
  }
}

