#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoopAlgs/AlgSelect.h>

#include "MyAnalysis/MyxAODAnalysis.h"
#include "MyAnalysis/OBJ_Base.h"

#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticlexAODHelpers.h"
#include "xAODJet/JetContainer.h"
#include <xAODMuon/Muon.h>
#include <xAODMuon/MuonContainer.h>
#include <xAODMuon/MuonAuxContainer.h>
#include "xAODEventInfo/EventInfo.h"

#include "xAODTruth/TruthEvent.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTruth/TruthVertexContainer.h"

#include <QuickAna/QuickAna.h>

#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"

#include "Initialize.h"
#include "Counting.h"


using namespace std;

// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)



MyxAODAnalysis :: MyxAODAnalysis ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  //m_setSysList=new std::vector<std::string>();
}

MyxAODAnalysis :: MyxAODAnalysis (string treename="physics")
{

  set = treename;
  //m_setSysList=new std::vector<std::string>();

}



EL::StatusCode MyxAODAnalysis :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  EL::OutputStream output1("tree_output");
  job.outputAdd (output1);

  EL::OutputStream output2("hist_output");
  job.outputAdd (output2);

  EL::OutputStream output3("cutflow");
  job.outputAdd (output3);




  job.useXAOD ();

  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init( "MyxAODAnalysis" ).ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}

// < move the "initialize ()" to a separate file "initialize.h"
/*
   EL::StatusCode MyxAODAnalysis :: initialize ()
   {
// Here you do everything that you need to do after the first input
// file has been connected and before the first event is processed,
// e.g. create additional histograms based on which variables are
// available in the input files.  You can also create all of your
// histograms and trees in here, but be aware that this method
// doesn't get called if no events are processed.  So any objects
// you create here won't be available in the output if you have no
// input events.

//  m_event = wk()->xaodEvent();

// as a check, let's see the number of events in our xAOD
//  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

Varinitialize();

return EL::StatusCode::SUCCESS;
}
*/

EL::StatusCode MyxAODAnalysis :: execute ()
{

  m_eventCounter++;
  if(m_eventCounter==1 || m_eventCounter%5000==0)
    cout << "event counter " << m_eventCounter << endl;

  const xAOD::EventInfo* eventInfo = 0;
  if( !m_event->retrieve( eventInfo, "EventInfo").isSuccess() )
  {
    Error("execute ()", "Failed to retrieve EventInfo. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //  float aveIntPC = eventInfo->actualInteractionsPerCrossing() + eventInfo->averageInteractionsPerCrossing();
  float ave_mu = eventInfo->averageInteractionsPerCrossing();
  uint32_t lb = eventInfo->lumiBlock();
  uint32_t run;
  unsigned long long event;
  vector<float> mcEvtWts;
  double pileWeight=1.0, gen_weight = 1.0;
  double weight=1.0;
  bool passJetCleaning = true;

  //cout << "ave_mu = " << ave_mu << endl;

  if(isMC) {
    run = eventInfo->mcChannelNumber();
    event = eventInfo->mcEventNumber();
    if(SETTING[SETNAME[0]]["doweight"]==1 && isMC) {
      mcEvtWts = eventInfo->mcEventWeights();
      gen_weight = mcEvtWts.size()>0?mcEvtWts[0]:1.0;
    }

  }else {
    run = eventInfo->runNumber();
    event = eventInfo->eventNumber();
  }

  if(!isMC){ // it's data!
    if(!m_grl->passRunLB(*eventInfo)){
      return EL::StatusCode::SUCCESS; // go to next event
    }
  } // end if not MC

  /// Change JetCleaning (a event-level cut) to the flag of eventClean in DAOD
  /// Use working point of "LooseBad"
  //passJetCleaning = eventInfo->DFCommonJets_eventClean_LooseBad();
  bool isEventCleanAvailable = eventInfo->isAvailable<char>("DFCommonJets_eventClean_LooseBad");
  passJetCleaning = eventInfo->auxdataConst<char>("DFCommonJets_eventClean_LooseBad");
  std::cout << "isEventCleanAvailable = " << isEventCleanAvailable << std::endl;
  std::cout << "passJetCleaning = " << passJetCleaning << std::endl;

  VOmuon goodm,goodmiso;
  VOelectron goode,goodeiso;
  VOjet goodj; 


  vector<TLorentzVector> p4muonp, p4muonm, p4elep, p4elem;
  p4muonp.clear(); p4muonm.clear(); p4elep.clear(); p4elem.clear();

  vector<TLorentzVector> sp4muonp, sp4muonm, sp4elep, sp4elem;
  sp4muonp.clear(); sp4muonm.clear(); sp4elep.clear(); sp4elem.clear();

  vector<TLorentzVector> dp4muonm, dp4muonp, dp4elem, dp4elep;
  dp4muonm.clear(); dp4muonp.clear(); dp4elem.clear(); dp4elep.clear();

  vector<TLorentzVector> fp4muonm, fp4muonp, fp4elem, fp4elep;
  fp4muonm.clear(); fp4muonp.clear(); fp4elem.clear(); fp4elep.clear();

  vector<TLorentzVector> vp4m, vp4e, vp4Z, vp4g; 
  vp4m.clear(); vp4e.clear(); 
  vp4Z.clear(); vp4g.clear(); 

  vector<TLorentzVector> Truth_l, Truth_Z;
  vector<int> Truth_pid;

  //  if(m_eventCounter!=45) return EL::StatusCode::SUCCESS;
  //  cout << endl;
  //  cout << "event number is " << event <<endl;
  double MTZZ=-999.0e3, ptZZ_truth=-999.e3;
  bool isTruth=false;
  vector<TLorentzVector> Truth_q;

  // TruthParticle
  if(isMC) {
    // OnShell ZZ
    bool onShell=false, onShell4m=false, onShell4e=false, onShell2e2m=false;
    bool fiducial1=true, fiducial2=true, fiducial3=true, fiducial4=false;

    bool WithTau=false;

    xAOD::TruthParticleContainer const * truth_particles = 0;
    if( !m_event->retrieve( truth_particles, "TruthParticles").isSuccess() )
    {
      Error("excute()", "Failed to retrieve Truth info. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    int nW = 0;
    xAOD::TruthParticleContainer::const_iterator trEvt_itr = truth_particles->begin();
    xAOD::TruthParticleContainer::const_iterator trEvt_end = truth_particles->end();
    for( ; trEvt_itr != trEvt_end; ++trEvt_itr ) {
      const xAOD::TruthParticle* trPart = (*trEvt_itr);
      if(!trPart) continue;
      int status = trPart->status();
      int PDG = trPart->pdgId();
      bool ProdVx = trPart->hasProdVtx();
      bool DecayVx = trPart->hasDecayVtx();
      int barcode = trPart->barcode();
      TLorentzVector p4 = trPart->p4();

      //cout << "Event: " << event;
      //cout << ", PDG=" << PDG << ", " << p4.Pt()/1000. << ", " << p4.Eta() << ", " << p4.Phi() << ", " << p4.M()/1000.;
      //cout << ", status=" << status << ", barcode=" << barcode << endl;

      const xAOD::TruthParticle* parPart = trPart->parent();
      int PDGpar=-1, statuspar=-1;
      TLorentzVector p4par(0,0,0,0);
      if(parPart) {
        PDGpar = parPart->pdgId();
        statuspar = parPart->status();
        p4par = parPart->p4();
      }

      if(fabs(PDG)==23 && status==62 && barcode<=200000) { // for MG signal sample
        Truth_Z.push_back(p4);
        //cout << "Nevt = " << m_eventCounter << ", MZ = " << p4.M()/1000. << endl;
      }

      if(fabs(PDG)==11 && status==1 && barcode<=200000) {
        bool IsParentHadron = CheckFromHadron(trPart, PDG);
        bool IsParentPhoton = false;
        if(!IsParentHadron && !IsParentPhoton) {
          if(PDG>0) p4elem.push_back(p4);
          if(PDG<0) p4elep.push_back(p4);
        }
      }

      if(fabs(PDG)==13 && status==1 && barcode<=200000) {
        bool IsParentHadron = CheckFromHadron(trPart, PDG);
        bool IsParentPhoton = false;
        if(!IsParentHadron && !IsParentPhoton) {
          if(PDG>0) p4muonm.push_back(p4);
          if(PDG<0) p4muonp.push_back(p4);
        }
      }

      if((fabs(PDG)==11 || fabs(PDG)==13) && status==1 && barcode<=200000) {
        bool IsParentHadron = CheckFromHadron(trPart, PDG);
        if(!IsParentHadron) {
          Truth_l.push_back(p4);
          Truth_pid.push_back(PDG);
          //cout << "Event: " << event;
          //cout << ", PDG=" << PDG << ", " << p4.Pt() << ", " << p4.Eta() << ", " << p4.Phi();
          //cout << ", status=" << status << ", barcode=" << barcode << endl;
        }
      }

      if(trPart->isTau() && (status==2 || status==10902)) {
        bool IsParentHadron = CheckFromHadron(trPart, PDG);
        bool IsParentPhoton = CheckFromPhoton(trPart,22);
        if(!IsParentHadron && !IsParentPhoton) {
          WithTau=true;
        }
      }

      if( fabs(PDG)==22 && status==1 && barcode<99999) {
        bool IsParentHadron = CheckFromHadron(trPart, PDG);
        if(!IsParentHadron) {
          vp4g.push_back(p4);
        }
      }

      // Truth quark
      // NOTE: Detailed selections will depend on generator setup
      // This is specified for mc15_13TeV.363998.PhPy8EG_NNPDF30NLO_AZNLOC6L1_QCD_WpWp.merge.DAOD_STDM3
      if((fabs(PDG)<=6 || fabs(PDG)==21) && status==23) {
        Truth_q.push_back(p4);
        //if(p4.Pt() > 10. && fabs(p4.Eta()) < 5.) {
        //  Truth_q.push_back(p4);
        //}
      }
    }
    //cout << "nQuark = " << Truth_q.size() << endl;

    int num_lep = (int)p4elep.size()+(int)p4elem.size()+(int)p4muonp.size()+(int)p4muonm.size();

    if(!WithTau && num_lep>=4) {
      m_filter++;
      m_nlep += num_lep;

      for(int i=0; i<(int)p4elep.size(); i++) {
        for(int j=0; j<(int)vp4g.size(); j++) {
          if(vp4g[j].DeltaR(p4elep[i])<0.1) {
            p4elep[i] = p4elep[i]+vp4g[j];
            vp4g.erase(vp4g.begin()+j);
            j--;
          }
        }
      }
      for(int i=0; i<(int)p4elem.size(); i++) {
        for(int j=0; j<(int)vp4g.size(); j++) {
          if(vp4g[j].DeltaR(p4elem[i])<0.1) {
            p4elem[i] = p4elem[i]+vp4g[j];
            vp4g.erase(vp4g.begin()+j);
            j--;
          }
        }
      }
      for(int i=0; i<(int)p4muonp.size(); i++) {
        for(int j=0; j<(int)vp4g.size(); j++) {
          if(vp4g[j].DeltaR(p4muonp[i])<0.1) {
            p4muonp[i] = p4muonp[i]+vp4g[j];
            vp4g.erase(vp4g.begin()+j);
            j--;
          }
        }
      }
      for(int i=0; i<(int)p4muonm.size(); i++) {
        for(int j=0; j<(int)vp4g.size(); j++) {
          if(vp4g[j].DeltaR(p4muonm[i])<0.1) {
            p4muonm[i] = p4muonm[i]+vp4g[j];
            vp4g.erase(vp4g.begin()+j);
            j--;
          }
        }
      }
      Quad Truth_Quad;
      bool passFid2e2m = CheckQuadTruth(p4muonp, p4muonm, p4elep, p4elem, CHN, "eemm", FLAG_cut_temp, Truth_Quad);
      if(passFid2e2m) m_2e2mfidu4 += 1.0;
      bool passFid4m = CheckQuadTruth(p4muonp, p4muonm, p4elep, p4elem, CHN, "mmmm", FLAG_cut_temp, Truth_Quad);
      if(passFid4m) m_4mfidu4 += 1.0;
      bool passFid4e = CheckQuadTruth(p4muonp, p4muonm, p4elep, p4elem, CHN, "eeee", FLAG_cut_temp, Truth_Quad);
      if(passFid4e) m_4efidu4 += 1.0;

      if(passFid2e2m || passFid4m || passFid4e) {
        MTZZ = (Truth_Quad.pair[0].Z+Truth_Quad.pair[1].Z).M(); 
        ptZZ_truth = (Truth_Quad.pair[0].Z+Truth_Quad.pair[1].Z).Pt(); 
        isTruth=true;
      }
    }   
  }

  //TruthJets
  vector<TLorentzVector> Truth_j;
  const xAOD::JetContainer* TruthJets = 0;
  if(isMC) {
    if( !m_event->retrieve( TruthJets, "AntiKt4TruthJets").isSuccess() ){
      Error("excute()", "Failed to retrieve Truth Jets info. Exiting." );
      return EL::StatusCode::FAILURE;
    }

    xAOD::JetContainer::const_iterator trJets_itr = TruthJets->begin();
    xAOD::JetContainer::const_iterator trJets_end = TruthJets->end();
    for( ; trJets_itr != trJets_end; ++trJets_itr ) {
      //if(((*trJets_itr)->pt()>25.e3 && fabs((*trJets_itr)->eta())<2.4) || ((*trJets_itr)->pt()>30.e3 && fabs((*trJets_itr)->eta())<4.5)) {}
      if((*trJets_itr)->pt()>20.e3 && fabs((*trJets_itr)->eta())<5.) {
        TLorentzVector tempj;
        tempj.SetPxPyPzE((*trJets_itr)->px(), (*trJets_itr)->py(), (*trJets_itr)->pz(), (*trJets_itr)->e());
        Truth_j.push_back(tempj);
      }
    }

    TLorentzVector tmp_jet;

    // sort PS level jets
    for (int i=0; i<(int)Truth_j.size()-1; i++) {
      for (int j=i+1; j<(int)Truth_j.size();j++){
        if(Truth_j[i].Pt()<Truth_j[j].Pt()){
          tmp_jet = Truth_j[i];
          Truth_j[i] = Truth_j[j];
          Truth_j[j] = tmp_jet;
        }
      }
    }
    //cout << "Jet Pt: ";
    //for(int i=0; i<Truth_j.size(); i++) cout << Truth_j[i].Pt() << ", ";
    //cout << endl;
  }


  const xAOD::Vertex *pv(0);
  const xAOD::VertexContainer* vertices = 0;
  if( !m_event->retrieve( vertices, "PrimaryVertices" ).isSuccess() )
  {
    Error("execute ()", "Failed to retrieve verteices. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  for ( const auto* const vtx_itr : *vertices )
  {
    if (vtx_itr->vertexType() != xAOD::VxType::VertexType::PriVtx) continue;
    else { pv = vtx_itr; break;}
  }
  if(!pv) return EL::StatusCode::SUCCESS;
  double pz0 = pv->z();

  //  cout << endl;
  //  cout << "start new event " << m_eventCounter << endl;

  for (auto sysListItr : m_sysList){
    string sysname = sysListItr.name();

    if(SETTING[SETNAME[0]]["dosys"]==0) {
      if(sysname != "") continue;
    }
    if(sysname == "") sysname = "NOMINAL";
    if(SETTING[SETNAME[0]]["docorr"]==0) sysname="NOCORR";
    // cout << "name sys is : " << sysListItr.name() << endl;
    if (quickAna->applySystematicVariation (sysListItr) == CP::SystematicCode::Ok) {


      ClearFlags(FLAG_cut_temp);
      ClearFlags(FLAG_cut);
      ClearWeight(Evt_Weight);

      ClearVariables(HistVar);
      ClearVariables(VVar);
      ClearVariables(Hist2DVar);
      ClearVariables(TreeIntVar);
      ClearVariables(TreeIntVVar);
      ClearVariables(TreeDouVar);
      ClearVariables(TreeLngVar);
      ClearVariables(TreeFltVar);
      ClearVariables(TreeFltVVar);
      ClearVariables(TreeTLVVar);
      ClearVariables(TreeTLVVVar);

      goodm.clear(); goodmiso.clear();
      goode.clear(); goodeiso.clear();

      quickAna->process(*m_event).ignore();

      // starting point

      //double pileWeight=1.0;
      auto evtInfo = quickAna->eventinfo();
      int RadNum = -1;
      if(SETTING[SETNAME[0]]["doweight"]==1 && isMC) {
        //pileWeight = evtInfo->auxdata<float>("PileupWeight");
        //RadNum = evtInfo->auxdata<unsigned int>("RandomRunNumber");
      }
      weight = pileWeight*gen_weight;

      if(evtInfo->auxdata<bool>("HLT_e24_lhmedium_L1EM18VH_passTrig")) {
        m_tag1++;
        if(evtInfo->auxdata<bool>("HLT_e24_lhmedium_L1EM20VH_passTrig")) m_tag2++;
      }

      //if(evtInfo->auxdata<bool>("HLT_e24_lhmedium_L1EM20VH_passTrig")) {
      //   m_tag1++;
      //   if(evtInfo->auxdata<bool>("HLT_e24_lhmedium_L1EM18VH_passTrig")) m_tag2++;
      //}

      SetFlag(FLAG_cut_temp,"xAOD","All", 1);
      SetWeight(Evt_Weight,"xAOD", "All", weight);

      bool passTrig = false;
      MapType_String::iterator it;
      for(it=TRIG_vec.begin(); it!=TRIG_vec.end(); it++) {
        string hltname = (*it).first;
        string status = (*it).second;

        if(isMC) {
          if(RadNum<290000) {
            if(status.find("2015")!=string::npos && status.find("MC")!=string::npos) {
              passTrig = passTrig || evtInfo->auxdata<bool>(hltname+"_passTrig");
            }
          }else if(RadNum>290000 && RadNum<302919) {
            if(RadNum==298687 && status.find("2016_298687")!=string::npos)
              passTrig = passTrig || evtInfo->auxdata<bool>(hltname+"_passTrig");

            if(status.find("2016")!=string::npos && status.find("MC")!=string::npos &&
                (status.find("ALL")!=string::npos || status.find("ABC")!=string::npos)) {
              passTrig = passTrig || evtInfo->auxdata<bool>(hltname+"_passTrig");
            }
          }else if(RadNum>=302919) {
            if(status.find("2016")!=string::npos && status.find("MC")!=string::npos &&
                (status.find("ALL")!=string::npos || status.find("D_")!=string::npos)) {
              passTrig = passTrig || evtInfo->auxdata<bool>(hltname+"_passTrig");
            }
          }
        }else {
          if(run<290000) {
            if(status.find("2015")!=string::npos && status.find("Data")!=string::npos) {
              passTrig = passTrig || evtInfo->auxdata<bool>(hltname+"_passTrig");
            }
          }else if(run>290000 && run<302919) {
            if(run==298687 && status.find("2016_298687")!=string::npos)
              passTrig = passTrig || evtInfo->auxdata<bool>(hltname+"_passTrig");

            if(status.find("2016")!=string::npos && status.find("Data")!=string::npos &&
                (status.find("ALL")!=string::npos || status.find("ABC")!=string::npos)) {
              passTrig = passTrig || evtInfo->auxdata<bool>(hltname+"_passTrig");
            }
          }else if(run>=302919) {
            if(status.find("2016")!=string::npos && status.find("Data")!=string::npos &&
                (status.find("ALL")!=string::npos || status.find("D_")!=string::npos)) {
              passTrig = passTrig || evtInfo->auxdata<bool>(hltname+"_passTrig");
            }
          }
        }
      }

      SetFlag(FLAG_cut_temp,"Trigger","All", passTrig);
      SetWeight(Evt_Weight, "Trigger", "All", weight);


      //if(!passTrig) continue;

      int index_muon = 0;
      int index_ele = 0;
      int index_jet = 0;
      TLorentzVector p4jet1, p4jet2;
      p4jet1.SetPtEtaPhiM(-999., 0., 0., 0.);
      p4jet2.SetPtEtaPhiM(-999., 0., 0., 0.);
      double leadJet=-999., secJet=-999.;

      for(auto muon : *quickAna->muons()) {

        // without correction
        //typedef ElementLink<xAOD::IParticleContainer> LinkType;
        //static const char* linkName = "originalObjectLink";
        //LinkType& auxLink = muon->auxdata<LinkType> (linkName);
        //const xAOD::Muon* origMuon = dynamic_cast<const xAOD::Muon*>(*auxLink.cptr());

        OBJ_MUON muonInfo;

        muonInfo.author = muon->author();
        muonInfo.charge = muon->charge();
        muonInfo.type =   muon->muonType();

        muonInfo.L.SetPtEtaPhiM( muon->pt(), muon->eta(), muon->phi(), m_mass ); 

        // isolation reconstruction is not good, kept for future

        float etcone20=0, ptcone20=0;
        if(muon->isolation(etcone20, xAOD::Iso::etcone20))
          muonInfo.etcone20 = etcone20;
        if(muon->isolation(ptcone20, xAOD::Iso::ptcone20))
          muonInfo.ptcone20 = ptcone20; 

        muon->auxdata< char >( "All" ) = true;

        if(muon->auxdata<char> ("ana_select_ID")) muon->auxdata< char >( "Tool" ) = true;
        if(muon->auxdata<char> ("ana_select_Pt")) muon->auxdata< char >( "Pt" ) = true;
        if(muon->auxdata<char> ("ana_select_Pt_Calo")) muon->auxdata< char >( "Pt_Calo" ) = true;
        if(muon->auxdata<char>("ana_select_D0")) muon->auxdata< char >( "D0" ) = true;
        if(muon->auxdata<char>("ana_select_Z0")) muon->auxdata< char >( "Z0" ) = true;
        if(muon->auxdata<char> ("ana_select")) muon->auxdata< char >( "OverLap" ) = true;

        bool passMuon=true;
        CountMuObj(muon, STEP_obj, CNT_obj, sysname, passMuon);
        if(passMuon) {
          muonInfo.trigM=false;

          MapType_String::iterator it;
          for(it=TRIG_vec.begin(); it!=TRIG_vec.end(); it++) {
            string hltname = (*it).first;
            string status = (*it).second;

            if(status.find("Muon")==string::npos) continue;

            if(isMC) {
              if(RadNum<290000) {
                if(status.find("2015")!=string::npos && status.find("MC")!=string::npos) {
                  muonInfo.trigM = (muonInfo.trigM || muon->auxdata<bool>(hltname+"_trigMatch"));
                }
              }else if(RadNum>290000 && RadNum<302919) {
                if(status.find("2016")!=string::npos && status.find("MC")!=string::npos &&
                    (status.find("ALL")!=string::npos || status.find("ABC")!=string::npos)) {
                  muonInfo.trigM = (muonInfo.trigM || muon->auxdata<bool>(hltname+"_trigMatch"));
                }
              }else if(RadNum>=302919) {
                if(status.find("2016")!=string::npos && status.find("MC")!=string::npos &&
                    (status.find("ALL")!=string::npos || status.find("D_")!=string::npos)) {
                  muonInfo.trigM = (muonInfo.trigM || muon->auxdata<bool>(hltname+"_trigMatch"));
                }
              }
            }else {
              if(run<290000) {
                if(status.find("2015")!=string::npos && status.find("Data")!=string::npos) {
                  muonInfo.trigM = (muonInfo.trigM || muon->auxdata<bool>(hltname+"_trigMatch"));
                }
              }else if(run>290000 && run<302919) {
                if(run==298687 && status.find("2016_298687")!=string::npos)
                  muonInfo.trigM = (muonInfo.trigM || muon->auxdata<bool>(hltname+"_trigMatch"));

                if(status.find("2016")!=string::npos && status.find("Data")!=string::npos &&
                    (status.find("ALL")!=string::npos || status.find("ABC")!=string::npos)) {
                  muonInfo.trigM = (muonInfo.trigM || muon->auxdata<bool>(hltname+"_trigMatch"));
                }
              }else if(run>=302919) {
                if(status.find("2016")!=string::npos && status.find("Data")!=string::npos &&
                    (status.find("ALL")!=string::npos || status.find("D_")!=string::npos)) {
                  muonInfo.trigM = (muonInfo.trigM || muon->auxdata<bool>(hltname+"_trigMatch"));
                }
              }
            }
          }

          muonInfo.sf=muon->auxdata<float>("ana_weight");
          muonInfo.index = index_muon;
          muonInfo.ptrMuon = muon;
          muon->auxdata<int>("index") = index_muon;
          goodm.push_back( muonInfo );

          index_muon++;
        }
      }    

      //cout << event << ":";
      for(auto electron : *quickAna->electrons()) {

        typedef ElementLink<xAOD::IParticleContainer> LinkType;
        static const char* linkName = "originalObjectLink";
        LinkType& auxLink = electron->auxdata<LinkType> (linkName);
        const xAOD::Electron* origEle = dynamic_cast<const xAOD::Electron*>(*auxLink.cptr());

        //cout << " Pt " << origEle->pt() << ", Eta " << origEle->eta() << ", Phi " << origEle->phi() << "; ";

        OBJ_ELECTRON eleInfo;

        eleInfo.author = electron->author();
        eleInfo.charge = electron->charge();

        eleInfo.L.SetPtEtaPhiM( electron->pt(), electron->eta(), electron->phi(), e_mass ); 

        float etcone20=0, ptcone20=0;
        if(electron->isolationValue(etcone20, xAOD::Iso::topoetcone20))
          eleInfo.etcone20 = etcone20;
        if(electron->isolationValue(ptcone20, xAOD::Iso::ptvarcone20))
          eleInfo.ptcone20 = ptcone20;



        uint8_t numberOfSCTHits = 0;
        electron->trackParticleSummaryValue(numberOfSCTHits, xAOD::numberOfSCTHits);
        eleInfo.Et = (numberOfSCTHits>=4) ? eleInfo.clE/cosh(eleInfo.trketa) : eleInfo.clE/cosh(eleInfo.cleta);

        eleInfo.passOQ = electron->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON);

        electron->auxdata< char >( "All" ) = true;

        if(electron->auxdata<char> ("ana_select_smzz4l_veryloose_Pt")) electron->auxdata< char >( "Pt" ) = true;
        if(electron->auxdata<char> ("ana_select_smzz4l_veryloose_Eta")) electron->auxdata< char >( "Eta" ) = true;
        if(electron->auxdata<char> ("ana_select_smzz4l_veryloose_OQ")) electron->auxdata< char >( "ObjQ" ) = true;
        if(electron->auxdata<char> ("ana_select_smzz4l_veryloose_selectionTool")) electron->auxdata< char >( "ID" ) = true;
        //if(electron->auxdata<char> ("ana_select_D0")) electron->auxdata< char >( "D0" ) = true;
        if(electron->auxdata<char> ("ana_select_smzz4l_veryloose_Z0")) electron->auxdata< char >( "Z0" ) = true;
        if(electron->auxdata<char> ("ana_select")) electron->auxdata< char >( "OverLap" ) = true;


        bool passEle=true;
        CountEleObj(electron, STEP_obj, CNT_obj, sysname, passEle);
        if(passEle) {

          eleInfo.ismedium = electron->auxdata<char> ("ana_select_hzhinv_medium_selectionTool");
          eleInfo.d0sig = electron->auxdata<double> ("d0Sig");
          eleInfo.d0 = electron->auxdata<double> ("d0value");
          eleInfo.z0 = electron->auxdata<double> ("z0value");
          //const xAOD::TrackParticle * IDtrkPart = electron->trackParticle();
          //eleInfo.vtxtype = -1;
          //if(IDtrkPart) {
          //  if(IDtrkPart->vertex()) {
          //    eleInfo.vtxtype = IDtrkPart->vertex()->vertexType();
          //    //cout << IDtrkPart->vertex()->chiSquared() << " " << IDtrkPart->vertex()->numberDoF() << endl;
          //  }                        
          //}


          eleInfo.trigM=false;

          MapType_String::iterator it;
          for(it=TRIG_vec.begin(); it!=TRIG_vec.end(); it++) {
            string hltname = (*it).first;
            string status = (*it).second;

            if(status.find("Electron")==string::npos) continue;

            if(isMC) {
              if(RadNum<290000) {
                if(status.find("2015")!=string::npos && status.find("MC")!=string::npos) {
                  eleInfo.trigM = (eleInfo.trigM || electron->auxdata<bool>(hltname+"_trigMatch"));
                }
              }else if(RadNum>290000 && RadNum<302919) {
                if(RadNum==298687 && status.find("2016_298687")!=string::npos)
                  eleInfo.trigM = (eleInfo.trigM || electron->auxdata<bool>(hltname+"_trigMatch"));

                if(status.find("2016")!=string::npos && status.find("MC")!=string::npos &&
                    (status.find("ALL")!=string::npos || status.find("ABC")!=string::npos)) {
                  eleInfo.trigM = (eleInfo.trigM || electron->auxdata<bool>(hltname+"_trigMatch"));
                }
              }else if(RadNum>=302919) {
                if(status.find("2016")!=string::npos && status.find("MC")!=string::npos &&
                    (status.find("ALL")!=string::npos || status.find("D_")!=string::npos)) {
                  eleInfo.trigM = (eleInfo.trigM || electron->auxdata<bool>(hltname+"_trigMatch"));
                }
              }
            }else {
              if(run<290000) {
                if(status.find("2015")!=string::npos && status.find("Data")!=string::npos) {
                  eleInfo.trigM = (eleInfo.trigM || electron->auxdata<bool>(hltname+"_trigMatch"));
                }
              }else if(run>290000 && run<302919) {
                if(run==298687 && status.find("2016_298687")!=string::npos)
                  eleInfo.trigM = (eleInfo.trigM || electron->auxdata<bool>(hltname+"_trigMatch"));

                if(status.find("2016")!=string::npos && status.find("Data")!=string::npos &&
                    (status.find("ALL")!=string::npos || status.find("ABC")!=string::npos)) {
                  eleInfo.trigM = (eleInfo.trigM || electron->auxdata<bool>(hltname+"_trigMatch"));
                }
              }else if(run>=302919) {
                if(status.find("2016")!=string::npos && status.find("Data")!=string::npos &&
                    (status.find("ALL")!=string::npos || status.find("D_")!=string::npos)) {
                  eleInfo.trigM = (eleInfo.trigM || electron->auxdata<bool>(hltname+"_trigMatch"));
                }
              }
            }
          }
          eleInfo.trigM = eleInfo.trigM && eleInfo.L.Pt()>25.e3;

          eleInfo.sf=electron->auxdata<float>("ana_weight_smzz4l");
          eleInfo.index = index_ele;
          eleInfo.ptrElectron = electron;
          electron->auxdata<int>("index") = index_ele;
          goode.push_back( eleInfo );

          index_ele++;
        }
      }
      //cout << endl;

      bool hasBadJet = false;
      for(auto jet : *quickAna->jets()) {

        typedef ElementLink<xAOD::IParticleContainer> LinkType;
        static const char* linkName = "originalObjectLink";
        LinkType& auxLink = jet->auxdata<LinkType> (linkName);
        //const xAOD::Jet* origJet = dynamic_cast<const xAOD::Jet*>(*auxLink.cptr());

        OBJ_JET jetInfo;
        jetInfo.pt = jet->pt();
        jetInfo.eta= jet->eta();
        jetInfo.L  = jet->p4();

        jet->auxdata< char >( "All" ) = true;

        if(fabs(jetInfo.eta) < 4.5) jet->auxdata< char >( "Eta") = true;

        bool pass_pteta = false;
        if(fabs(jetInfo.eta) < 2.4) jet->auxdata< char >( "PtEta") = jetInfo.pt/1000. > 30.;
        else if(fabs(jetInfo.eta) < 4.5) jet->auxdata< char >( "PtEta") = jetInfo.pt/1000. > 40.;
        //if(jet->auxdata<char>("ana_select_kin_select") ) jet->auxdata< char >( "PtEta") = true;

        if(jet->auxdata<char> ("ana_select")) jet->auxdata< char >( "OverLap") = true;

        //if(jetInfo.pt/1000. > 60.) jet->auxdata< char >( "JVT" ) = true;
        //else if(jet->auxdata<char>("ana_select_jvt")) jet->auxdata< char >( "JVT" ) = true;

        if(jet->auxdata<char>("ana_select_jvt")) jet->auxdata< char >( "JVT" ) = true;

        bool passJet=true;
        CountJetObj(jet, STEP_obj, CNT_obj, sysname, passJet);
        if(passJet) {
          jetInfo.index = index_jet;
          jet->auxdata<int>("index") = index_jet;

          //if(jet->auxdata<char>("ana_select_cleaning_tool"))
          //  DoCounting(sysname, CNT_obj, "jet", "Clean");
          //else {
          //  passJetCleaning=false;
          //  hasBadJet = true;
          //}

          goodj.push_back( jetInfo );
          index_jet++;

        }
      }

      int nMuons = goodm.size(); 
      int nElectrons = goode.size(); 
      int nJets = goodj.size();

      SetFlag(FLAG_cut_temp,"Cleaning","All", passJetCleaning==true);
      SetWeight(Evt_Weight, "Cleaning","All", weight);

      SetFlag(FLAG_cut_temp,"FourMore","eeee", nElectrons>=4);
      SetFlag(FLAG_cut_temp,"FourMore","mmmm", nMuons>=4);
      SetFlag(FLAG_cut_temp,"FourMore","eemm", nElectrons>=2 && nMuons>=2);
      SetFlag(FLAG_cut_temp,"FourMore","incl", nElectrons>=4 || nMuons>=4 || (nElectrons>=2 && nMuons>=2));
      SetWeight(Evt_Weight, "FourMore","All", weight);


      vector<Pair> pairs;
      vector<Quad> quads;

      for(int i=0; i<(int)goode.size(); i++) {
        for(int j=i+1; j<(int)goode.size(); j++) { // two electrons
          if(goode[i].charge*goode[j].charge != -1) continue;

          int index1, index2;
          if(goode[i].L.Pt()>=goode[j].L.Pt()) {index1=i; index2=j;}
          else {index1=j; index2=i;}

          Pair temp;
          temp.flavor=0;
          temp.index.push_back(index1); temp.index.push_back(index2);
          temp.lepton.push_back(goode[index1].L); temp.lepton.push_back(goode[index2].L);
          temp.lepton_trk.push_back(goode[index1].L_trk); temp.lepton_trk.push_back(goode[index2].L_trk);
          temp.charge.push_back((int)goode[index1].charge); temp.charge.push_back((int)goode[index2].charge);
          temp.ptrVElectron.push_back(goode[index1].ptrElectron); temp.ptrVElectron.push_back(goode[index2].ptrElectron);
          //temp.passIso.push_back(goode[index1].passIso); temp.passIso.push_back(goode[index2].passIso);
          temp.Z = goode[index1].L + goode[index2].L;
          pairs.push_back(temp);
        }
      }
      for(int i=0; i<(int)goodm.size(); i++) {
        for(int j=i+1; j<(int)goodm.size(); j++) {
          if(goodm[i].charge*goodm[j].charge != -1) continue;

          int index1, index2;
          if(goodm[i].L.Pt()>=goodm[j].L.Pt()) {index1=i; index2=j;}
          else {index1=j; index2=i;}

          Pair temp;
          temp.flavor=1;
          temp.index.push_back(index1); temp.index.push_back(index2);
          temp.lepton.push_back(goodm[index1].L); temp.lepton.push_back(goodm[index2].L);
          temp.lepton_trk.push_back(goodm[index1].L_id); temp.lepton_trk.push_back(goodm[index2].L_id);
          temp.charge.push_back((int)goodm[index1].charge); temp.charge.push_back((int)goodm[index2].charge);
          temp.ptrVMuon.push_back(goodm[index1].ptrMuon); temp.ptrVMuon.push_back(goodm[index2].ptrMuon);
          //temp.passIso.push_back(goodm[index1].passIso); temp.passIso.push_back(goodm[index2].passIso);
          temp.Z = temp.lepton[0] + temp.lepton[1];
          pairs.push_back(temp);
        }
      }


      bool vetopair=false;
      // veto events with any pair |m2l-mZ|<5GeV
      for(int i=0; i<(int)pairs.size();i++) {
        if( fabs(pairs[i].Z.M()-ZMass)<5.e3 ){
            vetopair=true;
            break;
        }
      }  

      SetFlag(FLAG_cut_temp,"VetoPair","All", vetopair==false);
      SetWeight(Evt_Weight, "VetoPair","All", weight);
      //if(vetopair) continue;

      for(int i=0; i<(int)pairs.size();i++) {
        for(int j=i+1; j<(int)pairs.size(); j++) {
          if(pairs[i].flavor==pairs[j].flavor)
            if(pairs[i].index[0]==pairs[j].index[0] || pairs[i].index[1]==pairs[j].index[0] || pairs[i].index[0]==pairs[j].index[1] || pairs[i].index[1]==pairs[j].index[1]) continue;

          Quad temp;

          temp.pair.push_back(pairs[i]);
          temp.index.push_back(i);
          temp.flavor_lep.push_back(pairs[i].flavor);
          temp.flavor_lep.push_back(pairs[i].flavor);
          temp.index_lep.push_back(pairs[i].index[0]);
          temp.index_lep.push_back(pairs[i].index[1]);
          temp.lepton.push_back(pairs[i].lepton[0]);
          temp.lepton.push_back(pairs[i].lepton[1]);
          temp.charge_lep.push_back(pairs[i].charge[0]);
          temp.charge_lep.push_back(pairs[i].charge[1]);

          temp.pair.push_back(pairs[j]);
          temp.index.push_back(j);
          temp.flavor_lep.push_back(pairs[j].flavor);
          temp.flavor_lep.push_back(pairs[j].flavor);
          temp.index_lep.push_back(pairs[j].index[0]);
          temp.index_lep.push_back(pairs[j].index[1]);
          temp.lepton.push_back(pairs[j].lepton[0]);
          temp.lepton.push_back(pairs[j].lepton[1]);
          temp.charge_lep.push_back(pairs[j].charge[0]);
          temp.charge_lep.push_back(pairs[j].charge[1]);

          temp.ZZ = pairs[i].Z + pairs[j].Z;
          if(pairs[i].flavor==0&&pairs[j].flavor==0) temp.type="eeee";
          if(pairs[i].flavor==1&&pairs[j].flavor==0) temp.type="eemm";
          if(pairs[i].flavor==0&&pairs[j].flavor==1) temp.type="eemm";
          if(pairs[i].flavor==1&&pairs[j].flavor==1) temp.type="mmmm";

          temp.weight = 1.0;
          for(int k=0;k<4;k++) {
            int ind=temp.index_lep[k];
            int fla=temp.flavor_lep[k];
            if(fla==0) temp.medium.push_back(goode[ind].ismedium?1:0);
            double w = (fla==0)? goode[ind].sf : goodm[ind].sf;
            if(SETTING[SETNAME[0]]["doweight"]==1 && isMC) temp.weight *= w;

            temp.etcone20.push_back(fla==0? goode[ind].etcone20 : goodm[ind].etcone20);
            temp.ptcone20.push_back(fla==0? goode[ind].ptcone20 : goodm[ind].ptcone20);
            temp.d0sig.push_back(fla==0? goode[ind].d0sig : goodm[ind].d0sig);
            temp.z0.push_back(fla==0? goode[ind].z0 : goodm[ind].z0);
            temp.d0.push_back(fla==0? goode[ind].d0 : goodm[ind].d0);
          }

          if(temp.pair[0].flavor==temp.pair[1].flavor) {
            if(temp.charge_lep[0] + temp.charge_lep[3] == 0) {
              temp.alter_pairs.push_back(temp.lepton[0]+temp.lepton[3]);
              temp.alter_pairs.push_back(temp.lepton[2]+temp.lepton[1]);
            }
            else {
              temp.alter_pairs.push_back(temp.lepton[0]+temp.lepton[2]);
              temp.alter_pairs.push_back(temp.lepton[1]+temp.lepton[3]);
            }
          }
          quads.push_back(temp);
        }
      }

      bool vetoquad=false;
      // veto events with any m4l+5GeV>mZ
      for(int i=0; i<(int)quads.size();i++) {
        if( quads[i].ZZ.M()+5.e3>ZMass ){
            vetoquad=true;
            break;
        }
      }  
      SetFlag(FLAG_cut_temp,"VetoQuad","All", vetoquad==false);
      SetWeight(Evt_Weight, "VetoQuad","All", weight);
      //if(vetoquad) continue;


      bool has4e=false, has4mu=false, has2e2mu=false;
      double Minpull=999.e3;
      Quad best_Quad;
      int best_i = -1, event_type=-99;
      for(int i=0; i<(int)quads.size(); i++) {

        SetFlag(FLAG_cut_temp,"Quad", quads[i].type, 1);
        SetFlag(FLAG_cut_temp,"Quad", "incl", 1);
        if((quads[i].type=="mmmm" && !has4mu) || (quads[i].type=="eeee" && !has4e) || (quads[i].type=="eemm" && !has2e2mu)) SetWeight(Evt_Weight,"Quad",quads[i].type, weight*quads[i].weight);

        if(quads[i].type=="mmmm") has4mu=true;
        else if(quads[i].type=="eemm") has2e2mu=true;
        else has4e=true;
      
        //bool pass_pt=false;
        //double lep_pt1=0., lep_pt2=0., lep_pt3=0., lep_pt4=0.;
        //for(int j=0; j<4; j++) {
        //  if(quads[i].lepton[j].Pt()>lep_pt1) lep_pt1=quads[i].lepton[j].Pt();
        //}
        //for(int j=0; j<4; j++) {
        //  if(quads[i].lepton[j].Pt()>lep_pt2 && quads[i].lepton[j].Pt()<lep_pt1)
        //    lep_pt2=quads[i].lepton[j].Pt();
        //}
        //for(int j=0; j<4; j++) {
        //  if(quads[i].lepton[j].Pt()>lep_pt3 && quads[i].lepton[j].Pt()<lep_pt2)
        //    lep_pt3=quads[i].lepton[j].Pt();
        //}
        //for(int j=0; j<4; j++) {
        //  if(quads[i].lepton[j].Pt()<lep_pt3)
        //    lep_pt4=quads[i].lepton[j].Pt();
        //}

        //pass_pt = lep_pt1>15.e3 && lep_pt2>8.e3 && lep_pt3>7.e3 && lep_pt4>5.e3;

        //if(!pass_pt) continue;

        //SetFlag(FLAG_cut_temp,"Kine", quads[i].type, 1);
        //SetFlag(FLAG_cut_temp,"Kine", "incl", 1);

        //int num_sa=0, num_ca=0;
        //for(int j=0; j<4; j++) {
        //  int ind=quads[i].index_lep[j];
        //  int fla=quads[i].flavor_lep[j];
        //  if(fla==0) continue;
        //  if(goodm[ind].type==1) num_sa++;
        //  if(goodm[ind].type==3) num_ca++;
        //}
        //if((num_sa+num_ca)>1) continue;

        //SetFlag(FLAG_cut_temp,"nCalo", quads[i].type, 1);
        //SetFlag(FLAG_cut_temp,"nCalo", "incl", 1);

        double massZ1 = quads[i].pair[0].Z.M();
        double massZ2 = quads[i].pair[1].Z.M();
        double pull = fabs(massZ1-massZ2);
        if(pull<Minpull) {
          Minpull = pull;
          best_Quad = quads[i];
          best_i = i;
          if(quads[i].type=="eeee")  event_type = quadType::_4e;
          if(quads[i].type=="eemm")  event_type = quadType::_2e2mu;
          if(quads[i].type=="mmmm")  event_type = quadType::_4mu;
        } 
        //        good_quads.push_back(quads[i]);
      }


      Quad Final_Quad;
      double sf_exp = 1.0;
      if(best_i!=-1) {

        bool passdR = false;
        double mindR1=99.0, mindR2=99.0;
        for(int j=0; j<4; j++) {
          for(int k=j+1; k<4; k++) {
            if(best_Quad.flavor_lep[j]==best_Quad.flavor_lep[k])
              if(best_Quad.lepton[j].DeltaR(best_Quad.lepton[k])<mindR1)
                mindR1=best_Quad.lepton[j].DeltaR(best_Quad.lepton[k]);

            if(best_Quad.flavor_lep[j]!=best_Quad.flavor_lep[k])
              if(best_Quad.lepton[j].DeltaR(best_Quad.lepton[k])<mindR2)
                mindR2=best_Quad.lepton[j].DeltaR(best_Quad.lepton[k]);
          }
        }

        passdR = mindR1>0.2 && mindR2>0.2;
        SetFlag(FLAG_cut_temp,"dR4l", best_Quad.type, passdR);
        SetFlag(FLAG_cut_temp,"dR4l", "incl", passdR);
        SetWeight(Evt_Weight,"dR4l",best_Quad.type, weight*best_Quad.weight);


        bool passJPsiVeto = true;
        if(best_Quad.pair[0].flavor==best_Quad.pair[1].flavor) {
          for(int j=0; j<(int)best_Quad.alter_pairs.size(); j++)
            if(best_Quad.alter_pairs[j].M()<5.e3) passJPsiVeto=false;
        }
        if(best_Quad.pair[0].Z.M()<5.e3) passJPsiVeto=false;
        if(best_Quad.pair[1].Z.M()<5.e3) passJPsiVeto=false;

        SetFlag(FLAG_cut_temp,"JPsiVeto", best_Quad.type, passJPsiVeto);
        SetFlag(FLAG_cut_temp,"JPsiVeto", "incl", passJPsiVeto);
        SetWeight(Evt_Weight,"JPsiVeto",best_Quad.type, weight*best_Quad.weight);

        bool pass_Iso=true;
        bool pass_d0sig=true;
        bool pass_Loose=true;
        for(int j=0; j<4; j++) {
          int ind=best_Quad.index_lep[j];
          int fla=best_Quad.flavor_lep[j];

          if(fla==1) {
            xAOD::Muon* muon = goodm[ind].ptrMuon;
            if(!muon->auxdata<char>("PassIso")) pass_Iso=false;
            // IP cut for non-SA muons
            if(muon->muonType() != xAOD::Muon::MuonType::MuonStandAlone) {
              auto const *Track = muon->primaryTrackParticle();
              double d0sig = xAOD::TrackingHelpers::d0significance(Track,eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY());
              if(fabs(d0sig) > 3.0) pass_d0sig=false;
            }
            //if(fabs(muon->auxdata<double>("d0Sig"))>3.0) pass_d0sig=false;          
          }
          if(fla==0) {
            xAOD::Electron* electron = goode[ind].ptrElectron;
            if(!electron->auxdata<char>("ana_select_smzz4l_selectionTool")) pass_Loose=false;
            if(!electron->auxdata<char>("PassIso")) pass_Iso=false;
            auto const *Track = electron->trackParticle();
            double d0sig = xAOD::TrackingHelpers::d0significance(Track,eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY());
            if(fabs(d0sig) > 5.0) pass_d0sig=false;
          }
        }

        if(pass_Loose) {
          SetFlag(FLAG_cut_temp,"EleLoose", best_Quad.type, 1);
          SetWeight(Evt_Weight, "EleLoose", best_Quad.type, weight*best_Quad.weight);
          SetFlag(FLAG_cut_temp,"EleLoose", "incl", 1);
          SetWeight(Evt_Weight, "EleLoose", "incl", weight*best_Quad.weight);
        }
        if(pass_Iso) { 
          SetFlag(FLAG_cut_temp,"Iso", best_Quad.type, 1);
          SetWeight(Evt_Weight, "Iso", best_Quad.type, weight*best_Quad.weight);
          SetFlag(FLAG_cut_temp,"Iso", "incl", 1);
          SetWeight(Evt_Weight, "Iso", "incl", weight*best_Quad.weight);
        }
        if(pass_d0sig) {
          SetFlag(FLAG_cut_temp,"IP", best_Quad.type, 1);
          SetWeight(Evt_Weight, "IP", best_Quad.type, weight*best_Quad.weight);
          SetFlag(FLAG_cut_temp,"IP", "incl", 1);
          SetWeight(Evt_Weight, "IP", "incl", weight*best_Quad.weight);
        }

        // m4l filter at 100 GeV for Sherpa 4l sample
        SetFlag(FLAG_cut_temp,"Filter", "All", 1);
        SetWeight(Evt_Weight, "Filter", "All", weight*best_Quad.weight);
      }

      CountEvt(sysname, CHN, STEP_cut, FLAG_cut_temp, FLAG_cut, CNT_cut, Evt_Weight);

      // Histograms
      // Event level
      HistVar["mu"]["Value"] = ave_mu;

      if(FLAG_cut["eeee"]["JPsiVeto"] || FLAG_cut["eemm"]["JPsiVeto"] || FLAG_cut["mmmm"]["JPsiVeto"]) {
        // Lepton
        HistVar["PtL1"]["Value"] = best_Quad.pair[0].lepton[0].Pt()/1000.;
        HistVar["PtL2"]["Value"] = best_Quad.pair[0].lepton[1].Pt()/1000.;
        HistVar["PtL3"]["Value"] = best_Quad.pair[1].lepton[0].Pt()/1000.;
        HistVar["PtL4"]["Value"] = best_Quad.pair[1].lepton[1].Pt()/1000.;
        HistVar["EtaL1"]["Value"] = best_Quad.pair[0].lepton[0].Eta();
        HistVar["EtaL2"]["Value"] = best_Quad.pair[0].lepton[1].Eta();
        HistVar["EtaL3"]["Value"] = best_Quad.pair[1].lepton[0].Eta();
        HistVar["EtaL4"]["Value"] = best_Quad.pair[1].lepton[1].Eta();
        HistVar["PhiL1"]["Value"] = best_Quad.pair[0].lepton[0].Phi();
        HistVar["PhiL2"]["Value"] = best_Quad.pair[0].lepton[1].Phi();
        HistVar["PhiL3"]["Value"] = best_Quad.pair[1].lepton[0].Phi();
        HistVar["PhiL4"]["Value"] = best_Quad.pair[1].lepton[1].Phi();

        // Di-lepton
        HistVar["MZ1"]["Value"] = best_Quad.pair[0].Z.M()/1000.;
        HistVar["PtZ1"]["Value"] = best_Quad.pair[0].Z.Pt()/1000.;
        HistVar["MZ2"]["Value"] = best_Quad.pair[1].Z.M()/1000.;
        HistVar["PtZ2"]["Value"] = best_Quad.pair[1].Z.Pt()/1000.;

        // 4lepton
        HistVar["MZZ"]["Value"] = best_Quad.ZZ.M()/1000.;
        HistVar["PtZZ"]["Value"] = best_Quad.ZZ.Pt()/1000.;
      }

      // Fill hist
      FillHistograms(sysname);

      // Fill tree
      bool fillReco = FLAG_cut["eeee"]["JPsiVeto"] || FLAG_cut["eemm"]["JPsiVeto"] || FLAG_cut["mmmm"]["JPsiVeto"];
      //bool fillReco = FLAG_cut["eeee"]["xAOD"] || FLAG_cut["eemm"]["xAOD"] || FLAG_cut["mmmm"]["xAOD"];
      double m4l=-9999.e3, ptZZ=-9999.e3, m2l1=-9999.e3, m2l2=-9999.e3, ptZ=-9999.e3;
      weight = 1.0;

      if(fillReco) { 
        Final_Quad = best_Quad;
        weight = pileWeight*gen_weight*Final_Quad.weight;
        m4l = (Final_Quad.pair[0].Z+Final_Quad.pair[1].Z).M();
        ptZZ = (Final_Quad.pair[0].Z+Final_Quad.pair[1].Z).Pt();
        
        // MiniTree
        // Event level
        TreeIntVar["run"]["Value"] = run;
        TreeLngVar["event"]["Value"] = event;
        TreeDouVar["weight"]["Value"] = weight;

        // Lepton
        TreeIntVar["nLep"]["Value"] = nMuons + nElectrons;
        TreeDouVar["PtL1"]["Value"] = best_Quad.pair[0].lepton[0].Pt()/1000.;
        TreeDouVar["PtL2"]["Value"] = best_Quad.pair[0].lepton[1].Pt()/1000.;
        TreeDouVar["PtL3"]["Value"] = best_Quad.pair[1].lepton[0].Pt()/1000.;
        TreeDouVar["PtL4"]["Value"] = best_Quad.pair[1].lepton[1].Pt()/1000.;
        TreeDouVar["EtaL1"]["Value"] = best_Quad.pair[0].lepton[0].Eta();
        TreeDouVar["EtaL2"]["Value"] = best_Quad.pair[0].lepton[1].Eta();
        TreeDouVar["EtaL3"]["Value"] = best_Quad.pair[1].lepton[0].Eta();
        TreeDouVar["EtaL4"]["Value"] = best_Quad.pair[1].lepton[1].Eta();
        TreeDouVar["PhiL1"]["Value"] = best_Quad.pair[0].lepton[0].Phi();
        TreeDouVar["PhiL2"]["Value"] = best_Quad.pair[0].lepton[1].Phi();
        TreeDouVar["PhiL3"]["Value"] = best_Quad.pair[1].lepton[0].Phi();
        TreeDouVar["PhiL4"]["Value"] = best_Quad.pair[1].lepton[1].Phi();
        int pid_tmp[4] = {11, 13, -11, -13};
        TreeIntVar["PIDL1"]["Value"] = pid_tmp[best_Quad.flavor_lep[0]+best_Quad.charge_lep[0]+1];
        TreeIntVar["PIDL2"]["Value"] = pid_tmp[best_Quad.flavor_lep[1]+best_Quad.charge_lep[1]+1];
        TreeIntVar["PIDL3"]["Value"] = pid_tmp[best_Quad.flavor_lep[2]+best_Quad.charge_lep[2]+1];
        TreeIntVar["PIDL4"]["Value"] = pid_tmp[best_Quad.flavor_lep[3]+best_Quad.charge_lep[3]+1];

        // All good leptons
        for(int i=0; i<goode.size(); i++) {
          TreeTLVVVar["v_tlv_L"]["Value"].push_back((goode[i].L));
          if(goode[i].charge < 0) TreeIntVVar["v_PID"]["Value"].push_back(11);
          else TreeIntVVar["v_PID"]["Value"].push_back(-11);
        }
        for(int i=0; i<goodm.size(); i++) {
          TreeTLVVVar["v_tlv_L"]["Value"].push_back((goodm[i].L));
          if(goodm[i].charge < 0) TreeIntVVar["v_PID"]["Value"].push_back(13);
          else TreeIntVVar["v_PID"]["Value"].push_back(-13);
        }
        for(int i=0; i<goodj.size(); i++) {
          TreeTLVVVar["v_tlv_J"]["Value"].push_back((goodj[i].L));
        }

        // More variables for BDT
        TreeDouVar["MZ1"]["Value"] = best_Quad.pair[0].Z.M()/1000.;
        TreeDouVar["PtZ1"]["Value"] = best_Quad.pair[0].Z.Pt()/1000.;
        TreeDouVar["MZ2"]["Value"] = best_Quad.pair[1].Z.M()/1000.;
        TreeDouVar["PtZ2"]["Value"] = best_Quad.pair[1].Z.Pt()/1000.;
        TreeDouVar["MZZ"]["Value"] = best_Quad.ZZ.M()/1000.;
        TreeDouVar["PtZZ"]["Value"] = best_Quad.ZZ.Pt()/1000.;

        // Triggers
        TreeIntVar["HLT_e26_lhtight_nod0_ivarloose"]["Value"]                            = evtInfo->auxdata<bool>("HLT_e26_lhtight_nod0_ivarloose_passTrig");
        TreeIntVar["HLT_e60_lhmedium_nod0"]["Value"]                                     = evtInfo->auxdata<bool>("HLT_e60_lhmedium_nod0_passTrig");
        TreeIntVar["HLT_e140_lhloose_nod0"]["Value"]                                     = evtInfo->auxdata<bool>("HLT_e140_lhloose_nod0_passTrig");
        TreeIntVar["HLT_e60_lhmedium_nod0"]["Value"]                                     = evtInfo->auxdata<bool>("HLT_e60_lhmedium_nod0_passTrig");
        TreeIntVar["HLT_e140_lhloose_nod0"]["Value"]                                     = evtInfo->auxdata<bool>("HLT_e140_lhloose_nod0_passTrig");
        TreeIntVar["HLT_e300_etcut"]["Value"]                                            = evtInfo->auxdata<bool>("HLT_e300_etcut_passTrig");
        TreeIntVar["HLT_2e17_lhvloose_nod0_L12EM15VHI"]["Value"]                         = evtInfo->auxdata<bool>("HLT_2e17_lhvloose_nod0_L12EM15VHI_passTrig");
        TreeIntVar["HLT_2e24_lhvloose_nod0"]["Value"]                                    = evtInfo->auxdata<bool>("HLT_2e24_lhvloose_nod0_passTrig");
        TreeIntVar["HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH"]["Value"] = evtInfo->auxdata<bool>("HLT_e24_lhvloose_nod0_2e12_lhvloose_nod0_L1EM20VH_3EM10VH_passTrig");
        TreeIntVar["HLT_mu26_ivarmedium"]["Value"]                                       = evtInfo->auxdata<bool>("HLT_mu26_ivarmedium_passTrig");
        TreeIntVar["HLT_mu50"]["Value"]                                                  = evtInfo->auxdata<bool>("HLT_mu50_passTrig");
        TreeIntVar["HLT_mu60_0eta105_msonly"]["Value"]                                   = evtInfo->auxdata<bool>("HLT_mu60_0eta105_msonly_passTrig");
        TreeIntVar["HLT_2mu14"]["Value"]                                                 = evtInfo->auxdata<bool>("HLT_2mu14_passTrig");
        TreeIntVar["HLT_mu22_mu8noL1"]["Value"]                                          = evtInfo->auxdata<bool>("HLT_mu22_mu8noL1_passTrig");
        TreeIntVar["HLT_mu22_mu8noL1_calotag_0eta010"]["Value"]                          = evtInfo->auxdata<bool>("HLT_mu22_mu8noL1_calotag_0eta010_passTrig");
        TreeIntVar["HLT_mu20_2mu4noL1"]["Value"]                                         = evtInfo->auxdata<bool>("HLT_mu20_2mu4noL1_passTrig");
        TreeIntVar["HLT_3mu6"]["Value"]                                                  = evtInfo->auxdata<bool>("HLT_3mu6_passTrig");
        TreeIntVar["HLT_3mu6_msonly"]["Value"]                                           = evtInfo->auxdata<bool>("HLT_3mu6_msonly_passTrig");
        TreeIntVar["HLT_4mu4"]["Value"]                                                  = evtInfo->auxdata<bool>("HLT_4mu4_passTrig");
        TreeIntVar["HLT_e17_lhloose_nod0_mu14"]["Value"]                                 = evtInfo->auxdata<bool>("HLT_e17_lhloose_nod0_mu14_passTrig");
        TreeIntVar["HLT_e26_lhmedium_nod0_mu8noL1"]["Value"]                             = evtInfo->auxdata<bool>("HLT_e26_lhmedium_nod0_mu8noL1_passTrig");
        TreeIntVar["HLT_e7_lhmedium_nod0_mu24"]["Value"]                                 = evtInfo->auxdata<bool>("HLT_e7_lhmedium_nod0_mu24_passTrig");
        TreeIntVar["HLT_e12_lhloose_nod0_2mu10"]["Value"]                                = evtInfo->auxdata<bool>("HLT_e12_lhloose_nod0_2mu10_passTrig");
        TreeIntVar["HLT_2e12_lhloose_nod0_mu10"]["Value"]                                = evtInfo->auxdata<bool>("HLT_2e12_lhloose_nod0_mu10_passTrig");

      }

      if(fillReco) Tree[sysname]->Fill();   
    }

  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: finalize ()
{

  delete quickAna; 
  //  delete configTool;
  //  delete trigDecTool;

  system("rm -rf cutflow"); 
  system("mkdir -p cutflow"); 


  TFile *file3 = wk()->getOutputFile ("cutflow");
  file3->ReOpen("Update");

  //<Cut Flow Output>
  for(int i=0; i<(int)CHN.size(); i++) {
    string chn=CHN[i];
    string filename = "log_eff_" + chn + "_physics.txt";
    ofstream file(filename.c_str());

    if(file.is_open()) {
      for(int n=0; n<(int)SYSNAME.size(); n++) {
        string sys = SYSNAME[n];

        if(sys!="NOCORR" && (!SETTING["physics"]["docorr"])) continue;

        if(sys!="NOMINAL" && SETTING["physics"]["docorr"] && (!SETTING["physics"]["dosys"])) continue;

        if(sys=="NOCORR" && SETTING["physics"]["docorr"] && SETTING["physics"]["dosys"]) continue;


        file << "### <Channel : " << chn << "> ; Systematic type is : " 
          << sys << " ###" << endl;
        file <<endl;
        file << "WeightSum = " << sumOfWeights << endl;


        for(int j=0; j<(int)STEP_cut.size(); j++) {
          string cut = STEP_cut[j];
          file<<cut<<" = "<<CNT_cut[sys][chn][cut].num<<
            " ; "<<CNT_cut[sys][chn][cut].wnum<<
            " +/- "<<CNT_cut[sys][chn][cut].err<<endl;
        }
        file << endl;
        file << " ======== border ======== " << endl;
        file << endl;
      }
    }
    else {
      cout<<"Can not open file "<<filename<<endl;
      exit(-1);
    }
    file.close();

    TMacro* m = new TMacro(filename.c_str());
    m->Write();
    delete m;

    string command1 = "mv " + filename + " cutflow";
    system(command1.c_str());
  }

  string filename_obj = "log_eff_obj_physics.txt";
  ofstream file_obj(filename_obj.c_str());
  for(int n=0; n<(int)SYSNAME.size(); n++) {
    string sys = SYSNAME[n];

    if(sys!="NOCORR" && (!SETTING["physics"]["docorr"])) continue;

    if(sys!="NOMINAL" && SETTING["physics"]["docorr"] && (!SETTING["physics"]["dosys"])) continue;

    if(sys=="NOCORR" && SETTING["physics"]["docorr"] && SETTING["physics"]["dosys"]) continue;

    file_obj << "### <Object Selection> Systematic type is : " 
      << sys << " ###" << endl;
    file_obj <<endl;

    MapType_VString::iterator it;
    for(it=STEP_obj.begin(); it!=STEP_obj.end(); it++) {
      string obj=(*it).first;
      for(int i=0; i<(int)STEP_obj[obj].size(); i++) {
        string cut=STEP_obj[obj][i];
        string cutname=obj+"_"+cut;
        file_obj<<cutname<<" = "<<CNT_obj[sys][obj][cut].num
          <<" ; "<<CNT_obj[sys][obj][cut].wnum
          <<" +/- "<<CNT_obj[sys][obj][cut].err<<endl;
      }
      file_obj<< endl;
    }
    file_obj << endl;
    file_obj << " ======== border ======== " << endl;
    file_obj << endl;

  }  
  file_obj.close();


  TMacro* m = new TMacro(filename_obj.c_str());
  m->Write();
  delete m;

  string command2= "mv " + filename_obj + " cutflow";
  system(command2.c_str());

  cout << endl;
  cout << "####" << endl;
  //double A_2e2m = double(m_2e2mfidu4)/double(m_2e2mfidu0);
  //double A_4e = double(m_4efidu4)/double(m_4efidu0);
  //double A_4m = double(m_4mfidu4)/double(m_4mfidu0);
  printf("Finalize : %i total events have been processed !\n", m_eventCounter);
  printf("Finalize : %i events rejecting tau   !\n", m_filter);
  printf("Finalize : %i leptons   !\n\n", m_nlep);

  printf(" >>> Finalize : %f 4mu    events have been processed !\n", m_4mfilter);
  printf("Finalize : %f events fiducial0 at truth  !\n", m_4mfidu0);
  printf("Finalize : %f events fiducial1 at truth  !\n", m_4mfidu1);
  printf("Finalize : %f events fiducial2 at truth  !\n", m_4mfidu2);
  printf("Finalize : %f events fiducial3 at truth  !\n", m_4mfidu3);
  printf("Finalize : %f events fiducial4 at truth  !\n", m_4mfidu4);
  //  printf("Final A factor : %f at truth  !\n", A_4m);
  printf(" >>> Finalize : %f 4e     events have been processed !\n", m_4efilter);
  printf("Finalize : %f events fiducial0 at truth  !\n", m_4efidu0);
  printf("Finalize : %f events fiducial1 at truth  !\n", m_4efidu1);
  printf("Finalize : %f events fiducial2 at truth  !\n", m_4efidu2);
  printf("Finalize : %f events fiducial3 at truth  !\n", m_4efidu3);
  printf("Finalize : %f events fiducial4 at truth  !\n", m_4efidu4);
  //  printf("Final A factor : %f at truth  !\n", A_4e);
  printf(" >>> Finalize : %f 2e2mu  events have been processed !\n", m_2e2mfilter);
  printf("Finalize : %f events fiducial0 at truth  !\n", m_2e2mfidu0);
  printf("Finalize : %f events fiducial1 at truth  !\n", m_2e2mfidu1);
  printf("Finalize : %f events fiducial2 at truth  !\n", m_2e2mfidu2);
  printf("Finalize : %f events fiducial3 at truth  !\n", m_2e2mfidu3);
  printf("Finalize : %f events fiducial4 at truth  !\n", m_2e2mfidu4);
  //  printf("Final A factor : %f at truth  !\n", A_2e2m);
  printf("Finalize MyxAODAnalysis !");
  cout << "####" << endl;
  cout << endl;

  printf("! Trigger information !");
  printf("! %i events pass HLT_e24_lhmedium_L1EM18VH !\n", m_tag1);
  printf("! %i events pass HLT_e24_lhmedium_L1EM20VH !\n", m_tag2);


  time(&end);
  double timedif = difftime(end,start);
  if(timedif>3600) { printf("Finalize : Time Cost: %f hours\n", timedif/3600.); }
  else if(timedif>60) { printf("Finalize : Time Cost: %f minutes\n", timedif/60.); }
  else { printf("Finalize : Time Cost: %f second\n", timedif); }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histFinalize ()
{
  return EL::StatusCode::SUCCESS;
}
