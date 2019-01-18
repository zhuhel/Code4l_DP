#include <string>
#include <vector>
#include <math.h>
#include <iostream>

using namespace std;

#include "xAODMuon/Muon.h"
#include "xAODEgamma/Electron.h"
#include "MyAnalysis/AnalysisVar.h"

void SetWeight(MapType2_Double& map, string cut, string chn, double value) {

  if(chn=="All" || chn=="ALL" || chn=="all") {
    MapType2_Double::iterator it;
    for(it=map.begin(); it!=map.end(); it++) {
      string chn1=(*it).first;
      map[chn1][cut] *= value;
    }
  }
  else {
    map[chn][cut] *= value;
  }
}


void SetFlag(MapType2_Int& map, string cut, string chn, int value) {

  if(chn=="All" || chn=="ALL" || chn=="all") {
    MapType2_Int::iterator it;
    for(it=map.begin(); it!=map.end(); it++) {
      string chn1=(*it).first;
      map[chn1][cut]=value;
    }
  }
  else {
    map[chn][cut]=value;
  }
}


void DoCounting(string sysname, MapType3_Counting& CNT, string chn, string cut, double w=1.0) {
  double num = 1.0;
  double wnum = w*1.0;
  CNT[sysname][chn][cut].num += num;
  CNT[sysname][chn][cut].wnum += wnum;
  CNT[sysname][chn][cut].err = sqrt(pow(CNT[sysname][chn][cut].err,2)+wnum);
}



void CountMuObj(xAOD::Muon* muon, MapType_VString STEP_obj, MapType3_Counting& CNT, string sysname, bool& passAll) {

  vector<string> objstr = STEP_obj["mu"];

  for(int i=0; i<(int)objstr.size(); i++) {
    //if(objstr[i].find("OverLap") != string::npos) break;
    passAll = passAll && muon->auxdata< char >( objstr[i].c_str() );
    if(passAll) DoCounting(sysname, CNT, "mu", objstr[i]); 
  }
}

void CountMuIso(xAOD::Muon* muon, MapType_VString STEP_obj, MapType3_Counting& CNT, string sysname, bool& passAll) {

  vector<string> objstr = STEP_obj["mu2"];

  for(int i=0; i<(int)objstr.size(); i++) {
    //if(objstr[i].find("OverLap") != string::npos) break;
    passAll = passAll && muon->auxdata< char >( objstr[i].c_str() );
    if(passAll) DoCounting(sysname, CNT, "mu2", objstr[i]);
  }
}


void CountEleObj(xAOD::Electron* electron, MapType_VString STEP_obj, MapType3_Counting& CNT, string sysname, bool& passAll) {

  vector<string> objstr = STEP_obj["ele"];

  for(int i=0; i<(int)objstr.size(); i++) {
    //if(objstr[i].find("OverLap") != string::npos) break;
    passAll = passAll && electron->auxdata< char >( objstr[i].c_str() );
    if(passAll) DoCounting(sysname, CNT, "ele", objstr[i]);
  }
}

void CountEleIso(xAOD::Electron* electron, MapType_VString STEP_obj, MapType3_Counting& CNT, string sysname, bool& passAll) {

  vector<string> objstr = STEP_obj["ele2"];

  for(int i=0; i<(int)objstr.size(); i++) {
    //if(objstr[i].find("OverLap") != string::npos) break;
    passAll = passAll && electron->auxdata< char >( objstr[i].c_str() );
    if(passAll) DoCounting(sysname, CNT, "ele2", objstr[i]);
  }
}


void CountJetObj(xAOD::Jet* jet, MapType_VString STEP_obj, MapType3_Counting& CNT, string sysname, bool& passAll) {

  vector<string> objstr = STEP_obj["jet"];
  for(int i=0; i<(int)objstr.size(); i++) {
    if(objstr[i].find("Clean") != string::npos) break;
    passAll = passAll && jet->auxdata< char >( objstr[i].c_str() );
    if(passAll) DoCounting(sysname, CNT, "jet", objstr[i]);
  }
}

void CountD0Out(xAOD::Muon* muon, MapType_VString STEP_obj, MapType3_Counting& CNT, string sysname, bool& passAll) {
  vector<string> objstr = STEP_obj["mufd0"];

  for(int i=0; i<(int)objstr.size(); i++) {
    passAll = passAll && muon->auxdata< char >( objstr[i].c_str() )==1 ;
    if(passAll) DoCounting(sysname, CNT, "mufd0", objstr[i]);
  }
}

void CountIsoOut(xAOD::Muon* muon, MapType_VString STEP_obj, MapType3_Counting& CNT, string sysname, bool& passAll) {
  vector<string> objstr = STEP_obj["mufiso"];

  for(int i=0; i<(int)objstr.size(); i++) {
    passAll = passAll && muon->auxdata< char >( objstr[i].c_str() )==1 ;
    if(passAll) DoCounting(sysname, CNT, "mufiso", objstr[i]);
  }
}


void CountEvt(string sysname, vector<string> CHN, vector<string> STEP_cut, MapType2_Int& FLAG_cut_temp, MapType2_Int& FLAG_cut, MapType3_Counting& CNT, MapType2_Double& Evt_Weight) {

  MapType_VString::iterator it;
  bool passAll;
  for(int i=0; i<(int)CHN.size(); i++) {
    string chl =CHN[i];
    passAll=true;
    for(int j=0; j<(int)STEP_cut.size(); j++) {
      string cut=STEP_cut[j];
      passAll = passAll && (bool)FLAG_cut_temp[chl][cut];
      if(passAll) {
        FLAG_cut[chl][cut]=1;
        double w = Evt_Weight[chl][cut];
        DoCounting(sysname, CNT, chl, cut, w);
      }
    }
  }
}




bool FindZ(const xAOD::TruthParticle* trhPart, int dauPDG) {

  if(!trhPart) return false;

  int PDG = trhPart->pdgId();
  int status = trhPart->status();
  if(PDG!=dauPDG && fabs(PDG)!=23) return false;

  bool result = false;
  bool ProdVx = trhPart->hasProdVtx();
  if(!ProdVx) return false;
  const xAOD::TruthParticle* parPart = trhPart->parent();

  if(fabs(PDG)==23&&status==62) result=true;
  else result=FindZ(parPart, dauPDG);

  return result;
}

bool SameChild(const xAOD::TruthParticle* trPart) {

  int PDG = trPart->pdgId();
  int nChildren = trPart->nChildren();

  bool result=false;
  for(int i=0; i<nChildren; i++) {
    const xAOD::TruthParticle* dauPart = trPart->child(i);
    if(!dauPart) continue;
    int PDGdau = dauPart->pdgId();

    if(PDGdau==PDG) result=true;
  }

  return result;
}

bool CheckTruthOri(const xAOD::TruthParticle* trPart, int PDG) {
  const xAOD::TruthParticle * parent = trPart->parent();
  if (!parent){return false;}
  if (parent->pdgId()==PDG){return true;}
  return CheckTruthOri(parent, PDG); 
}

bool CheckFromPhoton(const xAOD::TruthParticle* trPart, int PDG) {
  const xAOD::TruthParticle * parent = trPart->parent();
  if (!parent){return false;}
  if (parent->pdgId()==PDG){return true;}
  return CheckFromPhoton(parent, PDG); 
}

bool CheckFromHadron(const xAOD::TruthParticle* trPart, int PDG) {

  //if(trPart->pdgId()==PDG) {
  //  bool ProdVx = trPart->hasProdVtx();
  //  if(!ProdVx) return false;
  //  
  //  CheckFromHadron(trPart->parent(), PDG);
  //}
  //else {
  //
  if(!trPart->parent()) return false;

  //if(trPart->parent()->pdgId() == PDG) return CheckFromHadron(trPart->parent(), PDG);

  bool const parent_NOT_in_initial_state = (trPart->parent()->status() != 4 );
  bool const parent_NOT_in_hard_process  = (trPart->parent()->status() != 11);

  if(fabs(PDG)==15) {
    if(trPart->parent()->isHadron() && 
        parent_NOT_in_initial_state && parent_NOT_in_hard_process) return true;
    //if(trPart->parent()->isHadron()) return true;
    //else return false;

    return CheckFromHadron(trPart->parent(), PDG);
  }else {
    if(trPart->parent()->isHadron() &&
        parent_NOT_in_initial_state && parent_NOT_in_hard_process) return true;
    //if(trPart->parent()->isHadron()) return true;
    //else return false;

    return CheckFromHadron(trPart->parent(), PDG);
  }
  //}
  //if(result) return result;
  //else {

  //  if(trPart->pdgId()==23) return false;

  //  bool ProdVx = trPart->hasProdVtx();
  //  if(!ProdVx) return result;
  //  const xAOD::TruthParticle* parPart = trPart->parent();

  //  result=CheckFromHadron(parPart);

  //  return result;
  //}
}

bool CheckQuadTruth(vector<TLorentzVector> p4muonp, vector<TLorentzVector> p4muonm,
    vector<TLorentzVector> p4elep,  vector<TLorentzVector> p4elem,
    vector<string> CHN, string type, MapType2_Int& FLAG_cut_temp,
    Quad& Final_Quad) {

  if( (p4muonp.size()+p4muonm.size()+p4elep.size()+p4elem.size()) <4) return false;

  vector<Pair> pairs;
  vector<Quad> quads,good_quads;

  for(int i=0; i<(int)p4elep.size(); i++) {
    if(p4elep[i].Pt()<6.e3) continue;
    if(fabs(p4elep[i].Eta())>2.7) continue;

    for(int j=0; j<(int)p4elem.size(); j++) { // two electrons
      if(p4elem[j].Pt()<6.e3) continue;
      if(fabs(p4elem[j].Eta())>2.7) continue;

      Pair temp;
      temp.flavor=0;
      temp.index.push_back(i); temp.index.push_back(j);
      temp.lepton.push_back(p4elep[i]); temp.lepton.push_back(p4elem[j]);
      temp.Z = p4elep[i] + p4elem[j];
      pairs.push_back(temp);
    }     
  }
  for(int i=0; i<(int)p4muonp.size(); i++) {
    if(p4muonp[i].Pt()<6.e3) continue;
    if(fabs(p4muonp[i].Eta())>2.7) continue;

    for(int j=0; j<(int)p4muonm.size(); j++) { // two muons
      if(p4muonm[j].Pt()<6.e3) continue;
      if(fabs(p4muonm[j].Eta())>2.7) continue;

      Pair temp;
      temp.flavor=1;
      temp.index.push_back(i); temp.index.push_back(j);
      temp.lepton.push_back(p4muonp[i]); temp.lepton.push_back(p4muonm[j]);
      temp.Z = p4muonp[i] + p4muonm[j];
      pairs.push_back(temp);
    } 
  }

  Quad best_Quad;
  double Maxipull=999.e3;
  int best_i = -1, event_type=-99;
  int tagDr=0, tagJpsi=0, tagQuad=0;
  for(int i=0; i<(int)pairs.size();i++) {
    for(int j=i+1; j<(int)pairs.size(); j++) {
      if(pairs[i].flavor==pairs[j].flavor)
        if(pairs[i].index[0]==pairs[j].index[0] || pairs[i].index[1]==pairs[j].index[1]) continue;

      Quad temp;

      temp.pair.push_back(pairs[i]);
      temp.index.push_back(i);
      temp.lepton.push_back(pairs[i].lepton[0]);
      temp.lepton.push_back(pairs[i].lepton[1]);

      temp.pair.push_back(pairs[j]);
      temp.index.push_back(j);
      temp.lepton.push_back(pairs[j].lepton[0]);
      temp.lepton.push_back(pairs[j].lepton[1]);

      double lep_pt1=0., lep_pt2=0., lep_pt3=0.;
      for(int m=0; m<4; m++) {
        if(temp.lepton[m].Pt()>lep_pt1) lep_pt1=temp.lepton[m].Pt();
      }
      for(int m=0; m<4; m++) {
        if(temp.lepton[m].Pt()>lep_pt2 && temp.lepton[m].Pt()<lep_pt1)
          lep_pt2=temp.lepton[m].Pt();
      }
      for(int m=0; m<4; m++) {
        if(temp.lepton[m].Pt()>lep_pt3 && temp.lepton[m].Pt()<lep_pt2)
          lep_pt3=temp.lepton[m].Pt();
      }
      bool pass_pt = lep_pt1>20.e3 && lep_pt2>20.e3 && lep_pt3>10.e3;
      if(!pass_pt) continue;

      temp.flavor_lep.push_back(pairs[i].flavor);
      temp.flavor_lep.push_back(pairs[i].flavor);
      temp.flavor_lep.push_back(pairs[j].flavor);
      temp.flavor_lep.push_back(pairs[j].flavor);

      temp.ZZ = pairs[i].Z + pairs[j].Z;
      if(pairs[i].flavor==0&&pairs[j].flavor==0) temp.type="eeee";
      if(pairs[i].flavor==1&&pairs[j].flavor==0) temp.type="eemm";
      if(pairs[i].flavor==0&&pairs[j].flavor==1) temp.type="eemm";
      if(pairs[i].flavor==1&&pairs[j].flavor==1) temp.type="mmmm";

      if(temp.pair[0].flavor==temp.pair[1].flavor) {
        temp.alter_pairs.push_back(temp.lepton[0]+temp.lepton[3]);
        temp.alter_pairs.push_back(temp.lepton[2]+temp.lepton[1]);
      }

      double mindR1=999., mindR2=999.;
      for(int n=0; n<4; n++) {
        for(int m=n+1; m<4; m++) {
          if(temp.flavor_lep[n]==temp.flavor_lep[m])
            if(temp.lepton[n].DeltaR(temp.lepton[m])<mindR1)
              mindR1=temp.lepton[n].DeltaR(temp.lepton[m]);

          if(temp.flavor_lep[n]!=temp.flavor_lep[m])
            if(temp.lepton[n].DeltaR(temp.lepton[m])<mindR2)
              mindR2=temp.lepton[n].DeltaR(temp.lepton[m]);
        }
      }
      if(mindR1<0.1 || mindR2<0.2) continue;

      bool passJPsiVeto = true;
      if(temp.pair[0].flavor==temp.pair[1].flavor) {
        for(int n=0; n<(int)temp.alter_pairs.size(); n++)
          if(temp.alter_pairs[n].M()<5.e3) passJPsiVeto=false;
      }
      if(temp.pair[0].Z.M()<5.e3) passJPsiVeto=false;
      if(temp.pair[1].Z.M()<5.e3) passJPsiVeto=false;
      if(!passJPsiVeto) continue;

      double massZ1 = temp.pair[0].Z.M();
      double massZ2 = temp.pair[1].Z.M();
      double pull = fabs(massZ1-ZMass) + fabs(massZ2-ZMass);
      if(pull<Maxipull) {
        Maxipull = pull;
        best_Quad = temp;
        best_i = i;
        if(temp.type=="eeee")  event_type = quadType::_4e;
        if(temp.type=="eemm")  event_type = quadType::_2e2mu;
        if(temp.type=="mmmm")  event_type = quadType::_4mu;
      }

      //           good_quads.push_back(temp);
    }
  }


  //     int best_evt_type=99;
  if(best_i!=-1) {

    double massZ1 = best_Quad.pair[0].Z.M();
    double massZ2 = best_Quad.pair[1].Z.M();
    bool passOnShell = massZ1>66.e3 && massZ1<116.e3;
    passOnShell = (massZ2>66.e3 && massZ2<116.e3) && passOnShell;
    if(passOnShell && best_Quad.type==type) {
      Final_Quad = best_Quad;
      return true;
    }else return false;
    //if(event_type<best_evt_type) {
    //  best_evt_type=event_type;
    //  Final_Quad = best_Quad;
    //}
  }else return false;


  //if(best_i!=-1 && Final_Quad.type==type) {
  //  SetFlag(FLAG_cut_temp,"Final",type, 1); 
  //  return true;
  //}
  //else return false;    

}





