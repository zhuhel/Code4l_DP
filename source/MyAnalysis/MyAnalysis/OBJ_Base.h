#ifndef MyAnalysis_OBJ_H
#define MyAnalysis_OBJ_H
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <xAODMuon/Muon.h>
#include <xAODEgamma/Electron.h>
#include <xAODJet/Jet.h>

//<Common Header>
#include <TLorentzVector.h>
#include <string>
#include <vector>
using namespace std;

//<Predefined Constants>
const double PI=3.1415926;
const double ZMass=91.1876e3; //MeV
const double DWMass=160.77e3; //MeV
//const double m_mass=105.658367; //MeV
//const double e_mass = 0.510998910; //MeV
const double m_mass=105.658; //MeV
const double e_mass = 0.511; //MeV
const double Unit_GeV = 1000.; //MeV
int countinggood=0;
int countingbad=0;
int goodeee=0;
//<Object Class>
//<General, Muon, Electron, Jet, MET>
class OBJ{
  public:
    double E,Et;
    double p,pt,px,py,pz;
    double eta,phi,m;
    double d0,z0,d0err,z0err,d0sig,z0sintheta;
    double sf;
    int index, vtxtype;
    int passIso;
    bool trigM; 
    TLorentzVector L,Lori,Lcorr,LF;       
};    

class OBJ_MUON : public OBJ{
  public:
    int author, charge, type, quality; // type is added by Cong
    int iscombined, isloose, ismedium, istight, islowpt, istag, issa; // classified into type and quality
    double id; /* 1: staco, 2: ST, 3: SA, 4: calo */
    double caloId, caloLikelihood;
    int passhits_staco, passhits_calo, passhits_sa, passhits_3rdChain;
    int passhits_blayer, passhits_pix, passhits_sct, passhits_hole, passhits_trt;
    int isgood, iscr; //too define control region. isgood==1 for std good muon
    double ptme, ptid, ptms, pt_corr, ptme_corr, ptid_corr;
    double etcone20, etcone30, etcone40, ptcone20, ptcone30, ptcone40, etcone20_corr, ptcone20_corr;
    double topoetcone20;
    double mu_calo_eta, mu_calo_phi, mu_staco_eta, mu_staco_phi, id_theta;
    double met_corrx, met_corry;
    double sf_loose, sf_calo;
    double npv;
    vector<int> fsr_index;
    vector<double> fsr_dR;
    vector<string> fsr_type;
    TLorentzVector L_me, L_ms, L_id, Lcorr_me, Lcorr_ms, Lcorr_id, LF_me, LF_id;
    xAOD::Muon* ptrMuon;
};

class OBJ_ELECTRON: public OBJ{
  public:
    int source; //1: from Z, 2: from b/c, 3: from hadron, 4: from photon conversion, 5: unknown, 6: from W
    int iscr, isgood; //too define control region. isgood==1 for std good electron, iscr==1 for cr
    double charge;
    uint16_t author;
    int isloose, ismedium, istight, isloosepp, ismediumpp, istightpp, islooseppzz, ismultilep, islikelihood;
    bool passOQ;
    double clE, clpt, cleta, clphi, trkpt, trketa, trktheta, trkphi, qoverp;
    double trkpt2, trketa2, trkphi2;
    double clE_corr, pt2_corr, clE_noEP;
    double etcone20, etcone20_corr, etcone30, etcone40, ptcone20, ptcone30, ptcone40, ptcone20_corr;
    double topoetcone20;
    double el_etas2, el_phis2, el_rawcl_E;
    double met_corrx, met_corry;
    double sf_loosepp, sf_mediumpp, sf_tightpp, sf_id, sf_reco;
    int npv;
    TLorentzVector Lcl,L_trk,LF_trk,LF_tri; 
    xAOD::Electron* ptrElectron;
};       

class OBJ_JET : public OBJ{
  public:
    double pt_em, eta_em, phi_em, m_em;
    double isbad, jes, isugly, islooserbad;
    int isbjet77, isbjet85;
    double vtxf, fmax, smax;
    TLorentzVector LEM;
    TVector3 L3;
    xAOD::Jet* ptrJet;
};

class OBJ_MET : public OBJ{
  public:
    double pt_corr,px_corr,py_corr;
    TLorentzVector L_RefFinal,L_Topo;
    TLorentzVector L_RefFinal_Corr,L_Topo_Corr;
};

//class for lepton pairs
class Pair {
  public:
    int flavor; //0-electron, 1-muon
    vector<int> index;
    vector<int> charge;
    vector<TLorentzVector> lepton;
    vector<TLorentzVector> lepton_trk;
    TLorentzVector Z;
    double mass;
    vector<xAOD::Muon*> ptrVMuon;
    vector<xAOD::Electron*> ptrVElectron;
    vector<int> passIso;
};

//class for quads, i.e. Z pairs
class Quad {
  public:
    vector<int> index, medium;
    vector<Pair> pair;
    vector<TLorentzVector> alter_pairs;
    vector<TLorentzVector> lepton;
    vector<TLorentzVector> lepton_trk;
    vector<int> index_lep;
    vector<int> flavor_lep;
    vector<int> charge_lep;
    vector<float> ptcone20, etcone20, ptcone20_corr, etcone20_corr, d0sig, z0sintheta, d0, z0;
    vector<double> topoetcone20;
    TLorentzVector ZZ;
    double weight, trigSF;
    // for storage
    string type;
    int event_type;
    int prod_type;
    int mass_type;
    double mz1_fsr;
    double mz2_fsr;
    double mzz_fsr;
    double mz1_constrained;
    double mz2_constrained;
    double m4l_constrained;
    double mass2jet;
    double deta2jet;
    double jmax_pt;
    double jmax_eta;
    double jsec_pt;
};


typedef vector<OBJ_MUON> VOmuon;
typedef vector<OBJ_ELECTRON> VOelectron;
typedef vector<OBJ_JET> VOjet;
typedef vector<OBJ_MET> VOmet;


#endif
