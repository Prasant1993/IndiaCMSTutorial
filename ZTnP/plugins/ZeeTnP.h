#ifndef __ZeeTnP_h
#define __ZeeTnP_h
// -*- C++ -*-
//
// Package:    IndiaCMSTutorial/ZeeTnP
// Class:      ZeeTnP
// 
/**\class ZeeTnP ZeeTnP.cc IndiaCMSTutorial/ZeeTnP/plugins/ZeeTnP.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
// system include files
#include <memory>
#include <TDirectory.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TVector2.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

//
// class declaration
//
class ZeeTnP : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
 public:
  explicit ZeeTnP(const edm::ParameterSet&);
  ~ZeeTnP();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  template<class T>
    TLorentzVector getP4(const T& obj) {
    TLorentzVector P;
    P.SetPtEtaPhiE(obj.pt(),obj.eta(),obj.phi(),obj.energy());
    return P;
  }
  double combinedRelativeIso(const pat::Electron&);
  void selectZee();
  
 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  //TFileDirectory* zee;
  
  //edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  
  std::vector<pat::Electron>  selectedEle_;     
  const static Int_t kMaxTnP = 8;
  TTree* outTree_;
  int nselectedEle;
  int nTnP;
  float         TnP_pt[kMaxTnP];   
  float         TnP_eta[kMaxTnP];   
  float         TnP_phi[kMaxTnP];   
  float         TnP_mass[kMaxTnP];   
  //tag properties 
  int           TnP_l1_pdgId[kMaxTnP];   
  float         TnP_l1_pt[kMaxTnP];   
  float         TnP_l1_eta[kMaxTnP];   
  float         TnP_l1_phi[kMaxTnP];   
  float         TnP_l1_mass[kMaxTnP];   
  int           TnP_l1_charge[kMaxTnP];   
  float         TnP_l1_relIso[kMaxTnP];
  //probe properties
  int           TnP_l2_pdgId[kMaxTnP];   
  float         TnP_l2_pt[kMaxTnP];   
  float         TnP_l2_eta[kMaxTnP];   
  float         TnP_l2_phi[kMaxTnP];   
  float         TnP_l2_mass[kMaxTnP];   
  int           TnP_l2_charge[kMaxTnP];   
  float         TnP_l2_relIso[kMaxTnP];
  //Histogram booking
  TFileDirectory* histoDir;
  TH1D* mZeeAll_ptl50_barrel;
  TH1D* mZeePass_ptl50_barrel;
  TH1D* mZeeFail_ptl50_barrel;
  TH1D* mZeeAll_ptl50_endcap;
  TH1D* mZeePass_ptl50_endcap;
  TH1D* mZeeFail_ptl50_endcap;
  
  TH1D* mZeeAll_ptg50_barrel;
  TH1D* mZeePass_ptg50_barrel;
  TH1D* mZeeFail_ptg50_barrel;
  TH1D* mZeeAll_ptg50_endcap;
  TH1D* mZeePass_ptg50_endcap;
  TH1D* mZeeFail_ptg50_endcap;
};
#endif
