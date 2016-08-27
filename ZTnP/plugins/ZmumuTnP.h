#ifndef __ZmumuTnP_h
#define __ZmumuTnP_h
// -*- C++ -*-
//
// Package:    IndiaCMSTutorial/ZmumuTnP
// Class:      ZmumuTnP
// 
/**\class ZmumuTnP ZmumuTnP.cc IndiaCMSTutorial/ZmumuTnP/plugins/ZmumuTnP.cc

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

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//
// class declaration
//
class ZmumuTnP : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ZmumuTnP(const edm::ParameterSet&);
      ~ZmumuTnP();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      template<class T>
      TLorentzVector getP4(const T& obj) {
        TLorentzVector P;
        P.SetPtEtaPhiE(obj.pt(),obj.eta(),obj.phi(),obj.energy());
        return P;
      }
      void selectZmumu();
      
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      double mupfiso(const pat::Muon& mu);
      bool isMatchedtoTrigger (const pat::Muon& mu, double hlt2reco_deltaRmax);
      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
	  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
	  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
       
      std::vector<pat::Muon>  selectedMu_;
      std::vector<pat::TriggerObjectStandAlone> triggerObj_;
      const static Int_t kMaxTnP = 8;
      TTree* outTree_;
      int nselectedMu;   
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
      
      std::string hltPathName;
      std::string hltFilterName;
	  double DeltaR_, dEta, dPhi, dR;
	  bool isPath, isFilter, isMatched;
	  edm::TriggerNames names;	
      //Histogram booking
      TFileDirectory* histoDir;
      TH1D* mZmmAll_ptl50_barrel;
      TH1D* mZmmPass_ptl50_barrel;
      TH1D* mZmmFail_ptl50_barrel;
      TH1D* mZmmAll_ptl50_endcap;
      TH1D* mZmmPass_ptl50_endcap;
      TH1D* mZmmFail_ptl50_endcap;

      TH1D* mZmmAll_ptg50_barrel;
      TH1D* mZmmPass_ptg50_barrel;
      TH1D* mZmmFail_ptg50_barrel;
      TH1D* mZmmAll_ptg50_endcap;
      TH1D* mZmmPass_ptg50_endcap;
      TH1D* mZmmFail_ptg50_endcap;
};
#endif
