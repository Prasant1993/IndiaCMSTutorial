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
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
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
      void selectZmumu();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      //TFileDirectory* zee;

      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
       
      //std::vector<pat::Electron>  selectedEle_;     
};
#endif
