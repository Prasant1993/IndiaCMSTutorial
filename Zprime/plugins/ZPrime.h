#ifndef __ZPrime_h
#define __ZPrime_h
// system include files
#include <memory>

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
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

struct Dimuon {
  TLorentzVector mu1P4;
  TLorentzVector mu2P4;
  int mu1Charge;
  int mu2Charge;
  double mZp;
  double avgPt;
}

class ZPrime : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ZPrime(const edm::ParameterSet&);
      ~ZPrime();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      template<class T>
      TLorentzVector getP4(const T& obj) {
        TLorentzVector P;
        P.SetPtEtaPhiE(obj.pt(),obj.eta(),obj.phi(),obj.energy());
        return P;
      }
      void selectrecoDimuons();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
       
      std::vector<pat::Muon>  selectedMu_;
      Dimuon bestZpcand_;
      TH1D* mZp_reco;
      TH1D* mu1Pt_reco;
      TH1D* mu2Pt_reco;
};
#endif
