#ifndef __ExampleAnalysis_h
#define __ExampleAnalysis_h
// system include files
#include <memory>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class ExampleAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ExampleAnalysis(const edm::ParameterSet&);
      ~ExampleAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      TDirectory* muDir;
      TH1D* muPt;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
};

#endif
