// -*- C++ -*-
//
// Package:    IndiaCMSTutorial/ExampleAnalysis
// Class:      ExampleAnalysis
// 
/**\class ExampleAnalysis ExampleAnalysis.cc IndiaCMSTutorial/ExampleAnalysis/plugins/ExampleAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Suvankar Roy Chowdhury
//         Created:  Thu, 18 Aug 2016 06:25:43 GMT
//
//

#include "IndiaCMSTutorial/ExampleAnalysis/plugins/ExampleAnalysis.h"

//
// constructors and destructor
//
ExampleAnalysis::ExampleAnalysis(const edm::ParameterSet& iConfig):
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> outFile;
   muDir = new TFileDirectory(outFile->mkdir("Muon"));
   muPt = muDir->make<TH1D>("muonPt","MuonPT;p_T;#entries",500,0.,1000.);
}


ExampleAnalysis::~ExampleAnalysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ExampleAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  for (const pat::Muon &mu : *muons) {
    if (mu.pt() < 5 || !mu.isLooseMuon()) continue;
    muPt->Fill(mu.pt());
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
ExampleAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExampleAnalysis::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ExampleAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExampleAnalysis);
