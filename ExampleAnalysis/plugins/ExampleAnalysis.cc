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
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> outFile;
   vtxDir = new TFileDirectory(outFile->mkdir("Vertex"));
   muDir = new TFileDirectory(outFile->mkdir("Muon"));
   eleDir = new TFileDirectory(outFile->mkdir("Electron"));
   tauDir = new TFileDirectory(outFile->mkdir("Tau");
   photonDir = new TFileDirectory(outFile->mkdir("Phton"));
   jetDir = new TFileDirectory(outFile->mkdir("Jet"));
   metDir = new TFileDirectory(outFile->mkdir("MET"));

   hnVtx = vtxDir->make<TH1I>("nVtx","NVertices;nVertex;#entries",40,-0.5,39.5);
   hmuPt = muDir->make<TH1D>("pt","Muon p_T;p_T;#entries",500,0.,1000.); 
   helePt = eleDir->make<TH1D>("pt","Electron p_T;p_T;#entries",500,0.,1000.);
   htauPt = tau->make<TH1D>("pt","Tau p_T;p_T;#entries",500,0.,1000.);
   hphotonuPt = photonDir->make<TH1D>("pt","Photon p_T;p_T;#entries",500,0.,1000.);
   hjetPt = jetDir->make<TH1D>("pt","Jet p_T;p_T;#entries",500,0.,1000.);
   hmetEt = metDir->make<TH1D>("et","MET E_T;E_T;#entries",500,0.,1000.);

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
  //vertices
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();
  hnVtx->Fill(vertices->size());

  //muons
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  for (const pat::Muon &mu : *muons) {
    if (mu.pt() < 5 || !mu.isLooseMuon()) continue;
    hmuPt->Fill(mu.pt());
  }
  
  //electrons
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  for (const pat::Electron &el : *electrons) {
    if (el.pt() < 5) continue;
    helePt->Fill(ele.pt());
  }
  
  //photons
  edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken(photonToken_, photons);
  for(const pat::Photon &pho : *photons) {
    if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
    hphotonPt->Fill(pho.pt());
  }
  
  //tau
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauToken_, taus);
  for (const pat::Tau &tau : *taus) {
    if (tau.pt() < 20) continue;
    htauPt->Fill(tau.pt());
  }

  //Jets
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  for (const pat::Jet &j : *jets) {
    if (j.pt() < 20) continue;
    hjetPt->Fill(j.pt());
  }
  
  //MET
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();
  hmetPt->Fill(met.pt());
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
