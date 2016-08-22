#include "IndiaCMSTutorial/ZmumuTnP/plugins/ZmumuTnP.h"
ZmumuTnP::ZmumuTnP(const edm::ParameterSet& iConfig) :
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")))

{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


ZmumuTnP::~ZmumuTnP()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZmumuTnP::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  selectedMu_.clear();
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
    if(mu.pt() <= 5.) continue;
    if(std::fabs(mu.eta()) > 2.4) continue;
    reco::TrackRef tk = mu.muonBestTrack();
    double dxyWrtPV = tk->dxy(PV.position());
    double dzWrtPV = tk->dz(PV.position());
    if(std::fabs(dxyWrtPV) >= 0.5 )      continue;
    if(std::fabs(dzWrtPV) >= 1.) continue;
    bool quality = (mu.isGlobalMuon()
                    || ( mu.isTrackerMuon() && mu.numberOfMatches() >= 0)) 
                    && mu.muonBestTrackType()!=2 ;
    if(!quality) continue;
    if(std::fabs(mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D)) >= 4.)    continue;
    if(!mu.isPFMuon())    continue;
    selectedMu_.push_back(mu);
  }
  if(selectedMu_.size() > 2)   selectZmumu();
  

}

void ZmumuTnP::selectZmumu() {
  for(unsigned int i = 0; i < selectedMu_.size(); i++) {
      
    if()  
    for(unsigned int j = 0; j < selectedMu_.size(); j++) {
    }
  }
}
// ------------ method called once each job just before starting event loop  ------------
void 
ZmumuTnP::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZmumuTnP::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZmumuTnP::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZmumuTnP);
