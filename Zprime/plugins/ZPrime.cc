// -*- C++ -*-
//
// Package:    IndiaCMSTutorial/ZPrime
// Class:      ZPrime
// 
/**\class ZPrime ZPrime.cc IndiaCMSTutorial/ZPrime/plugins/ZPrime.cc
*/
// constructors and destructor
//
#include "IndiaCMSTutorial/Zprime/plugins/ZPrime.h"
ZPrime::ZPrime(const edm::ParameterSet& iConfig) :
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  prunedGenToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("packed")))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  //edm::Service<TFileService> outFile;
  //**Book Histograms here
}


ZPrime::~ZPrime()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZPrime::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  selectedMu_.clear();
  genZpdaughter_.clear();
  //vertices
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();
  
  // Pruned particles are the one containing "important" stuff
  edm::Handle<std::vector<reco::GenParticle>> pruned;
  iEvent.getByToken(prunedGenToken_,pruned);
  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  edm::Handle<std::vector<pat::PackedGenParticle>> packed;
  iEvent.getByToken(packedGenToken_,packed);
  //do GenLevel Analysis
  checkGenlevel(pruned,packed);

  //muons
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  
  for (const pat::Muon &mu : *muons) {
    //Apply muon cuts here
    selectedMu_.push_back(mu);
  }
  
  //call selectrecoDimuons() to select the best dimuon pair
  if(selectedMu_.size() >= 2)   selectrecoDimuons();
  //Fill reco level histograms if dimuon found

  //Gen level Analysis
  checkGenlevel(pruned, packed);
  //Fill GenLevel histos
}
//dimuon formation
//out of all selected muons, select the OS pair with highest sun pt
//save the dimuon properties in the Dimuon struct             
void ZPrime::selectrecoDimuons() {
  bestZpcand_.mu1P4.SetPtEtaPhiE(0.,0.,0.,0.);
  bestZpcand_.mu2P4.SetPtEtaPhiE(0.,0.,0.,0.);
  bestZpcand_.mZp = 0.;
  bestZpcand_.avgPt = 0.;
  for(unsigned int i = 0; i<selectedMu_.size(); i++) {
    TLorentzVector mu1P4 = getP4(selectedMu_[i]);
    for(unsigned int j = i+1; j<selectedMu_.size(); j++) { 
      TLorentzVector mu2P4 = getP4(selectedMu_[j]);
    }
  }
}

//Check recursively if any ancestor of particle is the given one
bool ZPrime::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;
  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++) {
    if(isAncestor(ancestor,particle->mother(i))) return true;
  }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}

//check Zp to mumu at gen level
void ZPrime::checkGenlevel(const edm::Handle<std::vector<reco::GenParticle>>& pruned,              
                           const edm::Handle<std::vector<pat::PackedGenParticle>>& packed) {
  
  //let's try to find all status1 originating directly ZPrime decay 
  for(const auto& gp : *pruned) { 
    //std::cout << "PdgID: " << gp.pdgId() 
    //          << " pt " << gp.pt() 
    //          << " eta: " << gp.eta() 
    //          << " phi: " << gp.phi() 
    //          << " phi: " << gp.mass() 
    //          << std::endl;
    
    if(std::abs(gp.pdgId()) == 23 && gp.status() == 22) {
      //std::cout << "  found daugthers for Z with status: " << gp.status() << std::endl;
      for(size_t j=0; j<packed->size();j++){
        //get the pointer to the first survied ancestor of a given 
	//packed GenParticle in the prunedCollection 
        if(std::fabs(packed->at(j).pdgId()) != 13)     continue;
	const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
	if(motherInPrunedCollection != nullptr && isAncestor( &gp , motherInPrunedCollection)){
          genZpdaughter_.push_back(packed->at(j));
	  //std::cout << "     PdgID: " << (*packed)[j].pdgId() 
	  //	      << " pt " << (*packed)[j].pt() 
	  //	      << " eta: " << (*packed)[j].eta() 
	  //	      << " phi: " << (*packed)[j].phi() << std::endl;
	}
      }
    }
    if(genZpdaughter_.size()==2)  break;
  }
}
// ------------ method called once each job just before starting event loop  ------------
void 
ZPrime::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZPrime::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZPrime::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZPrime);
