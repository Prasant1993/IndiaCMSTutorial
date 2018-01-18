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
  edm::Service<TFileService> outFile;
  mZp_reco = outFile->make<TH1D>("mZp","ZRrime Mass;m_Zp;#entries",2000,0.,4000.); 
  mu1Pt_reco = outFile->make<TH1D>("muon1Pt","Leading lepton Pt;mu1_pT;#entries",2000,0.,4000.);
  mu2Pt_reco = outFile->make<TH1D>("muon2Pt","Sub-Leading lepton Pt;mu2_pT;#entries",2000,0.,4000.);
  mZp_genDau = outFile->make<TH1D>("mZp_gen","ZRrime Mass from Daughters at Gen Level;m_Zp(Gen);#entries",2000,0.,4000.); 
  mu1Pt_gen = outFile->make<TH1D>("muon1Pt_gen","Leading lepton Pt at Gen Level;mu1_pT(Gen);#entries",2000,0.,4000.);
  mu2Pt_gen = outFile->make<TH1D>("muon2Pt_gen","Sub-Leading lepton Pt at Gen Level;mu2_pT(Gen);#entries",2000,0.,4000.);
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
    if (!mu.isGlobalMuon())   continue;
    if (!mu.isTrackerMuon())  continue;
    const reco::TrackRef& muTrack = mu.muonBestTrack();
    if (muTrack->ptError()/muTrack->pt() > 0.3)  continue;
    if (muTrack->pt() <= 53.)  continue;
    if (std::fabs(muTrack->dxy(PV.position())) > 0.2)  continue;
    if (std::fabs(mu.dB()) > 0.2)    continue;
    if (mu.isolationR03().sumPt / mu.innerTrack()->pt() > 0.10 )  continue;
    if (mu.globalTrack()->hitPattern().trackerLayersWithMeasurement() <= 5)  continue;
    if (mu.globalTrack()->hitPattern().numberOfValidPixelHits() == 0)  continue;
    if (mu.globalTrack()->hitPattern().numberOfValidMuonHits() == 0)   continue;
    if (mu.numberOfMatchedStations() < 2)  continue;
    selectedMu_.push_back(mu);
  }
  //std::cout << "Muon Size=" << muons->size() << "::Selected>>" << selectedMu_.size() << std::endl;
  if(selectedMu_.size() >= 2)   selectrecoDimuons();
  //std::cout << "Zprime Candidate=" << bestZpcand_.mZp << std::endl;
  if(bestZpcand_.mZp > 0.) {
    mZp_reco->Fill(bestZpcand_.mZp);
    if(bestZpcand_.mu1P4.Pt() > bestZpcand_.mu2P4.Pt()) {
      mu1Pt_reco->Fill(bestZpcand_.mu1P4.Pt());
      mu2Pt_reco->Fill(bestZpcand_.mu2P4.Pt());
    } else {
      mu2Pt_reco->Fill(bestZpcand_.mu1P4.Pt());
      mu1Pt_reco->Fill(bestZpcand_.mu2P4.Pt());
    }
  }
  //Gen level Analysis
  checkGenlevel(pruned, packed);
  if(genZpdaughter_.size() > 1) {
    TLorentzVector gmu1 = getP4(genZpdaughter_[0]);
    TLorentzVector gmu2 = getP4(genZpdaughter_[1]);
    mZp_genDau->Fill((gmu1+gmu2).M());
    if(gmu1.Pt() > gmu2.Pt()) { 
      mu1Pt_gen->Fill(gmu1.Pt());
      mu2Pt_gen->Fill(gmu2.Pt());
    } else {
      mu1Pt_gen->Fill(gmu2.Pt());
      mu2Pt_gen->Fill(gmu1.Pt());
    }   
  }
}
void ZPrime::selectrecoDimuons() {
  bestZpcand_.mu1P4.SetPtEtaPhiE(0.,0.,0.,0.);
  bestZpcand_.mu2P4.SetPtEtaPhiE(0.,0.,0.,0.);
  bestZpcand_.mZp = 0.;
  bestZpcand_.avgPt = 0.;
  for(unsigned int i = 0; i<selectedMu_.size(); i++) {
    TLorentzVector mu1P4 = getP4(selectedMu_[i]);
    for(unsigned int j = i+1; j<selectedMu_.size(); j++) { 
      if(selectedMu_[i].charge() + selectedMu_[j].charge() != 0)   continue;
      TLorentzVector mu2P4 = getP4(selectedMu_[j]);
      if(bestZpcand_.avgPt < (mu1P4.Pt() + mu2P4.Pt())/2.) {
        bestZpcand_.mu1P4 = mu1P4;
        bestZpcand_.mu2P4 = mu2P4;
        bestZpcand_.mZp = (mu1P4 + mu2P4).M();
        bestZpcand_.avgPt = (mu1P4.Pt() + mu2P4.Pt())/2.;
      }
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
      std::cout << "  found daugthers for Z with status: " << gp.status() << std::endl;
      for(size_t j=0; j<packed->size();j++){
        //get the pointer to the first survied ancestor of a given 
	//packed GenParticle in the prunedCollection 
        if(std::fabs(packed->at(j).pdgId()) != 13)     continue;
	const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
	if(motherInPrunedCollection != nullptr && isAncestor( &gp , motherInPrunedCollection)){
          genZpdaughter_.push_back(packed->at(j));
	  std::cout << "     PdgID: " << (*packed)[j].pdgId() 
		    << " pt " << (*packed)[j].pt() 
		    << " eta: " << (*packed)[j].eta() 
		    << " phi: " << (*packed)[j].phi() << std::endl;
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
