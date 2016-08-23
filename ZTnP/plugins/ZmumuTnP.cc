#include "TLorentzVector.h"
#include "IndiaCMSTutorial/ZTnP/plugins/ZmumuTnP.h"

ZmumuTnP::ZmumuTnP(const edm::ParameterSet& iConfig) :
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> outTreeFile;
   outTree_ = outTreeFile->make<TTree>("zmmtree","muTnP tree for iso");
   outTree_->Branch("TnP_pt",TnP_pt,"TnP_pt[nTnP]/F");
   outTree_->Branch("TnP_eta",TnP_eta,"TnP_eta[nTnP]/F");
   outTree_->Branch("TnP_phi",TnP_phi,"TnP_phi[nTnP]/F");
   outTree_->Branch("TnP_mass",TnP_mass,"TnP_mass[nTnP]/F");
   //tag properties
   outTree_->Branch("TnP_l1_pdgId",TnP_l1_pdgId,"TnP_l1_pdgId[nTnP]/I");
   outTree_->Branch("TnP_l1_pt",TnP_l1_pt,"TnP_l1_pt[nTnP]/F");
   outTree_->Branch("TnP_l1_eta",TnP_l1_eta,"TnP_l1_eta[nTnP]/F");
   outTree_->Branch("TnP_l1_phi",TnP_l1_phi,"TnP_l1_phi[nTnP]/F");
   outTree_->Branch("TnP_l1_mass",TnP_l1_mass,"TnP_l1_mass[nTnP]/F");
   outTree_->Branch("TnP_l1_charge",TnP_l1_charge,"TnP_l1_charge[nTnP]/I");
   outTree_->Branch("TnP_l1_relIso",TnP_l1_relIso,"TnP_l1_relIso[nTnP]/F");
   //probe properties
   outTree_->Branch("TnP_l2_pdgId",TnP_l2_pdgId,"TnP_l2_pdgId[nTnP]/I");
   outTree_->Branch("TnP_l2_pt",TnP_l2_pt,"TnP_l2_pt[nTnP]/F");
   outTree_->Branch("TnP_l2_eta",TnP_l2_eta,"TnP_l2_eta[nTnP]/F");
   outTree_->Branch("TnP_l2_phi",TnP_l2_phi,"TnP_l2_phi[nTnP]/F");
   outTree_->Branch("TnP_l2_mass",TnP_l2_mass,"TnP_eta[nTnP]/F");
   outTree_->Branch("TnP_l2_charge",TnP_l2_charge,"TnP_l2_charge[nTnP]/I");
   outTree_->Branch("TnP_l2_relIso",TnP_l2_relIso,"TnP_l2_relIso[nTnP]/F");
   
   histoDir = new TFileDirectory(outFile->mkdir("Histos_zmumu"));
   mZmmAll_ptl50_barrel = histoDir->make<TH1D>("mZmmAll_ptl50_barrel","Z Mass;m_Z;#entries",25,70.,120.); 
   mZmmPass_ptl50_barrel  = histoDir->make<TH1D>("mZmmPass_ptl50_barrel","Z Mass;m_Z;#entries",25,70.,120.);
   mZmmFail_ptl50_barrel = histoDir->make<TH1D>("mZmmFail_ptl50_barrel","Z Mass;m_Z;#entries",25,70.,120.);
   mZmmAll_ptl50_endcap = histoDir->make<TH1D>("mZmmAll_ptl50_endcap","Z Mass;m_Z;#entries",25,70.,120.);
   mZmmPass_ptl50_endcap = histoDir->make<TH1D>("mZmmPass_ptl50_endcap","Z Mass;m_Z;#entries",25,70.,120.);
   mZmmFail_ptl50_endcap = histoDir->make<TH1D>("mZmmFail_ptl50_endcap","Z Mass;m_Z;#entries",25,70.,120.);
   mZmmAll_ptg50_barrel = histoDir->make<TH1D>("mZmmAll_ptg50_barrel","Z Mass;m_Z;#entries",25,70.,120.);
   mZmmPass_ptg50_barrel = histoDir->make<TH1D>("mZmmPass_ptg50_barrel","Z Mass;m_Z;#entries",25,70.,120.);
   mZmmFail_ptg50_barrel = histoDir->make<TH1D>("mZmmFail_ptg50_barrel","Z Mass;m_Z;#entries",25,70.,120.);
   mZmmAll_ptg50_endcap = histoDir->make<TH1D>("mZmmAll_ptg50_endcap","Z Mass;m_Z;#entries",25,70.,120.);
   mZmmPass_ptg50_endcap = histoDir->make<TH1D>("mZmmPass_ptg50_endcap","Z Mass;m_Z;#entries",25,70.,120.);
   mZmmFail_ptg50_endcap = histoDir->make<TH1D>("mZmmFail_ptg50_endcap","Z Mass;m_Z;#entries",25,70.,120.);
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
//tag selection
//pt > 20. and rel_iso03 < 0.35
//tag + probe pair = OS + 80 < mass < 100
void ZmumuTnP::selectZmumu() {
  int nTnP = 0;
  for(unsigned int i = 0; i < selectedMu_.size(); i++) {
    double tagrelIso = pfiso(selectedMu_[i])/selectedMu_[i].pt();  
    if(tagrelIso >= 0.35)     continue;
    if(selectedMu_[i].pt() <= 20.)    continue;
    TLorentzVector tagP4 = getP4(selectedMu_[i]);
    int tagcharge =  selectedMu_[i].charge();  
    for(unsigned int j = 0; j < selectedMu_.size(); j++) {
      int probecharge = selectedMu_[j].charge();
      if(tagcharge + probecharge != 0)      continue;
      if(nTnP >= 8)    continue; 
      TLorentzVector probeP4 = getP4(selectedMu_[j]);
      TLorentzVector TnP = (tagP4 + probeP4);
      if( TnP.M() <= 80. || TnP.M() >= 100. )   continue;
      double proberelIso = pfiso(selectedMu_[j])/selectedMu_[j].pt();
      TnP_pt[nTnP] = TnP.Pt();   
      TnP_eta[nTnP] = TnP.Eta();   
      TnP_phi[nTnP] = TnP.Phi();   
      TnP_mass[nTnP] = TnP.M();   
      //tag properties 
      TnP_l1_pdgId[nTnP] = (tagcharge > 0) ? -13 : 13;   
      TnP_l1_pt[nTnP] = tagP4.Pt();   
      TnP_l1_eta[nTnP] = tagP4.Eta();   
      TnP_l1_phi[nTnP] = tagP4.Phi();   
      TnP_l1_mass[nTnP] = tagP4.M();   
      TnP_l1_charge[nTnP] = tagcharge;   
      TnP_l1_relIso[nTnP] = tagrelIso;
      //probe properties
      TnP_l2_pdgId[nTnP] = (probecharge > 0) ? -13 : 13;   
      TnP_l2_pt[nTnP] = probeP4.Pt();   
      TnP_l2_eta[nTnP] = probeP4.Eta();   
      TnP_l2_phi[nTnP]= probeP4.Phi();   
      TnP_l2_mass[nTnP]= probeP4.M();   
      TnP_l2_charge[nTnP]= probecharge;   
      TnP_l2_relIso[nTnP]= proberelIso;
      nTnP++;
      if(probeP4.Pt() <= 50.) {
        if(std::fabs(probeP4.Eta()) < 1.2 ) {
          mZmmAll_ptl50_barrel->Fill(TnP.M());
          mZmmPass_ptl50_barrel->Fill(TnP.M());
          mZmmFail_ptl50_barrel->Fill(TnP.M());
        } else {
          mZmmAll_ptl50_endcap->Fill(TnP.M());
          mZmmPass_ptl50_endcap->Fill(TnP.M());
          mZmmFail_ptl50_endcap->Fill(TnP.M());
        }
      } else {
        if(std::fabs(probeP4.Eta()) < 1.2 ) {
          mZmmAll_ptg50_barrel->Fill(TnP.M());
          mZmmPass_ptg50_barrel->Fill(TnP.M());
          mZmmFail_ptg50_barrel->Fill(TnP.M());
        } else {
          mZmmAll_ptg50_endcap->Fill(TnP.M());
          mZmmPass_ptg50_endcap->Fill(TnP.M());
          mZmmFail_ptg50_endcap->Fill(TnP.M());
        }
      }
    }
  }
}

//function to compute muon isolation

double ZmumuTnP::pfiso(const pat::Muon& mu) {
  return (mu.pfIsolationR03().sumChargedHadronPt 
          + std::max(0., mu.pfIsolationR03().sumNeutralHadronEt
                       + mu.pfIsolationR03().sumPhotonEt - 0.5 * mu.pfIsolationR03().sumPUPt));
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
