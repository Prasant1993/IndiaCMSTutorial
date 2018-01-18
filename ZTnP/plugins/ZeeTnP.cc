#include "IndiaCMSTutorial/ZTnP/plugins/ZeeTnP.h"
ZeeTnP::ZeeTnP(const edm::ParameterSet& iConfig) :
  //vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")))
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> outFile;
  outTree_ = outFile->make<TTree>("zeetree","eleTnP tree");
  outTree_->Branch("nselectedEle",&nselectedEle,"nselectedEle/I");
  outTree_->Branch("nTnP",&nTnP,"nTnP/I");
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
  outTree_->Branch("TnP_l1_sigmaIetaIeta",TnP_l1_sigmaIetaIeta,"TnP_l1_sigmaIetaIeta[nTnP]/F");
  //probe properties
  outTree_->Branch("TnP_l2_pdgId",TnP_l2_pdgId,"TnP_l2_pdgId[nTnP]/I");
  outTree_->Branch("TnP_l2_pt",TnP_l2_pt,"TnP_l2_pt[nTnP]/F");
  outTree_->Branch("TnP_l2_eta",TnP_l2_eta,"TnP_l2_eta[nTnP]/F");
  outTree_->Branch("TnP_l2_phi",TnP_l2_phi,"TnP_l2_phi[nTnP]/F");
  outTree_->Branch("TnP_l2_mass",TnP_l2_mass,"TnP_l2_mass[nTnP]/F");
  outTree_->Branch("TnP_l2_charge",TnP_l2_charge,"TnP_l2_charge[nTnP]/I");
  outTree_->Branch("TnP_l2_relIso",TnP_l2_relIso,"TnP_l2_relIso[nTnP]/F");
  outTree_->Branch("TnP_l2_sigmaIetaIeta",TnP_l2_sigmaIetaIeta,"TnP_l2_sigmaIetaIeta[nTnP]/F");
  
  histoDir = new TFileDirectory(outFile->mkdir("Histos_zee"));
  mZeeAll_ptl50_barrel = histoDir->make<TH1D>("mZeeAll_ptl50_barrel","Z Mass;m_Z;#entries",80,70.,110.); 
  mZeePass_ptl50_barrel  = histoDir->make<TH1D>("mZeePass_ptl50_barrel","Z Mass;m_Z;#entries",80,70.,110.);
  mZeeFail_ptl50_barrel = histoDir->make<TH1D>("mZeeFail_ptl50_barrel","Z Mass;m_Z;#entries",80,70.,110.);
  mZeeAll_ptl50_endcap = histoDir->make<TH1D>("mZeeAll_ptl50_endcap","Z Mass;m_Z;#entries",80,70.,110.);
  mZeePass_ptl50_endcap = histoDir->make<TH1D>("mZeePass_ptl50_endcap","Z Mass;m_Z;#entries",80,70.,110.);
  mZeeFail_ptl50_endcap = histoDir->make<TH1D>("mZeeFail_ptl50_endcap","Z Mass;m_Z;#entries",80,70.,110.);
  mZeeAll_ptg50_barrel = histoDir->make<TH1D>("mZeeAll_ptg50_barrel","Z Mass;m_Z;#entries",80,70.,110.);
  mZeePass_ptg50_barrel = histoDir->make<TH1D>("mZeePass_ptg50_barrel","Z Mass;m_Z;#entries",80,70.,110.);
  mZeeFail_ptg50_barrel = histoDir->make<TH1D>("mZeeFail_ptg50_barrel","Z Mass;m_Z;#entries",80,70.,110.);
  mZeeAll_ptg50_endcap = histoDir->make<TH1D>("mZeeAll_ptg50_endcap","Z Mass;m_Z;#entries",80,70.,110.);
  mZeePass_ptg50_endcap = histoDir->make<TH1D>("mZeePass_ptg50_endcap","Z Mass;m_Z;#entries",80,70.,110.);
  mZeeFail_ptg50_endcap = histoDir->make<TH1D>("mZeeFail_ptg50_endcap","Z Mass;m_Z;#entries",80,70.,110.);
}


ZeeTnP::~ZeeTnP()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZeeTnP::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  selectedEle_.clear();
  nTnP =0;
  //vertices
  //edm::Handle<reco::VertexCollection> vertices;
  //iEvent.getByToken(vtxToken_, vertices);
  //if (vertices->empty()) return; // skip the event if no PV found
  //const reco::Vertex &PV = vertices->front();
  //hnVtx->Fill(vertices->size());
  
  //electrons
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  for (const pat::Electron &ele : *electrons) {
    if(!ele.isPF())    continue;
    if(ele.pt() <= 20.) continue;
    reco::SuperClusterRef ele_sc = ele.superCluster();
    float eta = std::fabs(ele_sc->eta());
    if ((eta > 1.479 && eta < 1.556) || eta > 2.5) continue;

    reco::GsfTrackRef trk = ele.gsfTrack();
    double d0 = trk->d0();
    double dz = trk->dz();
    
    double ele_relIso = combinedRelativeIso(ele);
    if (eta <= 1.479) {
      //if (ele.full5x5_sigmaIetaIeta() >= 0.012) continue;
      if (ele.deltaEtaSuperClusterTrackAtVtx() >= 0.00311) continue;
      if (ele.deltaPhiSuperClusterTrackAtVtx() >= 0.103) continue;
      if (ele.hadronicOverEm() >= 0.253) continue;
      if (ele_relIso < 0.0695) continue;
      if (trk->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1) continue;
      if (d0 > 0.05) continue;
      if (dz > 0.1)  continue;
    } 
    else if (eta >= 1.556) {
      //if (ele.full5x5_sigmaIetaIeta() >= 0.04) continue;
      if (ele.deltaEtaSuperClusterTrackAtVtx() >= 0.00609) continue;
      if (ele.deltaPhiSuperClusterTrackAtVtx() >= 0.045) continue;
      if (ele.hadronicOverEm() >= 0.0878) continue;
      if (ele_relIso < 0.0821) continue;
      if (trk->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1) continue;
      if (d0 > 0.1)  continue;
      if (dz > 0.2)  continue;
    }
    selectedEle_.push_back(ele);
  }
  if(selectedEle_.size() >= 2)   selectZee();
  nselectedEle = selectedEle_.size();
  outTree_->Fill();
}
double ZeeTnP::combinedRelativeIso(const pat::Electron& electron) {
  float effArea = 0.0;
  reco::SuperClusterRef electron_sc = electron.superCluster();
  float ele_eta = std::fabs(electron_sc->eta());

  if (ele_eta <= 1.0)                     effArea = 0.1752;
  else if (ele_eta > 1.0 && ele_eta <= 1.479) effArea = 0.1862;
  else if (ele_eta > 1.479 && ele_eta <= 2.0) effArea = 0.1411;
  else if (ele_eta > 2.0 && ele_eta <= 2.2)   effArea = 0.1534;
  else if (ele_eta > 2.2 && ele_eta <= 2.3)   effArea = 0.1903;
  else if (ele_eta > 2.3 && ele_eta <= 2.4)   effArea = 0.2243;
  else if (ele_eta > 2.4)                     effArea = 0.2687;
  
  double scx = electron_sc->x();
  double scy = electron_sc->y();
  double scz = electron_sc->z();
  float rho = sqrt(scx*scx + scy*scy + scz*scz);
  
  reco::GsfElectron::PflowIsolationVariables pfIso = electron.pfIsolationVariables();
  double neuiso = pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - effArea * rho;
  double combinedRelIso = pfIso.sumChargedHadronPt + std::max(neuiso, 0.0) / (electron.pt());


  return combinedRelIso;
}

void ZeeTnP::selectZee() {
  
  for (uint i =0; i < selectedEle_.size(); i++) {
    double tagrelIso = combinedRelativeIso(selectedEle_[i]);
    if(selectedEle_[i].pt() <= 35.)    continue;
    double tagsigmaIetaIeta = selectedEle_[i].full5x5_sigmaIetaIeta();
    if(tagsigmaIetaIeta >= 0.00998) continue;
    TLorentzVector tagP4 = getP4(selectedEle_[i]);
    int tagcharge =  selectedEle_[i].charge();
    for(unsigned int j = 0; j < selectedEle_.size(); j++) {
      int probecharge = selectedEle_[j].charge();
      double probesigmaIetaIeta = selectedEle_[j].full5x5_sigmaIetaIeta();
      if(tagcharge + probecharge != 0)      continue;
      if(nTnP >= 8)    continue; 
      TLorentzVector probeP4 = getP4(selectedEle_[j]);
      TLorentzVector TnP = (tagP4 + probeP4);
      if( TnP.M() <= 80. || TnP.M() >= 110. )   continue;
      double proberelIso = combinedRelativeIso(selectedEle_[j]);
      TnP_pt[nTnP] = TnP.Pt();   
      TnP_eta[nTnP] = TnP.Eta();   
      TnP_phi[nTnP] = TnP.Phi();   
      TnP_mass[nTnP] = TnP.M();   
      //tag properties 
      TnP_l1_pdgId[nTnP] = (tagcharge > 0) ? -11 : 11;   
      TnP_l1_pt[nTnP] = tagP4.Pt();   
      TnP_l1_eta[nTnP] = tagP4.Eta();   
      TnP_l1_phi[nTnP] = tagP4.Phi();   
      TnP_l1_mass[nTnP] = tagP4.M();   
      TnP_l1_charge[nTnP] = tagcharge;   
      TnP_l1_relIso[nTnP] = tagrelIso;
      TnP_l1_sigmaIetaIeta[nTnP] = tagsigmaIetaIeta;
      //probe properties
      TnP_l2_pdgId[nTnP] = (probecharge > 0) ? -11 : 11;   
      TnP_l2_pt[nTnP] = probeP4.Pt();   
      TnP_l2_eta[nTnP] = probeP4.Eta();   
      TnP_l2_phi[nTnP]= probeP4.Phi();   
      TnP_l2_mass[nTnP]= probeP4.M();   
      TnP_l2_charge[nTnP]= probecharge;   
      TnP_l2_relIso[nTnP]= proberelIso;
      TnP_l2_sigmaIetaIeta[nTnP] = probesigmaIetaIeta;
      nTnP++;
      if(probeP4.Pt() <= 50.) {
        if(std::fabs(probeP4.Eta()) <= 1.479) {
          mZeeAll_ptl50_barrel->Fill(TnP.M());
	  if (selectedEle_[j].full5x5_sigmaIetaIeta() < 0.0115)      mZeePass_ptl50_barrel->Fill(TnP.M());
	  else                                                       mZeeFail_ptl50_barrel->Fill(TnP.M());
        } 
	else if (std::fabs(probeP4.Eta()) > 1.479 && std::fabs(probeP4.Eta()) < 2.5){
          mZeeAll_ptl50_endcap->Fill(TnP.M());
          if (selectedEle_[j].full5x5_sigmaIetaIeta() < 0.037)       mZeePass_ptl50_endcap->Fill(TnP.M());
          else                                                       mZeeFail_ptl50_endcap->Fill(TnP.M());
        }
      } 
      else {
        if(std::fabs(probeP4.Eta()) <= 1.479 ) {
          mZeeAll_ptg50_barrel->Fill(TnP.M());
          if (selectedEle_[j].full5x5_sigmaIetaIeta() < 0.0115)      mZeePass_ptg50_barrel->Fill(TnP.M());
          else                                                       mZeeFail_ptg50_barrel->Fill(TnP.M());
        } 
	else if (std::fabs(probeP4.Eta()) > 1.479 && std::fabs(probeP4.Eta()) < 2.5){
          mZeeAll_ptg50_endcap->Fill(TnP.M());
          if (selectedEle_[j].full5x5_sigmaIetaIeta() < 0.037)       mZeePass_ptg50_endcap->Fill(TnP.M());
          else                                                       mZeeFail_ptg50_endcap->Fill(TnP.M());
        }
      }
    }
    
  }
}
// ------------ method called once each job just before starting event loop  ------------
void 
ZeeTnP::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZeeTnP::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZeeTnP::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZeeTnP);
