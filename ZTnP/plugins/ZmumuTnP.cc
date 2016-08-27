#include "ZmumuTnP.h"

ZmumuTnP::ZmumuTnP(const edm::ParameterSet& iConfig) :
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects")))
{

   hltPathName=iConfig.getParameter<std::string>("hltPath");
   hltFilterName=iConfig.getParameter<std::string>("hltFilter"); 
   DeltaR_=iConfig.getParameter<double>("DeltaR");
	
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> outFile;
   outTree_ = outFile->make<TTree>("zmmtree","muTnP tree for iso");
   outTree_->Branch("nselectedMu",&nselectedMu,"nselectedMu/I");
   outTree_->Branch("nTnP",&nTnP,"nTnP/I");
   outTree_->Branch("TnP_pt",TnP_pt,"TnP_pt[nTnP]/F");
   //**Create the other branches in the tree here   

   histoDir = new TFileDirectory(outFile->mkdir("Histos_zmumu"));
   //**Book the mass histograms here
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
  triggerObj_.clear();
  nTnP = 0;
  
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  //vertices
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();
  //hnVtx->Fill(vertices->size());

  //muons
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  for (const pat::Muon &mu : *muons) {
    //**Apply muon cuts here and select muons
    selectedMu_.push_back(mu);
  }
  
  names = iEvent.triggerNames(*triggerBits);
  unsigned int index=names.triggerIndex(hltPathName);
  bool isPassed=triggerBits->accept(index);
  //std::cout<<hltPathName<<"\tindex= "<<index<<"\tPassed="<<isPassed<<std::endl;

  if(!isPassed){
	  //std::cout<<hltPathName<<" is not fired"<<std::endl;
	  return;
  }	  	
  
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
		
		triggerObj_.push_back(obj);			
  }
  //Select Zmumu candidates in selectZmumu function	
  if(selectedMu_.size() > 2)   selectZmumu();
  nselectedMu = selectedMu_.size();
  outTree_->Fill();

}
//tag selection
//pt > 20. and rel_iso03 < 0.35
//tag + probe pair = OS + 70 < mass < 110
void ZmumuTnP::selectZmumu() {
  for(unsigned int i = 0; i < selectedMu_.size(); i++) {
    double tagrelIso = mupfiso(selectedMu_[i])/selectedMu_[i].pt();  
    //**Apply tag selection here
    TLorentzVector tagP4 = getP4(selectedMu_[i]);
    int tagcharge =  selectedMu_[i].charge();  
    for(unsigned int j = 0; j < selectedMu_.size(); j++) {
      int probecharge = selectedMu_[j].charge();
      //**Apply probe selection 
      if(nTnP >= 8)    continue; 
      TLorentzVector probeP4 = getP4(selectedMu_[j]);
      TLorentzVector TnP = (tagP4 + probeP4);
      //TnP pair mass window cut
      if( TnP.M() <= 70. || TnP.M() >= 110. )   continue;
      double proberelIso = mupfiso(selectedMu_[j])/selectedMu_[j].pt();
      //Fill out the tree variables
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
          if(proberelIso < 0.35)   mZmmPass_ptl50_barrel->Fill(TnP.M());
          else                     mZmmFail_ptl50_barrel->Fill(TnP.M());
        } else {
          mZmmAll_ptl50_endcap->Fill(TnP.M());
          if(proberelIso < 0.35)   mZmmPass_ptl50_endcap->Fill(TnP.M());
          else                     mZmmFail_ptl50_endcap->Fill(TnP.M());
        }
      } else {
        if(std::fabs(probeP4.Eta()) < 1.2 ) {
          mZmmAll_ptg50_barrel->Fill(TnP.M());
          if(proberelIso < 0.35)   mZmmPass_ptg50_barrel->Fill(TnP.M());
          else                     mZmmFail_ptg50_barrel->Fill(TnP.M());
        } else {
          mZmmAll_ptg50_endcap->Fill(TnP.M());
          if(proberelIso < 0.35)   mZmmPass_ptg50_endcap->Fill(TnP.M());
          else                     mZmmFail_ptg50_endcap->Fill(TnP.M());
        }
      }
    }
  }
}

//function to compute muon isolation
//deltaB correction
double ZmumuTnP::mupfiso(const pat::Muon& mu) {
  return (mu.pfIsolationR03().sumChargedHadronPt 
          + std::max(0., mu.pfIsolationR03().sumNeutralHadronEt
                       + mu.pfIsolationR03().sumPhotonEt - 0.5 * mu.pfIsolationR03().sumPUPt));
}

// function for matching pat(reco) muon with its trigger object
bool ZmumuTnP::isMatchedtoTrigger (const pat::Muon& mu, double hlt2reco_deltaRmax){
	bool isMatch=false;
    
	for(unsigned int i = 0; i < triggerObj_.size(); i++){            
		triggerObj_[i].unpackPathNames(names);	
		isPath = triggerObj_[i].hasPathName(hltPathName);
		isFilter =  triggerObj_[i].hasFilterLabel(hltFilterName);      
		if(!isPath || !isFilter) continue;
		dEta= mu.eta() - triggerObj_[i].eta();
		dPhi=TVector2::Phi_mpi_pi(mu.phi() - triggerObj_[i].phi());
		dR = sqrt(pow(dEta,2)+pow(dPhi,2));
//		ptRatio= obj.pt()/mu.pt() ; 
		if(dR < hlt2reco_deltaRmax){
			isMatch=true;
			//std::cout<<"===== Kinematic info of Matched Trigger object ===="<<std::endl;
			//std::cout<<"pt= "<<triggerObj_[i].pt()<<"\tabs(eta)= "<<fabs(triggerObj_[i].eta())<<"\tphi= "<<triggerObj_[i].phi()<<std::endl;
			break;		
		}
	}
	return isMatch;
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
