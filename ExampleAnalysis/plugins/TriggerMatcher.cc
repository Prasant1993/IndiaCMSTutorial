#include <memory>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector2.h>
#include <TLorentzVector.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


class TriggerMatcher : public edm::EDAnalyzer {
   public:
      explicit TriggerMatcher(const edm::ParameterSet&);
      ~TriggerMatcher();
	  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
	  virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	  virtual void endJob() override;

	edm::EDGetTokenT<pat::ElectronCollection> electronToken_;	      
	edm::EDGetTokenT<pat::MuonCollection> muonToken_;
	edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
	edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
	std::string hltPathName;
	std::string hltFilterName;
	double hlt2reco_DeltaR, dEta, dPhi, dR, ptRatio;
//	double hlt2reco_ptRatio;
	TFile* outFile;
	TTree* trigObjInfo_;

	unsigned int index, ii;
	bool isPassed, isPath, isFilter;
	const static Int_t maxSize = 10;

	int ele_size, mu_size;
	double trigEle_Pt[maxSize];
	double trigEle_Eta[maxSize];
	double trigEle_Phi[maxSize];
//	double trigEle_E[maxSize];
	double trigEle_DeltaR[maxSize];
	double trigEle_PtRatio[maxSize];

	double trigMu_Pt[maxSize];
	double trigMu_Eta[maxSize];
	double trigMu_Phi[maxSize];
	double trigMu_DeltaR[maxSize];
	double trigMu_PtRatio[maxSize];
//	double trigMu_E[maxSize];

};

TriggerMatcher::TriggerMatcher(const edm::ParameterSet& iConfig):
electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects")))
{
	hltPathName=iConfig.getParameter<std::string>("hltPath");
	hltFilterName=iConfig.getParameter<std::string>("hltFilter"); 
	hlt2reco_DeltaR=iConfig.getParameter<double>("DeltaR");
//	hlt2reco_ptRatio=iConfig.getParameter<double>("PtRatio");

//	usesResource("TFileService"); 
	outFile=new TFile("TriggerTree.root","RECREATE");
	outFile->cd();
	trigObjInfo_=new TTree("matchedTrigObjects","tree for matched Trigger Objects");

	trigObjInfo_->Branch("ele_size",&ele_size);	
	trigObjInfo_->Branch("trigEle_Pt",trigEle_Pt,"trigEle_Pt[ele_size]/D");
	trigObjInfo_->Branch("trigEle_Eta",trigEle_Eta,"trigEle_Eta[ele_size]/D");
	trigObjInfo_->Branch("trigEle_Phi",trigEle_Phi,"trigEle_Phi[ele_size]/D");
//	trigObjInfo_->Branch("Ele_E",trigEle_E[maxSize],"trigEle_E[maxSize]/D");
	trigObjInfo_->Branch("trigEle_DeltaR",trigEle_DeltaR,"trigEle_DeltaR[ele_size]/D");
	trigObjInfo_->Branch("trigEle_PtRatio",trigEle_PtRatio,"trigEle_PtRatio[ele_size]/D");

	trigObjInfo_->Branch("mu_size",&mu_size);	
	trigObjInfo_->Branch("trigMu_Pt",trigMu_Pt,"trigMu_Pt[mu_size]/D");
	trigObjInfo_->Branch("trigMu_Eta",trigMu_Eta,"trigMu_Eta[mu_size]/D");
	trigObjInfo_->Branch("trigMu_Phi",trigMu_Phi,"trigMu_Phi[mu_size]/D");
//	trigObjInfo_->Branch("Mu_E",trigMu_E[maxSize],"trigMu_E[maxSize]/D");
	trigObjInfo_->Branch("trigMu_DeltaR",trigMu_DeltaR,"trigMu_DeltaR[mu_size]/D");
	trigObjInfo_->Branch("trigMu_PtRatio",trigMu_PtRatio,"trigMu_PtRatio[mu_size]/D");
		
}

TriggerMatcher::~TriggerMatcher()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void
TriggerMatcher::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{	
	edm::Handle<pat::MuonCollection> muons;
	iEvent.getByToken(muonToken_, muons);

	edm::Handle<pat::ElectronCollection> electrons;
	iEvent.getByToken(electronToken_, electrons);

	edm::Handle<edm::TriggerResults> triggerBits;
	iEvent.getByToken(triggerBits_, triggerBits);

	edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
	iEvent.getByToken(triggerObjects_, triggerObjects);
	
	const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
	
	index=names.triggerIndex(hltPathName);
	isPassed=triggerBits->accept(index);

	std::cout<<hltPathName<<"\tindex= "<<index<<"\tPassed="<<isPassed<<std::endl;

	if(!isPassed){
		std::cout<<hltPathName<<" is not fired"<<std::endl; 
		return;
	}	
	
	ele_size=0;	
	for (const pat::Electron &el : *electrons){
		if (el.pt()<5.0) continue;		
		for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
			obj.unpackPathNames(names);
			isPath = obj.hasPathName(hltPathName);
			isFilter =  obj.hasFilterLabel(hltFilterName);	
			if(!isPath || !isFilter) continue;
			dEta=el.eta() - obj.eta();
			dPhi=TVector2:: Phi_mpi_pi(el.phi() - obj.phi());
			dR = sqrt(pow(dEta,2)+pow(dPhi,2));
			ptRatio= obj.pt()/el.pt() ; 
			if(dR < hlt2reco_DeltaR){
				trigEle_Pt[ele_size]= obj.pt(); trigEle_Eta[ele_size]=obj.eta(); trigEle_Phi[ele_size]=obj.phi(); //trigEle_E[ii]=obj.E();
				trigEle_DeltaR[ele_size]=dR; trigEle_PtRatio[ele_size]=ptRatio;
				ele_size++;
				break;
			} 		
		}
	}	 
	
	mu_size=0;		
	for (const pat::Muon &mu : *muons) {
		if(mu.pt()<5.0) continue; 		  
        for (pat::TriggerObjectStandAlone obj : *triggerObjects) {            
			obj.unpackPathNames(names);	
            isPath = obj.hasPathName(hltPathName);
            isFilter =  obj.hasFilterLabel(hltFilterName);      
            if(!isPath || !isFilter) continue;
            dEta=mu.eta() - obj.eta();
            dPhi=TVector2::Phi_mpi_pi(mu.phi() - obj.phi());
            dR = sqrt(pow(dEta,2)+pow(dPhi,2));
            ptRatio= obj.pt()/mu.pt() ; 
            if(dR < hlt2reco_DeltaR){
				trigMu_Pt[mu_size]= obj.pt(); trigMu_Eta[mu_size]=obj.eta(); trigMu_Phi[mu_size]=obj.phi(); //trigMu_E[ii]=obj.E();
				trigMu_DeltaR[mu_size]=dR; trigMu_PtRatio[mu_size]=ptRatio;
				mu_size++;
				break;	
			}	        
        }
	}
	
	trigObjInfo_->Fill();
}		
	
void 
TriggerMatcher::beginJob()
{
}

void 
TriggerMatcher::endJob() 
{
	outFile->cd();
	trigObjInfo_->Write();
	outFile->Close();
}

void
TriggerMatcher::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
// The following says we do not know what parameters are allowed so do no validation
// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerMatcher);
