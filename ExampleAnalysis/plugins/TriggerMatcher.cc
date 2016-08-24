#include <memory>
#include <TDirectory.h>
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
      ~TriggerMatcher() {}
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
	double hlt2reco_DeltaR;
	double hlt2reco_ptRatio;
	TTree* trigObjInfo_;
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
	hlt2reco_ptRatio=iConfig.getParameter<double>("PtRatio");

	useResouce("TFileService"); 
	edm::Service<TFileService> outFile;
	trigObjInfo_=outFile->make<TTree>("matchedTrigObjects","matchedTrigObjects");
	trigObjInfo_->Branch("Pt",trigObj_Pt[maxSize],"trigObj_Pt[maxSize]/D");
	trigObjInfo_->Branch("Eta",trigObj_Eta[maxSize],"trigObj_Eta[maxSize]/D");
	trigObjInfo_->Branch("Phi",trigObj_Phi[maxSize],"trigObj_Phi[maxSize]/D");
	trigObjInfo_->Branch("E",trigObj_E[maxSize],"trigObj_E[maxSize]/D");
	trigObjInfo_->Branch("DeltaR",trigObj_DeltaR[maxSize],"trigObj_DeltaR[maxSize]/D");
	trigObjInfo_->Branch("PtRatio",trigObj_PtRatio[maxSize],"trigObj_PtRatio[maxSize]/D");	
}

TriggerMatcher::~TriggerMatcher()
{

}

void
TriggerMatcher::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
		
	edm::Handle<pat::MuonCollection> muons_;
	iEvent.getByToken(muonToken_, muons)

	edm::Handle<pat::ElectronCollection> electrons;
	iEvent.getByToken(electronToken_, electrons);

	edm::Handle<edm::TriggerResults> triggerBits;
	iEvent.getByToken(triggerBits_, triggerBits);

	edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
	iEvent.getByToken(triggerObjects_, triggerObjects);
	
	const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
	
	index=names.triggerIndex(hltPathName);
	isPassed=triggerBits->accept(i);
	
	if(!isPassed){
		std::cout<<names.triggerNames(index)<<" is not fired"<<std::endl; 
		return;
	}	
	
	for (const pat::Electron &el : *electrons) {
		for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
			count++;	
			obj.unpackPathNames(names);
			ispath = obj.hasPathName(hltPathName);
			isFilter =  obj.hasFilterLabel(hltFilterName);		
			if(!isPath || !isFilter) continue;
			else{
				dEta=el.eta() - obj.eta();
				dPhi=TVector2::mpi_pi(el.phi() - obj.phi());
				dR = sqrt(pow(dEta,2)+pow(dPhi,2));
				ptRatio= obj.pt()/el.pt() ; 
				if(dR < hlt2reco_DeltaR)
				 trigObj_Pt[ii]= obj.pt(); trigObj_Eta[ii]=obj.eta(); trigObj_Phi[ii]=obj.phi(); trigObj_E[ii]=obj.E();
				 trigObj_DeltaR[ii]=dR; trigObj_PtRatio=ptRatio;		
			}
		}
	}	 
			
	for (const pat::muon &mu : *muons) {
        for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
            count++;    
            obj.unpackPathNames(names);
            ispath = obj.hasPathName(hltPathName);
            isFilter =  obj.hasFilterLabel(hltFilterName);      
            if(!isPath || !isFilter) continue;
            else{
                dEta=mu.eta() - obj.eta();
                dPhi=TVector2::mpi_pi(mu.phi() - obj.phi());
                dR = sqrt(pow(dEta,2)+pow(dPhi,2));
                ptRatio= obj.pt()/mu.pt() ; 
                if(dR < hlt2reco_DeltaR)
                 trigObj_Pt[ii]= obj.pt(); trigObj_Eta[ii]=obj.eta(); trigObj_Phi[ii]=obj.phi(); trigObj_E[ii]=obj.E();
                 trigObj_DeltaR[ii]=dR; trigObj_PtRatio=ptRatio;        
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
