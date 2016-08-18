#ifndef __ZTnP_h
#define __ZTnP_h
// -*- C++ -*-
//
// Package:    IndiaCMSTutorial/ZTnP
// Class:      ZTnP
// 
/**\class ZTnP ZTnP.cc IndiaCMSTutorial/ZTnP/plugins/ZTnP.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//
class ZTnP : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ZTnP(const edm::ParameterSet&);
      ~ZTnP();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
};
#endif
