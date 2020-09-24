#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TSystem.h>
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexerTracks.h"
#include "ReadJson.C"
#include "AliFJWrapper.h"
#include "FJ_includes.h"
#endif

Bool_t SingleTrkCuts(AliESDtrack* trk, AliESDtrackCuts* esdTrackCuts, AliESDVertex* fV1, Double_t fBzkG)
{
  if (!trk->PropagateToDCA(fV1, fBzkG, kVeryBig))
    return kFALSE;
  trk->RelateToVertex(fV1, fBzkG, kVeryBig);
  return esdTrackCuts->AcceptTrack(trk);
}

Bool_t SingleTrkCutsSimple(AliESDtrack* trk, Int_t minclutpc, int ptmintrack, double dcatoprimxymin, AliESDVertex* fV1, Double_t fBzkG)
{
  Int_t status = trk->GetStatus();
  bool sel_track = status & AliESDtrack::kITSrefit && (trk->HasPointOnITSLayer(0) || trk->HasPointOnITSLayer(1));
  sel_track = sel_track && trk->GetNcls(1) >= minclutpc;
  sel_track = sel_track && trk->Pt() > ptmintrack;
  //AliExternalTrackParam* track = (AliExternalTrackParam*)trk;
  //double b[2];
  //double bCov[3];
  //track->PropagateToDCA(fV1, fBzkG, 100., b, bCov);
  //sel_track = sel_track && abs(b[0]) > dcatoprimxymin;
  return sel_track;
}



AliAODVertex* ConvertToAODVertex(AliESDVertex* trkv)
{
  Double_t pos_[3], cov_[6], chi2perNDF_;
  trkv->GetXYZ(pos_);       // position
  trkv->GetCovMatrix(cov_); // covariance matrix
  chi2perNDF_ = trkv->GetChi2toNDF();
  double dispersion_ = trkv->GetDispersion();
  //  printf(" pos_ %f %f %f \n", pos_[0], pos_[1], pos_[2]);
  AliAODVertex* vertexAOD = new AliAODVertex(pos_, cov_, chi2perNDF_, 0x0, -1, AliAODVertex::kUndef, 2);
  return vertexAOD;
}

Float_t RelativePhi(Float_t phi1, Float_t phi2){
  if(phi1 < -TMath::Pi()) {
    phi1 += (2*TMath::Pi());
  }
  else if (phi1 > TMath::Pi()){
    phi1 -= (2*TMath::Pi());
  }
  if(phi2 < -TMath::Pi()) {
    phi2 += (2*TMath::Pi());
  }
  else if(phi2 > TMath::Pi()){
    phi2 -= (2*TMath::Pi());
  }
  Float_t deltaPhi=phi2-phi1;
  if(deltaPhi < -TMath::Pi()){
    deltaPhi += (2*TMath::Pi());
  }
  else if (deltaPhi > TMath::Pi()) {
    deltaPhi -= (2*TMath::Pi());
  }
  return deltaPhi;
}


Int_t JetFindingRun1(TString esdfile = "AliESDs.root",
                          TString output = "JetRun1.root",
                          TString jsonconfig = "dpl-config_std.json",
                          double ptmintrack = 0.15,
                          int do3Prongs = 0,
                          TString triggerstring = "")
{

  TFile* esdFile = TFile::Open(esdfile.Data());
  if (!esdFile || !esdFile->IsOpen()) {
    printf("Error in opening ESD file");
    return 1;
  }

  AliESDEvent* esd = new AliESDEvent;
  TTree* tree = (TTree*)esdFile->Get("esdTree");
  if (!tree) {
    printf("Error: no ESD tree found");
    return 1;
  }
  esd->ReadFromTree(tree);

  //int minncluTPC = 71;
  float minTTPt,maxTTPt;
  // read configuration from json file
  if (jsonconfig != "" && gSystem->Exec(Form("ls %s > /dev/null", jsonconfig.Data())) == 0) {
    printf("Read configuration from JSON file\n");
    minTTPt = GetJsonFloat(jsonconfig.Data(), "f_trackTTMin");
    printf("Min pt TT = %f\n", minTTPt);
    maxTTPt = GetJsonInteger(jsonconfig.Data(), "f_trackTTMax");
    printf("Max pt TT = %f\n", maxTTPt);
  }

  TH1F* hjet_pt = new TH1F("hjet_pt", " ; pt jet (#GeV) ; Entries", 100, 0, 100);
  TH1F* hjet_phi = new TH1F("hjet_phi", " ; #phi jet  ; Entries", 80, -1.0, 7.0);
  TH1F* hjet_eta = new TH1F("hjet_eta", " ; #eta jet ; Entries", 70, -0.7, 0.7);
  TH1F* hjet_n = new TH1F("hjet_n", " ; n jet constituents ; Entries", 30, 0, 30);
  TH1F* hjet_zg = new TH1F("hjet_zg", " ; zg ; Entries", 10, 0.0, 0.5);
  TH1F* hjet_rg = new TH1F("hjet_rg", " ; rg ; Entries", 10, 0.0, 0.5);
  TH1F* hjet_nsd = new TH1F("hjet_nsd", " ; nsd ; Entries", 7, -0.5, 6.5);
  TH1F* hjet_TT_pt = new TH1F("hjet_TT_pt", " ; pt jet recoil (#GeV) ; Entries", 100, 0., 100.);
  TH1F* hhadron_TT_pt = new TH1F("hhadron_TT_pt", " ; pt trigger hadron (#GeV) ; Entries", 120, 0., 60.);
  TH1F* hhadronjet_TT_phi = new TH1F("hhadronjet_TT_phi", " ; #delta #phi jet-trigger hadron; Entries", 40, 0.0, 4.);

  /*
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "default");
  esdTrackCuts->SetPtRange(ptmintrack, 1.e10);
  esdTrackCuts->SetEtaRange(-0.9, +0.9); //why?
  esdTrackCuts->SetMinNClustersTPC(minncluTPC);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                         AliESDtrackCuts::kAny);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetMaxDCAToVertexZ(3.2);
  esdTrackCuts->SetMaxDCAToVertexXY(2.4);
  esdTrackCuts->SetDCAToVertex2D(kTRUE);
  */
  

  for (Int_t iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
    tree->GetEvent(iEvent);
    if (!esd) {
      printf("Error: no ESD object found for event %d", iEvent);
      return 1;
    }
    TString trClass = esd->GetFiredTriggerClasses();
    if (triggerstring != "" && !trClass.Contains(triggerstring))
      //continue;
    printf("\n------------ Run %d Event: %d  Tracks %d ------------------\n",
           esd->GetRunNumber(), iEvent, esd->GetNumberOfTracks());
    printf("      Fired Trigger Classes %s\n", trClass.Data());

    Int_t maxTracksToProcess = 999999999; /// temporary to limit the time duration of tests
    Int_t totTracks = TMath::Min(maxTracksToProcess, esd->GetNumberOfTracks());

    AliESDVertex* primvtx = (AliESDVertex*)esd->GetPrimaryVertex();
    if (!primvtx)
      //return 1;
    TString title = primvtx->GetTitle();
    if (primvtx->IsFromVertexer3D() || primvtx->IsFromVertexerZ())
      //continue;
    if (primvtx->GetNContributors() < 2)
      //continue;
    AliAODVertex *vertexAODp = ConvertToAODVertex(primvtx);
    if (triggerstring != "" && !trClass.Contains(triggerstring))
      //continue;

    Double_t fBzkG = (Double_t)esd->GetMagneticField();

    AliFJWrapper *fFastJetWrapper;
    fFastJetWrapper = new AliFJWrapper("fFastJetWrapper","fFastJetWrapper");
    fFastJetWrapper->Clear();
    fFastJetWrapper->SetR(0.4); 
    fFastJetWrapper->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
    fFastJetWrapper->SetRecombScheme(fastjet::RecombinationScheme::E_scheme);
    fFastJetWrapper->SetStrategy(fastjet::Strategy::Best);
    fFastJetWrapper->SetGhostArea(0.005); 
    fFastJetWrapper->SetAreaType(fastjet::AreaType::passive_area);

    std::vector<int> trackTTIndex;
    trackTTIndex.clear();

    Float_t pTComparison=0.0;
    Float_t PhiTT;
    Float_t PtTT;
    bool isTT=false;
    for (Int_t iTrack = 0; iTrack < totTracks; iTrack++) {
      AliESDtrack* track = esd->GetTrack(iTrack);
      //if (!SingleTrkCutsSimple(track, minncluTPC, ptmintrack, 0.03, primvtx, fBzkG)) continue;
      if (track->Pt() < 0.15 || TMath::Abs(track->Eta()) >= 0.9) continue;
      fFastJetWrapper->AddInputVector(track->Px(), track->Py(), track->Pz(), TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack+100);
      if (track->Pt() >= minTTPt && track->Pt() < maxTTPt){
	isTT=true;
	if (track->Pt() >= pTComparison){
	  pTComparison=track->Pt();
	  PtTT=track->Pt();
	  PhiTT=track->Phi();
	}
      }
    }
    if(isTT) hhadron_TT_pt->Fill(PtTT);
    fFastJetWrapper->Run();
    std::vector<fastjet::PseudoJet> jets = fFastJetWrapper->GetInclusiveJets();
    for (Int_t ijet=0; ijet<jets.size(); ijet++){
      fastjet::PseudoJet jet = jets[ijet];
      if (jet.pt() <= 0.15 || jet.perp() > 1000.0 || TMath::Abs(jet.eta()) >= 0.5) continue;
      hjet_pt->Fill(jet.pt());
      hjet_phi->Fill(jet.phi());
      hjet_eta->Fill(jet.eta());
      if (isTT){
	Float_t deltaPhi=TMath::Abs(RelativePhi(jet.phi(),PhiTT));
	if (deltaPhi >= (TMath::Pi() - 0.6)){
	hjet_TT_pt->Fill(jet.pt());
	}
	if(deltaPhi >= TMath::Pi()/2.0 && deltaPhi <= TMath::Pi()){
	  hhadronjet_TT_phi->Fill(deltaPhi);
	}
      }
      std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(ijet));
      hjet_n->Fill(constituents.size());
      fastjet::JetDefinition subJetDef(fastjet::JetAlgorithm::cambridge_algorithm , 0.4*2.5,fastjet::RecombinationScheme::E_scheme, fastjet::Best);
      try{
	fastjet::ClusterSequence reclusterSeq(constituents, subJetDef);
	std::vector<fastjet::PseudoJet> reclusteredJet =  reclusterSeq.inclusive_jets(0.0);
	reclusteredJet = sorted_by_pt(reclusteredJet);
         
	fastjet::PseudoJet daughterSubJet = reclusteredJet[0];
	fastjet::PseudoJet parentSubJet1; 
	fastjet::PseudoJet parentSubJet2;

	Float_t zg=-1.0,rg=-1.0;
	Int_t nsd=0;
	bool softDropped=false;
	while(daughterSubJet.has_parents(parentSubJet1,parentSubJet2)){
	  if(parentSubJet1.perp() < parentSubJet2.perp()) std::swap(parentSubJet1,parentSubJet2);
	  zg=parentSubJet2.perp()/(parentSubJet1.perp()+parentSubJet2.perp());
	  rg=parentSubJet1.delta_R(parentSubJet2);

	  if (zg >= 0.1*TMath::Power(rg/0.4,0.0)){
	    if(!softDropped){
	      hjet_zg->Fill(zg);
	      hjet_rg->Fill(rg);
	      softDropped=true;
	    }
	    nsd++;
	  }
	  daughterSubJet=parentSubJet1;
	}
	hjet_nsd->Fill(nsd);
      }catch (fastjet::Error) {}
    }
    delete fFastJetWrapper;

  }
  //delete esdTrackCuts;
  
  TFile* fout = new TFile(output.Data(), "recreate");
  hjet_pt->Write();
  hjet_phi->Write();
  hjet_eta->Write();
  hjet_n->Write();
  hjet_zg->Write();
  hjet_rg->Write();
  hjet_nsd->Write();
  hjet_TT_pt->Write();
  hhadron_TT_pt->Write();
  hhadronjet_TT_phi->Write();
  

  fout->Close();
  return 0;
}
