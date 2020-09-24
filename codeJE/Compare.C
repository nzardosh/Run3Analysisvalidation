Int_t Compare(TString filerun3 = "AnalysisResults.root", TString filerun1 = "Vertices2prong-ITS1.root", bool donorm = false)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);

  TFile* fRun3 = new TFile(filerun3.Data());
  TFile* fRun1 = new TFile(filerun1.Data());

  const int nhisto = 10;
  TString histonameRun1[nhisto] = {"hjet_pt",
				   "hjet_phi",
				   "hjet_eta",
				   "hjet_n",
				   "hjet_zg",
				   "hjet_rg",
				   "hjet_nsd",
				   "hjet_TT_pt",
				   "hhadron_TT_pt",
				   "hhadronjet_TT_phi"};
                                   //"hImpParErr",
                                   //"hDecLenErr",
                                   //"hDecLenXYErr",
                                   //"hCovPVXX",
                                   //"hCovSVXX"
                                   //"hvx3",
                                   //"hvy3",
                                   //"hvz3"};
  TString histonameRun3[nhisto] = {"jet-finder/h_jet_pt",
				   "jet-finder/h_jet_phi",
				   "jet-finder/h_jet_eta",
				   "jet-finder/h_jet_n",
				   "jet-substructure/h_jet_zg",
				   "jet-substructure/h_jet_rg",
				   "jet-substructure/h_jet_nsd",
				   "jet-finder-hadron-recoil/h_jet_pt",
				   "jet-finder-hadron-recoil/h_hadron_pt",
				   "jet-finder-hadron-recoil/h_jet_hadron_deltaphi"};
                                   //"hf-task-d0/hImpParErr",
                                   //"hf-task-d0/hDecLenErr",
                                   //"hf-task-d0/hDecLenXYErr",
                                   //"hf-cand-creator-2prong/hCovPVXX",
                                   //"hf-cand-creator-2prong/hCovSVXX"
                                   //"hf-track-index-skims-creator/hvtx3_x",
                                   //"hf-track-index-skims-creator/hvtx3_y",
                                   //"hf-track-index-skims-creator/hvtx3_z"};
  TString xaxis[nhisto] = {"jet pT GeV/c",
			   "#phi",
			   "#eta",
			   "n constituents",
			   "zg",
			   "rg",
			   "nsd",
			   "jet pT GeV/c",
			   "trigger hadron pT GeV/c",
			   "#Delta #varphi"};
                           //"impact parameter error",
                           //"decay length error",
                           //"decay length XY error",
                           //"XX element of PV cov. matrix",
                           //"XX element of SV cov. matrix"
                           //"secondary vtx x - 3prong",
                           //"secondary vtx y - 3prong",
                           //"secondary vtx z - 3prong"};
  int rebin[nhisto] = {1,
		       1,
		       1,
		       1,
		       1,
		       1,
		       1,
		       1,
		       1,
		       1};
  bool LogScale[nhisto] = {true,
			  false,
			  false,
			   true,
			   false,
			   false,
			   false,
			   true,
			   true,
			   false};
  TH1F* hRun1[nhisto];
  TH1F* hRun3[nhisto];
  TH1F* hRatio[nhisto];
  for (int index = 0; index < nhisto; index++) {
    hRun1[index] = (TH1F*)fRun1->Get(histonameRun1[index].Data());
    if (!hRun1[index]) {
      printf("Failed to load %s from %s\n", histonameRun1[index].Data(), filerun1.Data());
      return 1;
    }
    hRun3[index] = (TH1F*)fRun3->Get(histonameRun3[index].Data());
    if (!hRun3[index]) {
      printf("Failed to load %s from %s\n", histonameRun3[index].Data(), filerun3.Data());
      return 1;
    }
  }

  // Histogram plot vertical margins
  Float_t marginHigh = 0.05;
  Float_t marginLow = 0.05;
  Float_t k = 1. - marginHigh - marginLow;
  // Ratio plot vertical margins
  Float_t marginRHigh = 0.05;
  Float_t marginRLow = 0.05;
  Float_t kR = 1. - marginRHigh - marginRLow;
  bool LogScaleR = false;
  Float_t yMin, yMax, yRange;
  Int_t nRun1, nRun3;

  TCanvas* cv = new TCanvas("cv", "Histos", 3000, 3000);
  cv->Divide(3, 5);
  TCanvas* cr = new TCanvas("cr", "Ratios", 3000, 3000);
  cr->Divide(3, 5);
 
  Int_t i=0;
  for (int index = 0; index < nhisto; index++) {
    nRun1 = hRun1[index]->GetEntries();
    nRun3 = hRun3[index]->GetEntries();
    auto pad = cv->cd(i + 1);
    if (donorm) {
      hRun1[index]->Scale(1./nRun1);
      hRun3[index]->Scale(1./nRun3);
    }
    hRun1[index]->Rebin(rebin[index]);
    hRun3[index]->Rebin(rebin[index]);
    hRun1[index]->SetLineColor(1);
    hRun1[index]->SetLineWidth(2);
    hRun3[index]->SetLineColor(2);
    hRun3[index]->SetLineWidth(1);
    hRun1[index]->SetTitle(Form("Entries: Run1: %d, Run3: %d;%s;Entries", nRun1, nRun3, xaxis[index].Data()));
    hRun1[index]->GetYaxis()->SetMaxDigits(3);
    yMin = TMath::Min(hRun3[index]->GetMinimum(0), hRun1[index]->GetMinimum(0));
    yMax = TMath::Max(hRun3[index]->GetMaximum(), hRun1[index]->GetMaximum());
    if (LogScale[index] && yMin > 0 && yMax > 0) {
      yRange = yMax/yMin;
      hRun1[index]->GetYaxis()->SetRangeUser(yMin/std::pow(yRange, marginLow/k), yMax * std::pow(yRange, marginHigh/k));
      pad->SetLogy();
    } else {
      yRange = yMax - yMin;
      hRun1[index]->GetYaxis()->SetRangeUser(yMin - marginLow/k * yRange, yMax + marginHigh/k * yRange);
    }
    hRun1[index]->Draw();
    hRun3[index]->Draw("same");
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hRun1[index], "Run1", "L");
    legend->AddEntry(hRun3[index], "Run3", "L");
    legend->Draw();
    auto padR = cr->cd(i + 1);
    hRatio[index] = (TH1F*)hRun3[index]->Clone(Form("hRatio%d", index));
    hRatio[index]->Divide(hRun1[index]);
    hRatio[index]->SetTitle(Form("Entries ratio: %g;%s;Run3/Run1", (double)nRun3/(double)nRun1, xaxis[index].Data()));
    yMin = hRatio[index]->GetMinimum(0);
    yMax = hRatio[index]->GetMaximum();
    if (LogScaleR && yMin > 0 && yMax > 0) {
      yRange = yMax/yMin;
      hRatio[index]->GetYaxis()->SetRangeUser(yMin/std::pow(yRange, marginLow/k), yMax * std::pow(yRange, marginHigh/k));
      padR->SetLogy();
    } else {
      yRange = yMax - yMin;
      hRatio[index]->GetYaxis()->SetRangeUser(yMin - marginRLow/kR * yRange, yMax + marginRHigh/kR * yRange);
    }
    hRatio[index]->Draw();
    i=i+1;
  }

  
  cv->SaveAs("comparison_histos_jets.pdf");
  cr->SaveAs("comparison_ratios_jets.pdf");
  return 0;
}
