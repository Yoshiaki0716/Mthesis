#include<iostream>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <dirent.h>
#include <vector>
#include <string>
#include <cstring>

void DiffScan(){
    TFile *file = TFile::Open("./20UPGM23610078/INITIAL_WARM/hists_20UPGM23610078_warm_021044_std_thresholdscan_hd.root","READ");
    TFile *file2 = TFile::Open("./20UPGM23610078/FINAL_WARM/hists_20UPGM23610078_warm_022890_std_thresholdscan_hd.root","READ");

    //TH1F*hist = (TH1F*)file->Get("h1d_0078_021043_OccupancyMap;1");
    //TH1F*hist2 = (TH1F*)file2->Get("h1d_0078_022889_OccupancyMap;1");
    TH2F*hist = (TH2F*)file->Get("map_0078_021044_ThresholdMap_0_0_300_0_0_0;1");
    TH2F*hist2 = (TH2F*)file2->Get("map_0078_022890_ThresholdMap_0_0_300_0_0_0;1");
    // 差分ヒストグラムを作成
    //TH1F* diffHist = (TH1F*)hist->Clone("diffHist");
    TH2F* diffHist = (TH2F*)hist->Clone("diffHist");

    diffHist->Add(hist2, -1);

    diffHist->SetStats(0); 
    //diffHist->GetXaxis()->SetRangeUser(-200, 200);

     // 1Dヒストグラムを作成（bin数や範囲は適宜調整）
    int nbins = 100;
    double min = diffHist->GetMinimum();
    double max = diffHist->GetMaximum();
    TH1F* projHist = new TH1F("projHist", "Difference Distribution", nbins, min, max);

    // 2Dヒストグラムの全ビンを走査して1DヒストグラムにFill
    for(int ix = 1; ix <= diffHist->GetNbinsX(); ++ix){
        for(int iy = 1; iy <= diffHist->GetNbinsY(); ++iy){
            double val = diffHist->GetBinContent(ix, iy);
            projHist->Fill(val);
        }
    }

    // 1Dヒストグラムを描画
    TCanvas* c2 = new TCanvas("c2", "Diff Value Histogram", 800, 600);
    c2->SetLogy();
    projHist->SetStats(0); 
    projHist->Draw();
    c2->SaveAs("DiffThrHist.pdf");

    TCanvas* c1 = new TCanvas("c1", "Occupancy Map Diff", 1000, 600);
    //diffHist->Draw();
    diffHist->Draw("colz");
    c1->Update();
    c1->SaveAs("DiffThrScan.pdf");

    file->Close();
    file2->Close();

}