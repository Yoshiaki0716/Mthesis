#include<iostream>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <dirent.h>
#include <vector>
#include <string>
#include <cstring>

// 部分一致でディレクトリまたはファイル名を探す関数
std::string FindNameWithKeyword(const char* dir, const char* keyword, bool searchDir = false) {
    std::vector<std::string> matches;
    DIR* dp = opendir(dir);
    if (!dp) {
        std::cerr << "ディレクトリが開けません: " << dir << std::endl;
        return "";
    }
    struct dirent* entry;
    while ((entry = readdir(dp)) != NULL) {
        if (strstr(entry->d_name, keyword)) {
            if (searchDir) {
                if (entry->d_type == DT_DIR && strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0)
                    matches.push_back(entry->d_name);
            } else {
                if (entry->d_type == DT_REG)
                    matches.push_back(entry->d_name);
            }
        }
    }
    closedir(dp);

    if (matches.size() == 0) {
        std::cerr << "エラー: " << dir << " に " << keyword << " を含む" << (searchDir ? "ディレクトリ" : "ファイル") << "が見つかりません。" << std::endl;
        return "";
    }
    if (matches.size() > 1) {
        std::cerr << "エラー: " << dir << " に " << keyword << " を含む" << (searchDir ? "ディレクトリ" : "ファイル") << "が複数見つかりました。" << std::endl;
        for (auto& m : matches) std::cerr << "  " << m << std::endl;
        return "";
    }
    return matches[0];
}

// メイン関数
void DiffScan(const char* moduleKeyword, const char* stage1Keyword, const char* stage2Keyword, const char* fileKeyword, const char* histKeyword1, const char* histKeyword2) {
    // 1. モジュールディレクトリを探す
    std::string moduleDir = FindNameWithKeyword(".", moduleKeyword, true);
    if (moduleDir.empty()) return;

    // 2. ステージディレクトリを探す
    std::string stage1Dir = std::string("./") + moduleDir + "/" + FindNameWithKeyword((std::string("./") + moduleDir).c_str(), stage1Keyword, true);
    std::string stage2Dir = std::string("./") + moduleDir + "/" + FindNameWithKeyword((std::string("./") + moduleDir).c_str(), stage2Keyword, true);
    if (stage1Dir.empty() || stage2Dir.empty()) return;

    // 3. ファイル名を探す
    std::string file1 = FindNameWithKeyword(stage1Dir.c_str(), fileKeyword, false);
    std::string file2 = FindNameWithKeyword(stage2Dir.c_str(), fileKeyword, false);
    if (file1.empty() || file2.empty()) return;

    std::string path1 = stage1Dir + "/" + file1;
    std::string path2 = stage2Dir + "/" + file2;

    // 4. ファイルを開く
    TFile *f1 = TFile::Open(path1.c_str(), "READ");
    TFile *f2 = TFile::Open(path2.c_str(), "READ");
    if (!f1 || !f2) {
        std::cerr << "ファイルが開けません。" << std::endl;
        return;
    }

    // 5. ヒストグラム名を部分一致で探す（2単語両方を含むものを探す）
    std::string hist1Name = "";
    std::string hist2Name = "";
    TIter next1(f1->GetListOfKeys());
    TKey *key1;
    while ((key1 = (TKey*)next1())) {
        std::string name = key1->GetName();
        if (name.find(histKeyword1) != std::string::npos && name.find(histKeyword2) != std::string::npos) {
            if (!hist1Name.empty()) {
                std::cerr << "エラー: " << path1 << " に " << histKeyword1 << " と " << histKeyword2 << " を含むヒストグラムが複数見つかりました。" << std::endl;
                return;
            }
            hist1Name = name;
        }
    }
    TIter next2(f2->GetListOfKeys());
    TKey *key2;
    while ((key2 = (TKey*)next2())) {
        std::string name = key2->GetName();
        if (name.find(histKeyword1) != std::string::npos && name.find(histKeyword2) != std::string::npos) {
            if (!hist2Name.empty()) {
                std::cerr << "エラー: " << path2 << " に " << histKeyword1 << " と " << histKeyword2 << " を含むヒストグラムが複数見つかりました。" << std::endl;
                return;
            }
            hist2Name = name;
        }
    }
    if (hist1Name.empty() || hist2Name.empty()) {
        std::cerr << "ヒストグラムが見つかりません。" << std::endl;
        return;
    }

    // 6. ヒストグラムを取得
    TH2F* hist1 = (TH2F*)f1->Get(hist1Name.c_str());
    TH2F* hist2 = (TH2F*)f2->Get(hist2Name.c_str());
    if (!hist1 || !hist2) {
        std::cerr << "ヒストグラムの取得に失敗しました。" << std::endl;
        return;
    }

    // 7. 差分ヒストグラム作成
    TH2F* diffHist = (TH2F*)hist1->Clone("diffHist");
    diffHist->Add(hist2, -1);
    diffHist->SetStats(0);

    // 8. 1Dヒストグラム作成
    int nbins = 100;
    double min = diffHist->GetMinimum();
    double max = diffHist->GetMaximum();
    TH1F* projHist = new TH1F("projHist", "Difference Distribution", nbins, min, max);

    for(int ix = 1; ix <= diffHist->GetNbinsX(); ++ix){
        for(int iy = 1; iy <= diffHist->GetNbinsY(); ++iy){
            double val = diffHist->GetBinContent(ix, iy);
            projHist->Fill(val);
        }
    }

    // 9. プロット
    std::string outDir = "./Output/";
    std::string pdfHistName = outDir + std::string(moduleKeyword) + std::string(fileKeyword) + "DiffHist.pdf";
    std::string pdfScanName = outDir + std::string(moduleKeyword) + std::string(fileKeyword) + "DiffScan.pdf";

    // Outputディレクトリがなければ作成
    system(("mkdir -p " + outDir).c_str());

    TCanvas* c2 = new TCanvas("c2", "Diff Value Histogram", 800, 600);
    c2->SetLogy();
    projHist->SetStats(0); 
    projHist->Draw();
    c2->SaveAs(pdfHistName.c_str());

    TCanvas* c1 = new TCanvas("c1", "Occupancy Map Diff", 1000, 600);
    diffHist->Draw("colz");
    c1->Update();
    c1->SaveAs(pdfScanName.c_str());

    f1->Close();
    f2->Close();
}

//使用例：root -l 'DiffScan.C("23610078", "INITIAL_WARM", "FINAL_WARM", "threshold", "map", "ThresholdMap")'