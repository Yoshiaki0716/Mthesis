
////_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
//
// Understanding primitve variables.
//
//_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
#include "utils.h"
#include "LandauGauss.h"

void wakuATLAS(double xmin, double xmax, double ymin, double ymax, std::string xtit, std::string ytit) {
  TH1F *waku = new TH1F("","",1,xmin,xmax);
  waku -> SetMaximum(ymax);
  waku -> SetMinimum(ymin);
  waku -> GetXaxis() -> SetTitle(xtit.c_str());
  waku -> GetYaxis() -> SetTitle(ytit.c_str());
  waku -> DrawCopy();
  waku -> Delete();
}

void Hist1DPlot(int imode, int icl1, int ist1, double rnorm1, TH1F *histA, double xmin, double xmax, std::string xtit, std::string ytit) {

  if (imode>1) { histA -> Rebin(imode); }

  int NbinX    = histA->GetNbinsX();
  if (xmin == xmax) { 
    xmin = histA->GetXaxis()->GetXmin();
    xmax = histA->GetXaxis()->GetXmax(); 
  }
  double xWidth = histA->GetBinWidth(1);
  double totA = 0.0;
  double ymxA = 0.0;
  double ymnA = 1.0e10;

  double xminc = histA->GetXaxis()->GetXmin();
  for (int i = 0; i < NbinX; i++) {
    if (xminc+xWidth*i >= xmin && xminc+xWidth*i <= xmax) {
      double hn1 = histA->GetBinContent(i+1);
      totA += hn1;
      if (ymxA < hn1) { ymxA = hn1; }
      if (ymnA > hn1 && hn1 > 0.0) { ymnA = hn1; }
    }
  }

  double sumA = totA;
  if (sumA == 0.0) { sumA = 1.0; }

  double ymax  = ymxA/sumA*rnorm1;
  double ymin  = ymnA/sumA*rnorm1;

  TH1F *waku = new TH1F("","",1,xmin,xmax);
  if (imode%2 == 0) { waku -> SetMaximum(1.2*ymax); }
  if (imode%2 == 1) { waku -> SetMaximum(2.0*ymax); }
  if (imode%2 == 0) { waku -> SetMinimum(0.0); }
  if (imode%2 == 1) { waku -> SetMinimum(ymin); }
  waku -> GetXaxis() -> SetTitle(xtit.c_str());
  waku -> GetYaxis() -> SetTitle(ytit.c_str());
  waku -> DrawCopy(); 

  histA -> SetNormFactor(rnorm1);
  histA -> SetLineColor(icl1);
  histA -> SetLineStyle(ist1);
  histA -> SetLineWidth(2);
  histA -> Draw("HIST same");
  waku -> Delete();
}

void Hist1DPlot(int imode, int icl1, int icl2, int ist1, int ist2,
		double rnorm1, double rnorm2, TH1F *histA, TH1F *histB, 
		double xmin, double xmax, double Percentage, std::string xtit, std::string ytit) {

  if (imode>1) { histA -> Rebin(imode); }
  if (imode>1) { histB -> Rebin(imode); }

  int NbinX    = histA->GetNbinsX();
  if (xmin == xmax) { 
    xmin = histA->GetXaxis()->GetXmin();
    xmax = histA->GetXaxis()->GetXmax(); 
  }
  double xWidth = histA->GetBinWidth(1);
  double totA = 0.0;
  double totB = 0.0;
  double ymxA = 0.0;
  double ymxB = 0.0;
  double ymnA = 1.0e10;
  double ymnB = 1.0e10;

  double xminc = histA->GetXaxis()->GetXmin();
  for (int i = 0; i < NbinX; i++) {
    if (xminc+xWidth*i >= xmin && xminc+xWidth*i <= xmax) {
      double hn1 = histA->GetBinContent(i+1);
      double hn2 = histB->GetBinContent(i+1);
      totA += hn1;
      totB += hn2;
      if (ymxA < hn1) { ymxA = hn1; }
      if (ymxB < hn2) { ymxB = hn2; }
      if (ymnA > hn1 && hn1 > 0.0) { ymnA = hn1; }
      if (ymnB > hn2 && hn2 > 0.0) { ymnB = hn2; }
    }
  }

  double sumA = totA;
  double sumB = totB;
  if (sumA == 0.0) { sumA = 1.0; }
  if (sumB == 0.0) { sumB = 1.0; }

  double ymax  = TMath::Max(ymxA/sumA*rnorm1,ymxB/sumB*rnorm2);
  double ymin  = TMath::Min(ymnA/sumA*rnorm1,ymnB/sumB*rnorm2);

  TH1F *waku = new TH1F("","",1,xmin,xmax);
  if (imode%2 == 0) { waku -> SetMaximum(1.2*ymax); }
  if (imode%2 == 1) { waku -> SetMaximum(2.0*ymax); }
  if (imode%2 == 0) { waku -> SetMinimum(0.0); }
  if (imode%2 == 1) { waku -> SetMinimum(ymin); }
  waku -> GetXaxis() -> SetTitle(xtit.c_str());
  waku -> GetYaxis() -> SetTitle(ytit.c_str());
  waku -> DrawCopy(); 

  //histA -> SetNormFactor(rnorm1);
  //histB -> SetNormFactor(rnorm2);//見た目だけ変わる
  histA -> Scale(1/sumA);
  histB -> Scale(1/sumB);//値も変える？


  histA -> SetLineColor(icl1);
  histB -> SetLineColor(icl2);

  histA -> SetLineStyle(ist1);
  histB -> SetLineStyle(ist2);

  histA -> SetLineWidth(2);
  histB -> SetLineWidth(2);

  histB -> Draw("HIST same");
  histA -> Draw("HIST same");

  histA->SetName("histA");
  histB->SetName("histB");

  LangauFit fitBase_B;
  fitBase_B.set_initialParameter_landauMPV(0.6);
  fitBase_B.set_initialParameter_landauSigma(0.12);
  fitBase_B.set_initialParameter_totalArea(histB->GetEntries()/100);
  fitBase_B.set_fitPercentage(Percentage);
  fitBase_B.execute(histB,0);
  fitBase_B.get_fitFunc()->SetLineColor(kGreen);
  fitBase_B.get_fitFunc()->SetLineWidth(3);
  fitBase_B.get_fitFunc()->Draw("same");


  LangauFit fitBase_A;
  fitBase_A.set_initialParameter_landauMPV(0.6);
  fitBase_A.set_initialParameter_landauSigma(0.12);
  fitBase_A.set_initialParameter_totalArea(histA->GetEntries()/100);
  fitBase_A.set_fitPercentage(Percentage);
  fitBase_A.execute(histA,0);
  fitBase_A.get_fitFunc()->SetLineColor(kRed);
  fitBase_A.get_fitFunc()->SetLineWidth(3);
  fitBase_A.get_fitFunc()->Draw("same");



  waku -> Delete();
}

void Hist1DPlot(int imode, 
		int icl1, int icl2, int icl3,
		int ist1, int ist2, int ist3,
		double rnorm1, double rnorm2, double rnorm3,
		TH1F *histA, TH1F *histB, TH1F *histC,
		double xmin, double xmax,
		std::string xtit, std::string ytit) {

  if (imode>1) { histA -> Rebin(imode); }
  if (imode>1) { histB -> Rebin(imode); }
  if (imode>1) { histC -> Rebin(imode); }

  int NbinX    = histA->GetNbinsX();
  if (xmin == xmax) { 
    xmin = histA->GetXaxis()->GetXmin();
    xmax = histA->GetXaxis()->GetXmax(); 
  }
  double xWidth = histA->GetBinWidth(1);
  double totA = 0.0;
  double totB = 0.0;
  double totC = 0.0;
  double ymxA = 0.0;
  double ymxB = 0.0;
  double ymxC = 0.0;
  double ymnA = 1.0e10;
  double ymnB = 1.0e10;
  double ymnC = 1.0e10;

  double xminc = histA->GetXaxis()->GetXmin();
  for (int i = 0; i < NbinX; i++) {
    if (xminc+xWidth*i >= xmin && xminc+xWidth*i <= xmax) {
      double hn1 = histA->GetBinContent(i+1);
      double hn2 = histB->GetBinContent(i+1);
      double hn3 = histC->GetBinContent(i+1);
      totA += hn1;
      totB += hn2;
      totC += hn3;
      if (ymxA < hn1) { ymxA = hn1; }
      if (ymxB < hn2) { ymxB = hn2; }
      if (ymxC < hn3) { ymxC = hn3; }
      if (ymnA > hn1 && hn1 > 0.0) { ymnA = hn1; }
      if (ymnB > hn2 && hn2 > 0.0) { ymnB = hn2; }
      if (ymnC > hn3 && hn3 > 0.0) { ymnC = hn3; }
    }
  }

  double sumA = totA;
  double sumB = totB;
  double sumC = totC;
  if (sumA == 0.0) { sumA = 1.0; }
  if (sumB == 0.0) { sumB = 1.0; }
  if (sumC == 0.0) { sumC = 1.0; }

  double ymax  = TMath::Max(ymxA/sumA*rnorm1,ymxB/sumB*rnorm2);
  ymax = TMath::Max(ymax,ymxC/sumC*rnorm3);

  double ymin  = TMath::Min(ymnA/sumA*rnorm1,ymnB/sumB*rnorm2);
  ymin = TMath::Min(ymin,ymnC/sumC*rnorm3);

  TH1F *waku = new TH1F("","",1,xmin,xmax);
  if (imode%2 == 0) { waku -> SetMaximum(1.2*ymax); }
  if (imode%2 == 1) { waku -> SetMaximum(2.0*ymax); }
  if (imode%2 == 0) { waku -> SetMinimum(0.0); }
  if (imode%2 == 1) { waku -> SetMinimum(ymin); }
  waku -> GetXaxis() -> SetTitle(xtit.c_str());
  waku -> GetYaxis() -> SetTitle(ytit.c_str());
  waku -> DrawCopy(); 

  histA -> SetNormFactor(rnorm1);
  histB -> SetNormFactor(rnorm2);
  histC -> SetNormFactor(rnorm3);

  histA -> SetLineColor(icl1);
  histB -> SetLineColor(icl2);
  histC -> SetLineColor(icl3);

  histA -> SetLineStyle(ist1);
  histB -> SetLineStyle(ist2);
  histC -> SetLineStyle(ist3);

  histA -> SetLineWidth(2);
  histB -> SetLineWidth(2);
  histC -> SetLineWidth(2);

  histC -> Draw("HIST same");
  histB -> Draw("HIST same");
  histA -> Draw("HIST same");
  waku -> Delete();
}

void Hist1DPlot(int imode, 
		int icl1, int icl2, int icl3, int icl4,
		int ist1, int ist2, int ist3, int ist4,
		double rnorm1, double rnorm2, 
		double rnorm3, double rnorm4, 
		TH1F *histA, TH1F *histB, TH1F *histC, TH1F *histD, 
		double xmin, double xmax,
		std::string xtit, std::string ytit) {

  if (imode>1) { histA -> Rebin(imode); }
  if (imode>1) { histB -> Rebin(imode); }
  if (imode>1) { histC -> Rebin(imode); }
  if (imode>1) { histD -> Rebin(imode); }

  int NbinX    = histA->GetNbinsX();
  if (xmin == xmax) { 
    xmin = histA->GetXaxis()->GetXmin();
    xmax = histA->GetXaxis()->GetXmax(); 
  }
  double xWidth = histA->GetBinWidth(1);
  double totA = 0.0;
  double totB = 0.0;
  double totC = 0.0;
  double totD = 0.0;
  double ymxA = 0.0;
  double ymxB = 0.0;
  double ymxC = 0.0;
  double ymxD = 0.0;
  double ymnA = 1.0e10;
  double ymnB = 1.0e10;
  double ymnC = 1.0e10;
  double ymnD = 1.0e10;

  double xminc = histA->GetXaxis()->GetXmin();
  for (int i = 0; i < NbinX; i++) {
    if (xminc+xWidth*i >= xmin && xminc+xWidth*i <= xmax) {
      double hn1 = histA->GetBinContent(i+1);
      double hn2 = histB->GetBinContent(i+1);
      double hn3 = histC->GetBinContent(i+1);
      double hn4 = histD->GetBinContent(i+1);
      totA += hn1;
      totB += hn2;
      totC += hn3;
      totD += hn4;
      if (ymxA < hn1) { ymxA = hn1; }
      if (ymxB < hn2) { ymxB = hn2; }
      if (ymxC < hn3) { ymxC = hn3; }
      if (ymxD < hn4) { ymxD = hn4; }
      if (ymnA > hn1 && hn1 > 0.0) { ymnA = hn1; }
      if (ymnB > hn2 && hn2 > 0.0) { ymnB = hn2; }
      if (ymnC > hn3 && hn3 > 0.0) { ymnC = hn3; }
      if (ymnD > hn4 && hn4 > 0.0) { ymnD = hn4; }
    }
  }

  double sumA = totA;
  double sumB = totB;
  double sumC = totC;
  double sumD = totD;
  if (sumA == 0.0) { sumA = 1.0; }
  if (sumB == 0.0) { sumB = 1.0; }
  if (sumC == 0.0) { sumC = 1.0; }
  if (sumD == 0.0) { sumD = 1.0; }

  double ymax  = TMath::Max(ymxA/sumA*rnorm1,ymxB/sumB*rnorm2);
  ymax = TMath::Max(ymax,ymxC/sumC*rnorm3);
  ymax = TMath::Max(ymax,ymxD/sumD*rnorm4);

  double ymin  = TMath::Min(ymnA/sumA*rnorm1,ymnB/sumB*rnorm2);
  ymin = TMath::Min(ymin,ymnC/sumC*rnorm3);
  ymin = TMath::Min(ymin,ymnD/sumD*rnorm4);

  TH1F *waku = new TH1F("","",1,xmin,xmax);
  if (imode%2 == 0) { waku -> SetMaximum(1.2*ymax); }
  if (imode%2 == 1) { waku -> SetMaximum(2.0*ymax); }
  if (imode%2 == 0) { waku -> SetMinimum(0.0); }
  if (imode%2 == 1) { waku -> SetMinimum(ymin); }
  waku -> GetXaxis() -> SetTitle(xtit.c_str());
  waku -> GetYaxis() -> SetTitle(ytit.c_str());
  waku -> DrawCopy(); 

  histA -> SetNormFactor(rnorm1);
  histB -> SetNormFactor(rnorm2);
  histC -> SetNormFactor(rnorm3);
  histD -> SetNormFactor(rnorm4);

  histA -> SetLineColor(icl1);
  histB -> SetLineColor(icl2);
  histC -> SetLineColor(icl3);
  histD -> SetLineColor(icl4);

  histA -> SetLineStyle(ist1);
  histB -> SetLineStyle(ist2);
  histC -> SetLineStyle(ist3);
  histD -> SetLineStyle(ist4);

  histA -> SetLineWidth(2);
  histB -> SetLineWidth(2);
  histC -> SetLineWidth(2);
  histD -> SetLineWidth(2);

  histD -> Draw("HIST same");
  histC -> Draw("HIST same");
  histB -> Draw("HIST same");
  histA -> Draw("HIST same");
  waku -> Delete();
}

void Hist1DPlot(int imode, 
		int icl1, int icl2, int icl3, 
		int icl4, int icl5, 
		int ist1, int ist2, int ist3, 
		int ist4, int ist5, 
		double rnorm1, double rnorm2, double rnorm3, 
		double rnorm4, double rnorm5, 
		TH1F *histA, TH1F *histB, TH1F *histC, 
		TH1F *histD, TH1F *histE, 
		double xmin, double xmax,
		std::string xtit, std::string ytit) {

  if (imode>1) { histA -> Rebin(imode); }
  if (imode>1) { histB -> Rebin(imode); }
  if (imode>1) { histC -> Rebin(imode); }
  if (imode>1) { histD -> Rebin(imode); }
  if (imode>1) { histE -> Rebin(imode); }

  int NbinX    = histA->GetNbinsX();
  if (xmin == xmax) { 
    xmin = histA->GetXaxis()->GetXmin();
    xmax = histA->GetXaxis()->GetXmax(); 
  }
  double xWidth = histA->GetBinWidth(1);
  double totA = 0.0;
  double totB = 0.0;
  double totC = 0.0;
  double totD = 0.0;
  double totE = 0.0;
  double ymxA = 0.0;
  double ymxB = 0.0;
  double ymxC = 0.0;
  double ymxD = 0.0;
  double ymxE = 0.0;
  double ymnA = 1.0e10;
  double ymnB = 1.0e10;
  double ymnC = 1.0e10;
  double ymnD = 1.0e10;
  double ymnE = 1.0e10;

  double xminc = histA->GetXaxis()->GetXmin();
  for (int i = 0; i < NbinX; i++) {
    if (xminc+xWidth*i >= xmin && xminc+xWidth*i <= xmax) {
      double hn1 = histA->GetBinContent(i+1);
      double hn2 = histB->GetBinContent(i+1);
      double hn3 = histC->GetBinContent(i+1);
      double hn4 = histD->GetBinContent(i+1);
      double hn5 = histE->GetBinContent(i+1);
      totA += hn1;
      totB += hn2;
      totC += hn3;
      totD += hn4;
      totE += hn5;
      if (ymxA < hn1) { ymxA = hn1; }
      if (ymxB < hn2) { ymxB = hn2; }
      if (ymxC < hn3) { ymxC = hn3; }
      if (ymxD < hn4) { ymxD = hn4; }
      if (ymxE < hn5) { ymxE = hn5; }
      if (ymnA > hn1 && hn1 > 0.0) { ymnA = hn1; }
      if (ymnB > hn2 && hn2 > 0.0) { ymnB = hn2; }
      if (ymnC > hn3 && hn3 > 0.0) { ymnC = hn3; }
      if (ymnD > hn4 && hn4 > 0.0) { ymnD = hn4; }
      if (ymnE > hn5 && hn5 > 0.0) { ymnE = hn5; }
    }
  }

  double sumA = totA;
  double sumB = totB;
  double sumC = totC;
  double sumD = totD;
  double sumE = totE;
  if (sumA == 0.0) { sumA = 1.0; }
  if (sumB == 0.0) { sumB = 1.0; }
  if (sumC == 0.0) { sumC = 1.0; }
  if (sumD == 0.0) { sumD = 1.0; }
  if (sumE == 0.0) { sumE = 1.0; }

  double ymax  = TMath::Max(ymxA/sumA*rnorm1,ymxB/sumB*rnorm2);
  ymax = TMath::Max(ymax,ymxC/sumC*rnorm3);
  ymax = TMath::Max(ymax,ymxD/sumD*rnorm4);
  ymax = TMath::Max(ymax,ymxE/sumE*rnorm5);

  double ymin  = TMath::Min(ymnA/sumA*rnorm1,ymnB/sumB*rnorm2);
  ymin = TMath::Min(ymin,ymnC/sumC*rnorm3);
  ymin = TMath::Min(ymin,ymnD/sumD*rnorm4);
  ymin = TMath::Min(ymin,ymnE/sumE*rnorm5);

  TH1F *waku = new TH1F("","",1,xmin,xmax);
  if (imode%2 == 0) { waku -> SetMaximum(1.2*ymax); }
  if (imode%2 == 1) { waku -> SetMaximum(2.0*ymax); }
  if (imode%2 == 0) { waku -> SetMinimum(0.0); }
  if (imode%2 == 1) { waku -> SetMinimum(ymin); }
  waku -> GetXaxis() -> SetTitle(xtit.c_str());
  waku -> GetYaxis() -> SetTitle(ytit.c_str());
  waku -> DrawCopy(); 

  histA -> SetNormFactor(rnorm1);
  histB -> SetNormFactor(rnorm2);
  histC -> SetNormFactor(rnorm3);
  histD -> SetNormFactor(rnorm4);
  histE -> SetNormFactor(rnorm5);

  histA -> SetLineColor(icl1);
  histB -> SetLineColor(icl2);
  histC -> SetLineColor(icl3);
  histD -> SetLineColor(icl4);
  histE -> SetLineColor(icl5);

  histA -> SetLineStyle(ist1);
  histB -> SetLineStyle(ist2);
  histC -> SetLineStyle(ist3);
  histD -> SetLineStyle(ist4);
  histE -> SetLineStyle(ist5);

  histA -> SetLineWidth(2);
  histB -> SetLineWidth(2);
  histC -> SetLineWidth(2);
  histD -> SetLineWidth(2);
  histE -> SetLineWidth(2);

  histE -> Draw("HIST same");
  histD -> Draw("HIST same");
  histC -> Draw("HIST same");
  histB -> Draw("HIST same");
  histA -> Draw("HIST same");
  waku -> Delete();
}

void Hist1DPlot(int imode, 
		int icl1, int icl2, int icl3, 
		int icl4, int icl5, int icl6, 
		int ist1, int ist2, int ist3, 
		int ist4, int ist5, int ist6, 
		double rnorm1, double rnorm2, double rnorm3, 
		double rnorm4, double rnorm5, double rnorm6, 
		TH1F *histA, TH1F *histB, TH1F *histC, 
		TH1F *histD, TH1F *histE, TH1F *histF, 
		double xmin, double xmax,
		std::string xtit, std::string ytit) {

  if (imode>1) { histA -> Rebin(imode); }
  if (imode>1) { histB -> Rebin(imode); }
  if (imode>1) { histC -> Rebin(imode); }
  if (imode>1) { histD -> Rebin(imode); }
  if (imode>1) { histE -> Rebin(imode); }
  if (imode>1) { histF -> Rebin(imode); }

  int NbinX    = histA->GetNbinsX();
  if (xmin == xmax) { 
    xmin = histA->GetXaxis()->GetXmin();
    xmax = histA->GetXaxis()->GetXmax(); 
  }
  double xWidth = histA->GetBinWidth(1);
  double totA = 0.0;
  double totB = 0.0;
  double totC = 0.0;
  double totD = 0.0;
  double totE = 0.0;
  double totF = 0.0;
  double ymxA = 0.0;
  double ymxB = 0.0;
  double ymxC = 0.0;
  double ymxD = 0.0;
  double ymxE = 0.0;
  double ymxF = 0.0;
  double ymnA = 1.0e10;
  double ymnB = 1.0e10;
  double ymnC = 1.0e10;
  double ymnD = 1.0e10;
  double ymnE = 1.0e10;
  double ymnF = 1.0e10;

  double xminc = histA->GetXaxis()->GetXmin();
  for (int i = 0; i < NbinX; i++) {
    if (xminc+xWidth*i >= xmin && xminc+xWidth*i <= xmax) {
      double hn1 = histA->GetBinContent(i+1);
      double hn2 = histB->GetBinContent(i+1);
      double hn3 = histC->GetBinContent(i+1);
      double hn4 = histD->GetBinContent(i+1);
      double hn5 = histE->GetBinContent(i+1);
      double hn6 = histF->GetBinContent(i+1);
      totA += hn1;
      totB += hn2;
      totC += hn3;
      totD += hn4;
      totE += hn5;
      totF += hn6;
      if (ymxA < hn1) { ymxA = hn1; }
      if (ymxB < hn2) { ymxB = hn2; }
      if (ymxC < hn3) { ymxC = hn3; }
      if (ymxD < hn4) { ymxD = hn4; }
      if (ymxE < hn5) { ymxE = hn5; }
      if (ymxF < hn6) { ymxF = hn6; }
      if (ymnA > hn1 && hn1 > 0.0) { ymnA = hn1; }
      if (ymnB > hn2 && hn2 > 0.0) { ymnB = hn2; }
      if (ymnC > hn3 && hn3 > 0.0) { ymnC = hn3; }
      if (ymnD > hn4 && hn4 > 0.0) { ymnD = hn4; }
      if (ymnE > hn5 && hn5 > 0.0) { ymnE = hn5; }
      if (ymnF > hn6 && hn6 > 0.0) { ymnF = hn6; }
    }
  }

  double sumA = totA;
  double sumB = totB;
  double sumC = totC;
  double sumD = totD;
  double sumE = totE;
  double sumF = totF;
  if (sumA == 0.0) { sumA = 1.0; }
  if (sumB == 0.0) { sumB = 1.0; }
  if (sumC == 0.0) { sumC = 1.0; }
  if (sumD == 0.0) { sumD = 1.0; }
  if (sumE == 0.0) { sumE = 1.0; }
  if (sumF == 0.0) { sumF = 1.0; }

  double ymax  = TMath::Max(ymxA/sumA*rnorm1,ymxB/sumB*rnorm2);
  ymax = TMath::Max(ymax,ymxC/sumC*rnorm3);
  ymax = TMath::Max(ymax,ymxD/sumD*rnorm4);
  ymax = TMath::Max(ymax,ymxE/sumE*rnorm5);
  ymax = TMath::Max(ymax,ymxF/sumF*rnorm6);

  double ymin  = TMath::Min(ymnA/sumA*rnorm1,ymnB/sumB*rnorm2);
  ymin = TMath::Min(ymin,ymnC/sumC*rnorm3);
  ymin = TMath::Min(ymin,ymnD/sumD*rnorm4);
  ymin = TMath::Min(ymin,ymnE/sumE*rnorm5);
  ymin = TMath::Min(ymin,ymnF/sumF*rnorm6);

  TH1F *waku = new TH1F("","",1,xmin,xmax);
  if (imode%2 == 0) { waku -> SetMaximum(1.2*ymax); }
  if (imode%2 == 1) { waku -> SetMaximum(2.0*ymax); }
  if (imode%2 == 0) { waku -> SetMinimum(0.0); }
  if (imode%2 == 1) { waku -> SetMinimum(ymin); }
  waku -> GetXaxis() -> SetTitle(xtit.c_str());
  waku -> GetYaxis() -> SetTitle(ytit.c_str());
  waku -> DrawCopy(); 

  histA -> SetNormFactor(rnorm1);
  histB -> SetNormFactor(rnorm2);
  histC -> SetNormFactor(rnorm3);
  histD -> SetNormFactor(rnorm4);
  histE -> SetNormFactor(rnorm5);
  histF -> SetNormFactor(rnorm6);

  histA -> SetLineColor(icl1);
  histB -> SetLineColor(icl2);
  histC -> SetLineColor(icl3);
  histD -> SetLineColor(icl4);
  histE -> SetLineColor(icl5);
  histF -> SetLineColor(icl6);

  histA -> SetLineStyle(ist1);
  histB -> SetLineStyle(ist2);
  histC -> SetLineStyle(ist3);
  histD -> SetLineStyle(ist4);
  histE -> SetLineStyle(ist5);
  histF -> SetLineStyle(ist6);

  histA -> SetLineWidth(2);
  histB -> SetLineWidth(2);
  histC -> SetLineWidth(2);
  histD -> SetLineWidth(2);
  histE -> SetLineWidth(2);
  histF -> SetLineWidth(2);

  histF -> Draw("HIST same");
  histE -> Draw("HIST same");
  histD -> Draw("HIST same");
  histC -> Draw("HIST same");
  histB -> Draw("HIST same");
  histA -> Draw("HIST same");
  waku -> Delete();
}

void alabel(int itype) {
  if (itype == 1) {
    TPad *tp = new TPad("","",0.5,0.85,0.95,0.95,0,0,0);
    tp -> Draw("A");
    tp -> cd();
    TLatex *text = new TLatex(0.3,0.7,"pretag");
    text -> SetTextSize(0.5);
    text -> DrawLatex(0.25,0.5,"(A)");
    TLine *line = new TLine(0.05,0.75,0.25,0.75);
    line->SetLineWidth(2);
    line->SetLineColor(1); line->SetLineStyle(1);
    line->DrawLine(0.05,0.55,0.2,0.55);
  }
  if (itype == 2) {
    TPad *tp = new TPad("","",0.5,0.8,0.95,0.95,0,0,0);
    tp -> Draw("A");
    tp -> cd();
    TLatex *text = new TLatex(0.3,0.7,"pretag");
    text -> SetTextSize(0.3);
    text -> DrawLatex(0.25,0.75,"(A)");
    text -> DrawLatex(0.25,0.25,"(B)");
    TLine *line = new TLine(0.05,0.75,0.25,0.75);
    line->SetLineWidth(2);
    line->SetLineColor(1); line->SetLineStyle(1); 
    line->DrawLine(0.05,0.8,0.2,0.8);
    line->SetLineColor(4); line->SetLineStyle(1); 
    line->DrawLine(0.05,0.3,0.2,0.3);
  }
  if (itype == 3) {
    TPad *tp = new TPad("","",0.5,0.78,0.95,0.95,0,0,0);
    tp -> Draw("A");
    tp -> cd();
    TLatex *text = new TLatex(0.3,0.7,"pretag");
    text -> SetTextSize(0.26);
    text -> DrawLatex(0.25,0.8,"(A)");
    text -> DrawLatex(0.25,0.5,"(B)");
    text -> DrawLatex(0.25,0.2,"(C)");

    TLine *line = new TLine(0.05,0.75,0.25,0.75);
    line->SetLineWidth(2);
    line->SetLineColor(1); line->SetLineStyle(1); 
    line->DrawLine(0.05,0.85,0.2,0.85);
    line->SetLineColor(2); line->SetLineStyle(1); 
    line->DrawLine(0.05,0.55,0.2,0.55);
    line->SetLineColor(4); line->SetLineStyle(1); 
    line->DrawLine(0.05,0.25,0.2,0.25);
  }
  if (itype == 4) {
    TPad *tp = new TPad("","",0.5,0.75,0.95,0.95,0,0,0);
    tp -> Draw("A");
    tp -> cd();
    TLatex *text = new TLatex(0.3,0.7,"pretag");
    text -> SetTextSize(0.26);
    text -> DrawLatex(0.25,0.88,"VBF H#rightarrow #tau^{+}#tau^{-}");
    text -> DrawLatex(0.25,0.62,"Z#rightarrow e^{+}e^{-}+N jets");
    text -> DrawLatex(0.25,0.38,"Z#rightarrow #mu^{+}#mu^{-}+N jets");
    text -> DrawLatex(0.25,0.12,"W#rightarrow e#nu+N jets");
 
    TLine *line = new TLine(0.05,0.75,0.25,0.75);
    line->SetLineWidth(2);
    line->SetLineColor(1); line->SetLineStyle(1);
    line->DrawLine(0.05,0.94,0.2,0.94);
    line->SetLineColor(2); line->SetLineStyle(1);
    line->DrawLine(0.05,0.68,0.2,0.68);
    line->SetLineColor(4); line->SetLineStyle(1);
    line->DrawLine(0.05,0.44,0.2,0.44);
    line->SetLineColor(6); line->SetLineStyle(2);
    line->DrawLine(0.05,0.18,0.2,0.18);
  }
  if (itype == 5) {
    TPad *tp = new TPad("","",0.5,0.72,0.95,0.95,0,0,0);
    tp -> Draw("A");
    tp -> cd();
    TLatex *text = new TLatex(0.3,0.7,"pretag");
    text -> SetTextSize(0.22);
    text -> DrawLatex(0.25,0.90,"VBF H#rightarrow #tau^{+}#tau^{-}");
    text -> DrawLatex(0.25,0.70,"t#bar{t}");
    text -> DrawLatex(0.25,0.50,"Z#rightarrow e^{+}e^{-}+N jets");
    text -> DrawLatex(0.25,0.30,"Z#rightarrow #mu^{+}#mu^{-}+N jets");
    text -> DrawLatex(0.25,0.10,"W#rightarrow e#nu+N jets");
 
    TLine *line = new TLine(0.05,0.75,0.25,0.75);
    line->SetLineWidth(2);
    line->SetLineColor(1); line->SetLineStyle(1);
    line->DrawLine(0.05,0.95,0.2,0.95);
    line->SetLineColor(2); line->SetLineStyle(1);
    line->DrawLine(0.05,0.75,0.2,0.75);
    line->SetLineColor(4); line->SetLineStyle(1);
    line->DrawLine(0.05,0.55,0.2,0.55);
    line->SetLineColor(6); line->SetLineStyle(2);
    line->DrawLine(0.05,0.35,0.2,0.35);
    line->SetLineColor(7); line->SetLineStyle(2);
    line->DrawLine(0.05,0.15,0.2,0.15);
  }
}

void GraphPlot(const int N,double GX[N],double GY[N],double GXe[N],
	       double GYe[N],int mSyl,int mCol,double mSiz) {
  TGraphErrors *TG = new TGraphErrors(N,GX,GY,GXe,GYe);
  TG -> SetMarkerStyle(mSyl);
  TG -> SetMarkerSize(mSiz);
  TG -> SetMarkerColor(mCol);
  TG -> SetLineColor(mCol);
  TG -> Draw("PE1");
  TG -> Clear();
}

void makePlot_origin(){

  char xtit[50];
  char ytit[50];
  const int nfil = 2;
  std::string filename[nfil];

  std::string plottit = "plot.ps";
  filename[0] = "hist_data18_13TeV.00359441.physics_Main.merge.EVENT_NT.f964_m1831_r11266_p3670._0005.pool.root";
  //filename[0] = "hist_trig1_misal1_mc12.005334.HerwigVBFH120tautaulh.recon.BaseNT.v12000601_signal.root";
  //filename[1] = "hist_AlpgenJimmyZtautauNpVBFCut.root";
  filename[1] = "hist_data16_13TeV.00306448.physics_Main.recon.EVENT_NT.r9264_r11492._0008.pool.root";
  //filename[2] = "hist_trig1_misal1_mc12.005200.T1_McAtNlo_Jimmy.recon.BaseNT.v12000604_tid008037._00241.pool.root";
  //filename[3] = "hist_AlpgenJimmyZeeNpVBFCut.root";
  //filename[4] = "hist_AlpgenJimmyZmumuNpVBFCut.root";
  //filename[5] = "hist_AlpgenJimmyWlnuNpVBFCut.root";
 /*
  TH1F* hist_muValue[nfil];
  TH1F* hist_lumiBlk[nfil];
  TH1F* hist_numVtx[nfil];
  TH1F* hist_numJet[nfil];
  TH1F* hist_numJetjvt[nfil];
  TH1F* hist_muon_pt1[nfil];
  TH1F* hist_muon_eta1[nfil];
  TH1F* hist_muon_pt2[nfil];
  TH1F* hist_muon_eta2[nfil];
  TH1F* hist_dimuon_mass[nfil];
  TH1F* hist_dimuon_pt[nfil];
  TH1F* hist_dimuon_eta[nfil];
  TH1F* hist_trk_num[nfil];
  TH1F* hist_trk_pt[nfil];
  TH1F* hist_trk_eta[nfil];
  TH1F* hist_trk_phi[nfil];
  TH1F* hist_trk_qoverp[nfil];
  TH1F* hist_trk_d0[nfil];
  TH1F* hist_trk_z0[nfil];
  TH1F* hist_trk_deltaZ[nfil];
  TH2F* hist_trk_etaVSphi[nfil];
  TH1F* hist_truth_pt[nfil];
  TH1F* hist_truth_eta[nfil];
  TH1F* hist_truth_phi[nfil];
  TH1F* hist_truth_d0[nfil];
  TH1F* hist_truth_z0[nfil];
  TH1F* hist_trk_dpt[nfil];
  TH1F* hist_trk_dphi[nfil];
  TH1F* hist_trk_deta[nfil];
  TH1F* hist_trk_dd0[nfil];
  TH1F* hist_trk_dz0[nfil];
  TH1F* hist_trk_nPixHits[nfil];
  TH1F* hist_trk_nSCTHits[nfil];
  TH1F* hist_trk_nSiHits[nfil];
  TH1F* hist_trk_nGangedPix[nfil];
  TH1F* hist_trk_nPixLay[nfil];
  TH1F* hist_trk_nPixSharedHits[nfil];
  TH1F* hist_trk_nPixSplitHits[nfil];
  TH1F* hist_trk_nPixOutliers[nfil];
  TH1F* hist_trk_nPixHoles[nfil];
  TH1F* hist_trk_nPixelDeadSensors[nfil];
  TH1F* hist_trk_nSCTSharedHits[nfil];
  TH1F* hist_trk_nSCTOutliers[nfil];
  TH1F* hist_trk_nSCTHoles[nfil];
  TH1F* hist_trk_nSCTDeadSensors[nfil];
  TH1F* hist_trk_nTRTHits[nfil];
  TH1F* hist_trk_nTRTOutliers[nfil];
  TH1F* hist_trk_nTRTHoles[nfil];
  TH1F* hist_trk_nTRTHTHits[nfil];
  TH1F* hist_trk_chiSqPerDof[nfil];
  TH1F* hist_trk_nOutliers[nfil];
  TH1F* hist_trk_stdDevChi2OS[nfil];
  TH1F* hist_trk_truthMatchProb[nfil];
  TH2F* hist_trk_weight[nfil];
  TH1F* hist_IBL_IsEdge[nfil];
  TH1F* hist_IBL_IsOverflow[nfil];
  TH1F* hist_IBL_IsSplit[nfil];
  TH1F* hist_IBL_L1A[nfil];
  TH1F* hist_IBL_ToT[nfil];
  TH1F* hist_IBL_Charge[nfil];*/
  TH1F* hist_IBL_dEdx[nfil];
  TH2F* hist_IBL_dEdxVsP[nfil];
  TH1F* hist_IBL_HitSize[nfil];
  TH1F* hist_IBL_HitSizePhi[nfil];
  /*TH1F* hist_IBL_HitSizeZ[nfil];
  TH1F* hist_IBL_HitSizeRatio[nfil];
  TH1F* hist_IBL_unbiasedResidualX[nfil];
  TH1F* hist_IBL_unbiasedResidualY[nfil];
  TH1F* hist_IBL_Isolation10x2[nfil];
  TH1F* hist_IBL_Isolation20x4[nfil];
  TH1F* hist_IBL_numTotalClustersPerModule[nfil];
  TH1F* hist_IBL_numTotalPixelsPerModule[nfil];
  TProfile* profile_IBL_LorentzAngle[nfil];
  TH2F* hist_IBL_Map[nfil];
  TH2F* hist_IBL_MapHit[nfil];
  TH1F* hist_IBL_MapEta[nfil];
  TH1F* hist_IBL_MapHitEta[nfil];
  TH1F* hist_IBL_dEdx_true[nfil];
  TH2F* hist_IBL_dEdxVsP_true[nfil];
  TH1F* hist_IBL_ChargeFraction[nfil];
  TH2F* hist_IBL_ClusterShape[nfil];
  TH2F* hist_IBL_ModuleMap[nfil];
  TH1F* hist_BLY_IsEdge[nfil];
  TH1F* hist_BLY_IsOverflow[nfil];
  TH1F* hist_BLY_IsSplit[nfil];
  TH1F* hist_BLY_L1A[nfil];
  TH1F* hist_BLY_ToT[nfil];
  TH1F* hist_BLY_Charge[nfil];
  TH1F* hist_BLY_dEdx[nfil];
  TH2F* hist_BLY_dEdxVsP[nfil];
  TH1F* hist_BLY_HitSize[nfil];
  TH1F* hist_BLY_HitSizePhi[nfil];
  TH1F* hist_BLY_HitSizeZ[nfil];
  TH1F* hist_BLY_unbiasedResidualX[nfil];
  TH1F* hist_BLY_unbiasedResidualY[nfil];
  TH1F* hist_BLY_Isolation10x2[nfil];
  TH1F* hist_BLY_Isolation20x4[nfil];
  TH1F* hist_BLY_numTotalClustersPerModule[nfil];
  TH1F* hist_BLY_numTotalPixelsPerModule[nfil];
  TProfile* profile_BLY_LorentzAngle[nfil];
  TH2F* hist_BLY_Map[nfil];
  TH2F* hist_BLY_MapHit[nfil];
  TH1F* hist_BLY_MapEta[nfil];
  TH1F* hist_BLY_MapHitEta[nfil];
  TH1F* hist_BLY_dEdx_true[nfil];
  TH2F* hist_BLY_dEdxVsP_true[nfil];
  TH1F* hist_BLY_ChargeFraction[nfil];
  TH2F* hist_BLY_ClusterShape[nfil];
  TH1F* hist_LY1_IsEdge[nfil];
  TH1F* hist_LY1_IsOverflow[nfil];
  TH1F* hist_LY1_IsSplit[nfil];
  TH1F* hist_LY1_L1A[nfil];
  TH1F* hist_LY1_ToT[nfil];
  TH1F* hist_LY1_Charge[nfil];
  TH1F* hist_LY1_dEdx[nfil];
  TH2F* hist_LY1_dEdxVsP[nfil];
  TH1F* hist_LY1_HitSize[nfil];
  TH1F* hist_LY1_HitSizePhi[nfil];
  TH1F* hist_LY1_HitSizeZ[nfil];
  TH1F* hist_LY1_unbiasedResidualX[nfil];
  TH1F* hist_LY1_unbiasedResidualY[nfil];
  TH1F* hist_LY1_Isolation10x2[nfil];
  TH1F* hist_LY1_Isolation20x4[nfil];
  TH1F* hist_LY1_numTotalClustersPerModule[nfil];
  TH1F* hist_LY1_numTotalPixelsPerModule[nfil];
  TProfile* profile_LY1_LorentzAngle[nfil];
  TH2F* hist_LY1_Map[nfil];
  TH2F* hist_LY1_MapHit[nfil];
  TH1F* hist_LY1_MapEta[nfil];
  TH1F* hist_LY1_MapHitEta[nfil];
  TH1F* hist_LY1_dEdx_true[nfil];
  TH2F* hist_LY1_dEdxVsP_true[nfil];
  TH1F* hist_LY2_IsEdge[nfil];
  TH1F* hist_LY2_IsOverflow[nfil];
  TH1F* hist_LY2_IsSplit[nfil];
  TH1F* hist_LY2_L1A[nfil];
  TH1F* hist_LY2_ToT[nfil];
  TH1F* hist_LY2_Charge[nfil];
  TH1F* hist_LY2_dEdx[nfil];
  TH2F* hist_LY2_dEdxVsP[nfil];
  TH1F* hist_LY2_HitSize[nfil];
  TH1F* hist_LY2_HitSizePhi[nfil];
  TH1F* hist_LY2_HitSizeZ[nfil];
  TH1F* hist_LY2_unbiasedResidualX[nfil];
  TH1F* hist_LY2_unbiasedResidualY[nfil];
  TH1F* hist_LY2_Isolation10x2[nfil];
  TH1F* hist_LY2_Isolation20x4[nfil];
  TH1F* hist_LY2_numTotalClustersPerModule[nfil];
  TH1F* hist_LY2_numTotalPixelsPerModule[nfil];
  TProfile* profile_LY2_LorentzAngle[nfil];
  TH2F* hist_LY2_Map[nfil];
  TH2F* hist_LY2_MapHit[nfil];
  TH1F* hist_LY2_MapEta[nfil];
  TH1F* hist_LY2_MapHitEta[nfil];
  TH1F* hist_LY2_dEdx_true[nfil];
  TH2F* hist_LY2_dEdxVsP_true[nfil];
  TH1F* hist_END_IsEdge[nfil];
  TH1F* hist_END_IsOverflow[nfil];
  TH1F* hist_END_IsSplit[nfil];
  TH1F* hist_END_L1A[nfil];
  TH1F* hist_END_ToT[nfil];
  TH1F* hist_END_Charge[nfil];
  TH1F* hist_END_dEdx[nfil];
  TH2F* hist_END_dEdxVsP[nfil];
  TH1F* hist_END_HitSize[nfil];
  TH1F* hist_END_HitSizePhi[nfil];
  TH1F* hist_END_HitSizeZ[nfil];
  TH1F* hist_END_unbiasedResidualX[nfil];
  TH1F* hist_END_unbiasedResidualY[nfil];
  TH1F* hist_END_Isolation10x2[nfil];
  TH1F* hist_END_Isolation20x4[nfil];
  TH1F* hist_END_numTotalClustersPerModule[nfil];
  TH1F* hist_END_numTotalPixelsPerModule[nfil];
  TProfile* profile_END_LorentzAngle[nfil];
  TH2F* hist_ED1_Map[nfil];
  TH2F* hist_ED1_MapHit[nfil];
  TH1F* hist_ED1_MapEta[nfil];
  TH1F* hist_ED1_MapHitEta[nfil];
  TH2F* hist_ED2_3Map[nfil];
  TH2F* hist_ED2_3MapHit[nfil];
  TH1F* hist_ED2_3MapEta[nfil];
  TH1F* hist_ED2_3MapHitEta[nfil];
  TH1F* hist_END_dEdx_true[nfil];
  TH2F* hist_END_dEdxVsP_true[nfil];
  TH1F* hist_ALL_dEdx[nfil];
  TH2F* hist_ALL_dEdxVsP[nfil];
  TH1F* hist_ALL_dEdx_true[nfil];
  TH2F* hist_ALL_dEdxVsP_true[nfil];*/

  TFile *ffl[nfil];
  for (int i = 0; i < nfil; i++) {
    ffl[i] = new TFile(filename[i].c_str(),"read");

    /*hist_muValue[i] = (TH1F*) ffl[i]->Get("hist_muValue");
    hist_lumiBlk[i] = (TH1F*) ffl[i]->Get("hist_lumiBlk");
    hist_numVtx[i] = (TH1F*) ffl[i]->Get("hist_numVtx");
    hist_numJet[i] = (TH1F*) ffl[i]->Get("hist_numJet");
    hist_numJetjvt[i] = (TH1F*) ffl[i]->Get("hist_numJetjvt");
    hist_muon_pt1[i] = (TH1F*) ffl[i]->Get("hist_muon_pt1");
    hist_muon_eta1[i] = (TH1F*) ffl[i]->Get("hist_muon_eta1");
    hist_muon_pt2[i] = (TH1F*) ffl[i]->Get("hist_muon_pt2");
    hist_muon_eta2[i] = (TH1F*) ffl[i]->Get("hist_muon_eta2");
    hist_dimuon_mass[i] = (TH1F*) ffl[i]->Get("hist_dimuon_mass");
    hist_dimuon_pt[i] = (TH1F*) ffl[i]->Get("hist_dimuon_pt");
    hist_dimuon_eta[i] = (TH1F*) ffl[i]->Get("hist_dimuon_eta");
    hist_trk_num[i] = (TH1F*) ffl[i]->Get("hist_trk_num");
    hist_trk_pt[i] = (TH1F*) ffl[i]->Get("hist_trk_pt");
    hist_trk_eta[i] = (TH1F*) ffl[i]->Get("hist_trk_eta");
    hist_trk_phi[i] = (TH1F*) ffl[i]->Get("hist_trk_phi");
    hist_trk_qoverp[i] = (TH1F*) ffl[i]->Get("hist_trk_qoverp");
    hist_trk_d0[i] = (TH1F*) ffl[i]->Get("hist_trk_d0");
    hist_trk_z0[i] = (TH1F*) ffl[i]->Get("hist_trk_z0");
    hist_trk_deltaZ[i] = (TH1F*) ffl[i]->Get("hist_trk_deltaZ");
    hist_trk_etaVSphi[i] = (TH2F*) ffl[i]->Get("hist_trk_etaVSphi");
    hist_truth_pt[i] = (TH1F*) ffl[i]->Get("hist_truth_pt");
    hist_truth_eta[i] = (TH1F*) ffl[i]->Get("hist_truth_eta");
    hist_truth_phi[i] = (TH1F*) ffl[i]->Get("hist_truth_phi");
    hist_truth_d0[i] = (TH1F*) ffl[i]->Get("hist_truth_d0");
    hist_truth_z0[i] = (TH1F*) ffl[i]->Get("hist_truth_z0");
    hist_trk_dpt[i] = (TH1F*) ffl[i]->Get("hist_trk_dpt");
    hist_trk_dphi[i] = (TH1F*) ffl[i]->Get("hist_trk_dphi");
    hist_trk_deta[i] = (TH1F*) ffl[i]->Get("hist_trk_deta");
    hist_trk_dd0[i] = (TH1F*) ffl[i]->Get("hist_trk_dd0");
    hist_trk_dz0[i] = (TH1F*) ffl[i]->Get("hist_trk_dz0");
    hist_trk_nPixHits[i] = (TH1F*) ffl[i]->Get("hist_trk_nPixHits");
    hist_trk_nSCTHits[i] = (TH1F*) ffl[i]->Get("hist_trk_nSCTHits");
    hist_trk_nSiHits[i] = (TH1F*) ffl[i]->Get("hist_trk_nSiHits");
    hist_trk_nGangedPix[i] = (TH1F*) ffl[i]->Get("hist_trk_nGangedPix");
    hist_trk_nPixLay[i] = (TH1F*) ffl[i]->Get("hist_trk_nPixLay");
    hist_trk_nPixSharedHits[i] = (TH1F*) ffl[i]->Get("hist_trk_nPixSharedHits");
    hist_trk_nPixSplitHits[i] = (TH1F*) ffl[i]->Get("hist_trk_nPixSplitHits");
    hist_trk_nPixOutliers[i] = (TH1F*) ffl[i]->Get("hist_trk_nPixOutliers");
    hist_trk_nPixHoles[i] = (TH1F*) ffl[i]->Get("hist_trk_nPixHoles");
    hist_trk_nPixelDeadSensors[i] = (TH1F*) ffl[i]->Get("hist_trk_nPixelDeadSensors");
    hist_trk_nSCTSharedHits[i] = (TH1F*) ffl[i]->Get("hist_trk_nSCTSharedHits");
    hist_trk_nSCTOutliers[i] = (TH1F*) ffl[i]->Get("hist_trk_nSCTOutliers");
    hist_trk_nSCTHoles[i] = (TH1F*) ffl[i]->Get("hist_trk_nSCTHoles");
    hist_trk_nSCTDeadSensors[i] = (TH1F*) ffl[i]->Get("hist_trk_nSCTDeadSensors");
    hist_trk_nTRTHits[i] = (TH1F*) ffl[i]->Get("hist_trk_nTRTHits");
    hist_trk_nTRTOutliers[i] = (TH1F*) ffl[i]->Get("hist_trk_nTRTOutliers");
    hist_trk_nTRTHoles[i] = (TH1F*) ffl[i]->Get("hist_trk_nTRTHoles");
    hist_trk_nTRTHTHits[i] = (TH1F*) ffl[i]->Get("hist_trk_nTRTHTHits");
    hist_trk_chiSqPerDof[i] = (TH1F*) ffl[i]->Get("hist_trk_chiSqPerDof");
    hist_trk_nOutliers[i] = (TH1F*) ffl[i]->Get("hist_trk_nOutliers");
    hist_trk_stdDevChi2OS[i] = (TH1F*) ffl[i]->Get("hist_trk_stdDevChi2OS");
    hist_trk_truthMatchProb[i] = (TH1F*) ffl[i]->Get("hist_trk_truthMatchProb");
    hist_trk_weight[i] = (TH2F*) ffl[i]->Get("hist_trk_weight");
    hist_IBL_IsEdge[i] = (TH1F*) ffl[i]->Get("hist_IBL_IsEdge");
    hist_IBL_IsOverflow[i] = (TH1F*) ffl[i]->Get("hist_IBL_IsOverflow");
    hist_IBL_IsSplit[i] = (TH1F*) ffl[i]->Get("hist_IBL_IsSplit");
    hist_IBL_L1A[i] = (TH1F*) ffl[i]->Get("hist_IBL_L1A");
    hist_IBL_ToT[i] = (TH1F*) ffl[i]->Get("hist_IBL_ToT");
    hist_IBL_Charge[i] = (TH1F*) ffl[i]->Get("hist_IBL_Charge");*/
    hist_IBL_dEdx[i] = (TH1F*) ffl[i]->Get("hist_IBL_dEdx");
    hist_IBL_dEdxVsP[i] = (TH2F*) ffl[i]->Get("hist_IBL_dEdxVsP");
    hist_IBL_HitSize[i] = (TH1F*) ffl[i]->Get("hist_IBL_HitSize");
    hist_IBL_HitSizePhi[i] = (TH1F*) ffl[i]->Get("hist_IBL_HitSizePhi");
    /*hist_IBL_HitSizeZ[i] = (TH1F*) ffl[i]->Get("hist_IBL_HitSizeZ");
    hist_IBL_HitSizeRatio[i] = (TH1F*) ffl[i]->Get("hist_IBL_HitSizeRatio");
    hist_IBL_unbiasedResidualX[i] = (TH1F*) ffl[i]->Get("hist_IBL_unbiasedResidualX");
    hist_IBL_unbiasedResidualY[i] = (TH1F*) ffl[i]->Get("hist_IBL_unbiasedResidualY");
    hist_IBL_Isolation10x2[i] = (TH1F*) ffl[i]->Get("hist_IBL_Isolation10x2");
    hist_IBL_Isolation20x4[i] = (TH1F*) ffl[i]->Get("hist_IBL_Isolation20x4");
    hist_IBL_numTotalClustersPerModule[i] = (TH1F*) ffl[i]->Get("hist_IBL_numTotalClustersPerModule");
    hist_IBL_numTotalPixelsPerModule[i] = (TH1F*) ffl[i]->Get("hist_IBL_numTotalPixelsPerModule");
    profile_IBL_LorentzAngle[i] = (TProfile*) ffl[i]->Get("profile_IBL_LorentzAngle");
    hist_IBL_Map[i] = (TH2F*) ffl[i]->Get("hist_IBL_Map");
    hist_IBL_MapHit[i] = (TH2F*) ffl[i]->Get("hist_IBL_MapHit");
    hist_IBL_MapEta[i] = (TH1F*) ffl[i]->Get("hist_IBL_MapEta");
    hist_IBL_MapHitEta[i] = (TH1F*) ffl[i]->Get("hist_IBL_MapHitEta");
    hist_IBL_dEdx_true[i] = (TH1F*) ffl[i]->Get("hist_IBL_dEdx_true");
    hist_IBL_dEdxVsP_true[i] = (TH2F*) ffl[i]->Get("hist_IBL_dEdxVsP_true");
    hist_IBL_ChargeFraction[i] = (TH1F*) ffl[i]->Get("hist_IBL_ChargeFraction");
    hist_IBL_ClusterShape[i] = (TH2F*) ffl[i]->Get("hist_IBL_ClusterShape");
    hist_IBL_ModuleMap[i] = (TH2F*) ffl[i]->Get("hist_IBL_ModuleMap");
    hist_BLY_IsEdge[i] = (TH1F*) ffl[i]->Get("hist_BLY_IsEdge");
    hist_BLY_IsOverflow[i] = (TH1F*) ffl[i]->Get("hist_BLY_IsOverflow");
    hist_BLY_IsSplit[i] = (TH1F*) ffl[i]->Get("hist_BLY_IsSplit");
    hist_BLY_L1A[i] = (TH1F*) ffl[i]->Get("hist_BLY_L1A");
    hist_BLY_ToT[i] = (TH1F*) ffl[i]->Get("hist_BLY_ToT");
    hist_BLY_Charge[i] = (TH1F*) ffl[i]->Get("hist_BLY_Charge");
    hist_BLY_dEdx[i] = (TH1F*) ffl[i]->Get("hist_BLY_dEdx");
    hist_BLY_dEdxVsP[i] = (TH2F*) ffl[i]->Get("hist_BLY_dEdxVsP");
    hist_BLY_HitSize[i] = (TH1F*) ffl[i]->Get("hist_BLY_HitSize");
    hist_BLY_HitSizePhi[i] = (TH1F*) ffl[i]->Get("hist_BLY_HitSizePhi");
    hist_BLY_HitSizeZ[i] = (TH1F*) ffl[i]->Get("hist_BLY_HitSizeZ");
    hist_BLY_unbiasedResidualX[i] = (TH1F*) ffl[i]->Get("hist_BLY_unbiasedResidualX");
    hist_BLY_unbiasedResidualY[i] = (TH1F*) ffl[i]->Get("hist_BLY_unbiasedResidualY");
    hist_BLY_Isolation10x2[i] = (TH1F*) ffl[i]->Get("hist_BLY_Isolation10x2");
    hist_BLY_Isolation20x4[i] = (TH1F*) ffl[i]->Get("hist_BLY_Isolation20x4");
    hist_BLY_numTotalClustersPerModule[i] = (TH1F*) ffl[i]->Get("hist_BLY_numTotalClustersPerModule");
    hist_BLY_numTotalPixelsPerModule[i] = (TH1F*) ffl[i]->Get("hist_BLY_numTotalPixelsPerModule");
    profile_BLY_LorentzAngle[i] = (TProfile*) ffl[i]->Get("profile_BLY_LorentzAngle");
    hist_BLY_Map[i] = (TH2F*) ffl[i]->Get("hist_BLY_Map");
    hist_BLY_MapHit[i] = (TH2F*) ffl[i]->Get("hist_BLY_MapHit");
    hist_BLY_MapEta[i] = (TH1F*) ffl[i]->Get("hist_BLY_MapEta");
    hist_BLY_MapHitEta[i] = (TH1F*) ffl[i]->Get("hist_BLY_MapHitEta");
    hist_BLY_dEdx_true[i] = (TH1F*) ffl[i]->Get("hist_BLY_dEdx_true");
    hist_BLY_dEdxVsP_true[i] = (TH2F*) ffl[i]->Get("hist_BLY_dEdxVsP_true");
    hist_BLY_ChargeFraction[i] = (TH1F*) ffl[i]->Get("hist_BLY_ChargeFraction");
    hist_BLY_ClusterShape[i] = (TH2F*) ffl[i]->Get("hist_BLY_ClusterShape");
    hist_LY1_IsEdge[i] = (TH1F*) ffl[i]->Get("hist_LY1_IsEdge");
    hist_LY1_IsOverflow[i] = (TH1F*) ffl[i]->Get("hist_LY1_IsOverflow");
    hist_LY1_IsSplit[i] = (TH1F*) ffl[i]->Get("hist_LY1_IsSplit");
    hist_LY1_L1A[i] = (TH1F*) ffl[i]->Get("hist_LY1_L1A");
    hist_LY1_ToT[i] = (TH1F*) ffl[i]->Get("hist_LY1_ToT");
    hist_LY1_Charge[i] = (TH1F*) ffl[i]->Get("hist_LY1_Charge");
    hist_LY1_dEdx[i] = (TH1F*) ffl[i]->Get("hist_LY1_dEdx");
    hist_LY1_dEdxVsP[i] = (TH2F*) ffl[i]->Get("hist_LY1_dEdxVsP");
    hist_LY1_HitSize[i] = (TH1F*) ffl[i]->Get("hist_LY1_HitSize");
    hist_LY1_HitSizePhi[i] = (TH1F*) ffl[i]->Get("hist_LY1_HitSizePhi");
    hist_LY1_HitSizeZ[i] = (TH1F*) ffl[i]->Get("hist_LY1_HitSizeZ");
    hist_LY1_unbiasedResidualX[i] = (TH1F*) ffl[i]->Get("hist_LY1_unbiasedResidualX");
    hist_LY1_unbiasedResidualY[i] = (TH1F*) ffl[i]->Get("hist_LY1_unbiasedResidualY");
    hist_LY1_Isolation10x2[i] = (TH1F*) ffl[i]->Get("hist_LY1_Isolation10x2");
    hist_LY1_Isolation20x4[i] = (TH1F*) ffl[i]->Get("hist_LY1_Isolation20x4");
    hist_LY1_numTotalClustersPerModule[i] = (TH1F*) ffl[i]->Get("hist_LY1_numTotalClustersPerModule");
    hist_LY1_numTotalPixelsPerModule[i] = (TH1F*) ffl[i]->Get("hist_LY1_numTotalPixelsPerModule");
    profile_LY1_LorentzAngle[i] = (TProfile*) ffl[i]->Get("profile_LY1_LorentzAngle");
    hist_LY1_Map[i] = (TH2F*) ffl[i]->Get("hist_LY1_Map");
    hist_LY1_MapHit[i] = (TH2F*) ffl[i]->Get("hist_LY1_MapHit");
    hist_LY1_MapEta[i] = (TH1F*) ffl[i]->Get("hist_LY1_MapEta");
    hist_LY1_MapHitEta[i] = (TH1F*) ffl[i]->Get("hist_LY1_MapHitEta");
    hist_LY1_dEdx_true[i] = (TH1F*) ffl[i]->Get("hist_LY1_dEdx_true");
    hist_LY1_dEdxVsP_true[i] = (TH2F*) ffl[i]->Get("hist_LY1_dEdxVsP_true");
    hist_LY2_IsEdge[i] = (TH1F*) ffl[i]->Get("hist_LY2_IsEdge");
    hist_LY2_IsOverflow[i] = (TH1F*) ffl[i]->Get("hist_LY2_IsOverflow");
    hist_LY2_IsSplit[i] = (TH1F*) ffl[i]->Get("hist_LY2_IsSplit");
    hist_LY2_L1A[i] = (TH1F*) ffl[i]->Get("hist_LY2_L1A");
    hist_LY2_ToT[i] = (TH1F*) ffl[i]->Get("hist_LY2_ToT");
    hist_LY2_Charge[i] = (TH1F*) ffl[i]->Get("hist_LY2_Charge");
    hist_LY2_dEdx[i] = (TH1F*) ffl[i]->Get("hist_LY2_dEdx");
    hist_LY2_dEdxVsP[i] = (TH2F*) ffl[i]->Get("hist_LY2_dEdxVsP");
    hist_LY2_HitSize[i] = (TH1F*) ffl[i]->Get("hist_LY2_HitSize");
    hist_LY2_HitSizePhi[i] = (TH1F*) ffl[i]->Get("hist_LY2_HitSizePhi");
    hist_LY2_HitSizeZ[i] = (TH1F*) ffl[i]->Get("hist_LY2_HitSizeZ");
    hist_LY2_unbiasedResidualX[i] = (TH1F*) ffl[i]->Get("hist_LY2_unbiasedResidualX");
    hist_LY2_unbiasedResidualY[i] = (TH1F*) ffl[i]->Get("hist_LY2_unbiasedResidualY");
    hist_LY2_Isolation10x2[i] = (TH1F*) ffl[i]->Get("hist_LY2_Isolation10x2");
    hist_LY2_Isolation20x4[i] = (TH1F*) ffl[i]->Get("hist_LY2_Isolation20x4");
    hist_LY2_numTotalClustersPerModule[i] = (TH1F*) ffl[i]->Get("hist_LY2_numTotalClustersPerModule");
    hist_LY2_numTotalPixelsPerModule[i] = (TH1F*) ffl[i]->Get("hist_LY2_numTotalPixelsPerModule");
    profile_LY2_LorentzAngle[i] = (TProfile*) ffl[i]->Get("profile_LY2_LorentzAngle");
    hist_LY2_Map[i] = (TH2F*) ffl[i]->Get("hist_LY2_Map");
    hist_LY2_MapHit[i] = (TH2F*) ffl[i]->Get("hist_LY2_MapHit");
    hist_LY2_MapEta[i] = (TH1F*) ffl[i]->Get("hist_LY2_MapEta");
    hist_LY2_MapHitEta[i] = (TH1F*) ffl[i]->Get("hist_LY2_MapHitEta");
    hist_LY2_dEdx_true[i] = (TH1F*) ffl[i]->Get("hist_LY2_dEdx_true");
    hist_LY2_dEdxVsP_true[i] = (TH2F*) ffl[i]->Get("hist_LY2_dEdxVsP_true");
    hist_END_IsEdge[i] = (TH1F*) ffl[i]->Get("hist_END_IsEdge");
    hist_END_IsOverflow[i] = (TH1F*) ffl[i]->Get("hist_END_IsOverflow");
    hist_END_IsSplit[i] = (TH1F*) ffl[i]->Get("hist_END_IsSplit");
    hist_END_L1A[i] = (TH1F*) ffl[i]->Get("hist_END_L1A");
    hist_END_ToT[i] = (TH1F*) ffl[i]->Get("hist_END_ToT");
    hist_END_Charge[i] = (TH1F*) ffl[i]->Get("hist_END_Charge");
    hist_END_dEdx[i] = (TH1F*) ffl[i]->Get("hist_END_dEdx");
    hist_END_dEdxVsP[i] = (TH2F*) ffl[i]->Get("hist_END_dEdxVsP");
    hist_END_HitSize[i] = (TH1F*) ffl[i]->Get("hist_END_HitSize");
    hist_END_HitSizePhi[i] = (TH1F*) ffl[i]->Get("hist_END_HitSizePhi");
    hist_END_HitSizeZ[i] = (TH1F*) ffl[i]->Get("hist_END_HitSizeZ");
    hist_END_unbiasedResidualX[i] = (TH1F*) ffl[i]->Get("hist_END_unbiasedResidualX");
    hist_END_unbiasedResidualY[i] = (TH1F*) ffl[i]->Get("hist_END_unbiasedResidualY");
    hist_END_Isolation10x2[i] = (TH1F*) ffl[i]->Get("hist_END_Isolation10x2");
    hist_END_Isolation20x4[i] = (TH1F*) ffl[i]->Get("hist_END_Isolation20x4");
    hist_END_numTotalClustersPerModule[i] = (TH1F*) ffl[i]->Get("hist_END_numTotalClustersPerModule");
    hist_END_numTotalPixelsPerModule[i] = (TH1F*) ffl[i]->Get("hist_END_numTotalPixelsPerModule");
    profile_END_LorentzAngle[i] = (TProfile*) ffl[i]->Get("profile_END_LorentzAngle");
    hist_ED1_Map[i] = (TH2F*) ffl[i]->Get("hist_ED1_Map");
    hist_ED1_MapHit[i] = (TH2F*) ffl[i]->Get("hist_ED1_MapHit");
    hist_ED1_MapEta[i] = (TH1F*) ffl[i]->Get("hist_ED1_MapEta");
    hist_ED1_MapHitEta[i] = (TH1F*) ffl[i]->Get("hist_ED1_MapHitEta");
    hist_ED2_3Map[i] = (TH2F*) ffl[i]->Get("hist_ED2_3Map");
    hist_ED2_3MapHit[i] = (TH2F*) ffl[i]->Get("hist_ED2_3MapHit");
    hist_ED2_3MapEta[i] = (TH1F*) ffl[i]->Get("hist_ED2_3MapEta");
    hist_ED2_3MapHitEta[i] = (TH1F*) ffl[i]->Get("hist_ED2_3MapHitEta");
    hist_END_dEdx_true[i] = (TH1F*) ffl[i]->Get("hist_END_dEdx_true");
    hist_END_dEdxVsP_true[i] = (TH2F*) ffl[i]->Get("hist_END_dEdxVsP_true");
    hist_ALL_dEdx[i] = (TH1F*) ffl[i]->Get("hist_ALL_dEdx");
    hist_ALL_dEdxVsP[i] = (TH2F*) ffl[i]->Get("hist_ALL_dEdxVsP");
    hist_ALL_dEdx_true[i] = (TH1F*) ffl[i]->Get("hist_ALL_dEdx_true");
    hist_ALL_dEdxVsP_true[i] = (TH2F*) ffl[i]->Get("hist_ALL_dEdxVsP_true");
*/
  }

  TPostScript *ps = new TPostScript(plottit.c_str(),-111);

  TCanvas *c1 = new TCanvas("c1","c1",0,0,500,500);
  c1->SetBorderSize(0);
  gStyle -> SetOptStat(0);
  gStyle -> SetPadBorderMode(0);
  gStyle -> SetPadBorderSize(0);
  gStyle -> SetCanvasBorderMode(0);

  gStyle -> SetPadLeftMargin(0.15);
  //*gStyle -> SetPadRightMargin(0.05);
  //*gStyle -> SetPadTopMargin(0.05);
  gStyle -> SetPadBottomMargin(0.15);

  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.5);
  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetLabelOffset(0.003,"X");
  gStyle->SetLabelOffset(0.003,"Y");

  gStyle->SetPalette(1);
  gStyle->SetNdivisions(510);

  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();

/*
//  gStyle -> SetPadLeftMargin(0.17);
  gStyle -> SetPadLeftMargin(0.18);
//  gStyle->SetTitleYOffset(1.5);
//  gStyle->SetTitleYOffset(1.9);
//  gStyle -> SetPadRightMargin(0.114);
*/

  ps -> NewPage();
/*
  //========
  // page.1
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_muValue[0],
             hist_muValue[1],
             0,0,"hist_muValue","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_lumiBlk[0],
             hist_lumiBlk[1],
             0,0,"hist_lumiBlk","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_numVtx[0],
             hist_numVtx[1],
             0,0,"hist_numVtx","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_numJet[0],
             hist_numJet[1],
             0,0,"hist_numJet","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.2
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_numJetjvt[0],
             hist_numJetjvt[1],
             0,0,"hist_numJetjvt","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_muon_pt1[0],
             hist_muon_pt1[1],
             0,0,"hist_muon_pt1","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_muon_eta1[0],
             hist_muon_eta1[1],
             0,0,"hist_muon_eta1","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_muon_pt2[0],
             hist_muon_pt2[1],
             0,0,"hist_muon_pt2","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.3
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_muon_eta2[0],
             hist_muon_eta2[1],
             0,0,"hist_muon_eta2","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_dimuon_mass[0],
             hist_dimuon_mass[1],
             0,0,"hist_dimuon_mass","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_dimuon_pt[0],
             hist_dimuon_pt[1],
             0,0,"hist_dimuon_pt","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_dimuon_eta[0],
             hist_dimuon_eta[1],
             0,0,"hist_dimuon_eta","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.4
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_num[0],
             hist_trk_num[1],
             0,0,"hist_trk_num","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_pt[0],
             hist_trk_pt[1],
             0,0,"hist_trk_pt","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_eta[0],
             hist_trk_eta[1],
             0,0,"hist_trk_eta","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_phi[0],
             hist_trk_phi[1],
             0,0,"hist_trk_phi","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.5
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_qoverp[0],
             hist_trk_qoverp[1],
             0,0,"hist_trk_qoverp","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_d0[0],
             hist_trk_d0[1],
             0,0,"hist_trk_d0","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_z0[0],
             hist_trk_z0[1],
             0,0,"hist_trk_z0","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_deltaZ[0],
             hist_trk_deltaZ[1],
             0,0,"hist_trk_deltaZ","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.6
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  hist_trk_etaVSphi[0] -> Draw("Box");
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_truth_pt[0],
             hist_truth_pt[1],
             0,0,"hist_truth_pt","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_truth_eta[0],
             hist_truth_eta[1],
             0,0,"hist_truth_eta","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_truth_phi[0],
             hist_truth_phi[1],
             0,0,"hist_truth_phi","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.7
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_truth_d0[0],
             hist_truth_d0[1],
             0,0,"hist_truth_d0","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_truth_z0[0],
             hist_truth_z0[1],
             0,0,"hist_truth_z0","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_dpt[0],
             hist_trk_dpt[1],
             0,0,"hist_trk_dpt","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_dphi[0],
             hist_trk_dphi[1],
             0,0,"hist_trk_dphi","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.8
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_deta[0],
             hist_trk_deta[1],
             0,0,"hist_trk_deta","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_dd0[0],
             hist_trk_dd0[1],
             0,0,"hist_trk_dd0","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_dz0[0],
             hist_trk_dz0[1],
             0,0,"hist_trk_dz0","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nPixHits[0],
             hist_trk_nPixHits[1],
             0,0,"hist_trk_nPixHits","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.9
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nSCTHits[0],
             hist_trk_nSCTHits[1],
             0,0,"hist_trk_nSCTHits","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nSiHits[0],
             hist_trk_nSiHits[1],
             0,0,"hist_trk_nSiHits","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nGangedPix[0],
             hist_trk_nGangedPix[1],
             0,0,"hist_trk_nGangedPix","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nPixLay[0],
             hist_trk_nPixLay[1],
             0,0,"hist_trk_nPixLay","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.10
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nPixSharedHits[0],
             hist_trk_nPixSharedHits[1],
             0,0,"hist_trk_nPixSharedHits","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nPixSplitHits[0],
             hist_trk_nPixSplitHits[1],
             0,0,"hist_trk_nPixSplitHits","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nPixOutliers[0],
             hist_trk_nPixOutliers[1],
             0,0,"hist_trk_nPixOutliers","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nPixHoles[0],
             hist_trk_nPixHoles[1],
             0,0,"hist_trk_nPixHoles","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.11
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nPixelDeadSensors[0],
             hist_trk_nPixelDeadSensors[1],
             0,0,"hist_trk_nPixelDeadSensors","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nSCTSharedHits[0],
             hist_trk_nSCTSharedHits[1],
             0,0,"hist_trk_nSCTSharedHits","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nSCTOutliers[0],
             hist_trk_nSCTOutliers[1],
             0,0,"hist_trk_nSCTOutliers","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nSCTHoles[0],
             hist_trk_nSCTHoles[1],
             0,0,"hist_trk_nSCTHoles","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.12
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nSCTDeadSensors[0],
             hist_trk_nSCTDeadSensors[1],
             0,0,"hist_trk_nSCTDeadSensors","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nTRTHits[0],
             hist_trk_nTRTHits[1],
             0,0,"hist_trk_nTRTHits","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nTRTOutliers[0],
             hist_trk_nTRTOutliers[1],
             0,0,"hist_trk_nTRTOutliers","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nTRTHoles[0],
             hist_trk_nTRTHoles[1],
             0,0,"hist_trk_nTRTHoles","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.13
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nTRTHTHits[0],
             hist_trk_nTRTHTHits[1],
             0,0,"hist_trk_nTRTHTHits","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_chiSqPerDof[0],
             hist_trk_chiSqPerDof[1],
             0,0,"hist_trk_chiSqPerDof","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_nOutliers[0],
             hist_trk_nOutliers[1],
             0,0,"hist_trk_nOutliers","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_stdDevChi2OS[0],
             hist_trk_stdDevChi2OS[1],
             0,0,"hist_trk_stdDevChi2OS","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.14
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_trk_truthMatchProb[0],
             hist_trk_truthMatchProb[1],
             0,0,"hist_trk_truthMatchProb","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  hist_trk_weight[0] -> Draw("Box");
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_IsEdge[0],
             hist_IBL_IsEdge[1],
             0,0,"hist_IBL_IsEdge","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_IsOverflow[0],
             hist_IBL_IsOverflow[1],
             0,0,"hist_IBL_IsOverflow","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.15
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_IsSplit[0],
             hist_IBL_IsSplit[1],
             0,0,"hist_IBL_IsSplit","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_L1A[0],
             hist_IBL_L1A[1],
             0,0,"hist_IBL_L1A","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_ToT[0],
             hist_IBL_ToT[1],
             0,0,"hist_IBL_ToT","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_Charge[0],
             hist_IBL_Charge[1],
             0,0,"hist_IBL_Charge","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  */
  //========
  // page.16
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_dEdx[0],
             hist_IBL_dEdx[1],
             0,0,0.5,"hist_IBL_dEdx","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  hist_IBL_dEdxVsP[0]->GetXaxis()->SetTitle("p"); // X軸ラベルを設定
  hist_IBL_dEdxVsP[0]->GetYaxis()->SetTitle("dE/dx"); // Y軸ラベルを設定
  hist_IBL_dEdxVsP[0]->Draw("colz"); // colzオプションを追加
  c1->cd(3);   c1->SetLogy(0);
  /*Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_HitSize[0],
             hist_IBL_HitSize[1],
             0,0,0.1,"hist_IBL_HitSize","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_HitSizePhi[0],
             hist_IBL_HitSizePhi[1],
             0,0,0.1,"hist_IBL_HitSizePhi","");
  c1->cd(4); alabel(2);*/
  c1->Update(); ps->NewPage(); c1->Clear();
  /*
  //========
  // page.17
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_HitSizeZ[0],
             hist_IBL_HitSizeZ[1],
             0,0,"hist_IBL_HitSizeZ","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_HitSizeRatio[0],
             hist_IBL_HitSizeRatio[1],
             0,0,"hist_IBL_HitSizeRatio","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_unbiasedResidualX[0],
             hist_IBL_unbiasedResidualX[1],
             0,0,"hist_IBL_unbiasedResidualX","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_unbiasedResidualY[0],
             hist_IBL_unbiasedResidualY[1],
             0,0,"hist_IBL_unbiasedResidualY","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.18
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_Isolation10x2[0],
             hist_IBL_Isolation10x2[1],
             0,0,"hist_IBL_Isolation10x2","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_Isolation20x4[0],
             hist_IBL_Isolation20x4[1],
             0,0,"hist_IBL_Isolation20x4","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_numTotalClustersPerModule[0],
             hist_IBL_numTotalClustersPerModule[1],
             0,0,"hist_IBL_numTotalClustersPerModule","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_numTotalPixelsPerModule[0],
             hist_IBL_numTotalPixelsPerModule[1],
             0,0,"hist_IBL_numTotalPixelsPerModule","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.19
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  c1->cd(2);   c1->SetLogy(0);
  hist_IBL_Map[0] -> Draw("colz");
  c1->cd(3);   c1->SetLogy(0);
  hist_IBL_MapHit[0] -> Draw("Box");
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_MapEta[0],
             hist_IBL_MapEta[1],
             0,0,"hist_IBL_MapEta","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.20
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_MapHitEta[0],
             hist_IBL_MapHitEta[1],
             0,0,"hist_IBL_MapHitEta","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_dEdx_true[0],
             hist_IBL_dEdx_true[1],
             0,0,"hist_IBL_dEdx_true","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  hist_IBL_dEdxVsP_true[0] -> Draw("Box");
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_IBL_ChargeFraction[0],
             hist_IBL_ChargeFraction[1],
             0,0,"hist_IBL_ChargeFraction","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.21
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  hist_IBL_ClusterShape[0] -> Draw("Box");
  c1->cd(2);   c1->SetLogy(0);
  hist_IBL_ModuleMap[0] -> Draw("Box");
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_IsEdge[0],
             hist_BLY_IsEdge[1],
             0,0,"hist_BLY_IsEdge","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_IsOverflow[0],
             hist_BLY_IsOverflow[1],
             0,0,"hist_BLY_IsOverflow","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.22
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_IsSplit[0],
             hist_BLY_IsSplit[1],
             0,0,"hist_BLY_IsSplit","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_L1A[0],
             hist_BLY_L1A[1],
             0,0,"hist_BLY_L1A","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_ToT[0],
             hist_BLY_ToT[1],
             0,0,"hist_BLY_ToT","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_Charge[0],
             hist_BLY_Charge[1],
             0,0,"hist_BLY_Charge","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.23
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_dEdx[0],
             hist_BLY_dEdx[1],
             0,0,"hist_BLY_dEdx","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  hist_BLY_dEdxVsP[0] -> Draw("Box");
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_HitSize[0],
             hist_BLY_HitSize[1],
             0,0,"hist_BLY_HitSize","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_HitSizePhi[0],
             hist_BLY_HitSizePhi[1],
             0,0,"hist_BLY_HitSizePhi","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.24
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_HitSizeZ[0],
             hist_BLY_HitSizeZ[1],
             0,0,"hist_BLY_HitSizeZ","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_unbiasedResidualX[0],
             hist_BLY_unbiasedResidualX[1],
             0,0,"hist_BLY_unbiasedResidualX","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_unbiasedResidualY[0],
             hist_BLY_unbiasedResidualY[1],
             0,0,"hist_BLY_unbiasedResidualY","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_Isolation10x2[0],
             hist_BLY_Isolation10x2[1],
             0,0,"hist_BLY_Isolation10x2","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.25
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_Isolation20x4[0],
             hist_BLY_Isolation20x4[1],
             0,0,"hist_BLY_Isolation20x4","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_numTotalClustersPerModule[0],
             hist_BLY_numTotalClustersPerModule[1],
             0,0,"hist_BLY_numTotalClustersPerModule","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_numTotalPixelsPerModule[0],
             hist_BLY_numTotalPixelsPerModule[1],
             0,0,"hist_BLY_numTotalPixelsPerModule","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.26
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  hist_BLY_Map[0] -> Draw("Box");
  c1->cd(2);   c1->SetLogy(0);
  hist_BLY_MapHit[0] -> Draw("Box");
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_MapEta[0],
             hist_BLY_MapEta[1],
             0,0,"hist_BLY_MapEta","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_MapHitEta[0],
             hist_BLY_MapHitEta[1],
             0,0,"hist_BLY_MapHitEta","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.27
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_dEdx_true[0],
             hist_BLY_dEdx_true[1],
             0,0,"hist_BLY_dEdx_true","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  hist_BLY_dEdxVsP_true[0] -> Draw("Box");
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_BLY_ChargeFraction[0],
             hist_BLY_ChargeFraction[1],
             0,0,"hist_BLY_ChargeFraction","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  hist_BLY_ClusterShape[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.28
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_IsEdge[0],
             hist_LY1_IsEdge[1],
             0,0,"hist_LY1_IsEdge","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_IsOverflow[0],
             hist_LY1_IsOverflow[1],
             0,0,"hist_LY1_IsOverflow","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_IsSplit[0],
             hist_LY1_IsSplit[1],
             0,0,"hist_LY1_IsSplit","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_L1A[0],
             hist_LY1_L1A[1],
             0,0,"hist_LY1_L1A","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.29
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_ToT[0],
             hist_LY1_ToT[1],
             0,0,"hist_LY1_ToT","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_Charge[0],
             hist_LY1_Charge[1],
             0,0,"hist_LY1_Charge","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_dEdx[0],
             hist_LY1_dEdx[1],
             0,0,"hist_LY1_dEdx","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  hist_LY1_dEdxVsP[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.30
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_HitSize[0],
             hist_LY1_HitSize[1],
             0,0,"hist_LY1_HitSize","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_HitSizePhi[0],
             hist_LY1_HitSizePhi[1],
             0,0,"hist_LY1_HitSizePhi","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_HitSizeZ[0],
             hist_LY1_HitSizeZ[1],
             0,0,"hist_LY1_HitSizeZ","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_unbiasedResidualX[0],
             hist_LY1_unbiasedResidualX[1],
             0,0,"hist_LY1_unbiasedResidualX","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.31
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_unbiasedResidualY[0],
             hist_LY1_unbiasedResidualY[1],
             0,0,"hist_LY1_unbiasedResidualY","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_Isolation10x2[0],
             hist_LY1_Isolation10x2[1],
             0,0,"hist_LY1_Isolation10x2","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_Isolation20x4[0],
             hist_LY1_Isolation20x4[1],
             0,0,"hist_LY1_Isolation20x4","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_numTotalClustersPerModule[0],
             hist_LY1_numTotalClustersPerModule[1],
             0,0,"hist_LY1_numTotalClustersPerModule","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.32
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_numTotalPixelsPerModule[0],
             hist_LY1_numTotalPixelsPerModule[1],
             0,0,"hist_LY1_numTotalPixelsPerModule","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  c1->cd(3);   c1->SetLogy(0);
  hist_LY1_Map[0] -> Draw("Box");
  c1->cd(4);   c1->SetLogy(0);
  hist_LY1_MapHit[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.33
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_MapEta[0],
             hist_LY1_MapEta[1],
             0,0,"hist_LY1_MapEta","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_MapHitEta[0],
             hist_LY1_MapHitEta[1],
             0,0,"hist_LY1_MapHitEta","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY1_dEdx_true[0],
             hist_LY1_dEdx_true[1],
             0,0,"hist_LY1_dEdx_true","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  hist_LY1_dEdxVsP_true[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.34
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_IsEdge[0],
             hist_LY2_IsEdge[1],
             0,0,"hist_LY2_IsEdge","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_IsOverflow[0],
             hist_LY2_IsOverflow[1],
             0,0,"hist_LY2_IsOverflow","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_IsSplit[0],
             hist_LY2_IsSplit[1],
             0,0,"hist_LY2_IsSplit","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_L1A[0],
             hist_LY2_L1A[1],
             0,0,"hist_LY2_L1A","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.35
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_ToT[0],
             hist_LY2_ToT[1],
             0,0,"hist_LY2_ToT","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_Charge[0],
             hist_LY2_Charge[1],
             0,0,"hist_LY2_Charge","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_dEdx[0],
             hist_LY2_dEdx[1],
             0,0,"hist_LY2_dEdx","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  hist_LY2_dEdxVsP[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.36
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_HitSize[0],
             hist_LY2_HitSize[1],
             0,0,"hist_LY2_HitSize","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_HitSizePhi[0],
             hist_LY2_HitSizePhi[1],
             0,0,"hist_LY2_HitSizePhi","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_HitSizeZ[0],
             hist_LY2_HitSizeZ[1],
             0,0,"hist_LY2_HitSizeZ","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_unbiasedResidualX[0],
             hist_LY2_unbiasedResidualX[1],
             0,0,"hist_LY2_unbiasedResidualX","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.37
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_unbiasedResidualY[0],
             hist_LY2_unbiasedResidualY[1],
             0,0,"hist_LY2_unbiasedResidualY","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_Isolation10x2[0],
             hist_LY2_Isolation10x2[1],
             0,0,"hist_LY2_Isolation10x2","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_Isolation20x4[0],
             hist_LY2_Isolation20x4[1],
             0,0,"hist_LY2_Isolation20x4","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_numTotalClustersPerModule[0],
             hist_LY2_numTotalClustersPerModule[1],
             0,0,"hist_LY2_numTotalClustersPerModule","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.38
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_numTotalPixelsPerModule[0],
             hist_LY2_numTotalPixelsPerModule[1],
             0,0,"hist_LY2_numTotalPixelsPerModule","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  c1->cd(3);   c1->SetLogy(0);
  hist_LY2_Map[0] -> Draw("Box");
  c1->cd(4);   c1->SetLogy(0);
  hist_LY2_MapHit[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.39
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_MapEta[0],
             hist_LY2_MapEta[1],
             0,0,"hist_LY2_MapEta","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_MapHitEta[0],
             hist_LY2_MapHitEta[1],
             0,0,"hist_LY2_MapHitEta","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_LY2_dEdx_true[0],
             hist_LY2_dEdx_true[1],
             0,0,"hist_LY2_dEdx_true","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  hist_LY2_dEdxVsP_true[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.40
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_IsEdge[0],
             hist_END_IsEdge[1],
             0,0,"hist_END_IsEdge","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_IsOverflow[0],
             hist_END_IsOverflow[1],
             0,0,"hist_END_IsOverflow","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_IsSplit[0],
             hist_END_IsSplit[1],
             0,0,"hist_END_IsSplit","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_L1A[0],
             hist_END_L1A[1],
             0,0,"hist_END_L1A","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.41
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_ToT[0],
             hist_END_ToT[1],
             0,0,"hist_END_ToT","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_Charge[0],
             hist_END_Charge[1],
             0,0,"hist_END_Charge","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_dEdx[0],
             hist_END_dEdx[1],
             0,0,"hist_END_dEdx","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  hist_END_dEdxVsP[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.42
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_HitSize[0],
             hist_END_HitSize[1],
             0,0,"hist_END_HitSize","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_HitSizePhi[0],
             hist_END_HitSizePhi[1],
             0,0,"hist_END_HitSizePhi","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_HitSizeZ[0],
             hist_END_HitSizeZ[1],
             0,0,"hist_END_HitSizeZ","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_unbiasedResidualX[0],
             hist_END_unbiasedResidualX[1],
             0,0,"hist_END_unbiasedResidualX","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.43
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_unbiasedResidualY[0],
             hist_END_unbiasedResidualY[1],
             0,0,"hist_END_unbiasedResidualY","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_Isolation10x2[0],
             hist_END_Isolation10x2[1],
             0,0,"hist_END_Isolation10x2","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_Isolation20x4[0],
             hist_END_Isolation20x4[1],
             0,0,"hist_END_Isolation20x4","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_numTotalClustersPerModule[0],
             hist_END_numTotalClustersPerModule[1],
             0,0,"hist_END_numTotalClustersPerModule","");
  c1->cd(4); alabel(2);
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.44
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_numTotalPixelsPerModule[0],
             hist_END_numTotalPixelsPerModule[1],
             0,0,"hist_END_numTotalPixelsPerModule","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  c1->cd(3);   c1->SetLogy(0);
  hist_ED1_Map[0] -> Draw("Box");
  c1->cd(4);   c1->SetLogy(0);
  hist_ED1_MapHit[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.45
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_ED1_MapEta[0],
             hist_ED1_MapEta[1],
             0,0,"hist_ED1_MapEta","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_ED1_MapHitEta[0],
             hist_ED1_MapHitEta[1],
             0,0,"hist_ED1_MapHitEta","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  hist_ED2_3Map[0] -> Draw("Box");
  c1->cd(4);   c1->SetLogy(0);
  hist_ED2_3MapHit[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.46
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_ED2_3MapEta[0],
             hist_ED2_3MapEta[1],
             0,0,"hist_ED2_3MapEta","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_ED2_3MapHitEta[0],
             hist_ED2_3MapHitEta[1],
             0,0,"hist_ED2_3MapHitEta","");
  c1->cd(2); alabel(2);
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_END_dEdx_true[0],
             hist_END_dEdx_true[1],
             0,0,"hist_END_dEdx_true","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  hist_END_dEdxVsP_true[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
  //========
  // page.47
  //========
  c1 -> Divide(2,2);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_ALL_dEdx[0],
             hist_ALL_dEdx[1],
             0,0,"hist_ALL_dEdx","");
  c1->cd(1); alabel(2);
  c1->cd(2);   c1->SetLogy(0);
  hist_ALL_dEdxVsP[0] -> Draw("Box");
  c1->cd(3);   c1->SetLogy(0);
  Hist1DPlot(0,1,4,1,1,1,1,
             hist_ALL_dEdx_true[0],
             hist_ALL_dEdx_true[1],
             0,0,"hist_ALL_dEdx_true","");
  c1->cd(3); alabel(2);
  c1->cd(4);   c1->SetLogy(0);
  hist_ALL_dEdxVsP_true[0] -> Draw("Box");
  c1->Update(); ps->NewPage(); c1->Clear();
*/
  for (int i = 0; i < nfil; i++) {

    /*hist_muValue[i] -> Clear();
    hist_lumiBlk[i] -> Clear();
    hist_numVtx[i] -> Clear();
    hist_numJet[i] -> Clear();
    hist_numJetjvt[i] -> Clear();
    hist_muon_pt1[i] -> Clear();
    hist_muon_eta1[i] -> Clear();
    hist_muon_pt2[i] -> Clear();
    hist_muon_eta2[i] -> Clear();
    hist_dimuon_mass[i] -> Clear();
    hist_dimuon_pt[i] -> Clear();
    hist_dimuon_eta[i] -> Clear();
    hist_trk_num[i] -> Clear();
    hist_trk_pt[i] -> Clear();
    hist_trk_eta[i] -> Clear();
    hist_trk_phi[i] -> Clear();
    hist_trk_qoverp[i] -> Clear();
    hist_trk_d0[i] -> Clear();
    hist_trk_z0[i] -> Clear();
    hist_trk_deltaZ[i] -> Clear();
    hist_trk_etaVSphi[i] -> Clear();
    hist_truth_pt[i] -> Clear();
    hist_truth_eta[i] -> Clear();
    hist_truth_phi[i] -> Clear();
    hist_truth_d0[i] -> Clear();
    hist_truth_z0[i] -> Clear();
    hist_trk_dpt[i] -> Clear();
    hist_trk_dphi[i] -> Clear();
    hist_trk_deta[i] -> Clear();
    hist_trk_dd0[i] -> Clear();
    hist_trk_dz0[i] -> Clear();
    hist_trk_nPixHits[i] -> Clear();
    hist_trk_nSCTHits[i] -> Clear();
    hist_trk_nSiHits[i] -> Clear();
    hist_trk_nGangedPix[i] -> Clear();
    hist_trk_nPixLay[i] -> Clear();
    hist_trk_nPixSharedHits[i] -> Clear();
    hist_trk_nPixSplitHits[i] -> Clear();
    hist_trk_nPixOutliers[i] -> Clear();
    hist_trk_nPixHoles[i] -> Clear();
    hist_trk_nPixelDeadSensors[i] -> Clear();
    hist_trk_nSCTSharedHits[i] -> Clear();
    hist_trk_nSCTOutliers[i] -> Clear();
    hist_trk_nSCTHoles[i] -> Clear();
    hist_trk_nSCTDeadSensors[i] -> Clear();
    hist_trk_nTRTHits[i] -> Clear();
    hist_trk_nTRTOutliers[i] -> Clear();
    hist_trk_nTRTHoles[i] -> Clear();
    hist_trk_nTRTHTHits[i] -> Clear();
    hist_trk_chiSqPerDof[i] -> Clear();
    hist_trk_nOutliers[i] -> Clear();
    hist_trk_stdDevChi2OS[i] -> Clear();
    hist_trk_truthMatchProb[i] -> Clear();
    hist_trk_weight[i] -> Clear();
    hist_IBL_IsEdge[i] -> Clear();
    hist_IBL_IsOverflow[i] -> Clear();
    hist_IBL_IsSplit[i] -> Clear();
    hist_IBL_L1A[i] -> Clear();
    hist_IBL_ToT[i] -> Clear();
    hist_IBL_Charge[i] -> Clear();*/
    hist_IBL_dEdx[i] -> Clear();
    hist_IBL_dEdxVsP[i] -> Clear();
    hist_IBL_HitSize[i] -> Clear();
    hist_IBL_HitSizePhi[i] -> Clear();
    /*hist_IBL_HitSizeZ[i] -> Clear();
    hist_IBL_HitSizeRatio[i] -> Clear();
    hist_IBL_unbiasedResidualX[i] -> Clear();
    hist_IBL_unbiasedResidualY[i] -> Clear();
    hist_IBL_Isolation10x2[i] -> Clear();
    hist_IBL_Isolation20x4[i] -> Clear();
    hist_IBL_numTotalClustersPerModule[i] -> Clear();
    hist_IBL_numTotalPixelsPerModule[i] -> Clear();
    profile_IBL_LorentzAngle[i] -> Clear();
    hist_IBL_Map[i] -> Clear();
    hist_IBL_MapHit[i] -> Clear();
    hist_IBL_MapEta[i] -> Clear();
    hist_IBL_MapHitEta[i] -> Clear();
    hist_IBL_dEdx_true[i] -> Clear();
    hist_IBL_dEdxVsP_true[i] -> Clear();
    hist_IBL_ChargeFraction[i] -> Clear();
    hist_IBL_ClusterShape[i] -> Clear();
    hist_IBL_ModuleMap[i] -> Clear();
    hist_BLY_IsEdge[i] -> Clear();
    hist_BLY_IsOverflow[i] -> Clear();
    hist_BLY_IsSplit[i] -> Clear();
    hist_BLY_L1A[i] -> Clear();
    hist_BLY_ToT[i] -> Clear();
    hist_BLY_Charge[i] -> Clear();
    hist_BLY_dEdx[i] -> Clear();
    hist_BLY_dEdxVsP[i] -> Clear();
    hist_BLY_HitSize[i] -> Clear();
    hist_BLY_HitSizePhi[i] -> Clear();
    hist_BLY_HitSizeZ[i] -> Clear();
    hist_BLY_unbiasedResidualX[i] -> Clear();
    hist_BLY_unbiasedResidualY[i] -> Clear();
    hist_BLY_Isolation10x2[i] -> Clear();
    hist_BLY_Isolation20x4[i] -> Clear();
    hist_BLY_numTotalClustersPerModule[i] -> Clear();
    hist_BLY_numTotalPixelsPerModule[i] -> Clear();
    profile_BLY_LorentzAngle[i] -> Clear();
    hist_BLY_Map[i] -> Clear();
    hist_BLY_MapHit[i] -> Clear();
    hist_BLY_MapEta[i] -> Clear();
    hist_BLY_MapHitEta[i] -> Clear();
    hist_BLY_dEdx_true[i] -> Clear();
    hist_BLY_dEdxVsP_true[i] -> Clear();
    hist_BLY_ChargeFraction[i] -> Clear();
    hist_BLY_ClusterShape[i] -> Clear();
    hist_LY1_IsEdge[i] -> Clear();
    hist_LY1_IsOverflow[i] -> Clear();
    hist_LY1_IsSplit[i] -> Clear();
    hist_LY1_L1A[i] -> Clear();
    hist_LY1_ToT[i] -> Clear();
    hist_LY1_Charge[i] -> Clear();
    hist_LY1_dEdx[i] -> Clear();
    hist_LY1_dEdxVsP[i] -> Clear();
    hist_LY1_HitSize[i] -> Clear();
    hist_LY1_HitSizePhi[i] -> Clear();
    hist_LY1_HitSizeZ[i] -> Clear();
    hist_LY1_unbiasedResidualX[i] -> Clear();
    hist_LY1_unbiasedResidualY[i] -> Clear();
    hist_LY1_Isolation10x2[i] -> Clear();
    hist_LY1_Isolation20x4[i] -> Clear();
    hist_LY1_numTotalClustersPerModule[i] -> Clear();
    hist_LY1_numTotalPixelsPerModule[i] -> Clear();
    profile_LY1_LorentzAngle[i] -> Clear();
    hist_LY1_Map[i] -> Clear();
    hist_LY1_MapHit[i] -> Clear();
    hist_LY1_MapEta[i] -> Clear();
    hist_LY1_MapHitEta[i] -> Clear();
    hist_LY1_dEdx_true[i] -> Clear();
    hist_LY1_dEdxVsP_true[i] -> Clear();
    hist_LY2_IsEdge[i] -> Clear();
    hist_LY2_IsOverflow[i] -> Clear();
    hist_LY2_IsSplit[i] -> Clear();
    hist_LY2_L1A[i] -> Clear();
    hist_LY2_ToT[i] -> Clear();
    hist_LY2_Charge[i] -> Clear();
    hist_LY2_dEdx[i] -> Clear();
    hist_LY2_dEdxVsP[i] -> Clear();
    hist_LY2_HitSize[i] -> Clear();
    hist_LY2_HitSizePhi[i] -> Clear();
    hist_LY2_HitSizeZ[i] -> Clear();
    hist_LY2_unbiasedResidualX[i] -> Clear();
    hist_LY2_unbiasedResidualY[i] -> Clear();
    hist_LY2_Isolation10x2[i] -> Clear();
    hist_LY2_Isolation20x4[i] -> Clear();
    hist_LY2_numTotalClustersPerModule[i] -> Clear();
    hist_LY2_numTotalPixelsPerModule[i] -> Clear();
    profile_LY2_LorentzAngle[i] -> Clear();
    hist_LY2_Map[i] -> Clear();
    hist_LY2_MapHit[i] -> Clear();
    hist_LY2_MapEta[i] -> Clear();
    hist_LY2_MapHitEta[i] -> Clear();
    hist_LY2_dEdx_true[i] -> Clear();
    hist_LY2_dEdxVsP_true[i] -> Clear();
    hist_END_IsEdge[i] -> Clear();
    hist_END_IsOverflow[i] -> Clear();
    hist_END_IsSplit[i] -> Clear();
    hist_END_L1A[i] -> Clear();
    hist_END_ToT[i] -> Clear();
    hist_END_Charge[i] -> Clear();
    hist_END_dEdx[i] -> Clear();
    hist_END_dEdxVsP[i] -> Clear();
    hist_END_HitSize[i] -> Clear();
    hist_END_HitSizePhi[i] -> Clear();
    hist_END_HitSizeZ[i] -> Clear();
    hist_END_unbiasedResidualX[i] -> Clear();
    hist_END_unbiasedResidualY[i] -> Clear();
    hist_END_Isolation10x2[i] -> Clear();
    hist_END_Isolation20x4[i] -> Clear();
    hist_END_numTotalClustersPerModule[i] -> Clear();
    hist_END_numTotalPixelsPerModule[i] -> Clear();
    profile_END_LorentzAngle[i] -> Clear();
    hist_ED1_Map[i] -> Clear();
    hist_ED1_MapHit[i] -> Clear();
    hist_ED1_MapEta[i] -> Clear();
    hist_ED1_MapHitEta[i] -> Clear();
    hist_ED2_3Map[i] -> Clear();
    hist_ED2_3MapHit[i] -> Clear();
    hist_ED2_3MapEta[i] -> Clear();
    hist_ED2_3MapHitEta[i] -> Clear();
    hist_END_dEdx_true[i] -> Clear();
    hist_END_dEdxVsP_true[i] -> Clear();
    hist_ALL_dEdx[i] -> Clear();
    hist_ALL_dEdxVsP[i] -> Clear();
    hist_ALL_dEdx_true[i] -> Clear();
    hist_ALL_dEdxVsP_true[i] -> Clear();*/

    ffl[i] -> Close();
  }

  ps -> Close();
}

