#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TFitResult.h>

// Gaussian
double fitf(double* x, double* par) {
    return par[0] * TMath::Exp(-TMath::Power(*x - par[1], 2) / (2 * TMath::Power(par[2], 2)));
}

// Gaussian + quadratic background.
// The background is written in t = (x - 4000)/1000 so all coefficients are O(counts).
//   par[0..2] = Gaussian constant, mean, sigma
//   par[3..5] = background constant, linear, quadratic (in t)
double myfitbg(double* x, double* par) {
    double t = (*x - 4000.) / 1000.;
    double bckg = par[3] + par[4] * t + par[5] * t * t;
    return fitf(x, par) + bckg;
}

int fit_calibration() {

    ROOT::RDataFrame df("AnalysisTree", "../../data/R01855/R01855_0*root");
    auto h = df.Histo1D({"h", "calibration", 250, 0, 10000}, "tckAna_MaxTrackEnergy");

    TCanvas* c = new TCanvas();
    h->GetXaxis()->SetTitle("energy [ADC]");
    h->GetYaxis()->SetTitle("counts");

    // 1) Gaussian only
    TF1* func = new TF1("fit", fitf, 0, 10000, 3);
    func->SetParameters(200, 4000, 800);
    func->SetParLimits(1, 0, 10000);
    func->SetParNames("Constant", "Mean_value", "Sigma");
    func->SetLineColor(kRed);
    auto r1 = h->Fit(func, "LS", "", 2000, 6000);

    // 2) Gaussian + background, starting from the Gaussian result
    TF1* myfitf = new TF1("myfit", myfitbg, 0, 10000, 6);
    myfitf->SetParameters(func->GetParameter(0), func->GetParameter(1), func->GetParameter(2), 0, 0, 0);
    myfitf->SetParLimits(1, 0, 10000);
    myfitf->SetParNames("Constant", "Mean_value", "Sigma", "Bkg_const", "Bkg_lin", "Bkg_quad");
    myfitf->SetLineColor(kBlue);
    auto r2 = h->Fit(myfitf, "LS+", "", 0, 10000);   // "+" keeps the first function too

    // Compare
    printf("Gauss only : mean = %.1f +- %.1f, sigma = %.1f, chi2/ndf = %.2f\n",
           r1->Parameter(1), r1->ParError(1), r1->Parameter(2), r1->Chi2() / r1->Ndf());
    printf("Gauss + bkg: mean = %.1f +- %.1f, sigma = %.1f, chi2/ndf = %.2f\n",
           r2->Parameter(1), r2->ParError(1), r2->Parameter(2), r2->Chi2() / r2->Ndf());

    h->DrawClone();
    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(func, "Gaussian", "l");
    leg->AddEntry(myfitf, "Gaussian + pol2 bkg", "l");
    leg->Draw();

    c->SaveAs("fit_calibration_bkg.png");
    return 0;
}