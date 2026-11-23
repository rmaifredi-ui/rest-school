#include <ROOT/RDataFrame.hxx>
#include <TH2D.h>
#include <TEllipse.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <iostream>
#include <TSpectrum.h>


using namespace std;
using namespace ROOT;

double gauss_2bg(double* x, double* pars) {
    double t = (*x - 4000.) / 1000.;
    double bckg = pars[3] + pars[4] * t + pars[5] * t * t;
    double gauss = pars[0] * TMath::Exp(-TMath::Power(*x - pars[1], 2) / (2 * TMath::Power(pars[2], 2)));
    return gauss + bckg;
}

void calibrate() {

    // 1. Load the data and keep only central events (r < 10 mm)
    RDataFrame df("AnalysisTree", "../../data/R01855/R01855_0*root");
    auto df_filtered = df.Define("radius", "sqrt(hitsAna_xMean*hitsAna_xMean + hitsAna_yMean*hitsAna_yMean)")
                         .Filter("radius < 10");

    // 2. Energy histogram (model: name, title, nbins, xmin, xmax)
    auto h1 = df_filtered.Histo1D({"h", "...", 250, 0, 10000}, "tckAna_MaxTrackEnergy");
    TH1D h = h1.GetValue();      // h is now a real TH1D (a copy)
    
    //3. Peak search with TSpectrum
    TSpectrum *s = new TSpectrum();
    int nfound = s->Search(&h, 10, "", 0.5);    
    printf("Found %d candidate peaks to fit\n", nfound);
    if (nfound == 0) {
        cout << "No peaks found, exiting" << endl;
        return;
    }
    double* xpeaks = s->GetPositionX(); // returns the x position of the found peaks
    double peak = xpeaks[0]; // take the first peak as the main peak    
    for (int i = 1; i < nfound; i++)
        if (h.GetBinContent(h.FindBin(xpeaks[i])) > h.GetBinContent(h.FindBin(peak))) peak = xpeaks[i];
    cout << "Fe55 peak candidate at " << peak << " ADC" << endl;

    // --- Fit with start parameters taken from the peak ---
    TF1* f = new TF1("f", gauss_2bg, 0, 10000, 6);
    TF1* f_gauss = new TF1("f_gauss", "gaus", 0, 10000);
    f->SetParNames("Constant", "Mean", "Sigma", "Bkg_const", "Bkg_lin", "Bkg_quad");
    f->SetParameters(h.GetBinContent(h.FindBin(peak)),   // amplitude = height at the peak
                     peak,                                 // mean
                     0.1 * peak,                           // sigma ~ 10% of the position
                     0, 0, 0);                             // background starts at zero


    f->SetLineColor(kRed);
    f->SetLineWidth(2);

    f_gauss->SetLineColor(kAzure + 1);
    f_gauss->SetLineWidth(2);
    f_gauss->SetLineStyle(2);          // dashed: still distinguishable in black & white

    auto fitres = h.Fit(f, "LS", "", 0, 10000);
    auto fitres2 = h.Fit(f_gauss, "LS+", "", 0.5*peak, 1.5*peak);
    double mean = fitres->Parameter(1);
    cout << "Mean = " << mean << " +- " << fitres->ParError(1) << " ADC" << endl;
    cout << "Mean_gauss = " << fitres2->Parameter(1) << " +- " << fitres2->ParError(1) << " ADC" << endl;

    double calibration_factor = 5.899 / fitres2->Parameter(1); // 5.899 keV is the energy of the Fe55 peak
    cout << "Calibration factor = " << calibration_factor << " keV/ADC" << endl;

    // 4. Draw
    TCanvas* c = new TCanvas("c", "Calibration", 1000, 800);
    h.GetXaxis()->SetTitle("Energy [ADC]");
    h.GetYaxis()->SetTitle("pdf(E) [1/ADC]");
    h.SetMarkerStyle(23);
    h.SetMarkerColor(TColor::GetColor("#007700"));
    h.SetMarkerSize(0.8);
    h.SetStats(0);


    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(f, "Gauss + pol2 bkg", "l");
    leg->AddEntry(f_gauss, "Gauss only", "l");
    leg->Draw();
    
    h.DrawClone("E1");
    c->SaveAs("calibration_pdf.png");
    
    auto df_calibrated = df_filtered.Define("tckAna_MaxTrackEnergy_keV", [calibration_factor](double energy_adc) {
        return energy_adc * calibration_factor;
    }, {"tckAna_MaxTrackEnergy"});

    auto h2 = df_calibrated.Histo1D({"h2", "Calibrated data", 100, 0, 10}, "tckAna_MaxTrackEnergy_keV");
    TH1D h_calibrated = h2.GetValue();
    TCanvas* c2 = new TCanvas("c2", "Calibrated data", 1000, 800);
    h_calibrated.GetXaxis()->SetTitle("Energy [keV]");
    h_calibrated.GetYaxis()->SetTitle("pdf(E) [1/keV]");
    h_calibrated.SetLineWidth(2);
    h_calibrated.SetLineColor(kBlack);
    h_calibrated.SetFillColor(kCyan - 7);
    h_calibrated.DrawClone();
    c2->SaveAs("calibrated_data.png");
}


void heatmap(double calibration_factor) {
    RDataFrame df("AnalysisTree", "../../data/R01855/R01855_0*root");

    auto df_cal = df.Define("radius", "sqrt(hitsAna_xMean*hitsAna_xMean + hitsAna_yMean*hitsAna_yMean)")
                    .Define("energy_keV", [calibration_factor](double e) { return e * calibration_factor; },
                            {"tckAna_MaxTrackEnergy"});

    // Photopeak events only, for the gain map (window around 5.9 keV)
    auto df_peak = df_cal.Filter("energy_keV > 4 && energy_keV < 8");

    // Left: hit map of all events. Right: mean energy per 2 mm pixel (gain map)
    auto h1 = df_cal.Histo2D({"h1", "Hit map, all events", 60, -30, 30, 60, -30, 30},
                             "hitsAna_xMean", "hitsAna_yMean");
    auto p2 = df_peak.Profile2D({"p2", "Mean energy, 4-8 keV events", 30, -30, 30, 30, -30, 30},
                                "hitsAna_xMean", "hitsAna_yMean", "energy_keV");

    TH2D       hitmap  = h1.GetValue();
    TProfile2D gainmap = p2.GetValue();

    hitmap.SetTitle("Hit map, all events;x [mm];y [mm];counts");
    gainmap.SetTitle("Mean energy, 4-8 keV events;x [mm];y [mm];<E> [keV]");
    hitmap.SetStats(0);
    gainmap.SetStats(0);
    gainmap.SetMinimum(5.0);       // fix the colour scale around the peak
    gainmap.SetMaximum(6.8);

    gStyle->SetPalette(kViridis);

    TCanvas* c = new TCanvas("c_heat", "Hit and gain maps", 1600, 700);
    c->Divide(2, 1);

    c->cd(1);
    gPad->SetRightMargin(0.15);
    hitmap.DrawClone("COLZ");
    TEllipse* cut = new TEllipse(0, 0, 10, 10);
    cut->SetFillStyle(0);
    cut->SetLineColor(kRed);
    cut->SetLineWidth(2);
    cut->Draw();

    c->cd(2);
    gPad->SetRightMargin(0.15);
    gainmap.DrawClone("COLZ");
    cut->Draw();

    c->Modified();
    c->Update();
    c->SaveAs("heatmap.png");
}