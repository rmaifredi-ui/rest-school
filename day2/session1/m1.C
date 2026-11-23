#include <stdio.h>
#include <vector>
#include <TRestTools.h>
#include <TH1D.h>
#include <string>
#include <TCanvas.h>
#include <TString.h>

Int_t GenerateHistogram( std::string fname ) {

    /// Reading the table
    std::vector <std::vector <float> > myTable;
    TRestTools::ReadBinaryTable( "../../data/XenonNeon_25Pct_1bar.N150f", myTable );

    /// Creating the histogram to be filled
    TH1D *h = new TH1D( "Response", "XenonNeon_25Pct", myTable.size(), 0, 15);

    int n = 0;
    for( const auto &d : myTable ) {

        /// We add all the values in a row to set the value of each histogram bin
        float result = 0;
        for( const auto &elem : d )
            result += elem;

        n++; // First bin is 1

		/// We fill the histogram
        h->SetBinContent( n, result );
    }

    std::string fnameRoot = TRestTools::GetFileNameRoot( fname );

	/// We must create the canvas where to draw before drawing!!
    TCanvas c;

    h->Draw();
    h->GetXaxis()->SetTitle("Energy root [keV]");
    h->GetYaxis()->SetTitle("Efficiency");

    c.Print( (TString) fnameRoot + ".png" );

    return 0;
}
