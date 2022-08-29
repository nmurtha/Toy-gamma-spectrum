// NATHAN MURTHA, Jan 2017
//  Updated: June 17, 2018
//
// This code generates a toy spectra of incoming gamma rays of a specified
// energy. The code does not consider pair production, so input energies larger
// than 1.022 MeV will miss this effect. The code also assumes an infinitely
// large interaction medium, always deciding between Compton scatter or
// photoelectric absorption.

// Declare internal function headers.
double applySmearSPM(double E);

// Start main program.
void ToySpectrum()
{
    // Set seed for consistency.
    gRandom->SetSeed(12345);

    // Set number of photons to simulate and their max energy, reading these
    // values in from the user. Also read in and set a threshold energy below
    // which tracking of the photon completely ceases (whatever interaction
    // the photon last underwent is considered it's last interaction).
    int nPhotons;
    double Emax;    // [MeV]
    double Ecut;    // [MeV]
    cout << "How many photons will be simulated?" << endl;
    cin >> nPhotons;
    cout << "What energy [MeV] do the photons have?" << endl;
    cin >> Emax;
    cout << "To what energy [MeV] will photons be considered before calculations are terminated?" << endl;
    cin >> Ecut;
    if (Emax > 1.022)
    {
        printf(" WARNING: Toy simulation does not consider pair production,"
            "and you've set Emax > 1.022 MeV.\n");
    }

    // Begin assigning either absorptions or Compton scatters. Follow photon
    // until energy falls below threshold.
    TH1D* hTransSmear = new TH1D("","hTransSmear",200,0,Emax+0.1); // energy transfer with empirical energy smearing
    TH1D* hTrans      = new TH1D("","hTrans",200,0,Emax+0.1); // energy transfer
    TH1D* hLDcompton  = new TH1D("","hLDcompton",10,0,10); // number of photons with last interaction as Compton scatter
    TH1D* hLDdeposit  = new TH1D("","hLDdeposit",10,0,10); // number of photons with last interaction as photoelectric absorption
    double theta, Ephoton, Etrans;
    int nInters = 0, nHistory, flag;
    for (int iPhoton = 0; iPhoton < nPhotons; iPhoton++)
    {
        Ephoton  = Emax; // new photon starts with full energy
        nHistory = 0; // new photon has 0 interactions
        while (Ephoton > Ecut)
        {
            if (gRandom->Uniform(0,1) <= 0.5) // total deposition
            {
                nHistory++; // individual photon interaction index goes up by one
                flag = 0;
                Etrans = Ephoton;
                hTransSmear->Fill(applySmearSPM(Etrans)); // [MeV]
                hTrans->Fill(Etrans); // [MeV]
                Ephoton -= Etrans;
            }
            else // Compton scatter
            {
                nHistory++; // individual photon interaction index goes up by one
                flag = 1;
                theta  = gRandom->Uniform(0,TMath::Pi());
                Etrans = ( Ephoton*Ephoton*(1-cos(theta))  )
                    / ( Ephoton*(1-cos(theta)) + 0.511 ); // [MeV]
                hTransSmear->Fill(applySmearSPM(Etrans));
                hTrans->Fill(Etrans);
                Ephoton -= Etrans;
            }
            nInters += 1; // overall interaction index goes up by one
        }

        // Check status of flag to determine what photon last did. Add number of
        // interactions to appropriate histogram (either last interaction was
        // photoelectric absorption or was Compton scatter).
        if (flag == 0)
        {
            hLDdeposit->Fill(nHistory);
        }
        else if (flag == 1)
        {
            hLDcompton->Fill(nHistory);
        }
    }

    // Display the results.
    TCanvas* c1 = new TCanvas();
    c1->cd(1);
    c1->SetGrid();
    //c1->SetLogy();
    hTrans->SetTitle(Form("Spectra from %d %.4f MeV Photons (%d Interactions)",nPhotons,Emax,nInters));
    hTrans->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    hTrans->GetYaxis()->SetTitle("Counts");
    hTrans->SetStats(0);
    hTrans->Draw();
    hTransSmear->SetLineColor(kRed);
    hTransSmear->Draw("Same");

    TLegend *legend = new TLegend(0.15, 0.85, 0.4, 0.65, "", "brNDC");
    legend->AddEntry(hTrans, Form("No Energy Smearing"), "l");
    legend->AddEntry(hTransSmear, Form("Energy Smearing Applied"), "l");
    legend->Draw("same");

    TCanvas* c4 = new TCanvas();
    c4->cd(1);
    hLDdeposit->SetTitle("Number of Interactions Before Termination (Last Interaction as Photoelectric Absorption)");
    hLDdeposit->Draw();

    TCanvas* c5 = new TCanvas();
    c5->cd(1);
    hLDcompton->SetTitle("Number of Interactions Before Termination (Last Interaction as Compton Scatter)");
    hLDcompton->Draw();


}


// -----------------------------------------------------------------------------
// -----                        INTERNAL FUNCTIONS                         -----
// -----------------------------------------------------------------------------

// Empirical smearing function that takes an input energy in MeV and randomly
// resamples a Gaussian distribution centered at the nominal energy to model the
// effect of imperfect energy resolution.
double applySmearSPM(double E)
{
    // Set experimental coefficients for the energy smearing.
    double C_SPM[3] = {6.49750e+00,2.39530e-01, 2.00423e-01};

    // Limit the smearing at low, E, esp for E=0.000000 !
    if (E<0.001)
    {
        E=0.001;
    }

    double resFWHMperc = (C_SPM[0]+C_SPM[1]*E+C_SPM[2]*E*E) / sqrt(E);
    double sig = resFWHMperc/2.35/100*E;
    double Esmear = gRandom->Gaus(E,sig);

    // Catch negative energies.
    if (Esmear<=0)
    {
        Esmear=0.001;
    }

    // Return smeared energy.
    return Esmear;
}
