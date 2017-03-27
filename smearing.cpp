/*
 *--------------------------------------------------------------
 * An exercise to smear and attempt to unfold a simple histrogram 
 *
 *--------------------------------------------------------------
 * Author: Rhiannon Jones
 * Date  : March 2017
 *
 * For now, this is a simple macro for neutrino-mode only
 * I will amend this entire system to be a class structure with 
 * constructors and destructors for neutrino and antineutrino-modes
 *
*/

#include "smearing.h"

using namespace std; 

int smearing() {

    // -------------------------------------------------------------------------
    //                              Open event files
    // -------------------------------------------------------------------------
    
    // Open Default
    TFile f1("gntp.10000.gst.root");
    if (f1.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " Default event file is open " << endl;
    }
    //-----------------------------------------------------------------------
    //                            Parameters to consider
    // -------------------------------------------------------------------------
    
    // Number of MiniBooNE bins (data) 
    // costhetamu : 20
    // Tmu        : 18
    // 
    // Amount to smear 
    // costhetamu : 5 degrees
    // Tmu        : 10%
    //
    // Kinetic energy cuts
    // mu         : > 50MeV
    // charged pi : > 50MeV
   
    // -------------------------------------------------------------------------
    //                  Get the tree, normalisation and the branches
    // -------------------------------------------------------------------------
    TTree *def_tree = (TTree*) f1.Get( "gst" );
 
    // Define the log and normal histograms for the unsmeared distributions
    // El - m_mu : kinetic energy of the final state primary lepton (fspl : muon == 13 )
    // cthl      : costheta of the final state primary lepton
    // h_un      : unsmeared histogram
    // h_sm      : smeared histogram
    //

    TFile *f = new TFile("output.root","RECREATE");

    TH2D *h_un = new TH2D("h_un"," T_{#mu} - cos#theta_{#mu} distribution before smearing ",20,-1,1,18,0,2);
    def_tree->Draw("( El - 0.10566 ):cthl>>h_un","fspl == 13 && cc","colz"); 

    h_un->SetStats(kFALSE);
    h_un->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un->GetYaxis()->SetTitle("T_{#mu}");
    h_un->Write("distribution_before");

    std::vector<TH2*> v_un;
    std::vector<TH2*> v_sm;
    
    // The same histogram definitions for the smeared distributions
    TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} distribution after smearing, impurity ",20,-1,1,18,0,2);
    
    // Take h_un and smear it
    Smear(def_tree, h_un, h_sm, v_un, v_sm, f);

    // The the 2 2D histograms and draw slices in Tmu and cos theta mu
    Slices( h_un, h_sm );

    SliceStack( v_un, v_sm, f );

    f->Close();

    return 0;
}

// -------------------------------------------------------------------------
//                            Function definitions
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
//                             Smearing function
// -------------------------------------------------------------------------

void Smear ( TTree *tree, 
             TH2D  *h_unsmeared,
             TH2D  *h_smeared, 
             std::vector<TH2*> &v_un,
             std::vector<TH2*> &v_sm,
             TFile *f ){

    // Total momentum and angle distributions before smearing
    TH1D *h_Tmu_un = new TH1D("h_Tmu_un","T_{#mu} before smearing",18,0,2);
    TH1D *h_costhetamu_un = new TH1D("h_costheta_mu","cos(#theta_{#mu}) before smearing",20,-1,1);

    // Total momentum and angle distributions after smearing
    TH1D *h_Tmu_sm = new TH1D("h_Tmu_sm","T_{#mu} after smearing",18,0,2);
    TH1D *h_costhetamu_sm = new TH1D("h_costheta_sm","cos(#theta_{#mu}) after smearing",20,-1,1);

    // Stacked unsmeared histogram
    THStack *hs_un = new THStack("hs_un","Stacked Histograms");
    TH2D *h_ccqe_un = new TH2D("h_ccqe_un","CCQE",20,-1,1,18,0,2);
    TH2D *h_ccres_un = new TH2D("h_ccres_un","CCRES",20,-1,1,18,0,2);
    TH2D *h_ccdis_un = new TH2D("h_ccdis_un","CCDIS",20,-1,1,18,0,2);
    TH2D *h_cccoh_un = new TH2D("h_cccoh_un","CCCOH",20,-1,1,18,0,2);

    // Stacked smeared histogram
    THStack *hs_sm = new THStack("hs_sm","Stacked Histograms");
    TH2D *h_ccqe_sm = new TH2D("h_ccqe_sm","CCQE",20,-1,1,18,0,2);
    TH2D *h_ccres_sm = new TH2D("h_ccres_sm","CCRES",20,-1,1,18,0,2);
    TH2D *h_ccdis_sm = new TH2D("h_ccdis_sm","CCDIS",20,-1,1,18,0,2);
    TH2D *h_cccoh_sm = new TH2D("h_cccoh_sm","CCCOH",20,-1,1,18,0,2);
    TH2D *h_nc_sm = new TH2D("h_nc_sm","NC1pi",20,-1,1,18,0,2);

    ofstream file1;
    file1.open("unsmeared_results.tex"); 
    file1 << " \\begin{tabular}{ | * {20}{c} | } " << endl;
    
    ofstream file2;
    file2.open("impure_smeared_results.tex"); 
    file2 << " \\begin{tabular}{ | * {20}{c} | } " << endl;
     
    ofstream file3;
    file3.open("impure_results_difference.tex"); 
    file3 << " \\begin{tabular}{ | * {20}{c} | } " << endl;

    // Get the branches
    TBranch *b_nf    = tree->GetBranch( "nf" );
    TBranch *b_cthl  = tree->GetBranch( "cthl" );
    TBranch *b_El    = tree->GetBranch( "El" );
    TBranch *b_pl    = tree->GetBranch( "pl" );
    TBranch *b_fspl  = tree->GetBranch( "fspl" );
    TBranch *b_cthf  = tree->GetBranch( "cthf" );
    TBranch *b_Ef    = tree->GetBranch( "Ef" );
    TBranch *b_pf    = tree->GetBranch( "pf" );
    TBranch *b_pdgf  = tree->GetBranch( "pdgf" );
    TBranch *b_cc    = tree->GetBranch( "cc" );
    TBranch *b_nc    = tree->GetBranch( "nc" );
    TBranch *b_nfpip = tree->GetBranch( "nfpip" );
    TBranch *b_nfpim = tree->GetBranch( "nfpim" );
    TBranch *b_nfpi0 = tree->GetBranch( "nfpi0" );
    TBranch *b_qel   = tree->GetBranch( "qel" );
    TBranch *b_res   = tree->GetBranch( "res" );
    TBranch *b_dis   = tree->GetBranch( "dis" );
    TBranch *b_coh   = tree->GetBranch( "coh" );
    
    // Number of events in the TTree
    int n_values = tree->GetEntries();
   
    double m_mu = 0.10566; // Muon mass, GeV
    double m_pi = 0.13957; // Charged pion mass, GeV
    
    // Initiate some counters for when things went wrong
    int e_count     = 0;
    int other_count = 0;
    int pion_count  = 0;
    int muon_count  = 0;

    // Counters for different types of neutrino interaction
    int ccqe_count = 0;
    int ccres_count = 0;
    int ccdis_count = 0;
    int cccoh_count = 0;
    int other2_count = 0;
    int nc_count     = 0;

    // Vectors to fill for the impurity implementation
    vector< double > T_mu_vect;
    vector< double > T_pi_vect;
    vector< double > cos_pi_vect;
    
    vector< int > Impurity;
   
    // Counters to check what is going on
    int nc_1pi_count      = 0;

    // Loop over the entries of the tree and calculate the kinetic energies 
    // of the muons and pions and define the impurity
    for ( int i = 0; i < n_values; ++i ){
        
        tree->GetEntry(i);
     
        double T_mu, e_mu;

        int nfpi0 = b_nfpi0->GetLeaf("nfpi0")->GetValue();
        int nfpip = b_nfpip->GetLeaf("nfpip")->GetValue();
        int nfpim = b_nfpim->GetLeaf("nfpim")->GetValue();
        int fspl = b_fspl->GetLeaf("fspl")->GetValue();
        int nf   = b_nf->GetLeaf("nf")->GetValue();

        // Calculate the kinetic energy for muons
        if ( fspl == 13 ){  
         
            // Energy of the final state primary lepton
            e_mu = b_El->GetLeaf("El")->GetValue();
            T_mu = e_mu - m_mu;

            T_mu_vect.push_back(T_mu);

            muon_count++;
        }
        // If the final state primary is a lepton, push back a number that will
        // be removed in the cuts later
        else if ( fspl == 11 ){
            T_mu_vect.push_back(-99999);
            e_count++;
        }
        else{
            T_mu_vect.push_back(-99999);
            other_count++;
        }

        // Calculate the kinetic energy of the pions
        if ( nfpip + nfpim == 1 && nfpi0 == 0){ // Get rid of this for cc inclusive

          //For all the final state hadronic particles, get their pdg code
          for (int j = 0; j<nf; ++j) {

            b_pdgf->GetEntry(i);
            b_cthf->GetEntry(i);
            b_Ef->GetEntry(i);

            int pdgf = b_pdgf->GetLeaf("pdgf")->GetValue(j);
            double e_pi = b_Ef->GetLeaf("Ef")->GetValue(j);
            double cos_pi = b_cthf->GetLeaf("cthf")->GetValue(j);

            if (pdgf == 211 || pdgf == -211){// Could have more than 1 pion per event, so should this be a vector?

              double T_pi = e_pi - m_pi;

              T_pi_vect.push_back(T_pi);
              cos_pi_vect.push_back(cos_pi);

              pion_count++;   
            }
          }
        }
        else{
            T_pi_vect.push_back(-99999);
            cos_pi_vect.push_back(-99999);
        }

        // Apply the impurity cut
        if ( b_nc->GetLeaf("nc")->GetValue() != 0 && nfpip + nfpim == 1 && nfpi0 == 0){// get rid of pion requirement for cc inc 
            int random;
            random = rand() % 5 + 1;
            nc_1pi_count++;

            // If the random number = 5 then push back a 1 onto the vector
            // If the random number < 5 push back a 0
            if( random == 5 ){
                Impurity.push_back(1);       
            }
            else{
                Impurity.push_back(0);
            }
        }
        else{
            Impurity.push_back(0);
        }
    }
    std::cout<<"Number of muons = "<<muon_count<<endl<<"Number of electrons = "<<e_count<<endl<<"Number of pions = "<<pion_count<<endl<<"Number of others = "<<other_count<<endl<<"Number of NC1pi = "<<nc_1pi_count<<endl;
    // Vectors to hold the smeared values of T and cos theta
    vector< double > T_mu_prime;
    vector< double > cos_mu_prime;

    // Initiate the random number generation
    ROOT::Math::GSLRngMT *_random_gen = new ROOT::Math::GSLRngMT;
    _random_gen->Initialize();    
    _random_gen->SetSeed( time( NULL ) ); 
    
    // Event by event, generate Tmu_prime and Tpi_prime: lognormal
    // Then find thetamu_prime and thetapi_prime: gaussian
    for ( int i = 0; i < n_values; ++i ){

        tree->GetEntry(i);

        double El    = b_El->GetLeaf( "El" )->GetValue();
        double cthl  = b_cthl->GetLeaf( "cthl" )->GetValue();

        // -------------------------------------------------------
        //                      Kinetic energy
        // -------------------------------------------------------
        // Calculate the mean and sigma for the LogNormal function
        //      zeta  = TMath::Log( m * ( 1 / sqrt( 1 + ( var / pow( m, 2 ) ) ) ) );
        //      sigma = sqrt( log( 1 + ( var / pow( m, 2 ) ) ) );
        //      m     = expectation value = El - m_mu
        //      var   = variance = s.d.^2 = ( El - m_mu * 0.1 ) ^ 2
    
        double var_mu     = TMath::Power( ( El - m_mu ) * 0.1, 2 ); 
        double sigma_mu   = TMath::Sqrt( TMath::Log( 1 + ( var_mu / TMath::Power( ( El - m_mu ), 2 ) ) ) );
        double zeta_mu    = TMath::Log( ( El - m_mu ) * ( 1 / TMath::Sqrt( 1 + ( var_mu / TMath::Power( ( El - m_mu ), 2 ) ) ) ) );
        double lognorm_mu = _random_gen->LogNormal( zeta_mu, sigma_mu );

        T_mu_prime.push_back(lognorm_mu);
            
        // -------------------------------------------------------
        //                      Cos theta
        // -------------------------------------------------------
        
        // Calculate the mean and sigma for the LogNormal function
        //      theta = acos(costtheta)
        //      var   = 5 degrees
    
        double sd_thetamu = TMath::Pi() / 36; // 5 degrees
        double gaus_theta = TMath::ACos( cthl ) + _random_gen->Gaussian( sd_thetamu );
        double gaus_costheta = TMath::Cos( gaus_theta ); 

        cos_mu_prime.push_back(gaus_costheta);

    } 

    // Implement the smearing factors and plot the histograms both before and after to observe the effect
    // Let's start with Gaussian smearing for now and move onto log normal smearing later
    // Now try the smearing
    TCanvas *c2 = new TCanvas("c2","",1000,800);

    int ccinc_count = 0;
    int imp_count = 0;

    // Fill h_smeared normally and with the impurities and THEN clone it
    for ( int i = 0; i < n_values; ++i ){        
        tree->GetEntry(i);

        double fspl  = b_fspl->GetLeaf( "fspl" )->GetValue(); 
        double pdgf  = b_pdgf->GetLeaf( "pdgf" )->GetValue(); 
        double cc    = b_cc->GetLeaf( "cc" )->GetValue(); 
        double nc    = b_nc->GetLeaf( "nc" )->GetValue(); 
        double nfpip = b_nfpip->GetLeaf( "nfpip" )->GetValue();
        double nfpim = b_nfpim->GetLeaf( "nfpim" )->GetValue();
        double nfpi0 = b_nfpi0->GetLeaf( "nfpi0" )->GetValue();
        double Ef    = b_Ef->GetLeaf( "Ef" )->GetValue();
        double pf    = b_pf->GetLeaf( "pf" )->GetValue();
        double cthf  = b_cthf->GetLeaf( "cthf" )->GetValue();
        double El    = b_El->GetLeaf( "El" )->GetValue();
        double pl    = b_pl->GetLeaf( "pl" )->GetValue();
        double cthl  = b_cthl->GetLeaf( "cthl" )->GetValue();
        double qel   = b_qel->GetLeaf( "qel" )->GetValue();
        double res   = b_res->GetLeaf( "res" )->GetValue();
        double dis   = b_dis->GetLeaf( "dis" )->GetValue();
        double coh   = b_coh->GetLeaf( "coh" )->GetValue();

        //Count the event types
        if ( fspl == 13 && cc == 1 ) {
          if ( qel == 1 ) { h_ccqe_un->Fill(cthl,T_mu_vect[i]); ccqe_count++; }
          else if ( res == 1 ) { h_ccres_un->Fill(cthl,T_mu_vect[i]); ccres_count++; }
          else if ( dis == 1 ) { h_ccdis_un->Fill(cthl,T_mu_vect[i]); ccdis_count++; }
          else if ( coh == 1 ) { h_cccoh_un->Fill(cthl,T_mu_vect[i]); cccoh_count++; }
          else { other2_count++;}
        }
        // Set the energy and impurity cuts
        if ( fspl == 13 
            && cc != 0  
            && T_mu_vect[i] > 0.05 ){

            if ( qel == 1 ) { h_ccqe_sm->Fill(cos_mu_prime[i],T_mu_prime[i]); }
            else if ( res == 1 ) { h_ccres_sm->Fill(cos_mu_prime[i],T_mu_prime[i]); }
            else if ( dis == 1 ) { h_ccdis_sm->Fill(cos_mu_prime[i],T_mu_prime[i]); }
            else if ( coh == 1 ) { h_cccoh_sm->Fill(cos_mu_prime[i],T_mu_prime[i]); }

            h_Tmu_un->Fill(T_mu_vect[i]);
            h_Tmu_sm->Fill(T_mu_prime[i]);
            h_costhetamu_un->Fill(cthl);
            h_costhetamu_sm->Fill(cos_mu_prime[i]);
            h_smeared->Fill(cos_mu_prime[i], T_mu_prime[i]);
            ccinc_count++;
        }
        else if ( Impurity[i] == 1 
            && T_pi_vect[i] > 0.05 ){

            h_nc_sm->Fill(cos_pi_vect[i], T_pi_vect[i]);
            h_smeared->Fill(cos_pi_vect[i], T_pi_vect[i]);
            imp_count++;
            nc_count++;
        }
    }
    std::cout<<"Number of true CC events = "<<ccinc_count<<endl<<"Number of impurity events = "<<imp_count<<endl;
    std::cout<<endl<<"Signal event types:"<<endl<<"Number of CCQE = "<<ccqe_count<<endl<<"Number of CCRES = "<<ccres_count<<endl<<"Number of CCDIS = "<<ccdis_count<<endl<<"Number of CCCOH = "<<cccoh_count<<endl<<"Number of others = "<<other2_count<<endl;

    h_cccoh_un->SetFillColor(kMagenta+2);
    h_ccdis_un->SetFillColor(kBlue-4);
    h_ccres_un->SetFillColor(kGreen-4);
    h_ccqe_un->SetFillColor(kRed-4);
    TLegend *leg = new TLegend(0.2,0.8,0.4,0.4);
    leg->AddEntry(h_cccoh_un,"CCCOH","f");
    leg->AddEntry(h_ccdis_un,"CCDIS","f");
    leg->AddEntry(h_ccres_un,"CCRES","f");
    leg->AddEntry(h_ccqe_un,"CCQE","f");
    hs_un->Add(h_cccoh_un);
    hs_un->Add(h_ccdis_un);
    hs_un->Add(h_ccres_un);
    hs_un->Add(h_ccqe_un);

    v_un.push_back(h_cccoh_un);
    v_un.push_back(h_ccdis_un);
    v_un.push_back(h_ccres_un);
    v_un.push_back(h_ccqe_un);

    h_cccoh_sm->SetFillColor(kMagenta+2);
    h_ccdis_sm->SetFillColor(kBlue-4);
    h_ccres_sm->SetFillColor(kGreen-4);
    h_ccqe_sm->SetFillColor(kRed-4);
    h_nc_sm->SetFillColor(kOrange);
    TLegend *leg2 = new TLegend(0.2,0.8,0.4,0.4);
    leg2->AddEntry(h_cccoh_sm,"CCCOH","f");
    leg2->AddEntry(h_ccdis_sm,"CCDIS","f");
    leg2->AddEntry(h_ccres_sm,"CCRES","f");
    leg2->AddEntry(h_ccqe_sm,"CCQE","f");
    leg2->AddEntry(h_nc_sm,"NC1pi","f");
    hs_sm->Add(h_nc_sm);
    hs_sm->Add(h_cccoh_sm);
    hs_sm->Add(h_ccdis_sm);
    hs_sm->Add(h_ccres_sm);
    hs_sm->Add(h_ccqe_sm);

    v_sm.push_back(h_nc_sm);
    v_sm.push_back(h_cccoh_sm);
    v_sm.push_back(h_ccdis_sm);
    v_sm.push_back(h_ccres_sm);
    v_sm.push_back(h_ccqe_sm);

    double nX = h_smeared->GetNbinsX();
    double nY = h_smeared->GetNbinsY();
 
    // Write the bin contents to .tex files
    for ( int i = nY; i >= 1; --i){
        for ( int j = 1; j <= nX - 1; ++j ){
            
            file1 << h_unsmeared->GetBinContent(j,i) << " & ";
            file2 << h_smeared->GetBinContent(j,i) << " & "; 
            file3 << ( h_smeared->GetBinContent(j,i) - h_unsmeared->GetBinContent(j,i) ) << " & ";
        }

        file1 << h_unsmeared->GetBinContent(nX, i) << " \\\\ " << endl;
        file2 << h_smeared->GetBinContent(nX, i) << " \\\\ " << endl;
        file3 << ( h_smeared->GetBinContent(nX, i) - h_unsmeared->GetBinContent(nX, i) ) << " \\\\ " << endl;
        
    }
    file1 << " \\end{tabular} " << endl;
    file2 << " \\end{tabular} " << endl;
    file3 << " \\end{tabular} " << endl;

    f->cd();

    c2->SetRightMargin(0.13);
    
    h_smeared->Draw("colz");
    h_smeared->SetStats(kFALSE);
    h_smeared->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_smeared->GetYaxis()->SetTitle("T_{#mu}");
    //c2->SetLogz();
    c2->SaveAs("distribution_after_impure.png");
    //c2->SaveAs("distribution_after_log_impure.png");

    h_Tmu_un->Draw();
    h_Tmu_un->GetXaxis()->SetTitle("T_{#mu}");
    h_Tmu_un->Write("Tmu_before_smearing");

    h_Tmu_sm->Draw();
    h_Tmu_sm->GetXaxis()->SetTitle("T_{#mu}");
    h_Tmu_sm->Write("Tmu_after_smearing");

    h_costhetamu_un->Draw();
    h_costhetamu_un->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_costhetamu_un->Write("costhetamu_before_smearing");

    h_costhetamu_sm->Draw();
    h_costhetamu_sm->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_costhetamu_sm->Write("costhetamu_after_smearing");

    hs_un->Draw();
    leg->Draw("same");
    hs_un->Write("Stacked_hist_before_smearing");

    hs_sm->Draw();
    leg2->Draw("same");
    hs_sm->Write("Stacked_hist_after_smearing");

    THStack *hs_Tmu_un = new THStack("hs_Tmu_un","");
    hs_Tmu_un->Add(h_cccoh_un->ProjectionY());
    hs_Tmu_un->Add(h_ccdis_un->ProjectionY());
    hs_Tmu_un->Add(h_ccres_un->ProjectionY());
    hs_Tmu_un->Add(h_ccqe_un->ProjectionY());
    hs_Tmu_un->Write("Tmu unsmeared Stack");

    THStack *hs_ctmu_un = new THStack("hs_ctmu_un","");
    hs_ctmu_un->Add(h_cccoh_un->ProjectionX());
    hs_ctmu_un->Add(h_ccdis_un->ProjectionX());
    hs_ctmu_un->Add(h_ccres_un->ProjectionX());
    hs_ctmu_un->Add(h_ccqe_un->ProjectionX());
    hs_ctmu_un->Write("costhetamu unsmeared Stack");

    THStack *hs_Tmu_sm = new THStack("hs_Tmu_sm","");
    hs_Tmu_sm->Add(h_nc_sm->ProjectionY());
    hs_Tmu_sm->Add(h_cccoh_sm->ProjectionY());
    hs_Tmu_sm->Add(h_ccdis_sm->ProjectionY());
    hs_Tmu_sm->Add(h_ccres_sm->ProjectionY());
    hs_Tmu_sm->Add(h_ccqe_sm->ProjectionY());
    hs_Tmu_sm->Write("Tmu smeared Stack");

    THStack *hs_ctmu_sm = new THStack("hs_ctmu_sm","");
    hs_ctmu_sm->Add(h_nc_sm->ProjectionX());
    hs_ctmu_sm->Add(h_cccoh_sm->ProjectionX());
    hs_ctmu_sm->Add(h_ccdis_sm->ProjectionX());
    hs_ctmu_sm->Add(h_ccres_sm->ProjectionX());
    hs_ctmu_sm->Add(h_ccqe_sm->ProjectionX());
    hs_ctmu_sm->Write("costhetamu smeared Stack");


}

void Slices ( TH2D *h_unsmeared, TH2D *h_smeared ){

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = h_unsmeared->GetNbinsX(); // Cos theta
    int y_bins = h_unsmeared->GetNbinsY(); // Tmu

    TCanvas *c_Tmu   = new TCanvas ( "c_Tmu", "", 800, 600 );
   
    TLegend *leg_T = new TLegend( 0.12, 0.78, 0.28, 0.88 );

    TH1D *h_Tmu      = new TH1D ( "h_Tmu", "", x_bins, -1, 1 );
    TH1D *h_Tmu_sm   = new TH1D ( "h_Tmu_sm", "", x_bins, -1, 1 );
    
    leg_T->AddEntry( h_Tmu, " Unsmeared ", "l" );
    leg_T->AddEntry( h_Tmu_sm, " Smeared ", "l" );
    
    for ( int i = 1; i <= y_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_T;
        double upp_edge_T;

        low_edge_T = h_unsmeared->GetYaxis()->GetBinLowEdge(i);
        upp_edge_T = h_unsmeared->GetYaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream conv;
        conv.clear();

        string title;
        title.clear();

        char file_name[1024];

        conv << setprecision(4) << "Tmu_slice_" << low_edge_T << ";" << upp_edge_T << ".root";
        title = conv.str();
        
        strcpy( file_name, title.c_str() );

        // For the title of the histogram
        stringstream hist;
        hist.clear();

        char hist_name[1024];
        string temp1;

        hist << setprecision(4) << "T_{#mu} slice: " << low_edge_T << "-" << upp_edge_T;
        temp1 = hist.str();

        strcpy( hist_name, temp1.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= x_bins; ++j ){
            h_Tmu->SetBinContent( j, h_unsmeared->GetBinContent(j, i) );
            h_Tmu_sm->SetBinContent( j, h_smeared->GetBinContent(j, i) );
        }

        h_Tmu->Draw();
        h_Tmu_sm->Draw("same");
        h_Tmu->SetTitle(hist_name);
        h_Tmu->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        h_Tmu->GetYaxis()->SetTitle("Number of SBND events");   
        h_Tmu->SetLineColor( kRed + 2 );
        h_Tmu_sm->SetLineColor( kGreen + 2 );
        h_Tmu->SetTitleOffset(1.5, "Y");
        h_Tmu->SetStats(kFALSE);

        leg_T->Draw();
        c_Tmu->Print(file_name);

    } 
   
    delete h_Tmu;
    delete h_Tmu_sm;

    delete c_Tmu;
    
    delete leg_T;
    
    TCanvas *c_cosmu = new TCanvas ( "c_cosmu", "", 800, 600 );
    
    TLegend *leg_c = new TLegend( 0.72, 0.78, 0.88, 0.88 );

    TH1D *h_cosmu    = new TH1D ( "h_cosmu", "", y_bins, 0, 2 );
    TH1D *h_cosmu_sm = new TH1D ( "h_cosmu_sm", "", y_bins, 0, 2 );
     
    leg_c->AddEntry( h_cosmu, " Unsmeared ", "l" );
    leg_c->AddEntry( h_cosmu_sm, " Smeared ", "l" );
    
    // Cos theta mu slices
    for ( int i = 1; i <= x_bins; ++i ){

        // Define the histograms
        // Get the lower and upper bin edges of the y axis
        double low_edge_cos;
        double upp_edge_cos;

        low_edge_cos = h_unsmeared->GetXaxis()->GetBinLowEdge(i);
        upp_edge_cos = h_unsmeared->GetXaxis()->GetBinLowEdge(i + 1);
   
        // Make the title for the current histogram and the file name
        // Clear the title for the new loop
        stringstream conv1;
        conv1.clear();

        string title1;
        title1.clear();

        char file_name1[1024];

        conv1 << setprecision(4) << "cos_thetamu_slice_" << low_edge_cos << ";" << upp_edge_cos << ".root";
        title1 = conv1.str();
        
        strcpy( file_name1, title1.c_str() );

        // For the title of the histogram
        stringstream hist1;
        hist1.clear();

        char hist_name1[1024];
        string temp2;

        hist1 << setprecision(4) << "cos#theta_{#mu} slice: " << low_edge_cos << "," << upp_edge_cos;
        temp2 = hist1.str();

        strcpy( hist_name1, temp2.c_str() );

        // Fill the histogram
        for ( int j = 1; j <= y_bins; ++j ){
            h_cosmu->SetBinContent( j, h_unsmeared->GetBinContent(i, j) );
            h_cosmu_sm->SetBinContent( j, h_smeared->GetBinContent(i, j) );
        }

        h_cosmu->Draw();
        h_cosmu_sm->Draw("same");
        h_cosmu->SetTitle(hist_name1);
        h_cosmu->GetXaxis()->SetTitle("T_{#mu}");   
        h_cosmu->GetYaxis()->SetTitle("Number of SBND events");   
        h_cosmu->SetLineColor( kRed + 2 );
        h_cosmu_sm->SetLineColor( kGreen + 2 );
        h_cosmu->SetTitleOffset(1.5, "Y");
        h_cosmu->SetStats(kFALSE);

        leg_c->Draw();
        c_cosmu->SaveAs(file_name1);

    } 
    
    delete h_cosmu;
    delete h_cosmu_sm;

    delete c_cosmu;

    delete leg_c;

}

void SliceStack ( std::vector<TH2*> v_un, std::vector<TH2*> v_sm, TFile *f ){

    f->cd();

    for ( size_t j = 1; j < 21; j++ ){
   
      // Make the title for the current histogram and the file name
      // Clear the title for the new loop
      stringstream conv;
      conv.clear();

      string title;
      title.clear();

      conv << setprecision(4) << "Tmu_un_cosmu:" << v_un[0]->GetXaxis()->GetBinLowEdge(j) << ";" << v_un[0]->GetXaxis()->GetBinLowEdge(j+1);;
      title = conv.str();

      // Make the title for the current histogram and the file name
      // Clear the title for the new loop
      stringstream conv2;
      conv2.clear();

      string title2;
      title2.clear();

      conv2 << setprecision(4) << "Tmu_sm_cosmu:" << v_sm[0]->GetXaxis()->GetBinLowEdge(j) << ";" << v_sm[0]->GetXaxis()->GetBinLowEdge(j+1);;
      title2 = conv2.str();

      THStack *hs_Tmu_un_tmp = new THStack("hs_tmu_un_tmp","");
      for ( size_t i = 0; i < v_un.size(); i++){
        //Fill THStack with projection Y of v_un for the jth bin
        hs_Tmu_un_tmp->Add(v_un[i]->ProjectionY(std::to_string(i).c_str(),j,j));
      }
      hs_Tmu_un_tmp->Write(title.c_str());
      delete hs_Tmu_un_tmp;

      THStack *hs_Tmu_sm_tmp = new THStack("hs_tmu_sm_tmp","");
      for ( size_t i = 0; i < v_sm.size(); i++){
        //Fill THStack with projection Y of v_un for the jth bin
        hs_Tmu_sm_tmp->Add(v_sm[i]->ProjectionY(std::to_string(i).c_str(),j,j));
      }
      hs_Tmu_sm_tmp->Write(title2.c_str());
      delete hs_Tmu_sm_tmp;

    }

   for ( size_t j = 1; j < 19; j++ ){
   
      // Make the title for the current histogram and the file name
      // Clear the title for the new loop
      stringstream conv;
      conv.clear();

      string title;
      title.clear();

      conv << setprecision(4) << "cosmu_un_Tmu:" << v_un[0]->GetYaxis()->GetBinLowEdge(j) << ";" << v_un[0]->GetYaxis()->GetBinLowEdge(j+1);;
      title = conv.str();

      // Make the title for the current histogram and the file name
      // Clear the title for the new loop
      stringstream conv2;
      conv2.clear();

      string title2;
      title2.clear();

      conv2 << setprecision(4) << "cosmu_sm_Tmu:" << v_sm[0]->GetYaxis()->GetBinLowEdge(j) << ";" << v_sm[0]->GetYaxis()->GetBinLowEdge(j+1);;
      title2 = conv2.str();

      THStack *hs_cosmu_un_tmp = new THStack("hs_cosmu_un_tmp","");
      for ( size_t i = 0; i < v_un.size(); i++){
        //Fill THStack with projection Y of v_un for the jth bin
        hs_cosmu_un_tmp->Add(v_un[i]->ProjectionX(std::to_string(i).c_str(),j,j));
      }
      hs_cosmu_un_tmp->Write(title.c_str());
      delete hs_cosmu_un_tmp;

      THStack *hs_cosmu_sm_tmp = new THStack("hs_cosmu_sm_tmp","");
      for ( size_t i = 0; i < v_sm.size(); i++){
        //Fill THStack with projection Y of v_un for the jth bin
        hs_cosmu_sm_tmp->Add(v_sm[i]->ProjectionX(std::to_string(i).c_str(),j,j));
      }
      hs_cosmu_sm_tmp->Write(title2.c_str());
      delete hs_cosmu_sm_tmp;

    }

    
}
