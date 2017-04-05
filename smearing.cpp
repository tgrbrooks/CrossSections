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

double SmearKE(double ke, ROOT::Math::GSLRngMT *_random_gen){

  // -------------------------------------------------------
  //                      Kinetic energy
  // -------------------------------------------------------
  // Calculate the mean and sigma for the LogNormal function
  //      zeta  = TMath::Log( m * ( 1 / sqrt( 1 + ( var / pow( m, 2 ) ) ) ) );
  //      sigma = sqrt( log( 1 + ( var / pow( m, 2 ) ) ) );
  //      m     = expectation value = El - m_mu
  //      var   = variance = s.d.^2 = ( El - m_mu * 0.1 ) ^ 2

  double var     = TMath::Power( ( ke ) * 0.1, 2 ); 
  double sigma   = TMath::Sqrt( TMath::Log( 1 + ( var / TMath::Power( ( ke ), 2 ) ) ) );
  double zeta    = TMath::Log( ( ke ) * ( 1 / TMath::Sqrt( 1 + ( var / TMath::Power( ( ke ), 2 ) ) ) ) );
  double lognorm = _random_gen->LogNormal( zeta, sigma );

  return lognorm;
}

double SmearCosTheta(double costheta, ROOT::Math::GSLRngMT *_random_gen){
    
  // -------------------------------------------------------
  //                      Cos theta
  // -------------------------------------------------------

  // Calculate the mean and sigma for the LogNormal function
  //      theta = acos(costtheta)
  //      var   = 5 degrees

  double sd_theta = TMath::Pi() / 36; // 5 degrees
  double gaus_theta = TMath::ACos( costheta ) + _random_gen->Gaussian( sd_theta );
  double gaus_costheta = TMath::Cos( gaus_theta ); 

  return gaus_costheta;
}

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

    h_un->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un->GetYaxis()->SetTitle("T_{#mu}");
    h_un->Write("distribution_before");

    std::vector<TH2*> v_un;
    std::vector<TH2*> v_sm;
    
    // The same histogram definitions for the smeared distributions
    TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} distribution after smearing, impurity ",20,-1,1,18,0,2);
    
    // Take h_un and smear it
    Smear(def_tree, h_sm, v_un, v_sm, f);

    // The the 2 2D histograms and draw slices in Tmu and cos theta mu
    Slices( h_un, h_sm );

    SliceStack( v_un, v_sm, f );

    TH2D *h_eff = new TH2D("h_eff","h_eff;cos#theta_{#mu};T_{#mu}",20,-1,1,18,0,2);
    TH2D *h_pur = new TH2D("h_pur","h_pur;cos#theta_{#mu};T_{#mu}",20,-1,1,18,0,2);
    TH2D *h_effpur = new TH2D("h_effpur","h_effpur;cos#theta_{#mu};T_{#mu}",20,-1,1,18,0,2);
    TH2D *h_impur = new TH2D("h_impur","h_impur;cos#theta_{#mu};T_{#mu}",20,-1,1,18,0,2);

    //Fill efficiency histogram with Nsig/Nmcsig
    h_eff->Add(h_sm,v_sm[0],1.,-1.);
    h_eff->Divide(h_un);
    //Fill purity histogram with Nsig/(Nsig+Nbkg)
    h_pur->Add(h_sm,v_sm[0],1.,-1.);
    h_pur->Divide(h_sm);
    //Fill efficiency*purity histogram with h_eff*h_pur
    h_effpur->Multiply(h_eff,h_pur,1.,1.);
    //Fill impurity histogram with Nbkg/(Nsig+Nbkg)
    h_impur->Divide(v_sm[0],h_sm,1.,1.);

    h_eff->SetMinimum(0.5);
    h_eff->SetMaximum(1.0);
    h_eff->Write("h_eff");
    h_pur->SetMinimum(0.7);
    h_pur->SetMaximum(1.0);
    h_pur->Write("h_pur");
    h_effpur->SetMinimum(0.5);
    h_effpur->SetMaximum(1.0);
    h_effpur->Write("h_effpur");
    h_impur->Write("h_impur");

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
             TH2D  *h_smeared, 
             std::vector<TH2*> &v_un,
             std::vector<TH2*> &v_sm,
             TFile *f ){

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
    TH2D *h_nc_sm = new TH2D("h_nc_sm","NCnpi",20,-1,1,18,0,2);

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

    // Counters for different types of neutrino interaction
    int ccqe_count_un = 0; int ccqe_count_sm = 0;
    int ccres_count_un = 0; int ccres_count_sm = 0;
    int ccdis_count_un = 0; int ccdis_count_sm = 0;
    int cccoh_count_un = 0; int cccoh_count_sm = 0;
    int other_count_un = 0; int other_count_sm = 0;

    // Total counters
    int ccinc_count_un = 0;
    int ccinc_count_sm = 0;
    int imp_count = 0;
    int nc_npi_count = 0;

    // Initiate the random number generation
    ROOT::Math::GSLRngMT *_random_gen = new ROOT::Math::GSLRngMT;
    _random_gen->Initialize();    
    _random_gen->SetSeed( time( NULL ) ); 

    // Loop over the entries of the tree and calculate the kinetic energies 
    // of the muons and pions and define the impurity
    for ( int i = 0; i < n_values; ++i ){
        
        tree->GetEntry(i);
     
        double T_mu, e_mu;

        // Get values from the branches
        int nf   = b_nf->GetLeaf("nf")->GetValue();
        double fspl  = b_fspl->GetLeaf( "fspl" )->GetValue(); 
        double cc    = b_cc->GetLeaf( "cc" )->GetValue(); 
        double nc    = b_nc->GetLeaf( "nc" )->GetValue(); 
        double nfpip = b_nfpip->GetLeaf( "nfpip" )->GetValue();
        double nfpim = b_nfpim->GetLeaf( "nfpim" )->GetValue();
        double nfpi0 = b_nfpi0->GetLeaf( "nfpi0" )->GetValue();
        double El    = b_El->GetLeaf( "El" )->GetValue();
        double pl    = b_pl->GetLeaf( "pl" )->GetValue();
        double cthl  = b_cthl->GetLeaf( "cthl" )->GetValue();
        double qel   = b_qel->GetLeaf( "qel" )->GetValue();
        double res   = b_res->GetLeaf( "res" )->GetValue();
        double dis   = b_dis->GetLeaf( "dis" )->GetValue();
        double coh   = b_coh->GetLeaf( "coh" )->GetValue();

        // Kinetic energy of the final state primary lepton
        T_mu = El - m_mu;

        //Count the CC numu event types before smearing
        if ( fspl == 13 && cc == 1 ) {

          if ( qel == 1 ) { h_ccqe_un->Fill(cthl,T_mu); ccqe_count_un++; }
          else if ( res == 1 ) { h_ccres_un->Fill(cthl,T_mu); ccres_count_un++; }
          else if ( dis == 1 ) { h_ccdis_un->Fill(cthl,T_mu); ccdis_count_un++; }
          else if ( coh == 1 ) { h_cccoh_un->Fill(cthl,T_mu); cccoh_count_un++; }
          else { other_count_un++; }

          ccinc_count_un++;
        }

        //Count the CC numu event types after smearing
        if ( fspl == 13 && cc == 1  && T_mu > 0.05 ){

          //Smear the muon energy and angle
          if ( qel == 1 ) { h_ccqe_sm->Fill(SmearCosTheta(cthl,_random_gen),SmearKE(T_mu,_random_gen)); ccqe_count_sm++; }
          else if ( res == 1 ) { h_ccres_sm->Fill(SmearCosTheta(cthl,_random_gen),SmearKE(T_mu,_random_gen)); ccres_count_sm++; }
          else if ( dis == 1 ) { h_ccdis_sm->Fill(SmearCosTheta(cthl,_random_gen),SmearKE(T_mu,_random_gen)); ccdis_count_sm++; }
          else if ( coh == 1 ) { h_cccoh_sm->Fill(SmearCosTheta(cthl,_random_gen),SmearKE(T_mu,_random_gen)); cccoh_count_sm++; }
          else { other_count_sm++; }

          h_smeared->Fill(SmearCosTheta(cthl,_random_gen),SmearKE(T_mu,_random_gen));
          ccinc_count_sm++;

        }

        // Find NCnpi events
        if ( nc == 1 ){ 
          // Count total number
          if ( nfpip + nfpip > 0 ) nc_npi_count++;

          double prev_e_pi = 0;
          double prev_cos_pi = 0;

          //Loop over all the final state hadronic particles
          for (int j = 0; j<nf; ++j) {

            b_pdgf->GetEntry(i);
            b_cthf->GetEntry(i);
            b_Ef->GetEntry(i);

            int pdgf = b_pdgf->GetLeaf("pdgf")->GetValue(j);
            double e_pi = b_Ef->GetLeaf("Ef")->GetValue(j);
            double cos_pi = b_cthf->GetLeaf("cthf")->GetValue(j);

            // If one of the hadrons is a pion use a random number to check if it is misidentified
            if (pdgf == 211 || pdgf == -211){
              int random;
              random = rand() % 5 + 1;
              // If more than one is misidentified assume the higher energy one is muon
              if ( random == 5 && e_pi > prev_e_pi ) {
                prev_e_pi = e_pi;
                prev_cos_pi = cos_pi;
              }
            }
          }
          //Kinetic energy
          double T_pi = prev_e_pi - m_pi;
          // If at least one pion is misidentified and it has a high enough energy to be detected record it as a muon
          if( prev_e_pi != 0 && T_pi > 0.05){
            //Smear the pion energy and angle
            h_nc_sm->Fill(SmearCosTheta(prev_cos_pi,_random_gen), SmearKE(T_pi,_random_gen));
            h_smeared->Fill(SmearCosTheta(prev_cos_pi,_random_gen), SmearKE(T_pi,_random_gen));
            imp_count++;
          }
        }
    }
     
    //Output the total numbers of events before and after smearing
    cout<<" ---------------------------------------------------------- "<<endl
        <<"| Before Smearing        "          <<" | After Smearing          "          <<"       |"<<endl
        <<"|----------------------------------------------------------|"<<endl
        <<"| Total CCInc =    "<<setfill(' ')<<setw(6)<<ccinc_count_un<<" | Total CCInc =           "<<ccinc_count_sm<<" |"<<endl
        <<"| Total NCnpi+/- = "<<setfill(' ')<<setw(6)<<nc_npi_count  <<" | Total Impurity events = "<<setfill(' ')<<setw(6)<<imp_count     <<" |"<<endl
        <<"|----------------------------------------------------------|"<<endl
        <<"|----------------- Signal Event Types: --------------------|"<<endl
        <<"|----------------------------------------------------------|"<<endl
        <<"| CCQE =           "<<setfill(' ')<<setw(6)<<ccqe_count_un <<" | CCQE =                  "<<setfill(' ')<<setw(6)<<ccqe_count_sm <<" |"<<endl
        <<"| CCRES =          "<<setfill(' ')<<setw(6)<<ccres_count_un<<" | CCRES =                 "<<setfill(' ')<<setw(6)<<ccres_count_sm<<" |"<<endl
        <<"| CCDIS =          "<<setfill(' ')<<setw(6)<<ccdis_count_un<<" | CCDIS =                 "<<setfill(' ')<<setw(6)<<ccdis_count_sm<<" |"<<endl
        <<"| CCCOH =          "<<setfill(' ')<<setw(6)<<cccoh_count_un<<" | CCCOH =                 "<<setfill(' ')<<setw(6)<<cccoh_count_sm<<" |"<<endl
        <<"| Others =         "<<setfill(' ')<<setw(6)<<other_count_un<<" | Others =                "<<setfill(' ')<<setw(6)<<other_count_sm<<" |"<<endl
        <<" ---------------------------------------------------------- "<<endl;

    //Store the unsmeared component histograms
    v_un.push_back(h_cccoh_un);
    v_un.push_back(h_ccdis_un);
    v_un.push_back(h_ccres_un);
    v_un.push_back(h_ccqe_un);

    //Store the smeared component histograms
    v_sm.push_back(h_nc_sm);
    v_sm.push_back(h_cccoh_sm);
    v_sm.push_back(h_ccdis_sm);
    v_sm.push_back(h_ccres_sm);
    v_sm.push_back(h_ccqe_sm);

    f->cd();
    
    //Format the smeared total 2D histogram and write to file
    h_smeared->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_smeared->GetYaxis()->SetTitle("T_{#mu}");
    h_smeared->Write("distribution_after_impure");

    //Format the unsmeared stacked histogram and write to file
    h_cccoh_un->SetFillColor(kMagenta+2);
    h_ccdis_un->SetFillColor(kBlue-4);
    h_ccres_un->SetFillColor(kGreen-4);
    h_ccqe_un->SetFillColor(kRed-4);
    hs_un->Add(h_cccoh_un);
    hs_un->Add(h_ccdis_un);
    hs_un->Add(h_ccres_un);
    hs_un->Add(h_ccqe_un);
    hs_un->Write("Stacked_hist_before_smearing");

    //Format the smeared stacked histograms and write to file
    h_cccoh_sm->SetFillColor(kMagenta+2);
    h_ccdis_sm->SetFillColor(kBlue-4);
    h_ccres_sm->SetFillColor(kGreen-4);
    h_ccqe_sm->SetFillColor(kRed-4);
    h_nc_sm->SetFillColor(kOrange);
    hs_sm->Add(h_nc_sm);
    hs_sm->Add(h_cccoh_sm);
    hs_sm->Add(h_ccdis_sm);
    hs_sm->Add(h_ccres_sm);
    hs_sm->Add(h_ccqe_sm);
    hs_sm->Write("Stacked_hist_after_smearing");

    //Create stacked histogram for just Tmu dependence before smearing
    THStack *hs_Tmu_un = new THStack("hs_Tmu_un","");
    hs_Tmu_un->Add(h_cccoh_un->ProjectionY());
    hs_Tmu_un->Add(h_ccdis_un->ProjectionY());
    hs_Tmu_un->Add(h_ccres_un->ProjectionY());
    hs_Tmu_un->Add(h_ccqe_un->ProjectionY());
    hs_Tmu_un->Write("Tmu unsmeared Stack");

    //Create stacked histogram for cos theta mu before smearing
    THStack *hs_ctmu_un = new THStack("hs_ctmu_un","");
    hs_ctmu_un->Add(h_cccoh_un->ProjectionX());
    hs_ctmu_un->Add(h_ccdis_un->ProjectionX());
    hs_ctmu_un->Add(h_ccres_un->ProjectionX());
    hs_ctmu_un->Add(h_ccqe_un->ProjectionX());
    hs_ctmu_un->Write("costhetamu unsmeared Stack");

    //Create stacked histogram for Tmu after smearing
    THStack *hs_Tmu_sm = new THStack("hs_Tmu_sm","");
    hs_Tmu_sm->Add(h_nc_sm->ProjectionY());
    hs_Tmu_sm->Add(h_cccoh_sm->ProjectionY());
    hs_Tmu_sm->Add(h_ccdis_sm->ProjectionY());
    hs_Tmu_sm->Add(h_ccres_sm->ProjectionY());
    hs_Tmu_sm->Add(h_ccqe_sm->ProjectionY());
    hs_Tmu_sm->Write("Tmu smeared Stack");

    //Create stacked histogram for cos theta mu after smearing
    THStack *hs_ctmu_sm = new THStack("hs_ctmu_sm","");
    hs_ctmu_sm->Add(h_nc_sm->ProjectionX());
    hs_ctmu_sm->Add(h_cccoh_sm->ProjectionX());
    hs_ctmu_sm->Add(h_ccdis_sm->ProjectionX());
    hs_ctmu_sm->Add(h_ccres_sm->ProjectionX());
    hs_ctmu_sm->Add(h_ccqe_sm->ProjectionX());
    hs_ctmu_sm->Write("costhetamu smeared Stack");

}

// -------------------------------------------------------------------------
//                             Slicing function
// -------------------------------------------------------------------------
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

// -------------------------------------------------------------------------
//                  Stacked histogram slicing function
// -------------------------------------------------------------------------
void SliceStack ( std::vector<TH2*> v_un, std::vector<TH2*> v_sm, TFile *f ){

    f->cd();

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = v_un[0]->GetNbinsX(); // Cos theta
    int y_bins = v_un[0]->GetNbinsY(); // Tmu

    for ( int j = 1; j < x_bins+1; j++ ){
   
      // Make the title for the current histogram and the file name
      // Clear the title for the new loop
      stringstream conv;
      conv.clear();

      string title;
      title.clear();

      conv << setprecision(4) << "un_Tmu_cosmuRange:" << v_un[0]->GetXaxis()->GetBinLowEdge(j) << ";" << v_un[0]->GetXaxis()->GetBinLowEdge(j+1);;
      title = conv.str();

      // Make the title for the current histogram and the file name
      // Clear the title for the new loop
      stringstream conv2;
      conv2.clear();

      string title2;
      title2.clear();

      conv2 << setprecision(4) << "sm_Tmu_cosmuRange:" << v_sm[0]->GetXaxis()->GetBinLowEdge(j) << ";" << v_sm[0]->GetXaxis()->GetBinLowEdge(j+1);;
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

   for ( int j = 1; j < y_bins+1; j++ ){
   
      // Make the title for the current histogram and the file name
      // Clear the title for the new loop
      stringstream conv;
      conv.clear();

      string title;
      title.clear();

      conv << setprecision(4) << "un_cosmu_TmuRange:" << v_un[0]->GetYaxis()->GetBinLowEdge(j) << ";" << v_un[0]->GetYaxis()->GetBinLowEdge(j+1);;
      title = conv.str();

      // Make the title for the current histogram and the file name
      // Clear the title for the new loop
      stringstream conv2;
      conv2.clear();

      string title2;
      title2.clear();

      conv2 << setprecision(4) << "sm_cosmu_TmuRange:" << v_sm[0]->GetYaxis()->GetBinLowEdge(j) << ";" << v_sm[0]->GetYaxis()->GetBinLowEdge(j+1);;
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
