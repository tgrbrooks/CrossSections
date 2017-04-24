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

  double var     = TMath::Power( ( ke ) * 0.1, 2 ); //WAS 0.1 
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

void StackProjectionY(std::vector<TH2*> histograms, string legendlabels[], string filename, int bin, TCanvas *canvas){
  canvas->Clear();
  // Make the title for the current histogram and the file name
  // Clear the title for the new loop
  stringstream conv;
  conv.clear();

  string title;
  title.clear();

  conv << setprecision(4) << "~/Documents/PhD/CrossSections/output/cosslicestack/" << filename 
       << histograms[0]->GetXaxis()->GetBinLowEdge(bin) << ";" << histograms[0]->GetXaxis()->GetBinLowEdge(bin+1) << ".root";
  title = conv.str();

  THStack *hs_Tmu_un_tmp = new THStack("hs_tmu_un_tmp","");
  TLegend *leg_Tmu_un = new TLegend( 0.584, 0.614, 0.901, 0.904 );
  for ( size_t i = 0; i < histograms.size(); i++){
    //Fill THStack with projection Y of v_un for the jth bin
    TH1D *temphist = histograms[i]->ProjectionY(std::to_string(i).c_str(),bin,bin);
    temphist->SetLineWidth(1);
    temphist->SetLineColor(temphist->GetFillColor());
    temphist->SetFillStyle(1001);
    if(legendlabels[i]=="NC 1#pi"||legendlabels[i]=="NC Other"||legendlabels[i]=="NC n#pi"){
      temphist->SetFillStyle(3004);
    }
    hs_Tmu_un_tmp->Add(temphist);
    leg_Tmu_un->AddEntry(temphist,legendlabels[i].c_str(),"f");
  }
  hs_Tmu_un_tmp->Draw();
  hs_Tmu_un_tmp->GetXaxis()->SetTitle("T_{#mu} (GeV)");
  hs_Tmu_un_tmp->GetYaxis()->SetTitle("Number of SBND Events");
  canvas->Modified();
  leg_Tmu_un->Draw();
  canvas->SaveAs(title.c_str());
  delete hs_Tmu_un_tmp;
  delete leg_Tmu_un;
}

void StackProjectionX(std::vector<TH2*> histograms, string legendlabels[], string filename, int bin, TCanvas *canvas){
  canvas->Clear();
  // Make the title for the current histogram and the file name
  // Clear the title for the new loop
  stringstream conv;
  conv.clear();

  string title;
  title.clear();

  conv << setprecision(4) << "~/Documents/PhD/CrossSections/output/tmuslicestack/" << filename 
       << histograms[0]->GetYaxis()->GetBinLowEdge(bin) << ";" << histograms[0]->GetYaxis()->GetBinLowEdge(bin+1) << ".root";
  title = conv.str();

  THStack *hs_cosmu_un_tmp = new THStack("hs_cosmu_un_tmp","");
  TLegend *leg_cosmu_un = new TLegend( 0.152, 0.598, 0.332, 0.887 );
  for ( size_t i = 0; i < histograms.size(); i++){
    //Fill THStack with projection Y of v_un for the jth bin
    TH1D *temphist = histograms[i]->ProjectionX(std::to_string(i).c_str(),bin,bin);
    temphist->SetLineWidth(1);
    temphist->SetLineColor(temphist->GetFillColor());
    temphist->SetFillStyle(1001);
    if(legendlabels[i]=="NC 1#pi"||legendlabels[i]=="NC Other"||legendlabels[i]=="NC n#pi"){
      temphist->SetFillStyle(3004);
    }
    hs_cosmu_un_tmp->Add(temphist);
    leg_cosmu_un->AddEntry(temphist,legendlabels[i].c_str(),"f");
  }
  hs_cosmu_un_tmp->Draw();
  hs_cosmu_un_tmp->GetXaxis()->SetTitle("cos#theta_{#mu}");
  hs_cosmu_un_tmp->GetYaxis()->SetTitle("Number of SBND Events");
  canvas->Modified();
  leg_cosmu_un->Draw();
  canvas->SaveAs(title.c_str());
  delete hs_cosmu_un_tmp;
  delete leg_cosmu_un;
}



int smearing() {

    gErrorIgnoreLevel = kWarning;

    // -------------------------------------------------------------------------
    //                              Open event files
    // -------------------------------------------------------------------------
    
    // Open Default
    TFile f1("~/Documents/PhD/CrossSections/mec/gntp.10000.gst.root");
    if (f1.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " Default event file is open " << endl;
    }
    // Open training file
    TFile f2("~/Documents/PhD/CrossSections/nonmec/gntp.10000.gst.root");
    if (f2.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " Training event file is open " << endl;
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
    TTree *train_tree = (TTree*) f2.Get( "gst" );
 
    // Define the log and normal histograms for the unsmeared distributions
    // El - m_mu : kinetic energy of the final state primary lepton (fspl : muon == 13 )
    // cthl      : costheta of the final state primary lepton
    // h_un      : unsmeared histogram
    // h_sm      : smeared histogram
    //

    TH2D *h_un = new TH2D("h_un"," T_{#mu} - cos#theta_{#mu} distribution before smearing ",20,-1,1,18,0,2);
    def_tree->Draw("( El - 0.10566 ):cthl>>h_un","fspl == 13 && cc","colz"); 

    TCanvas *can = new TCanvas("can","",800,600);
    can->SetLeftMargin(0.12);
    can->SetRightMargin(0.15);
    h_un->SetTitleOffset(0.75,"Y");
    h_un->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un->GetYaxis()->SetTitle("T_{#mu} (GeV)");
    h_un->Draw("COLZ");
    can->SaveAs("~/Documents/PhD/CrossSections/output/total/unsmeared.root");
    delete can;

    std::vector<TH2*> v_un;
    std::vector<TH2*> v_sm;
    std::vector<TH2*> v_sm_rec;
    
    // The same histogram definitions for the smeared distributions
    TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} distribution after smearing, impurity ",20,-1,1,18,0,2);
    
    // Take h_un and smear it
    Smear(def_tree, h_sm, v_un, v_sm, v_sm_rec);

    // The the 2 2D histograms and draw slices in Tmu and cos theta mu
    Slices( h_un, h_sm );

    SliceStack( v_un, v_sm, v_sm_rec);

    TCanvas *canvas = new TCanvas("canvas","",800,600);
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.15);

    TH2D *h_eff = new TH2D("h_eff","h_eff;cos#theta_{#mu};T_{#mu} (GeV)",20,-1,1,18,0,2);
    TH2D *h_pur = new TH2D("h_pur","h_pur;cos#theta_{#mu};T_{#mu} (GeV)",20,-1,1,18,0,2);
    TH2D *h_effpur = new TH2D("h_effpur","h_effpur;cos#theta_{#mu};T_{#mu} (GeV)",20,-1,1,18,0,2);
    TH2D *h_impur = new TH2D("h_impur","h_impur;cos#theta_{#mu};T_{#mu} (GeV)",20,-1,1,18,0,2);

    //Fill efficiency histogram with Nsig/Nmcsig
    h_eff->Add(h_sm,v_sm[5],1.,-1.);
    h_eff->Divide(h_un);
    //Fill purity histogram with Nsig/(Nsig+Nbkg)
    h_pur->Add(h_sm,v_sm[5],1.,-1.);
    h_pur->Divide(h_sm);
    //Fill efficiency*purity histogram with h_eff*h_pur
    h_effpur->Multiply(h_eff,h_pur,1.,1.);
    //Fill impurity histogram with Nbkg/(Nsig+Nbkg)
    h_impur->Divide(v_sm[5],h_sm,1.,1.);

    h_eff->SetMinimum(0.5);
    h_eff->SetMaximum(1.0);
    h_eff->SetTitleOffset(0.75,"Y");
    h_eff->Draw("COLZ");
    canvas->SaveAs("~/Documents/PhD/CrossSections/output/effpur/eff.root");
    canvas->Clear();

    h_pur->SetMinimum(0.7);
    h_pur->SetMaximum(1.0);
    h_pur->SetTitleOffset(0.75,"Y");
    h_pur->Draw("COLZ");
    canvas->SaveAs("~/Documents/PhD/CrossSections/output/effpur/pur.root");
    canvas->Clear();

    h_effpur->SetMinimum(0.5);
    h_effpur->SetMaximum(1.0);
    h_effpur->SetTitleOffset(0.75,"Y");
    h_effpur->Draw("COLZ");
    canvas->SaveAs("~/Documents/PhD/CrossSections/output/effpur/effpur.root");
    canvas->Clear();

    h_impur->SetTitleOffset(0.75,"Y");
    h_impur->Draw("COLZ");
    canvas->SaveAs("~/Documents/PhD/CrossSections/output/effpur/impur.root");
    delete canvas;

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
             std::vector<TH2*> &v_sm_rec){

    TCanvas *canv = new TCanvas("canv","",800,600);
    canv->SetLeftMargin(0.12);
    canv->SetRightMargin(0.15);

    // Initiate the random number generation
    ROOT::Math::GSLRngMT *_random_gen = new ROOT::Math::GSLRngMT;
    _random_gen->Initialize();    
    _random_gen->SetSeed( time( NULL ) );
/*
    TRandom *_rand = new TRandom( time( NULL ) );

    TH2D *hTrainTrue= new TH2D ("traintrue", "Training Truth", 20, -1, 1, 18, 0, 2);
    hTrainTrue->SetLineColor(kBlue);
    TH2D *hTrain= new TH2D ("train", "Training Measured", 20, -1, 1, 18, 0, 2);
    hTrain->SetLineColor(kRed);
    RooUnfoldResponse response(hTrain, hTrainTrue);
    TH2D *hTrue = new TH2D("true", "Test Truth", 20, -1, 1, 18, 0, 2);
    TH2D *hMeas = new TH2D("meas", "Test Measured", 20, -1 ,1, 18, 0, 2);

    //Train on a flat distribution
    for (int i = 0; i<1000000; i++){
      double costheta = _rand->Uniform(-1,1);
      double tmu = _rand->Uniform(0,2);
      double costheta_sm = SmearCosTheta(costheta,_random_gen);
      double tmu_sm = SmearKE(tmu,_random_gen);
      hTrainTrue->Fill(costheta,tmu);
      if(tmu>0.05){
        hTrain->Fill(costheta_sm,tmu_sm);
        response.Fill(costheta_sm,tmu_sm,costheta,tmu);
      }
      else response.Miss(costheta,tmu);
    }
*/
    // Stacked unsmeared histogram - truth variables
    THStack *hs_un = new THStack("hs_un","Stacked Histograms");
    TH2D *h_ccqe_un = new TH2D("h_ccqe_un","CCQE",20,-1,1,18,0,2);
    TH2D *h_ccres_un = new TH2D("h_ccres_un","CCRES",20,-1,1,18,0,2);
    TH2D *h_ccdis_un = new TH2D("h_ccdis_un","CCDIS",20,-1,1,18,0,2);
    TH2D *h_cccoh_un = new TH2D("h_cccoh_un","CCCOH",20,-1,1,18,0,2);
    TH2D *h_ccmec_un = new TH2D("h_ccmec_un","CCMEC",20,-1,1,18,0,2);

    // Stacked smeared histogram - truth variables
    THStack *hs_sm = new THStack("hs_sm","Stacked Histograms");
    TH2D *h_ccqe_sm = new TH2D("h_ccqe_sm","CCQE",20,-1,1,18,0,2);
    TH2D *h_ccres_sm = new TH2D("h_ccres_sm","CCRES",20,-1,1,18,0,2);
    TH2D *h_ccdis_sm = new TH2D("h_ccdis_sm","CCDIS",20,-1,1,18,0,2);
    TH2D *h_cccoh_sm = new TH2D("h_cccoh_sm","CCCOH",20,-1,1,18,0,2);
    TH2D *h_nc_sm = new TH2D("h_nc_sm","NCnpi",20,-1,1,18,0,2);
    TH2D *h_ccmec_sm = new TH2D("h_ccmec_sm","CCMEC",20,-1,1,18,0,2);

    // Stacked smeared histogram - reco variables
    THStack *hs_sm_rec = new THStack("hs_sm_rec","Stacked Histograms");
    TH2D *h_cc0pi0p_sm = new TH2D("h_cc0pi0p_sm","CC 0pi 0p",20,-1,1,18,0,2);
    TH2D *h_cc0pi1p_sm = new TH2D("h_cc0pi1p_sm","CC 0pi 1p",20,-1,1,18,0,2);
    TH2D *h_cc0pi2p_sm = new TH2D("h_cc0pi2p_sm","CC 0pi 2p",20,-1,1,18,0,2);
    TH2D *h_cc0pi3p_sm = new TH2D("h_cc0pi3p_sm","CC 0pi 3p",20,-1,1,18,0,2);
    TH2D *h_cc0pinp_sm = new TH2D("h_cc0pinp_sm","CC 0pi >3p",20,-1,1,18,0,2);
    TH2D *h_cc1pip_sm = new TH2D("h_cc1pip_sm","CC 1pi+",20,-1,1,18,0,2);
    TH2D *h_cc1pim_sm = new TH2D("h_cc1pim_sm","CC 1pi-",20,-1,1,18,0,2);
    TH2D *h_cc1pi0_sm = new TH2D("h_cc1pi0_sm","CC 1pi0",20,-1,1,18,0,2);
    TH2D *h_cc2pipm_sm = new TH2D("h_cc2pipm_sm","CC 2pi+/-",20,-1,1,18,0,2);
    TH2D *h_ccother_sm = new TH2D("h_ccother_sm","CCOther",20,-1,1,18,0,2);
    TH2D *h_nc1pi_sm = new TH2D("h_nc1pi_sm","NC 1pi",20,-1,1,18,0,2);
    TH2D *h_ncother_sm = new TH2D("h_ncother_sm","NCOther",20,-1,1,18,0,2);

    // Get the branches
    TBranch *b_nf    = tree->GetBranch( "nf" );
    TBranch *b_nfp   = tree->GetBranch( "nfp" );
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
    TBranch *b_mec   = tree->GetBranch( "mec" );
    
    // Number of events in the TTree
    int n_values = tree->GetEntries();
   
    double m_mu = 0.10566; // Muon mass, GeV
    double m_pi = 0.13957; // Charged pion mass, GeV

    // Counters for different types of neutrino interaction
    int ccqe_count_un = 0; int ccqe_count_sm = 0;
    int ccres_count_un = 0; int ccres_count_sm = 0;
    int ccdis_count_un = 0; int ccdis_count_sm = 0;
    int cccoh_count_un = 0; int cccoh_count_sm = 0;
    int ccmec_count_un = 0; int ccmec_count_sm = 0;
    int other_count_un = 0; int other_count_sm = 0;

    // Total counters
    int ccinc_count_un = 0;
    int ccinc_count_sm = 0;
    int imp_count = 0;
    int nc_npi_count = 0;

    // Loop over the entries of the tree and calculate the kinetic energies 
    // of the muons and pions and define the impurity
    for ( int i = 0; i < n_values; ++i ){
        
        tree->GetEntry(i);
     
        double T_mu, e_mu;

        // Get values from the branches
        int nf   = b_nf->GetLeaf("nf")->GetValue();
        int nfp  = b_nfp->GetLeaf("nfp")->GetValue();
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
        double mec   = b_mec->GetLeaf( "mec" )->GetValue();

        // Kinetic energy of the final state primary lepton
        T_mu = El - m_mu;

        //Count the CC numu event types before smearing
        if ( fspl == 13 && cc == 1 ) {

          if ( qel == 1 ) { h_ccqe_un->Fill(cthl,T_mu); ccqe_count_un++; }
          else if ( res == 1 ) { h_ccres_un->Fill(cthl,T_mu); ccres_count_un++; }
          else if ( dis == 1 ) { h_ccdis_un->Fill(cthl,T_mu); ccdis_count_un++; }
          else if ( coh == 1 ) { h_cccoh_un->Fill(cthl,T_mu); cccoh_count_un++; }
          else if ( mec == 1 ) { h_ccmec_un->Fill(cthl,T_mu); ccmec_count_un++; }
          else { other_count_un++; }

          ccinc_count_un++;

          double cthl_smeared = SmearCosTheta(cthl,_random_gen);
          double T_mu_smeared = SmearKE(T_mu,_random_gen);

//          hTrue->Fill(cthl,T_mu);

          //Count the CC numu event types after smearing
          if ( T_mu > 0.05 ){
            //Smear the muon energy and angle
            if ( qel == 1 ) { h_ccqe_sm->Fill(cthl_smeared,T_mu_smeared); ccqe_count_sm++; }
            else if ( res == 1 ) { h_ccres_sm->Fill(cthl_smeared,T_mu_smeared); ccres_count_sm++; }
            else if ( dis == 1 ) { h_ccdis_sm->Fill(cthl_smeared,T_mu_smeared); ccdis_count_sm++; }
            else if ( coh == 1 ) { h_cccoh_sm->Fill(cthl_smeared,T_mu_smeared); cccoh_count_sm++; }
            else if ( mec == 1 ) { h_ccmec_sm->Fill(cthl_smeared,T_mu_smeared); ccmec_count_sm++; }
            else { other_count_sm++; }

            //Reco variables
            if ( nfpip+nfpim+nfpi0 == 0 ){
              if ( nfp == 0 ) { h_cc0pi0p_sm->Fill(cthl_smeared,T_mu_smeared); }
              else if ( nfp == 1 ) { h_cc0pi1p_sm->Fill(cthl_smeared,T_mu_smeared); }
              else if ( nfp == 2 ) { h_cc0pi2p_sm->Fill(cthl_smeared,T_mu_smeared); }
              else if ( nfp == 3 ) { h_cc0pi3p_sm->Fill(cthl_smeared,T_mu_smeared); }
              else { h_cc0pinp_sm->Fill(cthl_smeared,T_mu_smeared); }
            }
            else if ( nfpip == 1 && nfpim+nfpi0 == 0 ) { h_cc1pip_sm->Fill(cthl_smeared,T_mu_smeared); }
            else if ( nfpim == 1 && nfpip+nfpi0 == 0 ) { h_cc1pim_sm->Fill(cthl_smeared,T_mu_smeared); }
            else if ( nfpi0 == 1 && nfpip+nfpim == 0 ) { h_cc1pi0_sm->Fill(cthl_smeared,T_mu_smeared); }
            else if ( nfpip+nfpim == 2 ) { h_cc2pipm_sm->Fill(cthl_smeared,T_mu_smeared); }
            else { h_ccother_sm->Fill(cthl_smeared,T_mu_smeared); }

            h_smeared->Fill(SmearCosTheta(cthl,_random_gen),SmearKE(T_mu,_random_gen));
            ccinc_count_sm++;

  //          hMeas->Fill(cthl_smeared,T_mu_smeared);
          }
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
              random = rand() % 5 + 1; // CHANGE FROM 20% TO 2%
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
            double cos_pi_smeared = SmearCosTheta(prev_cos_pi,_random_gen);
            double T_pi_smeared = SmearKE(T_pi,_random_gen);
            h_nc_sm->Fill(cos_pi_smeared, T_pi_smeared);
            h_smeared->Fill(cos_pi_smeared, T_pi_smeared);
            imp_count++;

            if ( nfpip+nfpim==1 ){ h_nc1pi_sm->Fill(cos_pi_smeared,T_pi_smeared); }
            else { h_ncother_sm->Fill(cos_pi_smeared,T_pi_smeared); }

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
        <<"| CCMEC =          "<<setfill(' ')<<setw(6)<<ccmec_count_un<<" | CCMEC =                 "<<setfill(' ')<<setw(6)<<ccmec_count_sm<<" |"<<endl
        <<"| Others =         "<<setfill(' ')<<setw(6)<<other_count_un<<" | Others =                "<<setfill(' ')<<setw(6)<<other_count_sm<<" |"<<endl
        <<" ---------------------------------------------------------- "<<endl;

    //Store the unsmeared component histograms
    v_un.push_back(h_ccqe_un);
    v_un.push_back(h_ccmec_un);
    v_un.push_back(h_ccres_un);
    v_un.push_back(h_ccdis_un);
    v_un.push_back(h_cccoh_un);

    //Store the smeared component histograms
    v_sm.push_back(h_nc_sm);
    v_sm.push_back(h_ccqe_sm);
    v_sm.push_back(h_ccmec_sm);
    v_sm.push_back(h_ccres_sm);
    v_sm.push_back(h_ccdis_sm);
    v_sm.push_back(h_cccoh_sm);

    //Store the smeared reco component histograms
    v_sm_rec.push_back(h_nc1pi_sm);
    v_sm_rec.push_back(h_ncother_sm);
    v_sm_rec.push_back(h_cc1pim_sm);
    v_sm_rec.push_back(h_cc1pip_sm);
    v_sm_rec.push_back(h_cc1pi0_sm);
    v_sm_rec.push_back(h_cc2pipm_sm);
    v_sm_rec.push_back(h_ccother_sm);
    v_sm_rec.push_back(h_cc0pi0p_sm);
    v_sm_rec.push_back(h_cc0pi1p_sm);
    v_sm_rec.push_back(h_cc0pi2p_sm);
    v_sm_rec.push_back(h_cc0pi3p_sm);
    v_sm_rec.push_back(h_cc0pinp_sm);
    
    //Format the smeared total 2D histogram and write to file
    h_smeared->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_smeared->GetYaxis()->SetTitle("T_{#mu} (GeV)");
    h_smeared->SetTitleOffset(0.75,"Y");
    h_smeared->Draw("COLZ");
    canv->SaveAs("~/Documents/PhD/CrossSections/output/total/smeared.root");
    canv->Clear();

    //Format the unsmeared stacked histogram and write to file
    h_cccoh_un->SetFillColor(kViolet+2);
    h_ccdis_un->SetFillColor(kYellow+1);
    h_ccres_un->SetFillColor(kOrange+7);
    h_ccqe_un->SetFillColor(kRed);
    h_ccmec_un->SetFillColor(kGreen+2);
    hs_un->Add(h_cccoh_un);
    hs_un->Add(h_ccdis_un);
    hs_un->Add(h_ccres_un);
    hs_un->Add(h_ccqe_un);
    hs_un->Add(h_ccmec_un);
    hs_un->Draw();
    canv->SaveAs("~/Documents/PhD/CrossSections/output/total/unsmeared_stacked.root");
    canv->Clear();

    //Format the smeared stacked histograms and write to file
    h_cccoh_sm->SetFillColor(kViolet+2);
    h_ccdis_sm->SetFillColor(kYellow+1);
    h_ccres_sm->SetFillColor(kOrange+7);
    h_ccqe_sm->SetFillColor(kRed);
    h_nc_sm->SetFillColor(kTeal-1);
    h_ccmec_sm->SetFillColor(kGreen+2);
    hs_sm->Add(h_nc_sm);
    hs_sm->Add(h_cccoh_sm);
    hs_sm->Add(h_ccdis_sm);
    hs_sm->Add(h_ccres_sm);
    hs_sm->Add(h_ccqe_sm);
    hs_sm->Add(h_ccmec_sm);
    hs_sm->Draw();
    canv->SaveAs("~/Documents/PhD/CrossSections/output/total/smeared_stacked.root");
    canv->Clear();

    //Format the smeared reco etacked histograms and write to file
    h_nc1pi_sm->SetFillColor(kTeal-1);
    h_ncother_sm->SetFillColor(kSpring+2);
    h_cc1pim_sm->SetFillColor(kCyan+2);
    h_cc1pip_sm->SetFillColor(kYellow+1);
    h_cc1pi0_sm->SetFillColor(kAzure+2);
    h_cc2pipm_sm->SetFillColor(kViolet+2);
    h_ccother_sm->SetFillColor(kPink+2);
    h_cc0pi0p_sm->SetFillColor(kRed);
    h_cc0pi1p_sm->SetFillColor(kGreen+2);
    h_cc0pi2p_sm->SetFillColor(kOrange+7); 
    h_cc0pi3p_sm->SetFillColor(kBlue);
    h_cc0pinp_sm->SetFillColor(kMagenta+1);
    hs_sm_rec->Add(h_nc1pi_sm);
    hs_sm_rec->Add(h_ncother_sm);
    hs_sm_rec->Add(h_cc1pim_sm);
    hs_sm_rec->Add(h_cc1pip_sm);
    hs_sm_rec->Add(h_cc1pi0_sm);
    hs_sm_rec->Add(h_cc2pipm_sm);
    hs_sm_rec->Add(h_ccother_sm);
    hs_sm_rec->Add(h_cc0pi0p_sm);
    hs_sm_rec->Add(h_cc0pi1p_sm);
    hs_sm_rec->Add(h_cc0pi2p_sm);
    hs_sm_rec->Add(h_cc0pi3p_sm);
    hs_sm_rec->Add(h_cc0pinp_sm);
    hs_sm_rec->Draw("LEGO1");
    canv->SaveAs("~/Documents/PhD/CrossSections/output/total/reco_stacked.root");
    canv->Clear();
/*
    //Create stacked histogram for just Tmu dependence before smearing
    THStack *hs_Tmu_un = new THStack("hs_Tmu_un","");
    hs_Tmu_un->Add(h_cccoh_un->ProjectionY());
    hs_Tmu_un->Add(h_ccdis_un->ProjectionY());
    hs_Tmu_un->Add(h_ccres_un->ProjectionY());
    hs_Tmu_un->Add(h_ccqe_un->ProjectionY());
    hs_Tmu_un->Add(h_ccmec_un->ProjectionY());
    hs_Tmu_un->GetXaxis()->SetTitle("T_{#mu} (GeV)");
    hs_Tmu_un->GetYaxis()->SetTitle("Number of SBND Events");
    hs_Tmu_un->Draw();

    //Create stacked histogram for Tmu after smearing
    THStack *hs_Tmu_sm = new THStack("hs_Tmu_sm","");
    hs_Tmu_sm->Add(h_nc_sm->ProjectionY());
    hs_Tmu_sm->Add(h_cccoh_sm->ProjectionY());
    hs_Tmu_sm->Add(h_ccdis_sm->ProjectionY());
    hs_Tmu_sm->Add(h_ccres_sm->ProjectionY());
    hs_Tmu_sm->Add(h_ccqe_sm->ProjectionY());
    hs_Tmu_sm->Add(h_ccmec_sm->ProjectionY());
    hs_Tmu_sm->Draw("SAME");
    canv->SaveAs("~/Documents/PhD/CrossSections/output/total/Tmu_total.root");
    canv->Clear();

    //Create stacked histogram for cos theta mu before smearing
    THStack *hs_ctmu_un = new THStack("hs_ctmu_un","");
    hs_ctmu_un->Add(h_cccoh_un->ProjectionX());
    hs_ctmu_un->Add(h_ccdis_un->ProjectionX());
    hs_ctmu_un->Add(h_ccres_un->ProjectionX());
    hs_ctmu_un->Add(h_ccqe_un->ProjectionX());
    hs_ctmu_un->Add(h_ccmec_un->ProjectionX());
    hs_ctmu_un->GetXaxis()->SetTitle("cos#theta_{#mu}");
    hs_ctmu_un->GetYaxis()->SetTitle("Number of SBND Events");
    hs_ctmu_un->Draw();

    //Create stacked histogram for cos theta mu after smearing
    THStack *hs_ctmu_sm = new THStack("hs_ctmu_sm","");
    hs_ctmu_sm->Add(h_nc_sm->ProjectionX());
    hs_ctmu_sm->Add(h_cccoh_sm->ProjectionX());
    hs_ctmu_sm->Add(h_ccdis_sm->ProjectionX());
    hs_ctmu_sm->Add(h_ccres_sm->ProjectionX());
    hs_ctmu_sm->Add(h_ccqe_sm->ProjectionX());
    hs_ctmu_sm->Add(h_ccmec_sm->ProjectionX());
    hs_ctmu_sm->Draw("SAME");
    canv->SaveAs("~/Documents/PhD/CrossSections/output/total/cosmu_total.root");
    canv->Clear();
*/
    delete canv;
/*
    RooUnfoldBayes unfold(&response, hMeas, 4);
    TH2D* hReco = (TH2D*)unfold.Hreco();
    unfold.PrintTable(cout,hTrue);

    TCanvas *c_Unf = new TCanvas("c_Unf", "", 800, 600);
    c_Unf->Divide(2,4,0,0);

    TH1D* hTrainX; TH1D* hTrainY;
    TH1D* hTrainTrueX; TH1D* hTrainTrueY;
    hTrainX = hTrain->ProjectionX(); hTrainY = hTrain->ProjectionY();
    hTrainTrueX = hTrainTrue->ProjectionX(); hTrainTrueY = hTrainTrue->ProjectionY();
    TH1D* hRecoX; TH1D* hRecoY;
    hRecoX= hReco->ProjectionX(); hRecoY= hReco->ProjectionY();
    hRecoX->SetMarkerStyle(kFullDotLarge); hRecoY->SetMarkerStyle(kFullDotLarge);
    TH1D* hTrueX; TH1D* hTrueY;
    TH1D* hMeasX; TH1D* hMeasY;
    hTrueX = hTrue->ProjectionX(); hTrueY = hTrue->ProjectionY();
    hMeasX = hMeas->ProjectionX(); hMeasY = hMeas->ProjectionY();

    c_Unf->cd(1);
    hTrainTrueX->Draw();
    hTrainX->Draw("SAME");
    c_Unf->cd(2);
    hTrainTrueY->Draw();
    hTrainY->Draw("SAME");
    c_Unf->cd(3);
    hTrueX->Draw();
    hMeasX->Draw("SAME");
    hRecoX->Draw("SAME");
    c_Unf->cd(4);
    hTrueY->Draw();
    hMeasY->Draw("SAME");
    hRecoY->Draw("SAME");
    c_Unf->cd(5);
    hTrue->Draw("COLZ");
    c_Unf->cd(6);
    hMeas->Draw("COLZ");
    c_Unf->cd(7);
    hReco->Draw("COLZ");

    c_Unf->SaveAs("unfolded.root");
    delete c_Unf;
*/
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
   
    TLegend *leg_T = new TLegend( 0.152, 0.704, 0.469, 0.88 );

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

        conv << setprecision(4) << "~/Documents/PhD/CrossSections/output/tmuslicecomp/Tmu_slice_" << low_edge_T << ";" << upp_edge_T << ".root";
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

        leg_T->Draw();
        c_Tmu->Print(file_name);

    } 
   
    delete h_Tmu;
    delete h_Tmu_sm;

    delete c_Tmu;
    
    delete leg_T;
    
    TCanvas *c_cosmu = new TCanvas ( "c_cosmu", "", 800, 600 );
    
    TLegend *leg_c = new TLegend( 0.584, 0.727, 0.902, 0.903 );

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

        conv1 << setprecision(4) << "~/Documents/PhD/CrossSections/output/cosslicecomp/cos_thetamu_slice_" << low_edge_cos << ";" << upp_edge_cos << ".root";
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

        h_cosmu->SetTitle(hist_name1);
        h_cosmu->GetXaxis()->SetTitle("T_{#mu} (GeV)");   
        h_cosmu->GetYaxis()->SetTitle("Number of SBND events");   
        h_cosmu->SetLineColor( kRed + 2 );
        h_cosmu_sm->SetLineColor( kGreen + 2 );
        h_cosmu->Draw();
        h_cosmu_sm->Draw("same");

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
void SliceStack ( std::vector<TH2*> v_un, std::vector<TH2*> v_sm, std::vector<TH2*> v_sm_rec){

    TCanvas *canvas1 = new TCanvas ( "canvas1", "", 800, 600 );

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = v_un[0]->GetNbinsX(); // Cos theta
    int y_bins = v_un[0]->GetNbinsY(); // Tmu

    string unsmeared_labs[5] = {"CCQE","CCMEC","CCRES","CCDIS","CCCOH"};
    string smeared_labs[6] = {"NC n#pi","CCQE","CCMEC","CCRES","CCDIS","CCCOH"};
    string reco_labs[12] = {"NC 1#pi","NC Other","CC1#pi^{-}","CC1#pi^{+}","CC1#pi^{0}","CC2#pi^{#pm}","CC Other","#nu_{#mu} CC0#pi, 0p","#nu_{#mu} CC0#pi, 1p","#nu_{#mu} CC0#pi, 2p","#nu_{#mu} CC0#pi, 3p","#nu_{#mu} CC0#pi, >3p"};

    for ( int j = 1; j < x_bins+1; j++ ){

      StackProjectionY(v_un, unsmeared_labs, "un_Tmu_cosmuRange:", j, canvas1);
      StackProjectionY(v_sm, smeared_labs, "sm_Tmu_cosmuRange:", j, canvas1);
      StackProjectionY(v_sm_rec, reco_labs, "rec_sm_Tmu_cosmuRange:", j, canvas1);

    }

    for ( int j = 1; j < y_bins+1; j++ ){

      StackProjectionX(v_un, unsmeared_labs, "un_cosmu_TmuRange:", j, canvas1);
      StackProjectionX(v_sm, smeared_labs, "sm_cosmu_TmuRange:", j, canvas1);
      StackProjectionX(v_sm_rec, reco_labs, "rec_sm_cosmu_TmuRange:", j, canvas1);

    }

    delete canvas1;
    
}
