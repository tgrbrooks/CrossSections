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

#include "smearing_1D.h"

using namespace std;

// -------------------------------------------------------
//              Kinetic energy smearing function
// -------------------------------------------------------
double SmearKE(double ke, ROOT::Math::GSLRngMT *_random_gen){

  // Calculate the mean and sigma for the LogNormal function
  //      zeta  = TMath::Log( m * ( 1 / sqrt( 1 + ( var / pow( m, 2 ) ) ) ) );
  //      sigma = sqrt( log( 1 + ( var / pow( m, 2 ) ) ) );
  //      m     = expectation value = El - m_mu
  //      var   = variance = s.d.^2 = ( El - m_mu * 0.1 ) ^ 2

  double var     = TMath::Power( ( ke ) * 0.1, 2 ); //WAS 0.1 
  double sigma   = TMath::Sqrt( TMath::Log( 1 + ( var / TMath::Power( ( ke ), 2 ) ) ) );
  double zeta    = TMath::Log( ( ke ) * ( 1 / TMath::Sqrt( 1 + ( var / TMath::Power( ( ke ), 2 ) ) ) ) );
  double lognorm = _random_gen->LogNormal( zeta, sigma );
  /*double lognorm = ke + _random_gen->Gaussian(0.1);
  if (lognorm < 0) lognorm = 0;*/

  return lognorm;
}


// -------------------------------------------------------
//              Cos theta smearing function
// -------------------------------------------------------
double SmearCosTheta(double costheta, ROOT::Math::GSLRngMT *_random_gen){
    
  // Calculate the mean and sigma for the LogNormal function
  //      theta = acos(costtheta)
  //      var   = 5 degrees

  double sd_theta = TMath::Pi() / 36; // 5 degrees
  double gaus_theta = TMath::ACos( costheta ) + _random_gen->Gaussian( sd_theta );
  double gaus_costheta = TMath::Cos( gaus_theta );

  return gaus_costheta;
}

// -------------------------------------------------------
//  Function to get Y projection of stacked histograms
// -------------------------------------------------------
void StackProjectionY(std::vector<TH2*> histograms, string legendlabels[], string filename, int bin, TCanvas *canvas, TH2* h_unf, bool isTruth){
  
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

  //Create a temporary stacked histogram for the Y projection
  THStack *hs_Tmu_un_tmp = new THStack("hs_tmu_un_tmp","");
  TLegend *leg_Tmu_un = new TLegend( 0.584, 0.614, 0.901, 0.904 );
  for ( size_t i = 0; i < histograms.size(); i++){
    //Fill THStack with projection Y of ith histogram for the jth bin
    TH1D *temphist = histograms[i]->ProjectionY(std::to_string(i).c_str(),bin,bin);
    temphist->SetLineWidth(1);
    temphist->SetLineColor(temphist->GetFillColor());
    temphist->SetFillStyle(1001);
    //Change fill style for background events
    if(legendlabels[i]=="NC 1#pi"||legendlabels[i]=="NC Other"||legendlabels[i]=="NC n#pi"){
      temphist->SetFillStyle(3004);
    }
    hs_Tmu_un_tmp->Add(temphist);
    leg_Tmu_un->AddEntry(temphist,legendlabels[i].c_str(),"f");
  }
  hs_Tmu_un_tmp->Draw();

  // If the histograms are of unsmeared truth variables plot the unfolding results on top of them
  if (isTruth){
    TH1D* hUnf = h_unf->ProjectionY("unf",bin,bin);
    hUnf->SetLineWidth(1);
    hUnf->SetLineColor(kBlack);
    hUnf->SetMarkerStyle(20);
    hUnf->SetMarkerSize(0.6);
    hUnf->Draw("E1 SAME");
  }
  hs_Tmu_un_tmp->GetXaxis()->SetTitle("T_{#mu} (GeV)");
  hs_Tmu_un_tmp->GetYaxis()->SetTitle("Number of SBND Events");
  canvas->Modified();
  leg_Tmu_un->Draw();
  canvas->SaveAs(title.c_str());

  //Clean up
  delete hs_Tmu_un_tmp;
  delete leg_Tmu_un;

}

// -------------------------------------------------------
//  Function to get X projection of stacked histograms
// -------------------------------------------------------
void StackProjectionX(std::vector<TH2*> histograms, string legendlabels[], string filename, int bin, TCanvas *canvas, TH2* h_unf, bool isTruth){

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

  //Create a temporary stacked histogram for the X projection
  THStack *hs_cosmu_un_tmp = new THStack("hs_cosmu_un_tmp","");
  TLegend *leg_cosmu_un = new TLegend( 0.152, 0.598, 0.332, 0.887 );
  for ( size_t i = 0; i < histograms.size(); i++){
    //Fill THStack with projection X of ith histogram for the jth bin
    TH1D *temphist = histograms[i]->ProjectionX(std::to_string(i).c_str(),bin,bin);
    temphist->SetLineWidth(1);
    temphist->SetLineColor(temphist->GetFillColor());
    temphist->SetFillStyle(1001);
    //Change fill style for background events
    if(legendlabels[i]=="NC 1#pi"||legendlabels[i]=="NC Other"||legendlabels[i]=="NC n#pi"){
      temphist->SetFillStyle(3004);
    }
    hs_cosmu_un_tmp->Add(temphist);
    leg_cosmu_un->AddEntry(temphist,legendlabels[i].c_str(),"f");
  }
  hs_cosmu_un_tmp->Draw();

  // If the histograms are of unsmeared truth variables plot the unfolding results on top of them
  if (isTruth){
    TH1D* hUnf = h_unf->ProjectionX("unf",bin,bin);
    hUnf->SetLineWidth(1);
    hUnf->SetLineColor(kBlack);
    hUnf->SetMarkerStyle(20);
    hUnf->SetMarkerSize(0.6);
    hUnf->Draw("E1 SAME");
  }
  hs_cosmu_un_tmp->GetXaxis()->SetTitle("cos#theta_{#mu}");
  hs_cosmu_un_tmp->GetYaxis()->SetTitle("Number of SBND Events");
  canvas->Modified();
  leg_cosmu_un->Draw();
  canvas->SaveAs(title.c_str());

  // Clean up
  delete hs_cosmu_un_tmp;
  delete leg_cosmu_un;
}

// -------------------------------------------------------
//  Function to get correlation histogram from unfolding
// -------------------------------------------------------
TH2D* CorrelationHist (const TMatrixD& cov,
                       const char* name, const char* title,
                       Double_t lo, Double_t hi) 
{
  //Pulled from RooUnfold unit tests
  Int_t nb= cov.GetNrows();
  TH2D* h= new TH2D (name, title, nb, lo, hi, nb, lo, hi);
  h->SetAxisRange (-1.0, 1.0, "Z");
  for(int i=0; i < nb; i++)
    for(int j=0; j < nb; j++) {
      Double_t Viijj= cov(i,i)*cov(j,j);
      if (Viijj>0.0) h->SetBinContent (i+1, j+1, cov(i,j)/sqrt(Viijj));
    }   
  return h;
}

// -------------------------------------------------------
//  Function to test smearing on 2D delta function
// -------------------------------------------------------
void TestSmearing(double cos, double tmu) {

  //Create canvases and histograms
  TCanvas *ctest = new TCanvas("ctest","",800,600);
  TLegend *ltest = new TLegend( 0.152, 0.704, 0.469, 0.88 );
  TH2D* h_before = new TH2D("h_before","before;cos#theta_{#mu};T_{#mu} (GeV)",200,-1,1,180,0,2);
  TH2D* h_after = new TH2D("h_after","after;cos#theta_{#mu};T_{#mu} (GeV)",200,-1,1,180,0,2);

  //Initialise random number generator
  ROOT::Math::GSLRngMT *_rand_gen = new ROOT::Math::GSLRngMT;
  _rand_gen->Initialize();    
  _rand_gen->SetSeed( time( NULL ) );

  //Populate histograms
  for (int i=0; i<10000; ++i){
    double smear_cos = SmearCosTheta(cos,_rand_gen);
    h_before->Fill(cos,tmu);
    h_after->Fill(smear_cos,SmearKE(tmu,_rand_gen));
  }

  //Get X and Y projections
  TH1D* h_beforeX = h_before->ProjectionX("px6",1,180);
  TH1D* h_beforeY = h_before->ProjectionY("py6",1,200);
  TH1D* h_afterX = h_after->ProjectionX("px7",1,180);
  TH1D* h_afterY = h_after->ProjectionY("py7",1,200);

  ctest->SetLeftMargin(0.12);
  ctest->SetRightMargin(0.15);

  //Plot 2D histograms
  h_before->SetTitleOffset(0.75,"Y");
  h_before->Draw("COLZ");
  ctest->SaveAs("~/Documents/PhD/CrossSections/output/smeartest/before2D.root");
  ctest->Clear();

  h_after->SetTitleOffset(0.75,"Y");
  h_after->Draw("COLZ");
  ctest->SaveAs("~/Documents/PhD/CrossSections/output/smeartest/after2D.root");
  ctest->Clear();

  ctest->SetLeftMargin(0.12);
  ctest->SetRightMargin(0.1);

  //Plot X and Y projections
  //Cos theta mu
  h_beforeX->SetLineColor(kRed + 2);
  h_beforeX->GetYaxis()->SetTitle("Events");
  h_afterX->SetLineColor(kGreen + 2);
  ltest->AddEntry(h_beforeX,"Unsmeared","l");
  ltest->AddEntry(h_afterX,"Smeared","l");
  h_beforeX->Draw();
  h_afterX->Draw("same");
  ltest->Draw();
  ctest->SaveAs("~/Documents/PhD/CrossSections/output/smeartest/costheta.root");
  ctest->Clear();
  ltest->Clear();

  //Tmu
  h_beforeY->SetLineColor(kRed + 2);
  h_beforeY->GetYaxis()->SetTitle("Events");
  h_afterY->SetLineColor(kGreen + 2);
  ltest->AddEntry(h_beforeY,"Unsmeared","l");
  ltest->AddEntry(h_afterY,"Smeared","l");
  h_beforeY->Draw();
  h_afterY->Draw("same");
  ltest->Draw();
  ctest->SaveAs("~/Documents/PhD/CrossSections/output/smeartest/tmu.root");
  ctest->Clear();

  delete ctest;
  delete ltest;

}

// Function to transform a 2D histogram into a 1D hist
TH1D* ReBin(const char *name,TH2D* twoDhist){

  int xbins = twoDhist->GetNbinsX();
  int ybins = twoDhist->GetNbinsY();
  int newbins = xbins*ybins;
  TH1D* newhist = new TH1D(name,"",newbins,0,newbins);
  // Start from 1, as 0 is underflow bin
  for(int i=1;i<=xbins;++i){
    for(int j=1;j<=ybins;++j){
      newhist->SetBinContent(18*(i-1)+j,twoDhist->GetBinContent(i,j));
    }
  }
  return newhist;

}

// Function to transform 1D hist back into a 2D hist
TH2D* UnReBin(const char *name, TH1D* oneDhist){

  TH2D* newhist = new TH2D(name,"",20,-1,1,18,0,2);
  for(int i=1; i<=20; ++i){
    for(int j=1; j<=18; ++j){
      newhist->SetBinContent(i,j,oneDhist->GetBinContent(18*(i-1)+j));
    }
  }
  return newhist;

}

// Function to map pair of cos theta and Tmu to 1D bin
int MapTo1D(double cos, double tmu){

  double cos_width = 2./20.; 
  double tmu_width = 2./18.;
  int i = ceil((1.+cos)/cos_width); //Cos theta 2D bin
  int j = ceil(tmu/tmu_width); //Tmu 2D bin
  int bin = 18.*(i-1.)+j; //Transform to 1D bin
  if(tmu<=2.) return bin;
  else return 361;

}

int smearing_1D() {

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
    TFile f2("~/Documents/PhD/CrossSections/sortamec/gntp.10000.gst.root");
    if (f2.IsZombie()) {
       cout << " Error opening file " << endl;
       exit(-1);
    }
    else{
        cout << " Training event file is open " << endl;
    }
    //Open flux file
    TFile f3("~/Documents/PhD/CrossSections/flux/sbn_FHC_flux_hist.root");
    if (f3.IsZombie()) {
      cout << "Error opening flux file" << endl;
      exit(-1);
    }
    else{
      cout << " Flux file is open " <<endl;
    }
    TH1D *h_flux = (TH1D*)f3.Get("h_numu_110m");
    //-----------------------------------------------------------------------
    //                            Parameters to consider
    // -------------------------------------------------------------------------
    
    // Number of MiniBooNE bins (data) 
    // costhetamu : 20
    int cos_nbins = 20; double cos_min = -1.; double cos_max = 1.;
    // Tmu        : 18
    int tmu_nbins = 18; double tmu_min = 0.; double tmu_max = 2.;

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

    TH2D *h_un = new TH2D("h_un"," T_{#mu} - cos#theta_{#mu} distribution before smearing ",20,-1,1,18,0,2);
    def_tree->Draw("( El - 0.10566 ):cthl>>h_un","fspl == 13 && cc","colz"); 

    // Normalization, Rate ==> Cross section
    double cos_bins = h_un->GetXaxis()->GetBinWidth(1);
    double Tmu_bins = h_un->GetYaxis()->GetBinWidth(1);
    double flux_int = h_flux->Integral();
    double Na = 6.63e34;
    double M_fid = 1.016e8; //grams
    double A_ar = 39.948;
    double tot_tgt = (Na*M_fid)/A_ar;
    double barns = 1e38;
    double scalar_norm = barns/(cos_bins*Tmu_bins*flux_int*tot_tgt);

    TCanvas *can = new TCanvas("can","",800,600);
    can->SetLeftMargin(0.12);
    can->SetRightMargin(0.15);
    h_un->SetTitleOffset(0.75,"Y");
    h_un->GetXaxis()->SetTitle("cos#theta_{#mu}");
    h_un->GetYaxis()->SetTitle("T_{#mu} (GeV)");
    h_un->Draw("COLZ");
    can->SaveAs("~/Documents/PhD/CrossSections/output/total/unsmeared.root");

    can->Clear();
    TH1D* h_un_rb = ReBin("h_un_rb",h_un);
    h_un_rb->Draw();
    can->SaveAs("~/Documents/PhD/CrossSections/output/rbtest.root");
    can->Clear();
    TH2D* h_un_urb = UnReBin("h_un_urb",h_un_rb);
    h_un_urb->Draw("COLZ");
    can->SaveAs("~/Documents/PhD/CrossSections/output/urbtest.root");
    can->Clear();
    h_un_urb->Divide(h_un);
    h_un_urb->Draw("COLZ");
    can->SaveAs("~/Documents/PhD/CrossSections/output/urbtest_div.root");

    delete can;

    std::vector<TH2*> v_un;
    std::vector<TH2*> v_sm;
    std::vector<TH2*> v_sm_rec;
    std::vector<TH2*> v_unf;
    
    // The same histogram definitions for the smeared distributions
    TH2D *h_sm = new TH2D("h_sm"," T_{#mu} - cos#theta_{#mu} distribution after smearing, impurity ",20,-1,1,18,0,2);
    
    // Take h_un and smear it
    Smear(def_tree, h_sm, v_un, v_sm, v_sm_rec, scalar_norm);

    // Method can be "bayes", "svd" or "tunfold"
    Unfold(train_tree, h_sm, h_un, v_unf, "svd");

    // The the 2 2D histograms and draw slices in Tmu and cos theta mu
    Slices( h_un, h_sm, v_unf );

    SliceStack( v_un, v_sm, v_sm_rec, v_unf );

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

    TestSmearing(0.5,0.5);

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
             std::vector<TH2*> &v_sm_rec,
             double scalar_norm){

    TCanvas *canv = new TCanvas("canv","",800,600);
    canv->SetLeftMargin(0.12);
    canv->SetRightMargin(0.15);

    // Initiate the random number generation
    ROOT::Math::GSLRngMT *_random_gen = new ROOT::Math::GSLRngMT;
    _random_gen->Initialize();    
    _random_gen->SetSeed( time( NULL ) );

    TRandom *_rand = new TRandom( time( NULL ) );
   
    double m_mu = 0.10566; // Muon mass, GeV
    double m_pi = 0.13957; // Charged pion mass, GeV

    // Stacked unsmeared histogram - truth variables
    THStack *hs_un = new THStack("hs_un","Stacked Histograms;cos#theta_{#mu};T_{#mu} (GeV)");
    TH2D *h_ccqe_un = new TH2D("h_ccqe_un","CCQE",20,-1,1,18,0,2);
    TH2D *h_ccres_un = new TH2D("h_ccres_un","CCRES",20,-1,1,18,0,2);
    TH2D *h_ccdis_un = new TH2D("h_ccdis_un","CCDIS",20,-1,1,18,0,2);
    TH2D *h_cccoh_un = new TH2D("h_cccoh_un","CCCOH",20,-1,1,18,0,2);
    TH2D *h_ccmec_un = new TH2D("h_ccmec_un","CCMEC",20,-1,1,18,0,2);

    // Stacked smeared histogram - truth variables
    THStack *hs_sm = new THStack("hs_sm","Stacked Histograms;cos#theta_{#mu};T_{#mu} (GeV)");
    TH2D *h_ccqe_sm = new TH2D("h_ccqe_sm","CCQE",20,-1,1,18,0,2);
    TH2D *h_ccres_sm = new TH2D("h_ccres_sm","CCRES",20,-1,1,18,0,2);
    TH2D *h_ccdis_sm = new TH2D("h_ccdis_sm","CCDIS",20,-1,1,18,0,2);
    TH2D *h_cccoh_sm = new TH2D("h_cccoh_sm","CCCOH",20,-1,1,18,0,2);
    TH2D *h_nc_sm = new TH2D("h_nc_sm","NCnpi",20,-1,1,18,0,2);
    TH2D *h_ccmec_sm = new TH2D("h_ccmec_sm","CCMEC",20,-1,1,18,0,2);

    // Stacked smeared histogram - reco variables
    THStack *hs_sm_rec = new THStack("hs_sm_rec","Stacked Histograms;cos#theta_{#mu};T_{#mu} (GeV)");
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
    TBranch *b_fspl  = tree->GetBranch( "fspl" );
    TBranch *b_cthf  = tree->GetBranch( "cthf" );
    TBranch *b_Ef    = tree->GetBranch( "Ef" );
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

            h_smeared->Fill(cthl_smeared, T_mu_smeared);
            ccinc_count_sm++;

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
    hs_un->Draw("COLZ");
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
    hs_sm->Draw("COLZ");
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

}

void Unfold(TTree *traintree, TH2D *hMeas, TH2D *hTrue, std::vector<TH2*> &v_unf, std::string method){

    TCanvas *canv = new TCanvas("canv","",800,600);
    canv->SetLeftMargin(0.12);
    canv->SetRightMargin(0.15);

    // Initiate the random number generation
    ROOT::Math::GSLRngMT *_random_gen = new ROOT::Math::GSLRngMT;
    _random_gen->Initialize();    
    _random_gen->SetSeed( time( NULL ) );

    TRandom *_rand = new TRandom( time( NULL ) );

    TH2D *hTrainTrue= new TH2D ("traintrue", "Training Truth;cos#theta_{#mu};T_{#mu} (GeV)", 20, -1, 1, 18, 0, 2);
    hTrainTrue->SetLineColor(kRed + 2);
    TH2D *hTrain= new TH2D ("train", "Training Measured;cos#theta_{#mu};T_{#mu} (GeV)", 20, -1, 1, 18, 0, 2);
    hTrain->SetLineColor(kGreen + 2);
    TH2D *hTrainFake= new TH2D ("trainfake", "Training Fakes;cos#theta_{#mu};T_{#mu} (GeV)", 20, -1, 1, 18, 0, 2);
    TH1D *hSize = ReBin("hSize",hTrain);

    TH1D *hTrainTrue_1D= new TH1D ("traintrue1D", "Training Truth;cos#theta_{#mu};T_{#mu} (GeV)", 360, 0, 360);
    hTrainTrue_1D->SetLineColor(kRed + 2);
    TH1D *hTrain_1D= new TH1D ("train1D", "Training Measured;cos#theta_{#mu};T_{#mu} (GeV)", 360, 0, 360);
    //TH1D* hTrain_1D = ReBin("train1D",hTrain);
    hTrain_1D->SetLineColor(kGreen + 2);
    TH1D *hTrainFake_1D= new TH1D ("trainfake1D", "Training Fakes;cos#theta_{#mu};T_{#mu} (GeV)", 360, 0, 360);

    RooUnfoldResponse response(hSize,hSize);
    hTrue->SetLineColor(kRed + 2);
    hMeas->SetLineColor(kGreen + 2);

    //Train on events generated without MEC
    // Get the branches
    TBranch *bt_nf    = traintree->GetBranch( "nf" );
    TBranch *bt_nfp   = traintree->GetBranch( "nfp" );
    TBranch *bt_cthl  = traintree->GetBranch( "cthl" );
    TBranch *bt_El    = traintree->GetBranch( "El" );
    TBranch *bt_pl    = traintree->GetBranch( "pl" );
    TBranch *bt_fspl  = traintree->GetBranch( "fspl" );
    TBranch *bt_cthf  = traintree->GetBranch( "cthf" );
    TBranch *bt_Ef    = traintree->GetBranch( "Ef" );
    TBranch *bt_pdgf  = traintree->GetBranch( "pdgf" );
    TBranch *bt_cc    = traintree->GetBranch( "cc" );
    TBranch *bt_nc    = traintree->GetBranch( "nc" );
    
    // Number of events in the TTree
    int nt_values = traintree->GetEntries();
   
    double m_mu = 0.10566; // Muon mass, GeV
    double m_pi = 0.13957; // Charged pion mass, GeV


    // Loop over the entries of the tree and calculate the kinetic energies 
    // of the muons and pions and define the impurity
    for ( int i = 0; i < nt_values; ++i ){
        
        traintree->GetEntry(i);
     
        double T_mu, e_mu;

        // Get values from the branches
        int nf   = bt_nf->GetLeaf("nf")->GetValue();
        int nfp  = bt_nfp->GetLeaf("nfp")->GetValue();
        double fspl  = bt_fspl->GetLeaf( "fspl" )->GetValue(); 
        double cc    = bt_cc->GetLeaf( "cc" )->GetValue(); 
        double nc    = bt_nc->GetLeaf( "nc" )->GetValue(); 
        double El    = bt_El->GetLeaf( "El" )->GetValue();
        double pl    = bt_pl->GetLeaf( "pl" )->GetValue();
        double cthl  = bt_cthl->GetLeaf( "cthl" )->GetValue();

        // Kinetic energy of the final state primary lepton
        T_mu = El - m_mu;

        //Count the CC numu event types before smearing
        if ( fspl == 13 && cc == 1 ) {

          double cthl_smeared = SmearCosTheta(cthl,_random_gen);
          double T_mu_smeared = SmearKE(T_mu,_random_gen);

          hTrainTrue->Fill(cthl,T_mu);
          hTrainTrue_1D->Fill(MapTo1D(cthl,T_mu)-0.5);

          //Count the CC numu event types after smearing
          if ( T_mu > 0.05 ){
 
            hTrain->Fill(cthl_smeared, T_mu_smeared);
            hTrain_1D->Fill(MapTo1D(cthl_smeared,T_mu_smeared)-0.5);
            response.Fill(MapTo1D(cthl_smeared, T_mu_smeared)-0.5, MapTo1D(cthl, T_mu)-0.5);
          }
          else {
            response.Miss(MapTo1D(cthl, T_mu)-0.5);
          }
        }
        // Find NCnpi events
        if ( nc == 1 ){ 

          double prev_e_pi = 0;
          double prev_cos_pi = 0;

          //Loop over all the final state hadronic particles
          for (int j = 0; j<nf; ++j) {

            bt_pdgf->GetEntry(i);
            bt_cthf->GetEntry(i);
            bt_Ef->GetEntry(i);

            int pdgf = bt_pdgf->GetLeaf("pdgf")->GetValue(j);
            double e_pi = bt_Ef->GetLeaf("Ef")->GetValue(j);
            double cos_pi = bt_cthf->GetLeaf("cthf")->GetValue(j);

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
            hTrain->Fill(cos_pi_smeared, T_pi_smeared);
            hTrain_1D->Fill(MapTo1D(cos_pi_smeared, T_pi_smeared)-0.5);
            hTrainFake->Fill(cos_pi_smeared, T_pi_smeared);
            hTrainFake_1D->Fill(MapTo1D(cos_pi_smeared,T_pi_smeared)-0.5);
            response.Fake(MapTo1D(cos_pi_smeared, T_pi_smeared)-0.5);
          }
        }
    }

    TH1D *hMeas_1D = ReBin("hMeas_1D",hMeas);
    TH1D *hTrue_1D = ReBin("hTrue_1D",hTrue);

    hTrue_1D->Draw();
    hTrain_1D->Draw("SAME");
    canv->SaveAs("~/Documents/PhD/CrossSections/output/unfolding/bintest.root");
    canv->Clear();

    // Do the unfolding
    std::vector<TH1D*> vReco;
    if (method == "bayes"){
      cout<<"Using the iterative bayesian method"<<endl;
      RooUnfoldBayes unfold(&response, hMeas_1D, 10);
      TH1D* hReco_1D = (TH1D*)unfold.Hreco((RooUnfold::ErrorTreatment)2);
      vReco.push_back(hReco_1D);
    }
    else if (method == "svd"){
      cout<<"Using singular value decomposition"<<endl;
      RooUnfoldSvd unfold(&response, hMeas_1D, 2);
      TH1D* hReco_1D = (TH1D*)unfold.Hreco();
      vReco.push_back(hReco_1D);
    }
    else if (method == "tunfold"){
      cout<<"Using TUnfold"<<endl;
      RooUnfoldTUnfold unfold(&response, hMeas_1D);
      TH1D* hReco_1D = (TH1D*)unfold.Hreco((RooUnfold::ErrorTreatment)2);
      vReco.push_back(hReco_1D);
    }
    else{
      cout<<"Unfolding method not recognised: Using bin by bin method"<<endl;
      RooUnfoldBinByBin unfold(&response, hMeas_1D);
      TH1D* hReco_1D = (TH1D*)unfold.Hreco();
      vReco.push_back(hReco_1D);
    }
    
    //unfold.PrintTable(cout,hTrue);
    TH2D* hReco = UnReBin("hReco",vReco[0]);

    // Seperate the 2D histograms into their X (cos) and Y (Tmu) components
    TH1D* hTrainX; TH1D* hTrainY;
    TH1D* hTrainTrueX; TH1D* hTrainTrueY;
    TH1D* hTrainFakeX; TH1D* hTrainFakeY;
    hTrainX = hTrain->ProjectionX("px3",1,18); hTrainY = hTrain->ProjectionY("py3",1,20);
    hTrainTrueX = hTrainTrue->ProjectionX("px4",1,18); hTrainTrueY = hTrainTrue->ProjectionY("py4",1,20);
    hTrainFakeX = hTrainFake->ProjectionX("px5",1,18); hTrainFakeY = hTrainFake->ProjectionY("py5",1,20);
    TH1D* hRecoX; TH1D* hRecoY;
    hRecoX= hReco->ProjectionX(); hRecoY= hReco->ProjectionY();
    hRecoX->SetMarkerStyle(kFullDotLarge); hRecoY->SetMarkerStyle(kFullDotLarge);
    TH1D* hTrueX; TH1D* hTrueY;
    TH1D* hMeasX; TH1D* hMeasY;
    hTrueX = hTrue->ProjectionX("px1",1,18); hTrueY = hTrue->ProjectionY("py1",1,20);
    hMeasX = hMeas->ProjectionX("px2",1,18); hMeasY = hMeas->ProjectionY("py2",1,20);

    TLegend *leg = new TLegend( 0.172, 0.704, 0.489, 0.88 );

    canv->SetLeftMargin(0.15);
    canv->SetRightMargin(0.1);

    //Draw the cos theta mu distributions used for training
    leg->AddEntry( hTrainTrueX, "Unsmeared Training", "l" );
    leg->AddEntry( hTrainX, "Smeared Training", "l" );
    leg->AddEntry( hTrainFakeX, "Fakes Training", "l" );
    hTrainTrueX->GetYaxis()->SetTitle("Number of SBND Events");
    hTrainTrueX->SetTitleOffset(1.0,"Y");
    hTrainTrueX->SetTitleOffset(0.9,"X");
    hTrainTrueX->Draw();
    hTrainX->Draw("SAME");
    hTrainFakeX->Draw("SAME");
    leg->Draw();
    canv->SaveAs("~/Documents/PhD/CrossSections/output/unfolding/cosTrain.root");
    canv->Clear();
    leg->Clear();

    // Draw the unfolded distributions for cos theta mu
    leg->AddEntry( hTrueX, "Truth", "l" );
    leg->AddEntry( hMeasX, "Measured", "l" );
    leg->AddEntry( hRecoX, "Unfolded", "pl" );
    hTrueX->GetYaxis()->SetTitle("Number of SBND Events");
    hTrueX->SetTitleOffset(1.0,"Y");
    hTrueX->SetTitleOffset(0.9,"X");
    hTrueX->Draw();
    hMeasX->Draw("SAME");
    hRecoX->SetLineWidth(1);
    hRecoX->SetLineColor(kBlack);
    hRecoX->SetMarkerSize(0.6);
    hRecoX->Draw("E1 SAME");
    leg->Draw();
    canv->SaveAs("~/Documents/PhD/CrossSections/output/unfolding/cosUnfolded.root");
    canv->Clear();
    delete leg;

    TLegend *leg2 = new TLegend( 0.554, 0.727, 0.872, 0.903 );

    // Draw the Tmu distributions used for training
    leg2->AddEntry( hTrainTrueY, "Unsmeared Training", "l" );
    leg2->AddEntry( hTrainY, "Smeared Training", "l" );
    leg2->AddEntry( hTrainFakeY, "Fakes Training", "l" );
    hTrainTrueY->GetYaxis()->SetTitle("Number of SBND Events");
    hTrainTrueY->SetTitleOffset(1.0,"Y");
    hTrainTrueY->SetTitleOffset(0.9,"X");
    hTrainTrueY->Draw();
    hTrainY->Draw("SAME");
    hTrainFakeY->Draw("SAME");
    leg2->Draw();
    canv->SaveAs("~/Documents/PhD/CrossSections/output/unfolding/TmuTrain.root");
    canv->Clear();
    leg2->Clear();

    // Draw the unfolded distributions for Tmu
    leg2->AddEntry( hTrueY, "Truth", "l" );
    leg2->AddEntry( hMeasY, "Measured", "l" );
    leg2->AddEntry( hRecoY, "Unfolded", "pl" );
    hTrueY->GetYaxis()->SetTitle("Number of SBND Events");
    hTrueY->SetTitleOffset(1.0,"Y");
    hTrueY->SetTitleOffset(0.9,"X");
    hTrueY->Draw();
    hMeasY->Draw("SAME");
    hRecoY->SetLineWidth(1);
    hRecoY->SetLineColor(kBlack);
    hRecoY->SetMarkerSize(0.6);
    hRecoY->Draw("E1 SAME");
    leg2->Draw();
    canv->SaveAs("~/Documents/PhD/CrossSections/output/unfolding/TmuUnfolded.root");
    canv->Clear();
    delete leg2;

    canv->SetLeftMargin(0.12);
    canv->SetRightMargin(0.15);

    // Draw the full truth distribution
    hTrue->SetTitleOffset(0.75,"Y");
    hTrue->Draw("COLZ");
    canv->SaveAs("~/Documents/PhD/CrossSections/output/unfolding/Truth.root");
    canv->Clear();

    //Draw the full measured distribution
    hMeas->SetTitleOffset(0.75,"Y");
    hMeas->Draw("COLZ");
    canv->SaveAs("~/Documents/PhD/CrossSections/output/unfolding/Measured.root");
    canv->Clear();

    // Draw the full unfolded distribution
    hReco->SetTitleOffset(0.75,"Y");
    hReco->Draw("COLZ");
    canv->SaveAs("~/Documents/PhD/CrossSections/output/unfolding/Unfolded.root");
    canv->Clear();

    // Draw histogram of difference between truth and unfolded
    TH2D *h_diff = new TH2D("h_diff","h_diff;cos#theta_{#mu};T_{#mu} (GeV)",20,-1,1,18,0,2);
    h_diff->SetTitleOffset(0.75,"Y");
    h_diff->Add(hTrue,hReco,100.,-100.);
    h_diff->Divide(hTrue);
    h_diff->Draw("COLZ");
    canv->SaveAs("~/Documents/PhD/CrossSections/output/unfolding/Difference.root");
    canv->Clear();
/*
    TH2D* hCorr= CorrelationHist (unfold.Ereco((RooUnfold::ErrorTreatment)2),
                          "corr", "Unfolded correlation matrix",
                          response.Hresponse()->GetYaxis()->GetXmin(),
                          response.Hresponse()->GetYaxis()->GetXmax());
    hCorr->Draw("COLZ");
    canv->SaveAs("~/Documents/PhD/CrossSections/output/unfolding/Covariance.root");
    delete canv;
*/
    v_unf.push_back(hReco);

}

// -------------------------------------------------------------------------
//                             Slicing function
// -------------------------------------------------------------------------
void Slices ( TH2D *h_unsmeared, TH2D *h_smeared, std::vector<TH2*> v_unf){

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = h_unsmeared->GetNbinsX(); // Cos theta
    int y_bins = h_unsmeared->GetNbinsY(); // Tmu

    TCanvas *c_Tmu   = new TCanvas ( "c_Tmu", "", 800, 600 );
   
    TLegend *leg_T = new TLegend( 0.152, 0.704, 0.469, 0.88 );

    TH1D *h_Tmu      = new TH1D ( "h_Tmu", "", x_bins, -1, 1 );
    TH1D *h_Tmu_sm   = new TH1D ( "h_Tmu_sm", "", x_bins, -1, 1 );
    TH1D *h_Tmu_unf  = new TH1D ( "h_Tmu_unf", "", x_bins, -1, 1 );
    
    leg_T->AddEntry( h_Tmu, " Unsmeared ", "l" );
    leg_T->AddEntry( h_Tmu_sm, " Smeared ", "l" );
    leg_T->AddEntry( h_Tmu_unf, " Unfolded ", "pl" );
    
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
            h_Tmu_unf->SetBinContent( j, v_unf[0]->GetBinContent(j, i) );
        }

        h_Tmu->Draw();
        h_Tmu_sm->Draw("same");
        h_Tmu_unf->Draw("E1 same");
        h_Tmu->SetTitle(hist_name);
        h_Tmu->GetXaxis()->SetTitle("cos#theta_{#mu}");   
        h_Tmu->GetYaxis()->SetTitle("Number of SBND events");   
        h_Tmu->SetLineColor( kRed + 2 );
        h_Tmu_sm->SetLineColor( kGreen + 2 );
        h_Tmu_unf->SetLineWidth(1);
        h_Tmu_unf->SetLineColor(kBlack);
        h_Tmu_unf->SetMarkerStyle(20);
        h_Tmu_unf->SetMarkerSize(0.6);

        leg_T->Draw();
        c_Tmu->Print(file_name);

    } 
   
    delete h_Tmu;
    delete h_Tmu_sm;
    delete h_Tmu_unf;

    delete c_Tmu;
    
    delete leg_T;
    
    TCanvas *c_cosmu = new TCanvas ( "c_cosmu", "", 800, 600 );
    
    TLegend *leg_c = new TLegend( 0.584, 0.727, 0.902, 0.903 );

    TH1D *h_cosmu    = new TH1D ( "h_cosmu", "", y_bins, 0, 2 );
    TH1D *h_cosmu_sm = new TH1D ( "h_cosmu_sm", "", y_bins, 0, 2 );
    TH1D *h_cosmu_unf = new TH1D ( "h_cosmu_unf", "", y_bins, 0, 2 );
     
    leg_c->AddEntry( h_cosmu, " Unsmeared ", "l" );
    leg_c->AddEntry( h_cosmu_sm, " Smeared ", "l" );
    leg_c->AddEntry( h_cosmu_unf, " Unfolded ", "pl" );

    double cmu_tot = 0;
    
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
            h_cosmu_unf->SetBinContent( j, v_unf[0]->GetBinContent(i, j) );
        }

        h_cosmu->SetTitle(hist_name1);
        h_cosmu->GetXaxis()->SetTitle("T_{#mu} (GeV)");   
        h_cosmu->GetYaxis()->SetTitle("Number of SBND events");   
        h_cosmu->SetLineColor( kRed + 2 );
        h_cosmu_sm->SetLineColor( kGreen + 2 );
        h_cosmu->Draw();
        h_cosmu_sm->Draw("same");
        h_cosmu_unf->Draw("E1 same");
        h_cosmu_unf->SetLineWidth(1);
        h_cosmu_unf->SetLineColor(kBlack);
        h_cosmu_unf->SetMarkerStyle(20);
        h_cosmu_unf->SetMarkerSize(0.6);

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
void SliceStack ( std::vector<TH2*> v_un, std::vector<TH2*> v_sm, std::vector<TH2*> v_sm_rec, std::vector<TH2*> v_unf){

    TCanvas *canvas1 = new TCanvas ( "canvas1", "", 800, 600 );

    // Firstly, loop over all Tmu bins on draw slices in cos theta mu
    // Take the bin edges to be the title
    int x_bins = v_un[0]->GetNbinsX(); // Cos theta
    int y_bins = v_un[0]->GetNbinsY(); // Tmu

    string unsmeared_labs[5] = {"CCQE","CCMEC","CCRES","CCDIS","CCCOH"};
    string smeared_labs[6] = {"NC n#pi","CCQE","CCMEC","CCRES","CCDIS","CCCOH"};
    string reco_labs[12] = {"NC 1#pi","NC Other","CC1#pi^{-}","CC1#pi^{+}","CC1#pi^{0}","CC2#pi^{#pm}","CC Other","#nu_{#mu} CC0#pi, 0p","#nu_{#mu} CC0#pi, 1p","#nu_{#mu} CC0#pi, 2p","#nu_{#mu} CC0#pi, 3p","#nu_{#mu} CC0#pi, >3p"};

    for ( int j = 1; j < x_bins+1; j++ ){

      StackProjectionY(v_un, unsmeared_labs, "un_Tmu_cosmuRange:", j, canvas1, v_unf[0], 1);
      StackProjectionY(v_sm, smeared_labs, "sm_Tmu_cosmuRange:", j, canvas1, v_unf[0], 0);
      StackProjectionY(v_sm_rec, reco_labs, "rec_sm_Tmu_cosmuRange:", j, canvas1, v_unf[0], 0);

    }

    for ( int j = 1; j < y_bins+1; j++ ){

      StackProjectionX(v_un, unsmeared_labs, "un_cosmu_TmuRange:", j, canvas1, v_unf[0], 1);
      StackProjectionX(v_sm, smeared_labs, "sm_cosmu_TmuRange:", j, canvas1, v_unf[0], 0);
      StackProjectionX(v_sm_rec, reco_labs, "rec_sm_cosmu_TmuRange:", j, canvas1, v_unf[0], 0);

    }

    delete canvas1;
    
}
