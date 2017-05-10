/*
 *--------------------------------------------------------------
 *
 * Author: Rhiannon Jones
 * Date  : February 2017
 *
 *
*/

#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "THStack.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <Math/GSLRndmEngines.h>
#include "TLatex.h"
#include "TStyle.h"
#include "TPad.h"
#include "TColor.h"
#include "TObjArray.h"
#include "TGaxis.h"
#include "TError.h"

#include "TRandom.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfold.h"

// -------------------------------------------------------------------------
// Smearing function:
// For the 2D histogram, smear the cos theta bins by 5 deg and the Tmu bins
// by 10%
// -------------------------------------------------------------------------

void Smear ( TTree *tree,
             TTree *traintree,
             TH2D *h_smeared,
             std::vector<TH2*> &v_un,
             std::vector<TH2*> &v_sm,
             std::vector<TH2*> &v_sm_rec,
             std::vector<TH2*> &v_unf); 


// -------------------------------------------------------------------------
// Slices function:
// For the 2D histogram, draw 1D histograms slice-by-slice ( bin-by-bin )
// -------------------------------------------------------------------------

void Slices ( TH2D *h_unsmeared, 
              TH2D *h_smeared,
              std::vector<TH2*> v_unf);

void SliceStack ( std::vector<TH2*> v_un,
                  std::vector<TH2*> v_sm,
                  std::vector<TH2*> v_sm_rec,
                  std::vector<TH2*> v_unf);
