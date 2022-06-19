// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file HFLbToLcPiCandidateSelector.cxx
/// \brief Λb0 → Λc+ π- candidate selector
///
/// \author Mohammad Waris

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/Core/HFSelectorCuts.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "PWGHF/Core/HFSelectorCuts.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand_lb;
using namespace o2::analysis;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_cand_prong4;
using namespace o2::analysis::hf_cuts_xicc_tolcpikpi;

struct HfLbToLcPiCandidateSelector {
  Produces<aod::HFSelXiccToLcPiKPiCandidate> hFSelXiccToLcPiKPiCandidate;
  
  Configurable<double> pTCandMin{"pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> pTCandMax{"pTCandMax", 36., "Upper bound of candidate pT"};
  Configurable<bool> d_FilterPID{"d_FilterPID", true, "Bool to use or not the PID at filtering level"};
  //Track quality
  Configurable<double> TPCNClsFindablePIDCut{"TPCNClsFindablePIDCut", 70., "Lower bound of TPC findable clusters for good PID"};

  //TPC PID
  Configurable<double> pidTPCMinpT{"pidTPCMinpT", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> pidTPCMaxpT{"pidTPCMaxpT", 10., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTPC{"nSigmaTPC", 5., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTPCCombined{"nSigmaTPCCombined", 5., "Nsigma cut on TPC combined with TOF"};

  //TOF PID
  Configurable<double> pidTOFMinpT{"pidTOFMinpT", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> pidTOFMaxpT{"pidTOFMaxpT", 10., "Upper bound of track pT for TOF PID"};
  Configurable<double> nSigmaTOF{"nSigmaTOF", 5., "Nsigma cut on TOF only"};
  Configurable<double> nSigmaTOFCombined{"nSigmaTOFCombined", 5., "Nsigma cut on TOF combined with TPC"};
   
  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_xicc_tolcpikpi::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"Xicc_to_lcpikpi_cuts", {hf_cuts_xicc_tolcpikpi::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "Xicc candidate selection per pT bin"};
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc+"};
 
  
  // Apply topological cuts as defined in HFSelectorCuts.h; return true if candidate passes all cuts
  template <typename T1, typename T2, typename T3,typename T4 ,typename T5 >
  bool selectionTopol(const T1& hfCandXicc, const T2& hfCandLc, const T3& trackPi1, const T4& trackKa, const T5& trackPi2 )
  { 
    auto candpT = hfCandXicc.pt();
    int pTBin = findBin(pTBins, candpT);
    if (pTBin == -1) {
      // Printf("Xicc topol selection failed at getpTBin");
      return false;
    }
    
    //check that the candidate pT is within the analysis range
    //if (candpT < pTCandMin || candpT >= pTCandMax) {
    //  return false;
    //}
    
    
    // check candidate mass is within a defined mass window
    if (std::abs(InvMassXiccToLcPiKPi(hfCandXicc) - RecoDecay::getMassPDG(pdg::Code::kXiCCPlusPlus)) > cuts->get(pTBin, "m")) {
      return false;
    }

    //Xicc CPA cut
    if (hfCandXicc.cpa() < cuts->get(pTBin, "CPA")) {
      return false;
    }
    //Xicc CPA XY cut
    if (hfCandXicc.cpaXY() <= cuts->get(pTBin, "CPA XY")) {
      return false;
    }

 

    // candidate maximum decay length
    if (hfCandXicc.decayLength() > cuts->get(pTBin, "Max Xicc decLen")) {
      return false;
    }
    // candidate minimum decay length
    if (hfCandXicc.decayLength() < cuts->get(pTBin, "Min Xicc decLen")) {
      return false;
    }

 
    //Xicc Decay length XY
    if (hfCandXicc.decayLengthXY() < cuts->get(pTBin, "Min Xicc decLenXY")) {
      return false;
    }

    

    //Xicc chi2PCA cut
    if (hfCandXicc.chi2PCA() > cuts->get(pTBin, "Chi2PCA")) {
      //Printf("Lb selection failed at chi2PCA");
      return false;
    }



 
    //pt cuts on daughters

    //Pion  and Kaon pt cuts
    if (trackPi1.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    if (trackPi2.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    if (trackKa.pt() < cuts->get(pTBin, "pT Pi")) {
      return false;
    }
    //Lc pt cuts
    if (hfCandLc.pt() < cuts->get(pTBin, "pT Lc")) {
      return false;
    }


    /*
    //d0 
    if (std::abs(hfCandXicc.impactParameter1()) < cuts->get(pTBin, "d0 Pi")) {
      return false;
    }

     if (std::abs(hfCandXicc.impactParameter2()) < cuts->get(pTBin, "d0 Pi")) {
      return false;
    }

     if (std::abs(hfCandXicc.impactParameter3()) < cuts->get(pTBin, "d0 Pi")) {
      return false;
    }

    //d0 of Lc+
    if (std::abs(hfCandXicc.impactParameter0()) < cuts->get(pTBin, "d0 Lc")) {
      return false;
    }

    // impact parameter product
    if (hfCandXicc.impactParameterProduct() > cuts->get(pTBin, "d0d0")) {
      return false;
    }
    
    
     
    //Lc mass
    //if (trackPi.sign() < 0) {
    //if (std::abs(InvMassLcpKpi(hfCandLc) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus)) > cuts->get(pTBin, "DeltaMLc")) {
    //return false;
    //}
    //}

    //Xicc Decay length
    if (hfCandXicc.decayLength() < cuts->get(pTBin, "Xicc decLen")) {
      return false;
    }

    // candidate minimum decay length
    if (hfCandXicc.decayLength() <= cuts->get(pTBin, "min decay length")) {
      return false;
    }
    // candidate maximum decay length XY
    if (hfCandXicc.decayLengthXY() > cuts->get(pTBin, "max decay length XY")) {
      return false;
    }

    // candidate minimum decay length XY
    //if (hfCandXicc.decayLengthXY() < cuts->get(pTBin, "min decay length XY")) {
    //  return false;
    //}

    // minimum DCA of daughters
    if ((std::abs(hfCandXicc.impactParameter0()) <= cuts->get(pTBin, "min d0 Xic")) ||
        (std::abs(hfCandXicc.impactParameter1()) <= cuts->get(pTBin, "min d0 Pi"))) {
      return false;
    }

    // maximum DCA of daughters
    if ((std::abs(hfCandXicc.impactParameter0()) > cuts->get(pTBin, "max d0 Xic")) ||
        (std::abs(hfCandXicc.impactParameter1()) > cuts->get(pTBin, "max d0 Pi"))) {
      return false;
    } */

    return true;
  }

  void process(aod::HfCandXicctoLcPiKPi const& hfCandXicctoLcPiKPis, soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate>, aod::BigTracksPID const&)
  {
    for (auto& hfCandXicctoLcPiKPi  : hfCandXicctoLcPiKPis) { //looping over Xicc candidates

      int statusXicc = 0;

      // check if flagged as Λb --> Λc+ π-
      if (!(hfCandXicctoLcPiKPi.hfflag() & 1 << hf_cand_XicctoLcPiKPi::DecayType::XicctoLcPiKPi)) {
        hFSelXiccToLcPiKPiCandidate(statusXicc);
        //Printf("Xicc candidate selection failed at hfflag check");
        continue;
      }

      // Lc is always index0 and pi is index1 by default
      //auto candLc = hfCandLb.index0();
      auto hfCandLc = hfCandXicctoLcPiKPi.index0_as<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate>>();
      auto trackPi1 = hfCandXicctoLcPiKPi.index1_as<aod::BigTracksPID>();
      auto trackKa = hfCandXicctoLcPiKPi.index2_as<aod::BigTracksPID>();
      auto trackPi2 = hfCandXicctoLcPiKPi.index3_as<aod::BigTracksPID>();
      
      //topological cuts
      if (!selectionTopol(hfCandXicctoLcPiKPi, hfCandLc, trackPi1, trackKa, trackPi2)) {
      //if (!selectionTopol( hfCandLc, trackPi1, trackKa, trackPi2)) {
        hFSelXiccToLcPiKPiCandidate(statusXicc);
        // Printf("Xicc candidate selection failed at selection topology");
        continue;
      } 

      hFSelXiccToLcPiKPiCandidate(1);
      //Printf("Xicc candidate selection successful, candidate should be selected");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfLbToLcPiCandidateSelector>(cfgc));
  return workflow;
}