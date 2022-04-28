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
/// \author Panos Christakoglou <panos.christakoglou@cern.ch>, Nikhef

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
//using namespace o2::analysis::hf_cuts_lb_tolcpi;

struct HfLbToLcPiCandidateSelector {
  Produces<aod::HFSelXiccToLcPiKPiCandidate> hFSelXiccToLcPiKPiCandidate;

/*
  Configurable<double> pTCandMin{"pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> pTCandMax{"pTCandMax", 50., "Upper bound of candidate pT"};

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

  Configurable<std::vector<double>> pTBins{"pTBins", std::vector<double>{hf_cuts_lb_tolcpi::pTBins_v}, "pT bin limits"};
  Configurable<LabeledArray<double>> cuts{"Lb_to_lcpi_cuts", {hf_cuts_lb_tolcpi::cuts[0], npTBins, nCutVars, pTBinLabels, cutVarLabels}, "Lb0 candidate selection per pT bin"};
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc+"}; */

  // Apply topological cuts as defined in HFSelectorCuts.h; return true if candidate passes all cuts
  template <typename T1, typename T2, typename T3>
  bool selectionTopol(const T1& hfCandLb, const T2& hfCandLc, const T3& trackPi)
  {


 
  

    return true;
  }

  void process(aod::HfCandXicctoLcPiKPi const& hfCandXicctoLcPiKPis, soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate> )//, aod::BigTracksPID const&)
  {
    for (auto& hfCandXicctoLcPiKPi : hfCandXicctoLcPiKPis) { //looping over Lb candidates

      /*int statusLb = 0;

      // check if flagged as Λb --> Λc+ π-
      if (!(hfCandLb.hfflag() & 1 << hf_cand_lb::DecayType::LbToLcPi)) {
        hfSelLbToLcPiCandidate(statusLb);
        //Printf("Lb candidate selection failed at hfflag check");
        continue;
      }

      // Lc is always index0 and pi is index1 by default
      //auto candLc = hfCandLb.index0();
      auto candLc = hfCandLb.index0_as<soa::Join<aod::HfCandProng3, aod::HFSelLcCandidate>>();
      auto trackPi = hfCandLb.index1_as<aod::BigTracksPID>();

      //topological cuts
      if (!selectionTopol(hfCandLb, candLc, trackPi)) {
        hfSelLbToLcPiCandidate(statusLb);
        // Printf("Lb candidate selection failed at selection topology");
        continue;
      } */

      hFSelXiccToLcPiKPiCandidate(1);
      //Printf("Lb candidate selection successful, candidate should be selected");
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfLbToLcPiCandidateSelector>(cfgc));
  return workflow;
}
