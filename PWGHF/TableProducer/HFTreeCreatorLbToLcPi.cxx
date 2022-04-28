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

/// \file HFTreeCreatorLbToLcPi.cxx
/// \brief Writer of the 2 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
/// \note Extended from HFTreeCreatorD0ToKPi, HFTreeCreatorLcToPKPi, HFTreeCreatorXToJpsiPiPi
///
/// \author Panos Christakoglou <Panos.Christakoglou@cern.ch>, Nikhef

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/TrackSelectorPID.h"
#include "ALICE3/DataModel/RICH.h"
#include "ALICE3/DataModel/MID.h"
#include "Common/Core/PID/PIDResponse.h" //from D0 cand sel ALICE3 Barrel
#include "ReconstructionDataFormats/PID.h" //from D0 cand sel ALICE3 Barrel

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_XicctoLcPiKPi;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(PtProng2, ptProng2, float);
DECLARE_SOA_COLUMN(PProng2, pProng2, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised2, impactParameterNormalised2, float);
DECLARE_SOA_COLUMN(PtProng3, ptProng3, float);
DECLARE_SOA_COLUMN(PProng3, pProng3, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised3, impactParameterNormalised3, float);
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(CPA, cpa, float);
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(MCflag, mcflag, int8_t);
// DECLARE_SOA_COLUMN(NSigTOFPi0, nsigTOFPi0, float);
////Lc selection parameter
DECLARE_SOA_COLUMN(LcM, lcM, float);
DECLARE_SOA_COLUMN(LcCt, lcCt, float);
DECLARE_SOA_COLUMN(LcY, lcY, float);
DECLARE_SOA_COLUMN(LcE, lcE, float);
DECLARE_SOA_COLUMN(LcEta, lcEta, float);
DECLARE_SOA_COLUMN(LcCPA, lcCPA, float);
DECLARE_SOA_COLUMN(LcCPAXY, lcCPAXY, float);
DECLARE_SOA_COLUMN(LcChi2PCA, lcChi2PCA, float);
DECLARE_SOA_COLUMN(LcDecayLength, lcDecayLength, float);
DECLARE_SOA_COLUMN(LcDecayLengthXY, lcDecayLengthXY, float);
DECLARE_SOA_COLUMN(LcDecayLengthNormalised, lcDecayLengthNormalised, float);
/*DECLARE_SOA_COLUMN(NSigTOFTrk1Pi, nSigTOFTrk1Pi, float);
DECLARE_SOA_COLUMN(NSigTOFTrk2Ka, nSigTOFTrk2Ka, float);
DECLARE_SOA_COLUMN(NSigTOFTrk3Pi, nSigTOFTrk3Pi, float);
DECLARE_SOA_COLUMN(NSigTOFTrk4Pi, nSigTOFTrk4Pi, float);
DECLARE_SOA_COLUMN(NSigTOFTrk4Pr, nSigTOFTrk4Pr, float);
DECLARE_SOA_COLUMN(NSigTOFTrk5Ka, nSigTOFTrk5Ka, float);
DECLARE_SOA_COLUMN(NSigTOFTrk6Pi, nSigTOFTrk6Pi, float);
DECLARE_SOA_COLUMN(NSigTOFTrk6Pr, nSigTOFTrk6Pr, float); */
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

// put the arguments into the table
DECLARE_SOA_TABLE(HfXicc4Full, "AOD", "HFXicc4Full",
                  full::RSecondaryVertex,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  hf_cand::Chi2PCA,
                  full::ImpactParameterNormalised0,
                  full::PtProng0,
                  full::PProng0,
                  full::ImpactParameterNormalised1,
                  full::PtProng1,
                  full::PProng1,
                  full::ImpactParameterNormalised2,
                  full::PtProng2,
                  full::PProng2,
                  full::ImpactParameterNormalised3,
                  full::PtProng3,
                  full::PProng3,
                  hf_cand::PxProng0,
                  hf_cand::PyProng0,
                  hf_cand::PzProng0,
                  hf_cand::PxProng1,
                  hf_cand::PyProng1,
                  hf_cand::PzProng1,
                  hf_cand::PxProng2,
                  hf_cand::PyProng2,
                  hf_cand::PzProng2,
                  hf_cand::PxProng3,
                  hf_cand::PyProng3,
                  hf_cand::PzProng3,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ImpactParameter2,
                  hf_cand::ImpactParameter3,
                  hf_cand::ErrorImpactParameter0,
                  hf_cand::ErrorImpactParameter1,
                  hf_cand::ErrorImpactParameter2,
                  hf_cand::ErrorImpactParameter3,
                  /* full::NSigTOFTrk1Pi,
                  full::NSigTOFTrk2Ka,
                  full::NSigTOFTrk3Pi,
                  full::NSigTOFTrk4Pi,
                  full::NSigTOFTrk4Pr,
                  full::NSigTOFTrk5Ka,
                  full::NSigTOFTrk6Pi,
                  full::NSigTOFTrk6Pr, */
                  full::LcM,
                  full::LcCt,
                  full::LcY,
                  full::LcE,
                  full::LcEta,
                  full::LcCPA,
                  full::LcCPAXY,
                  full::LcChi2PCA,
                  full::LcDecayLength,
                  full::LcDecayLengthXY,
                  full::LcDecayLengthNormalised,
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::P,
                  full::CPA,
                  full::CPAXY,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::MCflag);

DECLARE_SOA_TABLE(HfXicc4FullEvents, "AOD", "HFXicc4FullE",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);
//
DECLARE_SOA_TABLE(HfXicc4FullParticles, "AOD", "HFXicc4FullP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::MCflag);

} // namespace o2::aod


/// Writes the full information in an output TTree
struct HfTreeCreatorXiccToLcpikpi {
  Produces<o2::aod::HfXicc4Full> rowCandidateFull;
  Produces<o2::aod::HfXicc4FullEvents> rowCandidateFullEvents;  // massive data tables keeping track of all particles. Not needed for ML
  Produces<o2::aod::HfXicc4FullParticles> rowCandidateFullParticles; // // massive data tables keeping track of all particles. Not needed for ML

  void init(InitContext const&)
  {
  }

  void process(aod::Collisions const& collisions,
               aod::McCollisions const& mccollisions,
               soa::Join<aod::HfCandXicctoLcPiKPi, aod::HfCandXicctoLcPiKPiMCRec, aod::HFSelXiccToLcPiKPiCandidate> const& candidates,
               soa::Join<aod::HfCandProng3, aod::HfCandProng3MCRec, aod::HFSelLcCandidate> const& Lccandidates,
               soa::Join<aod::McParticles, aod::HfCandXicctoLcPiKPiMCGen> const& particles,
               // aod::BigTracksPID const& tracks,
               // aod::BigTracksMC const& bigtracksmc)
               // , ExtendedTracksPID const&,
               //aod::RICHs const&)
  // aod::MIDs const&)
  {

    // Filling event properties
    rowCandidateFullEvents.reserve(collisions.size());
    for (auto& collision : collisions) {
      rowCandidateFullEvents(
        collision.bcId(),
        collision.numContrib(),
        collision.posX(),
        collision.posY(),
        collision.posZ(),
        0,
        1);
    }

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto& candidate : candidates) {
      auto fillTable = [&](int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY) {
        if (FunctionSelection >= 1) { // Set to true to keep unselected events as well  FunctionSelection >= 1
          auto LcCand = candidate.index0_as<soa::Join<aod::HfCandProng3, aod::HfCandProng3MCRec, aod::HFSelLcCandidate>>();
          /* auto track0 = candidate.index1_as<ExtendedTracksPID>(); //daughter pion track
          auto track1 = candidate.index2_as<ExtendedTracksPID>(); //daughter Kaon track
          auto track2 = candidate.index3_as<ExtendedTracksPID>(); //daughter pion track
          auto track3 = LcCand.index0_as<ExtendedTracksPID>(); //granddaughter tracks (lc decay particles)
          auto track4 = LcCand.index1_as<ExtendedTracksPID>();
          auto track5 = LcCand.index2_as<ExtendedTracksPID>(); */

          rowCandidateFull(
            candidate.rSecondaryVertex(),
            candidate.decayLength(),
            candidate.decayLengthXY(),
            candidate.decayLengthNormalised(),
            candidate.decayLengthXYNormalised(),
            candidate.chi2PCA(),
            candidate.impactParameterNormalised0(),
            candidate.ptProng0(),
            RecoDecay::P(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
            candidate.impactParameterNormalised1(),
            candidate.ptProng1(),
            RecoDecay::P(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
            candidate.impactParameterNormalised2(),
            candidate.ptProng2(),
            RecoDecay::P(candidate.pxProng2(), candidate.pyProng2(), candidate.pzProng2()),
            candidate.impactParameterNormalised3(),
            candidate.ptProng3(),
            RecoDecay::P(candidate.pxProng3(), candidate.pyProng3(), candidate.pzProng3()),
            candidate.pxProng0(),
            candidate.pyProng0(),
            candidate.pzProng0(),
            candidate.pxProng1(),
            candidate.pyProng1(),
            candidate.pzProng1(),
            candidate.pxProng2(),
            candidate.pyProng2(),
            candidate.pzProng2(),
            candidate.pxProng3(),
            candidate.pyProng3(),
            candidate.pzProng3(),
            candidate.impactParameter0(),
            candidate.impactParameter1(),
            candidate.impactParameter2(),
            candidate.impactParameter3(),
            candidate.errorImpactParameter0(),
            candidate.errorImpactParameter1(),
            candidate.errorImpactParameter2(),
            candidate.errorImpactParameter3(),
            /*candidate.index1_as<aod::BigTracksPID>().tofNSigmaPi(),
            candidate.index2_as<aod::BigTracksPID>().tofNSigmaKa(),
            candidate.index3_as<aod::BigTracksPID>().tofNSigmaPi(),
            LcCand.index0_as<aod::BigTracksPID>().tofNSigmaPi(),
            LcCand.index0_as<aod::BigTracksPID>().tofNSigmaPr(),
            LcCand.index1_as<aod::BigTracksPID>().tofNSigmaKa(),
            LcCand.index2_as<aod::BigTracksPID>().tofNSigmaPi(),
            LcCand.index2_as<aod::BigTracksPID>().tofNSigmaPr(), */
            o2::aod::hf_cand_prong3::InvMassLcpKpi(LcCand),
            o2::aod::hf_cand_prong3::CtLc(LcCand),
            o2::aod::hf_cand_prong3::YLc(LcCand),
            o2::aod::hf_cand_prong3::ELc(LcCand),
            LcCand.eta(),
            LcCand.cpa(),
            LcCand.cpaXY(),
            LcCand.chi2PCA(),
            LcCand.decayLength(),
            LcCand.decayLengthXY(),
            LcCand.decayLengthXYNormalised(),
            FunctionSelection,
            FunctionInvMass,
            candidate.pt(),
            candidate.p(),
            candidate.cpa(),
            candidate.cpaXY(),
            FunctionCt,
            candidate.eta(),
            candidate.phi(),
            FunctionY,
            candidate.flagMCMatchRec());
        }
      };

      fillTable(candidate.isSelXiccToLcPiKPi(), InvMassXiccToLcPiKPi(candidate), CtXicc(candidate), YXicc(candidate)); 
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (auto& particle : particles) {
      if (std::abs(particle.flagMCMatchGen()) == 1 << DecayType::XicctoLcPiKPi) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          particle.eta(),
          particle.phi(),
          RecoDecay::Y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode())),
          particle.flagMCMatchGen());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorXiccToLcpikpi>(cfgc));
  return workflow;
}