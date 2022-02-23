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

/// \file HFCandidateCreatorLb.cxx
/// \brief Reconstruction of Lb candidates
/// \note Adapted from HFCandidateCreatorXicc
///
/// \author Panos Christakoglou <panos.christakoglou@cern.ch>, Nikhef





// importing packages - functionally - header files - libraries

#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/V0.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"

//  namespace - used to label a function specifying which package it belongs
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_cand_prong3;
using namespace o2::aod::hf_cand_lb;
using namespace o2::framework::expressions;

// function below - void means that it doesnt return anything - within brackets are the arguments of the function
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
// function itself  - what it does - vector <o2::framework::ConfigParamSpec> is a template - workflowoptions is an input parameter of the format
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Perform MC matching."}};
  workflowOptions.push_back(optionDoMC);
}
 
// this function is enabling the MC matching that happens at near the end of the code
#include "Framework/runDataProcessing.h"

/// Reconstruction of Λb candidates
struct HFCandidateCreatorLb {
  /// this thing below: where we tell o2 that what comes out of the structure. its going to be different in our case.
  Produces<aod::HfCandLbBase> rowCandidateBase;
  /// strange way of assinging values to the variables, looks more like a dictionary       
  Configurable<double> magneticField{"magneticField", 20., "magnetic field"};
  Configurable<bool> b_propdca{"b_propdca", true, "create tracks version propagated to PCA"};
  /// maybe greater d_max than 200
  Configurable<double> d_maxr{"d_maxr", 200., "reject PCA's above this radius"};
  // no clue
  Configurable<double> d_maxdzini{"d_maxdzini", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  // what do you mean by change: when reconstructing you use a fitting function iteratively. we stop at some point
  Configurable<double> d_minparamchange{"d_minparamchange", 1.e-3, "stop iterations if largest change of any Lb is smaller than this"};
  //should remain the same
  Configurable<double> d_minrelchi2change{"d_minrelchi2change", 0.9, "stop iterations is chi2/chi2old > this"};
  //threshold is actually a candidiate selection parameter, check this.  Ones below are ignored, they will not be detected in the  real experiment.
  Configurable<double> ptPionMin{"ptPionMin", 0.5, "minimum pion pT threshold (GeV/c)"};



/// new particles require new entries in the histograms
  OutputObj<TH1F> hMassLcToPKPi{TH1F("hMassLcToPKPi", "#Lambda_{c}^{#plus} candidates;inv. mass (pK^{#minus} #pi^{#plus}) (GeV/#it{c}^{2});entries", 500, 0., 5.)};
  // histogram of mass of lc to pkpi
  OutputObj<TH1F> hPtLc{TH1F("hPtLc", "#Lambda_{c}^{#plus} candidates;#Lambda_{c}^{#plus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hPtPion{TH1F("hPtPion", "#pi^{#minus} candidates;#pi^{#minus} candidate #it{p}_{T} (GeV/#it{c});entries", 100, 0., 10.)};
  OutputObj<TH1F> hCPALc{TH1F("hCPALc", "#Lambda_{c}^{#plus} candidates;#Lambda_{c}^{#plus} cosine of pointing angle;entries", 110, -1.1, 1.1)};
  OutputObj<TH1F> hMassLbToLcPi{TH1F("hMassLbToLcPi", "2-prong candidates;inv. mass (#Lambda_{b}^{0} #rightarrow #Lambda_{c}^{#plus}#pi^{#minus} #rightarrow pK^{#minus}#pi^{#plus}#pi^{#minus}) (GeV/#it{c}^{2});entries", 500, 3., 8.)};
  // this needs to be changed to 4 prong vertex
  OutputObj<TH1F> hCovPVXX{TH1F("hCovPVXX", "2-prong candidates;XX element of cov. matrix of prim. vtx position (cm^{2});entries", 100, 0., 1.e-4)};
  // cov  =covariance matrix  (dont worry much about it)
  OutputObj<TH1F> hCovSVXX{TH1F("hCovSVXX", "2-prong candidates;XX element of cov. matrix of sec. vtx position (cm^{2});entries", 100, 0., 0.2)};


/// normal way of assigning value to variable
  double massPi = RecoDecay::getMassPDG(kPiMinus);
  double massLc = RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus);
  
/// this the actually the value to be calculated of lmbda b invariant masses
  double massLcPi = 0.;

  Configurable<int> d_selectionFlagLc{"d_selectionFlagLc", 1, "Selection Flag for Lc"};
/// selection flag =  applying flags to suitables candidates for Lc
/// why not apply to pions = they are simple to detect?
///is kaon a simple particle = also a simple particle doesnt req flagging?
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  /// This doesnt depend on the nature of particles, accepts everything
  Filter filterSelectCandidates = (aod::hf_selcandidate_lc::isSelLcpKpi >= d_selectionFlagLc || aod::hf_selcandidate_lc::isSelLcpiKp >= d_selectionFlagLc);
  /// assigning a variable again?. Filter is keeping only the flagged ones




  // this function does: 


  void process(aod::Collision const& collision,
               soa::Filtered<soa::Join<
                 aod::HfCandProng3,
                 aod::HFSelLcCandidate>> const& lcCands,
               aod::BigTracks const& tracks)   // check this format? what is this?
  {
    // 2-prong vertex fitter
    o2::vertexing::DCAFitterN<2> df2;
    df2.setBz(magneticField);
    df2.setPropagateToPCA(b_propdca);
    df2.setMaxR(d_maxr);
    df2.setMaxDZIni(d_maxdzini);
    df2.setMinParamChange(d_minparamchange);
    df2.setMinRelChi2Change(d_minrelchi2change);
    df2.setUseAbsDCA(true);

    // 3-prong vertex fitter (to rebuild Lc vertex)
    o2::vertexing::DCAFitterN<3> df3;
    df3.setBz(magneticField);
    df3.setPropagateToPCA(b_propdca);
    df3.setMaxR(d_maxr);
    df3.setMaxDZIni(d_maxdzini);
    df3.setMinParamChange(d_minparamchange);
    df3.setMinRelChi2Change(d_minrelchi2change);
    df3.setUseAbsDCA(true);


    /// my prongs
    o2::vertexing::DCAFitterN<4> df4;
    df4.setBz(magneticField);
    df4.setPropagateToPCA(b_propdca);
    df4.setMaxR(d_maxr);
    df4.setMaxDZIni(d_maxdzini);
    df4.setMinParamChange(d_minparamchange);
    df4.setMinRelChi2Change(d_minrelchi2change);
    df4.setUseAbsDCA(true);


    /// check the difference bw vertex fitter and reconstructor



    // loop over Lc candidates




    // auto& if you dont know its an integer or double
    //lc Cands is the list of candidates
    // iterators : pointer to first element of the container (something advanced than an array, its a collection of elements). Difference in array and vector in C++
    // lcCand  =  pointer to elements in container in lcCands.


    for (auto& lcCand : lcCands) {
      /// check if it has been labelled as a lc decaying into PKPi
      
      if (!(lcCand.hfflag() & 1 << o2::aod::hf_cand_prong3::DecayType::LcToPKPi)) {    /// check what this does.
        continue;
      }
      if (lcCand.isSelLcpKpi() >= d_selectionFlagLc) {         /// proton and pion are both positive. Kaon is negative. Taking two combaination togather + - + (three paths)
        hMassLcToPKPi->Fill(InvMassLcpKpi(lcCand), lcCand.pt());
      }
      if (lcCand.isSelLcpiKp() >= d_selectionFlagLc) {
        hMassLcToPKPi->Fill(InvMassLcpiKp(lcCand), lcCand.pt());
      }
      // adding values to histograms
      hPtLc->Fill(lcCand.pt());
      hCPALc->Fill(lcCand.cpa());
      /// declaring the tracks of daughters of lambda C
      auto track0 = lcCand.index0_as<aod::BigTracks>();
      auto track1 = lcCand.index1_as<aod::BigTracks>();
      auto track2 = lcCand.index2_as<aod::BigTracks>();
      auto trackParVar0 = getTrackParCov(track0);  /// 
      auto trackParVar1 = getTrackParCov(track1);
      auto trackParVar2 = getTrackParCov(track2);
      // this is collision for the tracks, all tracks result from same collision, so its same
      auto collision = track0.collision();

      // reconstruct the 3-prong secondary vertex

      // decay vertex is farther from original vertex * 
      if (df3.process(trackParVar0, trackParVar1, trackParVar2) == 0) {  // ?
        continue;
      }
      const auto& secondaryVertex = df3.getPCACandidate();
      trackParVar0.propagateTo(secondaryVertex[0], magneticField);
      trackParVar1.propagateTo(secondaryVertex[0], magneticField);
      trackParVar2.propagateTo(secondaryVertex[0], magneticField);

      array<float, 3> pvecpK = {track0.px() + track1.px(), track0.py() + track1.py(), track0.pz() + track1.pz()}; // computing the momentum vector by summing up momenta of daughters
      array<float, 3> pvecLc = {pvecpK[0] + track2.px(), pvecpK[1] + track2.py(), pvecpK[2] + track2.pz()};
      auto trackpK = o2::dataformats::V0(df3.getPCACandidatePos(), pvecpK, df3.calcPCACovMatrixFlat(),
                                         trackParVar0, trackParVar1, {0, 0}, {0, 0});
      auto trackLc = o2::dataformats::V0(df3.getPCACandidatePos(), pvecLc, df3.calcPCACovMatrixFlat(),
                                         trackpK, trackParVar2, {0, 0}, {0, 0});

      int index0Lc = track0.globalIndex();
      int index1Lc = track1.globalIndex();
      int index2Lc = track2.globalIndex();
      //int charge = track0.sign() + track1.sign() + track2.sign();


      /// look for two pions because i have two pions to look for aginst one that is generated by lambda c decay
      ///

      
      // pions = []  // initiate a dictionary/ container for trackpions
      // 
      for (auto& trackPion : tracks) {
        if (trackPion.pt() < ptPionMin) {    /// momenta being selection criteria? a safe slection criteria
          continue;                          /// ptPionMin change this for Xi cc /// limitations :  combinations and storage // maybe do more cuts (more selection criteria)
        }
        // check for my own decay
        if (trackPion.sign() > 0) {         // in my decay decay we have has positive pions /// change to trackPion.sign() < 0
          continue;
        }
        /// checks if the new pion is not the one coming Lc
        if (trackPion.globalIndex() == index0Lc || trackPion.globalIndex() == index1Lc || trackPion.globalIndex() == index2Lc) {
          continue;
        }
        // pions.append[trackPion]
      }
      
      //pairpion = []  // initiate a dictionary/ container for trackpions pairs
      for (auto& trackPion1 : range(len(pions))) {

        for (auto& trackPion2 : range(len(pions))) {
          if (trackPion1 == trackPion2){
            continue  
          }
          pairpion.append([pions[trackpion1], pions[trackpion2]])
        }
        
      }
    
        // get a list of filtered pions. loop over that list 

        for (auto& trackKaon : tracks){
          if (trackKaon.pt()< ptKaonMin){
            continue;
          }
          if(trackKaon.sign()>0){
            continue;
          }
          if (trackKaon.globalIndex() == index0Lc || trackPion.globalIndex() == index1Lc || trackPion.globalIndex() == index2Lc) {
          continue;
        }

        }




        //filling to histogram
        hPtPion->Fill(trackPion.pt());
        array<float, 3> pvecPion;
        auto trackParVarPi = getTrackParCov(trackPion);

        // reconstruct the 2-prong Lb vertex
        if (df2.process(trackLc, trackParVarPi) == 0) {  /// equals 0 check is if the process function was successfull
          continue;
        }

        // calculate relevant properties     (of potential lb candidate)
        const auto& secondaryVertexLb = df2.getPCACandidate();
        auto chi2PCA = df2.getChi2AtPCACandidate();
        auto covMatrixPCA = df2.calcPCACovMatrixFlat();
        // replace with df3 ?  df4????
        df2.propagateTracksToVertex();
        df2.getTrack(0).getPxPyPzGlo(pvecLc);
        df2.getTrack(1).getPxPyPzGlo(pvecPion);

        auto primaryVertex = getPrimaryVertex(collision);
        auto covMatrixPV = primaryVertex.getCov();
        o2::dataformats::DCA impactParameter0;
        o2::dataformats::DCA impactParameter1;
        trackLc.propagateToDCA(primaryVertex, magneticField, &impactParameter0);
        trackParVarPi.propagateToDCA(primaryVertex, magneticField, &impactParameter1);

        hCovSVXX->Fill(covMatrixPCA[0]);
        hCovPVXX->Fill(covMatrixPV[0]);

        // get uncertainty of the decay length
        double phi, theta;
        getPointDirection(array{collision.posX(), collision.posY(), collision.posZ()}, secondaryVertexLb, phi, theta);
        auto errorDecayLength = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixPCA, phi, theta));
        auto errorDecayLengthXY = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.) + getRotatedCovMatrixXX(covMatrixPCA, phi, 0.));

        int hfFlag = 1 << hf_cand_lb::DecayType::LbToLcPi;

        // fill the candidate table for the Lb here:
        rowCandidateBase(collision.globalIndex(),
                         collision.posX(), collision.posY(), collision.posZ(),
                         secondaryVertexLb[0], secondaryVertexLb[1], secondaryVertexLb[2],
                         errorDecayLength, errorDecayLengthXY,
                         chi2PCA,
                         pvecLc[0], pvecLc[1], pvecLc[2],
                         pvecPion[0], pvecPion[1], pvecPion[2],
                         impactParameter0.getY(), impactParameter1.getY(),
                         std::sqrt(impactParameter0.getSigmaY2()), std::sqrt(impactParameter1.getSigmaY2()),
                         lcCand.globalIndex(), trackPion.globalIndex(),
                         hfFlag);

        // calculate invariant mass
        auto arrayMomenta = array{pvecLc, pvecPion};
        massLcPi = RecoDecay::M(std::move(arrayMomenta), array{massLc, massPi});
        if (lcCand.isSelLcpKpi() > 0) {
          hMassLbToLcPi->Fill(massLcPi);
        }
        if (lcCand.isSelLcpiKp() > 0) {
          hMassLbToLcPi->Fill(massLcPi);
        }
      } // pi- loop
    }   // Lc loop
  }     // process
};      // struct

/// Extends the base table with expression columns.
struct HFCandidateCreatorLbExpressions {
  Spawns<aod::HfCandLbExt> rowCandidateLb;
  void init(InitContext const&) {}
};

/// Performs MC matching.
struct HFCandidateCreatorLbMC {
  Produces<aod::HfCandLbMCRec> rowMCMatchRec;
  Produces<aod::HfCandLbMCGen> rowMCMatchGen;

  void process(aod::HfCandLb const& candidates,
               aod::HfCandProng3,
               aod::BigTracksMC const& tracks,
               aod::McParticles_000 const& particlesMC)
  {
    int indexRec = -1;
    int8_t sign = 0;
    int8_t flag = 0;
    int8_t origin = 0;
    int8_t debug = 0;

    // Match reconstructed candidates.
    for (auto& candidate : candidates) {
      //Printf("New rec. candidate");
      flag = 0;
      origin = 0;
      debug = 0;
      auto lcCand = candidate.index0();
      auto arrayDaughters = array{lcCand.index0_as<aod::BigTracksMC>(),
                                  lcCand.index1_as<aod::BigTracksMC>(),
                                  lcCand.index2_as<aod::BigTracksMC>(),
                                  candidate.index1_as<aod::BigTracksMC>()};
      auto arrayDaughtersLc = array{lcCand.index0_as<aod::BigTracksMC>(),
                                    lcCand.index1_as<aod::BigTracksMC>(),
                                    lcCand.index2_as<aod::BigTracksMC>()};
      // Λb → Λc+ π-
      //Printf("Checking Λb → Λc+ π-");
      indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kLambdaB0, array{+kProton, -kKPlus, +kPiPlus, -kPiPlus}, true, &sign, 2);
      if (indexRec > -1) {
        // Λb → Λc+ π-
        //Printf("Checking Λb → Λc+ π-");
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughtersLc, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 1);
        if (indexRec > -1) {
          flag = 1 << hf_cand_lb::DecayType::LbToLcPi;
        } else {
          debug = 1;
          LOGF(info, "WARNING: Λb in decays in the expected final state but the condition on the intermediate state is not fulfilled");
        }
      }
      rowMCMatchRec(flag, origin, debug);
    }

    // Match generated particles.
    for (auto& particle : particlesMC) {
      //Printf("New gen. candidate");
      flag = 0;
      origin = 0;
      // Λb → Λc+ π-
      if (RecoDecay::isMatchedMCGen(particlesMC, particle, pdg::Code::kLambdaB0, array{int(pdg::Code::kLambdaCPlus), -kPiPlus}, true)) {
        // Match Λc+ -> pKπ
        auto LcCandMC = particlesMC.iteratorAt(particle.daughter0Id());
        //Printf("Checking Λc+ → p K- π+");
        if (RecoDecay::isMatchedMCGen(particlesMC, LcCandMC, int(pdg::Code::kLambdaCPlus), array{+kProton, -kKPlus, +kPiPlus}, true, &sign)) {
          flag = sign * (1 << hf_cand_lb::DecayType::LbToLcPi);
        }
      }
      rowMCMatchGen(flag, origin);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HFCandidateCreatorLb>(cfgc, TaskName{"hf-cand-creator-lb"}),
    adaptAnalysisTask<HFCandidateCreatorLbExpressions>(cfgc, TaskName{"hf-cand-creator-lb-expressions"})};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HFCandidateCreatorLbMC>(cfgc, TaskName{"hf-cand-creator-lb-mc"}));
  }
  return workflow;
}




// Particle identification part needs to be added/ actually not, see: PWGHF/TableProducer/HFLcK0sPCandidateSelector.cxx
// no pid for pion, pid for kaons? = 