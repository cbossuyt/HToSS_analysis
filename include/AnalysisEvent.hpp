//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 23 14:04:00 2014 by ROOT version 5.32/00
// from TChain tree/
//////////////////////////////////////////////////////////

#ifndef _AnalysisEvent_hpp_
#define _AnalysisEvent_hpp_

#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <iostream>
#include <string>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class AnalysisEvent
{
    public:
    TTree* fChain; //! pointer to the analyzed TTree or TChain
    Int_t fCurrent; //! current Tree number in a TChain

    // Declaration of leaf types
    static const size_t NELECTRONSMAX{30};
    Int_t numElePF2PAT;
    Float_t elePF2PATE[NELECTRONSMAX];
    Float_t elePF2PATET[NELECTRONSMAX];
    Float_t elePF2PATPX[NELECTRONSMAX];
    Float_t elePF2PATPY[NELECTRONSMAX];
    Float_t elePF2PATPZ[NELECTRONSMAX];
    Float_t elePF2PATPhi[NELECTRONSMAX];
    Float_t elePF2PATTheta[NELECTRONSMAX];
    Float_t elePF2PATEta[NELECTRONSMAX];
    Float_t elePF2PATPT[NELECTRONSMAX];
    Int_t elePF2PATCharge[NELECTRONSMAX];
    Int_t elePF2PATCutIdVeto[NELECTRONSMAX];
    Int_t elePF2PATCutIdLoose[NELECTRONSMAX];
    Int_t elePF2PATCutIdMedium[NELECTRONSMAX];
    Int_t elePF2PATCutIdTight[NELECTRONSMAX];
    Int_t elePF2PATImpactTransDist[NELECTRONSMAX];
    Float_t elePF2PATImpactTransError[NELECTRONSMAX];
    Float_t elePF2PATImpactTransSignificance[NELECTRONSMAX];
    Float_t elePF2PATImpact3DDist[NELECTRONSMAX];
    Float_t elePF2PATImpact3DError[NELECTRONSMAX];
    Float_t elePF2PATImpact3DSignificance[NELECTRONSMAX];
    Float_t elePF2PATChargedHadronIso[NELECTRONSMAX];
    Float_t elePF2PATNeutralHadronIso[NELECTRONSMAX];
    Float_t elePF2PATPhotonIso[NELECTRONSMAX];
    Float_t elePF2PATTrackPt[NELECTRONSMAX];
    Float_t elePF2PATTrackPhi[NELECTRONSMAX];
    Float_t elePF2PATTrackEta[NELECTRONSMAX];
    Float_t elePF2PATTrackChi2[NELECTRONSMAX];
    Float_t elePF2PATTrackNDOF[NELECTRONSMAX];
    Float_t elePF2PATTrackD0[NELECTRONSMAX];
    Float_t elePF2PATTrackDBD0[NELECTRONSMAX];
    Float_t elePF2PATD0PV[NELECTRONSMAX];
    Float_t elePF2PATDZPV[NELECTRONSMAX];
    Float_t elePF2PATBeamSpotCorrectedTrackD0[NELECTRONSMAX];
    Float_t elePF2PATTrackDz[NELECTRONSMAX];
    Float_t elePF2PATVtxZ[NELECTRONSMAX];
    Int_t elePF2PATIsGsf[NELECTRONSMAX];
    Float_t elePF2PATGsfPx[NELECTRONSMAX];
    Float_t elePF2PATGsfPy[NELECTRONSMAX];
    Float_t elePF2PATGsfPz[NELECTRONSMAX];
    Float_t elePF2PATGsfE[NELECTRONSMAX];
    Float_t elePF2PATEcalEnergy[NELECTRONSMAX];
    Float_t elePF2PATSCEta[NELECTRONSMAX];
    Float_t elePF2PATSCE[NELECTRONSMAX];
    Float_t elePF2PATSCPhi[NELECTRONSMAX];
    Float_t elePF2PATSCEoverP[NELECTRONSMAX];
    Float_t elePF2PATSCSigmaEtaEta[NELECTRONSMAX];
    Float_t elePF2PATSCSigmaIEtaIEta[NELECTRONSMAX];
    Float_t elePF2PATSCSigmaIEtaIEta5x5[NELECTRONSMAX];
    Float_t elePF2PATSCE1x5[NELECTRONSMAX];
    Float_t elePF2PATSCE5x5[NELECTRONSMAX];
    Float_t elePF2PATSCE2x5max[NELECTRONSMAX];
    Float_t elePF2PATTrackIso04[NELECTRONSMAX];
    Float_t elePF2PATEcalIso04[NELECTRONSMAX];
    Float_t elePF2PATHcalIso04[NELECTRONSMAX];
    Float_t elePF2PATTrackIso03[NELECTRONSMAX];
    Float_t elePF2PATEcalIso03[NELECTRONSMAX];
    Float_t elePF2PATHcalIso03[NELECTRONSMAX];
    Float_t elePF2PATdr04EcalRecHitSumEt[NELECTRONSMAX];
    Float_t elePF2PATdr03EcalRecHitSumEt[NELECTRONSMAX];
    Float_t elePF2PATEcalIsoDeposit[NELECTRONSMAX];
    Float_t elePF2PATHcalIsoDeposit[NELECTRONSMAX];
    Float_t elePF2PATComRelIso[NELECTRONSMAX];
    Float_t elePF2PATComRelIsodBeta[NELECTRONSMAX];
    Float_t elePF2PATComRelIsoRho[NELECTRONSMAX];
    Float_t elePF2PATChHadIso[NELECTRONSMAX];
    Float_t elePF2PATNtHadIso[NELECTRONSMAX];
    Float_t elePF2PATGammaIso[NELECTRONSMAX];
    Float_t elePF2PATRhoIso[NELECTRONSMAX];
    Float_t elePF2PATAEff03[NELECTRONSMAX];
    Int_t elePF2PATMissingInnerLayers[NELECTRONSMAX];
    Float_t elePF2PATHoverE[NELECTRONSMAX];
    Float_t elePF2PATDeltaPhiSC[NELECTRONSMAX];
    Float_t elePF2PATDeltaEtaSC[NELECTRONSMAX];
    Float_t elePF2PATDeltaEtaSeedSC[NELECTRONSMAX];
    Int_t elePF2PATIsBarrel[NELECTRONSMAX];
    Int_t elePF2PATPhotonConversionTag[NELECTRONSMAX];
    Float_t elePF2PATPhotonConversionDist[NELECTRONSMAX];
    Float_t elePF2PATPhotonConversionDcot[NELECTRONSMAX];
    Int_t elePF2PATPhotonConversionVeto[NELECTRONSMAX];
//    Int_t elePF2PATPhotonConversionTagCustom[NELECTRONSMAX];
    Float_t elePF2PATPhotonConversionDistCustom[NELECTRONSMAX];
    Float_t elePF2PATPhotonConversionDcotCustom[NELECTRONSMAX];
    Float_t elePF2PATTriggerMatch[NELECTRONSMAX];
    Float_t elePF2PATJetOverlap[NELECTRONSMAX];
    Float_t genElePF2PATPT[NELECTRONSMAX];
    Float_t genElePF2PATET[NELECTRONSMAX];
    Float_t genElePF2PATPX[NELECTRONSMAX];
    Float_t genElePF2PATPY[NELECTRONSMAX];
    Float_t genElePF2PATPZ[NELECTRONSMAX];
    Float_t genElePF2PATPhi[NELECTRONSMAX];
    Float_t genElePF2PATTheta[NELECTRONSMAX];
    Float_t genElePF2PATEta[NELECTRONSMAX];
    Int_t genElePF2PATCharge[NELECTRONSMAX];
    Int_t genElePF2PATPdgId[NELECTRONSMAX];
    Int_t genElePF2PATMotherId[NELECTRONSMAX];
    Int_t genElePF2PATPromptDecayed[NELECTRONSMAX];
    Int_t genElePF2PATPromptFinalState[NELECTRONSMAX];
    Int_t genElePF2PATHardProcess[NELECTRONSMAX];

    static constexpr size_t NMUONSMAX{20};
    Int_t numMuonPF2PAT;
    Float_t muonPF2PATE[NMUONSMAX];
    Float_t muonPF2PATET[NMUONSMAX];
    Float_t muonPF2PATPt[NMUONSMAX];
    Float_t muonPF2PATPX[NMUONSMAX];
    Float_t muonPF2PATPY[NMUONSMAX];
    Float_t muonPF2PATPZ[NMUONSMAX];
    Float_t muonPF2PATPhi[NMUONSMAX];
    Float_t muonPF2PATTheta[NMUONSMAX];
    Float_t muonPF2PATEta[NMUONSMAX];
    Int_t muonPF2PATCharge[NMUONSMAX];
    Int_t muonPF2PATLooseCutId[NMUONSMAX];
    Int_t muonPF2PATMediumCutId[NMUONSMAX];
    Int_t muonPF2PATTightCutId[NMUONSMAX];
    Int_t muonPF2PATPfIsoVeryLoose[NMUONSMAX];
    Int_t muonPF2PATPfIsoLoose[NMUONSMAX];
    Int_t muonPF2PATPfIsoMedium[NMUONSMAX];
    Int_t muonPF2PATPfIsoTight[NMUONSMAX];
    Int_t muonPF2PATPfIsoVeryTight[NMUONSMAX];
    Int_t muonPF2PATTkIsoLoose[NMUONSMAX];
    Int_t muonPF2PATTkIsoTight[NMUONSMAX];
    Int_t muonPF2PATMvaLoose[NMUONSMAX];
    Int_t muonPF2PATMvaMedium[NMUONSMAX];
    Int_t muonPF2PATMvaTight[NMUONSMAX];
    Float_t muonPF2PATGlobalID[NMUONSMAX];
    Float_t muonPF2PATTrackID[NMUONSMAX];
    Float_t muonPF2PATChi2[NMUONSMAX];
    Float_t muonPF2PATD0[NMUONSMAX];
    Float_t muonPF2PATTrackDBD0[NMUONSMAX];
    Float_t muonPF2PATDBInnerTrackD0[NMUONSMAX];
    Float_t muonPF2PATBeamSpotCorrectedD0[NMUONSMAX];
    Int_t muonPF2PATTrackNHits[NMUONSMAX];
    Int_t muonPF2PATMuonNHits[NMUONSMAX];
    Float_t muonPF2PATNDOF[NMUONSMAX];
    Float_t muonPF2PATVertX[NMUONSMAX];
    Float_t muonPF2PATVertY[NMUONSMAX];
    Float_t muonPF2PATVertZ[NMUONSMAX];
    Float_t muonPF2PATChargedHadronIso[NMUONSMAX];
    Float_t muonPF2PATNeutralHadronIso[NMUONSMAX];
    Float_t muonPF2PATPhotonIso[NMUONSMAX];
    Float_t muonPF2PATTrackIso[NMUONSMAX];
    Float_t muonPF2PATEcalIso[NMUONSMAX];
    Float_t muonPF2PATHcalIso[NMUONSMAX];
    Float_t muonPF2PATComRelIso[NMUONSMAX];
    Float_t muonPF2PATComRelIsodBeta[NMUONSMAX];
    Int_t muonPF2PATIsPFMuon[NMUONSMAX];
    Int_t muonPF2PATNChambers[NMUONSMAX];
    Int_t muonPF2PATNMatches[NMUONSMAX];
    Int_t muonPF2PATTkLysWithMeasurements[NMUONSMAX];
    Int_t muonPF2PATVldPixHits[NMUONSMAX];
    Int_t muonPF2PATMatchedStations[NMUONSMAX];
    Float_t muonPF2PATGlbTkNormChi2[NMUONSMAX];
    Float_t muonPF2PATValidFraction[NMUONSMAX];
    Float_t muonPF2PATChi2LocalPosition[NMUONSMAX];
    Float_t muonPF2PATTrkKick[NMUONSMAX];
    Float_t muonPF2PATSegmentCompatibility[NMUONSMAX];
    Float_t muonPF2PATDBPV[NMUONSMAX];
    Float_t muonPF2PATDZPV[NMUONSMAX];
    Float_t genMuonPF2PATPT[NMUONSMAX];
    Float_t genMuonPF2PATET[NMUONSMAX];
    Float_t genMuonPF2PATPX[NMUONSMAX];
    Float_t genMuonPF2PATPY[NMUONSMAX];
    Float_t genMuonPF2PATPZ[NMUONSMAX];
    Float_t genMuonPF2PATPhi[NMUONSMAX];
    Float_t genMuonPF2PATTheta[NMUONSMAX];
    Float_t genMuonPF2PATEta[NMUONSMAX];
    Int_t genMuonPF2PATCharge[NMUONSMAX];
    Int_t genMuonPF2PATPdgId[NMUONSMAX];
    Int_t genMuonPF2PATMotherId[NMUONSMAX];
    Int_t genMuonPF2PATPromptDecayed[NMUONSMAX];
    Int_t genMuonPF2PATPromptFinalState[NMUONSMAX];
    Int_t genMuonPF2PATHardProcess[NMUONSMAX];

    static constexpr size_t NJETSMAX{40};
    Int_t numJetPF2PAT;
    Double_t jetPF2PATE[NJETSMAX];
    Double_t jetPF2PATEt[NJETSMAX];
    Double_t jetPF2PATPt[NJETSMAX];
    Double_t jetPF2PATPtRaw[NJETSMAX];
    Double_t jetPF2PATUnCorEt[NJETSMAX];
    Double_t jetPF2PATUnCorPt[NJETSMAX];
    Double_t jetPF2PATEta[NJETSMAX];
    Double_t jetPF2PATTheta[NJETSMAX];
    Double_t jetPF2PATPhi[NJETSMAX];
    Double_t jetPF2PATPx[NJETSMAX];
    Double_t jetPF2PATPy[NJETSMAX];
    Double_t jetPF2PATPz[NJETSMAX];
    Double_t jetPF2PATdRClosestLepton[NJETSMAX];
    Int_t jetPF2PATNtracksInJet[NJETSMAX];
    Float_t jetPF2PATJetCharge[NJETSMAX];
    Float_t jetPF2PATfHPD[NJETSMAX];
    Float_t jetPF2PATBtagSoftMuonPtRel[NJETSMAX];
    Float_t jetPF2PATBtagSoftMuonQuality[NJETSMAX];
    Float_t jetPF2PATCorrFactor[NJETSMAX];
    Float_t jetPF2PATCorrResidual[NJETSMAX];
    Float_t jetPF2PATL2L3ResErr[NJETSMAX];
    Float_t jetPF2PATCorrErrLow[NJETSMAX];
    Float_t jetPF2PATCorrErrHi[NJETSMAX];
    Float_t jetPF2PATN90Hits[NJETSMAX];
    Float_t jetPF2PATTriggered[NJETSMAX];
    Float_t jetPF2PATSVX[NJETSMAX];
    Float_t jetPF2PATSVY[NJETSMAX];
    Float_t jetPF2PATSVZ[NJETSMAX];
    Float_t jetPF2PATSVDX[NJETSMAX];
    Float_t jetPF2PATSVDY[NJETSMAX];
    Float_t jetPF2PATSVDZ[NJETSMAX];
    Float_t jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags[NJETSMAX];
    Float_t jetPF2PATpfCombinedCvsLJetTags[NJETSMAX];
    Float_t jetPF2PATpfCombinedCvsBJetTags[NJETSMAX];
    Int_t jetPF2PATNConstituents[NJETSMAX];
    Int_t jetPF2PATPID[NJETSMAX];
    Float_t jetPF2PATClosestBPartonDeltaR[NJETSMAX];
    Float_t jetPF2PATClosestCPartonDeltaR[NJETSMAX];
    Float_t genJetPF2PATE[NJETSMAX];
    Float_t genJetPF2PATET[NJETSMAX];
    Float_t genJetPF2PATPT[NJETSMAX];
    Float_t genJetPF2PATPX[NJETSMAX];
    Float_t genJetPF2PATPY[NJETSMAX];
    Float_t genJetPF2PATPZ[NJETSMAX];
    Float_t genJetPF2PATMass[NJETSMAX];
    Float_t genJetPF2PATPhi[NJETSMAX];
    Float_t genJetPF2PATTheta[NJETSMAX];
    Float_t genJetPF2PATEta[NJETSMAX];
    Int_t genJetPF2PATPID[NJETSMAX];
    Int_t genJetPF2PATMotherPID[NJETSMAX];
    Int_t genJetPF2PATScalarAncestor[NJETSMAX];
    Float_t jetPF2PATMuEnergy[NJETSMAX];
    Float_t jetPF2PATMuEnergyFraction[NJETSMAX];
    Float_t jetPF2PATNeutralHadEnergy[NJETSMAX];
    Float_t jetPF2PATNeutralEmEnergy[NJETSMAX];
    Float_t jetPF2PATChargedHadronEnergyFraction[NJETSMAX];
    Float_t jetPF2PATNeutralHadronEnergyFraction[NJETSMAX];
    Float_t jetPF2PATChargedEmEnergyFraction[NJETSMAX];
    Float_t jetPF2PATNeutralEmEnergyFraction[NJETSMAX];
    Float_t jetPF2PATMuonFraction[NJETSMAX];
    Float_t jetPF2PATChargedHadronEnergyFractionCorr[NJETSMAX];
    Float_t jetPF2PATNeutralHadronEnergyFractionCorr[NJETSMAX];
    Float_t jetPF2PATChargedEmEnergyFractionCorr[NJETSMAX];
    Float_t jetPF2PATNeutralEmEnergyFractionCorr[NJETSMAX];
    Float_t jetPF2PATMuonFractionCorr[NJETSMAX];
    Int_t jetPF2PATNeutralMultiplicity[NJETSMAX];
    Int_t jetPF2PATChargedMultiplicity[NJETSMAX];
    Float_t fixedGridRhoFastjetAll;

    Double_t metPF2PATE;
    Double_t metPF2PATEt;
    Double_t metPF2PATEtRaw;
    Double_t metPF2PATPhi;
    Double_t metPF2PATPt;
    Double_t metPF2PATPx;
    Double_t metPF2PATPy;
    Double_t metPF2PATPz;
    Float_t metPF2PATScalarEt;
    Float_t metPF2PATEtUncorrected;
    Float_t metPF2PATPhiUncorrected;
    Float_t metPF2PATUnclusteredEnUp;
    Float_t metPF2PATUnclusteredEnDown;
    Float_t genMetPF2PATE;
    Float_t genMetPF2PATEt;
    Float_t genMetPF2PATPhi;
    Float_t genMetPF2PATPt;
    Float_t genMetPF2PATPx;
    Float_t genMetPF2PATPy;
    Float_t genMetPF2PATPz;

    static constexpr size_t NTAUSMAX{20};
    Int_t numTauPF2PAT;
    Float_t tauPF2PATE[NTAUSMAX];
    Float_t tauPF2PATPt[NTAUSMAX];
    Float_t tauPF2PATPhi[NTAUSMAX];
    Float_t tauPF2PATEta[NTAUSMAX];

    static constexpr size_t NPHOTONSMAX{30};
    Int_t numPhoPF2PAT;
    Float_t phoPF2PATE[NPHOTONSMAX];
    Float_t phoPF2PATSigmaE[NPHOTONSMAX];
    Float_t phoPF2PATET[NPHOTONSMAX];
    Float_t phoPF2PATPhi[NPHOTONSMAX];
    Float_t phoPF2PATEta[NPHOTONSMAX];
    Float_t phoPF2PATTheta[NPHOTONSMAX];
    Float_t phoPF2PATPt[NPHOTONSMAX];
    Float_t phoPF2PATPx[NPHOTONSMAX];
    Float_t phoPF2PATPy[NPHOTONSMAX];
    Float_t phoPF2PATPz[NPHOTONSMAX];
    Float_t phoPF2PATCalibE[NPHOTONSMAX];
    Float_t phoPF2PATCalibEt[NPHOTONSMAX];
    Float_t phoPF2PATSCE[NPHOTONSMAX];
    Float_t phoPF2PATSCRawE[NPHOTONSMAX];
    Float_t phoPF2PATESEnP1[NPHOTONSMAX];
    Float_t phoPF2PATESEnP2[NPHOTONSMAX];
    Float_t phoPF2PATSCEta[NPHOTONSMAX];
    Float_t phoPF2PATSCEtaWidth[NPHOTONSMAX];
    Float_t phoPF2PATSCPhi[NPHOTONSMAX];
    Float_t phoPF2PATSCPhiWidth[NPHOTONSMAX];
    Float_t phoPF2PATSCBrem[NPHOTONSMAX];
    Int_t phoPF2PATHasPixelSeed[NPHOTONSMAX];
    Int_t phoPF2PATEleVeto[NPHOTONSMAX];
    Float_t phoPF2PATR9[NPHOTONSMAX];
    Float_t phoPF2PATHoverE[NPHOTONSMAX];
    Float_t phoPF2PATESEffSigmaRR[NPHOTONSMAX];
    Float_t phoPF2PATSigmaIEtaIEtaFull5x5[NPHOTONSMAX];
    Float_t phoPF2PATSigmaIEtaIPhiFull5x5[NPHOTONSMAX];
    Float_t phoPF2PATSigmaIPhiIPhiFull5x5[NPHOTONSMAX];
    Float_t phoPF2PATE2x2Full5x5[NPHOTONSMAX];
    Float_t phoPF2PATE5x5Full5x5[NPHOTONSMAX];
    Float_t phoPF2PATR9Full5x5[NPHOTONSMAX];
    Float_t phoPF2PATPFChIso[NPHOTONSMAX];
    Float_t phoPF2PATPFPhoIso [NPHOTONSMAX];
    Float_t phoPF2PATPFNeuIso[NPHOTONSMAX];
    Float_t phoPF2PATPFChWorstIso[NPHOTONSMAX];
    Float_t phoPF2PATMIPTotEnergy[NPHOTONSMAX];
    Int_t phoPF2PATCutIdLoose[NPHOTONSMAX];
    Int_t phoPF2PATCutIdMedium[NPHOTONSMAX];
    Int_t phoPF2PATCutIdTight[NPHOTONSMAX];
    Int_t phoPF2PATMvaIdWp80[NPHOTONSMAX];
    Int_t phoPF2PATMvaIdWp90[NPHOTONSMAX];
    Float_t genPhoPF2PATPt[NPHOTONSMAX];
    Float_t genPhoPF2PATET[NPHOTONSMAX];
    Float_t genPhoPF2PATEta[NPHOTONSMAX];
    Float_t genPhoPF2PATTheta[NPHOTONSMAX];
    Float_t genPhoPF2PATPhi[NPHOTONSMAX];
    Float_t genPhoPF2PATPx[NPHOTONSMAX];
    Float_t genPhoPF2PATPy[NPHOTONSMAX];
    Float_t genPhoPF2PATPz[NPHOTONSMAX];
    Int_t genPhoPF2PATCharge[NPHOTONSMAX];
    Int_t genPhoPF2PATPdgId[NPHOTONSMAX];
    Int_t genPhoPF2PATMotherId[NPHOTONSMAX];
    Int_t genPhoPF2PATIsPhoton[NPHOTONSMAX];
    Int_t genPhoPF2PATIsConvertedPhoton[NPHOTONSMAX];
    Int_t genPhoPF2PATIsJet[NPHOTONSMAX];
/*
    static constexpr size_t NOOTPHOTONSMAX{30};
    Int_t numPhoOOT_PF2PAT;
    Float_t phoOOT_PF2PATE[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATSigmaE[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATET[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATPhi[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATEta[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATTheta[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATPt[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATPx[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATPy[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATPz[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATCalibE[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATCalibEt[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATSCE[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATSCRawE[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATESEnP1[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATESEnP2[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATSCEta[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATSCEtaWidth[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATSCPhi[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATSCPhiWidth[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATSCBrem[NOOTPHOTONSMAX];
    Int_t phoOOT_PF2PATHasPixelSeed[NOOTPHOTONSMAX];
    Int_t phoOOT_PF2PATEleVeto[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATR9[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATHoverE[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATESEffSigmaRR[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATSigmaIEtaIEtaFull5x5[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATSigmaIEtaIPhiFull5x5[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATSigmaIPhiIPhiFull5x5[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATE2x2Full5x5[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATE5x5Full5x5[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATR9Full5x5[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATPFChIso[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATPFPhoIso [NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATPFNeuIso[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATPFChWorstIso[NOOTPHOTONSMAX];
    Float_t phoOOT_PF2PATMIPTotEnergy[NOOTPHOTONSMAX];
    Int_t phoOOT_PF2PATCutIdLoose[NOOTPHOTONSMAX];
    Int_t phoOOT_PF2PATCutIdMedium[NOOTPHOTONSMAX];
    Int_t phoOOT_PF2PATCutIdTight[NOOTPHOTONSMAX];
    Int_t phoOOT_PF2PATMvaIdWp80[NOOTPHOTONSMAX];
    Int_t phoOOT_PF2PATMvaIdWp90[NOOTPHOTONSMAX];
    Float_t genPhoOOT_PF2PATPt[NOOTPHOTONSMAX];
    Float_t genPhoOOT_PF2PATET[NOOTPHOTONSMAX];
    Float_t genPhoOOT_PF2PATEta[NOOTPHOTONSMAX];
    Float_t genPhoOOT_PF2PATTheta[NOOTPHOTONSMAX];
    Float_t genPhoOOT_PF2PATPhi[NOOTPHOTONSMAX];
    Float_t genPhoOOT_PF2PATPx[NOOTPHOTONSMAX];
    Float_t genPhoOOT_PF2PATPy[NOOTPHOTONSMAX];
    Float_t genPhoOOT_PF2PATPz[NOOTPHOTONSMAX];
    Int_t genPhoOOT_PF2PATCharge[NOOTPHOTONSMAX];
    Int_t genPhoOOT_PF2PATPdgId[NOOTPHOTONSMAX];
    Int_t genPhoOOT_PF2PATMotherId[NOOTPHOTONSMAX];
    Int_t genPhoOOT_PF2PATIsPhoton[NOOTPHOTONSMAX];
    Int_t genPhoOOT_PF2PATIsConvertedPhoton[NOOTPHOTONSMAX];
    Int_t genPhoOOT_PF2PATIsJet[NOOTPHOTONSMAX];
*/
    static constexpr size_t NPACKEDCANDSMAX{700};
    Int_t numPackedCands;
//    Float_t packedCandsPt[NPACKEDCANDSMAX];
    Float_t packedCandsPx[NPACKEDCANDSMAX];
    Float_t packedCandsPy[NPACKEDCANDSMAX];
    Float_t packedCandsPz[NPACKEDCANDSMAX];
    Float_t packedCandsE[NPACKEDCANDSMAX];
//    Float_t packedCandsEta[NPACKEDCANDSMAX];
//    Float_t packedCandsTheta[NPACKEDCANDSMAX];
//    Float_t packedCandsPhi[NPACKEDCANDSMAX];
    Int_t packedCandsCharge[NPACKEDCANDSMAX];
    Int_t packedCandsPdgId[NPACKEDCANDSMAX];
//    Float_t packedCandsTime[NPACKEDCANDSMAX];
//    Int_t packedCandsFromPV[NPACKEDCANDSMAX];
//    Int_t packedCandsPVquality[NPACKEDCANDSMAX];
//    Float_t packedCandsVx[NPACKEDCANDSMAX];
//    Float_t packedCandsVy[NPACKEDCANDSMAX];
//    Float_t packedCandsVz[NPACKEDCANDSMAX];
//    Float_t packedCandsVEta[NPACKEDCANDSMAX];
//    Float_t packedCandsVPhi[NPACKEDCANDSMAX];
//    Float_t packedCandsBeamSpotCorrectedD0[NPACKEDCANDSMAX];
    Float_t packedCandsDz[NPACKEDCANDSMAX];
    Float_t packedCandsDxy[NPACKEDCANDSMAX];
//    Float_t packedCandsDzAssocPV[NPACKEDCANDSMAX];
//    Float_t packedCandsVtxChi2Norm[NPACKEDCANDSMAX];
    Int_t packedCandsHasTrackDetails[NPACKEDCANDSMAX];
//    Float_t packedCandsDzError[NPACKEDCANDSMAX];
//    Float_t packedCandsDxyError[NPACKEDCANDSMAX];
//    Float_t packedCandsTimeError[NPACKEDCANDSMAX];
//    Float_t packedCandsPseudoTrkPt[NPACKEDCANDSMAX];
    Float_t packedCandsPseudoTrkPx[NPACKEDCANDSMAX];
    Float_t packedCandsPseudoTrkPy[NPACKEDCANDSMAX];
    Float_t packedCandsPseudoTrkPz[NPACKEDCANDSMAX];
//    Float_t packedCandsPseudoTrkEta[NPACKEDCANDSMAX];
//    Float_t packedCandsPseudoTrkPhi[NPACKEDCANDSMAX];
    Int_t packedCandsPseudoTrkCharge[NPACKEDCANDSMAX];
    Float_t packedCandsPseudoTrkVx[NPACKEDCANDSMAX];
    Float_t packedCandsPseudoTrkVy[NPACKEDCANDSMAX];
    Float_t packedCandsPseudoTrkVz[NPACKEDCANDSMAX];
//    Float_t packedCandsPseudoTrkChi2Norm[NPACKEDCANDSMAX];
//    Int_t packedCandsPseudoTrkNumberOfHits[NPACKEDCANDSMAX];
//    Int_t packedCandsPseudoTrkNumberOfPixelHits[NPACKEDCANDSMAX];
//    Int_t packedCandsPseudoTrkPixelLayersWithMeasurement[NPACKEDCANDSMAX];
//    Int_t packedCandsPseudoTrkStripLayersWithMeasurement[NPACKEDCANDSMAX];
//    Int_t packedCandsPseudoTrkTrackerLayersWithMeasurement[NPACKEDCANDSMAX];
    Int_t packedCandsHighPurityTrack[NPACKEDCANDSMAX];

    static constexpr size_t NISOTRACKSMAX{40};
    Int_t numIsolatedTracks;
    Float_t isoTracksPt[NISOTRACKSMAX];
    Float_t isoTracksPx[NISOTRACKSMAX];
    Float_t isoTracksPy[NISOTRACKSMAX];
    Float_t isoTracksPz[NISOTRACKSMAX];
    Float_t isoTracksE[NISOTRACKSMAX];
    Float_t isoTracksEta[NISOTRACKSMAX];
    Float_t isoTracksTheta[NISOTRACKSMAX];
    Float_t isoTracksPhi[NISOTRACKSMAX];
    Int_t isoTracksCharge[NISOTRACKSMAX];
    Int_t isoTracksPdgId[NISOTRACKSMAX];
    Float_t isoTracksMatchedCaloJetEmEnergy[NISOTRACKSMAX];
    Float_t isoTracksMatchedCaloJetHadEnergy[NISOTRACKSMAX];
    Float_t isoTracksDz[NISOTRACKSMAX];
    Float_t isoTracksDxy[NISOTRACKSMAX];
    Float_t isoTracksDzError[NISOTRACKSMAX];
    Float_t isoTracksDxyError[NISOTRACKSMAX];
    Int_t isoTracksFromPV[NISOTRACKSMAX];
    Float_t isoTracksVx[NISOTRACKSMAX];
    Float_t isoTracksVy[NISOTRACKSMAX];
    Float_t isoTracksVz[NISOTRACKSMAX];
    Int_t isoTracksHighPurity[NISOTRACKSMAX];
    Int_t isoTracksTight[NISOTRACKSMAX];
    Int_t isoTracksLoose[NISOTRACKSMAX];
    Float_t isoTracksDeltaEta[NISOTRACKSMAX];
    Float_t isoTracksDeltaPhi[NISOTRACKSMAX];

    Int_t isElePlusJets;
    Float_t genPDFScale;
    Float_t genPDFx1;
    Float_t genPDFx2;
    Int_t genPDFf1;
    Int_t genPDFf2;
    Double_t topPtReweight;
    Int_t processId;
    Float_t processPtHat;
    Double_t processMCWeight;
    Float_t beamSpotX;
    Float_t beamSpotY;
    Float_t beamSpotZ;

    static constexpr size_t NPVSMAX{80};
    Int_t numPVs;
    Float_t pvX[NPVSMAX];
    Float_t pvY[NPVSMAX];
    Float_t pvZ[NPVSMAX];
    Float_t pvDX[NPVSMAX];
    Float_t pvDY[NPVSMAX];
    Float_t pvDZ[NPVSMAX];
    Float_t pvRho[NPVSMAX];
    Int_t pvIsFake[NPVSMAX];
    Float_t pvNdof[NPVSMAX];
    Float_t pvChi2[NPVSMAX];
    Float_t pvNtracks[NPVSMAX];
    Float_t pvNtracksW05[NPVSMAX];

    static constexpr size_t NSVSMAX{20};
    Int_t numSVs;
    Float_t svPt[NSVSMAX];
    Float_t svPx[NSVSMAX];
    Float_t svPy[NSVSMAX];
    Float_t svPz[NSVSMAX];
    Float_t svMass[NSVSMAX];
    Float_t svE[NSVSMAX];
    Float_t svEta[NSVSMAX];
    Float_t svTheta[NSVSMAX];
    Float_t svPhi[NSVSMAX];
    Float_t svX[NSVSMAX];
    Float_t svY[NSVSMAX];
    Float_t svZ[NSVSMAX];
    Float_t svVertexChi2[NSVSMAX];
    Float_t svVertexNdof[NSVSMAX];
    Int_t svNtracks[NSVSMAX];
    Float_t svDist3D[NSVSMAX];
    Float_t svDist3DSig[NSVSMAX];
    Float_t svDist3DError[NSVSMAX];
    Float_t svDistXY[NSVSMAX];
    Float_t svDistXYSig[NSVSMAX];
    Float_t svDistXYError[NSVSMAX];
    Float_t svAnglePV[NSVSMAX];
    Int_t svIsLambda[NSVSMAX];
    Int_t svIsKshort[NSVSMAX];

    Float_t mhtPt;
    Float_t mhtPy;
    Float_t mhtPx;
    Float_t mhtPhi;
    Float_t mhtSumEt;
    Float_t mhtSignif;

    static constexpr size_t NTRIGGERBITSMAX{700};
    Int_t nTriggerBits;
    Int_t TriggerBits[NTRIGGERBITSMAX];

    Double_t weight_muF0p5;
    Double_t weight_muF2;
    Double_t weight_muR0p5;
    Double_t weight_muR2;
    Double_t weight_muF0p5muR0p5;
    Double_t weight_muF2muR2;
    Double_t origWeightForNorm;
    Double_t weight_pdfMax;
    Double_t weight_pdfMin;
    Double_t weight_alphaMax;
    Double_t weight_alphaMin;
    Double_t isrRedHi{};
    Double_t fsrRedHi{};
    Double_t isrRedLo{};
    Double_t fsrRedLo{};
    Double_t isrDefHi{};
    Double_t fsrDefHi{};
    Double_t isrDefLo{};
    Double_t fsrDefLo{};
    Double_t isrConHi{};
    Double_t fsrConHi{};
    Double_t isrConLo{};
    Double_t fsrConLo{};
    //   Int_t           numVert;

    // 2016 lepton Triggers
/*
    Int_t HLT_Ele25_eta2p1_WPTight_Gsf_v1;
    Int_t HLT_Ele25_eta2p1_WPTight_Gsf_v2;
    Int_t HLT_Ele25_eta2p1_WPTight_Gsf_v3;
    Int_t HLT_Ele25_eta2p1_WPTight_Gsf_v4;
    Int_t HLT_Ele25_eta2p1_WPTight_Gsf_v5;
    Int_t HLT_Ele25_eta2p1_WPTight_Gsf_v6;
    Int_t HLT_Ele25_eta2p1_WPTight_Gsf_v7;
    Int_t HLT_Ele27_WPTight_Gsf_v1;
    Int_t HLT_Ele27_WPTight_Gsf_v2;
    Int_t HLT_Ele27_WPTight_Gsf_v3;
    Int_t HLT_Ele27_WPTight_Gsf_v4;
    Int_t HLT_Ele27_WPTight_Gsf_v5;
    Int_t HLT_Ele27_WPTight_Gsf_v6;
    Int_t HLT_Ele27_WPTight_Gsf_v7;
    Int_t HLT_Ele32_eta2p1_WPTight_Gsf_v2;
    Int_t HLT_Ele32_eta2p1_WPTight_Gsf_v3;
    Int_t HLT_Ele32_eta2p1_WPTight_Gsf_v4;
    Int_t HLT_Ele32_eta2p1_WPTight_Gsf_v5;
    Int_t HLT_Ele32_eta2p1_WPTight_Gsf_v6;
    Int_t HLT_Ele32_eta2p1_WPTight_Gsf_v7;
    Int_t HLT_Ele32_eta2p1_WPTight_Gsf_v8;
*/
    Int_t HLT_IsoMu24_v1;
    Int_t HLT_IsoMu24_v2;
    Int_t HLT_IsoMu24_v3;
    Int_t HLT_IsoMu24_v4;
    Int_t HLT_IsoTkMu24_v1;
    Int_t HLT_IsoTkMu24_v2;
    Int_t HLT_IsoTkMu24_v3;
    Int_t HLT_IsoTkMu24_v4;
/*
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9;
*/
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6;
    Int_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2;
    Int_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3;
    Int_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5;

    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4;
    // Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5;
    // Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v6;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7;
    Int_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2;
    Int_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3;
    // Int_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4;
    // Int_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5;
    Int_t HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6;
/*
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4;
*/
    // 2017 lepton triggers
/*
    Int_t HLT_Ele32_WPTight_Gsf_L1DoubleEG_v1;
    Int_t HLT_Ele32_WPTight_Gsf_L1DoubleEG_v2;
    Int_t HLT_Ele32_WPTight_Gsf_L1DoubleEG_v3;
    Int_t HLT_Ele32_WPTight_Gsf_L1DoubleEG_v4;
    Int_t HLT_Ele32_WPTight_Gsf_L1DoubleEG_v5;
    Int_t HLT_Ele32_WPTight_Gsf_L1DoubleEG_v6;
    Int_t HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7;
    Int_t HLT_Ele35_WPTight_Gsf_v1;
    Int_t HLT_Ele35_WPTight_Gsf_v2;
    Int_t HLT_Ele35_WPTight_Gsf_v3;
    Int_t HLT_Ele35_WPTight_Gsf_v4;
    Int_t HLT_Ele35_WPTight_Gsf_v5;
    Int_t HLT_Ele35_WPTight_Gsf_v6;
    Int_t HLT_Ele35_WPTight_Gsf_v7;
*/
    Int_t HLT_IsoMu27_v8;
    Int_t HLT_IsoMu27_v9;
    Int_t HLT_IsoMu27_v10;
    Int_t HLT_IsoMu27_v11;
    Int_t HLT_IsoMu27_v12;
    Int_t HLT_IsoMu27_v13;
    Int_t HLT_IsoMu27_v14;
/*
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v10;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v11;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v12;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v13;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v14;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v15;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v16;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v17;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v14;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v15;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v16;
    Int_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17;
*/
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v8;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v9;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v10;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v11;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v12;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v13;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v14;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v1;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v2;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v3;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v4;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v7;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v8;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v1;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v2;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v3;
    Int_t HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v4;
/*
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v5;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v12;
    Int_t HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v13;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v7;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10;
    Int_t HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11;
*/
    // 2016 + 2017 triggers
/*
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4;
    Int_t HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5;
*/
    // 2017 Displaced Jet Triggers
//    Int_t HLT_HT430_DisplacedDijet40_DisplacedTrack_v10; // Already in MET trigger list
//    Int_t HLT_HT430_DisplacedDijet60_DisplacedTrack_v10; // Already in MET trigger list
//    Int_t HLT_HT430_DisplacedDijet80_DisplacedTrack_v10; // Already in MET trigger list
    Int_t HLT_HT400_DisplacedDijet40_DisplacedTrack_v10;

//    Int_t HLT_HT650_DisplacedDijet60_Inclusive_v10; // Already in MET trigger list
    Int_t HLT_HT550_DisplacedDijet80_Inclusive_v8;
    Int_t HLT_HT550_DisplacedDijet60_Inclusive_v10;
//    Int_t HLT_HT650_DisplacedDijet80_Inclusive_v11; // Already in MET trigger list
//    Int_t HLT_HT750_DisplacedDijet80_Inclusive_v11; // Already in MET trigger list

    // 2017 Photon Triggers
/*
    Int_t HLT_DoublePhoton33_CaloIdL_v5;
    Int_t HLT_DoublePhoton70_v5;
    Int_t HLT_DoublePhoton85_v13;
    Int_t HLT_TriplePhoton_20_20_20_CaloIdLV2_v2;
    Int_t HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_v2;
    Int_t HLT_TriplePhoton_30_30_10_CaloIdLV2_v3;
    Int_t HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_v3;
    Int_t HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_v3;
    Int_t HLT_Photon25_v2;
    Int_t HLT_Photon33_v4;
    Int_t HLT_Photon50_v12;
    Int_t HLT_Photon75_v12;
    Int_t HLT_Photon90_v12;
    Int_t HLT_Photon120_v12;
    Int_t HLT_Photon150_v5;
    Int_t HLT_Photon175_v13;
    Int_t HLT_Photon200_v12;
    Int_t HLT_Photon50_R9Id90_HE10_IsoM_v13;
    Int_t HLT_Photon75_R9Id90_HE10_IsoM_v13;
    Int_t HLT_Photon90_R9Id90_HE10_IsoM_v13;
    Int_t HLT_Photon120_R9Id90_HE10_IsoM_v13;
    Int_t HLT_Photon165_R9Id90_HE10_IsoM_v14;
    Int_t HLT_Photon90_CaloIdL_PFHT700_v12;

    Int_t HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v12;
    Int_t HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v12;
    Int_t HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13;
    Int_t HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12;
    Int_t HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12;
    Int_t HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13;

    Int_t HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_v4;
    Int_t HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_v7;
    Int_t HLT_Photon20_HoverELoose_v9;
    Int_t HLT_Photon30_HoverELoose_v9;
    Int_t HLT_Photon40_HoverELoose_v9;
    Int_t HLT_Photon50_HoverELoose_v9;
    Int_t HLT_Photon60_HoverELoose_v9;
*/
    // 2016 MET Triggers
    Int_t HLT_MET200_v1;
    Int_t HLT_MET200_v2;
    Int_t HLT_MET200_v3;
    Int_t HLT_MET200_v4;
    Int_t HLT_MET200_v5;
    Int_t HLT_MET250_v1;
    Int_t HLT_MET250_v2;
    Int_t HLT_MET250_v3;
    Int_t HLT_MET250_v4;
    Int_t HLT_MET250_v5;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v2;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v3;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v4;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v5;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v6;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v7;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v8;
    Int_t HLT_PFMET170_HBHECleaned_v2;
    Int_t HLT_PFMET170_HBHECleaned_v3;
    Int_t HLT_PFMET170_HBHECleaned_v4;
    Int_t HLT_PFMET170_HBHECleaned_v5;
    Int_t HLT_PFMET170_HBHECleaned_v6;
    Int_t HLT_PFMET170_HBHECleaned_v7;
    Int_t HLT_PFMET170_HBHECleaned_v8;
    Int_t HLT_PFMET170_HBHECleaned_v9;
    Int_t HLT_PFHT800_v2;
    Int_t HLT_PFHT800_v3;
    Int_t HLT_PFHT800_v4;
    Int_t HLT_PFHT800_v5;
    Int_t HLT_PFHT900_v4;
    Int_t HLT_PFHT900_v5;
    Int_t HLT_PFHT900_v6;
    Int_t HLT_PFHT750_4JetPt50_v3;
    Int_t HLT_PFHT750_4JetPt50_v4;
    Int_t HLT_PFHT750_4JetPt50_v5;
    Int_t HLT_PFHT750_4JetPt50_v6;
    Int_t HLT_PFHT750_4JetPt70_v1;
    Int_t HLT_PFHT750_4JetPt70_v2;
    Int_t HLT_PFHT750_4JetPt80_v2;
    Int_t HLT_PFHT300_PFMET100_v1;
    Int_t HLT_PFHT300_PFMET100_v2;
    Int_t HLT_PFHT300_PFMET100_v3;
    Int_t HLT_PFHT300_PFMET100_v4;
    Int_t HLT_PFHT300_PFMET110_v4;
    Int_t HLT_PFHT300_PFMET110_v5;
    Int_t HLT_PFHT300_PFMET110_v6;

    // 2017 MET triggers
    Int_t HLT_MET105_IsoTrk50_v1;
    Int_t HLT_MET105_IsoTrk50_v3;
    Int_t HLT_MET105_IsoTrk50_v4;
    Int_t HLT_MET105_IsoTrk50_v5;
    Int_t HLT_MET105_IsoTrk50_v6;
    Int_t HLT_MET105_IsoTrk50_v7;
    Int_t HLT_MET105_IsoTrk50_v8;
    Int_t HLT_MET120_IsoTrk50_v1;
    Int_t HLT_MET120_IsoTrk50_v3;
    Int_t HLT_MET120_IsoTrk50_v4;
    Int_t HLT_MET120_IsoTrk50_v5;
    Int_t HLT_MET120_IsoTrk50_v6;
    Int_t HLT_MET120_IsoTrk50_v7;
    Int_t HLT_MET120_IsoTrk50_v8;
    Int_t HLT_HT430_DisplacedDijet40_DisplacedTrack_v10; // Already in Displaced Jets trigger list
    Int_t HLT_HT430_DisplacedDijet40_DisplacedTrack_v11;
    Int_t HLT_HT430_DisplacedDijet40_DisplacedTrack_v5;
    Int_t HLT_HT430_DisplacedDijet40_DisplacedTrack_v6;
    Int_t HLT_HT430_DisplacedDijet40_DisplacedTrack_v8;
    Int_t HLT_HT430_DisplacedDijet40_DisplacedTrack_v9;
    Int_t HLT_HT430_DisplacedDijet60_DisplacedTrack_v10; // Already in Displaced Jets trigger list
    Int_t HLT_HT430_DisplacedDijet60_DisplacedTrack_v11;
    Int_t HLT_HT430_DisplacedDijet60_DisplacedTrack_v5;
    Int_t HLT_HT430_DisplacedDijet60_DisplacedTrack_v6;
    Int_t HLT_HT430_DisplacedDijet60_DisplacedTrack_v8;
    Int_t HLT_HT430_DisplacedDijet60_DisplacedTrack_v9;
    Int_t HLT_HT430_DisplacedDijet80_DisplacedTrack_v10; // Already in Displaced Jets trigger list
    Int_t HLT_HT430_DisplacedDijet80_DisplacedTrack_v11;
    Int_t HLT_HT430_DisplacedDijet80_DisplacedTrack_v5;
    Int_t HLT_HT430_DisplacedDijet80_DisplacedTrack_v6;
    Int_t HLT_HT430_DisplacedDijet80_DisplacedTrack_v8;
    Int_t HLT_HT430_DisplacedDijet80_DisplacedTrack_v9;
    Int_t HLT_HT650_DisplacedDijet60_Inclusive_v10; // Already in Displaced Jets trigger list
    Int_t HLT_HT650_DisplacedDijet60_Inclusive_v11;
    Int_t HLT_HT650_DisplacedDijet60_Inclusive_v5;
    Int_t HLT_HT650_DisplacedDijet60_Inclusive_v6;
    Int_t HLT_HT650_DisplacedDijet60_Inclusive_v8;
    Int_t HLT_HT650_DisplacedDijet60_Inclusive_v9;
    Int_t HLT_HT650_DisplacedDijet80_Inclusive_v10;
    Int_t HLT_HT650_DisplacedDijet80_Inclusive_v11; // Already in Displaced Jets trigger list
    Int_t HLT_HT650_DisplacedDijet80_Inclusive_v12;
    Int_t HLT_HT650_DisplacedDijet80_Inclusive_v6;
    Int_t HLT_HT650_DisplacedDijet80_Inclusive_v7;
    Int_t HLT_HT650_DisplacedDijet80_Inclusive_v9;
    Int_t HLT_HT750_DisplacedDijet80_Inclusive_v10;
    Int_t HLT_HT750_DisplacedDijet80_Inclusive_v11; // Already in Displaced Jets trigger list
    Int_t HLT_HT750_DisplacedDijet80_Inclusive_v12;
    Int_t HLT_HT750_DisplacedDijet80_Inclusive_v6;
    Int_t HLT_HT750_DisplacedDijet80_Inclusive_v7;
    Int_t HLT_HT750_DisplacedDijet80_Inclusive_v9;
    Int_t HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v1;
    Int_t HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v2;
    Int_t HLT_PFMET120_PFMHT120_IDTight_L1ETMnoHF_v10;
    Int_t HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v1;
    Int_t HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v2;
    Int_t HLT_PFMET120_PFMHT120_IDTight_PFHT60_v2;
    Int_t HLT_PFMET120_PFMHT120_IDTight_PFHT60_v3;
    Int_t HLT_PFMET120_PFMHT120_IDTight_PFHT60_v4;
    Int_t HLT_PFMET120_PFMHT120_IDTight_PFHT60_v5;
    Int_t HLT_PFMET120_PFMHT120_IDTight_PFHT60_v6;
    Int_t HLT_PFMET120_PFMHT120_IDTight_PFHT60_v7;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v11;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v13;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v14;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v15;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v16;
    Int_t HLT_PFMET120_PFMHT120_IDTight_v9;
    Int_t HLT_PFMET130_PFMHT130_IDTight_v11;
    Int_t HLT_PFMET130_PFMHT130_IDTight_v13;
    Int_t HLT_PFMET130_PFMHT130_IDTight_v14;
    Int_t HLT_PFMET130_PFMHT130_IDTight_v15;
    Int_t HLT_PFMET130_PFMHT130_IDTight_v16;
    Int_t HLT_PFMET130_PFMHT130_IDTight_v9;
    Int_t HLT_PFMET140_PFMHT140_IDTight_v11;
    Int_t HLT_PFMET140_PFMHT140_IDTight_v13;
    Int_t HLT_PFMET140_PFMHT140_IDTight_v14;
    Int_t HLT_PFMET140_PFMHT140_IDTight_v15;
    Int_t HLT_PFMET140_PFMHT140_IDTight_v16;
    Int_t HLT_PFMET140_PFMHT140_IDTight_v17;
    Int_t HLT_PFMET140_PFMHT140_IDTight_v18;
    Int_t HLT_PFMET140_PFMHT140_IDTight_v9;
    Int_t HLT_PFMET200_HBHE_BeamHaloCleaned_v5;
    Int_t HLT_PFMET200_HBHE_BeamHaloCleaned_v6;
    Int_t HLT_PFMET200_HBHE_BeamHaloCleaned_v7;
    Int_t HLT_PFMET250_HBHECleaned_v2;
    Int_t HLT_PFMET250_HBHECleaned_v3;
    Int_t HLT_PFMET250_HBHECleaned_v4;
    Int_t HLT_PFMET250_HBHECleaned_v5;
    Int_t HLT_PFMET250_HBHECleaned_v6;
    Int_t HLT_PFMET250_HBHECleaned_v7;
    Int_t HLT_PFMET300_HBHECleaned_v2;
    Int_t HLT_PFMET300_HBHECleaned_v3;
    Int_t HLT_PFMET300_HBHECleaned_v4;
    Int_t HLT_PFMET300_HBHECleaned_v5;
    Int_t HLT_PFMET300_HBHECleaned_v6;
    Int_t HLT_PFMET300_HBHECleaned_v7;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v1;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v2;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_L1ETMnoHF_v10;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v2;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v3;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v4;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v5;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v6;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v7;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v11;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v13;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v14;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v15;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v16;
    Int_t HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v9;
    Int_t HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v10;
    Int_t HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v12;
    Int_t HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v13;
    Int_t HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v14;
    Int_t HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v15;
    Int_t HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v9;
    Int_t HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v10;
    Int_t HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v12;
    Int_t HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v13;
    Int_t HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v14;
    Int_t HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v15;
    Int_t HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v16;
    Int_t HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v17;
    Int_t HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v9;
    Int_t HLT_PFHT1050_v11;
    Int_t HLT_PFHT1050_v12;
    Int_t HLT_PFHT1050_v13;
    Int_t HLT_PFHT1050_v14;
    Int_t HLT_PFHT1050_v15;
    Int_t HLT_PFHT1050_v16;
    Int_t HLT_PFHT1050_v7;
    Int_t HLT_PFHT1050_v9;
    Int_t HLT_PFHT180_v7;
    Int_t HLT_PFHT500_PFMET100_PFMHT100_IDTight_v1;
    Int_t HLT_PFHT500_PFMET100_PFMHT100_IDTight_v10;
    Int_t HLT_PFHT500_PFMET100_PFMHT100_IDTight_v3;
    Int_t HLT_PFHT500_PFMET100_PFMHT100_IDTight_v5;
    Int_t HLT_PFHT500_PFMET100_PFMHT100_IDTight_v6;
    Int_t HLT_PFHT500_PFMET100_PFMHT100_IDTight_v7;
    Int_t HLT_PFHT500_PFMET100_PFMHT100_IDTight_v8;
    Int_t HLT_PFHT500_PFMET100_PFMHT100_IDTight_v9;
    Int_t HLT_PFHT500_PFMET110_PFMHT110_IDTight_v1;
    Int_t HLT_PFHT500_PFMET110_PFMHT110_IDTight_v10;
    Int_t HLT_PFHT500_PFMET110_PFMHT110_IDTight_v3;
    Int_t HLT_PFHT500_PFMET110_PFMHT110_IDTight_v5;
    Int_t HLT_PFHT500_PFMET110_PFMHT110_IDTight_v6;
    Int_t HLT_PFHT500_PFMET110_PFMHT110_IDTight_v7;
    Int_t HLT_PFHT500_PFMET110_PFMHT110_IDTight_v8;
    Int_t HLT_PFHT500_PFMET110_PFMHT110_IDTight_v9;
    Int_t HLT_PFHT700_PFMET85_PFMHT85_IDTight_v1;
    Int_t HLT_PFHT700_PFMET85_PFMHT85_IDTight_v10;
    Int_t HLT_PFHT700_PFMET85_PFMHT85_IDTight_v3;
    Int_t HLT_PFHT700_PFMET85_PFMHT85_IDTight_v5;
    Int_t HLT_PFHT700_PFMET85_PFMHT85_IDTight_v6;
    Int_t HLT_PFHT700_PFMET85_PFMHT85_IDTight_v7;
    Int_t HLT_PFHT700_PFMET85_PFMHT85_IDTight_v8;
    Int_t HLT_PFHT700_PFMET85_PFMHT85_IDTight_v9;
    Int_t HLT_PFHT700_PFMET95_PFMHT95_IDTight_v1;
    Int_t HLT_PFHT700_PFMET95_PFMHT95_IDTight_v10;
    Int_t HLT_PFHT700_PFMET95_PFMHT95_IDTight_v3;
    Int_t HLT_PFHT700_PFMET95_PFMHT95_IDTight_v5;
    Int_t HLT_PFHT700_PFMET95_PFMHT95_IDTight_v6;
    Int_t HLT_PFHT700_PFMET95_PFMHT95_IDTight_v7;
    Int_t HLT_PFHT700_PFMET95_PFMHT95_IDTight_v8;
    Int_t HLT_PFHT700_PFMET95_PFMHT95_IDTight_v9;
    Int_t HLT_PFHT800_PFMET75_PFMHT75_IDTight_v1;
    Int_t HLT_PFHT800_PFMET75_PFMHT75_IDTight_v10;
    Int_t HLT_PFHT800_PFMET75_PFMHT75_IDTight_v3;
    Int_t HLT_PFHT800_PFMET75_PFMHT75_IDTight_v5;
    Int_t HLT_PFHT800_PFMET75_PFMHT75_IDTight_v6;
    Int_t HLT_PFHT800_PFMET75_PFMHT75_IDTight_v7;
    Int_t HLT_PFHT800_PFMET75_PFMHT75_IDTight_v8;
    Int_t HLT_PFHT800_PFMET75_PFMHT75_IDTight_v9;
    Int_t HLT_PFHT800_PFMET85_PFMHT85_IDTight_v1;
    Int_t HLT_PFHT800_PFMET85_PFMHT85_IDTight_v10;
    Int_t HLT_PFHT800_PFMET85_PFMHT85_IDTight_v3;
    Int_t HLT_PFHT800_PFMET85_PFMHT85_IDTight_v5;
    Int_t HLT_PFHT800_PFMET85_PFMHT85_IDTight_v6;
    Int_t HLT_PFHT800_PFMET85_PFMHT85_IDTight_v7;
    Int_t HLT_PFHT800_PFMET85_PFMHT85_IDTight_v8;
    Int_t HLT_PFHT800_PFMET85_PFMHT85_IDTight_v9;

    // MET Filters
    // 2016
    Int_t Flag_ecalLaserCorrFilter;
    Int_t Flag_chargedHadronTrackResolutionFilter;
    Int_t Flag_muonBadTrackFilter;
    Int_t Flag_badMuons;
    Int_t Flag_duplicateMuons;
    Int_t Flag_noBadMuons;
    // 2017
    Int_t Flag_BadPFMuonFilter;
    Int_t Flag_BadChargedCandidateFilter;
    Int_t Flag_ecalBadCalibFilter;
    // 2016 + 2017
    Int_t Flag_HBHENoiseFilter;
    Int_t Flag_HBHENoiseIsoFilter;
    Int_t Flag_globalTightHalo2016Filter;
    Int_t Flag_EcalDeadCellTriggerPrimitiveFilter;
    Int_t Flag_goodVertices;
    Int_t Flag_eeBadScFilter;

    // Gen info
    static constexpr size_t NGENPARMAX{1000};
    Int_t nGenPar;
    Float_t genParEta[NGENPARMAX];
    Float_t genParPhi[NGENPARMAX];
    Float_t genParE[NGENPARMAX];
    Float_t genParPt[NGENPARMAX];
    Int_t genParId[NGENPARMAX];
    Int_t genParVx[NGENPARMAX];
    Int_t genParVy[NGENPARMAX];
    Int_t genParVz[NGENPARMAX];
    Int_t genParMotherId[NGENPARMAX];
    Int_t genParMotherIndex[NGENPARMAX];
    Int_t genParNumDaughters[NGENPARMAX];
    Int_t genParDaughterId1[NGENPARMAX];
    Int_t genParDaughterId2[NGENPARMAX];
    Int_t genParStatus[NGENPARMAX];
    Int_t genParCharge[NGENPARMAX];
    Int_t eventRun;
    Int_t eventNum;
    Float_t eventLumiblock;

    // List of branches
    TBranch* b_numElePF2PAT; //!
    TBranch* b_elePF2PATE; //!
    TBranch* b_elePF2PATET; //!
    TBranch* b_elePF2PATPX; //!
    TBranch* b_elePF2PATPY; //!
    TBranch* b_elePF2PATPZ; //!
    TBranch* b_elePF2PATPhi; //!
    TBranch* b_elePF2PATTheta; //!
    TBranch* b_elePF2PATEta; //!
    TBranch* b_elePF2PATPT; //!
    TBranch* b_elePF2PATCharge; //!
    TBranch* b_elePF2PATCutIdVeto; //!
    TBranch* b_elePF2PATCutIdLoose; //!
    TBranch* b_elePF2PATCutIdMedium; //!
    TBranch* b_elePF2PATCutIdTight; //!
    TBranch* b_elePF2PATImpactTransDist; //!
    TBranch* b_elePF2PATImpactTransError; //!
    TBranch* b_elePF2PATImpactTransSignificance; //!
    TBranch* b_elePF2PATImpact3DDist; //!
    TBranch* b_elePF2PATImpact3DError; //!
    TBranch* b_elePF2PATImpact3DSignificance; //!
    TBranch* b_elePF2PATChargedHadronIso; //!
    TBranch* b_elePF2PATNeutralHadronIso; //!
    TBranch* b_elePF2PATPhotonIso; //!
    TBranch* b_elePF2PATTrackPt; //!
    TBranch* b_elePF2PATTrackPhi; //!
    TBranch* b_elePF2PATTrackEta; //!
    TBranch* b_elePF2PATTrackChi2; //!
    TBranch* b_elePF2PATTrackNDOF; //!
    TBranch* b_elePF2PATTrackD0; //!
    TBranch* b_elePF2PATTrackDBD0; //!
    TBranch* b_elePF2PATD0PV; //!
    TBranch* b_elePF2PATDZPV; //!
    TBranch* b_elePF2PATBeamSpotCorrectedTrackD0; //!
    TBranch* b_elePF2PATTrackDz; //!
    TBranch* b_elePF2PATVtxZ; //!
    TBranch* b_elePF2PATIsGsf; //!
    TBranch* b_elePF2PATGsfPx; //!
    TBranch* b_elePF2PATGsfPy; //!
    TBranch* b_elePF2PATGsfPz; //!
    TBranch* b_elePF2PATGsfE; //!
    TBranch* b_elePF2PATEcalEnergy;
    TBranch* b_elePF2PATSCEta; //!
    TBranch* b_elePF2PATSCE; //!
    TBranch* b_elePF2PATSCPhi; //!
    TBranch* b_elePF2PATSCEoverP;
    TBranch* b_elePF2PATSCSigmaEtaEta; //!
    TBranch* b_elePF2PATSCSigmaIEtaIEta; //!
    TBranch* b_elePF2PATSCSigmaIEtaIEta5x5; //!
    TBranch* b_elePF2PATSCE1x5; //!
    TBranch* b_elePF2PATSCE5x5; //!
    TBranch* b_elePF2PATSCE2x5max; //!
    TBranch* b_elePF2PATTrackIso04; //!
    TBranch* b_elePF2PATEcalIso04; //!
    TBranch* b_elePF2PATHcalIso04; //!
    TBranch* b_elePF2PATTrackIso03; //!
    TBranch* b_elePF2PATEcalIso03; //!
    TBranch* b_elePF2PATHcalIso03; //!
    TBranch* b_elePF2PATdr04EcalRecHitSumEt; //!
    TBranch* b_elePF2PATdr03EcalRecHitSumEt; //!
    TBranch* b_elePF2PATEcalIsoDeposit; //!
    TBranch* b_elePF2PATHcalIsoDeposit; //!
    TBranch* b_elePF2PATComRelIso; //!
    TBranch* b_elePF2PATComRelIsodBeta; //!
    TBranch* b_elePF2PATComRelIsoRho; //!
    TBranch* b_elePF2PATChHadIso; //!
    TBranch* b_elePF2PATNtHadIso; //!
    TBranch* b_elePF2PATGammaIso; //!
    TBranch* b_elePF2PATRhoIso; //!
    TBranch* b_elePF2PATAEff03; //!
    TBranch* b_elePF2PATMissingInnerLayers; //!
    TBranch* b_elePF2PATHoverE; //!
    TBranch* b_elePF2PATDeltaPhiSC; //!
    TBranch* b_elePF2PATDeltaEtaSC; //!
    TBranch* b_elePF2PATDeltaEtaSeedSC; //!
    TBranch* b_elePF2PATIsBarrel; //!
    TBranch* b_elePF2PATPhotonConversionTag; //!
    TBranch* b_elePF2PATPhotonConversionDist; //!
    TBranch* b_elePF2PATPhotonConversionDcot; //!
    TBranch* b_elePF2PATPhotonConversionVeto; //!
//    TBranch* b_elePF2PATPhotonConversionTagCustom; //!
    TBranch* b_elePF2PATPhotonConversionDistCustom; //!
    TBranch* b_elePF2PATPhotonConversionDcotCustom; //!
    TBranch* b_elePF2PATTriggerMatch; //!
    TBranch* b_elePF2PATJetOverlap; //!
    TBranch* b_genElePF2PATPT; //!
    TBranch* b_genElePF2PATET; //!
    TBranch* b_genElePF2PATPX; //!
    TBranch* b_genElePF2PATPY; //!
    TBranch* b_genElePF2PATPZ; //!
    TBranch* b_genElePF2PATPhi; //!
    TBranch* b_genElePF2PATTheta; //!
    TBranch* b_genElePF2PATEta; //!
    TBranch* b_genElePF2PATCharge; //!
    TBranch* b_genElePF2PATPdgId; //!
    TBranch* b_genElePF2PATMotherId; //!
    TBranch* b_genElePF2PATPromptDecayed; //!
    TBranch* b_genElePF2PATPromptFinalState; //!
    TBranch* b_genElePF2PATHardProcess; //!
    TBranch* b_numMuonPF2PAT; //!
    TBranch* b_muonPF2PATE; //!
    TBranch* b_muonPF2PATET; //!
    TBranch* b_muonPF2PATPt; //!
    TBranch* b_muonPF2PATPX; //!
    TBranch* b_muonPF2PATPY; //!
    TBranch* b_muonPF2PATPZ; //!
    TBranch* b_muonPF2PATPhi; //!
    TBranch* b_muonPF2PATTheta; //!
    TBranch* b_muonPF2PATEta; //!
    TBranch* b_muonPF2PATCharge; //!
    TBranch* b_muonPF2PATLooseCutId; //!
    TBranch* b_muonPF2PATMediumCutId; //!
    TBranch* b_muonPF2PATTightCutId; //!
    TBranch* b_muonPF2PATPfIsoVeryLoose; //!
    TBranch* b_muonPF2PATPfIsoLoose; //!
    TBranch* b_muonPF2PATPfIsoMedium; //!
    TBranch* b_muonPF2PATPfIsoTight; //!
    TBranch* b_muonPF2PATPfIsoVeryTight; //!
    TBranch* b_muonPF2PATTkIsoLoose; //!
    TBranch* b_muonPF2PATTkIsoTight; //!
    TBranch* b_muonPF2PATMvaLoose; //!
    TBranch* b_muonPF2PATMvaMedium; //!
    TBranch* b_muonPF2PATMvaTight; //!
    TBranch* b_muonPF2PATGlobalID; //!
    TBranch* b_muonPF2PATTrackID; //!
    TBranch* b_muonPF2PATChi2; //!
    TBranch* b_muonPF2PATD0; //!
    TBranch* b_muonPF2PATTrackDBD0; //!
    TBranch* b_muonPF2PATDBInnerTrackD0; //!
    TBranch* b_muonPF2PATBeamSpotCorrectedD0; //!
    TBranch* b_muonPF2PATTrackNHits; //!
    TBranch* b_muonPF2PATMuonNHits; //!
    TBranch* b_muonPF2PATNDOF; //!
    TBranch* b_muonPF2PATVertX; //!
    TBranch* b_muonPF2PATVertY; //!
    TBranch* b_muonPF2PATVertZ; //!
    TBranch* b_muonPF2PATChargedHadronIso; //!
    TBranch* b_muonPF2PATNeutralHadronIso; //!
    TBranch* b_muonPF2PATPhotonIso; //!
    TBranch* b_muonPF2PATTrackIso; //!
    TBranch* b_muonPF2PATEcalIso; //!
    TBranch* b_muonPF2PATHcalIso; //!
    TBranch* b_muonPF2PATComRelIso; //!
    TBranch* b_muonPF2PATComRelIsodBeta; //!
    TBranch* b_muonPF2PATIsPFMuon; //!
    TBranch* b_muonPF2PATNChambers; //!
    TBranch* b_muonPF2PATNMatches; //!
    TBranch* b_muonPF2PATTkLysWithMeasurements; //!
    TBranch* b_muonPF2PATVldPixHits; //!
    TBranch* b_muonPF2PATMatchedStations; //!
    TBranch* b_muonPF2PATGlbTkNormChi2; //!
    TBranch* b_muonPF2PATValidFraction; //!
    TBranch* b_muonPF2PATChi2LocalPosition; //!
    TBranch* b_muonPF2PATTrkKick; //!
    TBranch* b_muonPF2PATSegmentCompatibility; //!
    TBranch* b_muonPF2PATDBPV; //!
    TBranch* b_muonPF2PATDZPV; //!
    TBranch* b_genMuonPF2PATPT; //!
    TBranch* b_genMuonPF2PATET; //!
    TBranch* b_genMuonPF2PATPX; //!
    TBranch* b_genMuonPF2PATPY; //!
    TBranch* b_genMuonPF2PATPZ; //!
    TBranch* b_genMuonPF2PATPhi; //!
    TBranch* b_genMuonPF2PATTheta; //!
    TBranch* b_genMuonPF2PATEta; //!
    TBranch* b_genMuonPF2PATCharge; //!
    TBranch* b_genMuonPF2PATPdgId; //!
    TBranch* b_genMuonPF2PATMotherId; //!
    TBranch* b_genMuonPF2PATPromptDecayed; //!
    TBranch* b_genMuonPF2PATPromptFinalState; //!
    TBranch* b_genMuonPF2PATHardProcess; //!
    TBranch* b_numJetPF2PAT; //!
    TBranch* b_jetPF2PATE; //!
    TBranch* b_jetPF2PATEt; //!
    TBranch* b_jetPF2PATPt; //!
    TBranch* b_jetPF2PATPtRaw; //!
    TBranch* b_jetPF2PATUnCorEt; //!
    TBranch* b_jetPF2PATUnCorPt; //!
    TBranch* b_jetPF2PATEta; //!
    TBranch* b_jetPF2PATTheta; //!
    TBranch* b_jetPF2PATPhi; //!
    TBranch* b_jetPF2PATPx; //!
    TBranch* b_jetPF2PATPy; //!
    TBranch* b_jetPF2PATPz; //!
    TBranch* b_jetPF2PATdRClosestLepton; //!
    TBranch* b_jetPF2PATNtracksInJet; //!
    TBranch* b_jetPF2PATJetCharge; //!
    TBranch* b_jetPF2PATfHPD; //!
    TBranch* b_jetPF2PATBtagSoftMuonPtRel; //!
    TBranch* b_jetPF2PATBtagSoftMuonQuality; //!
    TBranch* b_jetPF2PATCorrFactor; //!
    TBranch* b_jetPF2PATCorrResidual; //!
    TBranch* b_jetPF2PATL2L3ResErr; //!
    TBranch* b_jetPF2PATCorrErrLow; //!
    TBranch* b_jetPF2PATCorrErrHi; //!
    TBranch* b_jetPF2PATN90Hits; //!
    TBranch* b_jetPF2PATTriggered; //!
    TBranch* b_jetPF2PATSVX; //!
    TBranch* b_jetPF2PATSVY; //!
    TBranch* b_jetPF2PATSVZ; //!
    TBranch* b_jetPF2PATSVDX; //!
    TBranch* b_jetPF2PATSVDY; //!
    TBranch* b_jetPF2PATSVDZ; //!
    TBranch* b_jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags; //!
    TBranch* b_jetPF2PATpfCombinedCvsLJetTags; //!
    TBranch* b_jetPF2PATpfCombinedCvsBJetTags; //!
    TBranch* b_jetPF2PATNConstituents; //!
    TBranch* b_jetPF2PATPID; //!
    TBranch* b_jetPF2PATClosestBPartonDeltaR; //!
    TBranch* b_jetPF2PATClosestCPartonDeltaR; //!
    TBranch* b_genJetPF2PATE; //!
    TBranch* b_genJetPF2PATET; //!
    TBranch* b_genJetPF2PATPT; //!
    TBranch* b_genJetPF2PATPX; //!
    TBranch* b_genJetPF2PATPY; //!
    TBranch* b_genJetPF2PATPZ; //!
    TBranch* b_genJetPF2PATMass; //!
    TBranch* b_genJetPF2PATPhi; //!
    TBranch* b_genJetPF2PATTheta; //!
    TBranch* b_genJetPF2PATEta; //!
    TBranch* b_genJetPF2PATPID; //!
    TBranch* b_genJetPF2PATMotherPID; //!
    TBranch* b_genJetPF2PATScalarAncestor; //!
    TBranch* b_jetPF2PATMuEnergy; //!
    TBranch* b_jetPF2PATMuEnergyFraction; //!
    TBranch* b_jetPF2PATNeutralHadEnergy; //!
    TBranch* b_jetPF2PATNeutralEmEnergy; //!
    TBranch* b_jetPF2PATChargedHadronEnergyFraction; //!
    TBranch* b_jetPF2PATNeutralHadronEnergyFraction; //!
    TBranch* b_jetPF2PATChargedEmEnergyFraction; //!
    TBranch* b_jetPF2PATNeutralEmEnergyFraction; //!
    TBranch* b_jetPF2PATMuonFraction; //!
    TBranch* b_jetPF2PATChargedHadronEnergyFractionCorr; //!
    TBranch* b_jetPF2PATNeutralHadronEnergyFractionCorr; //!
    TBranch* b_jetPF2PATChargedEmEnergyFractionCorr; //!
    TBranch* b_jetPF2PATNeutralEmEnergyFractionCorr; //!
    TBranch* b_jetPF2PATMuonFractionCorr; //!
    TBranch* b_jetPF2PATNeutralMultiplicity; //!
    TBranch* b_jetPF2PATChargedMultiplicity; //!
    TBranch* b_fixedGridRhoFastjetAll; //!
    TBranch* b_metPF2PATE; //!
    TBranch* b_metPF2PATEt; //!
    TBranch* b_metPF2PATEtRaw; //!
    TBranch* b_metPF2PATPhi; //!
    TBranch* b_metPF2PATPt; //!
    TBranch* b_metPF2PATPx; //!
    TBranch* b_metPF2PATPy; //!
    TBranch* b_metPF2PATPz; //!
    TBranch* b_metPF2PATScalarEt; //!
    TBranch* b_metPF2PATEtUncorrected; //!
    TBranch* b_metPF2PATPhiUncorrected; //!
    TBranch* b_metPF2PATUnclusteredEnUp;
    TBranch* b_metPF2PATUnclusteredEnDown;
    TBranch* b_genMetPF2PATE; //!
    TBranch* b_genMetPF2PATEt; //!
    TBranch* b_genMetPF2PATPhi; //!
    TBranch* b_genMetPF2PATPt; //!
    TBranch* b_genMetPF2PATPx; //!
    TBranch* b_genMetPF2PATPy; //!
    TBranch* b_genMetPF2PATPz; //!
    TBranch* b_numTauPF2PAT; //!
    TBranch* b_tauPF2PATE; //!
    TBranch* b_tauPF2PATPt; //!
    TBranch* b_tauPF2PATPhi; //!
    TBranch* b_tauPF2PATEta; //!
    TBranch* b_numPhoPF2PAT; //!
    TBranch* b_phoPF2PATE; //!
    TBranch* b_phoPF2PATSigmaE; //!
    TBranch* b_phoPF2PATET; //!
    TBranch* b_phoPF2PATPhi; //!
    TBranch* b_phoPF2PATEta; //!
    TBranch* b_phoPF2PATTheta; //!
    TBranch* b_phoPF2PATPt; //!
    TBranch* b_phoPF2PATPx; //!
    TBranch* b_phoPF2PATPy; //!
    TBranch* b_phoPF2PATPz; //!
    TBranch* b_phoPF2PATCalibE; //!
    TBranch* b_phoPF2PATCalibEt; //!
    TBranch* b_phoPF2PATSCE; //!
    TBranch* b_phoPF2PATSCRawE; //!
    TBranch* b_phoPF2PATESEnP1; //!
    TBranch* b_phoPF2PATESEnP2; //!
    TBranch* b_phoPF2PATSCEta; //!
    TBranch* b_phoPF2PATSCEtaWidth; //!
    TBranch* b_phoPF2PATSCPhi; //!
    TBranch* b_phoPF2PATSCPhiWidth; //!
    TBranch* b_phoPF2PATSCBrem; //!
    TBranch* b_phoPF2PATHasPixelSeed; //!
    TBranch* b_phoPF2PATEleVeto; //!
    TBranch* b_phoPF2PATR9; //!
    TBranch* b_phoPF2PATHoverE; //!
    TBranch* b_phoPF2PATESEffSigmaRR; //!
    TBranch* b_phoPF2PATSigmaIEtaIEtaFull5x5; //!
    TBranch* b_phoPF2PATSigmaIEtaIPhiFull5x5; //!
    TBranch* b_phoPF2PATSigmaIPhiIPhiFull5x5; //!
    TBranch* b_phoPF2PATE2x2Full5x5; //!
    TBranch* b_phoPF2PATE5x5Full5x5; //!
    TBranch* b_phoPF2PATR9Full5x5; //!
    TBranch* b_phoPF2PATPFChIso; //!
    TBranch* b_phoPF2PATPFPhoIso ; //!
    TBranch* b_phoPF2PATPFNeuIso; //!
    TBranch* b_phoPF2PATPFChWorstIso; //!
    TBranch* b_phoPF2PATMIPTotEnergy; //!
    TBranch* b_phoPF2PATCutIdLoose; //!
    TBranch* b_phoPF2PATCutIdMedium; //!
    TBranch* b_phoPF2PATCutIdTight; //!
    TBranch* b_phoPF2PATMvaIdWp80; //!
    TBranch* b_phoPF2PATMvaIdWp90; //!
    TBranch* b_genPhoPF2PATPt; //!
    TBranch* b_genPhoPF2PATET; //!
    TBranch* b_genPhoPF2PATEta; //!
    TBranch* b_genPhoPF2PATTheta; //!
    TBranch* b_genPhoPF2PATPhi; //!
    TBranch* b_genPhoPF2PATPx; //!
    TBranch* b_genPhoPF2PATPy; //!
    TBranch* b_genPhoPF2PATPz; //!
    TBranch* b_genPhoPF2PATCharge; //!
    TBranch* b_genPhoPF2PATPdgId; //!
    TBranch* b_genPhoPF2PATMotherId; //!
    TBranch* b_genPhoPF2PATIsPhoton; //!
    TBranch* b_genPhoPF2PATIsConvertedPhoton; //!
    TBranch* b_genPhoPF2PATIsJet; //!
/*
    TBranch* b_numPhoOOT_PF2PAT; //!
    TBranch* b_phoOOT_PF2PATE; //!
    TBranch* b_phoOOT_PF2PATSigmaE; //!
    TBranch* b_phoOOT_PF2PATET; //!
    TBranch* b_phoOOT_PF2PATPhi; //!
    TBranch* b_phoOOT_PF2PATEta; //!
    TBranch* b_phoOOT_PF2PATTheta; //!
    TBranch* b_phoOOT_PF2PATPt; //!
    TBranch* b_phoOOT_PF2PATPx; //!
    TBranch* b_phoOOT_PF2PATPy; //!
    TBranch* b_phoOOT_PF2PATPz; //!
    TBranch* b_phoOOT_PF2PATCalibE; //!
    TBranch* b_phoOOT_PF2PATCalibEt; //!
    TBranch* b_phoOOT_PF2PATSCE; //!
    TBranch* b_phoOOT_PF2PATSCRawE; //!
    TBranch* b_phoOOT_PF2PATESEnP1; //!
    TBranch* b_phoOOT_PF2PATESEnP2; //!
    TBranch* b_phoOOT_PF2PATSCEta; //!
    TBranch* b_phoOOT_PF2PATSCEtaWidth; //!
    TBranch* b_phoOOT_PF2PATSCPhi; //!
    TBranch* b_phoOOT_PF2PATSCPhiWidth; //!
    TBranch* b_phoOOT_PF2PATSCBrem; //!
    TBranch* b_phoOOT_PF2PATHasPixelSeed; //!
    TBranch* b_phoOOT_PF2PATEleVeto; //!
    TBranch* b_phoOOT_PF2PATR9; //!
    TBranch* b_phoOOT_PF2PATHoverE; //!
    TBranch* b_phoOOT_PF2PATESEffSigmaRR; //!
    TBranch* b_phoOOT_PF2PATSigmaIEtaIEtaFull5x5; //!
    TBranch* b_phoOOT_PF2PATSigmaIEtaIPhiFull5x5; //!
    TBranch* b_phoOOT_PF2PATSigmaIPhiIPhiFull5x5; //!
    TBranch* b_phoOOT_PF2PATE2x2Full5x5; //!
    TBranch* b_phoOOT_PF2PATE5x5Full5x5; //!
    TBranch* b_phoOOT_PF2PATR9Full5x5; //!
    TBranch* b_phoOOT_PF2PATPFChIso; //!
    TBranch* b_phoOOT_PF2PATPFPhoIso ; //!
    TBranch* b_phoOOT_PF2PATPFNeuIso; //!
    TBranch* b_phoOOT_PF2PATPFChWorstIso; //!
    TBranch* b_phoOOT_PF2PATMIPTotEnergy; //!
    TBranch* b_phoOOT_PF2PATCutIdLoose; //!
    TBranch* b_phoOOT_PF2PATCutIdMedium; //!
    TBranch* b_phoOOT_PF2PATCutIdTight; //!
    TBranch* b_phoOOT_PF2PATMvaIdWp80; //!
    TBranch* b_phoOOT_PF2PATMvaIdWp90; //!
    TBranch* b_genPhoOOT_PF2PATPt; //!
    TBranch* b_genPhoOOT_PF2PATET; //!
    TBranch* b_genPhoOOT_PF2PATEta; //!
    TBranch* b_genPhoOOT_PF2PATTheta; //!
    TBranch* b_genPhoOOT_PF2PATPhi; //!
    TBranch* b_genPhoOOT_PF2PATPx; //!
    TBranch* b_genPhoOOT_PF2PATPy; //!
    TBranch* b_genPhoOOT_PF2PATPz; //!
    TBranch* b_genPhoOOT_PF2PATCharge; //!
    TBranch* b_genPhoOOT_PF2PATPdgId; //!
    TBranch* b_genPhoOOT_PF2PATMotherId; //!
    TBranch* b_genPhoOOT_PF2PATIsPhoton; //!
    TBranch* b_genPhoOOT_PF2PATIsConvertedPhoton; //!
    TBranch* b_genPhoOOT_PF2PATIsJet; //!
*/
    TBranch* b_numPackedCands; //!
//    TBranch* b_packedCandsPt; //!
    TBranch* b_packedCandsPx; //!
    TBranch* b_packedCandsPy; //!
    TBranch* b_packedCandsPz; //!
    TBranch* b_packedCandsE; //!
//    TBranch* b_packedCandsEta; //!
//    TBranch* b_packedCandsTheta; //!
//    TBranch* b_packedCandsPhi; //!
    TBranch* b_packedCandsCharge; //!
    TBranch* b_packedCandsPdgId; //!
//    TBranch* b_packedCandsTime; //!
//    TBranch* b_packedCandsFromPV; //!
//    TBranch* b_packedCandsPVquality; //!
//    TBranch* b_packedCandsVx; //!
//    TBranch* b_packedCandsVy; //!
//    TBranch* b_packedCandsVz; //!
//    TBranch* b_packedCandsVEta; //!
//    TBranch* b_packedCandsVPhi; //!
//    TBranch* b_packedCandsBeamSpotCorrectedD0; //!
    TBranch* b_packedCandsDz; //!
    TBranch* b_packedCandsDxy; //!
//    TBranch* b_packedCandsDzAssocPV; //!
//    TBranch* b_packedCandsVtxChi2Norm; //!
    TBranch* b_packedCandsHasTrackDetails; //!
//    TBranch* b_packedCandsDzError; //!
//    TBranch* b_packedCandsDxyError; //!
//    TBranch* b_packedCandsTimeError; //!
//    TBranch* b_packedCandsPseudoTrkPt; //!
    TBranch* b_packedCandsPseudoTrkPx; //!
    TBranch* b_packedCandsPseudoTrkPy; //!
    TBranch* b_packedCandsPseudoTrkPz; //!
//    TBranch* b_packedCandsPseudoTrkEta; //!
//    TBranch* b_packedCandsPseudoTrkPhi; //!
    TBranch* b_packedCandsPseudoTrkCharge; //!
    TBranch* b_packedCandsPseudoTrkVx; //!
    TBranch* b_packedCandsPseudoTrkVy; //!
    TBranch* b_packedCandsPseudoTrkVz; //!
//    TBranch* b_packedCandsPseudoTrkChi2Norm; //!
//    TBranch* b_packedCandsPseudoTrkNumberOfHits; //!
//    TBranch* b_packedCandsPseudoTrkNumberOfPixelHits; //!
//    TBranch* b_packedCandsPseudoTrkPixelLayersWithMeasurement; //!
//    TBranch* b_packedCandsPseudoTrkStripLayersWithMeasurement; //!
//    TBranch* b_packedCandsPseudoTrkTrackerLayersWithMeasurement; //!
    TBranch* b_packedCandsHighPurityTrack; //!

    TBranch* b_numIsolatedTracks; //!
    TBranch* b_isoTracksPt; //!
    TBranch* b_isoTracksPx; //!
    TBranch* b_isoTracksPy; //!
    TBranch* b_isoTracksPz; //!
    TBranch* b_isoTracksE; //!
    TBranch* b_isoTracksEta; //!
    TBranch* b_isoTracksTheta; //!
    TBranch* b_isoTracksPhi; //!
    TBranch* b_isoTracksCharge; //!
    TBranch* b_isoTracksPdgId; //!
    TBranch* b_isoTracksMatchedCaloJetEmEnergy; //!
    TBranch* b_isoTracksMatchedCaloJetHadEnergy; //!
    TBranch* b_isoTracksDz; //!
    TBranch* b_isoTracksDxy; //!
    TBranch* b_isoTracksDzError; //!
    TBranch* b_isoTracksDxyError; //!
    TBranch* b_isoTracksFromPV; //!
    TBranch* b_isoTracksVx; //!
    TBranch* b_isoTracksVy; //!
    TBranch* b_isoTracksVz; //!
    TBranch* b_isoTracksHighPurity; //!
    TBranch* b_isoTracksTight; //!
    TBranch* b_isoTracksLoose; //!
    TBranch* b_isoTracksDeltaEta; //!
    TBranch* b_isoTracksDeltaPhi; //!

    TBranch* b_isElePlusJets; //!
    TBranch* b_genPDFScale; //!
    TBranch* b_genPDFx1; //!
    TBranch* b_genPDFx2; //!
    TBranch* b_genPDFf1; //!
    TBranch* b_genPDFf2; //!
    TBranch* b_topPtReweight; //!
    TBranch* b_processId; //!
    TBranch* b_processPtHat; //!
    TBranch* b_processMCWeight; //!
    TBranch* b_beamSpotX; //!
    TBranch* b_beamSpotY; //!
    TBranch* b_beamSpotZ; //!
    TBranch* b_numPVs; //!
    TBranch* b_pvX; //!
    TBranch* b_pvY; //!
    TBranch* b_pvZ; //!
    TBranch* b_pvDX; //!
    TBranch* b_pvDY; //!
    TBranch* b_pvDZ; //!
    TBranch* b_pvRho; //!
    TBranch* b_pvIsFake; //!
    TBranch* b_pvNdof; //!
    TBranch* b_pvChi2; //!
    TBranch* b_numSVs; //!
    TBranch* b_svPt; //!
    TBranch* b_svPx; //!
    TBranch* b_svPy; //!
    TBranch* b_svPz; //!
    TBranch* b_svMass; //!
    TBranch* b_svE; //!
    TBranch* b_svEta; //!
    TBranch* b_svTheta; //!
    TBranch* b_svPhi; //!
    TBranch* b_svX; //!
    TBranch* b_svY; //!
    TBranch* b_svZ; //!
    TBranch* b_svVertexChi2; //!
    TBranch* b_svVertexNdof; //!
    TBranch* b_svNtracks; //!
    TBranch* b_svDist3D; //!
    TBranch* b_svDist3DSig; //!
    TBranch* b_svDist3DError; //!
    TBranch* b_svDistXY; //!
    TBranch* b_svDistXYSig; //!
    TBranch* b_svDistXYError; //!
    TBranch* b_svAnglePV; //!
    TBranch* b_svIsLambda; //!
    TBranch* b_svIsKshort; //!
    TBranch* b_mhtPt; //!
    TBranch* b_mhtPy; //!
    TBranch* b_mhtPx; //!
    TBranch* b_mhtPhi; //!
    TBranch* b_mhtSumEt; //!
    TBranch* b_mhtSignif; //!
    TBranch* b_nTriggerBits; //!
    TBranch* b_TriggerBits; //!
    TBranch* b_weight_muF0p5; //!
    TBranch* b_weight_muF2; //!
    TBranch* b_weight_muR0p5; //!
    TBranch* b_weight_muR2; //!
    TBranch* b_weight_muF0p5muR0p5; //!
    TBranch* b_weight_muF2muR2; //!
    TBranch* b_origWeightForNorm; //!
    TBranch* b_weight_pdfMax; //!
    TBranch* b_weight_pdfMin; //!
    TBranch* b_weight_alphaMax; //!
    TBranch* b_weight_alphaMin; //!
    TBranch* b_isrRedHi;
    TBranch* b_fsrRedHi;
    TBranch* b_isrRedLo;
    TBranch* b_fsrRedLo;
    TBranch* b_isrDefHi;
    TBranch* b_fsrDefHi;
    TBranch* b_isrDefLo;
    TBranch* b_fsrDefLo;
    TBranch* b_isrConHi;
    TBranch* b_fsrConHi;
    TBranch* b_isrConLo;
    TBranch* b_fsrConLo;
    //   TBranch        *b_numVert;    //!
/*
    TBranch* b_HLT_Ele25_eta2p1_WPTight_Gsf_v1;
    TBranch* b_HLT_Ele25_eta2p1_WPTight_Gsf_v2;
    TBranch* b_HLT_Ele25_eta2p1_WPTight_Gsf_v3;
    TBranch* b_HLT_Ele25_eta2p1_WPTight_Gsf_v4;
    TBranch* b_HLT_Ele25_eta2p1_WPTight_Gsf_v5;
    TBranch* b_HLT_Ele25_eta2p1_WPTight_Gsf_v6;
    TBranch* b_HLT_Ele25_eta2p1_WPTight_Gsf_v7;
    TBranch* b_HLT_Ele27_WPTight_Gsf_v1;
    TBranch* b_HLT_Ele27_WPTight_Gsf_v2;
    TBranch* b_HLT_Ele27_WPTight_Gsf_v3;
    TBranch* b_HLT_Ele27_WPTight_Gsf_v4;
    TBranch* b_HLT_Ele27_WPTight_Gsf_v5;
    TBranch* b_HLT_Ele27_WPTight_Gsf_v6;
    TBranch* b_HLT_Ele27_WPTight_Gsf_v7;
    TBranch* b_HLT_Ele32_eta2p1_WPTight_Gsf_v2;
    TBranch* b_HLT_Ele32_eta2p1_WPTight_Gsf_v3;
    TBranch* b_HLT_Ele32_eta2p1_WPTight_Gsf_v4;
    TBranch* b_HLT_Ele32_eta2p1_WPTight_Gsf_v5;
    TBranch* b_HLT_Ele32_eta2p1_WPTight_Gsf_v6;
    TBranch* b_HLT_Ele32_eta2p1_WPTight_Gsf_v7;
    TBranch* b_HLT_Ele32_eta2p1_WPTight_Gsf_v8;
*/
    TBranch* b_HLT_IsoMu24_v1;
    TBranch* b_HLT_IsoMu24_v2;
    TBranch* b_HLT_IsoMu24_v3;
    TBranch* b_HLT_IsoMu24_v4;
    TBranch* b_HLT_IsoTkMu24_v1;
    TBranch* b_HLT_IsoTkMu24_v2;
    TBranch* b_HLT_IsoTkMu24_v3;
    TBranch* b_HLT_IsoTkMu24_v4;
/*
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9;
*/
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4;
    // TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5;
    // TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v6;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6;
    TBranch* b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2;
    TBranch* b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3;
    // TBranch* b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4;
    // TBranch* b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5;
    TBranch* b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6;
    TBranch* b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2;
    TBranch* b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3;
    TBranch* b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5;
/*
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4;
*/
/*
    TBranch* b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v1;
    TBranch* b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v2;
    TBranch* b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v3;
    TBranch* b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v4;
    TBranch* b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v5;
    TBranch* b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v6;
    TBranch* b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7;
    TBranch* b_HLT_Ele35_WPTight_Gsf_v1;
    TBranch* b_HLT_Ele35_WPTight_Gsf_v2;
    TBranch* b_HLT_Ele35_WPTight_Gsf_v3;
    TBranch* b_HLT_Ele35_WPTight_Gsf_v4;
    TBranch* b_HLT_Ele35_WPTight_Gsf_v5;
    TBranch* b_HLT_Ele35_WPTight_Gsf_v6;
    TBranch* b_HLT_Ele35_WPTight_Gsf_v7;
*/
    TBranch* b_HLT_IsoMu27_v8;
    TBranch* b_HLT_IsoMu27_v9;
    TBranch* b_HLT_IsoMu27_v10;
    TBranch* b_HLT_IsoMu27_v11;
    TBranch* b_HLT_IsoMu27_v12;
    TBranch* b_HLT_IsoMu27_v13;
    TBranch* b_HLT_IsoMu27_v14;
/*
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v10;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v11;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v12;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v13;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v14;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v15;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v16;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v17;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v14;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v15;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v16;
    TBranch* b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17;
*/
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v8;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v9;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v10;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v11;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v12;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v13;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v14;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v1;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v2;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v3;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v4;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v7;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v8;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v1;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v2;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v3;
    TBranch* b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v4;
/*
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12;
    TBranch* b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v5;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v12;
    TBranch* b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v13;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v7;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10;
    TBranch* b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11;
*/
//    TBranch* b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v10; // Already in MET trigger list
//    TBranch* b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v10; // Already in MET trigger list
//    TBranch* b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v10; // Already in MET trigger list
    TBranch* b_HLT_HT400_DisplacedDijet40_DisplacedTrack_v10;
//    TBranch* b_HLT_HT650_DisplacedDijet60_Inclusive_v10; // Already in MET trigger list
    TBranch* b_HLT_HT550_DisplacedDijet80_Inclusive_v8;
    TBranch* b_HLT_HT550_DisplacedDijet60_Inclusive_v10;
//    TBranch* b_HLT_HT650_DisplacedDijet80_Inclusive_v11; // Already in MET trigger list
//    TBranch* b_HLT_HT750_DisplacedDijet80_Inclusive_v11; // Already in MET trigger list
/*
    TBranch* b_HLT_DoublePhoton33_CaloIdL_v5;
    TBranch* b_HLT_DoublePhoton70_v5;
    TBranch* b_HLT_DoublePhoton85_v13;
    TBranch* b_HLT_TriplePhoton_20_20_20_CaloIdLV2_v2;
    TBranch* b_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_v2;
    TBranch* b_HLT_TriplePhoton_30_30_10_CaloIdLV2_v3;
    TBranch* b_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_v3;
    TBranch* b_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_v3;
    TBranch* b_HLT_Photon25_v2;
    TBranch* b_HLT_Photon33_v4;
    TBranch* b_HLT_Photon50_v12;
    TBranch* b_HLT_Photon75_v12;
    TBranch* b_HLT_Photon90_v12;
    TBranch* b_HLT_Photon120_v12;
    TBranch* b_HLT_Photon150_v5;
    TBranch* b_HLT_Photon175_v13;
    TBranch* b_HLT_Photon200_v12;
    TBranch* b_HLT_Photon50_R9Id90_HE10_IsoM_v13;
    TBranch* b_HLT_Photon75_R9Id90_HE10_IsoM_v13;
    TBranch* b_HLT_Photon90_R9Id90_HE10_IsoM_v13;
    TBranch* b_HLT_Photon120_R9Id90_HE10_IsoM_v13;
    TBranch* b_HLT_Photon165_R9Id90_HE10_IsoM_v14;
    TBranch* b_HLT_Photon90_CaloIdL_PFHT700_v12;

    TBranch* b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v12;
    TBranch* b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v12;
    TBranch* b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13;
    TBranch* b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12;
    TBranch* b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12;
    TBranch* b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13;

    TBranch* b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_v4;
    TBranch* b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_v7;
    TBranch* b_HLT_Photon20_HoverELoose_v9;
    TBranch* b_HLT_Photon30_HoverELoose_v9;
    TBranch* b_HLT_Photon40_HoverELoose_v9;
    TBranch* b_HLT_Photon50_HoverELoose_v9;
    TBranch* b_HLT_Photon60_HoverELoose_v9;
*/

    TBranch* b_HLT_MET200_v1;
    TBranch* b_HLT_MET200_v2;
    TBranch* b_HLT_MET200_v3;
    TBranch* b_HLT_MET200_v4;
    TBranch* b_HLT_MET200_v5;
    TBranch* b_HLT_MET250_v1;
    TBranch* b_HLT_MET250_v2;
    TBranch* b_HLT_MET250_v3;
    TBranch* b_HLT_MET250_v4;
    TBranch* b_HLT_MET250_v5;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v2;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v3;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v4;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v5;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v6;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v7;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v8;
    TBranch* b_HLT_PFMET170_HBHECleaned_v2;
    TBranch* b_HLT_PFMET170_HBHECleaned_v3;
    TBranch* b_HLT_PFMET170_HBHECleaned_v4;
    TBranch* b_HLT_PFMET170_HBHECleaned_v5;
    TBranch* b_HLT_PFMET170_HBHECleaned_v6;
    TBranch* b_HLT_PFMET170_HBHECleaned_v7;
    TBranch* b_HLT_PFMET170_HBHECleaned_v8;
    TBranch* b_HLT_PFMET170_HBHECleaned_v9;
    TBranch* b_HLT_PFHT800_v2;
    TBranch* b_HLT_PFHT800_v3;
    TBranch* b_HLT_PFHT800_v4;
    TBranch* b_HLT_PFHT800_v5;
    TBranch* b_HLT_PFHT900_v4;
    TBranch* b_HLT_PFHT900_v5;
    TBranch* b_HLT_PFHT900_v6;
    TBranch* b_HLT_PFHT750_4JetPt50_v3;
    TBranch* b_HLT_PFHT750_4JetPt50_v4;
    TBranch* b_HLT_PFHT750_4JetPt50_v5;
    TBranch* b_HLT_PFHT750_4JetPt50_v6;
    TBranch* b_HLT_PFHT750_4JetPt70_v1;
    TBranch* b_HLT_PFHT750_4JetPt70_v2;
    TBranch* b_HLT_PFHT750_4JetPt80_v2;
    TBranch* b_HLT_PFHT300_PFMET100_v1;
    TBranch* b_HLT_PFHT300_PFMET100_v2;
    TBranch* b_HLT_PFHT300_PFMET100_v3;
    TBranch* b_HLT_PFHT300_PFMET100_v4;
    TBranch* b_HLT_PFHT300_PFMET110_v4;
    TBranch* b_HLT_PFHT300_PFMET110_v5;
    TBranch* b_HLT_PFHT300_PFMET110_v6;

    TBranch* b_HLT_MET105_IsoTrk50_v1;
    TBranch* b_HLT_MET105_IsoTrk50_v3;
    TBranch* b_HLT_MET105_IsoTrk50_v4;
    TBranch* b_HLT_MET105_IsoTrk50_v5;
    TBranch* b_HLT_MET105_IsoTrk50_v6;
    TBranch* b_HLT_MET105_IsoTrk50_v7;
    TBranch* b_HLT_MET105_IsoTrk50_v8;
    TBranch* b_HLT_MET120_IsoTrk50_v1;
    TBranch* b_HLT_MET120_IsoTrk50_v3;
    TBranch* b_HLT_MET120_IsoTrk50_v4;
    TBranch* b_HLT_MET120_IsoTrk50_v5;
    TBranch* b_HLT_MET120_IsoTrk50_v6;
    TBranch* b_HLT_MET120_IsoTrk50_v7;
    TBranch* b_HLT_MET120_IsoTrk50_v8;
    TBranch* b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v10; // Already in Displaced Jets trigger list
    TBranch* b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v11;
    TBranch* b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v5;
    TBranch* b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v6;
    TBranch* b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v8;
    TBranch* b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v9;
    TBranch* b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v10; // Already in Displaced Jets trigger list
    TBranch* b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v11;
    TBranch* b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v5;
    TBranch* b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v6;
    TBranch* b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v8;
    TBranch* b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v9;
    TBranch* b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v10; // Already in Displaced Jets trigger list
    TBranch* b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v11;
    TBranch* b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v5;
    TBranch* b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v6;
    TBranch* b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v8;
    TBranch* b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v9;
    TBranch* b_HLT_HT650_DisplacedDijet60_Inclusive_v10; // Already in Displaced Jets trigger list
    TBranch* b_HLT_HT650_DisplacedDijet60_Inclusive_v11;
    TBranch* b_HLT_HT650_DisplacedDijet60_Inclusive_v5;
    TBranch* b_HLT_HT650_DisplacedDijet60_Inclusive_v6;
    TBranch* b_HLT_HT650_DisplacedDijet60_Inclusive_v8;
    TBranch* b_HLT_HT650_DisplacedDijet60_Inclusive_v9;
    TBranch* b_HLT_HT650_DisplacedDijet80_Inclusive_v10;
    TBranch* b_HLT_HT650_DisplacedDijet80_Inclusive_v11; // Already in Displaced Jets trigger list
    TBranch* b_HLT_HT650_DisplacedDijet80_Inclusive_v12;
    TBranch* b_HLT_HT650_DisplacedDijet80_Inclusive_v6;
    TBranch* b_HLT_HT650_DisplacedDijet80_Inclusive_v7;
    TBranch* b_HLT_HT650_DisplacedDijet80_Inclusive_v9;
    TBranch* b_HLT_HT750_DisplacedDijet80_Inclusive_v10;
    TBranch* b_HLT_HT750_DisplacedDijet80_Inclusive_v11; // Already in Displaced Jets trigger list
    TBranch* b_HLT_HT750_DisplacedDijet80_Inclusive_v12;
    TBranch* b_HLT_HT750_DisplacedDijet80_Inclusive_v6;
    TBranch* b_HLT_HT750_DisplacedDijet80_Inclusive_v7;
    TBranch* b_HLT_HT750_DisplacedDijet80_Inclusive_v9;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v1;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v2;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_L1ETMnoHF_v10;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v1;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v2;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v2;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v3;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v4;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v5;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v6;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v7;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v11;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v13;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v14;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v15;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v16;
    TBranch* b_HLT_PFMET120_PFMHT120_IDTight_v9;
    TBranch* b_HLT_PFMET130_PFMHT130_IDTight_v11;
    TBranch* b_HLT_PFMET130_PFMHT130_IDTight_v13;
    TBranch* b_HLT_PFMET130_PFMHT130_IDTight_v14;
    TBranch* b_HLT_PFMET130_PFMHT130_IDTight_v15;
    TBranch* b_HLT_PFMET130_PFMHT130_IDTight_v16;
    TBranch* b_HLT_PFMET130_PFMHT130_IDTight_v9;
    TBranch* b_HLT_PFMET140_PFMHT140_IDTight_v11;
    TBranch* b_HLT_PFMET140_PFMHT140_IDTight_v13;
    TBranch* b_HLT_PFMET140_PFMHT140_IDTight_v14;
    TBranch* b_HLT_PFMET140_PFMHT140_IDTight_v15;
    TBranch* b_HLT_PFMET140_PFMHT140_IDTight_v16;
    TBranch* b_HLT_PFMET140_PFMHT140_IDTight_v17;
    TBranch* b_HLT_PFMET140_PFMHT140_IDTight_v18;
    TBranch* b_HLT_PFMET140_PFMHT140_IDTight_v9;
    TBranch* b_HLT_PFMET200_HBHE_BeamHaloCleaned_v5;
    TBranch* b_HLT_PFMET200_HBHE_BeamHaloCleaned_v6;
    TBranch* b_HLT_PFMET200_HBHE_BeamHaloCleaned_v7;
    TBranch* b_HLT_PFMET250_HBHECleaned_v2;
    TBranch* b_HLT_PFMET250_HBHECleaned_v3;
    TBranch* b_HLT_PFMET250_HBHECleaned_v4;
    TBranch* b_HLT_PFMET250_HBHECleaned_v5;
    TBranch* b_HLT_PFMET250_HBHECleaned_v6;
    TBranch* b_HLT_PFMET250_HBHECleaned_v7;
    TBranch* b_HLT_PFMET300_HBHECleaned_v2;
    TBranch* b_HLT_PFMET300_HBHECleaned_v3;
    TBranch* b_HLT_PFMET300_HBHECleaned_v4;
    TBranch* b_HLT_PFMET300_HBHECleaned_v5;
    TBranch* b_HLT_PFMET300_HBHECleaned_v6;
    TBranch* b_HLT_PFMET300_HBHECleaned_v7;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v1;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v2;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_L1ETMnoHF_v10;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v2;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v3;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v4;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v5;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v6;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v7;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v11;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v13;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v14;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v15;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v16;
    TBranch* b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v9;
    TBranch* b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v10;
    TBranch* b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v12;
    TBranch* b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v13;
    TBranch* b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v14;
    TBranch* b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v15;
    TBranch* b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v9;
    TBranch* b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v10;
    TBranch* b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v12;
    TBranch* b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v13;
    TBranch* b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v14;
    TBranch* b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v15;
    TBranch* b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v16;
    TBranch* b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v17;
    TBranch* b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v9;
    TBranch* b_HLT_PFHT1050_v11;
    TBranch* b_HLT_PFHT1050_v12;
    TBranch* b_HLT_PFHT1050_v13;
    TBranch* b_HLT_PFHT1050_v14;
    TBranch* b_HLT_PFHT1050_v15;
    TBranch* b_HLT_PFHT1050_v16;
    TBranch* b_HLT_PFHT1050_v7;
    TBranch* b_HLT_PFHT1050_v9;
    TBranch* b_HLT_PFHT180_v7;
    TBranch* b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v1;
    TBranch* b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v10;
    TBranch* b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v3;
    TBranch* b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v5;
    TBranch* b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v6;
    TBranch* b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v7;
    TBranch* b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v8;
    TBranch* b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v9;
    TBranch* b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v1;
    TBranch* b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v10;
    TBranch* b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v3;
    TBranch* b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v5;
    TBranch* b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v6;
    TBranch* b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v7;
    TBranch* b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v8;
    TBranch* b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v9;
    TBranch* b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v1;
    TBranch* b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v10;
    TBranch* b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v3;
    TBranch* b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v5;
    TBranch* b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v6;
    TBranch* b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v7;
    TBranch* b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v8;
    TBranch* b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v9;
    TBranch* b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v1;
    TBranch* b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v10;
    TBranch* b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v3;
    TBranch* b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v5;
    TBranch* b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v6;
    TBranch* b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v7;
    TBranch* b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v8;
    TBranch* b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v9;
    TBranch* b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v1;
    TBranch* b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v10;
    TBranch* b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v3;
    TBranch* b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v5;
    TBranch* b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v6;
    TBranch* b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v7;
    TBranch* b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v8;
    TBranch* b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v9;
    TBranch* b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v1;
    TBranch* b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v10;
    TBranch* b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v3;
    TBranch* b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v5;
    TBranch* b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v6;
    TBranch* b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v7;
    TBranch* b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v8;
    TBranch* b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v9;

    TBranch* b_Flag_HBHENoiseFilter;
    TBranch* b_Flag_HBHENoiseIsoFilter;
    TBranch* b_Flag_globalTightHalo2016Filter;
    TBranch* b_Flag_EcalDeadCellTriggerPrimitiveFilter;
    TBranch* b_Flag_goodVertices;
    TBranch* b_Flag_eeBadScFilter;
    TBranch* b_Flag_ecalLaserCorrFilter;
    TBranch* b_Flag_chargedHadronTrackResolutionFilter;
    TBranch* b_Flag_muonBadTrackFilter;
    TBranch* b_Flag_badMuons;
    TBranch* b_Flag_duplicateMuons;
    TBranch* b_Flag_noBadMuons;
    TBranch* b_Flag_BadPFMuonFilter;
    TBranch* b_Flag_BadChargedCandidateFilter;
    TBranch* b_Flag_ecalBadCalibFilter;

    TBranch* b_nGenPar; //!
    TBranch* b_genParEta; //!
    TBranch* b_genParPhi; //!
    TBranch* b_genParE; //!
    TBranch* b_genParPt; //!
    TBranch* b_genParId; //!
    TBranch* b_genParVx; //!
    TBranch* b_genParVy; //!
    TBranch* b_genParVz; //!
    TBranch* b_genParMotherId; //!
    TBranch* b_genParMotherIndex; //!
    TBranch* b_genParNumDaughters; //!
    TBranch* b_genParDaughterId1; //!
    TBranch* b_genParDaughterId2; //!
    TBranch* b_genParStatus; //!
    TBranch* b_genParCharge; //!
    TBranch* b_eventRun; //!
    TBranch* b_eventNum; //!
    TBranch* b_eventLumiblock; //!

    bool isMC_{};
    const bool is2016_{};

    std::vector<double> muonMomentumSF;
    std::vector<double> jetSmearValue;

    std::vector<int> electronIndexTight;
    std::vector<int> electronIndexLoose;
    std::vector<int> muonIndexTight;
    std::vector<int> muonIndexLoose;
    std::vector<int> photonIndexTight;
    std::vector<int> photonIndexLoose;
    std::vector<int> jetIndex;
    std::vector<int> displacedJetIndex;
    double totalJetHt;
    double selectedJetIncHt;
    double selectedJetTracksHt;
    std::vector<int> bTagIndex;

    std::pair<TLorentzVector, TLorentzVector> zPairLeptons;
    std::pair<float, float> zPairRelIso;
    std::pair<int, int> zPairIndex;
    bool muonLeads{};

    std::pair<TLorentzVector, TLorentzVector> wPairQuarks;
    std::pair<int, int> wPairIndex;

    TLorentzVector wLepton;
    int wLepIndex;
    float wLeptonRelIso;

    Int_t numVert;
    TBranch* b_numVert;

    AnalysisEvent(bool isMC, TTree* tree, bool is2016);
    virtual ~AnalysisEvent();
    Int_t GetEntry(const Long64_t entry);
    Long64_t LoadTree(const Long64_t entry);
    void Loop();
    void Show(const Long64_t entry = -1) const;
//    bool eTrig() const;
    bool muTrig() const;
//    bool eeTrig() const;
//    bool muEGTrig() const;
    bool mumuTrig() const;

//    bool displacedJetTrig() const;
//    bool photonTrig() const;
//    bool diphotonTrig() const;
//    bool triphotonTrig() const;
//    bool photonHETrig() const;
//    bool diphotonHETrig() const;
//    bool diphotonHEpixelVetoTrig() const;
//    bool diphotonHEnoPixelVetoTrig() const;
//    bool photonDisplacedIdTrig() const;
//    bool photonHoverETrig() const;

};

inline AnalysisEvent::AnalysisEvent(const bool isMC,
                                    TTree* tree,
                                    const bool is2016)
    : fChain{nullptr}, isMC_{isMC}, is2016_{is2016}
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == nullptr)
    {
#ifdef SINGLE_TREE
        // The following code should be used if you want this class to access
        // a single tree instead of a chain
        TFile* f = (TFile*)gROOT->GetListOfFiles()->FindObject(
            "/data1/tW2012/mc/ttbarInclusive/MC_Ntuple_out_9_0_MJP_skim.root");
        if (!f || !f->IsOpen())
        {
            f = new TFile("/data1/tW2012/mc/ttbarInclusive/"
                          "MC_Ntuple_out_9_0_MJP_skim.root");
        }
        f->GetObject("tree", tree);

#else // SINGLE_TREE

        // The following code should be used if you want this class to access a
        // chain of trees.
        TChain* chain = new TChain("tree", "");
        chain->Add("/data1/tW2012/mc/ttbarInclusive/"
                   "MC_Ntuple_out_100_0_Gu6_skim.root/tree");
        tree = chain;
#endif // SINGLE_TREE
    }

    // clang-format off
   if (!tree)
   {
       return;
   }
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("numElePF2PAT", &numElePF2PAT, &b_numElePF2PAT);
   fChain->SetBranchAddress("elePF2PATE", elePF2PATE, &b_elePF2PATE);
   fChain->SetBranchAddress("elePF2PATET", elePF2PATET, &b_elePF2PATET);
   fChain->SetBranchAddress("elePF2PATPX", elePF2PATPX, &b_elePF2PATPX);
   fChain->SetBranchAddress("elePF2PATPY", elePF2PATPY, &b_elePF2PATPY);
   fChain->SetBranchAddress("elePF2PATPZ", elePF2PATPZ, &b_elePF2PATPZ);
   fChain->SetBranchAddress("elePF2PATPhi", elePF2PATPhi, &b_elePF2PATPhi);
   fChain->SetBranchAddress("elePF2PATTheta", elePF2PATTheta, &b_elePF2PATTheta);
   fChain->SetBranchAddress("elePF2PATEta", elePF2PATEta, &b_elePF2PATEta);
   fChain->SetBranchAddress("elePF2PATPT", elePF2PATPT, &b_elePF2PATPT);
   fChain->SetBranchAddress("elePF2PATCharge", elePF2PATCharge, &b_elePF2PATCharge);
   fChain->SetBranchAddress("elePF2PATCutIdVeto", elePF2PATCutIdVeto, &b_elePF2PATCutIdVeto);
   fChain->SetBranchAddress("elePF2PATCutIdLoose", elePF2PATCutIdLoose, &b_elePF2PATCutIdLoose);
   fChain->SetBranchAddress("elePF2PATCutIdMedium", elePF2PATCutIdMedium, &b_elePF2PATCutIdMedium);
   fChain->SetBranchAddress("elePF2PATCutIdTight", elePF2PATCutIdTight, &b_elePF2PATCutIdTight);
   fChain->SetBranchAddress("elePF2PATImpactTransDist", elePF2PATImpactTransDist, &b_elePF2PATImpactTransDist);
   fChain->SetBranchAddress("elePF2PATImpactTransError", elePF2PATImpactTransError, &b_elePF2PATImpactTransError);
   fChain->SetBranchAddress("elePF2PATImpactTransSignificance", elePF2PATImpactTransSignificance, &b_elePF2PATImpactTransSignificance);
   fChain->SetBranchAddress("elePF2PATImpact3DDist", elePF2PATImpact3DDist, &b_elePF2PATImpact3DDist);
   fChain->SetBranchAddress("elePF2PATImpact3DError", elePF2PATImpact3DError, &b_elePF2PATImpact3DError);
   fChain->SetBranchAddress("elePF2PATImpact3DSignificance", elePF2PATImpact3DSignificance, &b_elePF2PATImpact3DSignificance);
   fChain->SetBranchAddress("elePF2PATChargedHadronIso", elePF2PATChargedHadronIso, &b_elePF2PATChargedHadronIso);
   fChain->SetBranchAddress("elePF2PATNeutralHadronIso", elePF2PATNeutralHadronIso, &b_elePF2PATNeutralHadronIso);
   fChain->SetBranchAddress("elePF2PATPhotonIso", elePF2PATPhotonIso, &b_elePF2PATPhotonIso);
   fChain->SetBranchAddress("elePF2PATTrackPt", elePF2PATTrackPt, &b_elePF2PATTrackPt);
   fChain->SetBranchAddress("elePF2PATTrackPhi", elePF2PATTrackPhi, &b_elePF2PATTrackPhi);
   fChain->SetBranchAddress("elePF2PATTrackEta", elePF2PATTrackEta, &b_elePF2PATTrackEta);
   fChain->SetBranchAddress("elePF2PATTrackChi2", elePF2PATTrackChi2, &b_elePF2PATTrackChi2);
   fChain->SetBranchAddress("elePF2PATTrackNDOF", elePF2PATTrackNDOF, &b_elePF2PATTrackNDOF);
   fChain->SetBranchAddress("elePF2PATTrackD0", elePF2PATTrackD0, &b_elePF2PATTrackD0);
   fChain->SetBranchAddress("elePF2PATTrackDBD0", elePF2PATTrackDBD0, &b_elePF2PATTrackDBD0);
   fChain->SetBranchAddress("elePF2PATD0PV", elePF2PATD0PV, &b_elePF2PATD0PV);
   fChain->SetBranchAddress("elePF2PATDZPV", elePF2PATDZPV, &b_elePF2PATDZPV);
   fChain->SetBranchAddress("elePF2PATBeamSpotCorrectedTrackD0", elePF2PATBeamSpotCorrectedTrackD0, &b_elePF2PATBeamSpotCorrectedTrackD0);
   fChain->SetBranchAddress("elePF2PATTrackDz", elePF2PATTrackDz, &b_elePF2PATTrackDz);
   fChain->SetBranchAddress("elePF2PATVtxZ", elePF2PATVtxZ, &b_elePF2PATVtxZ);
   fChain->SetBranchAddress("elePF2PATIsGsf", elePF2PATIsGsf, &b_elePF2PATIsGsf);
   fChain->SetBranchAddress("elePF2PATGsfPx", elePF2PATGsfPx, &b_elePF2PATGsfPx);
   fChain->SetBranchAddress("elePF2PATGsfPy", elePF2PATGsfPy, &b_elePF2PATGsfPy);
   fChain->SetBranchAddress("elePF2PATGsfPz", elePF2PATGsfPz, &b_elePF2PATGsfPz);
   fChain->SetBranchAddress("elePF2PATGsfE", elePF2PATGsfE, &b_elePF2PATGsfE);
   fChain->SetBranchAddress("elePF2PATEcalEnergy", elePF2PATEcalEnergy, &b_elePF2PATEcalEnergy);
   fChain->SetBranchAddress("elePF2PATSCEta", elePF2PATSCEta, &b_elePF2PATSCEta);
   fChain->SetBranchAddress("elePF2PATSCE", elePF2PATSCE, &b_elePF2PATSCE);
   fChain->SetBranchAddress("elePF2PATSCPhi", elePF2PATSCPhi, &b_elePF2PATSCPhi);
   fChain->SetBranchAddress("elePF2PATSCEoverP", elePF2PATSCEoverP, &b_elePF2PATSCEoverP);
   fChain->SetBranchAddress("elePF2PATSCSigmaEtaEta", elePF2PATSCSigmaEtaEta, &b_elePF2PATSCSigmaEtaEta);
   fChain->SetBranchAddress("elePF2PATSCSigmaIEtaIEta", elePF2PATSCSigmaIEtaIEta, &b_elePF2PATSCSigmaIEtaIEta);
   fChain->SetBranchAddress("elePF2PATSCSigmaIEtaIEta5x5", elePF2PATSCSigmaIEtaIEta5x5, &b_elePF2PATSCSigmaIEtaIEta5x5);
   fChain->SetBranchAddress("elePF2PATSCE1x5", elePF2PATSCE1x5, &b_elePF2PATSCE1x5);
   fChain->SetBranchAddress("elePF2PATSCE5x5", elePF2PATSCE5x5, &b_elePF2PATSCE5x5);
   fChain->SetBranchAddress("elePF2PATSCE2x5max", elePF2PATSCE2x5max, &b_elePF2PATSCE2x5max);
   fChain->SetBranchAddress("elePF2PATTrackIso04", elePF2PATTrackIso04, &b_elePF2PATTrackIso04);
   fChain->SetBranchAddress("elePF2PATEcalIso04", elePF2PATEcalIso04, &b_elePF2PATEcalIso04);
   fChain->SetBranchAddress("elePF2PATHcalIso04", elePF2PATHcalIso04, &b_elePF2PATHcalIso04);
   fChain->SetBranchAddress("elePF2PATTrackIso03", elePF2PATTrackIso03, &b_elePF2PATTrackIso03);
   fChain->SetBranchAddress("elePF2PATEcalIso03", elePF2PATEcalIso03, &b_elePF2PATEcalIso03);
   fChain->SetBranchAddress("elePF2PATHcalIso03", elePF2PATHcalIso03, &b_elePF2PATHcalIso03);
   fChain->SetBranchAddress("elePF2PATdr04EcalRecHitSumEt", elePF2PATdr04EcalRecHitSumEt, &b_elePF2PATdr04EcalRecHitSumEt);
   fChain->SetBranchAddress("elePF2PATdr03EcalRecHitSumEt", elePF2PATdr03EcalRecHitSumEt, &b_elePF2PATdr03EcalRecHitSumEt);
   fChain->SetBranchAddress("elePF2PATEcalIsoDeposit", elePF2PATEcalIsoDeposit, &b_elePF2PATEcalIsoDeposit);
   fChain->SetBranchAddress("elePF2PATHcalIsoDeposit", elePF2PATHcalIsoDeposit, &b_elePF2PATHcalIsoDeposit);
   fChain->SetBranchAddress("elePF2PATComRelIso", elePF2PATComRelIso, &b_elePF2PATComRelIso);
   fChain->SetBranchAddress("elePF2PATComRelIsodBeta", elePF2PATComRelIsodBeta, &b_elePF2PATComRelIsodBeta);
   fChain->SetBranchAddress("elePF2PATComRelIsoRho", elePF2PATComRelIsoRho, &b_elePF2PATComRelIsoRho);
   fChain->SetBranchAddress("elePF2PATChHadIso", elePF2PATChHadIso, &b_elePF2PATChHadIso);
   fChain->SetBranchAddress("elePF2PATNtHadIso", elePF2PATNtHadIso, &b_elePF2PATNtHadIso);
   fChain->SetBranchAddress("elePF2PATGammaIso", elePF2PATGammaIso, &b_elePF2PATGammaIso);
   fChain->SetBranchAddress("elePF2PATRhoIso", elePF2PATRhoIso, &b_elePF2PATRhoIso);
   fChain->SetBranchAddress("elePF2PATAEff03", elePF2PATAEff03, &b_elePF2PATAEff03);
   fChain->SetBranchAddress("elePF2PATMissingInnerLayers", elePF2PATMissingInnerLayers, &b_elePF2PATMissingInnerLayers);
   fChain->SetBranchAddress("elePF2PATHoverE", elePF2PATHoverE, &b_elePF2PATHoverE);
   fChain->SetBranchAddress("elePF2PATDeltaPhiSC", elePF2PATDeltaPhiSC, &b_elePF2PATDeltaPhiSC);
   fChain->SetBranchAddress("elePF2PATDeltaEtaSC", elePF2PATDeltaEtaSC, &b_elePF2PATDeltaEtaSC);
   fChain->SetBranchAddress("elePF2PATDeltaEtaSeedSC", elePF2PATDeltaEtaSeedSC, &b_elePF2PATDeltaEtaSeedSC);
   fChain->SetBranchAddress("elePF2PATIsBarrel", elePF2PATIsBarrel, &b_elePF2PATIsBarrel);
   fChain->SetBranchAddress("elePF2PATPhotonConversionTag", elePF2PATPhotonConversionTag, &b_elePF2PATPhotonConversionTag);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDist", elePF2PATPhotonConversionDist, &b_elePF2PATPhotonConversionDist);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDcot", elePF2PATPhotonConversionDcot, &b_elePF2PATPhotonConversionDcot);
   fChain->SetBranchAddress("elePF2PATPhotonConversionVeto", elePF2PATPhotonConversionVeto, &b_elePF2PATPhotonConversionVeto);
//   fChain->SetBranchAddress("elePF2PATPhotonConversionTagCustom", elePF2PATPhotonConversionTagCustom, &b_elePF2PATPhotonConversionTagCustom);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDistCustom", elePF2PATPhotonConversionDistCustom, &b_elePF2PATPhotonConversionDistCustom);
   fChain->SetBranchAddress("elePF2PATPhotonConversionDcotCustom", elePF2PATPhotonConversionDcotCustom, &b_elePF2PATPhotonConversionDcotCustom);
   fChain->SetBranchAddress("elePF2PATTriggerMatch", elePF2PATTriggerMatch, &b_elePF2PATTriggerMatch);
   fChain->SetBranchAddress("elePF2PATJetOverlap", elePF2PATJetOverlap, &b_elePF2PATJetOverlap);
   if (isMC)
   {
       fChain->SetBranchAddress("genElePF2PATPT", genElePF2PATPT, &b_genElePF2PATPT);
       fChain->SetBranchAddress("genElePF2PATET", genElePF2PATET, &b_genElePF2PATET);
       fChain->SetBranchAddress("genElePF2PATPX", genElePF2PATPX, &b_genElePF2PATPX);
       fChain->SetBranchAddress("genElePF2PATPY", genElePF2PATPY, &b_genElePF2PATPY);
       fChain->SetBranchAddress("genElePF2PATPZ", genElePF2PATPZ, &b_genElePF2PATPZ);
       fChain->SetBranchAddress("genElePF2PATPhi", genElePF2PATPhi, &b_genElePF2PATPhi);
       fChain->SetBranchAddress("genElePF2PATTheta", genElePF2PATTheta, &b_genElePF2PATTheta);
       fChain->SetBranchAddress("genElePF2PATEta", genElePF2PATEta, &b_genElePF2PATEta);
       fChain->SetBranchAddress("genElePF2PATCharge", genElePF2PATCharge, &b_genElePF2PATCharge);
       fChain->SetBranchAddress("genElePF2PATPdgId", genElePF2PATPdgId, &b_genElePF2PATPdgId);
       fChain->SetBranchAddress("genElePF2PATMotherId", genElePF2PATMotherId, &b_genElePF2PATMotherId);
       fChain->SetBranchAddress("genElePF2PATPromptDecayed", genElePF2PATPromptDecayed, &b_genElePF2PATPromptDecayed);
       fChain->SetBranchAddress("genElePF2PATPromptFinalState", genElePF2PATPromptFinalState, &b_genElePF2PATPromptFinalState);
       fChain->SetBranchAddress("genElePF2PATHardProcess", genElePF2PATHardProcess, &b_genElePF2PATHardProcess);
   }
   fChain->SetBranchAddress("numMuonPF2PAT", &numMuonPF2PAT, &b_numMuonPF2PAT);
   fChain->SetBranchAddress("muonPF2PATE", muonPF2PATE, &b_muonPF2PATE);
   fChain->SetBranchAddress("muonPF2PATET", muonPF2PATET, &b_muonPF2PATET);
   fChain->SetBranchAddress("muonPF2PATPt", muonPF2PATPt, &b_muonPF2PATPt);
   fChain->SetBranchAddress("muonPF2PATPX", muonPF2PATPX, &b_muonPF2PATPX);
   fChain->SetBranchAddress("muonPF2PATPY", muonPF2PATPY, &b_muonPF2PATPY);
   fChain->SetBranchAddress("muonPF2PATPZ", muonPF2PATPZ, &b_muonPF2PATPZ);
   fChain->SetBranchAddress("muonPF2PATPhi", muonPF2PATPhi, &b_muonPF2PATPhi);
   fChain->SetBranchAddress("muonPF2PATTheta", muonPF2PATTheta, &b_muonPF2PATTheta);
   fChain->SetBranchAddress("muonPF2PATEta", muonPF2PATEta, &b_muonPF2PATEta);
   fChain->SetBranchAddress("muonPF2PATCharge", muonPF2PATCharge, &b_muonPF2PATCharge);
   fChain->SetBranchAddress("muonPF2PATLooseCutId", muonPF2PATLooseCutId, &b_muonPF2PATLooseCutId);
   fChain->SetBranchAddress("muonPF2PATMediumCutId", muonPF2PATMediumCutId, &b_muonPF2PATMediumCutId);
   fChain->SetBranchAddress("muonPF2PATTightCutId", muonPF2PATTightCutId, &b_muonPF2PATTightCutId);
   if (!is2016_)
   {
       fChain->SetBranchAddress("muonPF2PATPfIsoVeryLoose", muonPF2PATPfIsoVeryLoose, &b_muonPF2PATPfIsoVeryLoose);
       fChain->SetBranchAddress("muonPF2PATPfIsoLoose", muonPF2PATPfIsoLoose, &b_muonPF2PATPfIsoLoose);
       fChain->SetBranchAddress("muonPF2PATPfIsoMedium", muonPF2PATPfIsoMedium, &b_muonPF2PATPfIsoMedium);
       fChain->SetBranchAddress("muonPF2PATPfIsoTight", muonPF2PATPfIsoTight, &b_muonPF2PATPfIsoTight);
       fChain->SetBranchAddress("muonPF2PATPfIsoVeryTight", muonPF2PATPfIsoVeryTight, &b_muonPF2PATPfIsoVeryTight);
       fChain->SetBranchAddress("muonPF2PATTkIsoLoose", muonPF2PATTkIsoLoose, &b_muonPF2PATTkIsoLoose);
       fChain->SetBranchAddress("muonPF2PATTkIsoTight", muonPF2PATTkIsoTight, &b_muonPF2PATTkIsoTight);
       fChain->SetBranchAddress("muonPF2PATMvaLoose", muonPF2PATMvaLoose, &b_muonPF2PATMvaLoose);
       fChain->SetBranchAddress("muonPF2PATMvaMedium", muonPF2PATMvaMedium, &b_muonPF2PATMvaMedium);
       fChain->SetBranchAddress("muonPF2PATMvaTight", muonPF2PATMvaTight, &b_muonPF2PATMvaTight);
   }
   fChain->SetBranchAddress("muonPF2PATGlobalID", muonPF2PATGlobalID, &b_muonPF2PATGlobalID);
   fChain->SetBranchAddress("muonPF2PATTrackID", muonPF2PATTrackID, &b_muonPF2PATTrackID);
   fChain->SetBranchAddress("muonPF2PATChi2", muonPF2PATChi2, &b_muonPF2PATChi2);
   fChain->SetBranchAddress("muonPF2PATD0", muonPF2PATD0, &b_muonPF2PATD0);
   fChain->SetBranchAddress("muonPF2PATTrackDBD0", muonPF2PATTrackDBD0, &b_muonPF2PATTrackDBD0);
   fChain->SetBranchAddress("muonPF2PATDBInnerTrackD0", muonPF2PATDBInnerTrackD0, &b_muonPF2PATDBInnerTrackD0);
   fChain->SetBranchAddress("muonPF2PATBeamSpotCorrectedD0", muonPF2PATBeamSpotCorrectedD0, &b_muonPF2PATBeamSpotCorrectedD0);
   fChain->SetBranchAddress("muonPF2PATTrackNHits", muonPF2PATTrackNHits, &b_muonPF2PATTrackNHits);
   fChain->SetBranchAddress("muonPF2PATMuonNHits", muonPF2PATMuonNHits, &b_muonPF2PATMuonNHits);
   fChain->SetBranchAddress("muonPF2PATNDOF", muonPF2PATNDOF, &b_muonPF2PATNDOF);
   fChain->SetBranchAddress("muonPF2PATVertX", muonPF2PATVertX, &b_muonPF2PATVertX);
   fChain->SetBranchAddress("muonPF2PATVertY", muonPF2PATVertY, &b_muonPF2PATVertY);
   fChain->SetBranchAddress("muonPF2PATVertZ", muonPF2PATVertZ, &b_muonPF2PATVertZ);
   fChain->SetBranchAddress("muonPF2PATChargedHadronIso", muonPF2PATChargedHadronIso, &b_muonPF2PATChargedHadronIso);
   fChain->SetBranchAddress("muonPF2PATNeutralHadronIso", muonPF2PATNeutralHadronIso, &b_muonPF2PATNeutralHadronIso);
   fChain->SetBranchAddress("muonPF2PATPhotonIso", muonPF2PATPhotonIso, &b_muonPF2PATPhotonIso);
   fChain->SetBranchAddress("muonPF2PATTrackIso", muonPF2PATTrackIso, &b_muonPF2PATTrackIso);
   fChain->SetBranchAddress("muonPF2PATEcalIso", muonPF2PATEcalIso, &b_muonPF2PATEcalIso);
   fChain->SetBranchAddress("muonPF2PATHcalIso", muonPF2PATHcalIso, &b_muonPF2PATHcalIso);
   fChain->SetBranchAddress("muonPF2PATComRelIso", muonPF2PATComRelIso, &b_muonPF2PATComRelIso);
   fChain->SetBranchAddress("muonPF2PATComRelIsodBeta", muonPF2PATComRelIsodBeta, &b_muonPF2PATComRelIsodBeta);
   fChain->SetBranchAddress("muonPF2PATIsPFMuon", muonPF2PATIsPFMuon, &b_muonPF2PATIsPFMuon);
   fChain->SetBranchAddress("muonPF2PATNChambers", muonPF2PATNChambers, &b_muonPF2PATNChambers);
   fChain->SetBranchAddress("muonPF2PATNMatches", muonPF2PATNMatches, &b_muonPF2PATNMatches);
   fChain->SetBranchAddress("muonPF2PATTkLysWithMeasurements", muonPF2PATTkLysWithMeasurements, &b_muonPF2PATTkLysWithMeasurements);
   fChain->SetBranchAddress("muonPF2PATGlbTkNormChi2", muonPF2PATGlbTkNormChi2, &b_muonPF2PATGlbTkNormChi2);
   fChain->SetBranchAddress("muonPF2PATValidFraction", muonPF2PATValidFraction, &b_muonPF2PATValidFraction);
   fChain->SetBranchAddress("muonPF2PATChi2LocalPosition", muonPF2PATChi2LocalPosition, &b_muonPF2PATChi2LocalPosition);
   fChain->SetBranchAddress("muonPF2PATTrkKick", muonPF2PATTrkKick, &b_muonPF2PATTrkKick);
   fChain->SetBranchAddress("muonPF2PATSegmentCompatibility", muonPF2PATSegmentCompatibility, &b_muonPF2PATSegmentCompatibility);
   fChain->SetBranchAddress("muonPF2PATDBPV", muonPF2PATDBPV, &b_muonPF2PATDBPV);
   fChain->SetBranchAddress("muonPF2PATDZPV", muonPF2PATDZPV, &b_muonPF2PATDZPV);
   fChain->SetBranchAddress("muonPF2PATVldPixHits", muonPF2PATVldPixHits, &b_muonPF2PATVldPixHits);
   fChain->SetBranchAddress("muonPF2PATMatchedStations", muonPF2PATMatchedStations, &b_muonPF2PATMatchedStations);
   if (isMC)
   {
       fChain->SetBranchAddress("genMuonPF2PATPT", genMuonPF2PATPT, &b_genMuonPF2PATPT);
       fChain->SetBranchAddress("genMuonPF2PATET", genMuonPF2PATET, &b_genMuonPF2PATET);
       fChain->SetBranchAddress("genMuonPF2PATPX", genMuonPF2PATPX, &b_genMuonPF2PATPX);
       fChain->SetBranchAddress("genMuonPF2PATPY", genMuonPF2PATPY, &b_genMuonPF2PATPY);
       fChain->SetBranchAddress("genMuonPF2PATPZ", genMuonPF2PATPZ, &b_genMuonPF2PATPZ);
       fChain->SetBranchAddress("genMuonPF2PATPhi", genMuonPF2PATPhi, &b_genMuonPF2PATPhi);
       fChain->SetBranchAddress("genMuonPF2PATTheta", genMuonPF2PATTheta, &b_genMuonPF2PATTheta);
       fChain->SetBranchAddress("genMuonPF2PATEta", genMuonPF2PATEta, &b_genMuonPF2PATEta);
       fChain->SetBranchAddress("genMuonPF2PATCharge", genMuonPF2PATCharge, &b_genMuonPF2PATCharge);
       fChain->SetBranchAddress("genMuonPF2PATPdgId", genMuonPF2PATPdgId, &b_genMuonPF2PATPdgId);
       fChain->SetBranchAddress("genMuonPF2PATMotherId", genMuonPF2PATMotherId, &b_genMuonPF2PATMotherId);
       fChain->SetBranchAddress("genMuonPF2PATPromptDecayed", genMuonPF2PATPromptDecayed, &b_genMuonPF2PATPromptDecayed);
       fChain->SetBranchAddress("genMuonPF2PATPromptFinalState", genMuonPF2PATPromptFinalState, &b_genMuonPF2PATPromptFinalState);
       fChain->SetBranchAddress("genMuonPF2PATHardProcess", genMuonPF2PATHardProcess, &b_genMuonPF2PATHardProcess);
   }
   fChain->SetBranchAddress("numJetPF2PAT", &numJetPF2PAT, &b_numJetPF2PAT);
   fChain->SetBranchAddress("jetPF2PATE", jetPF2PATE, &b_jetPF2PATE);
   fChain->SetBranchAddress("jetPF2PATEt", jetPF2PATEt, &b_jetPF2PATEt);
   fChain->SetBranchAddress("jetPF2PATPt", jetPF2PATPt, &b_jetPF2PATPt);
   fChain->SetBranchAddress("jetPF2PATPtRaw", jetPF2PATPtRaw, &b_jetPF2PATPtRaw);
   fChain->SetBranchAddress("jetPF2PATUnCorEt", jetPF2PATUnCorEt, &b_jetPF2PATUnCorEt);
   fChain->SetBranchAddress("jetPF2PATUnCorPt", jetPF2PATUnCorPt, &b_jetPF2PATUnCorPt);
   fChain->SetBranchAddress("jetPF2PATEta", jetPF2PATEta, &b_jetPF2PATEta);
   fChain->SetBranchAddress("jetPF2PATTheta", jetPF2PATTheta, &b_jetPF2PATTheta);
   fChain->SetBranchAddress("jetPF2PATPhi", jetPF2PATPhi, &b_jetPF2PATPhi);
   fChain->SetBranchAddress("jetPF2PATPx", jetPF2PATPx, &b_jetPF2PATPx);
   fChain->SetBranchAddress("jetPF2PATPy", jetPF2PATPy, &b_jetPF2PATPy);
   fChain->SetBranchAddress("jetPF2PATPz", jetPF2PATPz, &b_jetPF2PATPz);
   fChain->SetBranchAddress("jetPF2PATdRClosestLepton", jetPF2PATdRClosestLepton, &b_jetPF2PATdRClosestLepton);
   fChain->SetBranchAddress("jetPF2PATNtracksInJet", jetPF2PATNtracksInJet, &b_jetPF2PATNtracksInJet);
   fChain->SetBranchAddress("jetPF2PATJetCharge", jetPF2PATJetCharge, &b_jetPF2PATJetCharge);
   fChain->SetBranchAddress("jetPF2PATfHPD", jetPF2PATfHPD, &b_jetPF2PATfHPD);
   fChain->SetBranchAddress("jetPF2PATBtagSoftMuonPtRel", jetPF2PATBtagSoftMuonPtRel, &b_jetPF2PATBtagSoftMuonPtRel);
   fChain->SetBranchAddress("jetPF2PATBtagSoftMuonQuality", jetPF2PATBtagSoftMuonQuality, &b_jetPF2PATBtagSoftMuonQuality);
   fChain->SetBranchAddress("jetPF2PATCorrFactor", jetPF2PATCorrFactor, &b_jetPF2PATCorrFactor);
   fChain->SetBranchAddress("jetPF2PATCorrResidual", jetPF2PATCorrResidual, &b_jetPF2PATCorrResidual);
   fChain->SetBranchAddress("jetPF2PATL2L3ResErr", jetPF2PATL2L3ResErr, &b_jetPF2PATL2L3ResErr);
   fChain->SetBranchAddress("jetPF2PATCorrErrLow", jetPF2PATCorrErrLow, &b_jetPF2PATCorrErrLow);
   fChain->SetBranchAddress("jetPF2PATCorrErrHi", jetPF2PATCorrErrHi, &b_jetPF2PATCorrErrHi);
   fChain->SetBranchAddress("jetPF2PATN90Hits", jetPF2PATN90Hits, &b_jetPF2PATN90Hits);
   fChain->SetBranchAddress("jetPF2PATTriggered", jetPF2PATTriggered, &b_jetPF2PATTriggered);
   fChain->SetBranchAddress("jetPF2PATSVX", jetPF2PATSVX, &b_jetPF2PATSVX);
   fChain->SetBranchAddress("jetPF2PATSVY", jetPF2PATSVY, &b_jetPF2PATSVY);
   fChain->SetBranchAddress("jetPF2PATSVZ", jetPF2PATSVZ, &b_jetPF2PATSVZ);
   fChain->SetBranchAddress("jetPF2PATSVDX", jetPF2PATSVDX, &b_jetPF2PATSVDX);
   fChain->SetBranchAddress("jetPF2PATSVDY", jetPF2PATSVDY, &b_jetPF2PATSVDY);
   fChain->SetBranchAddress("jetPF2PATSVDZ", jetPF2PATSVDZ, &b_jetPF2PATSVDZ);
   if (is2016)
   {
       fChain->SetBranchAddress("jetPF2PATBDiscriminator", jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags, &b_jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags);
   }
   else
   {
       fChain->SetBranchAddress("jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags", jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags, &b_jetPF2PATpfCombinedInclusiveSecondaryVertexV2BJetTags);
       fChain->SetBranchAddress("jetPF2PATpfCombinedCvsLJetTags", jetPF2PATpfCombinedCvsLJetTags, &b_jetPF2PATpfCombinedCvsLJetTags);
       fChain->SetBranchAddress("jetPF2PATpfCombinedCvsBJetTags", jetPF2PATpfCombinedCvsBJetTags, &b_jetPF2PATpfCombinedCvsBJetTags);
   }
   fChain->SetBranchAddress("jetPF2PATNConstituents", jetPF2PATNConstituents, &b_jetPF2PATNConstituents);
   fChain->SetBranchAddress("jetPF2PATPID", jetPF2PATPID, &b_jetPF2PATPID);
   fChain->SetBranchAddress("jetPF2PATClosestBPartonDeltaR", jetPF2PATClosestBPartonDeltaR, &b_jetPF2PATClosestBPartonDeltaR);
   fChain->SetBranchAddress("jetPF2PATClosestCPartonDeltaR", jetPF2PATClosestCPartonDeltaR, &b_jetPF2PATClosestCPartonDeltaR);
   if (isMC)
   {
       fChain->SetBranchAddress("genJetPF2PATE", genJetPF2PATE, &b_genJetPF2PATE);
       fChain->SetBranchAddress("genJetPF2PATET", genJetPF2PATET, &b_genJetPF2PATET);
       fChain->SetBranchAddress("genJetPF2PATPT", genJetPF2PATPT, &b_genJetPF2PATPT);
       fChain->SetBranchAddress("genJetPF2PATPX", genJetPF2PATPX, &b_genJetPF2PATPX);
       fChain->SetBranchAddress("genJetPF2PATPY", genJetPF2PATPY, &b_genJetPF2PATPY);
       fChain->SetBranchAddress("genJetPF2PATPZ", genJetPF2PATPZ, &b_genJetPF2PATPZ);
       fChain->SetBranchAddress("genJetPF2PATMass", genJetPF2PATMass, &b_genJetPF2PATMass);
       fChain->SetBranchAddress("genJetPF2PATPhi", genJetPF2PATPhi, &b_genJetPF2PATPhi);
       fChain->SetBranchAddress("genJetPF2PATTheta", genJetPF2PATTheta, &b_genJetPF2PATTheta);
       fChain->SetBranchAddress("genJetPF2PATEta", genJetPF2PATEta, &b_genJetPF2PATEta);
       fChain->SetBranchAddress("genJetPF2PATPID", genJetPF2PATPID, &b_genJetPF2PATPID);
       fChain->SetBranchAddress("genJetPF2PATMotherPID", genJetPF2PATMotherPID, &b_genJetPF2PATMotherPID);
       fChain->SetBranchAddress("genJetPF2PATScalarAncestor", genJetPF2PATScalarAncestor, &b_genJetPF2PATScalarAncestor);
   }
   fChain->SetBranchAddress("jetPF2PATMuEnergy", jetPF2PATMuEnergy, &b_jetPF2PATMuEnergy);
   fChain->SetBranchAddress("jetPF2PATMuEnergyFraction", jetPF2PATMuEnergyFraction, &b_jetPF2PATMuEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATNeutralHadEnergy", jetPF2PATNeutralHadEnergy, &b_jetPF2PATNeutralHadEnergy);
   fChain->SetBranchAddress("jetPF2PATNeutralEmEnergy", jetPF2PATNeutralEmEnergy, &b_jetPF2PATNeutralEmEnergy);
   fChain->SetBranchAddress("jetPF2PATChargedHadronEnergyFraction", jetPF2PATChargedHadronEnergyFraction, &b_jetPF2PATChargedHadronEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATNeutralHadronEnergyFraction", jetPF2PATNeutralHadronEnergyFraction, &b_jetPF2PATNeutralHadronEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATChargedEmEnergyFraction", jetPF2PATChargedEmEnergyFraction, &b_jetPF2PATChargedEmEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATNeutralEmEnergyFraction", jetPF2PATNeutralEmEnergyFraction, &b_jetPF2PATNeutralEmEnergyFraction);
   fChain->SetBranchAddress("jetPF2PATMuonFraction", jetPF2PATMuonFraction, &b_jetPF2PATMuonFraction);
   fChain->SetBranchAddress("jetPF2PATChargedHadronEnergyFractionCorr", jetPF2PATChargedHadronEnergyFractionCorr, &b_jetPF2PATChargedHadronEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATNeutralHadronEnergyFractionCorr", jetPF2PATNeutralHadronEnergyFractionCorr, &b_jetPF2PATNeutralHadronEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATChargedEmEnergyFractionCorr", jetPF2PATChargedEmEnergyFractionCorr, &b_jetPF2PATChargedEmEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATNeutralEmEnergyFractionCorr", jetPF2PATNeutralEmEnergyFractionCorr, &b_jetPF2PATNeutralEmEnergyFractionCorr);
   fChain->SetBranchAddress("jetPF2PATMuonFractionCorr", jetPF2PATMuonFractionCorr, &b_jetPF2PATMuonFractionCorr);
   fChain->SetBranchAddress("jetPF2PATNeutralMultiplicity", jetPF2PATNeutralMultiplicity, &b_jetPF2PATNeutralMultiplicity);
   fChain->SetBranchAddress("jetPF2PATChargedMultiplicity", jetPF2PATChargedMultiplicity, &b_jetPF2PATChargedMultiplicity);
   if (!is2016)
   {
       fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   }
   fChain->SetBranchAddress("jetPF2PATChargedMultiplicity", jetPF2PATChargedMultiplicity, &b_jetPF2PATChargedMultiplicity);
   fChain->SetBranchAddress("metPF2PATE", &metPF2PATE, &b_metPF2PATE);
   fChain->SetBranchAddress("metPF2PATEt", &metPF2PATEt, &b_metPF2PATEt);
   fChain->SetBranchAddress("metPF2PATEtRaw", &metPF2PATEtRaw, &b_metPF2PATEtRaw);
   fChain->SetBranchAddress("metPF2PATPhi", &metPF2PATPhi, &b_metPF2PATPhi);
   fChain->SetBranchAddress("metPF2PATPt", &metPF2PATPt, &b_metPF2PATPt);
   fChain->SetBranchAddress("metPF2PATPx", &metPF2PATPx, &b_metPF2PATPx);
   fChain->SetBranchAddress("metPF2PATPy", &metPF2PATPy, &b_metPF2PATPy);
   fChain->SetBranchAddress("metPF2PATPz", &metPF2PATPz, &b_metPF2PATPz);
   fChain->SetBranchAddress("metPF2PATScalarEt", &metPF2PATScalarEt, &b_metPF2PATScalarEt);
   fChain->SetBranchAddress("metPF2PATEtUncorrected", &metPF2PATEtUncorrected, &b_metPF2PATEtUncorrected);
   fChain->SetBranchAddress("metPF2PATPhiUncorrected", &metPF2PATPhiUncorrected, &b_metPF2PATPhiUncorrected);
   fChain->SetBranchAddress("metPF2PATUnclusteredEnUp", &metPF2PATUnclusteredEnUp, &b_metPF2PATUnclusteredEnUp);
   fChain->SetBranchAddress("metPF2PATUnclusteredEnDown", &metPF2PATUnclusteredEnDown, &b_metPF2PATUnclusteredEnDown);
   if (isMC)
   {
       fChain->SetBranchAddress("genMetPF2PATE", &genMetPF2PATE, &b_genMetPF2PATE);
       fChain->SetBranchAddress("genMetPF2PATEt", &genMetPF2PATEt, &b_genMetPF2PATEt);
       fChain->SetBranchAddress("genMetPF2PATPhi", &genMetPF2PATPhi, &b_genMetPF2PATPhi);
       fChain->SetBranchAddress("genMetPF2PATPt", &genMetPF2PATPt, &b_genMetPF2PATPt);
       fChain->SetBranchAddress("genMetPF2PATPx", &genMetPF2PATPx, &b_genMetPF2PATPx);
       fChain->SetBranchAddress("genMetPF2PATPy", &genMetPF2PATPy, &b_genMetPF2PATPy);
       fChain->SetBranchAddress("genMetPF2PATPz", &genMetPF2PATPz, &b_genMetPF2PATPz);
   }
   fChain->SetBranchAddress("numTauPF2PAT", &numTauPF2PAT, &b_numTauPF2PAT);
   fChain->SetBranchAddress("tauPF2PATE", tauPF2PATE, &b_tauPF2PATE);
   fChain->SetBranchAddress("tauPF2PATPt", tauPF2PATPt, &b_tauPF2PATPt);
   fChain->SetBranchAddress("tauPF2PATPhi", tauPF2PATPhi, &b_tauPF2PATPhi);
   fChain->SetBranchAddress("tauPF2PATEta", tauPF2PATEta, &b_tauPF2PATEta);

   fChain->SetBranchAddress("numPhoPF2PAT", &numPhoPF2PAT, &b_numPhoPF2PAT);
   fChain->SetBranchAddress("phoPF2PATE", phoPF2PATE, &b_phoPF2PATE);
   fChain->SetBranchAddress("phoPF2PATSigmaE", phoPF2PATSigmaE, &b_phoPF2PATSigmaE);
   fChain->SetBranchAddress("phoPF2PATET", phoPF2PATET, &b_phoPF2PATET);
   fChain->SetBranchAddress("phoPF2PATPhi", phoPF2PATPhi, &b_phoPF2PATPhi);
   fChain->SetBranchAddress("phoPF2PATEta", phoPF2PATEta, &b_phoPF2PATEta);
   fChain->SetBranchAddress("phoPF2PATTheta", phoPF2PATTheta, &b_phoPF2PATTheta);
   fChain->SetBranchAddress("phoPF2PATPt", phoPF2PATPt, &b_phoPF2PATPt);
   fChain->SetBranchAddress("phoPF2PATPx", phoPF2PATPx, &b_phoPF2PATPx);
   fChain->SetBranchAddress("phoPF2PATPy", phoPF2PATPy, &b_phoPF2PATPy);
   fChain->SetBranchAddress("phoPF2PATPz", phoPF2PATPz, &b_phoPF2PATPz);
   fChain->SetBranchAddress("phoPF2PATCalibE", phoPF2PATCalibE, &b_phoPF2PATCalibE);
   fChain->SetBranchAddress("phoPF2PATCalibEt", phoPF2PATCalibEt, &b_phoPF2PATCalibEt);
   fChain->SetBranchAddress("phoPF2PATSCE", phoPF2PATSCE, &b_phoPF2PATSCE);
   fChain->SetBranchAddress("phoPF2PATSCRawE", phoPF2PATSCRawE, &b_phoPF2PATSCRawE);
   fChain->SetBranchAddress("phoPF2PATESEnP1", phoPF2PATESEnP1, &b_phoPF2PATESEnP1);
   fChain->SetBranchAddress("phoPF2PATESEnP2", phoPF2PATESEnP2, &b_phoPF2PATESEnP1);
   fChain->SetBranchAddress("phoPF2PATSCEta", phoPF2PATSCEta, &b_phoPF2PATSCEta);
   fChain->SetBranchAddress("phoPF2PATSCEtaWidth", phoPF2PATSCEtaWidth, &b_phoPF2PATSCEtaWidth);
   fChain->SetBranchAddress("phoPF2PATSCPhi", phoPF2PATSCPhi, &b_phoPF2PATSCPhi);
   fChain->SetBranchAddress("phoPF2PATSCPhiWidth", phoPF2PATSCPhiWidth, &b_phoPF2PATSCPhiWidth);
   fChain->SetBranchAddress("phoPF2PATSCBrem", phoPF2PATSCBrem, &b_phoPF2PATSCBrem);
   fChain->SetBranchAddress("phoPF2PATHasPixelSeed", phoPF2PATHasPixelSeed, &b_phoPF2PATHasPixelSeed);
   fChain->SetBranchAddress("phoPF2PATEleVeto", phoPF2PATEleVeto, &b_phoPF2PATEleVeto);
   fChain->SetBranchAddress("phoPF2PATR9", phoPF2PATR9, &b_phoPF2PATR9);
   fChain->SetBranchAddress("phoPF2PATHoverE", phoPF2PATHoverE, &b_phoPF2PATHoverE);
   fChain->SetBranchAddress("phoPF2PATESEffSigmaRR", phoPF2PATESEffSigmaRR, &b_phoPF2PATESEffSigmaRR);
   fChain->SetBranchAddress("phoPF2PATSigmaIEtaIEtaFull5x5", phoPF2PATSigmaIEtaIEtaFull5x5, &b_phoPF2PATSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoPF2PATSigmaIEtaIPhiFull5x5", phoPF2PATSigmaIEtaIPhiFull5x5, &b_phoPF2PATSigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("phoPF2PATSigmaIPhiIPhiFull5x5", phoPF2PATSigmaIPhiIPhiFull5x5, &b_phoPF2PATSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("phoPF2PATE2x2Full5x5", phoPF2PATE2x2Full5x5, &b_phoPF2PATE2x2Full5x5);
   fChain->SetBranchAddress("phoPF2PATE5x5Full5x5", phoPF2PATE5x5Full5x5, &b_phoPF2PATE5x5Full5x5);
   fChain->SetBranchAddress("phoPF2PATR9Full5x5", phoPF2PATR9Full5x5, &b_phoPF2PATR9Full5x5);
   fChain->SetBranchAddress("phoPF2PATPFChIso", phoPF2PATPFChIso, &b_phoPF2PATPFChIso);
   fChain->SetBranchAddress("phoPF2PATPFPhoIso", phoPF2PATPFPhoIso, &b_phoPF2PATPFPhoIso);
   fChain->SetBranchAddress("phoPF2PATPFNeuIso", phoPF2PATPFNeuIso, &b_phoPF2PATPFNeuIso);
   fChain->SetBranchAddress("phoPF2PATPFChWorstIso", phoPF2PATPFChWorstIso, &b_phoPF2PATPFChWorstIso);
   fChain->SetBranchAddress("phoPF2PATMIPTotEnergy", phoPF2PATMIPTotEnergy, &b_phoPF2PATMIPTotEnergy);
   fChain->SetBranchAddress("phoPF2PATCutIdLoose", phoPF2PATCutIdLoose, &b_phoPF2PATCutIdLoose);
   fChain->SetBranchAddress("phoPF2PATCutIdMedium", phoPF2PATCutIdMedium, &b_phoPF2PATCutIdMedium);
   fChain->SetBranchAddress("phoPF2PATCutIdTight", phoPF2PATCutIdTight, &b_phoPF2PATCutIdTight);
   fChain->SetBranchAddress("phoPF2PATMvaIdWp80", phoPF2PATMvaIdWp80, &b_phoPF2PATMvaIdWp80);
   fChain->SetBranchAddress("phoPF2PATMvaIdWp90", phoPF2PATMvaIdWp90, &b_phoPF2PATMvaIdWp90);
   if (isMC) {
      fChain->SetBranchAddress("genPhoPF2PATPt", genPhoPF2PATPt, &b_genPhoPF2PATPt);
      fChain->SetBranchAddress("genPhoPF2PATET", genPhoPF2PATET, &b_genPhoPF2PATET);
      fChain->SetBranchAddress("genPhoPF2PATEta", genPhoPF2PATEta, &b_genPhoPF2PATEta);
      fChain->SetBranchAddress("genPhoPF2PATTheta", genPhoPF2PATTheta, &b_genPhoPF2PATTheta);
      fChain->SetBranchAddress("genPhoPF2PATPhi", genPhoPF2PATPhi, &b_genPhoPF2PATPhi);
      fChain->SetBranchAddress("genPhoPF2PATPx", genPhoPF2PATPx, &b_genPhoPF2PATPx);
      fChain->SetBranchAddress("genPhoPF2PATPy", genPhoPF2PATPy, &b_genPhoPF2PATPy);
      fChain->SetBranchAddress("genPhoPF2PATPz", genPhoPF2PATPz, &b_genPhoPF2PATPz);
      fChain->SetBranchAddress("genPhoPF2PATCharge", genPhoPF2PATCharge, &b_genPhoPF2PATCharge);
      fChain->SetBranchAddress("genPhoPF2PATPdgId", genPhoPF2PATPdgId, &b_genPhoPF2PATPdgId);
      fChain->SetBranchAddress("genPhoPF2PATMotherId", genPhoPF2PATMotherId, &b_genPhoPF2PATMotherId);
      fChain->SetBranchAddress("genPhoPF2PATIsPhoton", genPhoPF2PATIsPhoton, &b_genPhoPF2PATIsPhoton);
      fChain->SetBranchAddress("genPhoPF2PATIsConvertedPhoton", genPhoPF2PATIsConvertedPhoton, &b_genPhoPF2PATIsConvertedPhoton);
      fChain->SetBranchAddress("genPhoPF2PATIsJet", genPhoPF2PATIsJet, &b_genPhoPF2PATIsJet);
   }
/*
   fChain->SetBranchAddress("numPhoOOT_PF2PAT", &numPhoOOT_PF2PAT, &b_numPhoOOT_PF2PAT);
   fChain->SetBranchAddress("phoOOT_PF2PATE", phoOOT_PF2PATE, &b_phoOOT_PF2PATE);
   fChain->SetBranchAddress("phoOOT_PF2PATSigmaE", phoOOT_PF2PATSigmaE, &b_phoOOT_PF2PATSigmaE);
   fChain->SetBranchAddress("phoOOT_PF2PATET", phoOOT_PF2PATET, &b_phoOOT_PF2PATET);
   fChain->SetBranchAddress("phoOOT_PF2PATPhi", phoOOT_PF2PATPhi, &b_phoOOT_PF2PATPhi);
   fChain->SetBranchAddress("phoOOT_PF2PATEta", phoOOT_PF2PATEta, &b_phoOOT_PF2PATEta);
   fChain->SetBranchAddress("phoOOT_PF2PATTheta", phoOOT_PF2PATTheta, &b_phoOOT_PF2PATTheta);
   fChain->SetBranchAddress("phoOOT_PF2PATPt", phoOOT_PF2PATPt, &b_phoOOT_PF2PATPt);
   fChain->SetBranchAddress("phoOOT_PF2PATPx", phoOOT_PF2PATPx, &b_phoOOT_PF2PATPx);
   fChain->SetBranchAddress("phoOOT_PF2PATPy", phoOOT_PF2PATPy, &b_phoOOT_PF2PATPy);
   fChain->SetBranchAddress("phoOOT_PF2PATPz", phoOOT_PF2PATPz, &b_phoOOT_PF2PATPz);
   fChain->SetBranchAddress("phoOOT_PF2PATCalibE", phoOOT_PF2PATCalibE, &b_phoOOT_PF2PATCalibE);
   fChain->SetBranchAddress("phoOOT_PF2PATCalibEt", phoOOT_PF2PATCalibEt, &b_phoOOT_PF2PATCalibEt);
   fChain->SetBranchAddress("phoOOT_PF2PATSCE", phoOOT_PF2PATSCE, &b_phoOOT_PF2PATSCE);
   fChain->SetBranchAddress("phoOOT_PF2PATSCRawE", phoOOT_PF2PATSCRawE, &b_phoOOT_PF2PATSCRawE);
   fChain->SetBranchAddress("phoOOT_PF2PATESEnP1", phoOOT_PF2PATESEnP1, &b_phoOOT_PF2PATESEnP1);
   fChain->SetBranchAddress("phoOOT_PF2PATESEnP2", phoOOT_PF2PATESEnP2, &b_phoOOT_PF2PATESEnP1);
   fChain->SetBranchAddress("phoOOT_PF2PATSCEta", phoOOT_PF2PATSCEta, &b_phoOOT_PF2PATSCEta);
   fChain->SetBranchAddress("phoOOT_PF2PATSCEtaWidth", phoOOT_PF2PATSCEtaWidth, &b_phoOOT_PF2PATSCEtaWidth);
   fChain->SetBranchAddress("phoOOT_PF2PATSCPhi", phoOOT_PF2PATSCPhi, &b_phoOOT_PF2PATSCPhi);
   fChain->SetBranchAddress("phoOOT_PF2PATSCPhiWidth", phoOOT_PF2PATSCPhiWidth, &b_phoOOT_PF2PATSCPhiWidth);
   fChain->SetBranchAddress("phoOOT_PF2PATSCBrem", phoOOT_PF2PATSCBrem, &b_phoOOT_PF2PATSCBrem);
   fChain->SetBranchAddress("phoOOT_PF2PATHasPixelSeed", phoOOT_PF2PATHasPixelSeed, &b_phoOOT_PF2PATHasPixelSeed);
   fChain->SetBranchAddress("phoOOT_PF2PATEleVeto", phoOOT_PF2PATEleVeto, &b_phoOOT_PF2PATEleVeto);
   fChain->SetBranchAddress("phoOOT_PF2PATR9", phoOOT_PF2PATR9, &b_phoOOT_PF2PATR9);
   fChain->SetBranchAddress("phoOOT_PF2PATHoverE", phoOOT_PF2PATHoverE, &b_phoOOT_PF2PATHoverE);
   fChain->SetBranchAddress("phoOOT_PF2PATESEffSigmaRR", phoOOT_PF2PATESEffSigmaRR, &b_phoOOT_PF2PATESEffSigmaRR);
   fChain->SetBranchAddress("phoOOT_PF2PATSigmaIEtaIEtaFull5x5", phoOOT_PF2PATSigmaIEtaIEtaFull5x5, &b_phoOOT_PF2PATSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoOOT_PF2PATSigmaIEtaIPhiFull5x5", phoOOT_PF2PATSigmaIEtaIPhiFull5x5, &b_phoOOT_PF2PATSigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("phoOOT_PF2PATSigmaIPhiIPhiFull5x5", phoOOT_PF2PATSigmaIPhiIPhiFull5x5, &b_phoOOT_PF2PATSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("phoOOT_PF2PATE2x2Full5x5", phoOOT_PF2PATE2x2Full5x5, &b_phoOOT_PF2PATE2x2Full5x5);
   fChain->SetBranchAddress("phoOOT_PF2PATE5x5Full5x5", phoOOT_PF2PATE5x5Full5x5, &b_phoOOT_PF2PATE5x5Full5x5);
   fChain->SetBranchAddress("phoOOT_PF2PATR9Full5x5", phoOOT_PF2PATR9Full5x5, &b_phoOOT_PF2PATR9Full5x5);
   fChain->SetBranchAddress("phoOOT_PF2PATPFChIso", phoOOT_PF2PATPFChIso, &b_phoOOT_PF2PATPFChIso);
   fChain->SetBranchAddress("phoOOT_PF2PATPFPhoIso", phoOOT_PF2PATPFPhoIso, &b_phoOOT_PF2PATPFPhoIso);
   fChain->SetBranchAddress("phoOOT_PF2PATPFNeuIso", phoOOT_PF2PATPFNeuIso, &b_phoOOT_PF2PATPFNeuIso);
   fChain->SetBranchAddress("phoOOT_PF2PATPFChWorstIso", phoOOT_PF2PATPFChWorstIso, &b_phoOOT_PF2PATPFChWorstIso);
   fChain->SetBranchAddress("phoOOT_PF2PATMIPTotEnergy", phoOOT_PF2PATMIPTotEnergy, &b_phoOOT_PF2PATMIPTotEnergy);
   fChain->SetBranchAddress("phoOOT_PF2PATCutIdLoose", phoOOT_PF2PATCutIdLoose, &b_phoOOT_PF2PATCutIdLoose);
   fChain->SetBranchAddress("phoOOT_PF2PATCutIdMedium", phoOOT_PF2PATCutIdMedium, &b_phoOOT_PF2PATCutIdMedium);
   fChain->SetBranchAddress("phoOOT_PF2PATCutIdTight", phoOOT_PF2PATCutIdTight, &b_phoOOT_PF2PATCutIdTight);
   fChain->SetBranchAddress("phoOOT_PF2PATMvaIdWp80", phoOOT_PF2PATMvaIdWp80, &b_phoOOT_PF2PATMvaIdWp80);
   fChain->SetBranchAddress("phoOOT_PF2PATMvaIdWp90", phoOOT_PF2PATMvaIdWp90, &b_phoOOT_PF2PATMvaIdWp90);
   if (isMC) {
      fChain->SetBranchAddress("genPhoOOT_PF2PATPt", genPhoOOT_PF2PATPt, &b_genPhoOOT_PF2PATPt);
      fChain->SetBranchAddress("genPhoOOT_PF2PATET", genPhoOOT_PF2PATET, &b_genPhoOOT_PF2PATET);
      fChain->SetBranchAddress("genPhoOOT_PF2PATEta", genPhoOOT_PF2PATEta, &b_genPhoOOT_PF2PATEta);
      fChain->SetBranchAddress("genPhoOOT_PF2PATTheta", genPhoOOT_PF2PATTheta, &b_genPhoOOT_PF2PATTheta);
      fChain->SetBranchAddress("genPhoOOT_PF2PATPhi", genPhoOOT_PF2PATPhi, &b_genPhoOOT_PF2PATPhi);
      fChain->SetBranchAddress("genPhoOOT_PF2PATPx", genPhoOOT_PF2PATPx, &b_genPhoOOT_PF2PATPx);
      fChain->SetBranchAddress("genPhoOOT_PF2PATPy", genPhoOOT_PF2PATPy, &b_genPhoOOT_PF2PATPy);
      fChain->SetBranchAddress("genPhoOOT_PF2PATPz", genPhoOOT_PF2PATPz, &b_genPhoOOT_PF2PATPz);
      fChain->SetBranchAddress("genPhoOOT_PF2PATCharge", genPhoOOT_PF2PATCharge, &b_genPhoOOT_PF2PATCharge);
      fChain->SetBranchAddress("genPhoOOT_PF2PATPdgId", genPhoOOT_PF2PATPdgId, &b_genPhoOOT_PF2PATPdgId);
      fChain->SetBranchAddress("genPhoOOT_PF2PATMotherId", genPhoOOT_PF2PATMotherId, &b_genPhoOOT_PF2PATMotherId);
      fChain->SetBranchAddress("genPhoOOT_PF2PATIsPhoton", genPhoOOT_PF2PATIsPhoton, &b_genPhoOOT_PF2PATIsPhoton);
      fChain->SetBranchAddress("genPhoOOT_PF2PATIsConvertedPhoton", genPhoOOT_PF2PATIsConvertedPhoton, &b_genPhoOOT_PF2PATIsConvertedPhoton);
      fChain->SetBranchAddress("genPhoOOT_PF2PATIsJet", genPhoOOT_PF2PATIsJet, &b_genPhoOOT_PF2PATIsJet);
   }
*/

   fChain->SetBranchAddress("numPackedCands", &numPackedCands, &b_numPackedCands);
//   fChain->SetBranchAddress("packedCandsPt", packedCandsPt, &b_packedCandsPt);
   fChain->SetBranchAddress("packedCandsPx", packedCandsPx, &b_packedCandsPx);
   fChain->SetBranchAddress("packedCandsPy", packedCandsPy, &b_packedCandsPy);
   fChain->SetBranchAddress("packedCandsPz", packedCandsPz, &b_packedCandsPz);
   fChain->SetBranchAddress("packedCandsE", packedCandsE, &b_packedCandsE);
//   fChain->SetBranchAddress("packedCandsEta", packedCandsEta, &b_packedCandsEta);
//   fChain->SetBranchAddress("packedCandsTheta", packedCandsTheta, &b_packedCandsTheta);
//   fChain->SetBranchAddress("packedCandsPhi", packedCandsPhi, &b_packedCandsPhi);
   fChain->SetBranchAddress("packedCandsCharge", packedCandsCharge, &b_packedCandsCharge);
   fChain->SetBranchAddress("packedCandsPdgId", packedCandsPdgId, &b_packedCandsPdgId);
//   fChain->SetBranchAddress("packedCandsTime", packedCandsTime, &b_packedCandsTime);
//   fChain->SetBranchAddress("packedCandsFromPV", packedCandsFromPV, &b_packedCandsFromPV);
//   fChain->SetBranchAddress("packedCandsPVquality", packedCandsPVquality, &b_packedCandsPVquality);
//   fChain->SetBranchAddress("packedCandsVx", packedCandsVx, &b_packedCandsVx);
//   fChain->SetBranchAddress("packedCandsVy", packedCandsVy, &b_packedCandsVy);
//   fChain->SetBranchAddress("packedCandsVz", packedCandsVz, &b_packedCandsVz);
//   fChain->SetBranchAddress("packedCandsVEta", packedCandsVEta, &b_packedCandsVEta);
//   fChain->SetBranchAddress("packedCandsVPhi", packedCandsVPhi, &b_packedCandsVPhi);
//   fChain->SetBranchAddress("packedCandsBeamSpotCorrectedD0", packedCandsBeamSpotCorrectedD0, &b_packedCandsBeamSpotCorrectedD0);
   fChain->SetBranchAddress("packedCandsDz", packedCandsDz, &b_packedCandsDz);
   fChain->SetBranchAddress("packedCandsDxy", packedCandsDxy, &b_packedCandsDxy);
//   fChain->SetBranchAddress("packedCandsDzAssocPV", packedCandsDzAssocPV, &b_packedCandsDzAssocPV);
//   fChain->SetBranchAddress("packedCandsVtxChi2Norm", packedCandsVtxChi2Norm, &b_packedCandsVtxChi2Norm);
   fChain->SetBranchAddress("packedCandsHasTrackDetails", packedCandsHasTrackDetails, &b_packedCandsHasTrackDetails);
//   fChain->SetBranchAddress("packedCandsDzError", packedCandsDzError, &b_packedCandsDzError);
//   fChain->SetBranchAddress("packedCandsDxyError", packedCandsDxyError, &b_packedCandsDxyError);
//   fChain->SetBranchAddress("packedCandsTimeError", packedCandsTimeError, &b_packedCandsTimeError);
//   fChain->SetBranchAddress("packedCandsPseudoTrkPt", packedCandsPseudoTrkPt, &b_packedCandsPseudoTrkPt);
   fChain->SetBranchAddress("packedCandsPseudoTrkPx", packedCandsPseudoTrkPx, &b_packedCandsPseudoTrkPx);
   fChain->SetBranchAddress("packedCandsPseudoTrkPy", packedCandsPseudoTrkPy, &b_packedCandsPseudoTrkPy);
   fChain->SetBranchAddress("packedCandsPseudoTrkPz", packedCandsPseudoTrkPz, &b_packedCandsPseudoTrkPz);
//   fChain->SetBranchAddress("packedCandsPseudoTrkEta", packedCandsPseudoTrkEta, &b_packedCandsPseudoTrkEta);
//   fChain->SetBranchAddress("packedCandsPseudoTrkPhi", packedCandsPseudoTrkPhi, &b_packedCandsPseudoTrkPhi);
   fChain->SetBranchAddress("packedCandsPseudoTrkCharge", packedCandsPseudoTrkCharge, &b_packedCandsPseudoTrkCharge);
   fChain->SetBranchAddress("packedCandsPseudoTrkVx", packedCandsPseudoTrkVx, &b_packedCandsPseudoTrkVx);
   fChain->SetBranchAddress("packedCandsPseudoTrkVy", packedCandsPseudoTrkVy, &b_packedCandsPseudoTrkVy);
   fChain->SetBranchAddress("packedCandsPseudoTrkVz", packedCandsPseudoTrkVz, &b_packedCandsPseudoTrkVz);
//   fChain->SetBranchAddress("packedCandsPseudoTrkChi2Norm", packedCandsPseudoTrkChi2Norm, &b_packedCandsPseudoTrkChi2Norm);
//   fChain->SetBranchAddress("packedCandsPseudoTrkNumberOfHits", packedCandsPseudoTrkNumberOfHits, &b_packedCandsPseudoTrkNumberOfHits);
//   fChain->SetBranchAddress("packedCandsPseudoTrkNumberOfPixelHits", packedCandsPseudoTrkNumberOfPixelHits, &b_packedCandsPseudoTrkNumberOfHits);
//   fChain->SetBranchAddress("packedCandsPseudoTrkPixelLayersWithMeasurement", packedCandsPseudoTrkPixelLayersWithMeasurement, &b_packedCandsPseudoTrkPixelLayersWithMeasurement);
//   fChain->SetBranchAddress("packedCandsPseudoTrkStripLayersWithMeasurement", packedCandsPseudoTrkStripLayersWithMeasurement, &b_packedCandsPseudoTrkStripLayersWithMeasurement);
//   fChain->SetBranchAddress("packedCandsPseudoTrkTrackerLayersWithMeasurement", packedCandsPseudoTrkTrackerLayersWithMeasurement, &b_packedCandsPseudoTrkTrackerLayersWithMeasurement);
   fChain->SetBranchAddress("packedCandsHighPurityTrack", packedCandsHighPurityTrack, &b_packedCandsHighPurityTrack);

   fChain->SetBranchAddress("numIsolatedTracks", &numIsolatedTracks, &b_numIsolatedTracks);
   fChain->SetBranchAddress("isoTracksPt", isoTracksPt, &b_isoTracksPt);
   fChain->SetBranchAddress("isoTracksPx", isoTracksPx, &b_isoTracksPx);
   fChain->SetBranchAddress("isoTracksPy", isoTracksPy, &b_isoTracksPy);
   fChain->SetBranchAddress("isoTracksPz", isoTracksPz, &b_isoTracksPz);
   fChain->SetBranchAddress("isoTracksE", isoTracksE, &b_isoTracksE);
   fChain->SetBranchAddress("isoTracksEta", isoTracksEta, &b_isoTracksEta);
   fChain->SetBranchAddress("isoTracksTheta", isoTracksTheta, &b_isoTracksTheta);
   fChain->SetBranchAddress("isoTracksPhi", isoTracksPhi, &b_isoTracksPhi);
   fChain->SetBranchAddress("isoTracksCharge", isoTracksCharge, &b_isoTracksCharge);
   fChain->SetBranchAddress("isoTracksPdgId", isoTracksPdgId, &b_isoTracksPdgId);
   fChain->SetBranchAddress("isoTracksMatchedCaloJetEmEnergy", isoTracksMatchedCaloJetEmEnergy, &b_isoTracksMatchedCaloJetEmEnergy);
   fChain->SetBranchAddress("isoTracksMatchedCaloJetHadEnergy", isoTracksMatchedCaloJetHadEnergy, &b_isoTracksMatchedCaloJetHadEnergy);
   fChain->SetBranchAddress("isoTracksDz", isoTracksDz, &b_isoTracksDz);
   fChain->SetBranchAddress("isoTracksDxy", isoTracksDxy, &b_isoTracksDxy);
   fChain->SetBranchAddress("isoTracksDzError", isoTracksDzError, &b_isoTracksDzError);
   fChain->SetBranchAddress("isoTracksDxyError", isoTracksDxyError, &b_isoTracksDxyError);
   fChain->SetBranchAddress("isoTracksFromPV", isoTracksFromPV, &b_isoTracksFromPV);
   fChain->SetBranchAddress("isoTracksVx", isoTracksVx, &b_isoTracksVx);
   fChain->SetBranchAddress("isoTracksVy", isoTracksVy, &b_isoTracksVy);
   fChain->SetBranchAddress("isoTracksVz", isoTracksVz, &b_isoTracksVz);
   fChain->SetBranchAddress("isoTracksHighPurity", isoTracksHighPurity, &b_isoTracksHighPurity);
   fChain->SetBranchAddress("isoTracksTight", isoTracksTight, &b_isoTracksTight);
   fChain->SetBranchAddress("isoTracksLoose", isoTracksLoose, &b_isoTracksLoose);
   fChain->SetBranchAddress("isoTracksDeltaEta", isoTracksDeltaEta, &b_isoTracksDeltaEta);
   fChain->SetBranchAddress("isoTracksDeltaPhi", isoTracksDeltaPhi, &b_isoTracksDeltaPhi);

   if (isMC)
   {
       fChain->SetBranchAddress("isElePlusJets", &isElePlusJets, &b_isElePlusJets);
       fChain->SetBranchAddress("genPDFScale", &genPDFScale, &b_genPDFScale);
       fChain->SetBranchAddress("genPDFx1", &genPDFx1, &b_genPDFx1);
       fChain->SetBranchAddress("genPDFx2", &genPDFx2, &b_genPDFx2);
       fChain->SetBranchAddress("genPDFf1", &genPDFf1, &b_genPDFf1);
       fChain->SetBranchAddress("genPDFf2", &genPDFf2, &b_genPDFf2);
       fChain->SetBranchAddress("topPtReweight", &topPtReweight, &b_topPtReweight);
   }
   fChain->SetBranchAddress("processId", &processId, &b_processId);
   fChain->SetBranchAddress("processPtHat", &processPtHat, &b_processPtHat);
   fChain->SetBranchAddress("processMCWeight", &processMCWeight, &b_processMCWeight);
   fChain->SetBranchAddress("beamSpotX", &beamSpotX, &b_beamSpotX);
   fChain->SetBranchAddress("beamSpotY", &beamSpotY, &b_beamSpotY);
   fChain->SetBranchAddress("beamSpotZ", &beamSpotZ, &b_beamSpotZ);
   fChain->SetBranchAddress("numPVs", &numPVs, &b_numPVs);
   fChain->SetBranchAddress("pvX", pvX, &b_pvX);
   fChain->SetBranchAddress("pvY", pvY, &b_pvY);
   fChain->SetBranchAddress("pvZ", pvZ, &b_pvZ);
   fChain->SetBranchAddress("pvDX", pvDX, &b_pvDX);
   fChain->SetBranchAddress("pvDY", pvDY, &b_pvDY);
   fChain->SetBranchAddress("pvDZ", pvDZ, &b_pvDZ);
   fChain->SetBranchAddress("pvRho", pvRho, &b_pvRho);
   fChain->SetBranchAddress("pvIsFake", pvIsFake, &b_pvIsFake);
   fChain->SetBranchAddress("pvNdof", pvNdof, &b_pvNdof);
   fChain->SetBranchAddress("pvChi2", pvChi2, &b_pvChi2);
   fChain->SetBranchAddress("numSVs", &numSVs, &b_numSVs);
   fChain->SetBranchAddress("svPt", svPt, &b_svPt);
   fChain->SetBranchAddress("svPx", svPx, &b_svPx);
   fChain->SetBranchAddress("svPy", svPy, &b_svPy);
   fChain->SetBranchAddress("svPz", svPz, &b_svPz);
   fChain->SetBranchAddress("svMass", svMass, &b_svMass);
   fChain->SetBranchAddress("svE", svE, &b_svE);
   fChain->SetBranchAddress("svEta", svEta, &b_svEta);
   fChain->SetBranchAddress("svTheta", svTheta, &b_svTheta);
   fChain->SetBranchAddress("svPhi", svPhi, &b_svPhi);
   fChain->SetBranchAddress("svX", svX, &b_svX);
   fChain->SetBranchAddress("svY", svY, &b_svY);
   fChain->SetBranchAddress("svZ", svZ, &b_svZ);
   fChain->SetBranchAddress("svVertexChi2", svVertexChi2, &b_svVertexChi2);
   fChain->SetBranchAddress("svVertexNdof", svVertexNdof, &b_svVertexNdof);
   fChain->SetBranchAddress("svNtracks", svNtracks, &b_svNtracks);
   fChain->SetBranchAddress("svDist3D", svDist3D, &b_svDist3D);
   fChain->SetBranchAddress("svDist3DSig", svDist3DSig, &b_svDist3DSig);
   fChain->SetBranchAddress("svDist3DError", svDist3DError, &b_svDist3DError);
   fChain->SetBranchAddress("svDistXY", svDistXY, &b_svDistXY);
   fChain->SetBranchAddress("svDistXYSig", svDistXYSig, &b_svDistXYSig);
   fChain->SetBranchAddress("svDistXYError", svDistXYError, &b_svDistXYError);
   fChain->SetBranchAddress("svAnglePV", svAnglePV, &b_svAnglePV);
   fChain->SetBranchAddress("svIsLambda", svIsLambda, &b_svIsLambda);
   fChain->SetBranchAddress("svIsKshort", svIsKshort, &b_svIsKshort);
   fChain->SetBranchAddress("mhtPt", &mhtPt, &b_mhtPt);
   fChain->SetBranchAddress("mhtPy", &mhtPy, &b_mhtPy);
   fChain->SetBranchAddress("mhtPx", &mhtPx, &b_mhtPx);
   fChain->SetBranchAddress("mhtPhi", &mhtPhi, &b_mhtPhi);
   fChain->SetBranchAddress("mhtSumEt", &mhtSumEt, &b_mhtSumEt);
   fChain->SetBranchAddress("mhtSignif", &mhtSignif, &b_mhtSignif);
   fChain->SetBranchAddress("nTriggerBits", &nTriggerBits, &b_nTriggerBits);
   fChain->SetBranchAddress("TriggerBits", &TriggerBits, &b_TriggerBits);
   if (isMC)
   {
       fChain->SetBranchAddress("weight_muF0p5", &weight_muF0p5, &b_weight_muF0p5);
       fChain->SetBranchAddress("weight_muF2", &weight_muF2, &b_weight_muF2);
       fChain->SetBranchAddress("weight_muR0p5", &weight_muR0p5, &b_weight_muR0p5);
       fChain->SetBranchAddress("weight_muR2", &weight_muR2, &b_weight_muR2);
       fChain->SetBranchAddress("weight_muF0p5muR0p5", &weight_muF0p5muR0p5, &b_weight_muF0p5muR0p5);
       fChain->SetBranchAddress("weight_muF2muR2", &weight_muF2muR2, &b_weight_muF2muR2);
       fChain->SetBranchAddress("origWeightForNorm", &origWeightForNorm, &b_origWeightForNorm);
       fChain->SetBranchAddress("weight_pdfMax", &weight_pdfMax, &b_weight_pdfMax);
       fChain->SetBranchAddress("weight_pdfMin", &weight_pdfMin, &b_weight_pdfMin);
       fChain->SetBranchAddress("weight_alphaMax", &weight_alphaMax, &b_weight_alphaMax);
       fChain->SetBranchAddress("weight_alphaMin", &weight_alphaMin, &b_weight_alphaMin);
       if (!is2016)
       {
           fChain->SetBranchAddress("isrRedHi", &isrRedHi, &b_isrRedHi);
           fChain->SetBranchAddress("fsrRedHi", &fsrRedHi, &b_fsrRedHi);
           fChain->SetBranchAddress("isrRedLo", &isrRedLo, &b_isrRedLo);
           fChain->SetBranchAddress("fsrRedLo", &fsrRedLo, &b_fsrRedLo);
           fChain->SetBranchAddress("isrDefHi", &isrDefHi, &b_isrDefHi);
           fChain->SetBranchAddress("fsrDefHi", &fsrDefHi, &b_fsrDefHi);
           fChain->SetBranchAddress("isrDefLo", &isrDefLo, &b_isrDefLo);
           fChain->SetBranchAddress("fsrDefLo", &fsrDefLo, &b_fsrDefLo);
           fChain->SetBranchAddress("isrConHi", &isrConHi, &b_isrConHi);
           fChain->SetBranchAddress("fsrConHi", &fsrConHi, &b_fsrConHi);
           fChain->SetBranchAddress("isrConLo", &isrConLo, &b_isrConLo);
           fChain->SetBranchAddress("fsrConLo", &fsrConLo, &b_fsrConLo);
       }
   }
   //Trigger branches
   if (is2016)
   {
/*
       fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v1", &HLT_Ele25_eta2p1_WPTight_Gsf_v1, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v1);
       fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v2", &HLT_Ele25_eta2p1_WPTight_Gsf_v2, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v2);
       fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v3", &HLT_Ele25_eta2p1_WPTight_Gsf_v3, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v3);
       fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v4", &HLT_Ele25_eta2p1_WPTight_Gsf_v4, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v4);
       fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v5", &HLT_Ele25_eta2p1_WPTight_Gsf_v5, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v5);
       fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v6", &HLT_Ele25_eta2p1_WPTight_Gsf_v6, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v6);
       fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf_v7", &HLT_Ele25_eta2p1_WPTight_Gsf_v7, &b_HLT_Ele25_eta2p1_WPTight_Gsf_v7);
       fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v1", &HLT_Ele27_WPTight_Gsf_v1, &b_HLT_Ele27_WPTight_Gsf_v1);
       fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v2", &HLT_Ele27_WPTight_Gsf_v2, &b_HLT_Ele27_WPTight_Gsf_v2);
       fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v3", &HLT_Ele27_WPTight_Gsf_v3, &b_HLT_Ele27_WPTight_Gsf_v3);
       fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v4", &HLT_Ele27_WPTight_Gsf_v4, &b_HLT_Ele27_WPTight_Gsf_v4);
       fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v5", &HLT_Ele27_WPTight_Gsf_v5, &b_HLT_Ele27_WPTight_Gsf_v5);
       fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v6", &HLT_Ele27_WPTight_Gsf_v6, &b_HLT_Ele27_WPTight_Gsf_v6);
       fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v7", &HLT_Ele27_WPTight_Gsf_v7, &b_HLT_Ele27_WPTight_Gsf_v7);
       fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v2", &HLT_Ele32_eta2p1_WPTight_Gsf_v2, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v2);
       fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v3", &HLT_Ele32_eta2p1_WPTight_Gsf_v3, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v3);
       fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v4", &HLT_Ele32_eta2p1_WPTight_Gsf_v4, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v4);
       fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v5", &HLT_Ele32_eta2p1_WPTight_Gsf_v5, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v5);
       fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v6", &HLT_Ele32_eta2p1_WPTight_Gsf_v6, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v6);
       fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v7", &HLT_Ele32_eta2p1_WPTight_Gsf_v7, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v7);
       fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_v8", &HLT_Ele32_eta2p1_WPTight_Gsf_v8, &b_HLT_Ele32_eta2p1_WPTight_Gsf_v8);
*/
       fChain->SetBranchAddress("HLT_IsoMu24_v1", &HLT_IsoMu24_v1, &b_HLT_IsoMu24_v1);
       fChain->SetBranchAddress("HLT_IsoMu24_v2", &HLT_IsoMu24_v2, &b_HLT_IsoMu24_v2);
       fChain->SetBranchAddress("HLT_IsoMu24_v3", &HLT_IsoMu24_v3, &b_HLT_IsoMu24_v3);
       fChain->SetBranchAddress("HLT_IsoMu24_v4", &HLT_IsoMu24_v4, &b_HLT_IsoMu24_v4);
       fChain->SetBranchAddress("HLT_IsoTkMu24_v1", &HLT_IsoTkMu24_v1, &b_HLT_IsoTkMu24_v1);
       fChain->SetBranchAddress("HLT_IsoTkMu24_v2", &HLT_IsoTkMu24_v2, &b_HLT_IsoTkMu24_v2);
       fChain->SetBranchAddress("HLT_IsoTkMu24_v3", &HLT_IsoTkMu24_v3, &b_HLT_IsoTkMu24_v3);
       fChain->SetBranchAddress("HLT_IsoTkMu24_v4", &HLT_IsoTkMu24_v4, &b_HLT_IsoTkMu24_v4);
/*
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9);
*/
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4);
       // fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5);
       // fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v6", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v6, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v6);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3);
       // fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4);
       // fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6);
/*
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4);
*/
   }
   else
   {
/*
       fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v1", &HLT_Ele32_WPTight_Gsf_L1DoubleEG_v1, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v1);
       fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v2", &HLT_Ele32_WPTight_Gsf_L1DoubleEG_v2, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v2);
       fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v3", &HLT_Ele32_WPTight_Gsf_L1DoubleEG_v3, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v3);
       fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v4", &HLT_Ele32_WPTight_Gsf_L1DoubleEG_v4, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v4);
       fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v5", &HLT_Ele32_WPTight_Gsf_L1DoubleEG_v5, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v5);
       fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v6", &HLT_Ele32_WPTight_Gsf_L1DoubleEG_v6, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v6);
       fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7", &HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7);
       fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_v1", &HLT_Ele35_WPTight_Gsf_v1, &b_HLT_Ele35_WPTight_Gsf_v1);
       fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_v2", &HLT_Ele35_WPTight_Gsf_v2, &b_HLT_Ele35_WPTight_Gsf_v2);
       fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_v3", &HLT_Ele35_WPTight_Gsf_v3, &b_HLT_Ele35_WPTight_Gsf_v3);
       fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_v4", &HLT_Ele35_WPTight_Gsf_v4, &b_HLT_Ele35_WPTight_Gsf_v4);
       fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_v5", &HLT_Ele35_WPTight_Gsf_v5, &b_HLT_Ele35_WPTight_Gsf_v5);
       fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_v6", &HLT_Ele35_WPTight_Gsf_v6, &b_HLT_Ele35_WPTight_Gsf_v6);
       fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_v7", &HLT_Ele35_WPTight_Gsf_v7, &b_HLT_Ele35_WPTight_Gsf_v7);
*/
       fChain->SetBranchAddress("HLT_IsoMu27_v8", &HLT_IsoMu27_v8, &b_HLT_IsoMu27_v8);
       fChain->SetBranchAddress("HLT_IsoMu27_v9", &HLT_IsoMu27_v9, &b_HLT_IsoMu27_v9);
       fChain->SetBranchAddress("HLT_IsoMu27_v10", &HLT_IsoMu27_v10, &b_HLT_IsoMu27_v10);
       fChain->SetBranchAddress("HLT_IsoMu27_v11", &HLT_IsoMu27_v11, &b_HLT_IsoMu27_v11);
       fChain->SetBranchAddress("HLT_IsoMu27_v12", &HLT_IsoMu27_v12, &b_HLT_IsoMu27_v12);
       fChain->SetBranchAddress("HLT_IsoMu27_v13", &HLT_IsoMu27_v13, &b_HLT_IsoMu27_v13);
       fChain->SetBranchAddress("HLT_IsoMu27_v14", &HLT_IsoMu27_v14, &b_HLT_IsoMu27_v14);
/*
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v10", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v10, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v10);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v11", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v11, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v11);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v12", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v12, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v12);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v13", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v13, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v13);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v14", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v14, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v14);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v15", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v15, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v15);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v16", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v16, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v16);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v17", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v17, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v17);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v14", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v14, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v14);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v15", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v15, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v15);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v16", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v16, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v16);
       fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17);
*/
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v8);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v9", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v9, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v9);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v10", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v10, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v10);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v11", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v11, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v11);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v12", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v12, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v12);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v13", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v13, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v13);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v14", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v14, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v14);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v1", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v1, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v1);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v2", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v2, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v2);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v3", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v3, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v3);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v4", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v4, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v4);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v7", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v7, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v7);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v8", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v8, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v8);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v1", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v1, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v1);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v2", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v2, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v2);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v3", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v3, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v3);
       fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v4", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v4, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v4);
/*
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12);
       fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v5", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v5, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v5);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v12", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v12, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v12);
       fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v13", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v13, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v13);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v7", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v7, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v7);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10);
       fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11);
*/
   }
/*
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5);
*/

   // Displacced jets trigger branches 2017
   if ( !is2016 ) {
//       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack_v10", &HLT_HT430_DisplacedDijet40_DisplacedTrack_v10, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v10); // Already in MET trigger list
//       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack_v10", &HLT_HT430_DisplacedDijet60_DisplacedTrack_v10, &b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v10); // Already in MET trigger list
//       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet80_DisplacedTrack_v10", &HLT_HT430_DisplacedDijet80_DisplacedTrack_v10, &b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v10); // Already in MET trigger list
       fChain->SetBranchAddress("HLT_HT400_DisplacedDijet40_DisplacedTrack_v10", &HLT_HT400_DisplacedDijet40_DisplacedTrack_v10, &b_HLT_HT400_DisplacedDijet40_DisplacedTrack_v10);
//       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive_v10", &HLT_HT650_DisplacedDijet60_Inclusive_v10, &b_HLT_HT650_DisplacedDijet60_Inclusive_v10); // Already in MET trigger list
       fChain->SetBranchAddress("HLT_HT550_DisplacedDijet80_Inclusive_v8", &HLT_HT550_DisplacedDijet80_Inclusive_v8, &b_HLT_HT550_DisplacedDijet80_Inclusive_v8);
       fChain->SetBranchAddress("HLT_HT550_DisplacedDijet60_Inclusive_v10", &HLT_HT550_DisplacedDijet60_Inclusive_v10, &b_HLT_HT550_DisplacedDijet60_Inclusive_v10);
//       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive_v11", &HLT_HT650_DisplacedDijet80_Inclusive_v11, &b_HLT_HT650_DisplacedDijet80_Inclusive_v11); // Already in MET trigger list
//       fChain->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive_v11", &HLT_HT750_DisplacedDijet80_Inclusive_v11, &b_HLT_HT750_DisplacedDijet80_Inclusive_v11); // Already in MET trigger list
   }

   // Photon trigger branches 2017
   if ( !is2016 ) {
/*
       fChain->SetBranchAddress("HLT_DoublePhoton33_CaloIdL_v5", &HLT_DoublePhoton33_CaloIdL_v5, &b_HLT_DoublePhoton33_CaloIdL_v5);
       fChain->SetBranchAddress("HLT_DoublePhoton70_v5", &HLT_DoublePhoton70_v5, &b_HLT_DoublePhoton70_v5);
       fChain->SetBranchAddress("HLT_DoublePhoton85_v13", &HLT_DoublePhoton85_v13, &b_HLT_DoublePhoton85_v13);
       fChain->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2_v2", &HLT_TriplePhoton_20_20_20_CaloIdLV2_v2, &b_HLT_TriplePhoton_20_20_20_CaloIdLV2_v2);
       fChain->SetBranchAddress("HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_v2", &HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_v2, &b_HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_v2);
       fChain->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2_v3", &HLT_TriplePhoton_30_30_10_CaloIdLV2_v3, &b_HLT_TriplePhoton_30_30_10_CaloIdLV2_v3);
       fChain->SetBranchAddress("HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_v3", &HLT_TriplePhoton_30_30_10_CaloIdLV2_v3, &b_HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_v3);
       fChain->SetBranchAddress("HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_v3", &HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_v3, &b_HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_v3);
       fChain->SetBranchAddress("HLT_Photon25_v2", &HLT_Photon25_v2, &b_HLT_Photon25_v2);
       fChain->SetBranchAddress("HLT_Photon33_v4", &HLT_Photon33_v4, &b_HLT_Photon33_v4);
       fChain->SetBranchAddress("HLT_Photon50_v12", &HLT_Photon50_v12, &b_HLT_Photon50_v12);
       fChain->SetBranchAddress("HLT_Photon75_v12", &HLT_Photon75_v12, &b_HLT_Photon75_v12);
       fChain->SetBranchAddress("HLT_Photon90_v12", &HLT_Photon90_v12, &b_HLT_Photon90_v12);
       fChain->SetBranchAddress("HLT_Photon120_v12", &HLT_Photon120_v12, &b_HLT_Photon120_v12);
       fChain->SetBranchAddress("HLT_Photon150_v5", &HLT_Photon150_v5, &b_HLT_Photon150_v5);
       fChain->SetBranchAddress("HLT_Photon175_v13", &HLT_Photon175_v13, &b_HLT_Photon175_v13);
       fChain->SetBranchAddress("HLT_Photon200_v12", &HLT_Photon200_v12, &b_HLT_Photon200_v12);
       fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM_v13", &HLT_Photon50_R9Id90_HE10_IsoM_v13, &b_HLT_Photon50_R9Id90_HE10_IsoM_v13);
       fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM_v13", &HLT_Photon75_R9Id90_HE10_IsoM_v13, &b_HLT_Photon75_R9Id90_HE10_IsoM_v13);
       fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM_v13", &HLT_Photon90_R9Id90_HE10_IsoM_v13, &b_HLT_Photon90_R9Id90_HE10_IsoM_v13);
       fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM_v13", &HLT_Photon120_R9Id90_HE10_IsoM_v13, &b_HLT_Photon120_R9Id90_HE10_IsoM_v13);
       fChain->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM_v14", &HLT_Photon165_R9Id90_HE10_IsoM_v14, &b_HLT_Photon165_R9Id90_HE10_IsoM_v14);
       fChain->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT700_v12", &HLT_Photon90_CaloIdL_PFHT700_v12, &b_HLT_Photon90_CaloIdL_PFHT700_v12);

       fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v12", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v12, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v12);
       fChain->SetBranchAddress("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v12", &HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v12, &b_HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v12);
       fChain->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13, &b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13);
       fChain->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12, &b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12);
       fChain->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12, &b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_NoPixelVeto_Mass55_v12);
       fChain->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13, &b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v13);
   
       fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_v4", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_v4, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_v4);
       fChain->SetBranchAddress("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_v7", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_v7, &b_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_v7);
       fChain->SetBranchAddress("HLT_Photon20_HoverELoose_v9", &HLT_Photon20_HoverELoose_v9, &b_HLT_Photon20_HoverELoose_v9);
       fChain->SetBranchAddress("HLT_Photon30_HoverELoose_v9", &HLT_Photon30_HoverELoose_v9, &b_HLT_Photon30_HoverELoose_v9);
       fChain->SetBranchAddress("HLT_Photon40_HoverELoose_v9", &HLT_Photon40_HoverELoose_v9, &b_HLT_Photon40_HoverELoose_v9);
       fChain->SetBranchAddress("HLT_Photon50_HoverELoose_v9", &HLT_Photon50_HoverELoose_v9, &b_HLT_Photon50_HoverELoose_v9);
       fChain->SetBranchAddress("HLT_Photon60_HoverELoose_v9", &HLT_Photon60_HoverELoose_v9, &b_HLT_Photon60_HoverELoose_v9);
*/
   }

   //MET trigger branches
   if (is2016)
   {
       fChain->SetBranchAddress("HLT_MET200_v1", &HLT_MET200_v1, &b_HLT_MET200_v1);
       fChain->SetBranchAddress("HLT_MET200_v2", &HLT_MET200_v2, &b_HLT_MET200_v2);
       fChain->SetBranchAddress("HLT_MET200_v3", &HLT_MET200_v3, &b_HLT_MET200_v3);
       fChain->SetBranchAddress("HLT_MET200_v4", &HLT_MET200_v4, &b_HLT_MET200_v4);
       fChain->SetBranchAddress("HLT_MET200_v5", &HLT_MET200_v5, &b_HLT_MET200_v5);
       fChain->SetBranchAddress("HLT_MET250_v1", &HLT_MET250_v1, &b_HLT_MET250_v1);
       fChain->SetBranchAddress("HLT_MET250_v2", &HLT_MET250_v2, &b_HLT_MET250_v2);
       fChain->SetBranchAddress("HLT_MET250_v3", &HLT_MET250_v3, &b_HLT_MET250_v3);
       fChain->SetBranchAddress("HLT_MET250_v4", &HLT_MET250_v4, &b_HLT_MET250_v4);
       fChain->SetBranchAddress("HLT_MET250_v5", &HLT_MET250_v5, &b_HLT_MET250_v5);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v2", &HLT_PFMET120_PFMHT120_IDTight_v2, &b_HLT_PFMET120_PFMHT120_IDTight_v2);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v3", &HLT_PFMET120_PFMHT120_IDTight_v3, &b_HLT_PFMET120_PFMHT120_IDTight_v3);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v4", &HLT_PFMET120_PFMHT120_IDTight_v4, &b_HLT_PFMET120_PFMHT120_IDTight_v4);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v5", &HLT_PFMET120_PFMHT120_IDTight_v5, &b_HLT_PFMET120_PFMHT120_IDTight_v5);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v6", &HLT_PFMET120_PFMHT120_IDTight_v6, &b_HLT_PFMET120_PFMHT120_IDTight_v6);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v7", &HLT_PFMET120_PFMHT120_IDTight_v7, &b_HLT_PFMET120_PFMHT120_IDTight_v7);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v8", &HLT_PFMET120_PFMHT120_IDTight_v8, &b_HLT_PFMET120_PFMHT120_IDTight_v8);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v2", &HLT_PFMET170_HBHECleaned_v2, &b_HLT_PFMET170_HBHECleaned_v2);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v3", &HLT_PFMET170_HBHECleaned_v3, &b_HLT_PFMET170_HBHECleaned_v3);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v4", &HLT_PFMET170_HBHECleaned_v4, &b_HLT_PFMET170_HBHECleaned_v4);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v5", &HLT_PFMET170_HBHECleaned_v5, &b_HLT_PFMET170_HBHECleaned_v5);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v6", &HLT_PFMET170_HBHECleaned_v6, &b_HLT_PFMET170_HBHECleaned_v6);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v7", &HLT_PFMET170_HBHECleaned_v7, &b_HLT_PFMET170_HBHECleaned_v7);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v8", &HLT_PFMET170_HBHECleaned_v8, &b_HLT_PFMET170_HBHECleaned_v8);
       fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned_v9", &HLT_PFMET170_HBHECleaned_v9, &b_HLT_PFMET170_HBHECleaned_v9);
       fChain->SetBranchAddress("HLT_PFHT800_v2", &HLT_PFHT800_v2, &b_HLT_PFHT800_v2);
       fChain->SetBranchAddress("HLT_PFHT800_v3", &HLT_PFHT800_v3, &b_HLT_PFHT800_v3);
       fChain->SetBranchAddress("HLT_PFHT800_v4", &HLT_PFHT800_v4, &b_HLT_PFHT800_v4);
       fChain->SetBranchAddress("HLT_PFHT800_v5", &HLT_PFHT800_v5, &b_HLT_PFHT800_v5);
       fChain->SetBranchAddress("HLT_PFHT900_v4", &HLT_PFHT900_v4, &b_HLT_PFHT900_v4);
       fChain->SetBranchAddress("HLT_PFHT900_v5", &HLT_PFHT900_v5, &b_HLT_PFHT900_v5);
       fChain->SetBranchAddress("HLT_PFHT900_v6", &HLT_PFHT900_v6, &b_HLT_PFHT900_v6);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt50_v3", &HLT_PFHT750_4JetPt50_v3, &b_HLT_PFHT750_4JetPt50_v3);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt50_v4", &HLT_PFHT750_4JetPt50_v4, &b_HLT_PFHT750_4JetPt50_v4);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt50_v5", &HLT_PFHT750_4JetPt50_v5, &b_HLT_PFHT750_4JetPt50_v5);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt50_v6", &HLT_PFHT750_4JetPt50_v6, &b_HLT_PFHT750_4JetPt50_v6);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt70_v1", &HLT_PFHT750_4JetPt70_v1, &b_HLT_PFHT750_4JetPt70_v1);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt70_v2", &HLT_PFHT750_4JetPt70_v2, &b_HLT_PFHT750_4JetPt70_v2);
       fChain->SetBranchAddress("HLT_PFHT750_4JetPt80_v2", &HLT_PFHT750_4JetPt80_v2, &b_HLT_PFHT750_4JetPt80_v2);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET100_v1", &HLT_PFHT300_PFMET100_v1, &b_HLT_PFHT300_PFMET100_v1);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET100_v2", &HLT_PFHT300_PFMET100_v2, &b_HLT_PFHT300_PFMET100_v2);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET100_v3", &HLT_PFHT300_PFMET100_v3, &b_HLT_PFHT300_PFMET100_v3);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET100_v4", &HLT_PFHT300_PFMET100_v4, &b_HLT_PFHT300_PFMET100_v4);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET110_v4", &HLT_PFHT300_PFMET110_v4, &b_HLT_PFHT300_PFMET110_v4);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET110_v5", &HLT_PFHT300_PFMET110_v5, &b_HLT_PFHT300_PFMET110_v5);
       fChain->SetBranchAddress("HLT_PFHT300_PFMET110_v6", &HLT_PFHT300_PFMET110_v6, &b_HLT_PFHT300_PFMET110_v6);
   }
   else
   { // 2017
       fChain->SetBranchAddress("HLT_MET105_IsoTrk50_v1", &HLT_MET105_IsoTrk50_v1, &b_HLT_MET105_IsoTrk50_v1);
       fChain->SetBranchAddress("HLT_MET105_IsoTrk50_v3", &HLT_MET105_IsoTrk50_v3, &b_HLT_MET105_IsoTrk50_v3);
       fChain->SetBranchAddress("HLT_MET105_IsoTrk50_v4", &HLT_MET105_IsoTrk50_v4, &b_HLT_MET105_IsoTrk50_v4);
       fChain->SetBranchAddress("HLT_MET105_IsoTrk50_v5", &HLT_MET105_IsoTrk50_v5, &b_HLT_MET105_IsoTrk50_v5);
       fChain->SetBranchAddress("HLT_MET105_IsoTrk50_v6", &HLT_MET105_IsoTrk50_v6, &b_HLT_MET105_IsoTrk50_v6);
       fChain->SetBranchAddress("HLT_MET105_IsoTrk50_v7", &HLT_MET105_IsoTrk50_v7, &b_HLT_MET105_IsoTrk50_v7);
       fChain->SetBranchAddress("HLT_MET105_IsoTrk50_v8", &HLT_MET105_IsoTrk50_v8, &b_HLT_MET105_IsoTrk50_v8);
       fChain->SetBranchAddress("HLT_MET120_IsoTrk50_v1", &HLT_MET120_IsoTrk50_v1, &b_HLT_MET120_IsoTrk50_v1);
       fChain->SetBranchAddress("HLT_MET120_IsoTrk50_v3", &HLT_MET120_IsoTrk50_v3, &b_HLT_MET120_IsoTrk50_v3);
       fChain->SetBranchAddress("HLT_MET120_IsoTrk50_v4", &HLT_MET120_IsoTrk50_v4, &b_HLT_MET120_IsoTrk50_v4);
       fChain->SetBranchAddress("HLT_MET120_IsoTrk50_v5", &HLT_MET120_IsoTrk50_v5, &b_HLT_MET120_IsoTrk50_v5);
       fChain->SetBranchAddress("HLT_MET120_IsoTrk50_v6", &HLT_MET120_IsoTrk50_v6, &b_HLT_MET120_IsoTrk50_v6);
       fChain->SetBranchAddress("HLT_MET120_IsoTrk50_v7", &HLT_MET120_IsoTrk50_v7, &b_HLT_MET120_IsoTrk50_v7);
       fChain->SetBranchAddress("HLT_MET120_IsoTrk50_v8", &HLT_MET120_IsoTrk50_v8, &b_HLT_MET120_IsoTrk50_v8);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack_v10", &HLT_HT430_DisplacedDijet40_DisplacedTrack_v10, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v10); // Already in Displaced Jets trigger list
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack_v11", &HLT_HT430_DisplacedDijet40_DisplacedTrack_v11, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v11);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack_v5", &HLT_HT430_DisplacedDijet40_DisplacedTrack_v5, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v5);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack_v6", &HLT_HT430_DisplacedDijet40_DisplacedTrack_v6, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v6);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack_v8", &HLT_HT430_DisplacedDijet40_DisplacedTrack_v8, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v8);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet40_DisplacedTrack_v9", &HLT_HT430_DisplacedDijet40_DisplacedTrack_v9, &b_HLT_HT430_DisplacedDijet40_DisplacedTrack_v9);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack_v10", &HLT_HT430_DisplacedDijet60_DisplacedTrack_v10, &b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v10); // Already in Displaced Jets trigger list
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack_v11", &HLT_HT430_DisplacedDijet60_DisplacedTrack_v11, &b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v11);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack_v5", &HLT_HT430_DisplacedDijet60_DisplacedTrack_v5, &b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v5);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack_v6", &HLT_HT430_DisplacedDijet60_DisplacedTrack_v6, &b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v6);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack_v8", &HLT_HT430_DisplacedDijet60_DisplacedTrack_v8, &b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v8);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet60_DisplacedTrack_v9", &HLT_HT430_DisplacedDijet60_DisplacedTrack_v9, &b_HLT_HT430_DisplacedDijet60_DisplacedTrack_v9);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet80_DisplacedTrack_v10", &HLT_HT430_DisplacedDijet80_DisplacedTrack_v10, &b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v10); // Already in Displaced Jet trigger list
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet80_DisplacedTrack_v11", &HLT_HT430_DisplacedDijet80_DisplacedTrack_v11, &b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v11);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet80_DisplacedTrack_v5", &HLT_HT430_DisplacedDijet80_DisplacedTrack_v5, &b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v5);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet80_DisplacedTrack_v6", &HLT_HT430_DisplacedDijet80_DisplacedTrack_v6, &b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v6);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet80_DisplacedTrack_v8", &HLT_HT430_DisplacedDijet80_DisplacedTrack_v8, &b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v8);
       fChain->SetBranchAddress("HLT_HT430_DisplacedDijet80_DisplacedTrack_v9", &HLT_HT430_DisplacedDijet80_DisplacedTrack_v9, &b_HLT_HT430_DisplacedDijet80_DisplacedTrack_v9);
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive_v10", &HLT_HT650_DisplacedDijet60_Inclusive_v10, &b_HLT_HT650_DisplacedDijet60_Inclusive_v10); // Already in Displaced Jet trigger list
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive_v11", &HLT_HT650_DisplacedDijet60_Inclusive_v11, &b_HLT_HT650_DisplacedDijet60_Inclusive_v11);
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive_v5", &HLT_HT650_DisplacedDijet60_Inclusive_v5, &b_HLT_HT650_DisplacedDijet60_Inclusive_v5);
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive_v6", &HLT_HT650_DisplacedDijet60_Inclusive_v6, &b_HLT_HT650_DisplacedDijet60_Inclusive_v6);
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive_v8", &HLT_HT650_DisplacedDijet60_Inclusive_v8, &b_HLT_HT650_DisplacedDijet60_Inclusive_v8);
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet60_Inclusive_v9", &HLT_HT650_DisplacedDijet60_Inclusive_v9, &b_HLT_HT650_DisplacedDijet60_Inclusive_v9);
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive_v10", &HLT_HT650_DisplacedDijet80_Inclusive_v10, &b_HLT_HT650_DisplacedDijet80_Inclusive_v10);
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive_v11", &HLT_HT650_DisplacedDijet80_Inclusive_v11, &b_HLT_HT650_DisplacedDijet80_Inclusive_v11); // Already in Displaced Jet trigger list
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive_v12", &HLT_HT650_DisplacedDijet80_Inclusive_v12, &b_HLT_HT650_DisplacedDijet80_Inclusive_v12);
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive_v6", &HLT_HT650_DisplacedDijet80_Inclusive_v6, &b_HLT_HT650_DisplacedDijet80_Inclusive_v6);
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive_v7", &HLT_HT650_DisplacedDijet80_Inclusive_v7, &b_HLT_HT650_DisplacedDijet80_Inclusive_v7);
       fChain->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive_v9", &HLT_HT650_DisplacedDijet80_Inclusive_v9, &b_HLT_HT650_DisplacedDijet80_Inclusive_v9);
       fChain->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive_v10", &HLT_HT750_DisplacedDijet80_Inclusive_v10, &b_HLT_HT750_DisplacedDijet80_Inclusive_v10);
       fChain->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive_v11", &HLT_HT750_DisplacedDijet80_Inclusive_v11, &b_HLT_HT750_DisplacedDijet80_Inclusive_v11); // Already in Displaced Jet trigger list
       fChain->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive_v12", &HLT_HT750_DisplacedDijet80_Inclusive_v12, &b_HLT_HT750_DisplacedDijet80_Inclusive_v12);
       fChain->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive_v6", &HLT_HT750_DisplacedDijet80_Inclusive_v6, &b_HLT_HT750_DisplacedDijet80_Inclusive_v6);
       fChain->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive_v7", &HLT_HT750_DisplacedDijet80_Inclusive_v7, &b_HLT_HT750_DisplacedDijet80_Inclusive_v7);
       fChain->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive_v9", &HLT_HT750_DisplacedDijet80_Inclusive_v9, &b_HLT_HT750_DisplacedDijet80_Inclusive_v9);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v1", &HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v1, &b_HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v1);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v2", &HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v2, &b_HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v2);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_L1ETMnoHF_v10", &HLT_PFMET120_PFMHT120_IDTight_L1ETMnoHF_v10, &b_HLT_PFMET120_PFMHT120_IDTight_L1ETMnoHF_v10);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v1", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v1, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v1);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v2", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v2, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v2);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v2", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_v2, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v2);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v3", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_v3, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v3);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v4", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_v4, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v4);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v5", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_v5, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v5);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v6", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_v6, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v6);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v7", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_v7, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v7);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v11", &HLT_PFMET120_PFMHT120_IDTight_v11, &b_HLT_PFMET120_PFMHT120_IDTight_v11);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v13", &HLT_PFMET120_PFMHT120_IDTight_v13, &b_HLT_PFMET120_PFMHT120_IDTight_v13);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v14", &HLT_PFMET120_PFMHT120_IDTight_v14, &b_HLT_PFMET120_PFMHT120_IDTight_v14);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v15", &HLT_PFMET120_PFMHT120_IDTight_v15, &b_HLT_PFMET120_PFMHT120_IDTight_v15);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v16", &HLT_PFMET120_PFMHT120_IDTight_v16, &b_HLT_PFMET120_PFMHT120_IDTight_v16);
       fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v9", &HLT_PFMET120_PFMHT120_IDTight_v9, &b_HLT_PFMET120_PFMHT120_IDTight_v9);
       fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight_v11", &HLT_PFMET130_PFMHT130_IDTight_v11, &b_HLT_PFMET130_PFMHT130_IDTight_v11);
       fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight_v13", &HLT_PFMET130_PFMHT130_IDTight_v13, &b_HLT_PFMET130_PFMHT130_IDTight_v13);
       fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight_v14", &HLT_PFMET130_PFMHT130_IDTight_v14, &b_HLT_PFMET130_PFMHT130_IDTight_v14);
       fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight_v15", &HLT_PFMET130_PFMHT130_IDTight_v15, &b_HLT_PFMET130_PFMHT130_IDTight_v15);
       fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight_v16", &HLT_PFMET130_PFMHT130_IDTight_v16, &b_HLT_PFMET130_PFMHT130_IDTight_v16);
       fChain->SetBranchAddress("HLT_PFMET130_PFMHT130_IDTight_v9", &HLT_PFMET130_PFMHT130_IDTight_v9, &b_HLT_PFMET130_PFMHT130_IDTight_v9);
       fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_v11", &HLT_PFMET140_PFMHT140_IDTight_v11, &b_HLT_PFMET140_PFMHT140_IDTight_v11);
       fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_v13", &HLT_PFMET140_PFMHT140_IDTight_v13, &b_HLT_PFMET140_PFMHT140_IDTight_v13);
       fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_v14", &HLT_PFMET140_PFMHT140_IDTight_v14, &b_HLT_PFMET140_PFMHT140_IDTight_v14);
       fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_v15", &HLT_PFMET140_PFMHT140_IDTight_v15, &b_HLT_PFMET140_PFMHT140_IDTight_v15);
       fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_v16", &HLT_PFMET140_PFMHT140_IDTight_v16, &b_HLT_PFMET140_PFMHT140_IDTight_v16);
       fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_v17", &HLT_PFMET140_PFMHT140_IDTight_v17, &b_HLT_PFMET140_PFMHT140_IDTight_v17);
       fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_v18", &HLT_PFMET140_PFMHT140_IDTight_v18, &b_HLT_PFMET140_PFMHT140_IDTight_v18);
       fChain->SetBranchAddress("HLT_PFMET140_PFMHT140_IDTight_v9", &HLT_PFMET140_PFMHT140_IDTight_v9, &b_HLT_PFMET140_PFMHT140_IDTight_v9);
       fChain->SetBranchAddress("HLT_PFMET200_HBHE_BeamHaloCleaned_v5", &HLT_PFMET200_HBHE_BeamHaloCleaned_v5, &b_HLT_PFMET200_HBHE_BeamHaloCleaned_v5);
       fChain->SetBranchAddress("HLT_PFMET200_HBHE_BeamHaloCleaned_v6", &HLT_PFMET200_HBHE_BeamHaloCleaned_v6, &b_HLT_PFMET200_HBHE_BeamHaloCleaned_v6);
       fChain->SetBranchAddress("HLT_PFMET200_HBHE_BeamHaloCleaned_v7", &HLT_PFMET200_HBHE_BeamHaloCleaned_v7, &b_HLT_PFMET200_HBHE_BeamHaloCleaned_v7);
       fChain->SetBranchAddress("HLT_PFMET250_HBHECleaned_v2", &HLT_PFMET250_HBHECleaned_v2, &b_HLT_PFMET250_HBHECleaned_v2);
       fChain->SetBranchAddress("HLT_PFMET250_HBHECleaned_v3", &HLT_PFMET250_HBHECleaned_v3, &b_HLT_PFMET250_HBHECleaned_v3);
       fChain->SetBranchAddress("HLT_PFMET250_HBHECleaned_v4", &HLT_PFMET250_HBHECleaned_v4, &b_HLT_PFMET250_HBHECleaned_v4);
       fChain->SetBranchAddress("HLT_PFMET250_HBHECleaned_v5", &HLT_PFMET250_HBHECleaned_v5, &b_HLT_PFMET250_HBHECleaned_v5);
       fChain->SetBranchAddress("HLT_PFMET250_HBHECleaned_v6", &HLT_PFMET250_HBHECleaned_v6, &b_HLT_PFMET250_HBHECleaned_v6);
       fChain->SetBranchAddress("HLT_PFMET250_HBHECleaned_v7", &HLT_PFMET250_HBHECleaned_v7, &b_HLT_PFMET250_HBHECleaned_v7);
       fChain->SetBranchAddress("HLT_PFMET300_HBHECleaned_v2", &HLT_PFMET300_HBHECleaned_v2, &b_HLT_PFMET300_HBHECleaned_v2);
       fChain->SetBranchAddress("HLT_PFMET300_HBHECleaned_v3", &HLT_PFMET300_HBHECleaned_v3, &b_HLT_PFMET300_HBHECleaned_v3);
       fChain->SetBranchAddress("HLT_PFMET300_HBHECleaned_v4", &HLT_PFMET300_HBHECleaned_v4, &b_HLT_PFMET300_HBHECleaned_v4);
       fChain->SetBranchAddress("HLT_PFMET300_HBHECleaned_v5", &HLT_PFMET300_HBHECleaned_v5, &b_HLT_PFMET300_HBHECleaned_v5);
       fChain->SetBranchAddress("HLT_PFMET300_HBHECleaned_v6", &HLT_PFMET300_HBHECleaned_v6, &b_HLT_PFMET300_HBHECleaned_v6);
       fChain->SetBranchAddress("HLT_PFMET300_HBHECleaned_v7", &HLT_PFMET300_HBHECleaned_v7, &b_HLT_PFMET300_HBHECleaned_v7);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v1", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v1, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v1);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v2", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v2, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v2);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_L1ETMnoHF_v10", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_L1ETMnoHF_v10, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_L1ETMnoHF_v10);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v2", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v2, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v2);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v3", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v3, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v3);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v4", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v4, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v4);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v5", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v5, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v5);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v6", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v6, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v6);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v7", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v7, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v7);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v11", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v11, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v11);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v13", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v13, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v13);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v14", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v14, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v14);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v15", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v15, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v15);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v16", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v16, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v16);
       fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v9", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v9, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v9);
       fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v10", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v10, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v10);
       fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v12", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v12, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v12);
       fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v13", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v13, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v13);
       fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v14", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v14, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v14);
       fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v15", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v15, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v15);
       fChain->SetBranchAddress("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v9", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v9, &b_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v9);
       fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v10", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v10, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v10);
       fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v12", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v12, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v12);
       fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v13", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v13, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v13);
       fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v14", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v14, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v14);
       fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v15", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v15, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v15);
       fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v16", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v16, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v16);
       fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v17", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v17, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v17);
       fChain->SetBranchAddress("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v9", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v9, &b_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v9);
       fChain->SetBranchAddress("HLT_PFHT1050_v11", &HLT_PFHT1050_v11, &b_HLT_PFHT1050_v11);
       fChain->SetBranchAddress("HLT_PFHT1050_v12", &HLT_PFHT1050_v12, &b_HLT_PFHT1050_v12);
       fChain->SetBranchAddress("HLT_PFHT1050_v13", &HLT_PFHT1050_v13, &b_HLT_PFHT1050_v13);
       fChain->SetBranchAddress("HLT_PFHT1050_v14", &HLT_PFHT1050_v14, &b_HLT_PFHT1050_v14);
       fChain->SetBranchAddress("HLT_PFHT1050_v15", &HLT_PFHT1050_v15, &b_HLT_PFHT1050_v15);
       fChain->SetBranchAddress("HLT_PFHT1050_v16", &HLT_PFHT1050_v16, &b_HLT_PFHT1050_v16);
       fChain->SetBranchAddress("HLT_PFHT1050_v7", &HLT_PFHT1050_v7, &b_HLT_PFHT1050_v7);
       fChain->SetBranchAddress("HLT_PFHT1050_v9", &HLT_PFHT1050_v9, &b_HLT_PFHT1050_v9);
       fChain->SetBranchAddress("HLT_PFHT180_v7", &HLT_PFHT180_v7, &b_HLT_PFHT180_v7);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v1", &HLT_PFHT500_PFMET100_PFMHT100_IDTight_v1, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v1);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v10", &HLT_PFHT500_PFMET100_PFMHT100_IDTight_v10, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v10);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v3", &HLT_PFHT500_PFMET100_PFMHT100_IDTight_v3, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v3);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v5", &HLT_PFHT500_PFMET100_PFMHT100_IDTight_v5, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v5);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v6", &HLT_PFHT500_PFMET100_PFMHT100_IDTight_v6, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v6);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v7", &HLT_PFHT500_PFMET100_PFMHT100_IDTight_v7, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v7);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v8", &HLT_PFHT500_PFMET100_PFMHT100_IDTight_v8, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v8);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v9", &HLT_PFHT500_PFMET100_PFMHT100_IDTight_v9, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v9);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight_v1", &HLT_PFHT500_PFMET110_PFMHT110_IDTight_v1, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v1);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight_v10", &HLT_PFHT500_PFMET110_PFMHT110_IDTight_v10, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v10);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight_v3", &HLT_PFHT500_PFMET110_PFMHT110_IDTight_v3, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v3);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight_v5", &HLT_PFHT500_PFMET110_PFMHT110_IDTight_v5, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v5);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight_v6", &HLT_PFHT500_PFMET110_PFMHT110_IDTight_v6, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v6);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight_v7", &HLT_PFHT500_PFMET110_PFMHT110_IDTight_v7, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v7);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight_v8", &HLT_PFHT500_PFMET110_PFMHT110_IDTight_v8, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v8);
       fChain->SetBranchAddress("HLT_PFHT500_PFMET110_PFMHT110_IDTight_v9", &HLT_PFHT500_PFMET110_PFMHT110_IDTight_v9, &b_HLT_PFHT500_PFMET110_PFMHT110_IDTight_v9);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v1", &HLT_PFHT700_PFMET85_PFMHT85_IDTight_v1, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v1);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v10", &HLT_PFHT700_PFMET85_PFMHT85_IDTight_v10, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v10);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v3", &HLT_PFHT700_PFMET85_PFMHT85_IDTight_v3, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v3);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v5", &HLT_PFHT700_PFMET85_PFMHT85_IDTight_v5, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v5);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v6", &HLT_PFHT700_PFMET85_PFMHT85_IDTight_v6, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v6);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v7", &HLT_PFHT700_PFMET85_PFMHT85_IDTight_v7, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v7);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v8", &HLT_PFHT700_PFMET85_PFMHT85_IDTight_v8, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v8);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v9", &HLT_PFHT700_PFMET85_PFMHT85_IDTight_v9, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v9);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight_v1", &HLT_PFHT700_PFMET95_PFMHT95_IDTight_v1, &b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v1);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight_v10", &HLT_PFHT700_PFMET95_PFMHT95_IDTight_v10, &b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v10);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight_v3", &HLT_PFHT700_PFMET95_PFMHT95_IDTight_v3, &b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v3);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight_v5", &HLT_PFHT700_PFMET95_PFMHT95_IDTight_v5, &b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v5);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight_v6", &HLT_PFHT700_PFMET95_PFMHT95_IDTight_v6, &b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v6);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight_v7", &HLT_PFHT700_PFMET95_PFMHT95_IDTight_v7, &b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v7);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight_v8", &HLT_PFHT700_PFMET95_PFMHT95_IDTight_v8, &b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v8);
       fChain->SetBranchAddress("HLT_PFHT700_PFMET95_PFMHT95_IDTight_v9", &HLT_PFHT700_PFMET95_PFMHT95_IDTight_v9, &b_HLT_PFHT700_PFMET95_PFMHT95_IDTight_v9);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v1", &HLT_PFHT800_PFMET75_PFMHT75_IDTight_v1, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v1);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v10", &HLT_PFHT800_PFMET75_PFMHT75_IDTight_v10, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v10);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v3", &HLT_PFHT800_PFMET75_PFMHT75_IDTight_v3, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v3);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v5", &HLT_PFHT800_PFMET75_PFMHT75_IDTight_v5, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v5);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v6", &HLT_PFHT800_PFMET75_PFMHT75_IDTight_v6, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v6);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v7", &HLT_PFHT800_PFMET75_PFMHT75_IDTight_v7, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v7);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v8", &HLT_PFHT800_PFMET75_PFMHT75_IDTight_v8, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v8);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v9", &HLT_PFHT800_PFMET75_PFMHT75_IDTight_v9, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v9);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight_v1", &HLT_PFHT800_PFMET85_PFMHT85_IDTight_v1, &b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v1);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight_v10", &HLT_PFHT800_PFMET85_PFMHT85_IDTight_v10, &b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v10);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight_v3", &HLT_PFHT800_PFMET85_PFMHT85_IDTight_v3, &b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v3);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight_v5", &HLT_PFHT800_PFMET85_PFMHT85_IDTight_v5, &b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v5);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight_v6", &HLT_PFHT800_PFMET85_PFMHT85_IDTight_v6, &b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v6);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight_v7", &HLT_PFHT800_PFMET85_PFMHT85_IDTight_v7, &b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v7);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight_v8", &HLT_PFHT800_PFMET85_PFMHT85_IDTight_v8, &b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v8);
       fChain->SetBranchAddress("HLT_PFHT800_PFMET85_PFMHT85_IDTight_v9", &HLT_PFHT800_PFMET85_PFMHT85_IDTight_v9, &b_HLT_PFHT800_PFMET85_PFMHT85_IDTight_v9);
   }
   //MET filter branches
   if (is2016)
   {
       fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
       fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
       fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
       fChain->SetBranchAddress("Flag_badMuons", &Flag_badMuons, &b_Flag_badMuons);
       fChain->SetBranchAddress("Flag_duplicateMuons", &Flag_duplicateMuons, &b_Flag_duplicateMuons);
       fChain->SetBranchAddress("Flag_noBadMuons", &Flag_noBadMuons, &b_Flag_noBadMuons);
   }
   else
   {
       fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
       fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
       fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   }
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   if (isMC)
   {
       fChain->SetBranchAddress("nGenPar", &nGenPar, &b_nGenPar);
       fChain->SetBranchAddress("genParEta", genParEta, &b_genParEta);
       fChain->SetBranchAddress("genParPhi", genParPhi, &b_genParPhi);
       fChain->SetBranchAddress("genParE", genParE, &b_genParE);
       fChain->SetBranchAddress("genParPt", genParPt, &b_genParPt);
       fChain->SetBranchAddress("genParId", genParId, &b_genParId);
       fChain->SetBranchAddress("genParVx", genParVx, &b_genParVx);
       fChain->SetBranchAddress("genParVy", genParVy, &b_genParVy);
       fChain->SetBranchAddress("genParVz", genParVz, &b_genParVz);
       fChain->SetBranchAddress("genParMotherId", genParMotherId, &b_genParMotherId);
       fChain->SetBranchAddress("genParMotherIndex", genParMotherIndex, &b_genParMotherIndex);
       fChain->SetBranchAddress("genParNumDaughters", genParNumDaughters, &b_genParNumDaughters);
       fChain->SetBranchAddress("genParDaughterId1", genParDaughterId1, &b_genParDaughterId1);
       fChain->SetBranchAddress("genParDaughterId2", genParDaughterId2, &b_genParDaughterId2);
       fChain->SetBranchAddress("genParStatus", genParStatus, &b_genParStatus);
       fChain->SetBranchAddress("genParCharge", genParCharge, &b_genParCharge);
   }
   fChain->SetBranchAddress("eventRun", &eventRun, &b_eventRun);
   fChain->SetBranchAddress("eventNum", &eventNum, &b_eventNum);
   fChain->SetBranchAddress("eventLumiblock", &eventLumiblock, &b_eventLumiblock);
   fChain->SetBranchAddress("numVert", &numVert, &b_numVert);
}
// clang-format on

inline AnalysisEvent::~AnalysisEvent()
{
    if (!fChain)
    {
        return;
    }
    delete fChain->GetCurrentFile();
}

inline Int_t AnalysisEvent::GetEntry(const Long64_t entry)
{
    // Read contents of entry.
    if (!fChain)
    {
        return 0;
    }
    return fChain->GetEntry(entry);
}

inline Long64_t AnalysisEvent::LoadTree(const Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain)
    {
        return -5;
    }
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0)
    {
        return centry;
    }
    if (fChain->GetTreeNumber() != fCurrent)
    {
        fCurrent = fChain->GetTreeNumber();
    }
    return centry;
}

inline void AnalysisEvent::Show(const Long64_t entry) const
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain)
    {
        return;
    }
    fChain->Show(entry);
}

inline void AnalysisEvent::Loop()
{
    if (fChain == nullptr)
    {
        return;
    }

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0;
    Long64_t nb = 0;
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
        {
            break;
        }
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
    }
}
/*
inline bool AnalysisEvent::eTrig() const
{
    return is2016_ ? HLT_Ele32_eta2p1_WPTight_Gsf_v2 > 0
                         || HLT_Ele32_eta2p1_WPTight_Gsf_v3 > 0
                         || HLT_Ele32_eta2p1_WPTight_Gsf_v4 > 0
                         || HLT_Ele32_eta2p1_WPTight_Gsf_v5 > 0
                         || HLT_Ele32_eta2p1_WPTight_Gsf_v6 > 0
                         || HLT_Ele32_eta2p1_WPTight_Gsf_v7 > 0
                         || HLT_Ele32_eta2p1_WPTight_Gsf_v8 > 0
                   : HLT_Ele32_WPTight_Gsf_L1DoubleEG_v1 > 0
                         || HLT_Ele32_WPTight_Gsf_L1DoubleEG_v2 > 0
                         || HLT_Ele32_WPTight_Gsf_L1DoubleEG_v3 > 0
                         || HLT_Ele32_WPTight_Gsf_L1DoubleEG_v4 > 0
                         || HLT_Ele32_WPTight_Gsf_L1DoubleEG_v5 > 0
                         || HLT_Ele32_WPTight_Gsf_L1DoubleEG_v6 > 0
                         || HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7 > 0
                         || HLT_Ele35_WPTight_Gsf_v1 > 0
                         || HLT_Ele35_WPTight_Gsf_v2 > 0
                         || HLT_Ele35_WPTight_Gsf_v3 > 0
                         || HLT_Ele35_WPTight_Gsf_v4 > 0
                         || HLT_Ele35_WPTight_Gsf_v5 > 0
                         || HLT_Ele35_WPTight_Gsf_v6 > 0
                         || HLT_Ele35_WPTight_Gsf_v7 > 0;
}
*/
inline bool AnalysisEvent::muTrig() const
{
    return is2016_
               ? HLT_IsoMu24_v1 > 0 || HLT_IsoMu24_v2 > 0 || HLT_IsoMu24_v3 > 0
                     || HLT_IsoMu24_v4 > 0 || HLT_IsoTkMu24_v1 > 0
                     || HLT_IsoTkMu24_v2 > 0 || HLT_IsoTkMu24_v3 > 0
                     || HLT_IsoTkMu24_v4 > 0
               : HLT_IsoMu27_v8 > 0 || HLT_IsoMu27_v9 > 0 || HLT_IsoMu27_v10 > 0
                     || HLT_IsoMu27_v11 > 0 || HLT_IsoMu27_v12 > 0
                     || HLT_IsoMu27_v13 > 0 || HLT_IsoMu27_v14 > 0;
}
/*
inline bool AnalysisEvent::eeTrig() const
{
    return is2016_ ? HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9 > 0
                   : HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v10 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v11 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v12 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v13 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v14 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v15 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v16 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v17 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v14 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v15 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v16 > 0
                         || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17 > 0;
}
*/
/*
inline bool AnalysisEvent::muEGTrig() const
{
    // clang-format off
    return
      is2016_ ? eventRun < 280919 // different triggers for run H
                    ? HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0
                          || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4 > 0
                          || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5 > 0
                          || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6 > 0
                          || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7 > 0
                          || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8 > 0
                          || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9 > 0
                          || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3 > 0
                          || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4 > 0
                          || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5 > 0
                          || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6 > 0
                          || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7 > 0
                          || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8 > 0
                          || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9 > 0
                    : HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0
                          || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0
                          || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0
                          || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0
                          || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1 > 0
                          || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2 > 0
                          || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3 > 0
                          || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0
              : HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v11 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v13 > 0
                    || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v5 > 0
                    || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6 > 0
                    || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8 > 0
                    || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9 > 0
                    || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10 > 0
                    || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11 > 0
                    || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v12 > 0
                    || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v13 > 0
                    || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4 > 0
                    || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6 > 0
                    || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v7 > 0
                    || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v8 > 0
                    || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v9 > 0
                    || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v10 > 0
                    || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v11 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4 > 0
                    || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5 > 0;
    // clang-format on
}
*/
inline bool AnalysisEvent::mumuTrig() const
{
    return is2016_ ? (eventRun < 280919
                      && (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2 > 0
                          || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3 > 0
                          || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4 > 0
                          || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6 > 0
                          || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2 > 0
                          || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3 > 0
                          || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5 > 0))
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7 > 0
                         || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2 > 0
                         || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3 > 0
                         || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6 > 0
                   : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v8 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v9 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v10 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v11 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v12 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v13 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v14 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v1 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v2 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v3 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v4 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v7 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v8 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v1 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v2 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v3 > 0
                         || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v4 > 0;
}

#endif
