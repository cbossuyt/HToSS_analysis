#include "AnalysisEvent.hpp"
#include "TChain.h"
#include "TFile.h"
#include "TH1I.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TASImage.h"
#include "TLatex.h"
#include "TMVA/Timer.h"
#include "TTree.h"
#include "TString.h"
#include "config_parser.hpp"

#include <boost/filesystem.hpp>
#include <boost/functional/hash.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <boost/range/iterator_range.hpp>

#include <algorithm> 
#include <chrono> 
#include <fstream>
#include <iostream>
#include <regex>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <map>

std::string pythiaStatus (const Int_t status);
std::string pdgIdCode (const Int_t status, const bool unicode = false);
bool scalarGrandparent(const AnalysisEvent event, const Int_t k, const Int_t pdgId_);

uint debugCounter;

namespace fs = boost::filesystem;

int main(int argc, char* argv[])
{
    auto timerStart = std::chrono::high_resolution_clock::now(); 

    std::string config;
    std::vector<Dataset> datasets;
    double totalLumi;
    double usePreLumi;

    std::map<int, int> pdgIdMap;
    std::map<int, int> pdgIdMapStatus1;
    std::map<int, int> pdgIdMapStatus2;
    std::map<int, int> pdgIdMapStatus23;
    std::map<int, int> pdgIdMapStatus33;
    std::map<int, int> pdgIdMapStatus6X;
    std::map<int, int> pdgIdMapStatus7X;

    std::map<int, int> pdgIdMapScalarMother;

    std::map<int, int> pdgIdMapFromScalar;
    std::map<int, int> pdgIdMapFromScalarStatus1;
    std::map<int, int> pdgIdMapFromScalarStatus23;
    std::map<int, int> pdgIdMapFromScalarStatus33;
    std::map<int, int> pdgIdMapFromScalarStatus71;

    std::map<int, int> pdgIdMapScalarDecayProducts;

    std::vector <int>  nOutgoingStatus;
    std::vector <int>  nOutgoingStatus23;
    std::vector <int>  nOutgoingStatus33;

    std::map<int, int> pdgIdMapJetsFromScalar;
    std::map<int, int> pdgIdMapGenJetsFromScalar;

    std::string outFileString{"plots/distributions/output.root"};
    bool is2016_;
    int numFiles;
    Long64_t nEvents;
    Long64_t totalEvents {0};
    const std::regex mask{".*\\.root"};

////

//    int maxGenPars {0};

    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()("help,h", "Print this message.")(
        "config,c",
        po::value<std::string>(&config)->required(),
        "The configuration file to be used.")(
        "nFiles,f",
        po::value<int>(&numFiles)->default_value(-1),
        "Number of files to run over. All if set to -1.")(
        "lumi,l",
        po::value<double>(&usePreLumi)->default_value(41528.0),
        "Lumi to scale MC plots to.")(
        "outfile,o",
        po::value<std::string>(&outFileString)->default_value(outFileString),
        "Output file for plots.")(
        ",n",
        po::value<Long64_t>(&nEvents)->default_value(0),
        "The number of events to be run over. All if set to 0.")(
        "2016", po::bool_switch(&is2016_), "Use 2016 conditions (SFs, et al.).");
    po::variables_map vm;

    try
    {
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            std::cout << desc;
            return 0;
        }

        po::notify(vm);
    }
    catch (po::error& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    // Some vectors that will be filled in the parsing
    totalLumi = 0;

    try
    {
        Parser::parse_config(config,
                             datasets,
                             totalLumi);
    }
    catch (const std::exception)
    {
        std::cerr << "ERROR Problem with a confugration file, see previous "
                     "errors for more details. If this is the only error, the "
                     "problem is with the main configuration file."
                  << std::endl;
        throw;
    }

    if (totalLumi == 0.)
    {
        totalLumi = usePreLumi;
    }
    std::cout << "Using lumi: " << totalLumi << std::endl;

    bool datasetFilled{false};

    // Begin to loop over all datasets
    for (auto dataset = datasets.begin(); dataset != datasets.end(); ++dataset)
    {
        datasetFilled = false;
        TChain* datasetChain{new TChain{dataset->treeName().c_str()}};
        datasetChain->SetAutoSave(0);

        if (!datasetFilled) {
            if (!dataset->fillChain(datasetChain, numFiles)) {
                std::cerr << "There was a problem constructing the chain for " << dataset->name() << " made of " << numFiles << " files. Continuing with next dataset.\n";
                continue;
            }
            datasetFilled=true;
        }

        // extract the dataset weight. MC = (lumi*crossSection)/(totalEvents), data = 1.0
        float datasetWeight{dataset->getDatasetWeight(totalLumi)};
        std::cout << datasetChain->GetEntries() << " number of items in tree. Dataset weight: " << datasetWeight << std::endl;
        if (datasetChain->GetEntries() == 0) {
            std::cout << "No entries in tree, skipping..." << std::endl;
            continue;
        }

        AnalysisEvent event{dataset->isMC(), datasetChain, is2016_};

        Long64_t numberOfEvents{datasetChain->GetEntries()};
        if (nEvents && nEvents < numberOfEvents) numberOfEvents = nEvents;

        TMVA::Timer* lEventTimer{ new TMVA::Timer{boost::numeric_cast<int>(numberOfEvents), "Running over dataset ...", false}}; 
        lEventTimer->DrawProgressBar(0, "");
    
        totalEvents += numberOfEvents;
        for (Long64_t i{0}; i < numberOfEvents; i++) {

            lEventTimer->DrawProgressBar(i,"");

            event.GetEntry(i);      

            int nOutgoingStatusCounter {0};
            int nOutgoingStatus23Counter {0};
            int nOutgoingStatus33Counter {0};

//            std::cout << "eventNum: " << event.eventNum << std::endl;


            //////// JET STUFF

            for (Int_t k{0}; k < event.numJetPF2PAT; k++) {
                const Int_t jetPid       {event.jetPF2PATPID[k]};
                const Int_t genJetPid    {event.genJetPF2PATPID[k]};
                const Int_t fromScalar   {event.genJetPF2PATScalarGrandmother[k]};

                if ( fromScalar ) pdgIdMapJetsFromScalar[std::abs(jetPid)]++;
                if ( fromScalar ) pdgIdMapGenJetsFromScalar[std::abs(genJetPid)]++;
            }

            //////// GENERATOR PARTICLE STUFF

///

            for (Int_t k{0}; k < event.nGenPar; k++) {
                const Int_t pdgId    { std::abs(event.genParId[k]) };
		const Int_t status   { event.genParStatus[k] };
		const Int_t motherId { std::abs(event.genParMotherId[k]) };
		const Int_t numDaughters { event.genParNumDaughters[k] };
		const bool isOwnParent { pdgId == motherId ? true : false };
                const Int_t motherIndex  { std::abs(event.genParMotherIndex[k]) };
                const bool isScalarGrandparent{ scalarGrandparent(event, k, 9000006) };

//                if ( motherId == 9000006 ) std::cout << "MOTHER IS SCALAR and has " << daughters << " daughters and status " << status << " and pdgId " << pdgId << std::endl;                
	  
		if ( numDaughters == 0 && (status == 1 || status == 2 || status == 71 || status == 72) ) pdgIdMap[pdgId]++;
		if (status == 1 && numDaughters == 0) pdgIdMapStatus1[pdgId]++;
		if (status == 2 && numDaughters == 0) pdgIdMapStatus2[pdgId]++;
                if (status == 23) {
                   uint localCounter {0};
                   for (Int_t m{0}; m < event.nGenPar; m++) {
                       if (event.genParStatus[m] == 23) localCounter++;
                   }
//                   if (localCounter == 4) pdgIdMapStatus23[pdgId]++;
                    pdgIdMapStatus23[pdgId]++;
                    nOutgoingStatusCounter++;
                    nOutgoingStatus23Counter++;
                }
                if (status == 33) {
                    pdgIdMapStatus33[pdgId]++;
                    nOutgoingStatusCounter++;
                    nOutgoingStatus33Counter++;
                }
		if ((status == 61 && status == 62 || status == 63) && numDaughters == 0) pdgIdMapStatus6X[pdgId]++;
		if ((status == 71 || status == 72 || status == 74) && numDaughters == 0) pdgIdMapStatus7X[pdgId]++;

                // particles directly outgoing from scalar decays
                if (motherId == 9000006 && ( status == 23 || status == 33) && !isOwnParent) pdgIdMapScalarMother[pdgId]++;

                // search ancestry of final state for final state particles from scalars
                if ( numDaughters == 0 ) {
                    if ( isScalarGrandparent == true ) {
                        pdgIdMapFromScalar[pdgId]++;
                        if (status == 1) pdgIdMapFromScalarStatus1[pdgId]++;
                        if (status == 23) pdgIdMapFromScalarStatus23[pdgId]++;
                        if (status == 33) pdgIdMapFromScalarStatus33[pdgId]++;
                        if (status == 71) pdgIdMapFromScalarStatus71[pdgId]++;
                    }
                }

                // search for physical particles from scalar decays before they decay
                if ( isScalarGrandparent == true ) {
                    const bool quarkParent = std::abs(motherId < 7);
                    const bool isQuark = std::abs(pdgId < 7);

                    if ( status == 1 && motherId == 9000006 ) pdgIdMapScalarDecayProducts[pdgId]++; // final state is quarks directly from scalar

                    // grouped colour singlet quarks
                    else if ( status == 23 && motherId == 9000006 && numDaughters == 0 ) pdgIdMapScalarDecayProducts[pdgId]++; // quarks that don't decay/get grouped into colour singlet
                    else if ( status == 71 && isOwnParent && numDaughters == 0 ) pdgIdMapScalarDecayProducts[pdgId]++; // quarks from status 23 quarks that don't decay further

                    // omega, eta, eta', omega, kshort
                    else if ( status == 2 && quarkParent && !isQuark ) pdgIdMapScalarDecayProducts[pdgId]++; // is from quark from scalar (i.e. colour grouping into physical particle)
                    else if ( status == 71 && quarkParent && !isQuark ) pdgIdMapScalarDecayProducts[pdgId]++; // is from quark from scalar (i.e. colour grouping into physical particle)

                    // k-short, pi+
                    else if ( status == 1 && quarkParent && !isQuark ) pdgIdMapScalarDecayProducts[pdgId]++; // is from quark from scalar and is final status
                 }

//		if ( !daughters ) {
//              if ( isScalarGrandparent == true ) {
//		std::cout << "index / pdgId / mother / motherIndex / nDaughers / status: " << std::endl;
//  	        std::cout << k << " / " << pdgIdCode( pdgId, false ) << " / " << pdgIdCode( motherId, false ) << " / " << motherIndex << " / " << numDaughters << " / " << pythiaStatus( status ) << std::endl;
//              }
	    } ///
            nOutgoingStatus.emplace_back(nOutgoingStatusCounter);
            nOutgoingStatus23.emplace_back(nOutgoingStatus23Counter);
            nOutgoingStatus33.emplace_back(nOutgoingStatus33Counter);
//	    std::cout << std::endl;
	} 
    }

    std::cout << std::endl;
    std::cout << "Total no. of events:\t\t\t" << totalEvents << std::endl;
    std::cout << std::endl;

    // Do scalable histograms
    int nPdgIds                   = pdgIdMap.size();                   // number of different final state pdgIds
    int nPdgIdsStatus1            = pdgIdMapStatus1.size();            // number of different final state pdgIds with status 1
    int nPdgIdsStatus2            = pdgIdMapStatus2.size();            // number of different final state pdgIds with status 2
    int nPdgIdsStatus23           = pdgIdMapStatus23.size();           // number of different final state pdgIds with status 23
    int nPdgIdsStatus33           = pdgIdMapStatus33.size();           // number of different final state pdgIds with status 33
    int nPdgIdsStatus6X           = pdgIdMapStatus6X.size();           // number of different final state pdgIds with status 6X
    int nPdgIdsStatus7X           = pdgIdMapStatus7X.size();           // number of different final state pdgIds with status 7X
    int nPdgIdsScalarMother       = pdgIdMapScalarMother.size();       // number of different particles produced from scalar decay
    int nPdgIdsFromScalar         = pdgIdMapFromScalar.size();         // number of different final state pdgIds that are descended from scalar particle
    int nPdgIdsFromScalarStatus1  = pdgIdMapFromScalarStatus1.size();  // number of different final state pdgIds that are descended from scalar particle
    int nPdgIdsFromScalarStatus23 = pdgIdMapFromScalarStatus23.size(); // number of different final state pdgIds that are descended from scalar particle
    int nPdgIdsFromScalarStatus33 = pdgIdMapFromScalarStatus33.size(); // number of different final state pdgIds that are descended from scalar particle
    int nPdgIdsFromScalarStatus71 = pdgIdMapFromScalarStatus71.size(); // number of different final state pdgIds that are descended from scalar particle

    int nPdgIdsScalarDecayProducts = pdgIdMapScalarDecayProducts.size(); // number of physical decay products from scalars (before decaying)

    int nOutgoingStatusMax   = *std::max_element(nOutgoingStatus.begin(),nOutgoingStatus.end());
    int nOutgoingStatus23Max = *std::max_element(nOutgoingStatus23.begin(),nOutgoingStatus23.end());
    int nOutgoingStatus33Max = *std::max_element(nOutgoingStatus33.begin(),nOutgoingStatus33.end());

    int nPdgIdsFromJetsFromScalar     = pdgIdMapJetsFromScalar.size();
    int nPdgIdsFromGenJetsFromScalar  = pdgIdMapGenJetsFromScalar.size();

    // status == 1 for final state particles
    // status == 2 for a decayed Standard Model hadron or tau or mu lepton, excepting virtual intermediate states thereof (i.e. the particle must undergo a normal decay, not e.g. a shower branching);
    // status == 61-63 for particles produced by beam-remnant treatment
    // status == 71 for partons in preparation of hadronization process and 72+74 (but exclude particles who are their own parent)

    TH1I* h_pdgId            {new TH1I{"h_pdgId",            "Final state content - all final state codes"   , nPdgIds,           0, Double_t(nPdgIds)           }};
    TH1I* h_pdgIdStatus1     {new TH1I{"h_pdgIdStatus1",     "Final state content - status code 1"           , nPdgIdsStatus1,    0, Double_t(nPdgIdsStatus1)    }};
    TH1I* h_pdgIdStatus2     {new TH1I{"h_pdgIdStatus2",     "Final state content - status code 2"           , nPdgIdsStatus2,    0, Double_t(nPdgIdsStatus2)    }};
    TH1I* h_pdgIdStatus23    {new TH1I{"h_pdgIdStatus23",    "Outgoing particles of hardest subprocess"      , nPdgIdsStatus23,   0, Double_t(nPdgIdsStatus23)   }};
    TH1I* h_pdgIdStatus33    {new TH1I{"h_pdgIdStatus33",    "Outgoing particles of subsequent subprocesses" , nPdgIdsStatus33,   0, Double_t(nPdgIdsStatus33)   }};
    TH1I* h_pdgIdStatus6X    {new TH1I{"h_pdgIdStatus6X",    "Final state content - status code 6X"          , nPdgIdsStatus6X,   0, Double_t(nPdgIdsStatus6X)   }};
    TH1I* h_pdgIdStatus7X    {new TH1I{"h_pdgIdStatus7X",    "Final state content - status code 7X"          , nPdgIdsStatus7X,   0, Double_t(nPdgIdsStatus7X)   }};

    TH1I* h_pdgIdScalarMother       {new TH1I{"h_pdgIdScalarMother",         "Direct decays from scalars"                      , nPdgIdsScalarMother, 0, Double_t(nPdgIdsScalarMother) }};

    TH1I* h_pdgIdFromScalar          {new TH1I{"h_pdgIdFromScalar",          "Final state content from scalars"                , nPdgIdsFromScalar, 0, Double_t(nPdgIdsFromScalar) }}; 
    TH1I* h_pdgIdFromScalarStatus1   {new TH1I{"h_pdgIdFromScalarStatus1",   "Final state content from scalars with status 1"  , nPdgIdsFromScalarStatus1, 0, Double_t(nPdgIdsFromScalarStatus1) }}; 
    TH1I* h_pdgIdFromScalarStatus23  {new TH1I{"h_pdgIdFromScalarStatus23",  "Final state content from scalars with status 23" , nPdgIdsFromScalarStatus23, 0, Double_t(nPdgIdsFromScalarStatus23) }}; 
    TH1I* h_pdgIdFromScalarStatus33  {new TH1I{"h_pdgIdFromScalarStatus33",  "Final state content from scalars with status 33" , nPdgIdsFromScalarStatus33, 0, Double_t(nPdgIdsFromScalarStatus33) }}; 
    TH1I* h_pdgIdFromScalarStatus71  {new TH1I{"h_pdgIdFromScalarStatus71",  "Final state content from scalars with status 71" , nPdgIdsFromScalarStatus71, 0, Double_t(nPdgIdsFromScalarStatus71) }}; 

    TH1* h_pdgIdMapScalarDecayProducts {new TH1I{"h_pdgIdMapScalarDecayProducts", "Physical decay products from scalar decays" , nPdgIdsScalarDecayProducts, 0, Double_t(nPdgIdsScalarDecayProducts) }};

    TH1I* h_outgoingStatus   {new TH1I{"h_outgoingStatus"  , "Number of outgoing particles - hardest subprocess"   , nOutgoingStatusMax+1   , -0.5, Double_t(nOutgoingStatusMax+0.5) }};
    TH1I* h_outgoingStatus23 {new TH1I{"h_outgoingStatus23", "Number of outgoing particles - hardest subprocess"   , nOutgoingStatus23Max+1 , -0.5, Double_t(nOutgoingStatus23Max+0.5) }};
    TH1I* h_outgoingStatus33 {new TH1I{"h_outgoingStatus33", "Number of outgoing particles - hardest subprocess"   , nOutgoingStatus33Max+1 , -0.5, Double_t(nOutgoingStatus33Max+0.5) }};

    TH1I* h_pdgIdMapJetsFromScalar {new TH1I{"h_pdgIdMapJetsFromScalar",       "reco jet which is descended from generator scalar", nPdgIdsFromJetsFromScalar+1,    -0.5, Double_t(nPdgIdsFromJetsFromScalar+0.5)    }};
    TH1I* h_pdgIdMapGenJetsFromScalar {new TH1I{"h_pdgIdMapGenJetsFromScalar", "gen jet which is descended from generator scalar",  nPdgIdsFromGenJetsFromScalar+1, -0.5, Double_t(nPdgIdsFromGenJetsFromScalar+0.5) }};

    uint binCounter {1};
    for (auto it = pdgIdMap.begin(); it != pdgIdMap.end(); ++it) {
//        std::cout << "Add " << it->second << " to bin " << binCounter << " for pdgId " << it->first << std::endl;
        h_pdgId->SetBinContent(binCounter, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
//        const char *label2 = ( pdgIdCode(it->first, true) ).c_str();
//        std::cout << "label : " << label2 << std::endl;
        h_pdgId->GetXaxis()->SetBinLabel(binCounter, label);
        binCounter++;
    }

    uint binCounterStatus1 {1};
    for (auto it = pdgIdMapStatus1.begin(); it != pdgIdMapStatus1.end(); ++it) {
        h_pdgIdStatus1->SetBinContent(binCounterStatus1, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdStatus1->GetXaxis()->SetBinLabel(binCounterStatus1, label);
        binCounterStatus1++;
    }
    uint binCounterStatus2 {1};
    for (auto it = pdgIdMapStatus2.begin(); it != pdgIdMapStatus2.end(); ++it) {
        h_pdgIdStatus2->SetBinContent(binCounterStatus2, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdStatus2->GetXaxis()->SetBinLabel(binCounterStatus2, label);
        binCounterStatus2++;
    }
    uint binCounterStatus23 {1};
    for (auto it = pdgIdMapStatus23.begin(); it != pdgIdMapStatus23.end(); ++it) {
        h_pdgIdStatus23->SetBinContent(binCounterStatus23, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdStatus23->GetXaxis()->SetBinLabel(binCounterStatus23, label);
        binCounterStatus23++;
    }
    uint binCounterStatus33 {1};
    for (auto it = pdgIdMapStatus33.begin(); it != pdgIdMapStatus33.end(); ++it) {
        h_pdgIdStatus33->SetBinContent(binCounterStatus33, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdStatus33->GetXaxis()->SetBinLabel(binCounterStatus33, label);
        binCounterStatus33++;
    }
    uint binCounterStatus6X {1};
    for (auto it = pdgIdMapStatus6X.begin(); it != pdgIdMapStatus6X.end(); ++it) {
        h_pdgIdStatus6X->SetBinContent(binCounterStatus6X, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdStatus6X->GetXaxis()->SetBinLabel(binCounterStatus6X, label);
        binCounterStatus6X++;
    }
    uint binCounterStatus7X {1};
    for (auto it = pdgIdMapStatus7X.begin(); it != pdgIdMapStatus7X.end(); ++it) {
        h_pdgIdStatus7X->SetBinContent(binCounterStatus7X, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdStatus7X->GetXaxis()->SetBinLabel(binCounterStatus7X, label);
        binCounterStatus7X++;
    }

    uint binCounterScalarMother {1};
    for (auto it = pdgIdMapScalarMother.begin(); it != pdgIdMapScalarMother.end(); ++it) {
        h_pdgIdScalarMother->SetBinContent(binCounterScalarMother, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdScalarMother->GetXaxis()->SetBinLabel(binCounterScalarMother, label);
        binCounterScalarMother++;
    }

    uint binCounterFromScalar {1};
    for (auto it = pdgIdMapFromScalar.begin(); it != pdgIdMapFromScalar.end(); ++it) {
        h_pdgIdFromScalar->SetBinContent(binCounterFromScalar, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdFromScalar->GetXaxis()->SetBinLabel(binCounterFromScalar, label);
        binCounterFromScalar++;
    }

    uint binCounterFromScalarStatus1 {1};
    for (auto it = pdgIdMapFromScalarStatus1.begin(); it != pdgIdMapFromScalarStatus1.end(); ++it) {
        h_pdgIdFromScalarStatus1->SetBinContent(binCounterFromScalarStatus1, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdFromScalarStatus1->GetXaxis()->SetBinLabel(binCounterFromScalarStatus1, label);
        binCounterFromScalarStatus1++;
    }

    uint binCounterFromScalarStatus23 {1};
    for (auto it = pdgIdMapFromScalarStatus23.begin(); it != pdgIdMapFromScalarStatus23.end(); ++it) {
        h_pdgIdFromScalarStatus23->SetBinContent(binCounterFromScalarStatus23, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdFromScalarStatus23->GetXaxis()->SetBinLabel(binCounterFromScalarStatus23, label);
        binCounterFromScalarStatus23++;
    }

    uint binCounterFromScalarStatus33 {1};
    for (auto it = pdgIdMapFromScalarStatus33.begin(); it != pdgIdMapFromScalarStatus33.end(); ++it) {
        h_pdgIdFromScalarStatus33->SetBinContent(binCounterFromScalarStatus33, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdFromScalarStatus33->GetXaxis()->SetBinLabel(binCounterFromScalarStatus33, label);
        binCounterFromScalarStatus33++;
    }

    uint binCounterFromScalarStatus71 {1};
    for (auto it = pdgIdMapFromScalarStatus71.begin(); it != pdgIdMapFromScalarStatus71.end(); ++it) {
        h_pdgIdFromScalarStatus71->SetBinContent(binCounterFromScalarStatus71, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdFromScalarStatus71->GetXaxis()->SetBinLabel(binCounterFromScalarStatus71, label);
        binCounterFromScalarStatus71++;
    }

    uint binCounterScalarDecayProducts {1};
    for (auto it = pdgIdMapScalarDecayProducts.begin(); it != pdgIdMapScalarDecayProducts.end(); ++it) {
        h_pdgIdMapScalarDecayProducts->SetBinContent(binCounterScalarDecayProducts, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdMapScalarDecayProducts->GetXaxis()->SetBinLabel(binCounterScalarDecayProducts, label);
        binCounterScalarDecayProducts++;        
    }

    for (auto it = nOutgoingStatus.begin(); it != nOutgoingStatus.end(); ++it) {
        h_outgoingStatus->Fill(*it);
    }
    for (auto it = nOutgoingStatus23.begin(); it != nOutgoingStatus23.end(); ++it) {
        h_outgoingStatus23->Fill(*it);
    }
    for (auto it = nOutgoingStatus33.begin(); it != nOutgoingStatus33.end(); ++it) {
        h_outgoingStatus33->Fill(*it);
    }

    uint binCounterJetsFromScalar {1};
    for (auto it = pdgIdMapJetsFromScalar.begin(); it != pdgIdMapJetsFromScalar.end(); ++it) {
        h_pdgIdMapJetsFromScalar->SetBinContent(binCounterJetsFromScalar, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdMapJetsFromScalar->GetXaxis()->SetBinLabel(binCounterJetsFromScalar, label);
        binCounterJetsFromScalar++;
    }
    uint binCounterGenJetsFromScalar {1};
    for (auto it = pdgIdMapGenJetsFromScalar.begin(); it != pdgIdMapGenJetsFromScalar.end(); ++it) {
        h_pdgIdMapGenJetsFromScalar->SetBinContent(binCounterGenJetsFromScalar, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        h_pdgIdMapGenJetsFromScalar->GetXaxis()->SetBinLabel(binCounterGenJetsFromScalar, label);
        binCounterGenJetsFromScalar++;
    }

    TFile* outFile{new TFile{outFileString.c_str(), "RECREATE"}};
    outFile->cd();

    h_pdgId->Write();
    h_pdgIdStatus1->Write();
    h_pdgIdStatus2->Write();
    h_pdgIdStatus23->Write();
    h_pdgIdStatus33->Write();
    h_pdgIdStatus6X->Write();
    h_pdgIdStatus7X->Write();
    h_pdgIdScalarMother->Write();
    h_pdgIdFromScalar->Write();
    h_pdgIdFromScalarStatus1->Write();
    h_pdgIdFromScalarStatus23->Write();
    h_pdgIdFromScalarStatus33->Write();
    h_pdgIdFromScalarStatus71->Write();
    h_pdgIdMapScalarDecayProducts->Write();
    h_outgoingStatus->Write();
    h_outgoingStatus23->Write();
    h_outgoingStatus33->Write();
    h_pdgIdMapJetsFromScalar->Write();
    h_pdgIdMapGenJetsFromScalar->Write();

    outFile->Close();

//    std::cout << "Max nGenPar: " << maxGenPars << std::endl;    
    auto timerStop = std::chrono::high_resolution_clock::now(); 
    auto duration  = std::chrono::duration_cast<std::chrono::seconds>(timerStop - timerStart);

    std::cout << "\nFinished. Took " << duration.count() << " seconds" <<std::endl;
}

std::string pythiaStatus (const Int_t status) {

   std::string message;
   message += std::to_string(status);

   switch (status) {

       case 0 : message += " : an empty entry, with no meaningful information"; break;
       case 1 : message += " : a final-state particle"; break;
       case 2 : message += " : a decayed Standard Model hadron or tau or mu lepton, excepting virtual intermediate states thereof"; break;
       case 3 : message += " : a documentation entry - not used in PYTHIA"; break;
       case 4 : message += " : an incoming beam particle"; break;

       case 11 : message += " : beam particles - the event as a whole";
           break;
       case 12 : message += " : beam particles - incoming beam";
           break;
       case 13 : message += " : beam particles - incoming beam-inside-beam (e.g. gamma inside e)";
           break;
       case 14 : message += " : beam particles - outgoing elastically scattered";
           break;
       case 15 : message += " : beam particles - outgoing diffractively scattered";
           break;
       case 21 : message += " : particles of the hardest subprocess - incoming";
           break;
       case 22 : message += " : particles of the hardest subprocess - intermediate (intended to have preserved mass)";
           break;
       case 23 : message += " : particles of the hardest subprocess - outgoing";
           break;
       case 24 : message += " : particles of the hardest subprocess - outgoing, nonperturbatively kicked out in diffraction";
           break;
       case 31 : message += " : particles of subsequent subprocesses - incoming";
           break;
       case 32 : message += " : particles of subsequent subprocesses - intermediate (intended to have preserved mass)";
           break;
       case 33 : message += " : particles of subsequent subprocesses - outgoing";
           break;
       case 34 : message += " : particles of subsequent subprocesses - incoming that has already scattered";
           break;
       case 41 : message += " : particles produced by initial-state-showers - incoming on spacelike main branch";
           break;
       case 42 : message += " : particles produced by initial-state-showers - incoming copy of recoiler";
           break;
       case 43 : message += " : particles produced by initial-state-showers - outgoing produced by a branching";
           break;
       case 44 : message += " : particles produced by initial-state-showers - outgoing shifted by a branching";
           break;
       case 45 : message += " : particles produced by initial-state-showers - incoming rescattered parton, with changed kinematics owing to ISR in the mother system (cf. status 34)";
           break;
       case 46 : message += " : particles produced by initial-state-showers - incoming copy of recoiler when this is a rescattered parton (cf. status 42)";
           break;
       case 47 : message += " : particles produced by initial-state-showers - a W or Z gauge boson produced in the shower evolution";
           break;
       case 51 : message += " : particles produced by final-state-showers - outgoing produced by parton branching";
           break;
       case 52 : message += " : particles produced by final-state-showers - outgoing copy of recoiler, with changed momentum";
           break;
       case 53 : message += " : particles produced by final-state-showers - copy of recoiler when this is incoming parton, with changed momentum";
           break;
       case 54 : message += " : particles produced by final-state-showers - copy of a recoiler, when in the initial state of a different system from the radiator";
           break;
       case 55 : message += " : particles produced by final-state-showers - copy of a recoiler, when in the final state of a different system from the radiator";
           break;
       case 56 : message += " : particles produced by final-state-showers - a W or Z gauge boson produced in a shower branching (special case of 51)";
           break;
       case 61 : message += " : particles produced by beam-remnant treatment - incoming subprocess particle with primordial kT included";
           break;
       case 62 : message += " : particles produced by beam-remnant treatment - outgoing subprocess particle with primordial kT included";
           break;
       case 63 : message += " : particles produced by beam-remnant treatment - outgoing beam remnant";
           break;
       case 71 : message += " : partons in preparation of hadronization process - copied partons to collect into contiguous colour singlet";
           break;
       case 72 : message += " : partons in preparation of hadronization process - copied recoiling singlet when ministring collapses to one hadron and momentum has to be reshuffled";
           break;
       case 73 : message += " : partons in preparation of hadronization process - combination of very nearby partons into one";
           break;
       case 74 : message += " : partons in preparation of hadronization process - combination of two junction quarks (+ nearby gluons) to a diquark";
           break;
       case 75 : message += " : partons in preparation of hadronization process - gluons split to decouple a junction-antijunction pair";
           break;
       case 76 : message += " : partons in preparation of hadronization process - partons with momentum shuffled to decouple a junction-antijunction pair";
           break;
       case 77 : message += " : partons in preparation of hadronization process - temporary opposing parton when fragmenting first two strings in to junction (should disappear again)";
           break;
       case 78 : message += " : partons in preparation of hadronization process - temporary combined diquark end when fragmenting last string in to junction (should disappear again)";
           break;
       case 81 : message += " : primary hadrons produced by hadronization process - from ministring into one hadron";
           break;
       case 82 : message += " : primary hadrons produced by hadronization process - from ministring into two hadron";
           break;
       case 83 : message += " : primary hadrons produced by hadronization process - from normal string (the difference between the two is technical, whether fragmented off from the top of the string system or from the bottom, useful for debug only)";
           break;
       case 84 : message += " : primary hadrons produced by hadronization process - from normal string (the difference between the two is technical, whether fragmented off from the top of the string system or from the bottom, useful for debug only)";
           break;
       case 85 : message += " : primary hadrons produced by hadronization process - primary produced hadrons in junction fragmentation of the first two string legs in to the junction, in order of treatment";
           break;
       case 86 : message += " : primary hadrons produced by hadronization process - primary produced hadrons in junction fragmentation of the first two string legs in to the junction, in order of treatment";
           break;
       case 91 : message += " : particles produced in decay process, or by Bose-Einstein effects - normal decay products";
           break;
       case 92 : message += " : particles produced in decay process, or by Bose-Einstein effects - decay products after oscillation";
           break;
       case 93 : message += " : particles produced in decay process, or by Bose-Einstein effects - decay handled by external program, normally or with oscillation";
           break;
       case 94 : message += " : particles produced in decay process, or by Bose-Einstein effects - decay handled by external program, normally or with oscillation";
           break;
       case 99 : message += " : particles produced in decay process, or by Bose-Einstein effects - particles with momenta shifted by Bose-Einstein effects";
           break;
       case 101 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - when a string system contains two such long-lived particles, the system is split up by the production of a new q-qbar pair";
           break;
       case 102 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - partons rearranged from the long-lived particle end to prepare for fragmentation from this end";
           break;
       case 103 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - intermediate half-R-hadron ormed when a colour octet particle (like the gluino) has been fragmented on one side, but not yet on the other";
           break;
       case 104 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - an R-hadron";
           break;
       case 105 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - partons or particles formed together with the R-hadron during the fragmentation treatment";
           break;
       case 106 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - subdivision of an R-hadron into its flavour content, with momentum split accordingly, in preparation of the decay of the heavy new particle, if it is unstable";
           break;
       case 107 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - two temporary leftover gluons joined into one in the formation of a gluino-gluon R-hadron";
           break;
//       default: message += " : unknown status code ";
   }

   return message;
}

std::string pdgIdCode (const Int_t parId, const bool unicode) {

   std::string particle;
   int id = std::abs(parId);

   switch (id) {
      case 0 : particle += "-"; break;
      case 1 : particle += "d"; break;
      case 2 : particle += "u"; break;
      case 3 : particle += "s"; break;
      case 4 : particle += "c"; break;
      case 5 : particle += "b"; break;
      case 6 : particle += "t"; break;
      case 9 : particle += "g"; break;

      case 11 : particle += "e"; break;
      case 12 : particle += unicode ? "\u03BD_e" : "#nu_{e}" ; break;
      case 13 : particle += unicode ? "\u03BC" : "#mu"; break;
      case 14 : particle += unicode ? "\u03BD_\u03BC" : "#nu_{#mu}"; break;
      case 15 : particle += unicode ? "\u03C4" : "#tau"; break;
      case 16 : particle += unicode ? "\u03BD_\u03C4" : "#nu_{#tau}"; break;

      case 21 : particle += "g"; break;
      case 22 : particle += unicode ? "Y" : "#gamma"; break;
      case 23 : particle += "Z_{0}"; break;
      case 24 : particle += "W"; break;
      case 25 : particle += "h_{0}"; break;

      case 111 : particle += unicode ? "\u03C00" : "#pi^{0}"; break;
      case 211 : particle += unicode ? "\u03C0" : "#pi^{+}"; break;
      case 113 : particle += unicode ? "\u03C1" : "#rho(770)^{0}"; break;
      case 213 : particle += unicode ? "\u03C1" : "#rho(770)^{+}"; break;
      case 225 : particle += unicode ? "f2(1270)" : "f_{2}(1270)"; break;

      case 130 : particle += unicode ? "K0_L" : "K_{L}^{0}"; break;
      case 310 : particle += unicode ? "K0_S" : "K_{S}^{0}"; break;
      case 311 : particle += unicode ? "K0" : "K^{0}"; break;
      case 321 : particle += unicode ? "K" : "K^{+}"; break;

      case 313 : particle += unicode ? "K*0" : "K^{*}(892)^{0}"; break;
      case 315 : particle += unicode ? "K*0_2" : "K^{*}_{2}(1430)"; break;
      case 323 : particle += unicode ? "K*" : "K^{*}(892)^{+}"; break;

      case 10311 : particle += unicode ? "K0*_0(1430)" : "K^{*}_{0}(1430)"; break;

      case 221 : particle += unicode ? "\u03B7" : "#eta"; break;
      case 331 : particle += unicode ? "\u03B7\'" : "#eta'"; break;
      case 223 : particle += unicode ? "\u03C9" : "#omega"; break;
      case 333 : particle += unicode ? "\u03C6" : "#phi"; break;

      case 411 : particle += unicode ? "D" : "D^{+}"; break;
      case 421 : particle += unicode ? "D0" : "D^{0}"; break;
      case 413 : particle += unicode ? "D*" : "D^{*} 2010)^{+}"; break;
      case 423 : particle += unicode ? "D*0" : "D^{*}(2007)^{0}"; break;
      case 431 : particle += unicode ? "D_S" : "D_{s}^{+}"; break;
      case 433 : particle += unicode ? "D*_S" : "D_{s}^{*+}"; break;
      case 443 : particle += unicode ? "J/\u03C8" : "J/#psi"; break;
      case 445 : particle += unicode ? "\u03C7_c2(1P)" : "#chi_{c2}(1P)"; break;
 
      case 511 : particle += unicode ? "B0" : "B^{0}"; break;
      case 521 : particle += unicode ? "B" : "B^{+}"; break;
      case 513 : particle += unicode ? "B*0" : "B^{*0}"; break;
      case 523 : particle += unicode ? "B*" : "B^{*+}"; break;
      case 531 : particle += unicode ? "B0_S" : "B^{0}_{s}"; break;
      case 533 : particle += unicode ? "B*0_S" : "B^{*0}_{s}"; break;

      case 1101 : particle += unicode ? "(dd)0" : "(dd)_{0}"; break;
      case 2101 : particle += unicode ? "(ud)0" : "(ud)_{0}"; break;
      case 2103 : particle += unicode ? "(ud)1" : "(ud)_{1}"; break;
      case 2203 : particle += unicode ? "(uu)1" : "(uu)_{1}"; break;

      case 2212 : particle += "p"; break;
      case 2112 : particle += "n"; break;
      case 2224 : particle += unicode ? "\u0394++" : "#Delta^{++}"; break;
      case 2214 : particle += unicode ? "\u0394+" : "#Delta^{+}"; break;
      case 2114 : particle += unicode ? "\u03940" : "#Delta^{0}"; break;
      case 1114 : particle += unicode ? "\u0394-" : "#Delta^{-}"; break;

      case 3122 : particle += unicode ? "\u0394" : "#Lambda"; break;
      case 3222 : particle += unicode ? "\u03A3" : "#Sigma^{+}"; break;
      case 3224 : particle += unicode ? "\u03A3*" : "#Sigma^{*+}"; break;
      case 3212 : particle += unicode ? "\u03A30" : "#Sigma^{0}"; break;
      case 3322 : particle += unicode ? "\u03A30" : "#Xi^{0}"; break;
      case 3312 : particle += unicode ? "\u03A3-" : "#Xi^{-}"; break;
      case 3112 : particle += unicode ? "\u03A3-" : "#Sigma{-}"; break;

      case 3324 : particle += unicode ? "\u039E*0" : "#Xi^{*0}"; break;
      case 3334 : particle += unicode ? "\u03A9-" : "#Omega^{-}"; break;

      case 4214 : particle += unicode ? "\u03A3*_C" : "#Sigma_{c}^{*+}"; break;
      case 4222 : particle += unicode ? "\u03A3_C" : "#Sigma_{c}^{++}"; break;
      case 4122 : particle += unicode ? "\u039BC" : "#Lambda_{c}^{+}"; break;
      case 4114 : particle += unicode ? "\u03A3*0_C" : "#Sigma_{c}^{*0}"; break;
      case 4224 : particle += unicode ? "\u03A3+C" : "#Sigma_{c}^{*++}"; break;

      case 5122 : particle += unicode ? "\u039B0_b" : "#Lambda^{0}_{b}"; break;
      case 5212 : particle += unicode ? "\u03A30_b" : "#Sigma_{b}^{0}"; break;
      case 5232 : particle += unicode ? "\u039E0_b" : "#Xi_{b}^{0}"; break;

      case 10313: particle += unicode ? "K0_1 (1270)" : "K_{1} (1270)^{0}"; break;
      case 10441: particle += unicode ? "\u03C7_c0(1P)" : "#chi_{c0}(1P)"; break;
      case 10551: particle += unicode ? "\u03C7_b0(1P)" : "#chi_{b0}(1P)"; break;
      case 20443: particle += unicode ? "\u03C7_c1(1P)" : "#chi_{c1}(1P)"; break;
      case 20313: particle += unicode ? "K0_1 (1400)" : "K_{1}(1400)^0}"; break;
      case 20213: particle += unicode ? "a1" : "a_{1} (1260)^{+}"; break;

      case 9000006 : particle += unicode ? "Scalar" : "SCALAR"; break;
      case 9010221 : particle += unicode ? "f0(980)" : "f_{0}(980)"; break;

      default : {particle += std::to_string(std::abs(parId)); /*std::cout << "UNKNOWN PID: " << parId << std::endl;*/}
   }

   if ( unicode && (parId == 211 || parId == 213 || parId == 321 || parId == 323 || parId == 411 || parId == 431 || parId == 433 || parId == 521 || parId == 4122 || parId == 20213 || parId == 4214 || parId == 523
   || parId == 3224 || parId == 3222
))  particle += "+";

   if ( unicode && parId < 0) particle += "-";
   if ( unicode && (parId == -4222 || parId == -4224) ) particle += "-";
   if ( unicode && (parId == 4222  || parId == 4222) ) particle += "+";

   return particle;
}

bool scalarGrandparent (const AnalysisEvent event, const Int_t k, const Int_t grandparentId) {

    const Int_t pdgId        { std::abs(event.genParId[k]) };
    const Int_t numDaughters { event.genParNumDaughters[k] };
    const Int_t motherId     { std::abs(event.genParMotherId[k]) };
    const Int_t motherIndex  { std::abs(event.genParMotherIndex[k]) };


    if (motherId == 0 || motherIndex == -1) return false; // if no parent, then mother Id is null and there's no index, quit search
    else if (motherId == std::abs(grandparentId)) return true; // if mother is granparent being searched for, return true
    else {
//        std::cout << "Going up the ladder ... pdgId = " << pdgId << " : motherIndex = " << motherIndex << " : motherId = " << motherId << std::endl;
//        debugCounter++;
//        std::cout << "debugCounter: " << debugCounter << std::endl;
        return scalarGrandparent(event, motherIndex, grandparentId); // otherwise check mother's mother ...
    }
}
