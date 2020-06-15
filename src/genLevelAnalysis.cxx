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

namespace fs = boost::filesystem;

int main(int argc, char* argv[])
{
    std::string config;
    std::vector<Dataset> datasets;
    double totalLumi;
    double usePreLumi;

    std::map<int, int> pdgIdMap;

    std::string outFileString{"plots/distributions/output.root"};
    bool is2016_;
    int numFiles;
    Long64_t nEvents;
    Long64_t totalEvents {0};
    const std::regex mask{".*\\.root"};

// status == 1 for final state particles
// status == 2 for a decayed Standard Model hadron or tau or mu lepton, excepting virtual intermediate states thereof (i.e. the particle must undergo a normal decay, not e.g. a shower branching);
// status == 61-63 for particles produced by beam-remnant treatment
// status == 71 for partons in preparation of hadronization process and 72+74 (but exclude particles who are their own parent)

    TH1I* histPdgId{new TH1I{"histPdgId", "Final state content", 5001, -.5, 5000.5}};
    TH1I* histPdgId_2{new TH1I{"histPdgId_2", "Final state content", 80, -.5, 79.5}};
    TH1I* histPdgIdStatus1{new TH1I{"histPdgIdStatus1", "Final state content", 501, -.5, 500.5}};
    TH1I* histPdgIdStatus2{new TH1I{"histPdgIdStatus2", "Decayed SM hadron or tau or mu", 501, -0.5, 500.5}};
    TH1I* histPdgIdStatus6X{new TH1I{"histPdgIdStatus6X", "Beam remnants", 501, -0.5, 500.5}};
    TH1I* histPdgIdStatus7X{new TH1I{"histPdgIdStatus7X", "Partons in preparation of hadronization process", 501, -0.5, 500.5}};    

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
        if (	nEvents && nEvents < numberOfEvents) numberOfEvents = nEvents;

        TMVA::Timer* lEventTimer{ new TMVA::Timer{boost::numeric_cast<int>(numberOfEvents), "Running over dataset ...", false}}; 
        lEventTimer->DrawProgressBar(0, "");
    
        totalEvents += numberOfEvents;
        for (Long64_t i{0}; i < numberOfEvents; i++) {

            lEventTimer->DrawProgressBar(i,"");

            event.GetEntry(i);      

            for (Int_t k{0}; k < event.nGenPar; k++) {
                const Int_t pdgId    { std::abs(event.genParId[k]) };
		const Int_t status   { event.genParStatus[k] };
		const Int_t motherId { event.genParMotherId[k] };
		const Int_t daughters { event.genParNumDaughters[k] };
		const bool isOwnParent { pdgId == motherId ? true : false };
	  
		if ( daughters == 0 && (status == 1 || status == 2 || status == 71 || status == 72) ) { 
                    histPdgId->Fill(pdgId);
                    histPdgId_2->Fill(pdgId);
                    pdgIdMap[pdgId]++;
                }
		if (status == 1 && daughters == 0) histPdgIdStatus1->Fill(pdgId);
		if (status == 2 && daughters == 0) histPdgIdStatus2->Fill(pdgId);
		if ((status == 61 && status == 62 || status == 63) && daughters == 0) histPdgIdStatus6X->Fill(pdgId);
		if ((status == 71 || status == 72 || status == 74) && daughters == 0) histPdgIdStatus7X->Fill(pdgId);


		//if ( !daughters ) {
		//std::cout << "pdgId / mother / nDaughers / status: " << std::endl;
	        //std::cout << pdgIdCode( pdgId, true ) << " / " << pdgIdCode( motherId, true ) << " / " << daughters << " / " << pythiaStatus( status ) << std::endl;}
	    }
//	    std::cout << std::endl;
	}
    }
    
    std::cout << std::endl;
    std::cout << "Total no. of events:\t\t\t" << totalEvents << std::endl;
    std::cout << std::endl;

    // Do scalable histograms
    int nPdgIds = pdgIdMap.size(); // number of different pdgIds

    TH1I* h_pdgId{new TH1I{"h_pdgId", "Final state content", nPdgIds, 0, Double_t(nPdgIds)}};

    uint binCounter {1};
    for (auto it = pdgIdMap.begin(); it != pdgIdMap.end(); ++it) {
//        std::cout << "Add " << it->second << " to bin " << binCounter << " for pdgId " << it->first << std::endl;
        h_pdgId->SetBinContent(binCounter, it->second);
        const char *label = ( pdgIdCode(it->first, false) ).c_str();
        const char *label2 = ( pdgIdCode(it->first, true) ).c_str();
//        std::cout << "label : " << label2 << std::endl;
        h_pdgId->GetXaxis()->SetBinLabel(binCounter, label);
        binCounter++;
    }

    TFile* outFile{new TFile{outFileString.c_str(), "RECREATE"}};
    outFile->cd();

    histPdgId->Write();
    histPdgId_2->Write();
    histPdgIdStatus1->Write();
    histPdgIdStatus2->Write();
    histPdgIdStatus6X->Write();
    histPdgIdStatus7X->Write();
    h_pdgId->Write();    

    outFile->Close();

//    std::cout << "Max nGenPar: " << maxGenPars << std::endl;    
 
    std::cout << "\nFinished." << std::endl;
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
      case 22 : particle += unicode ? "\u03B3" : "#gamma"; break;
      case 23 : particle += "Z_{0}"; break;
      case 24 : particle += "W"; break;
      case 25 : particle += "h_{0}"; break;

      case 111 : particle += unicode ? "\u03C00" : "#pi^{0}"; break;
      case 211 : particle += unicode ? "\u03C0" : "#pi^{+}"; break;
      case 113 : particle += unicode ? "\u03C10" : "#rho(770)^{0}"; break;
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
      case 10441: particle += unicode ? "\u03C7_C0(1P)" : "#chi_{c0}(1P)"; break;
      case 20313: particle += unicode ? "K0_1 (1400)" : "K_{1}(1400)^0}"; break;
      case 20213: particle += unicode ? "a1" : "a_{1} (1260)^{+}"; break;

      case 9000006 : particle += unicode ? "S" : "ABCD"; break;
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

