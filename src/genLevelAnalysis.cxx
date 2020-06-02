#include "AnalysisEvent.hpp"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMVA/Timer.h"
#include "TTree.h"

#include <boost/filesystem.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <boost/range/iterator_range.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

std::string pythiaStatus (const Int_t status);
std::string pdgIdCode (const Int_t status);

int main(int argc, char* argv[])
{
    std::string inputDir{};
    std::string outFileString{"plots/distributions/output.root"};
    bool is2016;
    Long64_t nEvents;

//    int maxGenPars {0};

    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()("help,h", "Print this message.")(
        "indir,i",
        po::value<std::string>(&inputDir)->required(),
        "Input folder for nTuples.")(
        "outfile,o",
        po::value<std::string>(&outFileString)->default_value(outFileString),
        "Output file for plots.")(
        ",n",
        po::value<Long64_t>(&nEvents)->default_value(0),
        "The number of events to be run over. All if set to 0.")(
        "2016", po::bool_switch(&is2016), "Use 2016 conditions (SFs, et al.).");
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

    std::vector<TTree*> inputTrees;

    if (boost::filesystem::is_directory(inputDir))
    {
        const long int count{
            std::distance(boost::filesystem::directory_iterator{inputDir},
                          boost::filesystem::directory_iterator())};

        TMVA::Timer lTimer{
            boost::numeric_cast<int>(count), "Attaching files to TTree", false};
        Int_t lCounter{1};

        lTimer.DrawProgressBar(0, "");

        for (const auto& file : boost::make_iterator_range(
                 boost::filesystem::directory_iterator{inputDir}, {}))
        {
            TFile* inputFile{new TFile{file.path().string().c_str()}};
            TTree* lTempTree{dynamic_cast<TTree*>(inputFile->Get("tree"))};
            inputTrees.emplace_back(lTempTree);

            lTimer.DrawProgressBar(lCounter++, "");
        }
    }
    else
    {
        std::cout << "ERROR: " << inputDir << "is not a valid directory"
                  << std::endl;
        return 1;
    }

    std::cout << std::endl;
    std::cout << "Attached all files to TTree!" << std::endl;

// status == 1 for final state particles
// status == 2 for a decayed Standard Model hadron or tau or mu lepton, excepting virtual intermediate states thereof (i.e. the particle must undergo a normal decay, not e.g. a shower branching);
// status == 61-63 for particles produced by beam-remnant treatment
// status == 71 for partons in preparation of hadronization process and 72+74 (but exclude particles who are their own parent)

    TH1D* histPdgId{new TH1D{"histPdgId", "Final state content", 501, 0., 500.5}};
    TH1D* histPdgIdStatus1{new TH1D{"histPdgIdStatus1", "Final state content", 501, -.5, 500.5}};
    TH1D* histPdgIdStatus2{new TH1D{"histPdgIdStatus2", "Decayed SM hadron or tau or mu", 501, -0.5, 500.5}};
    TH1D* histPdgIdStatus6X{new TH1D{"histPdgIdStatus6X", "Beam remnants", 501, -0.5, 500.5}};
    TH1D* histPdgIdStatus7X{new TH1D{"histPdgIdStatus7X", "Oartons in preparation of hadronization process", 501, -0.5, 500.5}};

    TMVA::Timer* lTimer{
        new TMVA::Timer{boost::numeric_cast<int>(inputTrees.size()),
                        "Running over trees",
                        false}};

    lTimer->DrawProgressBar(0, "");

    Int_t lCounter{1};

    // Event counters
    int totalEvents{0};

    for (std::vector<TTree*>::const_iterator lIt = inputTrees.begin(); lIt != inputTrees.end(); ++lIt) {
      AnalysisEvent* lEvent{new AnalysisEvent{true, *lIt, is2016}};

      Long64_t lNumEvents{(*lIt)->GetEntries()};
      if (nEvents && nEvents < lNumEvents) lNumEvents = nEvents;
      
      totalEvents += lNumEvents;
      
      for (Int_t j{0}; j < lNumEvents; j++) {
	(*lIt)->GetEvent(j);
	
	
//        if (maxGenPars < lEvent->nGenPar) maxGenPars = lEvent->nGenPar;
//        std::cout << "nGenPar: " << lEvent->nGenPar << std::endl;

	for (Int_t k{0}; k < lEvent->nGenPar; k++) {
	  
	  const Int_t pdgId    { lEvent->genParId[k] };
	  const Int_t status   { lEvent->genParStatus[k] };
	  const Int_t motherId { lEvent->genParMotherId[k] };
	  const Int_t daughters { lEvent->genParNumDaughters[k] };
	  const bool isOwnParent { pdgId == motherId ? true : false };
	  
	  if ( daughters == 0 && (status == 1 || status == 2 || status == 71 || status == 72) )  histPdgId->Fill(pdgId);
	  if (status == 1 && daughters == 0) histPdgIdStatus1->Fill(pdgId);
	  if (status == 2 && daughters == 0) histPdgIdStatus2->Fill(pdgId);
	  if ((status == 61 && status == 62 || status == 63) && daughters == 0) histPdgIdStatus6X->Fill(pdgId);
	  if ((status == 71 || status == 72 || status == 74) && daughters == 0) histPdgIdStatus7X->Fill(pdgId);

//          if ( !daughters ) {
//	     std::cout << "pdgId / mother / nDaughers / status: " << std::endl;
//	     std::cout << pdgIdCode( pdgId ) << " / " << pdgIdCode( motherId ) << " / " << daughters << " / " << pythiaStatus( status ) << std::endl;}
	}
      }
//        std::cout << std::endl;
        lTimer->DrawProgressBar(lCounter++, "");
    }
    std::cout << std::endl;
    std::cout << "Total no. of events:\t\t\t" << totalEvents << std::endl;
    std::cout << std::endl;

    TFile* outFile{new TFile{outFileString.c_str(), "RECREATE"}};

    histPdgId->Write();
    histPdgIdStatus1->Write();
    histPdgIdStatus2->Write();
    histPdgIdStatus6X->Write();
    histPdgIdStatus7X->Write();

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

std::string pdgIdCode (const Int_t parId) {

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
      case 7 : particle += "b\'"; break;
      case 8 : particle += "t\t"; break;
      case 9 : particle += "g"; break;

      case 11 : particle += "e"; break;
      case 12 : particle += "\u03BD_e"; break;
      case 13 : particle += "mu"; break;
      case 14 : particle += "\u03BD_mu"; break;
      case 15 : particle += "\u03C4"; break;
      case 16 : particle += "\u03BD_\u03C4"; break;
      case 17 : particle += "\u03C4\'"; break;
      case 18 : particle += "\u03BD_\u03C4\'"; break;

      case 21 : particle += "g"; break;
      case 22 : particle += "\u03B3"; break;
      case 23 : particle += "Z0"; break;
      case 24 : particle += "W"; break;
      case 25 : particle += "h0"; break;

      case 111 : particle += "\u03C00"; break;
      case 211 : particle += "\u03C0"; break;
      case 113 : particle += "\u03C10"; break;
      case 213 : particle += "\u03C1"; break;

      case 130 : particle += "K0_L"; break;
      case 310 : particle += "K0_S"; break;
      case 311 : particle += "K0"; break;
      case 321 : particle += "K"; break;

      case 313 : particle += "K*0"; break;
      case 323 : particle += "K*"; break;

      case 221 : particle += "\u03B7"; break;
      case 331 : particle += "\u03B7\'"; break;
      case 223 : particle += "\u03C9"; break;
      case 333 : particle += "\u03C6"; break;

      case 411 : particle += "D"; break;
      case 421 : particle += "D0"; break;
      case 413 : particle += "D*"; break;
      case 423 : particle += "D*0"; break;
      case 431 : particle += "D_S"; break;
      case 433 : particle += "D*_S"; break;
      case 443 : particle += "J/\u03C8"; break;
      case 445 : particle += "\u03C7_c2(1P)"; break;
 
      case 511 : particle += "B0"; break;
      case 521 : particle += "B"; break;
      case 513 : particle += "B*0"; break;
      case 523 : particle += "B*"; break;
      case 531 : particle += "B0_S"; break;
      case 533 : particle += "B*0_S"; break;

      case 1101 : particle += "(dd)0"; break;
      case 2101 : particle += "(ud)0"; break;
      case 2103 : particle += "(ud)1"; break;
      case 2203 : particle += "(uu)1"; break;

      case 2212 : particle += "p"; break;
      case 2112 : particle += "n"; break;
      case 2224 : particle += "\u0394++"; break;
      case 2214 : particle += "\u0394+"; break;
      case 2114 : particle += "\u03940"; break;
      case 1114 : particle += "\u0394-"; break;

      case 3122 : particle += "\u0394"; break;
      case 3222 : particle += "\u03A3"; break;
      case 3224 : particle += "\u03A3*"; break;
      case 3212 : particle += "\u03A30"; break;
      case 3322 : particle += "\u03A30"; break;
      case 3112 : particle += "\u03A3-"; break;

      case 3324 : particle += "\u039E*0"; break;

      case 4214 : particle += "\u03A3*_C"; break;
      case 4222 : particle += "\u03A3_C"; break;
      case 4122 : particle += "\u039BC"; break;
      case 4114 : particle += "\u03A3*0_C"; break;
      case 4224 : particle += "\u03A3+C"; break;

      case 5122 : particle += "\u039B0_b"; break;
      case 5212 : particle += "\u03A30_b"; break;
      case 5232 : particle += "\u039E0_b"; break;

      case 10313: particle += "K0_1 (1270)"; break;
      case 10441: particle += "\u03C7_C0(1P)"; break;
      case 20313: particle += "K0_1 (1400)"; break;
      case 20213: particle += "a1"; break;

      case 9000006 : particle += "S"; break;
      default : {particle += std::to_string(std::abs(parId)); std::cout << "UNKNOWN PID: " << parId << std::endl;}
   }

   if (parId == 211 || parId == 213 || parId == 321 || parId == 323 || parId == 411 || parId == 431 || parId == 433 || parId == 521 || parId == 4122 || parId == 20213 || parId == 4214 || parId == 523
   || parId == 3224 || parId == 3222
)  particle += "+";

   if (parId < 0) particle += "-";
   if (parId == -4222 || parId == -4224) particle += "-";
   if (parId == 4222  || parId == 4222) particle += "+";
   return particle;
}

