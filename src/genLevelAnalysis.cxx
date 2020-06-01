#include "AnalysisEvent.hpp"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
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

std::string pythiaStatus (int status);

int main(int argc, char* argv[])
{
    std::string inputDir{};
    std::string outFileString{"plots/distributions/output.root"};
    bool is2016;
    Long64_t nEvents;

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

//    TH1F* histElePt{new TH1F{"histEleGenEta", "Distribution of gen-electron #eta", 500, -2.5, 2.5}};


    TMVA::Timer* lTimer{
        new TMVA::Timer{boost::numeric_cast<int>(inputTrees.size()),
                        "Running over trees",
                        false}};

    lTimer->DrawProgressBar(0, "");

    Int_t lCounter{1};

    // Event counters
    int totalEvents{0};

    for (std::vector<TTree*>::const_iterator lIt = inputTrees.begin(); lIt != inputTrees.end(); ++lIt)
    {
        AnalysisEvent* lEvent{new AnalysisEvent{true, *lIt, is2016}};

        Long64_t lNumEvents{(*lIt)->GetEntries()};
        if (nEvents && nEvents < lNumEvents) lNumEvents = nEvents;

        totalEvents += lNumEvents;

        for (Int_t j{0}; j < lNumEvents; j++)
        {
            (*lIt)->GetEvent(j);
            std::cout << "nGenPar: " << lEvent->nGenPar << std::endl;
            std::cout << "pdgId / mother / status: " << std::endl;
            for (Int_t k{0}; k < lEvent->nGenPar; k++) {
                std::cout << lEvent->genParId[k] << " / " << lEvent->genParMotherId[k] << " / " << pythiaStatus( lEvent->genParStatus[k] ) << std::endl;
            }
            std::cout << std::endl;

//            for (Int_t k{0}; k < lEvent->numElePF2PAT; k++) {
//                histElePt->Fill(lEvent->elePF2PATPT[k]);
//            }
        }
        lTimer->DrawProgressBar(lCounter++, "");
    }

    std::cout << std::endl << std::endl;
    std::cout << "Total no. of events:\t\t\t" << totalEvents << std::endl;
    std::cout << std::endl;

//    TFile* outFile{new TFile{outFileString.c_str(), "RECREATE"}};

//    histElePt->Write();

//    outFile->Close();
    std::cout << "\n Finished." << std::endl;
}

std::string pythiaStatus (int status) {

   std::string message;
   message += std::to_string(status);

   switch (status) {
       case 11 : message += " : beam particles - the event as a whole";
           break;
   }
   switch (status) {
       case 12 : message += " : beam particles - incoming beam";
           break;
   }
   switch (status) {
       case 13 : message += " : beam particles - incoming beam-inside-beam (e.g. gamma inside e)";
           break;
   }
   switch (status) {
       case 14 : message += " : beam particles - outgoing elastically scattered";
           break;
   }
   switch (status) {
       case 15 : message += " : beam particles - outgoing diffractively scattered";
           break;
   }
   switch (status) {
       case 21 : message += " : particles of the hardest subprocess - incoming";
           break;
   }
   switch (status) {
       case 22 : message += " : particles of the hardest subprocess - intermediate (intended to have preserved mass)";
           break;
   }
   switch (status) {
       case 23 : message += " : particles of the hardest subprocess - outgoing";
           break;
   }
   switch (status) {
       case 24 : message += " : particles of the hardest subprocess - outgoing, nonperturbatively kicked out in diffraction";
           break;
   }
   switch (status) {
       case 31 : message += " : particles of subsequent subprocesses - incoming";
           break;
   }
   switch (status) {
       case 32 : message += " : particles of subsequent subprocesses - intermediate (intended to have preserved mass)";
           break;
   }
   switch (status) {
       case 33 : message += " : particles of subsequent subprocesses - outgoing";
           break;
   }
   switch (status) {
       case 34 : message += " : particles of subsequent subprocesses - incoming that has already scattered";
           break;
   }
   switch (status) {
       case 41 : message += " : particles produced by initial-state-showers - incoming on spacelike main branch";
           break;
   }
   switch (status) {
       case 42 : message += " : particles produced by initial-state-showers - incoming copy of recoiler";
           break;
   }
   switch (status) {
       case 43 : message += " : particles produced by initial-state-showers - outgoing produced by a branching";
           break;
   }
   switch (status) {
       case 44 : message += " : particles produced by initial-state-showers - outgoing shifted by a branching";
           break;
   }
   switch (status) {
       case 45 : message += " : particles produced by initial-state-showers - incoming rescattered parton, with changed kinematics owing to ISR in the mother system (cf. status 34)";
           break;
   }
   switch (status) {
       case 46 : message += " : particles produced by initial-state-showers - incoming copy of recoiler when this is a rescattered parton (cf. status 42)";
           break;
   }
   switch (status) {
       case 47 : message += " : particles produced by initial-state-showers - a W or Z gauge boson produced in the shower evolution";
           break;
   }
   switch (status) {
       case 51 : message += " : particles produced by final-state-showers - outgoing produced by parton branching";
           break;
   }
   switch (status) {
       case 52 : message += " : particles produced by final-state-showers - outgoing copy of recoiler, with changed momentum";
           break;
   }
   switch (status) {
       case 53 : message += " : particles produced by final-state-showers - copy of recoiler when this is incoming parton, with changed momentum";
           break;
   }
   switch (status) {
       case 54 : message += " : particles produced by final-state-showers - copy of a recoiler, when in the initial state of a different system from the radiator";
           break;
   }
   switch (status) {
       case 55 : message += " : particles produced by final-state-showers - copy of a recoiler, when in the final state of a different system from the radiator";
           break;
   }
   switch (status) {
       case 56 : message += " : particles produced by final-state-showers - a W or Z gauge boson produced in a shower branching (special case of 51)";
           break;
   }
   switch (status) {
       case 61 : message += " : particles produced by beam-remnant treatment - incoming subprocess particle with primordial kT included";
           break;
   }
   switch (status) {
       case 62 : message += " : particles produced by beam-remnant treatment - outgoing subprocess particle with primordial kT included";
           break;
   }
   switch (status) {
       case 63 : message += " : particles produced by beam-remnant treatment - outgoing beam remnant";
           break;
   }
   switch (status) {
       case 71 : message += " : partons in preparation of hadronization process - copied partons to collect into contiguous colour singlet";
           break;
   }
   switch (status) {
       case 72 : message += " : partons in preparation of hadronization process - copied recoiling singlet when ministring collapses to one hadron and momentum has to be reshuffled";
           break;
   }
   switch (status) {
       case 73 : message += " : partons in preparation of hadronization process - combination of very nearby partons into one";
           break;
   }
   switch (status) {
       case 74 : message += " : partons in preparation of hadronization process - combination of two junction quarks (+ nearby gluons) to a diquark";
           break;
   }
   switch (status) {
       case 75 : message += " : partons in preparation of hadronization process - gluons split to decouple a junction-antijunction pair";
           break;
   }
   switch (status) {
       case 76 : message += " : partons in preparation of hadronization process - partons with momentum shuffled to decouple a junction-antijunction pair";
           break;
   }
   switch (status) {
       case 77 : message += " : partons in preparation of hadronization process - temporary opposing parton when fragmenting first two strings in to junction (should disappear again)";
           break;
   }
   switch (status) {
       case 78 : message += " : partons in preparation of hadronization process - temporary combined diquark end when fragmenting last string in to junction (should disappear again)";
           break;
   }
   switch (status) {
       case 81 : message += " : primary hadrons produced by hadronization process - from ministring into one hadron";
           break;
   }
   switch (status) {
       case 82 : message += " : primary hadrons produced by hadronization process - from ministring into two hadron";
           break;
   }
   switch (status) {
       case 83 : message += " : primary hadrons produced by hadronization process - from normal string (the difference between the two is technical, whether fragmented off from the top of the string system or from the bottom, useful for debug only)";
           break;
   }
   switch (status) {
       case 84 : message += " : primary hadrons produced by hadronization process - from normal string (the difference between the two is technical, whether fragmented off from the top of the string system or from the bottom, useful for debug only)";
           break;
   }
   switch (status) {
       case 85 : message += " : primary hadrons produced by hadronization process - primary produced hadrons in junction fragmentation of the first two string legs in to the junction, in order of treatment";
           break;
   }
   switch (status) {
       case 86 : message += " : primary hadrons produced by hadronization process - primary produced hadrons in junction fragmentation of the first two string legs in to the junction, in order of treatment";
           break;
   }
   switch (status) {
       case 91 : message += " : particles produced in decay process, or by Bose-Einstein effects - normal decay products";
           break;
   }
   switch (status) {
       case 92 : message += " : particles produced in decay process, or by Bose-Einstein effects - decay products after oscillation";
           break;
   }
   switch (status) {
       case 93 : message += " : particles produced in decay process, or by Bose-Einstein effects - decay handled by external program, normally or with oscillation";
           break;
   }
   switch (status) {
       case 94 : message += " : particles produced in decay process, or by Bose-Einstein effects - decay handled by external program, normally or with oscillation";
           break;
   }
   switch (status) {
       case 99 : message += " : particles produced in decay process, or by Bose-Einstein effects - particles with momenta shifted by Bose-Einstein effects";
           break;
   }
   switch (status) {
       case 101 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - when a string system contains two such long-lived particles, the system is split up by the production of a new q-qbar pair";
           break;
   }
   switch (status) {
       case 102 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - partons rearranged from the long-lived particle end to prepare for fragmentation from this end";
           break;
   }
   switch (status) {
       case 103 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - intermediate half-R-hadron ormed when a colour octet particle (like the gluino) has been fragmented on one side, but not yet on the other";
           break;
   }
   switch (status) {
       case 104 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - an R-hadron";
           break;
   }
   switch (status) {
       case 105 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - partons or particles formed together with the R-hadron during the fragmentation treatment";
           break;
   }
   switch (status) {
       case 106 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - subdivision of an R-hadron into its flavour content, with momentum split accordingly, in preparation of the decay of the heavy new particle, if it is unstable";
           break;
   }
   switch (status) {
       case 107 : message += " : particles in the handling of R-hadron production and decay, i.e. long-lived (or stable) particles with very heavy flavour - two temporary leftover gluons joined into one in the formation of a gluino-gluon R-hadron";
           break;
   }

   if (status > 107 ) message += " :  : unknown status code ";

   return message;
}
