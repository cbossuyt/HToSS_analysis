#A tool to pull plots from an input ROOT file

import ROOT
import subprocess
import sys

def main():

  # input file - location set by user
  infile = ROOT.TFile.Open("exampleInputs/ggH/HToSS_SmuonHadronFiltered/MS2_ctauS10/skimFile0.root")

  # Histogram declarations
  genParPt  = ROOT.TH1F("h_genParPt", "genPar p_{T}", 1000, 0., 1000.)
  genParEta = ROOT.TH1F("h_genParEta", "genPar #eta", 200, -7., 7.)
  genParPhi = ROOT.TH1F("h_genParPhi", "genPar #phi", 100, -3.5, 3.5)
  genParE   = ROOT.TH1F("h_genParE", "genPar energy", 1000, 0., 1000.)

  #Loop over each event in the input file
  for event in infile.tree :

    #Fill histogram
    for i in range (event.nGenPar):
      genParPt.Fill(event.genParPt[i])
      genParEta.Fill(event.genParEta[i])
      genParPhi.Fill(event.genParPhi[i])
      genParE.Fill(event.genParE[i])

  # Save each histogram in its own ROOT file
  genParPt.SaveAs("genParPt.root")
  genParEta.SaveAs("genParEta.root")
  genParPhi.SaveAs("genParPhi.root")
  genParE.SaveAs("genParE.root")


if __name__ == "__main__":
    main()

