from  optparse  import OptionParser
import array, pickle
import numpy as np
import itertools, math
import ROOT as rt
#canvas style file                                                                                                     
import rootlogon
rootlogon.style()

parser = OptionParser()

parser.add_option("-s", "--signal", dest="sig_file",
                  help="input signal root file",default="signal.root",
                  action="store",type="string")

parser.add_option("-b", "--background", dest="bkg_file",
                  help="input background root file",default="qcd.root",
                  action="store",type="string")

parser.add_option("-t", "--tree", dest="tree",
                  help="tree inside the input file to use",default="jets",
                  action="store",type="string")

(options, args) = parser.parse_args()

class analysis:
    def __init__(self, signal_file, bkg_file):
        self.signal_file = signal_file
        self.bkg_file = bkg_file
    
    def build_jetcollections(self):

        print "-- building jet collections:"

        self.signal_jets = jet_collection(signal_file)
        self.bkg_jets = jet_collection(bkg_file)


    def run_disc_scan(self, disc_scan):
        
        
class disc_scan:

    # name of the discriminant
    # number of configurable parameters for the discriminant
    # par_range a list of [x^1_min,x^2 
    def __init__(self, name, npars, par_ranges):
        self.name = name
        self.npars = npars
        self.par_ranges = par_ranges

        #sanity check on the parameter_ranges

        #build a direct product of the parameters
        self.grid  = list(itertools.product(*par_ranges))

class jet_collection:

    def __init__(self, file):
        self.file = file
        self.jetlist = []
        self.build_jets() 

    def build_discriminants(self, pars)
        
    def build_jets(self):
        print "-- building individual jets"

        tree.Draw(">>iterlist", "", "entrylist")
        iev = 0

        while iev < tree.GetEntries():
            if iev % 500 == 0: print "Filling Jets...",iev

            iev += 1
            entry = itlist.Next()
            tree.GetEntry(entry)                                                
            
            for jj in range(tree.nCaloJets):
                jetvars = {}
                jetID = tree.jetID[jj]
                jetvars["jetIPSigLogSum2D"] = tree.jetIPSigLogSum2D[jj]
                jetvars["jetMedianIPSig2D"] = tree.jetMedianIPSig2D[jj]
                jetvars["jetELogIPSig2D"] = tree.jetELogIPSig2D[jj]

                thisJet = jet(jetid, jetvars)
                #add it to the list of jets
                self.jetlist.append(thisJet)
                
class jet:
    #jetid = unique id of the jet
    #jet vars = dictionary look up for each parameter used for the 
    def __Init__(self, jetid, jetvars):
        self.jetid = jetid 
        self.jetvars = jetvars
        self.disc_dict = {}

    #method of grabbing a specific discriminant
    # disc : string name corresponding to metric
    # pars : list of parameters for generating the value
    def fill_disc(self, disc, pars):
        if disc == "metric":
             disc_dict[disc] = get_metric_val1(pars)
        else:
            return "INVALID CHOICE OF METRIC NAME:", disc 

    #simple spatial metric score
    def get_metric_val1(self, pars):
        (w1, w2, w3) = pars

        x1 = jetvars["jetELogIPSig2D"]
        x2 = jetvars["jetMedianIPSig2D"]
        x3 = jetvars["jetIPSigLogSum2D"]

        return w1*x1 + w2*x2 + w3*x3
        
ana = analysis(options.sig_file, options.bkg_file)

ana.build_jetcollections()
