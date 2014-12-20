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


parser.add_option("-o", "--output", dest="output",
                  help="output file ",default="jets",
                  action="store",type="string")


(options, args) = parser.parse_args()

class analysis:
    def __init__(self, signal_file, bkg_file):
        self.signal_file = signal_file
        self.bkg_file = bkg_file
    
    def build_jetcollections(self):

        print "-- building jet collections -- "
        self.signal_jets = jet_collection(self.signal_file)
        self.bkg_jets = jet_collection(self.bkg_file)

    def run_disc_scan(self, disc_scan):
        
        print "-- running discriminant scan -- " 

#container for a range of weights to be scaned 
class weight_range:
    def __init__(self, name, wmin, wmax, npoints, scale):
        self.name = name
        self.wmin = float(wmin)
        self.wmax = float(wmax)
        self.npoints = int(npoints)

        self.vector = np.linspace(wmin / scale, wmax / scale, npoints)

#container for the outer product of weight ranges    
class weight_space:
    # name of the corresponding discriminant to the weight space
    # number of configurable parameters for the discriminant
    def __init__(self, discname, nweights, tuple_of_weight_ranges):
        self.discname = discname
        self.nweights = nweights
        self.tuple_of_weight_ranges = tuple_of_weight_ranges

        #sanity check on the weight ranges
        if int(nweights) != len(tuple_of_weight_ranges):
            print "[ERROR] number of weights does not equal number of weight ranges. "
              
        # build a direct product of the parameters
        print tuple_of_weight_ranges
        self.grid  = list(itertools.product(*tuple_of_weight_ranges))

class jet_collection:

    def __init__(self, filename):
        self.filename = filename
        self.jetlist = []
        self.build_jets()     
        self.disc_list = []
        
    def build_grid_tree(self, tree_name):
        
        tree = rt.TTree(tree_name, tree_name)    
    
            # Create a struct for the run information
        line_to_process = "struct MyStruct{ Int_t id;Int_t genmatch;Float_t pt;Float_t eta;Float_t phi;Float_t ipsiglog;Float_t ipelog;Float_t ipmed;Int_t genmatch; Int_t nGrid;"

        for disc in self.disc_list:
            line_to_process += " Float_t %s[1000];" % disc
            line_to_process += " Int_t %sID[1000];" % disc
        
        line_to_process += "};"

        print line_to_process

        rt.gROOT.ProcessLine(line_to_process)

        from ROOT import MyStruct
        
        s = MyStruct()
        
        tree.Branch("jetid",rt.AddressOf(s,"id"),"jetid/I")
        tree.Branch("nGrid",rt.AddressOf(s,"nGrid"),"nGrid/I")

        tree.Branch("pt",rt.AddressOf(s,"pt"),"pt/F")
        tree.Branch("eta",rt.AddressOf(s,"eta"),"eta/F")
        tree.Branch("phi",rt.AddressOf(s,"phi"),"phi/F")

        tree.Branch("ipsiglog",rt.AddressOf(s,"ipsiglog"),"ipsiglog/F")
        tree.Branch("ipmed",rt.AddressOf(s,"ipmed"),"ipmed/F")
        tree.Branch("ipelog",rt.AddressOf(s,"ipelog"),"ipelog/F")
        tree.Branch("genmatch",rt.AddressOf(s,"genmatch"),"genmatch/F")

        tree.Branch(disc, eval("s.%s" % disc) ,"%s[nGrid]/F" % disc )
        tree.Branch(disc+"ID", eval("s.%sID" % disc) ,"%sID[nGrid]/I" % disc )

        # fill  the information common to every discriminant
        s.id = s.genmatch = s.pt = s.eta = s.phi = 0
        s.ipsiglog = s.ipmed = s.ipelog = 0

        #fill the tree in terms of jets
        for jet in self.jetlist:
            s.ipsiglog = jet.jetvars["jetIPSigLogSum2D"] 
            s.ipmed = jet.jetvars["jetMedianIPSig2D"] 
            s.ipelog = jet.jetvars["jetELogIPSig2D"] 
            s.pt = jet.jetvars["caloJetPt"]
            s.eta = jet.jetvars["caloJetEta"]
            s.phi = jet.jetvars["caloJetPhi"]
            s.id = jet.jetvars["jetID"] 
            s.genmatch = jet.jetvars["genMatch"]

            for disc in self.disc_list:            

                jet_dict = jet.disc_dict
                gridID = 0
                
                s.nGrid = len(jet_dict)
                
                disc_list = []
                discID_array = array.array("i", range(1,s.nGrid+1))

                #loop over every calculated disc and assign it a unique ID
                for key in jet_dict.keys():
                    (disc, wvector) = key 
                    val = jet_dict[key] 

                    disc_list.append(val)
                
                #print "disc values ", disc_array

                disc_array = array.array("f",disc_list)
                exec("s.%s = disc_array" % disc) 
                exec("s.%sID = discID_array" % disc)
                                

            #fill the tree for each jet
            tree.Fill()

        return tree

    def fill_discriminant(self, disc, weight_space):

        self.disc_list.append(disc)

        for jet in self.jetlist:
            jet.fill_disc(disc, weight_space)

        
    def build_jets(self):

        print "-- building individual jets"

        input_file = rt.TFile(self.filename)
        tree = input_file.Get(options.tree)

        tree.Draw(">>iterlist", "", "entrylist")
        iev = 0
        
        itlist = rt.gDirectory.Get("iterlist")

        while iev < tree.GetEntries():
            if iev % 500 == 0: print self.filename, "Filling Jets...",iev

            if iev > 500: continue 

            iev += 1
            entry = itlist.Next()
            tree.GetEntry(entry)                                                
            
            for jj in range(tree.nCaloJets):
                jetvars = {}
                jetID = tree.jetID[jj]
                jetvars["jetIPSigLogSum2D"] = tree.jetIPSigLogSum2D[jj]
                jetvars["jetMedianIPSig2D"] = tree.jetMedianIPSig2D[jj]
                jetvars["jetELogIPSig2D"] = tree.jetELogIPSig2D[jj]
                jetvars["caloJetPt"] = tree.caloJetPt[jj]
                jetvars["caloJetEta"] = tree.caloJetEta[jj]
                jetvars["caloJetPhi"] = tree.caloJetPhi[jj]
                jetvars["jetID"] = tree.jetID[jj]
                jetvars["genMatch"] = tree.genMatch[jj]

                thisJet = jet(tree.jetID[jj], jetvars)

                #add it to the list of jets
                self.jetlist.append(thisJet)
                
class jet:
    #jetid = unique id of the jet
    #jet vars = dictionary look up for each parameter used for the 
    def __init__(self, jetid, jetvars):
        self.jetid = jetid 
        self.jetvars = jetvars #dictionary of thet variable values
        self.disc_dict = {} #dictionary of ["discriminant", weight_vector] -> calculated value
        
    # disc : string name corresponding to metric
    # weight_space: class of outer product of the ranges of the weights

    def get_disc_val(self, disc, wvector):

        return disc_dict[disc, wvector]

    # most naive descriminant 
    def calc_metric_val1(self, pars):
        (w1, w2, w3) = pars

        x1 = self.jetvars["jetELogIPSig2D"]
        x2 = self.jetvars["jetMedianIPSig2D"]
        x3 = self.jetvars["jetIPSigLogSum2D"]

        return w1*x1 + w2*x2 + w3*x3

    def fill_disc(self, disc, weight_space):
        
        wgrid = weight_space.grid

        #go over each point in the weight space            
        for wvector in wgrid:            
            if disc == "metric":
                self.disc_dict[disc, wvector] = self.calc_metric_val1(wvector)
            else:
                return "INVALID CHOICE OF DISCRIMINANT -- NAME:", disc 


output_file = rt.TFile(options.output, "RECREATE")
        
ana = analysis(options.sig_file, options.bkg_file)

ana.build_jetcollections()

#build the space for the metric
w1_range = weight_range("elogipsig", 0, 1, 5, 20.)
w2_range = weight_range("jetmedianipsig", 0, 1, 5, .05)
w3_range = weight_range("ipsiglogsum", 0, 1, 5, 10.)
metric_weight_space = weight_space("metric", 3, (w1_range.vector, w2_range.vector, w3_range.vector))

ana.signal_jets.fill_discriminant("metric", metric_weight_space)
ana.bkg_jets.fill_discriminant("metric", metric_weight_space)


sig_tree = ana.signal_jets.build_grid_tree("sig")
bkg_tree = ana.bkg_jets.build_grid_tree("bkg")
output_file.cd()
sig_tree.Write()
bkg_tree.Write()
