from  optparse  import OptionParser
import copy
import array, pickle
import numpy as np
import CMS_lumi
import itertools, math
import ROOT as rt
#canvas style file                                                                                                     
import rootlogon
rootlogon.style()

N_BEST_ROC_POINTS = 20
MIN_BEST = .6
MAX_BEST = 1.0
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
parser.add_option( "--roc", dest="do_roc",
                  help="build the roc curves ",default=False,
                  action="store_true")
parser.add_option( "--gidout", dest="output_GID",
                  help="name of output file to give combinations of weights ",
                  action="store", type="string", default="GID.txt")


(options, args) = parser.parse_args()

class analysis:
    def __init__(self, signal_file, bkg_file):
        self.signal_file = signal_file
        self.bkg_file = bkg_file
        self.disc_list = [] 
        self.all_roc_curves = []
        self.signal_jets = self.bkg_jets = None 
        self.nGID = 0

    def fill_discriminant(self, disc, metric_weight_space):
        self.disc_list.append(disc)
        self.signal_jets.fill_discriminant(disc, metric_weight_space)
        self.bkg_jets.fill_discriminant(disc, metric_weight_space)
        self.nGID = len(metric_weight_space.grid)
        
    def get_trees(self):
        sig_tree = self.signal_jets.build_grid_tree("sig")
        bkg_tree = self.bkg_jets.build_grid_tree("bkg")
        
        return (sig_tree, bkg_tree)

    def build_jetcollections(self):
        print "-- building jet collections -- "
        self.signal_jets = jet_collection(self.signal_file)
        self.bkg_jets = jet_collection(self.bkg_file)

    def build_best_roc(self):        
        # build a new histogram with the best possible points
        best_roc = rt.TH1F("best_roc", "best_roc", N_BEST_ROC_POINTS, MIN_BEST, MAX_BEST)        
        for curve in self.all_roc_curves:
            hist = curve.sig_hist             
            for bin in range(hist.GetXaxis().GetNbins()):
                if best_roc.GetBinContent(bin) < hist.GetBinContent(bin):
                    best_roc.SetBinContent(bin, hist.GetBinContent(bin))

        #turn it into a tgraph
        x = []
        y = []
        
        for bin in range(best_roc.GetXaxis().GetNbins()):
            x.append(best_roc.GetXaxis().GetBinCenter(bin))
            y.append(best_roc.GetBinContent(bin))
            
        x_ar = array.array("d", x)
        y_ar = array.array("d", y)
        
        tgraph = rt.TGraph(len(x), x_ar, y_ar)
        tgraph.SetTitle("combination ROC")
        tgraph.SetName("best_roc")
        tgraph.GetXaxis().SetTitle("Displaced Jet Efficiency")
        tgraph.GetYaxis().SetTitle("QCD Jet Rejection")

        return tgraph

            
    def build_roc_canvas(self):

        self.canvas = rt.TCanvas("roc_curves","roc curves")

        color = 0
        colors = [rt.kRed, rt.kBlue]

        first_graph = self.tgraphs[0]
        first_graph.SetLineColor(1)
        first_graph.SetFillColor(0)
        first_graph.SetLineWidth(2)

        first_graph.GetXaxis().SetTitle("Displaced Jet Efficiency")
        first_graph.GetXaxis().SetRangeUser(.2,1)
        first_graph.GetYaxis().SetTitle("QCD Jet Rejection")
        first_graph.GetYaxis().SetRangeUser(.1,1)
        first_graph.Draw()

        self.canvas.SetGridy(1)
        self.canvas.SetGridx(1)

        for graph in self.tgraphs[1:]:
            graph.SetLineColor(colors[color])
            graph.SetFillColor(0)
            graph.SetLineWidth(2)

            graph.Draw("same")
            color+=1
            if color > 100:
                color = 0

        leg = self.canvas.BuildLegend()
        leg.SetFillColor(0)
        
        CMS_lumi.CMS_lumi(self.canvas, 4, 0)        

        return self.canvas          

    def generate_tgraphs(self):

        self.tgraphs = []

        print "-- generating tgraphs --"
        for curve in self.all_roc_curves:

            # fill the arrays from the dictionary
            x_array = array.array("d", curve.x)
            y_array = array.array("d", curve.y)

            # identifying information
            name = str(curve.disc) + "_" + str(curve.gridID)  
            label = "Disc: " + str(curve.disc) + " gridID:" + str(curve.gridID) 
        
            # make the tgraphs recognizable in the output
            graph = rt.TGraph(len(curve.x), copy.copy(x_array), copy.copy(y_array) )
            graph.SetTitle(label)
            graph.SetName(name)

            self.tgraphs.append(graph)            

        print "-- finished generating tgraphs" 

        return self.tgraphs

    def build_roc_curves(self, sig_tree, bkg_tree, xmin, xmax, npoints):
        print "-- Building ROC curves" 

#        scan_points = list(np.linspace(xmin, xmax, npoints))
#        if explicit_array !=  None:
#            scan_points = explicit_array 
#        print scan_points

        #roc curve for each discriminant
        print "-- Looping Discriminants -- " 
        for disc in self.disc_list:

            #total number of jets
#            nsig_tot = sig_tree.GetEntries("genMatch > 0 && %sID == 1" % disc)
            nsig_tot = sig_tree.GetEntries("caloGenMatch == 0 && %sID == 1" % disc)
            nbkg_tot = bkg_tree.GetEntries(" %sID == 1" % disc)        

            print "-- Looping GID -- " 
            # one curve per grid ID
            for gid in range(1, self.nGID+1):
                print "-- %s, GID =  %i of %s" % (disc, gid, self.nGID)
                
                curve = roc_curve(disc, gid)
                scan_points = list(np.linspace(xmin, xmax, npoints))

                # one point per threshold
                for threshold in scan_points:
                    sig_cut = "caloGenMatch > 0 && %sID == %i && %s > %s" % (disc, gid, disc, threshold)
                    bkg_cut = "%sID == %i && %s > %s" % (disc, gid, disc, threshold)
                    
                    nsig_pass = sig_tree.GetEntries(sig_cut)
                    nbkg_pass = bkg_tree.GetEntries(bkg_cut) 

                    eff_sig = float(nsig_pass) / float(nsig_tot)
                    eff_bkg = float(nbkg_pass) / float(nbkg_tot)
                    


                    curve.add_point(threshold, eff_sig, 1 - eff_bkg) 

                self.all_roc_curves.append(curve)

    def get_fischer_weights(self, tuple_of_weight_ranges):
        print "-- building fischer weights --" 

        #get array indexed by (jet1, jet2, ---) with jet1 = ( var1, var2, var3...)
        signal_vars = self.signal_jets.get_jet_var_index(tuple_of_weight_ranges)
        bkg_vars = self.bkg_jets.get_jet_var_index(tuple_of_weight_ranges)
            
        class1 = np.array(signal_vars, dtype=np.float)
        class2 = np.array(bkg_vars, dtype=np.float)
        
        mean1 = np.mean(class1, axis=0)
        mean2 = np.mean(class2, axis=0)
    
        #calculate variance within class
        class1_dot = np.dot((class1-mean1).T , class1-mean1) 
        class2_dot = np.dot((class2-mean2).T , class2-mean2) 
        Sw = class1_dot + class2_dot

        #calculate weights which maximize linear separation
        w = np.dot(np.linalg.inv(Sw), (mean2-mean1))
    
        neg_w = -1 * w

        return neg_w
    
#container for a range of weights to be scaned 
class weight_range:
    def __init__(self, name, wmin, wmax, npoints, scale, explicit=None):
        self.name = name
        self.wmin = float(wmin)
        self.wmax = float(wmax)
        self.npoints = int(npoints)
        self.scale = scale

        self.vector = np.linspace(wmin / float(scale), wmax / float(scale), npoints)

        if explicit != None:
            print "EXPLICIT RANGE SPECIFIED -- ", name, explicit
            self.vector = explicit
            self.npoints = len(explicit)
            self.wmin = min(explicit)
            self.wmax = max(explicit)
            self.scale = 1.0

#container for the outer product of weight ranges    
class weight_space:
    # name of the corresponding discriminant to the weight space
    # number of configurable parameters for the discriminant
    def __init__(self, discname, nweights, tuple_of_weight_ranges):
        self.discname = discname
        self.nweights = nweights
        self.tuple_of_weight_ranges = tuple_of_weight_ranges 
        self.tuple_of_array_ranges = map( lambda x: x.vector, tuple_of_weight_ranges)

        #sanity check on the weight ranges
        if int(nweights) != len(tuple_of_weight_ranges):
            print "[ERROR] number of weights does not equal number of weight ranges. "
              
        # build a direct product of the parameters
        print tuple_of_weight_ranges
        self.grid  = list(itertools.product(*self.tuple_of_array_ranges))


    def write_output(self, name):
        
        outfile = open(name, "w")

        labels = ""
        scales = []

        print "SCALES: ", scales 

        for range_var in self.tuple_of_weight_ranges:
            labels += (range_var.name + "\t" )
            scales.append(range_var.scale)

        outfile.write(labels + "\n")

        gid = 1

        for point in self.grid:
            value_string = "GID: %i -- " % gid
            plist = list(point)

            for ii in range(len(plist)):
                val = plist[ii]
                value_string +=  str(val * scales[ii]) + "\t" 
            
            outfile.write(value_string + "\n")
            gid +=1           
            
class roc_curve:
    def __init__(self, disc, gridID):
        self.disc = disc
        self.eff_dict = {} 
        self.gridID = gridID 
        
        name = disc + "__" + str(gridID)

        print "building roc curve histogram"
        
        self.sig_hist = rt.TH1F(name, name, N_BEST_ROC_POINTS, MIN_BEST, MAX_BEST)        

        self.x = []
        self.y = []
        
    def add_point(self, thres, signal_eff, bkg_rej):
        self.x.append(signal_eff)
        self.y.append(bkg_rej)
        self.eff_dict[thres] = (signal_eff, bkg_rej)

        #fill the binned ROC curve with best of points
        axis = self.sig_hist.GetXaxis()
        bin = axis.FindBin(signal_eff)
        
        bincontent = self.sig_hist.GetBinContent(bin)

        if bkg_rej > bincontent:
            self.sig_hist.SetBinContent(bin, bkg_rej)
        
    def build_tgraph(self):
        print "building tgraph" 

class jet_collection:

    def __init__(self, filename):
        self.filename = filename
        self.jetlist = []
        self.disc_list = []
        self.njets = 0

        self.build_jets()

    def get_variable_array(self, var):
        temp = []

        for jet in self.jetlist:
            temp.append(jet.jetvars[var])

        return temp

    #indexed by (jet1, jet2, ...) with jet1 = (var1, var2, var3...)
    def get_jet_var_index(self, tuple_of_weight_ranges):
        names_to_fill = []
        result = []
        nvars = len(tuple_of_weight_ranges)

        #fill in the blanks
        for ii in range(self.njets):
            result.append([])
            for jj in range(nvars):
                result[ii].append([])

        for rr in tuple_of_weight_ranges:
            names_to_fill.append(rr.name)            

        for jj in range(len(self.jetlist)):
            thisJet = self.jetlist[jj]
            for name in names_to_fill:
                index = names_to_fill.index(name)                
                result[jj][index] = thisJet.jetvars[name]

        print result

        return result

    def build_grid_tree(self, tree_name):
        
        tree = rt.TTree(tree_name, tree_name)    
    
            # Create a struct for the run information
        line_to_process = "struct MyStruct{ Int_t id;Int_t caloGenMatch;Float_t pt;Float_t eta;Float_t phi;Float_t ipsiglog;Float_t ipelog;Float_t ipmed; Int_t nGrid; Int_t evNum; Int_t svntrack; Float_t svlxy; Float_t svlxysig; Float_t svmass; Int_t event; Int_t lumi; Int_t run;"

        for disc in self.disc_list:
            line_to_process += " Float_t %s[1000];" % disc
            line_to_process += " Int_t %sID[1000];" % disc
        
        line_to_process += "};"

        print line_to_process

        rt.gROOT.ProcessLine(line_to_process)
        from ROOT import MyStruct        
        s = MyStruct()
        
        # identification
        tree.Branch("jetid",rt.AddressOf(s,"id"),"jetid/I")
        tree.Branch("nGrid",rt.AddressOf(s,"nGrid"),"nGrid/I")
        tree.Branch("evNum",rt.AddressOf(s,"evNum"),"evNum/I")

        tree.Branch("event",rt.AddressOf(s,"event"),"event/I")
        tree.Branch("run",rt.AddressOf(s,"run"),"run/I")
        tree.Branch("lumi",rt.AddressOf(s,"lumi"),"lumi/I")

        # kinematics
        tree.Branch("pt",rt.AddressOf(s,"pt"),"pt/F")
        tree.Branch("eta",rt.AddressOf(s,"eta"),"eta/F")
        tree.Branch("phi",rt.AddressOf(s,"phi"),"phi/F")

        # ip information for the metric
        tree.Branch("ipsiglog",rt.AddressOf(s,"ipsiglog"),"ipsiglog/F")
        tree.Branch("ipmed",rt.AddressOf(s,"ipmed"),"ipmed/F")
        tree.Branch("ipelog",rt.AddressOf(s,"ipelog"),"ipelog/F")
        tree.Branch("caloGenMatch",rt.AddressOf(s,"caloGenMatch"),"caloGenMatch/I")

        # sv information
        tree.Branch("svlxy",rt.AddressOf(s,"svlxy"),"svlxy/F")
        tree.Branch("svlxysig",rt.AddressOf(s,"svlxysig"),"svlxysig/F")
        tree.Branch("svmass",rt.AddressOf(s,"svmass"),"svmass/F")
        tree.Branch("svntrack",rt.AddressOf(s,"svntrack"),"svntrack/I")

        # add a branch for each discriminator
        for disc in self.disc_list:
            tree.Branch(disc, eval("s.%s" % disc) ,"%s[nGrid]/F" % disc )
            tree.Branch(disc+"ID", eval("s.%sID" % disc) ,"%sID[nGrid]/I" % disc )

        # fill  the information common to every discriminant
        s.id = s.caloGenMatch = s.pt = s.eta = s.phi = 0
        s.ipsiglog = s.ipmed = s.ipelog = 0
        s.run = s.ls = s.event = -1

        #fill the tree in terms of jets
        for jet in self.jetlist:
            s.ipsiglog = jet.jetvars["jetIPSigLogSum2D"] 
            s.ipmed = jet.jetvars["jetMedianIPSig2D"] 
            s.ipelog = jet.jetvars["jetELogIPSig2D"] 
            s.pt = jet.jetvars["caloJetPt"]
            s.eta = jet.jetvars["caloJetEta"]
            s.phi = jet.jetvars["caloJetPhi"]
            s.id = jet.jetvars["jetID"] 
            s.caloGenMatch = jet.jetvars["caloGenMatch"]
            s.evNum = jet.jetvars["evNum"]
            s.svlxysig = jet.jetvars["jetSvLxySig"] 
            s.svlxy = jet.jetvars["jetSvLxy"] 
            s.svntrack = jet.jetvars["jetSvNTrack"]
            s.svmass = jet.jetvars["jetSvMass"] 

            s.run = jet.jetvars["run"]
            s.event = jet.jetvars["event"]
            s.lumi = jet.jetvars["lumi"]

            for disc in self.disc_list:            
                jet_dict = jet.disc_dict                
                s.nGrid = len(jet_dict)
                
                disc_list = []
                discID_array = array.array("i", range(1,s.nGrid+1))

                #loop over every calculated disc and assign it a unique ID
                for key in jet_dict.keys():
                    (thisdisc, wvector) = key 
                    val = jet_dict[key] 
                    
                    if thisdisc == disc: disc_list.append(val)

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

            iev += 1
            entry = itlist.Next()
            tree.GetEntry(entry)                                                            

            for jj in range(tree.nCaloJets):
                jetvars = {}

                fillZero = False
                #jet must have 1 track
                if tree.jetNTracks[jj] == 0 or tree.sumTrackPt[jj] == 0: continue #fillZero = True
                
                self.njets+= 1
                # if fillZero:
                #     jetvars["jetMedianIPLogSig2D"] = tree.jetMedianIPLogSig2D[jj]
                #     jetvars["caloJetAlphaMax"]     = tree.caloJetAlphaMax[jj] / tree.sumTrackPt[jj]
                #     jetvars["jetMedianJetDist"]    = tree.jetMedianJetDist[jj]

                jetID = tree.jetID[jj]
#                print "jetNTracks", tree.jetNTracks[jj], "sumtrackpt", tree.sumTrackPt[jj], "jetid", jetID
                jetvars["jetIVFLxySig"]        = tree.jetIVFLxySig[jj]
                jetvars["jetIVFNTrack"]        = tree.jetIVFNTrack[jj]
                jetvars["jetMedianIPLogSig2D"] = tree.jetMedianIPLogSig2D[jj]
                jetvars["jetIPSigLogSum2D"]    = tree.jetIPSigLogSum2D[jj]                    
                jetvars["caloJetAlphaMax"]     = tree.caloJetAlphaMax[jj] / tree.sumTrackPt[jj]
                jetvars["jetMedianJetDist"]    = tree.jetMedianJetDist[jj]
                jetvars["jetMedianIPSig2D"]    = tree.jetMedianIPSig2D[jj]
                jetvars["jetELogIPSig2D"]      = tree.jetELogIPSig2D[jj]
                jetvars["caloJetPt"]           = tree.caloJetPt[jj]
                jetvars["caloJetEta"]          = tree.caloJetEta[jj]
                jetvars["caloJetPhi"]          = tree.caloJetPhi[jj]
                jetvars["jetID"]               = tree.jetID[jj]
                jetvars["caloGenMatch"]        = tree.caloGenMatch[jj]
                jetvars["evNum"]               = tree.evNum
                jetvars["jetSvLxySig"]         = tree.jetSvLxySig[jj]
                jetvars["jetSvLxy"]            = tree.jetSvLxy[jj]
                jetvars["jetSvNTrack"]         = tree.jetSvNTrack[jj]
                jetvars["jetSvMass"]           = tree.jetSvMass[jj]
                jetvars["event"]               = tree.event
                jetvars["run"]                 = tree.run
                jetvars["lumi"]                = tree.lumi

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
    def calc_metric_vtx(self, pars):
        (w2, w3, w5) = pars #, w4, w5) = pars

        x1 = self.jetvars["jetELogIPSig2D"]
        x2 = self.jetvars["jetMedianIPSig2D"]
        x3 = self.jetvars["jetIPSigLogSum2D"]
        x4 = self.jetvars["jetSvMass"]
        x5 = self.jetvars["jetSvLxySig"]

        return  w2*(x2) + w3*(x3) + w5*x5

    def calc_metric_simple(self, pars):
        (w2, w3) = pars

        x2 = self.jetvars["jetMedianIPSig2D"]
        x3 = self.jetvars["jetIPSigLogSum2D"]

        val =   w2*(x2) + w3*(x3)

        return val



    def calc_metric_ivf(self, pars):
        #(w1, w2, w3) = pars
        (w1, w3) = pars

        x1 = self.jetvars["jetMedianIPLogSig2D"]
        x2 = self.jetvars["caloJetAlphaMax"] 
        x3 = self.jetvars["jetMedianJetDist"]

        val = w1*x1 + w2*x2 + w3*x3

        return val


    def calc_metric_alph(self, pars):
        (w1, w2, w3) = pars

        # x1 = self.jetvars["jetMedianIPLogSig2D"]
        # # x2 = self.jetvars["jetIVFLxySig"]
        # x3 = self.jetvars["jetIVFNTrack"]

        x1 = self.jetvars["jetMedianIPLogSig2D"]
        x2 = self.jetvars["caloJetAlphaMax"] 
        x3 = self.jetvars["jetMedianJetDist"]

        val = w1*x1 + w2*x2 + w3*x3

        return val

    def fill_disc(self, disc, weight_space):        
        wgrid = weight_space.grid
        #go over each point in the weight space            
        for wvector in wgrid:            
            #self.disc_dict[disc, wvector] = self.calc_metric_val1(wvector)
            if disc == "simple":
                self.disc_dict[disc, wvector] = self.calc_metric_simple(wvector)
            elif disc == "vtx":
                self.disc_dict[disc, wvector] = self.calc_metric_vtx(wvector)
            elif disc == "ivf":
                self.disc_dict[disc, wvector] = self.calc_metric_ivf(wvector)
            elif disc == "alph":
                self.disc_dict[disc, wvector] = self.calc_metric_alph(wvector)

            else:
                print "FALSE DISCRIMINANT: ", disc, "---EXITING----"
                exit(1)
            

output_file = rt.TFile(options.output, "RECREATE")
        
ana = analysis(options.sig_file, options.bkg_file)

ana.build_jetcollections()

######################SIMPLE DISCRIMINANT###################

#build a dummy range for the simple discriminant
#med_range        = weight_range("jetMedianIPSig2D", 0, 1, 2, 1, explicit=[0])
#ipsiglog_range   = weight_range("jetIPSigLogSum2D", 0, 1, 2, 1, explicit=[0])

mediplog_range   = weight_range("jetMedianIPLogSig2D", 0, 1, 2, 1, explicit=[0])
ivflxysig_range  = weight_range("jetIVFLxySig", 0, 1, 2, 1, explicit=[0])
alpha_range    = weight_range("caloJetAlphaMax", 0, 1, 2, 1,explicit=[0])
dist_range    = weight_range("jetMedianJetDist", 0, 1, 2, 1, explicit=[0])
                             
#ivfntracks_range = weight_range("jetIVFNTrack", 0, 1, 2, 1, explicit=[0])

#disc_simple_tuple = (med_range, ipsiglog_range)
#disc_simple_tuple = (mediplog_range, ivflxysig_range, ivfntracks_range ) 
#disc_simple_tuple = (mediplog_range, ivfntracks_range ) 
disc_simple_tuple = (mediplog_range, alpha_range, dist_range)

#calculate the fischer weights for this range
print "-- Building Weights --"
fweights = ana.get_fischer_weights(disc_simple_tuple)
print "Fischer Weights: ", fweights 

#build the new fischer ranges with the weights
#fischer_med_range = weight_range("jetMedianIPSig2D", 0, 1, 2, 1, explicit=[fweights[0]])
#fischer_ipsiglog_range = weight_range("jetIPSigLogSum2D", 0, 1, 2, 1, explicit=[fweights[1]])
#fischer_tuple = (fischer_med_range, fischer_ipsiglog_range)

fischer_mediplog_range    = weight_range("jetMedianIPLogSig2D", 0, 1, 2, 1, [-1 * fweights[0]])
fischer_alpha_range    = weight_range("caloJetAlphaMax", 0, 1, 2, 1, [-1 * fweights[0]])
fischer_dist_range    = weight_range("jetMedianJetDist", 0, 1, 2, 1, [-1 * fweights[0]])
#fischer_ivflxysig_range   = weight_range("jetIVFLxySig", 0, 1, 2, 1, [-1 * fweights[1]])
#fischer_ivfntracks_range = weight_range("jetIVFNTrack", 0, 1, 2, 1, [fweights[1]])
#fischer_tuple            = (fischer_mediplog_range, fischer_ivflxysig_range, fischer_ivfntracks_range)
#fischer_tuple            = (fischer_mediplog_range, fischer_ivfntracks_range)
fischer_tuple             = (fischer_mediplog_range, fischer_alpha_range, fischer_dist_range)#fischer_ivflxysig_range)

#build the weight space to fill the discriminant
disc_simple_weight_space = weight_space("alph", 3, fischer_tuple)
ana.fill_discriminant("alph", disc_simple_weight_space)

#disc_simple_weight_space.write_output(options.output_GID)

######################SIMPLE DISCRIMINANT END ###################

######################VTX DISCRIMINANT ###################
# med_range = weight_range("jetMedianIPSig2D", 0, 1, 2, 1, explicit=[0])
# ipsiglog_range = weight_range("jetIPSigLogSum2D", 0, 1, 2, 1, explicit=[0])
# svlxysig_range = weight_range("jetSvLxySig", 0, 1, 2, 1, explicit=[0])
# disc_vtx_tuple = (med_range, ipsiglog_range, svlxysig_range)

# vtx_fweights = ana.get_fischer_weights(disc_vtx_tuple)

# fischer_med_range = weight_range("jetMedianIPSig2D", 0, 1, 2, 1, explicit=[vtx_fweights[0]])
# fischer_ipsiglog_range = weight_range("jetIPSigLogSum2D", 0, 1, 2, 1, explicit=[vtx_fweights[1]])
# fischer_svlxysig_range = weight_range("jetSvLxySig", 0, 1, 2, 1, explicit=[vtx_fweights[2]])
# fischer_vtx_tuple  = (fischer_med_range, fischer_ipsiglog_range, fischer_svlxysig_range)

# disc_vtx_weight_space =  weight_space("vtx", 3, fischer_vtx_tuple)
# ana.fill_discriminant("vtx", disc_vtx_weight_space)

######################END DISCRIMINANT ###################

(sig_tree, bkg_tree) = ana.get_trees() 

output_file.cd()
sig_tree.Write()
bkg_tree.Write()

if options.do_roc:    
    ana.build_roc_curves(sig_tree, bkg_tree, -.0001, .0005 , 1000)
    tgraphs = ana.generate_tgraphs() 
#    roc_canvas = ana.build_roc_canvas()
#    roc_canvas.Write()

    print "-- Writing TGraphs--"
    
    for tgraph in tgraphs: 
        print tgraph
        tgraph.Write()    

output_file.Close()
