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
            graph.SetLineColor(color)
            graph.SetFillColor(0)
            graph.SetLineWidth(2)

            graph.Draw("same")
            color+=1
            if color > 100:
                color = 0

        #leg = self.canvas.BuildLegend()
        #leg.SetFillColor(0)
        
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

    def build_roc_curves(self, sig_tree, bkg_tree, xmin, xmax, npoints, explicit_array=None):
        print "-- Building ROC curves" 

        scan_points = list(np.linspace(xmin, xmax, npoints))

        if explicit_array !=  None:
            scan_points = explicit_array 

        #roc curve for each discriminant
        for disc in self.disc_list:

            #total number of jets
            nsig_tot = sig_tree.GetEntries("genMatch > 0 && %sID == 1" % disc)
            nbkg_tot = bkg_tree.GetEntries(" %sID == 1" % disc)        

            # one curve per grid ID
            for gid in range(1, self.nGID+1):
                if gid % 25 == 0:
                    print "-- %s, GID =  %i of %s" % (disc, gid, self.nGID)
                
                curve = roc_curve(disc, gid)

                # one point per threshold
                for threshold in scan_points:
                    sig_cut = "genMatch > 0 && %sID == %i && %s > %s" % (disc, gid, disc, threshold)
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
    
        return w
    
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
        

    # def build_arrays(self):
    #     for key in self.eff_dict.keys():
    #         (x, y) = self.eff_dict[key] 
    #         self.x.append(x)
    #         self.y.append(y)

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

        print result[:4]

        return result

    def build_grid_tree(self, tree_name):
        
        tree = rt.TTree(tree_name, tree_name)    
    
            # Create a struct for the run information
        line_to_process = "struct MyStruct{ Int_t id;Int_t genmatch;Float_t pt;Float_t eta;Float_t phi;Float_t ipsiglog;Float_t ipelog;Float_t ipmed;Int_t genmatch; Int_t nGrid; Int_t evNum; Int_t svntrack; Float_t svlxy; Float_t svlxysig; Float_t svmass; Int_t event; Int_t lumi; Int_t run;"

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
        tree.Branch("genMatch",rt.AddressOf(s,"genmatch"),"genMatch/I")

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
        s.id = s.genmatch = s.pt = s.eta = s.phi = 0
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
            s.genmatch = jet.jetvars["genMatch"]
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
            
            self.njets += tree.nCaloJets 

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
                jetvars["evNum"] = tree.evNum
                jetvars["jetSvLxySig"] = tree.jetSvLxySig[jj]
                jetvars["jetSvLxy"] = tree.jetSvLxy[jj]
                jetvars["jetSvNTrack"] = tree.jetSvNTrack[jj]
                jetvars["jetSvMass"] = tree.jetSvMass[jj]
                jetvars["event"] = tree.event
                jetvars["run"] = tree.run
                jetvars["lumi"] = tree.lumi

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
        (w1, w2, w3) = pars #, w4, w5) = pars

        x1 = self.jetvars["jetELogIPSig2D"]
        x2 = self.jetvars["jetMedianIPSig2D"]
        x3 = self.jetvars["jetIPSigLogSum2D"]
        x4 = self.jetvars["jetSvMass"]
        x5 = self.jetvars["jetSvLxySig"]

        return w1*(x1) + w2*(x2) + w3*(x3) #+ w4*x4 + w5*x5 #+ (self.jetvars["jetSvNTrack"] - 2)

    def fill_disc(self, disc, weight_space):        
        wgrid = weight_space.grid
        #go over each point in the weight space            
        for wvector in wgrid:            
            self.disc_dict[disc, wvector] = self.calc_metric_val1(wvector)

output_file = rt.TFile(options.output, "RECREATE")
        
ana = analysis(options.sig_file, options.bkg_file)

ana.build_jetcollections()

#build the space for the metric

#Fsicher Weights [  3.74524639e-10  -5.10701610e-05  -1.10238882e-05]
#Fischer Weights:  [  1.82221547e-10  -1.15562473e-04  -2.77406473e-05]
w1_range = weight_range("jetELogIPSig2D", 0, 1, 2, 1, explicit=[-3.74524639e-10])
w2_range = weight_range("jetMedianIPSig2D", 0, 1, 2, 1, explicit=[5.1070161e-05])
w3_range = weight_range("jetIPSigLogSum2D", 0, 1, 2, 1, explicit=[1.10238882e-05])
#w4_range = weight_range("jetSvMass", 0, 1, 2, 1, explicit=[2.16e-04])
#w5_range = weight_range("jetSvLxySig", 0, 1, 2, 1, explicit=[1.02e-06])

tuple_range = (w1_range, w2_range, w3_range) #, w4_range, w5_range)

metric_weight_space = weight_space("metric", 3, tuple_range)
ana.fill_discriminant("metric", metric_weight_space)

metric_weight_space.write_output(options.output_GID)

# w1_range2 = weight_range("elogipsig", 0, 1, 5, 4.)
# w2_range2 = weight_range("jetmedianipsig", 0, 1, 5, .05)
# w3_range2 = weight_range("ipsiglogsum", 0, 1, 5, 10.)
# metric2_weight_space = weight_space("metric2", 3, (w1_range2.vector, w2_range2.vector, w3_range2.vector))
# ana.fill_discriminant("metric2", metric2_weight_space)

print "-- Building Weights --"
fweights = ana.get_fischer_weights(tuple_range)
print "Fischer Weights: ", fweights 

(sig_tree, bkg_tree) = ana.get_trees() 

output_file.cd()
sig_tree.Write()
bkg_tree.Write()

#--xmin -.001 --xmax .005 
if options.do_roc:
    ana.build_roc_curves(sig_tree, bkg_tree, -.001, .005 , 1000)
    tgraphs = ana.generate_tgraphs() 
    output_file.cd()

    print "-- Building ROC Canvas -- "
    roc_canvas = ana.build_roc_canvas()
    output_file.cd()
    roc_canvas.Write()

    # print "-- Building Best ROC -- "
    # best_roc = ana.build_best_roc()
    # output_file.cd()
    # best_roc.Write()

    # print "-- Building Best ROC Cavnas -- "
    # best_roc_canvas = rt.TCanvas("best_roc_canvas","cavnas")
    # best_roc.SetLineWidth(3)
    # best_roc.Draw()
    # CMS_lumi.CMS_lumi(best_roc_canvas, 4, 0)            

    # output_file.cd()
    # best_roc_canvas.Write()
    
    #print "-- Writing tgraphs --"
    #for graph in tgraphs:
    #    graph.Write()

output_file.Close()
