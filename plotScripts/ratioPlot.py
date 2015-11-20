
import ROOT as rt
from  optparse  import OptionParser
import sys, os, array
import CMS_lumi 
import rootlogon
rootlogon.style()

rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptTitle(0)

#exponentially increasing bin sizes
def makebins(start_,end_,inc_,inc_inc_):
    bin = start_
    inc = inc_
    list = []
    while bin < end_:
        list.append(bin)
        bin+=inc
        inc*=(1+inc_inc_)

    list[-1] = end_
    return list

parser = OptionParser()


parser.add_option("-f","--file",dest="file",
		                    help="text file containing samples and configuration",
		                    action="store",type="string")

parser.add_option("-o", "--output", dest="output",
		                    help="output destination",
		                    action="store",type="string")

parser.add_option("-c", "--cut", dest="cut",
		                    help="cut to apply to variable",
		                    action="store",type="string")

parser.add_option("--numcut", dest="numcut",
		                    help="cut to apply to the numerator of the ratio",
		                    action="store",type="string",default="(1)")

parser.add_option( "--dencut", dest="dencut",
		                    help="cut to apply to the denominator of the ratio",
		                    action="store",type="string",default="(1)")

parser.add_option("-v", "--variable", dest="var",
		                    help="variable to plot",
		                    action="store",type="string")

parser.add_option("-t", "--tree", dest="tree",
		                    help="tree for variable to be plotted from",
		                    action="store",type="string",default="jets")

parser.add_option( "--varbins", dest="var_bin",
		                    help="variable binning. Configure within script",
		                    action="store", type="float",default=0)

parser.add_option( "--log", dest="islog",
		                    help="Y axis Log Scale",
		                    action="store_true", default=False)


parser.add_option( "--genmatch", dest="genmatch",
		                    help="Do generator matching for the signal MC",
		                    action="store_true", default=False)

parser.add_option( "--ivfvtxmatch", dest="ivfvtxmatch",
		                    help="Do generator vertex matching for the signal MC",
		                    action="store_true", default=False)

parser.add_option( "--svvtxmatch", dest="svvtxmatch",
		                    help="Do generator vertex matching for the signal MC",
		                    action="store_true", default=False)


parser.add_option( "--genmatchtrack", dest="genmatchtrack",
		                    help="Do generator matching for the signal MC track collections",
		                    action="store_true", default=False)


parser.add_option( "--xmin", dest="xmin",
		                    help="minimum x for variable",
		                    action="store",type="float",default=0)

parser.add_option( "--xmax", dest="xmax",
		                    help="maximum x for variable",
		                    action="store",type="float",default = 100)

parser.add_option( "--ymin", dest="ymin",
		                    help="minimum x for variable",
		                    action="store",type="float",default=0.001)

parser.add_option( "--ymax", dest="ymax",
		                    help="maximum x for variable",
		                    action="store",type="float",default = 1)



parser.add_option( "--xlabel", dest="xlabel",
		                    help="label for the x axis",
		                    action="store",type="string",default = "variable")

parser.add_option( "--ylabel", dest="ylabel",
		                    help="label for the y axis",
		                    action="store",type="string",default = "N Jets")

parser.add_option( "--nbins", dest="nbins",
		                    help="number of bins. Non-variable binning",
		                    action="store",type="int", default = 20)

parser.add_option("-l", "--lumi", dest="lumi",
		                    help="integrated luminosity for normalization",
		                    action="store",type="float", default = .20)

parser.add_option("--norm1", dest="norm1",
		                    help="The plotte histogram is normalized to 1",
		                    action="store_true", default = False)


(options, args) = parser.parse_args()

#output root file
output = rt.TFile(options.output, "RECREATE")

variable = options.var
cut = options.cut        
lumi = options.lumi
tree = options.tree
xmin = xmax = nbins = 0
varbins = None

if options.var_bin == 0:
    xmin = options.xmin
    xmax = options.xmax
    nbins = options.nbins
else:
    varbins = makebins(options.xmin, options.xmax, float(options.var_bin), .5)
    print varbins
class sample:
    def __init__(self, file_name, cut, tree_name, isSignal, isData, xsec, fillColor, fillStyle, lineStyle, lineWidth, stack, label):

        # set configuration
        self.file_name = file_name 
        self.cut = cut
        self.tree_name = options.tree 
        self.xsec = xsec 
        self.fillColor = fillColor
        self.fillStyle = fillStyle
        self.lineStyle = lineStyle
        self.lineWidth = lineWidth 
        self.stack = stack
        self.label = label
        self.hist_name = file_name.split(".")[0] + "_" + tree_name
        self.hist_name_den = file_name.split(".")[0] + "_" + tree_name + "_den"
        self.hist_name_num = file_name.split(".")[0] + "_" + tree_name + "_num"
        self.ratio_graph = None
        self.isSignal = int(isSignal)
        self.isData =  int(isData)
        self.total_events     = -1
        self.total_events_den = 1 
        self.total_events_num = 1 

        #build the corresponding histogram
        self.hist_den = None
        self.hist_num = None

        if options.var_bin == 0:
            self.hist_den = rt.TH1F(self.hist_name_den, self.label, nbins, xmin, xmax)
            self.hist_num = rt.TH1F(self.hist_name_num, self.label, nbins, xmin, xmax)
        else:
            self.hist_den = rt.TH1F(self.hist_name_den, self.label, len(varbins)-1, array.array("d",varbins))
            self.hist_num = rt.TH1F(self.hist_name_num, self.label, len(varbins)-1, array.array("d",varbins))
        
    def init_graph_style(self):
        
        if int(self.fillStyle) != 0:
            eval("self.ratio_graph.SetFillColor(rt.%s)" % self.fillColor)
            eval("self.ratio_graph.SetFillColor(rt.%s)" % self.fillColor)

        eval("self.ratio_graph.SetFillStyle(%s)" % self.fillStyle)
        eval("self.ratio_graph.SetLineColor(rt.%s)" % self.fillColor)
        eval("self.ratio_graph.SetLineStyle(%s)" % self.lineStyle)
        self.ratio_graph.SetLineWidth(int(self.lineWidth))        

config_file = open(options.file, "r")
lines = map(lambda x: x.rstrip("\n"), config_file.readlines())

#dictionary from filename to sample class
samples = {} 

#dictionary from stack string to array of samples corresponding to stack
stacks = {}

#ignore the first line for labels
for line in lines[1:]:
    if  "//" in line: continue
    (file_name, cut, tree_name, isSignal, isData,  xSec, fillColor, fillStyle, lineStyle, lineWidth, stack, label) = line.split("|")
    samples[file_name] = sample(*line.split("|"))

    #check if the stack already exists
    if stack not in stacks:
        stacks[stack] = [samples[file_name]]
    else:
        stacks[stack].append(samples[file_name])



#for stacked objects (assume only one object per stack for now)
for key in stacks.keys(): 
    for samp in stacks[key]:
        thisFile = rt.TFile(samp.file_name)
        thisTree = thisFile.Get(samp.tree_name)        
        output.cd()
        
        thisCut = options.cut
        denCut =  "(%s) && (%s) && (%s)" % (samp.cut, options.cut, options.dencut) 
        numCut =  "(%s) && (%s) && (%s)" % (samp.cut, options.cut, options.numcut) 
        #add in the gen amtching requirement to the cut
        # if samp.isSignal and options.genmatch:
        #     print samp.file_name, isSignal
        #     thisCut = "(" + options.cut + ") && ( caloGenMatch > 0 )" 
        # if samp.isSignal and options.genmatchtrack:
        #     print samp.file_name, isSignal
        #     thisCut = "(" + options.cut + ") && ( genMatchTrack > 0 )" 
        # if samp.isSignal and options.ivfvtxmatch:
        #     print samp.file_name, isSignal
        #     thisCut = "(" + options.cut + ") && ( jetIVFGenVertexMatched > 0 )" 
        # if samp.isSignal and options.svvtxmatch:
        #     print samp.file_name, isSignal
        #     thisCut = "(" + options.cut + ") && ( jetSvGenVertexMatched > 0 )" 

        print "Sample: ", samp.file_name,  "den Cut: ", denCut, " num Cut", numCut

        #fill each histogram 
        draw_string_den = "%s>>%s" % (options.var, samp.hist_name_den)        
        draw_string_num = "%s>>%s" % (options.var, samp.hist_name_num)        
        print draw_string_den
        print draw_string_num
        samp.total_events = thisTree.GetEntries() 

        nevents_den = thisTree.Draw(draw_string_den , denCut)
        nevents_num = thisTree.Draw(draw_string_num , numCut)

        #nevents = thisTree.GetEntries(thisCut)

        samp.hist_den.Sumw2()
        samp.hist_num.Sumw2()

        samp.ratio_graph = rt.TGraphAsymmErrors()
        samp.ratio_graph.BayesDivide(samp.hist_num, samp.hist_den)

        #set the style after theg graph is built
        samp.init_graph_style()

        if samp.isData:
            samp.ratio_graph.SetMarkerStyle(8)
        else:
            samp.ratio_graph.SetMarkerStyle(4)

            
        samp.ratio_graph.GetXaxis().SetTitle(options.xlabel)
        samp.ratio_graph.GetYaxis().SetTitle(options.ylabel)

        eval("samp.ratio_graph.SetMarkerColor(rt.%s)" % samp.fillColor)

output.cd()
#build the canvas
canvas = rt.TCanvas("plot","plot", 1024, 768)
canvas.cd()

draw_first = None
draw_rest = []
max = -99999999
for key in stacks.keys(): 
    for samp in stacks[key]:
        #draw each histogram
        thisMax = samp.ratio_graph.GetMaximum() 
        print "THIS IS THE MAX OF THE GRAPH: ", thisMax

        if thisMax > max: 
            if draw_first != None and samp.ratio_graph not in draw_rest:
                print "shifting down", samp.hist_name
                draw_rest.append(draw_first)

            max = thisMax
            draw_first = samp.ratio_graph
        else:
            if samp.ratio_graph not in draw_rest:
                print "appending", samp.hist_name
                draw_rest.append(samp.ratio_graph)


print "length of draw_rest", len(draw_rest)

print draw_first
draw_first.Draw("ap")


for ii in draw_rest: ii.Draw("psame")


if options.islog:
    canvas.SetLogy()

draw_first.GetYaxis().SetRangeUser(options.ymin,options.ymax)

leg = rt.TLegend(.4, .45 ,.75, .75)#canvas.BuildLegend()
for key in stacks.keys(): 
    for samp in stacks[key]:
        if samp.isData:
            leg.AddEntry(samp.ratio_graph, samp.label, 'pl')
        else:
            leg.AddEntry(samp.ratio_graph, samp.label, 'l')

leg.SetFillColor(0)
leg.SetLineColor(0)
leg.Draw("same")
CMS_lumi.CMS_lumi(canvas, 4, 0)


raw_input("RAW INPUT")

plot_name = "%s_%s" % ( options.file[:-5], options.var)
#canvas.SaveAs("%s.C" % plot_name)
#canvas.Print("%s.pdf" % plot_name)

