
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
    varbins = [300, 325, 350, 375, 400, 425, 450, 600]
#    varbins = [400, 425, 450, 500, 525, 550, 575, 600, 625, 650, 1000]

    #varbins = makebins(options.xmin, options.xmax, float(options.var_bin), .5)
    print varbins

class stack:
    def __init__(self,  fillColor, fillStyle, lineStyle, lineWidth, label):
        self.fillColor   = fillColor
        self.fillStyle   = fillStyle
        self.lineStyle   = lineStyle
        self.lineWidth   = lineWidth
        self.lineWidth   = lineWidth
        self.label       = label

        #list of samples and the summed histgraphs
        self.samples     = []        
        self.hist_den    = None
        self.hist_num    = None
        self.ratio_graph = None
        
class sample:
    def __init__(self, file_name, cut, tree_name, isSignal, isData, xsec, fillColor, fillStyle, lineStyle, lineWidth, stack, label):

        # set configuration
        self.file_name = file_name 
        self.cut = cut
        self.tree_name = options.tree 
        self.xsec = eval("float(%s)"%xsec) 
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
        
    def init_graph_style(self, graph):
        
        if int(self.fillStyle) != 0:
            eval("graph.SetFillColor(rt.%s)" % self.fillColor)
            eval("graph.SetFillColor(rt.%s)" % self.fillColor)

        eval("graph.SetFillStyle(%s)" % self.fillStyle)
        eval("graph.SetLineColor(rt.%s)" % self.fillColor)
        eval("graph.SetLineStyle(%s)" % self.lineStyle)
        graph.SetLineWidth(int(self.lineWidth))        

        eval("graph.SetMarkerColor(rt.%s)" % self.fillColor)

        if self.isData:
            graph.SetMarkerStyle(8)
            graph.SetMarkerSize(1.1)
        else:
            graph.SetMarkerStyle(4)
            graph.SetMarkerSize(1.1)


config_file = open(options.file, "r")
lines = map(lambda x: x.rstrip("\n"), config_file.readlines())

#dictionary from filename to sample class
samples         = {} 
#dictionary from stack string to array of samples corresponding to stack
stacks          = {}
stack_hists_den = {}
stack_hists_num = {}
global_key_order = []

#ignore the first line for labels
for line in lines[1:]:
    if  "//" in line: continue
    (file_name, cut, tree_name, isSignal, isData,  xSec, fillColor, fillStyle, lineStyle, lineWidth, stack_name, label) = line.split("|")
    samples[file_name] = sample(*line.split("|"))

    #check if the stack already exists
    if stack_name not in stacks.keys():
        global_key_order.append(stack_name)
        stacks[stack_name] = stack(fillColor, fillStyle, lineStyle, lineWidth, stack_name)
        stacks[stack_name].samples = [samples[file_name]]        
        stacks[stack_name].hist_den = samples[file_name].hist_den.Clone() 
        stacks[stack_name].hist_num = samples[file_name].hist_num.Clone() 
    else:
        stacks[stack_name].samples.append(samples[file_name])

print stacks

stack_ratios = []

#for stacked objects (assume only one object per stack for now)
for key in stacks.keys(): 
    #build the stack ratio

    for samp in stacks[key].samples:
        thisFile = rt.TFile(samp.file_name)
        thisTree = thisFile.Get(samp.tree_name)        
        output.cd()
        
        thisCut = options.cut
        denCut =  "((%s) && (%s) && (%s))" % ( samp.cut, options.cut, options.dencut) 
        numCut =  "((%s) && (%s) && (%s))" % ( samp.cut, options.cut, options.numcut) 

        print "Sample: ", samp.file_name,  "den Cut: ", denCut, " num Cut", numCut

        #fill each histogram 
        draw_string_den = "%s>>%s" % (options.var, samp.hist_name_den)        
        draw_string_num = "%s>>%s" % (options.var, samp.hist_name_num)        
        print draw_string_den
        print draw_string_num
        samp.total_events = thisTree.GetEntries() 

        nevents_den = thisTree.Draw(draw_string_den , denCut)
        nevents_num = thisTree.Draw(draw_string_num , numCut)

        #restrict the weight to the number of events passing the triggers
        total_events = thisTree.GetEntries("(%s)" % (denCut))

        print "contributing numerator", nevents_num, " jet weight: ", float(samp.xsec)#,   float(samp.xsec) / float(total_events)

        #nevents = thisTree.GetEntries(thisCut)

        samp.hist_den.Sumw2()
        samp.hist_num.Sumw2()

        if not samp.isSignal and not samp.isData:
            samp.hist_den.Scale(float(samp.xsec / float(total_events)))
            samp.hist_num.Scale(float(samp.xsec / float(total_events)))

        #samp.ratio_graph = rt.TGraphAsymmErrors()
        #samp.ratio_graph.BayesDivide(samp.hist_num, samp.hist_den)

        hist_den = stacks[key].hist_den

        if(hist_den.Integral() > 0):                        
            stacks[key].hist_den.Add(samp.hist_den)
            stacks[key].hist_num.Add(samp.hist_num)
        else:
            stacks[key].hist_den = samp.hist_den.Clone()
            stacks[key].hist_num = samp.hist_num.Clone()

    stacks[key].ratio_graph = rt.TGraphAsymmErrors()
    stacks[key].hist_num.Print()
    stacks[key].hist_den.Print()
    stacks[key].ratio_graph.BayesDivide(stacks[key].hist_num, stacks[key].hist_den)

    print "PRINT? -----"
    stacks[key].ratio_graph.Print()
    #use the first sample in the stack to set style 
    stacks[key].samples[0].init_graph_style(stacks[key].ratio_graph)

    
output.cd()
#build the canvas
canvas = rt.TCanvas("plot","plot", 1024, 768)
canvas.cd()

first_key = stacks.keys()[0]
print "drawing first: ", first_key
draw_first = stacks[first_key].ratio_graph
draw_first.Print()
draw_first.Draw("ap")


if len(stacks.keys()) > 1: 
    for key in stacks.keys()[1:]: stacks[key].ratio_graph.Draw("psame")

if options.islog: canvas.SetLogy()

draw_first.GetYaxis().SetRangeUser(options.ymin,options.ymax)
draw_first.GetXaxis().SetTitle(options.xlabel)
draw_first.GetYaxis().SetTitle(options.ylabel)

leg = rt.TLegend(.4, .45, .75, .75) #canvas.BuildLegend()
leg = rt.TLegend(.73, .41, .92, .83) 


for key in global_key_order:
    stack = stacks[key]
    leg.AddEntry(stack.ratio_graph, stack.label, 'pl')

# for key in stacks.keys(): 
#     first = stacks[key][0]
#     if first.isData:
#         leg.AddEntry(first.ratio_graph, first.stack, 'pl')
#     else:
#         leg.AddEntry(first.ratio_graph, first.stack, 'pl')

leg.SetFillColor(0)
leg.SetLineColor(0)
leg.Draw("same")
CMS_lumi.CMS_lumi(canvas, 4, 0)



raw_input("RAW INPUT")

plot_name = "%s_%s" % ( options.file[:-5], options.var)
#canvas.SaveAs("%s.C" % plot_name)
#canvas.Print("%s.pdf" % plot_name)

