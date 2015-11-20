
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

parser.add_option("--refstack", dest="refstack",
		                    help="stack to take ratios against",
		                    action="store",type="string",default="DATA")


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

class stack:
    def __init__(self,  fillColor, fillStyle, lineStyle, lineWidth, label, isSig, isData):
        self.fillColor = fillColor
        self.fillStyle = fillStyle
        self.lineStyle = lineStyle
        self.lineWidth = lineWidth
        self.lineWidth = lineWidth
        self.label     = label
        self.isSig     = isSig
        self.isData    = isData

        #list of samples and the summed histgraphs
        self.samples     = []  
        self.hist    = None      
        self.ratio_hist = None
        
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
        #self.hist_name_den = file_name.split(".")[0] + "_" + tree_name + "_den"
        #self.hist_name_num = file_name.split(".")[0] + "_" + tree_name + "_num"
        self.ratio_hist = None
        self.isSignal = int(isSignal)
        self.isData =  int(isData)
        self.total_events     = -1

        #build the corresponding histogram
        self.hist = None

        if options.var_bin == 0:
            self.hist = rt.TH1F(self.hist_name, self.label, nbins, xmin, xmax)
        else:
            self.hist = rt.TH1F(self.hist_name, self.label, len(varbins)-1, array.array("d",varbins))
        
        if int(self.fillStyle) != 0:
            eval("self.histSetFillColor(rt.%s)" % self.fillColor)
            eval("self.histSetFillColor(rt.%s)" % self.fillColor)

        eval("self.hist.SetFillStyle(%s)" % self.fillStyle)
        eval("self.hist.SetLineColor(rt.%s)" % self.fillColor)
        eval("self.hist.SetLineStyle(%s)" % self.lineStyle)
        self.hist.SetLineWidth(int(self.lineWidth))        

        eval("self.hist.SetMarkerColor(rt.%s)" % self.fillColor)

        if self.isData:
            self.hist.SetMarkerStyle(8)
            self.hist.SetMarkerSize(1.3)
        elif self.isSignal:
            self.hist.SetMarkerStyle(1)
            self.hist.SetMarkerSize(1)
        else:
            self.hist.SetFillStyle(3001)
            eval('self.hist.SetFillColorAlpha(rt.%s, 0.20)' % self.fillColor)
            self.hist.SetMarkerStyle(25)
            self.hist.SetMarkerSize(1.3)            
        
    # def init_hist_style(self, hist):
        
    #     if int(self.fillStyle) != 0:
    #         eval("hist.SetFillColor(rt.%s)" % self.fillColor)
    #         eval("hist.SetFillColor(rt.%s)" % self.fillColor)

    #     eval("hist.SetFillStyle(%s)" % self.fillStyle)
    #     eval("hist.SetLineColor(rt.%s)" % self.fillColor)
    #     eval("hist.SetLineStyle(%s)" % self.lineStyle)
    #     hist.SetLineWidth(int(self.lineWidth))        

    #     eval("hist.SetMarkerColor(rt.%s)" % self.fillColor)

    #     if self.isData:
    #         hist.SetMarkerStyle(8)
    #         hist.SetMarkerSize(2)
    #     else:
    #         hist.SetMarkerStyle(4)
    #         hist.SetMarkerSize(2)

config_file = open(options.file, "r")
lines = map(lambda x: x.rstrip("\n"), config_file.readlines())

#dictionary from filename to sample class
samples         = {} 
#dictionary from stack string to array of samples corresponding to stack
stacks          = {}
global_key_order = []

#ignore the first line for labels
for line in lines[1:]:
    if  "//" in line: continue
    (file_name, cut, tree_name, isSignal, isData,  xSec, fillColor, fillStyle, lineStyle, lineWidth, stack_name, label) = line.split("|")
    samples[file_name] = sample(*line.split("|"))

    #check if the stack already exists
    if stack_name not in stacks.keys():
        print "creating stack", label, " is signal? ", isSignal, " is data? ", isData
        global_key_order.append(stack_name)
        stacks[stack_name]         = stack(fillColor, fillStyle, lineStyle, lineWidth, stack_name, int(isSignal), int(isData))
        stacks[stack_name].samples = [samples[file_name]]        
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
        
        #        thisCut = options.cut
        thisCut =  "((%s) && (%s))" % (samp.cut, options.cut) 

        print "Sample: ", samp.file_name,  " Cut: ", thisCut

        #fill each histogram 
        draw_string = "%s>>%s" % (options.var, samp.hist_name)        
        print draw_string

        nevents = thisTree.Draw(draw_string , thisCut)

        samp.total_events = thisTree.GetEntries(options.cut) 

        #restrict the weight to the number of events passing the triggers
        total_events = thisTree.GetEntries(options.cut)

        print "contributing to hist", samp.hist

        #nevents = thisTree.GetEntries(thisCut)
        samp.hist.Sumw2()

        if not samp.isSignal and not samp.isData:
            samp.hist.Scale(float(samp.xsec))

        if(stacks[key].hist != None):                        
            stacks[key].hist.Add(samp.hist)
        else:
            stacks[key].hist = samp.hist.Clone()
            
#parse the ref stack histogram to take all the ratios relative to
if options.refstack not in stacks.keys():
    print "[ERROR] reference stack does not exist:", options.refstack
    print "available stacks", stacks.keys()
    exit(1)


if options.norm1 and stacks[options.refstack].hist.Integral() != 0:
    stacks[options.refstack].hist.Scale(1. / float(stacks[options.refstack].hist.Integral()))


ref_hist = stacks[options.refstack].hist.Clone()

#loop over the stacks and build the ratio relative to the reference sample
for key in stacks.keys(): 
    if options.norm1 and stacks[key].hist.Integral() != 0:
        weight = 1. / float(stacks[key].hist.Integral())
        stacks[key].hist.Scale(weight)
        print "normalized histogram", stacks[key].label,  "by factor", weight

    #build the difference and ratio from the reference hist
    stacks[key].ratio_hist = stacks[key].hist.Clone()
    # diff_hist = stack_hist.Clone()
    # diff_hist.Add(ref_hist, -1.0)
    stacks[key].ratio_hist.Add(ref_hist, -1)
    stacks[key].ratio_hist.Divide(stacks[key].hist)       
    
output.cd()
#build the canvas
canvas     = rt.TCanvas("plot","plot", 1150, 1300)
top_pad    = rt.TPad("pad1", "The pad 80% of the height",0.0, 0.4, 1.0, 1.0, 21)
top_pad.SetBottomMargin(0)
top_pad.SetBottomMargin(0)
top_pad.SetRightMargin(.05)
top_pad.SetLeftMargin(.2)
top_pad.SetTopMargin(.11)
top_pad.SetFillColor(0)
bottom_pad = rt.TPad("pad2", "The pad 20% of the height",0.0, 0.0, 1.0, 0.4, 22)
bottom_pad.SetTopMargin(0)
bottom_pad.SetRightMargin(0.05)
bottom_pad.SetBottomMargin(.4)
bottom_pad.SetLeftMargin(.2)
bottom_pad.SetFillColor(0)

#draw the two pads
top_pad.Draw()
bottom_pad.Draw()

draw_first      = None
draw_rest       = []
draw_rest_names = []
max             = -1

print "keys in stack", 
for ii in stacks.keys(): print ii, 
print "\n"

for key in stacks.keys():     
    #draw each histogram
    stack_hist = stacks[key].hist
    thisMax = stack_hist.GetMaximum() 
    if thisMax > max: 
        if draw_first != None and stacks[key].label not in draw_rest_names:
            print "shifting down", stacks[key].label
            draw_rest.append(draw_first)
            draw_rest_names.append(draw_first.label)

        max = thisMax
        draw_first = stacks[key]
    else:
        if stacks[key] not in draw_rest:
            print "appending", stacks[key].label
            draw_rest_names.append(stacks[key].label)
            draw_rest.append(stacks[key])

print "DRAW FIRST:", draw_first.label

print "Draw Rest", 
for ii in draw_rest:
    print ii.label, 


top_pad.cd()
rootlogon.style()
# rt.gStyle.SetPadBottomMargin(0)
# rt.gStyle.SetLabelSize(0.0,"x")
# rt.gStyle.SetOptStat(0)
# rt.gStyle.SetOptTitle(0)


# if draw_first.hist.Integral() != 0 and options.norm1:
#     print "DRAW FIRST NOT 0? Rescaling...", 
#     draw_first.hist.Scale(1. / float(draw_first.hist.Integral()))

#draw the first hist
if draw_first.isSig:
    draw_first.hist.Draw("hist")
elif draw_first.isData:
    draw_first.hist.Draw("p")
else:    
    draw_first.hist.Draw("hist e0 bar")


# draw the rest
for ii in draw_rest: 
    if ii.hist.Integral() != 1 and options.norm1:
        print "WTF", ii.hist.Integral()
    if ii.isSig:        
        ii.hist.Draw("hist same")
    elif ii.isData:
        ii.hist.Draw("p same")
    else:
        ii.hist.Draw("hist e0 bar same")


#draw_first.hist.GetYaxis().SetRangeUser(options.ymin,options.ymax)
draw_first.hist.GetXaxis().SetTitle("")
#draw_first.hist.GetXaxis().SetFontSize(0)
draw_first.hist.GetYaxis().SetLabelSize(.06)
draw_first.hist.GetYaxis().SetTitleSize(.1)
draw_first.hist.GetYaxis().SetTitleOffset(.8)
draw_first.hist.GetYaxis().SetTitle(options.ylabel)

leg = rt.TLegend(.73, .41, .92, .83) 

for key in global_key_order:
    stack = stacks[key]
    leg.AddEntry(stack.hist, stack.label, 'pl')

# for key in stacks.keys(): 
#     first = stacks[key][0]
#     if first.isData:
#         leg.AddEntry(first.ratio_graph, first.stack, 'pl')
#     else:
#         leg.AddEntry(first.ratio_graph, first.stack, 'pl')

leg.SetFillColor(0)
leg.SetLineColor(0)
leg.Draw("same")
CMS_lumi.CMS_lumi(top_pad, 4, 0)

#draw the ratio pad
bottom_pad.cd()
rootlogon.style()

# rt.gStyle.SetPadTopMargin(0)
# rt.gStyle.SetLabelSize(0.12,"xy")
# rt.gStyle.SetOptStat(0)
# rt.gStyle.SetOptTitle(0)

draw_order = []

for key in stacks.keys():
    print "checking stack", stacks[key].label, "...",
    if key == options.refstack: 
        print "is reference stack"
        continue
    
    if stacks[key].isSig: 
        print "is signal stack"
        continue

    print "adding stack...."
    draw_order.append(stacks[key])

print "draw order for ratio plot",
for stack in draw_order: print stack.label,
print "\n"


if len(draw_order) != 0:

    draw_first = draw_order[0].ratio_hist

    if len(draw_order) == 1:
        draw_first.SetLineColor(rt.kBlack)

    draw_first.Draw("e1")
    # y axis
    draw_first.GetYaxis().SetRangeUser(-1, 1)
    draw_first.GetYaxis().SetTitle("(MC-DATA)/MC")
    draw_first.GetYaxis().SetTitleSize(0.1)
    draw_first.GetYaxis().SetTitleOffset(.8)
    draw_first.GetYaxis().SetLabelSize(0.07)

    # x axis
    draw_first.GetXaxis().SetLabelSize(0.10)
    draw_first.GetXaxis().SetTitleSize(0.11)
    draw_first.GetXaxis().SetTitleOffset(1.15)
    draw_first.GetXaxis().SetTitle(options.xlabel)

    if len(draw_order) > 1:
        for stack in draw_order[1:]:
            stack.ratio_hist.Draw("e1 same")

bottom_pad.Update()
top_pad.Update()
raw_input("RAW INPUT")

plot_name = "%s_%s" % ( options.file[:-5], options.var)
#canvas.SaveAs("%s.C" % plot_name)
#canvas.Print("%s.pdf" % plot_name)

