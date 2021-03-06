
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
		                    action="store",type="float", default = 1.0)

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
        self.isSignal = int(isSignal)
        self.isData =  int(isData)
        self.total_events = 1 

        #build the corresponding histogram
        self.hist = None
        if options.var_bin == 0:
            self.hist = rt.TH1F(self.hist_name, self.label, nbins, xmin, xmax)
        else:
            self.hist = rt.TH1F(self.hist_name, self.label, len(varbins)-1, array.array("d",varbins))

        if int(fillStyle) != 0:
            eval("self.hist.SetFillColor(rt.%s)" % fillColor)

        eval("self.hist.SetFillStyle(%s)" % fillStyle)
        eval("self.hist.SetLineColor(rt.%s)" % fillColor)
        eval("self.hist.SetLineStyle(%s)" % lineStyle)
        self.hist.SetLineWidth(int(lineWidth))        

    
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
        #add in the gen amtching requirement to the cut
        if samp.isSignal and options.genmatch:
            print samp.file_name, isSignal
            thisCut = "(" + options.cut + ") && ( caloGenMatch > 0 )" 
        if samp.isSignal and options.genmatchtrack:
            print samp.file_name, isSignal
            thisCut = "(" + options.cut + ") && ( genMatchTrack > 0 )" 
        if samp.isSignal and options.ivfvtxmatch:
            print samp.file_name, isSignal
            thisCut = "(" + options.cut + ") && ( jetIVFGenVertexMatched > 0 )" 
        if samp.isSignal and options.svvtxmatch:
            print samp.file_name, isSignal
            thisCut = "(" + options.cut + ") && ( jetSvGenVertexMatched > 0 )" 

        thisCut = thisCut + "&& (%s)" % samp.cut 

        print "Sample: ", samp.file_name,  "Cut: ", thisCut

        #fill each histogram 
        draw_string = "%s>>%s" % (options.var, samp.hist_name)
        print draw_string
        samp.total_events = thisTree.GetEntries() 
        nevents = thisTree.Draw(draw_string , thisCut)

        #nevents = thisTree.GetEntries(thisCut)

        if samp.isData:
            print "KEEPING ERRORS FOR: ", samp.hist
            samp.hist.Sumw2()
            samp.hist.SetMarkerStyle(8)
        else:
            samp.hist.SetMarkerStyle(4)

        if samp.hist.Integral() > 0: # and not samp.isData:
            scale_factor = float(options.lumi) * float(samp.xsec) / float(samp.total_events)

            if options.norm1:
                scale_factor = 1 / samp.hist.Integral()

            print samp, "SCALE FACTOR", scale_factor, "INTEGRAL", samp.hist.Integral()
            samp.hist.Scale(scale_factor) #float(samp.hist.Integral()))

            
        samp.hist.GetXaxis().SetTitle(options.xlabel)
        samp.hist.GetYaxis().SetTitle(options.ylabel)

        eval("samp.hist.SetMarkerColor(rt.%s)" % samp.fillColor)

output.cd()
#build the canvas
canvas = rt.TCanvas("plot","plot", 1024, 768)
canvas.cd()

draw_first = None
draw_rest = []
max = -1
for key in stacks.keys(): 
    for samp in stacks[key]:
        #draw each histogram
        thisMax = samp.hist.GetMaximum() 

        if thisMax > max: 

            if draw_first != None and samp.hist not in draw_rest:
                print "shifting down", samp.hist_name
                draw_rest.append(draw_first)

            max = thisMax
            draw_first = samp.hist
        else:
            if samp.hist not in draw_rest:
                print "appending", samp.hist_name
                draw_rest.append(samp.hist)


print "length of draw_rest", len(draw_rest)

print draw_first
draw_first.Draw("")

for ii in draw_rest: ii.Draw("same")


if options.islog:
    canvas.SetLogy()
    #    canvas.SetMaximum(1.5)
    #    canvas.SetMinimum(.0001)

leg = rt.TLegend(.4, .45 ,.75, .75)#canvas.BuildLegend()
for key in stacks.keys(): 
    for samp in stacks[key]:
        if samp.isData:
            leg.AddEntry(samp.hist, samp.label, 'pl')
        else:
            leg.AddEntry(samp.hist, samp.label, 'l')

leg.SetFillColor(0)
leg.SetLineColor(0)
leg.Draw("same")
CMS_lumi.CMS_lumi(canvas, 4, 0)


raw_input("RAW INPUT")

plot_name = "%s_%s" % ( options.file[:-5], options.var)
#canvas.SaveAs("%s.C" % plot_name)
#canvas.Print("%s.pdf" % plot_name)

