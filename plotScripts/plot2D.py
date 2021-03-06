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
    return list

parser = OptionParser()


parser.add_option("-f","--file",dest="file",
		                    help="text file containing samples and configuration",
		                    action="store",type="string")

parser.add_option( "--ivfvtxmatch", dest="ivfvtxmatch",
		                    help="Do generator vertex matching for the signal MC",
		                    action="store_true", default=False)

parser.add_option( "--svvtxmatch", dest="svvtxmatch",
		                    help="Do generator vertex matching for the signal MC",
		                    action="store_true", default=False)

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
		                    action="store_true", default=False)

parser.add_option( "--genmatch", dest="genmatch",
		                    help="Do generator matching for the signal MC",
		                    action="store_true", default=False)


parser.add_option( "--square", dest="square",
		                    help="use a square canvas",
		                    action="store_true" )

parser.add_option( "--black", dest="black",
		                    help="draw black and white",
		                    action="store_true" )



parser.add_option( "--ymin", dest="ymin",
		                    help="minimum y for variable",
		                    action="store",type="float",default=0)

parser.add_option( "--ymax", dest="ymax",
		                    help="maximum y for variable",
		                    action="store",type="float",default = 100)


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

parser.add_option( "--nxbins", dest="nybins",
		                    help="number of x bins. Non-variable binning",
		                    action="store",type="int", default = 20)

parser.add_option( "--nybins", dest="nxbins",
		                    help="number of y bins. Non-variable binning",
		                    action="store",type="int", default = 20)

parser.add_option("-l", "--lumi", dest="lumi",
		                    help="integrated luminosity for normalization",
		                    action="store",type="float", default = 1.0)




(options, args) = parser.parse_args()

#output root file
output = rt.TFile(options.output, "RECREATE")

variable = options.var
cut = options.cut        
lumi = options.lumi
tree = options.tree
xmin = xmax = nxbins = 0
ymin = ymax = nybins = 0
varbins = None

if not options.var_bin:
    xmin = options.xmin
    xmax = options.xmax
    nxbins = options.nxbins

    ymin = options.ymin
    ymax = options.ymax
    nybins = options.nybins

else:
    varbins = makebins(options.xmin, options.xmax, .1, .3)

class sample:
    def __init__(self, file_name, isSignal, xsec, fillColor, fillStyle, lineWidth, stack, label):

        # set configuration
        self.file_name = file_name 
        self.xsec = xsec 
        self.fillColor = fillColor
        self.fillStyle = fillStyle
        self.lineWidth = lineWidth 
        self.stack = stack
        self.label = label
        self.hist_name = file_name.split(".")[0]
        self.isSignal = int(isSignal)
        
        #build the corresponding histogram
        self.hist = None
        if not options.var_bin:
            self.hist = rt.TH2F(self.hist_name, self.label, nxbins, xmin, xmax, nybins, ymin, ymax)

#no option for variable binning yet
#         else: 
#             self.hist = rt.TH2F(self.hist_name, self.label, array.array("d",varsbins))

#         eval("self.hist.SetFillColor(rt.%s)" % fillColor)
#         eval("self.hist.SetFillStyle(%s)" % fillStyle)
#         eval("self.hist.SetLineColor(rt.%s)" % fillColor)

        eval("self.hist.SetFillColor(0)")
        eval("self.hist.SetFillStyle(0)")
        eval("self.hist.SetLineColor(0)")

        self.hist.SetLineWidth(int(lineWidth))        

    
config_file = open(options.file, "r")
lines = map(lambda x: x.rstrip("\n"), config_file.readlines())

#dictionary from filename to sample class
samples = {} 

#dictionary from stack string to array of samples corresponding to stack
stacks = {}

#ignore the first line for labels
for line in lines[1:]:
    print line.split("|")
    (file_name, isSignal, xSec, fillColor, fillStyle, lineWidth, stack, label) = line.split("|")
    samples[file_name] = sample(file_name, isSignal, xSec, fillColor, fillStyle, lineWidth, stack, label)

    #check if the stack already exists
    if stack not in stacks:
        stacks[stack] = [samples[file_name]]
    else:
        stacks[stack].append(samples[file_name])



#for stacked objects (assume only one object per stack for now)
for key in stacks.keys(): 
    for samp in stacks[key]:
        thisFile = rt.TFile(samp.file_name)
        thisTree = thisFile.Get(options.tree)        
        output.cd()
        
        thisCut = options.cut
        #add in the gen amtching requirement to the cut
        if samp.isSignal and options.genmatch:
            print samp.file_name, isSignal
            thisCut = "(" + options.cut + ") && ( caloGenMatch > 0 )" 
        if samp.isSignal and options.ivfvtxmatch:
            print samp.file_name, isSignal
            thisCut = "(" + options.cut + ") && ( jetIVFGenVertexMatched > 0 )" 
        if samp.isSignal and options.svvtxmatch:
            print samp.file_name, isSignal
            thisCut = "(" + options.cut + ") && ( jetSvGenVertexMatched > 0 )" 


        print "Sample: ", samp.file_name,  "Cut: ", thisCut

        #fill each histogram 
        nevents =  thisTree.Draw("%s>>%s" % (options.var, samp.hist_name), thisCut,"colz")
        samp.hist.Scale( float(options.lumi) * float(samp.xsec) / float(nevents))
        samp.hist.GetXaxis().SetTitle(options.xlabel)
        samp.hist.GetYaxis().SetTitle(options.ylabel)


output.cd()
#build the canvas
canvas = None
if options.square:
    canvas = rt.TCanvas("plot","plot", 1024, 1024)
else:
    canvas = rt.TCanvas("plot","plot", 1024, 768)
canvas.cd()

draw_first = None
draw_rest = []
max = -1

# for now there should be only 1 stack and 1 sample per plot
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

draw_first.Draw("" if options.black else "colz")
draw_first.GetXaxis().SetLimits(options.xmin, options.xmax)
draw_first.GetYaxis().SetLimits(options.ymin, options.ymax)


leg = canvas.BuildLegend()
leg.SetFillColor(0)
leg.SetLineColor(0)
CMS_lumi.CMS_lumi(canvas, 4, 0)

raw_input("RAW INPUT")

split_var = options.var.split(":") 

plot_name = "%s_%s_%s.pdf" % (options.file, split_var[0], split_var[1])

#canvas.Print(plot_name)
