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

parser.add_option("-o", "--output", dest="output",
		                    help="output destination",
		                    action="store",type="string")

parser.add_option("-c", "--cut", dest="cut",
		                    help="cut to apply to variable",
		                    action="store",type="string")

parser.add_option("-t", "--tree", dest="tree",
		                    help="tree for variable to be plotted from",
		                    action="store",type="string",default="eventInfo")

parser.add_option("-x", "--threshold", dest="thres",
		                    help="which threhsold to use for taggs",
		                    action="store",type="int",default=1)



# parser.add_option( "--genmatch", dest="genmatch",
# 		                    help="Do generator matching for the signal MC",
# 		                    action="store_true", default=False)

# parser.add_option( "--ivfvtxmatch", dest="ivfvtxmatch",
# 		                    help="Do generator vertex matching for the signal MC",
# 		                    action="store_true", default=False)

# parser.add_option( "--svvtxmatch", dest="svvtxmatch",
# 		                    help="Do generator vertex matching for the signal MC",
# 		                    action="store_true", default=False)

parser.add_option( "--xlabel", dest="xlabel",
		                    help="label for the x axis",
		                    action="store",type="string",default = "variable")

parser.add_option( "--ylabel", dest="ylabel",
		                    help="label for the y axis",
		                    action="store",type="string",default = "N Jets")

parser.add_option("-l", "--lumi", dest="lumi",
		                    help="integrated luminosity for normalization",
		                    action="store",type="float", default = 20.0)

(options, args) = parser.parse_args()

#output root file
output = rt.TFile(options.output, "RECREATE")

#variable = options.var
#cut = options.cut        
lumi = options.lumi
#tree = options.tree

nbins = 16
xmin = 0
xmax = 16



class sample:
    def __init__(self, file_name, tree_name, isSignal, xsec, fillColor, fillStyle, lineStyle, lineWidth, stack, label):

        # set configuration
        self.file_name = file_name 
        self.tree_name = tree_name 
        print "XSEC", xsec, "EVAL(XSEC)", eval(xsec)
        self.xsec = eval(xsec)
        print xsec
        self.fillColor = fillColor
        self.fillStyle = fillStyle
        self.lineStyle = lineStyle
        self.lineWidth = lineWidth 
        self.stack = stack
        self.label = label
        self.hist_name = file_name.split(".")[0] + "_" + tree_name
        self.isSignal = int(isSignal)
        
        #build the corresponding histogram
        self.hist = rt.TH1F(self.hist_name, self.label, nbins, xmin, xmax)
        self.hist_total = rt.TH1F(self.hist_name+"total", self.label, 4, 0, 4)

        
        #setup the bins
        self.hist.Fill("1 No Vertex",0)
        self.hist.Fill("2 No Vertex",0)
        self.hist.Fill("3 No Vertex",0)
        self.hist.Fill("4+ No Vertex",0)
        self.hist.Fill("1 Short",0)
        self.hist.Fill("2 Short",0)
        self.hist.Fill("3 Short",0)
        self.hist.Fill("4+ Short",0)
        self.hist.Fill("1 Medium",0)
        self.hist.Fill("2 Medium",0)
        self.hist.Fill("3 Medium",0)
        self.hist.Fill("4+ Medium",0)
        self.hist.Fill("1 Long",0)
        self.hist.Fill("2 Long",0)
        self.hist.Fill("3 Long",0)
        self.hist.Fill("4+ Long",0)
        self.hist_total.Fill("1 Total",0)
        self.hist_total.Fill("2 Total",0)
        self.hist_total.Fill("3 Total",0)
        self.hist_total.Fill("4+ Total",0)

        # binned in number
        # self.hist.Fill("1 No Vertex",0)
        # self.hist.Fill("1 Short",0)
        # self.hist.Fill("1 Medium",0)
        # self.hist.Fill("1 Long",0)
        # self.hist.Fill("2 No Vertex",0)
        # self.hist.Fill("2 Short",0)
        # self.hist.Fill("2 Medium",0)
        # self.hist.Fill("2 Long",0)
        # self.hist.Fill("3 No Vertex",0)
        # self.hist.Fill("3 Short",0)
        # self.hist.Fill("3 Medium",0)
        # self.hist.Fill("3 Long",0)
        # self.hist.Fill("4+ No Vertex",0)
        # self.hist.Fill("4+ Short",0)
        # self.hist.Fill("4+ Medium",0)
        # self.hist.Fill("4+ Long",0)
        # self.hist_total.Fill("1 Total",0)
        # self.hist_total.Fill("2 Total",0)
        # self.hist_total.Fill("3 Total",0)
        # self.hist_total.Fill("4+ Total",0)



        if int(fillStyle) != 0:
            eval("self.hist.SetFillColor(rt.%s)" % fillColor)

        eval("self.hist.SetFillStyle(%s)" % fillStyle)
        eval("self.hist.SetLineColor(rt.%s)" % fillColor)
        eval("self.hist.SetLineStyle(%s)" % lineStyle)
        eval("self.hist_total.SetFillStyle(%s)" % fillStyle)
        eval("self.hist_total.SetLineColor(rt.%s)" % fillColor)
        eval("self.hist_total.SetLineStyle(%s)" % lineStyle)

        self.hist.SetLineWidth(int(lineWidth))        
        self.hist_total.SetLineWidth(int(lineWidth))        

    
config_file = open(options.file, "r")
lines = map(lambda x: x.rstrip("\n"), config_file.readlines())

#dictionary from filename to sample class
samples = {} 

#dictionary from stack string to array of samples corresponding to stack
stacks = {}

#ignore the first line for labels
for line in lines[1:]:
    (file_name, tree_name, isSignal, xSec, fillColor, fillStyle, lineStyle, lineWidth, stack, label) = line.split("|")
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
        eventTree = thisFile.Get("runStats")
        output.cd()
        
        thisCut = options.cut
        # #add in the gen amtching requirement to the cut
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

        print "Sample: ", samp.file_name,  "Cut: ", thisCut

        #fill the iterlist
        n_pass = 0

        thisTree.Draw(">>iterlist", thisCut, "entrylist")
        itlist = rt.gDirectory.Get("iterlist")
        iev = 0

        # loop over events in the tree
        while iev < thisTree.GetEntries():
            if iev % 1000 == 0: print "Filling Tree...", iev
            iev += 1
            entry = itlist.Next()
            thisTree.GetEntry(entry)

#            print list(thisTree.eventNNoVertexTags)
#            print list(thisTree.eventNShortTags)

            # parse out the tags
            n_novtx = thisTree.eventNNoVertexTags[options.thres]
            n_short = thisTree.eventNShortTags[options.thres]
            n_med   = thisTree.eventNMediumTags[options.thres]
            n_long  = thisTree.eventNLongTags[options.thres]
            n_total = n_novtx + n_short + n_med + n_long
            
            # there needs to be at least one tag to keep the event
            if n_total == 0: 
                continue
            else:
                n_pass+=1

            w_novtx = float(n_novtx) / float(n_total)
            w_short = float(n_short) / float(n_total)
            w_med   = float(n_med)   / float(n_total)
            w_long  = float(n_long)  / float(n_total)

            # individual tags
            if int(n_novtx) >= 4:  n_novtx = "4+"
            if int(n_short) >= 4:  n_short = "4+"
            if int(n_med) >= 4:    n_med   = "4+"
            if int(n_long) >= 4:   n_long  = "4+"
            # sum tags
            if int(n_total) >= 4:  n_total = "4+"
            
            # individual tags
            if n_novtx == "4+" or int(n_novtx) > 0: samp.hist.Fill("%s No Vertex" % str(n_novtx), w_novtx)
            if n_short == "4+" or int(n_short) > 0: samp.hist.Fill("%s Short" % str(n_short), w_short)
            if n_med == "4+" or int(n_med) > 0:     samp.hist.Fill("%s Medium" % str(n_med), w_med)
            if n_long == "4+" or int(n_long) > 0:   samp.hist.Fill("%s Long" % str(n_long), w_long)
            # sum tags
            if n_total == "4+" or int(n_total) > 0: samp.hist_total.Fill("%s Total" % str(n_total), 1)        

        #nevents = thisTree.GetEntries(thisCut)
        n_analyzed = eventTree.GetEntries()
        hist_scale =  float(options.lumi) * 1000 * float(samp.xsec) /  n_analyzed #* float(n_pass)
#        samp.hist.Sumw2()
#        samp.hist_total.Sumw2()

        samp.hist.Scale( hist_scale )
        samp.hist_total.Scale( hist_scale )
        samp.hist.GetXaxis().SetTitle(options.xlabel)
        samp.hist.GetYaxis().SetTitle(options.ylabel)

        samp.hist_total.GetXaxis().SetTitle(options.xlabel)
        samp.hist_total.GetYaxis().SetTitle(options.ylabel)

        samp.hist_total.SetMarkerStyle(4)
        samp.hist.SetMarkerStyle(4)

        samp.hist_total.SetMarkerStyle(4)

        eval("samp.hist.SetMarkerColor(rt.%s)" % samp.fillColor)
        eval("samp.hist_total.SetMarkerColor(rt.%s)" % samp.fillColor)

output.cd()
#build the canvas
canvas = rt.TCanvas("plot","plot", 1500, 768)
canvas.cd()

draw_first = None
draw_rest = []
max = -1
for key in stacks.keys(): 
    for samp in stacks[key]:
        #draw each histogram
        thisMax = samp.hist.GetMaximum() 
        print samp.hist_name, " INTEGRAL ", samp.hist.Integral()
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
draw_first.Draw("HIST")

for ii in draw_rest: ii.Draw("HISTsame")


leg = canvas.BuildLegend()
leg.SetFillColor(0)
leg.SetLineColor(0)
CMS_lumi.CMS_lumi(canvas, 4, 0)


raw_input("RAW INPUT")

canvas_tot = rt.TCanvas("plot","plot", 1000, 768)
canvas_tot.cd()
draw_first = None
draw_rest = []
max = -1
for key in stacks.keys(): 
    for samp in stacks[key]:
        #draw each histogram
        thisMax = samp.hist_total.GetMaximum() 
        print samp.hist_name, " INTEGRAL ", samp.hist_total.Integral()        
        if thisMax > max: 
            if draw_first != None and samp.hist_total not in draw_rest:
                print "shifting down", samp.hist_name
                draw_rest.append(draw_first)

            max = thisMax
            draw_first = samp.hist_total
        else:
            if samp.hist_total not in draw_rest:
                print "appending", samp.hist_name
                draw_rest.append(samp.hist_total)

draw_first.Draw("HIST")
for ii in draw_rest: ii.Draw("HISTsame")



leg_tot = canvas_tot.BuildLegend()
leg_tot.SetFillColor(0)
leg_tot.SetLineColor(0)
CMS_lumi.CMS_lumi(canvas_tot, 4, 0)

raw_input("RAW INPUT")
#plot_name = "%s_tags" % options.file[:-5]
#canvas.SaveAs("%s.C" % plot_name)
#canvas.Print("%s.pdf" % plot_name)

