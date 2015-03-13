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

parser.add_option( "--genmatch", dest="genmatch",
		                    help="Do generator matching for the signal MC",
		                    action="store_true", default=False)


parser.add_option( "--xlabel", dest="xlabel",
		                    help="label for the x axis",
		                    action="store",type="string",default = "Exactly N Jet Tags")

parser.add_option( "--ylabel", dest="ylabel",
		                    help="label for the y axis",
		                    action="store",type="string",default = "N Events")

parser.add_option("-l", "--lumi", dest="lumi",
		                    help="integrated luminosity for normalization",
		                    action="store",type="float", default = 1.0)

(options, args) = parser.parse_args()


class sample:
    def __init__(self, file_name, tree_name, isSignal, metricID, metricCut, xsec, fillColor, fillStyle, lineStyle, lineWidth, stack, label):

        # set configuration
        self.file_name = file_name 
        self.tree_name = tree_name 
        self.xsec = float(xsec)
        self.fillColor = fillColor
        self.fillStyle = fillStyle
        self.lineStyle = lineStyle
        self.lineWidth = lineWidth 
        self.stack = stack
        self.label = label
        self.hist_name = file_name.split(".")[0] + "_" + tree_name + metricCut
        self.isSignal = int(isSignal)
        self.metricID = int(metricID)
        self.metricCut = float(metricCut)


        self.eventDict = {}
        #build the corresponding histogram
        self.hist = rt.TH1F(self.hist_name, self.label, 8, -.5, 7.5)
        if int(fillStyle) != 0:
            eval("self.hist.SetFillColor(rt.%s)" % fillColor)
        eval("self.hist.SetFillStyle(%s)" % fillStyle)
        eval("self.hist.SetLineColor(rt.%s)" % fillColor)
        self.hist.SetLineWidth(int(lineWidth))        


#output root file
output = rt.TFile(options.output, "RECREATE")

config_file = open(options.file, "r")
lines = map(lambda x: x.rstrip("\n"), config_file.readlines())

#dictionary from filename to sample class
samples = {} 

#dictionary from stack string to array of samples corresponding to stack
stacks = {}

#ignore the first line for labels
for line in lines[1:]:
    print line.split("|")

    #fill all the configuration info and build the sample class
    (file_name, tree_name, isSignal, metricID, metricCut, xSec, fillColor, fillStyle, lineStyle, lineWidth, stack, label) = line.split("|")
    samples[file_name] = sample(file_name, tree_name, isSignal, metricID, metricCut, xSec, fillColor, fillStyle, lineStyle, lineWidth, stack, label)

    #check if the stack already exists
    if stack not in stacks:
        stacks[stack] = [samples[file_name]]
    else:
        stacks[stack].append(samples[file_name])

# for stacked objects (assume only one object per stack for now)
# determine the number of tag counts        
for key in stacks.keys(): 
    for samp in stacks[key]:
        
        thisFile = rt.TFile(samp.file_name)
        tree = thisFile.Get(samp.tree_name)                
        tree.Draw(">>iterlist", "", "entrylist")
        iev = 0

        itlist = rt.gDirectory.Get("iterlist")
        
        #loop over all events in the given tree
        while iev < tree.GetEntries():
            iev += 1
            entry = itlist.Next()
            tree.GetEntry(entry)

            event = tree.event
            ls = tree.lumi
            run = tree.run
            genMatch = tree.genMatch

            #get the metric we are looking for            
            exec_string = "metric = tree.metric[%i]" % (int(samp.metricID)-1) #zero indexing
            exec(exec_string) 

            # apply the cut and gen matching if asked for            
            pass_genmatch = genMatch > 0 or not options.genmatch
            evTuple = (run, ls, event)
            if metric > samp.metricCut and (not samp.isSignal or pass_genmatch):
                if evTuple not in samp.eventDict.keys():
                    samp.eventDict[evTuple] = 1
                else:
                    samp.eventDict[evTuple] += 1
            else:
                if evTuple not in samp.eventDict.keys():
                    samp.eventDict[evTuple] = 0
                                        
        #print samp.file_name, samp.tree_name, samp.eventDict
        #loop over the event dict and fill the histogram with number of tags
        for key in samp.eventDict.keys():
            ntags = samp.eventDict[key]
            samp.hist.Fill(ntags, 1)

        #        samp.hist.Sumw2()
        samp.hist.Scale(samp.xsec * options.lumi / float(len(samp.eventDict)))
        
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

draw_first.Draw("")
draw_first.GetXaxis().SetTitle(options.xlabel)
draw_first.GetYaxis().SetTitle(options.ylabel)
draw_first.GetYaxis().SetRangeUser(.0001,1)

for ii in draw_rest: ii.Draw("same")

leg = canvas.BuildLegend()
leg.SetFillColor(0)
leg.SetLineColor(0)
CMS_lumi.CMS_lumi(canvas, 4, 0)

raw_input("RAW INPUT")
canvas.Print("%s.pdf" % options.file)
