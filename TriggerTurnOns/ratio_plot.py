import ROOT as rt
from  optparse  import OptionParser
import sys, os, array
import CMS_lumi 
import rootlogon
rootlogon.style()

parser = OptionParser()

rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptTitle(0)

parser.add_option("-f","--file",dest="file",
		                    help="text file containing samples and configuration",
		                    action="store",type="string")

parser.add_option("-o", "--output", dest="output",
		                    help="output destination",
		                    action="store",type="string")

parser.add_option("--label1", dest="label1",
		                    help="first cut label to apply to variable",
		                    action="store",type="string")

parser.add_option("--label2", dest="label2",
		                    help="first cut label to apply to variable",
		                    action="store",type="string")


parser.add_option("--cut1", dest="cut1",
		                    help="first cut to apply to variable",
		                    action="store",type="string")

parser.add_option("--cut2", dest="cut2",
		                    help="second cut to apply to variable",
		                    action="store",type="string")


parser.add_option("-v", "--variable", dest="var",
		                    help="variable to plot",
		                    action="store",type="string")
parser.add_option("-t", "--tree", dest="tree",
		                    help="tree for variable to be plotted from",
		                    action="store",type="string",default="jets")


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
		                    action="store",type="string",default = "N Events")

parser.add_option( "--nbins", dest="nbins",
		                    help="number of bins. Non-variable binning",
		                    action="store",type="int", default = 20)

parser.add_option("-l", "--lumi", dest="lumi",
		                    help="integrated luminosity for normalization",
		                    action="store",type="float", default = 1.0)

parser.add_option( "--norm1", dest="norm1",
		                    help="normalize to 1",
		                    action="store_true",default = False)

(options, args) = parser.parse_args()

#output root file
output = rt.TFile(options.output, "RECREATE")

variable = options.var
cut1 = options.cut1        
cut2 = options.cut1        
lumi = options.lumi
tree = options.tree
xmin = xmax = nbins = 0
varbins = None

xmin = options.xmin
xmax = options.xmax
nbins = options.nbins

output.cd()
#build the canvas
canvas = rt.TCanvas("plot","plot", 800, 600)
canvas.cd()
canvas.SetGridx()


thisFile = rt.TFile(options.file)
thisTree = thisFile.Get(options.tree)        

hist1_name  = "h1"
hist2_name  = "h2"

hist1 = rt.TH1F(hist1_name, options.label1, nbins, xmin, xmax)
hist1.SetMarkerStyle(21)
hist2 = rt.TH1F(hist2_name, options.label2, nbins, xmin, xmax)
hist2.SetMarkerStyle(21)



n1 = thisTree.Draw("%s>>%s" % (options.var, hist1_name), options.cut1)
n2 = thisTree.Draw("%s>>%s" % (options.var, hist2_name), options.cut2)

print "passsing cut1", n1
print "passsing cut2", n2


hist1.Sumw2()
hist2.Sumw2()

hist1.SetLineWidth(2)
hist1.SetLineColor(rt.kRed)
hist1.SetMarkerColor(rt.kRed)

hist2.SetLineWidth(2)
hist2.SetLineColor(rt.kBlack)
hist2.SetMarkerColor(rt.kBlack)

hist1_norm = hist1.Clone()
hist2_norm = hist2.Clone()


hist1_norm.GetXaxis().SetTitle(options.xlabel)
hist1_norm.GetYaxis().SetTitle(options.ylabel)

if options.norm1:
    hist1_norm.Scale( 1. / n1)
    hist2_norm.Scale( 1. / n2)

if hist1_norm.GetMaximum() > hist2_norm.GetMaximum():
    hist1_norm.Draw("e")
    hist1_norm.SetMarkerSize(2)
    hist2_norm.Draw("same")
else:
    hist2_norm.Draw("e")
    hist2_norm.SetMarkerSize(2)
    hist1_norm.Draw("same")

hist1_norm.GetXaxis().SetTitle(options.xlabel)
hist1_norm.GetYaxis().SetTitle(options.ylabel)


leg1 = canvas.BuildLegend()

CMS_lumi.CMS_lumi(canvas, 4, 0)

leg1.SetFillColor(0)
#leg1.SetLineColor(0)

canvas2 = rt.TCanvas("plot2","plot2", 800, 600)
canvas2.SetGridy()
canvas2.SetGridx()
canvas2.cd()

hcopy2 = hist2.Clone()
hcopy1 = hist1.Clone()

#hcopy2.Scale(1. / n2)
#hcopy1.Scale(1. / n1)

#hcopy2.Divide(hcopy1)
graph = rt.TGraphAsymmErrors()
graph.BayesDivide(hcopy2, hcopy1)

graph.SetLineWidth(2)
graph.SetLineColor(rt.kBlack)
graph.SetFillColor(0)


graph.GetXaxis().SetTitle(options.xlabel)
graph.GetYaxis().SetTitle("Ratio")

graph.SetTitle(options.label2)
graph.SetMarkerStyle(21)
graph.SetMarkerColor(rt.kBlack)
graph.Draw("ap")

graph.GetXaxis().SetTitle(options.xlabel)

#leg2 = canvas2.BuildLegend()

CMS_lumi.CMS_lumi(canvas2, 4, 0)

#leg2.SetFillColor(0)
#leg2.SetLineColor(0)
raw_input("RAWINPUT")
