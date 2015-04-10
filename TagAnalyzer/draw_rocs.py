import ROOT as rt
from  optparse  import OptionParser
import rootlogon
rootlogon.style()
import CMS_lumi

rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptTitle(0)

parser = OptionParser()

parser.add_option("-f","--file",dest="file",
		                    help="text file containing samples and configuration",
		                    action="store",type="string")

(options, args) = parser.parse_args()


config_file = open(options.file, "r")
lines = map(lambda x: x.rstrip("\n"), config_file.readlines())

canvas = rt.TCanvas("c1")
roc_list = []
for line in lines[1:]:
    (file_name, metric_name, lineColor, lineStyle, lineWidth, label) = line.split("|")


    infile =  rt.TFile(file_name) 
    roc = infile.Get(metric_name)
    roc.SetLineColor(eval(lineColor))
    roc.SetLineStyle(int(lineStyle))
    roc.SetLineWidth(int(lineWidth))
    roc.SetFillColor(0)
    roc.SetTitle(label)
    
    roc_list.append(roc)

roc_list[0].Draw()
for roc in roc_list[1:]:
    roc.Draw("same")


leg = canvas.BuildLegend()    
leg.SetFillColor(0)
    
CMS_lumi.CMS_lumi(canvas, 4, 0)
    
raw_input("RAW_INPUT:")

canvas.SaveAs("%s.pdf" % options.file[:-5])
canvas.SaveAs("%s.C" % options.file[:-5])
