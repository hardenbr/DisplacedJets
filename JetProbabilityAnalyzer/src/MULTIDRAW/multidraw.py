import CMS_lumi, rootlogon
from  optparse  import OptionParser
import ROOT as rt
import os, sys, json
import array
parser = OptionParser()

parser.add_option("-l", "--list", dest="list",
                  help="list root files and objects to be draw",default="multidraw.json",
                  action="store",type="string")


parser.add_option( "--log", dest="log",
                  help="set canvas to be logscale",default=True,
                  action="store_true")


parser.add_option( "--xlabel", dest="xlabel",
                  help="label for the x axis",default="N Tracks",
                  action="store",type="string")


parser.add_option( "--ylabel", dest="ylabel",
                  help="label for the y axis",default="Fraction of Jets",
                  action="store",type="string")


parser.add_option( "--ymin", dest="ymin",
                  help="min y value",default=1E-5,
                  action="store",type="float")

parser.add_option( "--ymax", dest="ymax",
                  help="max y value",default=1,
                  action="store",type="float")


if __name__ == '__main__':
    rootlogon.style()

    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptTitle(0)

    canvas = rt.TCanvas("c1", "c1", 1024, 768)
    legend = rt.TLegend(.53, .41, .72, .63) 

    (options, args) = parser.parse_args()
    parser.print_help()

    json_data = open(options.list, "r").read()
    json      = json.loads(json_data)

    for sample in json["samples"]:
        index = json["samples"].index(sample)
        print "Processing sample: #", index, sample["path"] 
        
        # skip samples not meant to run
        if not sample["runSample"] : 
            print "Not Running sample...."
            continue

        #open the root file designated by the path
        thisFile = rt.TFile(sample["path"], "READ")        

        #global configuration parameters
        objectType = sample["type"]
        draw_options = sample["draw_options"]

        #parse the various modifiable parameters
        objects     = sample["objects"]
        labels       = sample["labels"]
        colors      = sample["colors"]
        line_styles = sample["line_styles"]        

        #parse the number of objects in the given root file to look at         
        nObj = len(objects)
        for ii in range(nObj):            
            #get the individual object
            obj = thisFile.Get(objects[ii])
            obj.Print()

            legend.AddEntry(obj, labels[ii], "l")
            
            # style changes
            obj.SetLineColor(eval("rt."+colors[ii]))
            obj.SetLineStyle(line_styles[ii])
            obj.SetLineWidth(2)
            obj.SetMarkerStyle(25)

            thisDrawOption = draw_options

            isFirst = index == 0 and ii == 0
            #append same to the draw options if necessary
            if isFirst: 
                thisDrawOption = "AP"                
            else:
            #only draw axis for the first sample first object
                thisDrawOption = "P" 

            canvas.cd()

            print "drawing", objects[ii], " with color ", colors[ii], "with options ", thisDrawOption
            obj.Draw(thisDrawOption)

            # special formatting for the first object in the first file
            if isFirst: 
                obj.GetXaxis().SetTitle(options.xlabel)
                obj.GetYaxis().SetTitle(options.ylabel)
                obj.GetYaxis().SetRangeUser(options.ymin, options.ymax)

            
    if options.log: canvas.SetLogy(1)

    legend.SetLineColor(0)
    legend.SetFillColor(0)
    legend.Draw("same")
    
    CMS_lumi.CMS_lumi(canvas, 4, 0)

    raw_input("RAW INPUT")
