import os, sys, math
import ROOT as rt
from  optparse  import OptionParser

class sample:
    def __init__(self, file_name, xsec, label, stack):
        self.file_name = file_name
        self.xsec      = xsec
        self.label     = label
        self.stack     = stack
class cut:
    def __init__(self, tree, cutString, label):
        self.tree  = tree.strip()
        self.cutString   = cutString.strip()
        self.label = label

class cut_flow:
    def __init__(self, name, lumi):
        self.name = name
        self.lumi = lumi
        self.cuts = []
        self.samples = []
        self.cutflow = {}
        self.stackflow = {}

    def add_cut(self, thisCut): 

        if len(self.cuts) == 0:
            self.cuts.append(thisCut)
        else:
            thisCut.cutString = "(%s) && (%s)" % (self.cuts[-1].cutString, thisCut.cutString)
            self.cuts.append(thisCut)

        print "adding cut...", thisCut.cutString
    def add_sample(self, thisSample): self.samples.append(thisSample)
           
    def get_pass(self, sample, thisCut):
        try:
            file_in = rt.TFile(sample.file_name)
            tree        = file_in.Get(thisCut.tree)

            runStatTree = file_in.Get("runStats")
            n_events = tree.GetEntries(thisCut.cutString)
            scale    = self.lumi * float(sample.xsec) / float(runStatTree.GetEntries())

            n_pass = n_events * scale  
            n_err  = math.sqrt(n_events) * scale
        except AttributeError:
            return (-999,-999)
        return (n_pass, n_err)

    def fill_cutflow(self):        
       for sample in self.samples:
           for cut in self.cuts:
               (n_pass, n_err) = self.get_pass(sample,cut)
               cut_result = (cut, n_pass, n_err)
               if sample in self.cutflow.keys():
                   self.cutflow[sample].append(cut_result)
               else:
                   self.cutflow[sample] = []
                   self.cutflow[sample].append(cut_result)
                   

    def fill_stackflow(self):        
       for sample in self.samples:
           for cut_obj in self.cuts:
               stack = sample.stack
               (n_pass, n_err) = self.get_pass(sample, cut_obj)

               # the file has problems/ doesnt exist
               if n_pass == -999 and n_err == -999: continue

               # build the result to append
               cut_result = (cut_obj, n_pass, n_err)               

               #check if the stack is already in the stackflow
               if stack in self.stackflow.keys():
                   found_cut = False

                   # find the cut this sample should be included in
                   for o_cut_result in self.stackflow[stack]:

                       # index for the old cut within the stack array
                       index = self.stackflow[stack].index(o_cut_result)
                       (o_cut, o_n_pass, o_n_err) = o_cut_result

                       # you found the cut in the stack
                       if o_cut.label == cut_obj.label:
                           # combine the entry with this sample 
                           new_cut_result = (cut_obj, n_pass + o_n_pass, math.sqrt(n_err*n_err + o_n_err*o_n_err))
                           self.stackflow[stack][index] = new_cut_result                   
                           
                           found_cut = True
                           break

                   # the cut is new, but the stack is not
                   if not found_cut:
                       self.stackflow[stack].append(cut_result)                       
               else: # the stack is new to the stackflow
                   self.stackflow[stack] = [cut_result]

parser = OptionParser()
            
parser.add_option("-f","--files",dest="files",
		                    help="text file containing sample names",
		                    action="store",type="string",default="spring15.cfg")

parser.add_option("-c","--cutflow",dest="cutflow",
		                    help="text file containing cuts to apply to samples",
		                    action="store",type="string",default="cuts.cfg")

parser.add_option("-o", "--output", dest="output",
		                    help="output destination",
		                    action="store",type="string",default="test.root")

(options, args) = parser.parse_args()

cutflow = cut_flow("Cut Flow", 20.0*1000.0)

#output root file
output = rt.TFile(options.output, "RECREATE")

#read in the config file containing the samples to produce cut flow    
config_file = open(options.files, "r")
lines = map(lambda x: x.rstrip("\n"), config_file.readlines())
for line in lines:
    print line.split("|")
    (file_name, xSec, label, stack) = line.split("|")
    thisSample = sample(*line.split("|"))
    cutflow.add_sample(thisSample)

#read in the config file containing the samples to produce cut flow    
cutflow_file = open(options.cutflow, "r")
lines = map(lambda x: x.rstrip("\n"), cutflow_file.readlines())

for line in lines:
    print line.split("|")
    (tree, cutstring, label) = line.split("|")
    thisCut = cut(tree, cutstring, label)
    cutflow.add_cut(thisCut)

# fill everything in 
cutflow.fill_cutflow()
cutflow.fill_stackflow()

# print out the cut flow per sample
for sample in cutflow.cutflow.keys():
    print "---------%s--------".center(50) % sample.label
    for cutResult in cutflow.cutflow[sample]:
        (cut, n_pass, n_err) = cutResult
        print cut.label.ljust(25),": " "%2.1f +/- %2.1f".rjust(25) % (n_pass, n_err)

print "\n\n\n"

# print out the cut flow per stack
for stack in cutflow.stackflow.keys():
    print "---------%s--------".center(50) % stack
    for cutResult in cutflow.stackflow[stack]:
        (cut, n_pass, n_err) = cutResult
        print  cut.label.ljust(25),": %2.1f +/- %2.1f".rjust(25) % (n_pass, n_err)


