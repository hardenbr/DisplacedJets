from  optparse  import OptionParser
#import ROOT as rt
import os, sys, json

parser = OptionParser()

parser.add_option("-j", "--json", dest="json",
                  help="json of samples to combined",default="composite.json",
                  action="store",type="string")


parser.add_option("-o", "--output", dest="output",
                  help="name of the composite json to output",default="compositeProbabilitiesResult.json",
                  action="store",type="string")


parser.add_option("--xsec", dest="useXsec",
                  help="""use the xsec to weight the efficiencies. 
               Otherwise the eventWeight of the sample will be used""",default=False,
                  action="store_true")

# sample corresponding to a container of categories
class cat_container:
    def __init__(self, name, weight, xsec, doSubtract):
        self.name        = name
        self.xsec        = xsec
        self.eventWeight = eventWeight
        self.doSubtract  = doSubtract
        self.categories  = []

#
class category:
    def __init__(self, name):
        self.name = name
        self.all            = []
        self.tagged         = []
        self.binning        = []
        self.binningVar     = None
        self.eff            = []
        self.effErrDn       = []
        self.effErrUp       = []
        self.eventTagString = ""
        self.jetTagString   = ""

def scale_list(list1, scale):
    scaled_list = []
    for ii in list1: scaled_list.append(ii * scale)
    return scaled_list

def add_lists(list1, list2, scale):
    if list1 == [] and len(list2) != 0: 
        list1 = [0]*len(list2)
    if list2 == [] and len(list1) != 0: 
        list2 = [0]*len(list1)

    if len(list1) != len(list2): 
        print "LIST LENGTH DOES NOT MATCH"
        exit(1)

    sum_list = []
    for ii in range(len(list1)):
        sum_list.append(list1[ii] + scale * list2[ii])

    return sum_list

if __name__ == '__main__':
    (options, args) = parser.parse_args()
    parser.print_help()
    
    # parse the json containing the probability jsons to use
    setup_json_data = open(options.json, "r").read()
    setup_json      = json.loads(setup_json_data)

    # parse the samples to add and subtract
    sum_samples      = setup_json["sum"]
    subtract_samples = setup_json["subtract"]





    # if len(subtract_samples) != 0 and options.useXsec:
    #     print "Cannot use xsec weighting if a sample is being subtracted?    

    allSamples = []
    total_xsec = 0
    total_weight = 0

    for sample in sum_samples + subtract_samples:
        prob_json_data = open(sample, "r").read() 
        sample_json    = json.loads(prob_json_data)

        # sample parameters
        eventWeight = sample_json["eventWeight"]
        xsec        = sample_json["xsec"]
        doSubtract  = sample in subtract_samples
        
        # keep track of the total xsec
        total_xsec   += xsec
        total_weight += eventWeight

        # build the sample
        thisSample  = cat_container(sample, eventWeight, xsec, doSubtract)

        # check for each category
        for item in sample_json:
            if "cat" not in item: continue # make sure it is a category
            
            cat = item
            #build the new category
            thisCat = category(item)

            print "\n----------------------"
            print "sample name", sample
            print "sample eff", sample_json[cat]["eff"]
            print "sample all", sample_json[cat]["binning"]
            print "sample tagged", sample_json[cat]["tagged"]

            #parse the necessary information
            thisCat.binning  = sample_json[cat]["binning"]
            thisCat.all      = sample_json[cat]["all"]
            thisCat.tagged   = sample_json[cat]["tagged"]
            thisCat.eff      = sample_json[cat]["eff"]
            thisCat.effErrUp = sample_json[cat]["effErrUp"]
            thisCat.effErrDn = sample_json[cat]["effErrDn"]
            
            # selection information
            thisCat.binningVar     = sample_json[cat]["binningVar"]
            thisCat.eventTagString = sample_json[cat]["eventTagString"]
            thisCat.jetTagString   = sample_json[cat]["jetTagString"]

            # add it to the sample
            thisSample.categories.append(thisCat)

        # add the sample
        allSamples.append(thisSample)

    # total weighted averages
    total_container = cat_container("total", 1, 1, False)
    nCats           = len(allSamples[0].categories)

    # loop over the categories and add the contribution from each sample
    for catNum in range(nCats):
        isFirstCat = (catNum  == 0)
        total_cat = category("cat%i" % catNum)

        # loop over categories within the sample
        for sample in allSamples:
            # check if this is the first sample for this category
            isFirstSample = allSamples.index(sample) == 0
            cat     = sample.categories[catNum]

            # if its the first sample for the category just take the nominal values
            if isFirstSample:
                total_cat.binning        = cat.binning 
                total_cat.binningVar     = cat.binningVar
                total_cat.eventTagString = cat.eventTagString
                total_cat.jetTagString   = cat.jetTagString
                
            else: # make sure this sample is compatible with the rest of them
                binDiff       = not total_cat.binning == cat.binning
                binVarDiff    = not total_cat.binningVar == cat.binningVar
                tagStringDiff = not total_cat.eventTagString == cat.eventTagString
                jetStringDiff = not total_cat.jetTagString == cat.jetTagString
                nCatsDiff     = not nCats == len(sample.categories)
                
                # make a check that all the sample probabilities are compatible
                if binDiff or binVarDiff or tagStringDiff or jetStringDiff or nCatsDiff:
                    print "DIFFERENCES IN GLOBAL VARIABLES"
                    exit(1)
            
            # if everything matchs add it
            sign = -1 if doSubtract else 1            
            scale = sample.xsec if options.useXsec else sample.eventWeight

            total_cat.eff      = add_lists(total_cat.eff, cat.eff, sign*scale)
            total_cat.effErrUp = add_lists(total_cat.effErrUp, cat.effErrUp, sign*scale)
            total_cat.effErrDn = add_lists(total_cat.effErrDn, cat.effErrDn, sign*scale)
            total_cat.all      = add_lists(total_cat.all, cat.all, sign*scale)
            total_cat.tagged   = add_lists(total_cat.tagged, cat.tagged, sign*scale)

        print "\n-----categorie totale devant scaling------- \n"
        print "total eff", total_cat.eff, "\n"
        print "total tagged", total_cat.tagged, "\n"
        print "total all", total_cat.all, "\n"
        print "\n-----categorie totale apres scaling------- \n"

        #determine which normalization to use
        norm_scale = 1./total_xsec if options.useXsec else 1./total_weight

        #scale out the total cross section given the weighting
        total_cat.eff    = scale_list(total_cat.eff, norm_scale)
        total_cat.effErrUp  = scale_list(total_cat.effErrUp, norm_scale)
        total_cat.effErrDn  = scale_list(total_cat.effErrDn, norm_scale)
        total_cat.all    = scale_list(total_cat.all, norm_scale)
        total_cat.tagged = scale_list(total_cat.tagged, norm_scale)

        # add the category to the total container
        total_container.categories.append(total_cat)
    
    #now output the total container        
    json_data_dict = {}
    json_data_dict["eventWeight"] = total_container.eventWeight
    json_data_dict["xsec"]        = total_container.xsec

    # add the information into the dictionary
    for cat in total_container.categories:
        json_data_dict[cat.name]                   = {}
        json_data_dict[cat.name]["binning"]        = cat.binning
        json_data_dict[cat.name]["binningVar"]     = cat.binningVar
        json_data_dict[cat.name]["all"]            = cat.all
        json_data_dict[cat.name]["tagged"]         = cat.tagged
        json_data_dict[cat.name]["eff"]            = cat.eff
        json_data_dict[cat.name]["effErrUp"]       = cat.effErrUp
        json_data_dict[cat.name]["effErrDn"]       = cat.effErrDn
        json_data_dict[cat.name]["eventTagString"] = cat.eventTagString
        json_data_dict[cat.name]["jetTagString"]   = cat.jetTagString

        print "total xsec", total_xsec, "\n"
        print "total eff", cat.eff, "\n"
        print "total effErrUp", cat.effErrUp, "\n"
        print "total effErrDn", cat.effErrDn, "\n"
        print "total tagged", cat.tagged, "\n"
        print "total all", cat.all, "\n"

    
    json_string = json.dumps(json_data_dict, sort_keys=True, indent=4, separators=(',', ': '))
    
    output_name = options.output
    if ".json" not in options.output: output_name += ".json"

    # write the file
    output_file = open(output_name, "w")
    output_file.write(json_string)
