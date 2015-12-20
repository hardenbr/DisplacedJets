#include "globalJetProbabilities.h"

globalJetProbabilities::globalJetProbabilities(const std::string& label, 
					       const std::string& stack,
					       const bool & isMC_ , 
					       const bool & isSig_ , 
					       const float& evWeight_, 
					       Json::Value probabilities,
					       const int & debug_) : 
  isMC(isMC_), 
  isSig(isSig_), 
  evWeight(evWeight_),
  debug(debug_) {

  // names for hte histograms generated by the label
  std::string	taggedJetHistName = "tagged_" + label;
  std::string	allJetHistName    = "allJets_" + label;
  std::string	effHistName	  = "effHist_" + label;
  std::string	effHistNameUp	  = "effHist_" + label + "Up";
  std::string	effHistNameDn	  = "effHist_" + label + "Dn";
  
  int nCat  = probabilities.size();
  if(nCat != 1) {
    std::cout << "[globalJetProbabilities] Only can handle one category for now" << std::endl;
    exit(1);
    }

  if(debug > 2) std::cout << "[globalJetProbabilities] Getting first category" << std::endl;

  // TODO!!!!: generalize to multiple cateogires
  std::string	catName		= "cat0";

  // get the first category
  Json::Value	cat1		  = probabilities[catName];
  // get the arrays for each of the histograms
  Json::Value	binning_json	  = cat1["binning"];
  Json::Value	allHist_json	  = cat1["all"];
  Json::Value	taggedHist_json	  = cat1["tagged"];
  Json::Value	effHist_json	  = cat1["eff"];
  Json::Value	effHistErrUp_json = cat1["effErrUp"];
  Json::Value	effHistErrDn_json = cat1["effErrDn"];

  if(debug > 2) std::cout << "[globalJetProbabilities] Building Binning : ";

  // make a vector of the binning
  std::vector<double> binningVec;
  nBins = binning_json.size() - 1 ;
  // build vectors from the arrays
  for(int bin = 0; bin < binning_json.size(); ++bin) {
    double binVal = binning_json[bin].asDouble();    
    if(debug > 2) std::cout << binVal << " ";
    binningVec.push_back(binVal);
  }

  if(debug > 2) std::cout << "\n[globalJetProbabilities] Building Temp Histograms" << std::endl;
  // build the jet histograms
  TH1D allJetHist_temp(allJetHistName.c_str(), "", nBins, &(binningVec[0]));
  TH1D taggedJetHist_temp(taggedJetHistName.c_str(), "", nBins, &(binningVec[0]));
  TH1D ratioHistEff_temp(effHistName.c_str(), "", nBins, &(binningVec[0]));
  TH1D ratioHistEffErrUp_temp(effHistNameUp.c_str(), "", nBins, &(binningVec[0]));
  TH1D ratioHistEffErrDn_temp(effHistNameDn.c_str(), "", nBins, &(binningVec[0]));

  if(debug > 2) std::cout << "[globalJetProbabilities] setting bin content of hists" << std::endl;
  // set the histogram bin content from the json
  for(int bin = 0; bin < nBins; ++bin) {
    allJetHist_temp.SetBinContent(bin+1, allHist_json[bin].asDouble());
    taggedJetHist_temp.SetBinContent(bin+1, taggedHist_json[bin].asDouble());
    // the efficiency histogram and the associated errors
    ratioHistEff_temp.SetBinContent(bin+1, effHist_json[bin].asDouble());
    ratioHistEffErrUp_temp.SetBinContent(bin+1, effHistErrUp_json[bin].asDouble());
    ratioHistEffErrDn_temp.SetBinContent(bin+1, effHistErrDn_json[bin].asDouble());
  }

  if(debug > 2) std::cout << "[globalJetProbabilities] Setting Pointers" << std::endl;

  // set the local hists to the temp hists
  taggedJetHist	    = taggedJetHist_temp;
  allJetHist	    = allJetHist_temp;
  ratioHistEff	    = ratioHistEff_temp;
  ratioHistEffErrUp = ratioHistEffErrUp_temp;
  ratioHistEffErrDn = ratioHistEffErrDn_temp;

  if(debug > 2) std::cout << "[globalJetProbabilities] Parsing cfg strings from config  " << std::endl;
  // set the string variables  from json
  binningVar	 = probabilities[catName].get("binningVar","1").asString();  
  categoryVar	 = probabilities[catName].get("categoryVar","caloJetPt").asString();  
  jetCutString	 = probabilities[catName].get("jetTagString","(1)").asString();
  eventCutString = probabilities[catName].get("eventTagString","(1)").asString();      
  if(debug > 2) {
    std::cout << "\t\t\t binvar: " << binningVar << std::endl;
    std::cout << "\t\t\t catVar: " << categoryVar << std::endl;
    std::cout << "\t\t\t cutstring: " << jetCutString << std::endl;
    std::cout << "\t\t\t eventCutString: " << eventCutString << std::endl;
  }  
}

// create a new global jet probability for outtputing to a JSON
globalJetProbabilities::globalJetProbabilities(const std::string& label, 
					       const std::string& stack,
					       const bool & isMC_ , 
					       const bool & isSig_ , 
					       const float& evWeight_, 
					       TTree*& tree,
					       const jetSelector& jetSel, 
					       const int & debug_) : 
  nBins(jetSel.getHistBinning().size() - 1),
  isMC(isMC_), 
  isSig(isSig_), 
  evWeight(evWeight_),
  debug(debug_) {
  
  // parse out the jet selection
  jetCutString	       = jetSel.getJetCutString();
  eventCutString       = jetSel.getEventCutString();
  baselineJetCutString = jetSel.getJetBaselineCutString();

  // build the fake rate histograms from the binning contained in the jetSelection
  histBinVals = jetSel.getHistBinning();
  catBinVals  = jetSel.getCatBinning();  

  // for now just make the single histogram for the fake Rate
  //  TH1D hist("genericHist", "genericHist", nBins, histBinVals); 
  std::string		    taggedJetHistName = "tagged_" + label;
  std::string		    allJetHistName    = "allJets_" + label;

  // get the names of the variables
  binningVar	   = jetSel.getBinningVarName();
  categoryVar	   = jetSel.getCategoryVarName();

  // form the draw string that will fill the histograms
  std::string	drawSelectedString = binningVar+">>" + taggedJetHistName+"(1000,-100,100)";
  std::string	drawAllString      = binningVar+">>" + allJetHistName + "(1000,-100,100)"; 

  if(debug > -1) { 
    std::cout << "-----------------------" << std::endl;
    std::cout << "binning variable string: " << binningVar << std::endl;
    std::cout << "jet selection string: " << jetCutString << std::endl;
    std::cout << "baseline jet  selection string: " << baselineJetCutString << std::endl;
    std::cout << "event selection string: " << eventCutString << std::endl;
  }

  // fill the histograms with the appropriate draw command
  // -- the tagged jets  
  tree->Draw(drawSelectedString.c_str(), ("(" + eventCutString + ") && (" + jetCutString + ")" + 
					  "&& (" + baselineJetCutString + ")").c_str(), "goff");
  // -- the not tagged jets
  tree->Draw(drawAllString.c_str(), ("(" + eventCutString + ") && (" + baselineJetCutString + ")").c_str(), "goff");

  // get the histograms from the pipe 
  taggedJetHist = (TH1D)*(TH1D*)gDirectory->Get(taggedJetHistName.c_str());
  allJetHist	= (TH1D)*(TH1D*)gDirectory->Get(allJetHistName.c_str());

  // check they are filled
  taggedJetHist.Print();
  allJetHist.Print();

  // rebin the histograms so they can be divided
  taggedJetHist = (TH1D)*(TH1D*)taggedJetHist.Rebin(nBins, "", &(histBinVals[0]));
  allJetHist	= (TH1D)*(TH1D*)allJetHist.Rebin(nBins, "", &(histBinVals[0]));

  // build the corresponding efficieny graph
  ratioGraph.BayesDivide(&taggedJetHist, &allJetHist);
  std::string graphName = "efficieny_" + label;
  ratioGraph.SetTitle(graphName.c_str());
  ratioGraph.SetName(graphName.c_str());

  if(debug > -1) std::cout << "Building a histogram based on the tgraph " << std::endl;         
  // build histographs based on the ratioGraph for the fake rate
  // copy the all jets histogram  and clear the values
  std::string effHistName = "effHist_" + label;

  // make copies off the histogram for the efficiency with the same binning
  ratioHistEff = (TH1D)*(TH1D*)allJetHist.Clone(effHistName.c_str());
  ratioHistEff.Reset();
  ratioHistEffErrUp = (TH1D)*(TH1D*)allJetHist.Clone((effHistName+"Up").c_str());
  ratioHistEffErrUp.Reset();
  ratioHistEffErrDn = (TH1D)*(TH1D*)allJetHist.Clone((effHistName+"Dn").c_str());
  ratioHistEffErrDn.Reset();

  // parse the central values and errors of the ratio Graph
  double *  xEffVals	    = ratioGraph.GetX();
  double *  yEffVals	    = ratioGraph.GetY();
  double *  yEffErrorUpVals = ratioGraph.GetEYhigh();
  double *  yEffErrorDnVals = ratioGraph.GetEYlow();
  int	    nPoints	    = ratioGraph.GetN();

  if(debug > -1) std::cout << "Filling the parsed values " << std::endl;         
  // fill the histogram with the values from the graph
  for(int ii = 0; ii < nPoints; ++ii) {
    // parse the efficiency 
    float   eff = yEffVals[ii];
    float   var = xEffVals[ii];

    // find and fill the bin in the histogram translation
    int	bin   = ratioHistEff.FindBin(var);    
    int	binUp = ratioHistEffErrUp.FindBin(var);    
    int	binDn = ratioHistEffErrDn.FindBin(var);    

    ratioHistEff.SetBinContent(bin, eff);    
    ratioHistEffErrDn.SetBinContent(bin, yEffErrorDnVals[ii]);    
    ratioHistEffErrUp.SetBinContent(bin, yEffErrorUpVals[ii]);    
  }   
  
  // check the histogram was printed
  ratioHistEff.Print();
}


Json::Value globalJetProbabilities::getProbabilitiesJSON() {
  // create the json 
  Json::Value event;     

  // loop over all categories for the fake rate
  for(int cat = 0; cat < catBinVals.size(); ++cat) {
    // build the arrays for the bins and values
    Json::Value binning(Json::arrayValue);
    Json::Value tagged(Json::arrayValue);
    Json::Value all(Json::arrayValue);
    Json::Value eff(Json::arrayValue);
    Json::Value effErrUp(Json::arrayValue);    
    Json::Value effErrDn(Json::arrayValue);    

    // loop over bins in the histogram for that category 
    for(int bin = 1; bin < histBinVals.size(); ++bin) {      

      // histbinvals is the array of the binning and begins with 0
      float binVal	= histBinVals[bin-1];
      // get the values from the histograms 
      // histograms first bin is numbered 1 (0 is the underflow bin)
      float tagVal	= taggedJetHist.GetBinContent(bin);
      float allJetVal	= allJetHist.GetBinContent(bin);
      float effVal	= ratioHistEff.GetBinContent(bin);
      float effValErrUp = ratioHistEffErrUp.GetBinContent(bin);
      float effValErrDn = ratioHistEffErrDn.GetBinContent(bin);
      
      // add the into the array 
      binning.append(Json::Value(binVal));
      tagged.append(Json::Value(tagVal));      
      all.append(Json::Value(allJetVal));      
      eff.append(Json::Value(effVal));
      effErrUp.append(Json::Value(effValErrUp));
      effErrDn.append(Json::Value(effValErrDn));
    }
    // add the last bin value for the binning array (which has nbins + 1 entries)
    binning.append(Json::Value(histBinVals[histBinVals.size()]));

    std::string catName = "cat" + std::to_string(cat); 
    
    // set info about how this was generated
    event[catName]["binningVar"]     = binningVar;
    event[catName]["binning"]	     = binning;
    event[catName]["jetTagString"]   = jetCutString;
    event[catName]["eventTagString"] = eventCutString;
    // encode the values of the histograms
    event[catName]["tagged"]	     = tagged;
    event[catName]["all"]	     = all;
    event[catName]["eff"]	     = eff;
    event[catName]["effErrUp"]	     = effErrUp;
    event[catName]["effErrDN"]	     = effErrDn;
  }
  
  return event;
} 

// does the lookup in the ratiograph for a specific jet configuration
double globalJetProbabilities::getJetFakeProbability(float binVariable, float catVar) {
  // find the bin and return the central value
  //  ratioHistEff.Print();
  if(debug > 20) std::cout << "[globalJetProb] binVariable: " << binVariable << " catVar " << catVar << std::endl;
  int bin = ratioHistEff.FindBin(binVariable);
  if(debug > 20) std::cout << "[globalJetProb] bin: " << bin << std::endl;
  double value = ratioHistEff.GetBinContent(bin);
  if(debug > 20) std::cout << "[globalJetProb] bin Content: " << value << std::endl;

  if(debug > 5 && binVariable == 1) std::cout << "JET HAS 1 TRACK -- binvar = " << 
				      binVariable << " fake probability = " << value << std::endl;

  return value;
}


// does the lookup in the ratiograph for a specific jet configuration
std::pair<double, double> globalJetProbabilities::getJetFakeProbabilityError(float binVariable, float catVar) {
  // find the bin and return the central value
  //  ratioHistEff.Print();
  if(debug > 5) std::cout << "[globalJetProb] binVariable: " << binVariable << " catVar " << catVar << std::endl;
  int bin = ratioHistEff.GetBin(binVariable);
  if(debug > 5) std::cout << "[globalJetProb] bin: " << bin << std::endl;
  double errUp = ratioHistEffErrUp.GetBinContent(bin);
  double errDn = ratioHistEffErrDn.GetBinContent(bin);
  if(debug > 5) std::cout << "[globalJetProb] bin errors up: " << errUp << " error down: "  <<  errDn << std::endl;

  std::pair<double, double> errors(errUp, errDn);
  return errors;
}

// does the lookup in the ratiograph for a specific jet configuration
// std::vector<float> globalJetProbabilities::getJetFakeProbabilityVector(TTree* tree, int evNum ) {
//   std::vector<float> probabilityVector;

//   // find the bin and return the central value
//   int bin = ratioHistEff->GetBin(binVariable);
//   // for now dont worry about the category variable
//   return ratioHistEff->GetBinContent(bin);
// }

void globalJetProbabilities::printHistStatus() {  
  // local variables
  std::cout << "----------global jet probabilities to apply object status----------" << std::endl;
  std::cout << "tagged histogram and all jet histogram" << std::endl;
  taggedJetHist.Print();
  allJetHist.Print();
  std::cout << "ratio graph for efficiency" << std::endl;
  ratioGraph.Print();
  std::cout << "efficiency histogram" << std::endl;
  ratioHistEff.Print();
  std::cout << "efficiency error histograms" << std::endl;
  ratioHistEffErrUp.Print();
  ratioHistEffErrDn.Print();
}


// accessors
TGraphAsymmErrors globalJetProbabilities::getRatioGraph() { return ratioGraph; }
TH1D globalJetProbabilities::getTaggedHist() { return taggedJetHist; }
TH1D globalJetProbabilities::getAllHist() { return allJetHist; }
TH1D globalJetProbabilities::getCentralEffHist() { return ratioHistEff;}
TH1D globalJetProbabilities::getCentralEffHistErrUp() { return ratioHistEffErrUp;}
TH1D globalJetProbabilities::getCentralEffHistErrDn() { return ratioHistEffErrDn;}

// build from preset global Probabilities
// globalJetProbabilities::globalJetProbabilities(const Json::Value & assignedProb) {  
// }
