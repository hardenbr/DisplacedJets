//#include "jetTagAnalyzer.h"
#include "json/json.h"
#include "json/json-forwards.h"
//#include "jetSelector.h"
//#include "globalJetProbabilities.h"
#include "jetProbabilityMasterComputer.h"
#include <iostream>
#include <fstream>
#include <string.h>  
#include "TCanvas.h"
#include "TColor.h"

int main(int argc, char* argv[]) {

  // names of the json files to be parsed
  std::string sample_json;
  std::string selection_json;
  std::string run_config_json;
  std::string globalProb_json;
  // name of the json containing the global probabilities
  std::string globalProbabilities;
  // name of the tree to use in each sample
  std::string treeName;

  // debug status
  int debug = 0;

  // give default values for the configurations if not specified
  if (argc == 2) {
    sample_json	    = "sample.json";
    selection_json  = "selection.json";
    run_config_json = "run_config.json";
    globalProb_json = argv[1];
  }
  else {
    sample_json	    = "sample.json";
    selection_json  = "selection.json";
    run_config_json = "run_config.json";
    globalProb_json = "none";
  }
  // else { //otherwise parse from the arguments
  //   sample_json	    = argv[1];
  //   selection_json  = argv[2];
  //   run_config_json = argv[3];
  //   globalProb_json = argv[4];
  // }

  // check if global probabilities were provided
  bool probProvided = (globalProb_json != "none");
  if(debug > -1){
    std::cout << "sample      json name: " << sample_json << std::endl;
    std::cout << "selection   json name: " << selection_json << std::endl;
    std::cout << "run config  json name: " << run_config_json << std::endl;
    std::cout << "globalProb  json name: " << globalProb_json << std::endl;
  }

  // json to pull from
  Json::Value	sample_root;
  Json::Value	selection_root;
  Json::Value	run_config_root;
  Json::Value   globalProb_root;

  // file to input into the json reader
  std::ifstream sample_doc(sample_json, std::ifstream::binary);
  std::ifstream selection_doc(selection_json, std::ifstream::binary);
  std::ifstream run_config_doc(run_config_json, std::ifstream::binary);
  std::ifstream globalProb_doc(globalProb_json, std::ifstream::binary);

  // pointer to the globalProbabilities
  globalJetProbabilities * globalJetProbToApply = NULL; 

  // input the files into the JSON values
  if(debug > -1) std::cout << "Reading in JSONs..." << std::endl;
  if(debug > -1) std::cout << "\t..." << sample_json << std::endl;
  sample_doc     >> sample_root;
  if(debug > -1) std::cout << "\t..." << selection_json << std::endl;
  selection_doc  >> selection_root;
  if(debug > -1) std::cout << "\t..." << run_config_json << std::endl;
  run_config_doc >> run_config_root;

  // read parameters from the run config
  if(debug > -1) std::cout << "Building Output File..." << std::endl;
  std::string	outputFileName = run_config_root.get("outputFileName","test.root").asString();
  std::string	outputDir      = run_config_root.get("outputDir","/tmp/hardenbr/").asString();
  debug			       = run_config_root.get("debug",0).asInt();
  long int	maxEvents      = run_config_root.get("maxEvents",-1).asInt();
  // limit on the number of ntags to compute
  int		maxJetTags     = run_config_root.get("maxJetTags",3).asInt(); 

  if(debug > -1) std::cout << "\t..." << outputFileName << std::endl;
  // read in the json if the config is provided 

  if(debug > -1) std::cout << "Checking for provided global probabilities configuration..." << std::endl;

  if(probProvided) {
    globalProb_doc >> globalProb_root;
    std::cout << "Global Jet Probabilities provided. Parsing JSON... " << std::endl;
  }
  else {
    std::cout << "No global jet probabilities provided. Will compute per sample... " << std::endl;
  }

  if(debug > -1) std::cout << "Building Jet Selection..." << std::endl;
  // build the jet selector
  jetSelector jetSel(selection_root, debug);

  // parse the name of the tree to use for the globalProbabilities
  treeName = selection_root.get("tree","jets").asString();  

  // access the sample parameterization
  const Json::Value samples = sample_root["samples"];
  std::cout << "Number of samples: " << samples.size() << std::endl;

  // output file containing combined information of each sample
  TFile outputFile((outputDir+"/"+outputFileName).c_str(), "RECREATE");

  if(debug > -1) std::cout << "Begin Sample Loop..." << std::endl;
  // loop over each sample and build the global jet probabilities 
  for(int ss = 0; ss < samples.size(); ++ss ) {
    bool runSample  = samples[ss].get("runSample", false).asBool();
    if(!runSample) continue;

    // extract sample informatino from json
    std::string label = samples[ss].get("label", "NO_LABEL").asString();
    std::string path  = samples[ss].get("path", "NO_PATH_PROVIDED").asString();
    std::string stack = samples[ss].get("stack", "NO_STACK_PROVIDED").asString();
    bool	isMC  = samples[ss].get("isMC", true).asBool();
    bool	isSig = samples[ss].get("isSig", false).asBool();
    float		xsec  = samples[ss].get("xsec", 1).asFloat();

    // build  names for tagging histograms  based on the label
    std::string nTagHistTrueName = label + "_nTagTrue";
    std::string nTagHistPredName = label + "_nTagPred";

    // initialize histograms
    const int maxNTags  = 5;
    TH1D nTagHistTrue(nTagHistTrueName.c_str(), "Number of true tags in event", maxNTags - 1, 1, maxNTags);
    TH1D nTagHistPred(nTagHistPredName.c_str(), "Number of predicted tags in event", maxNTags - 1, 1, maxNTags);
    TH1D nTagHistPredErrUp((nTagHistPredName+"ErrUp").c_str(), 
			   "Number of predicted tags in event with combinatorial error varied up", maxNTags - 1, 1, maxNTags);
    TH1D nTagHistPredErrDn((nTagHistPredName+"ErrDn").c_str(), 
			   "Number of predicted tags in event with combinatorial error varried dn", maxNTags - 1, 1, maxNTags);

    if(debug > -1) {
      std::cout << "label: " << label << std::endl;
      std::cout << "\t path: " << path << std::endl;
      std::cout << "\t stack: " << stack << std::endl;
      std::cout << "\t isMC: " << isMC << std::endl;
      std::cout << "\t isSig: " << isSig << std::endl;
      std::cout << "\t xsec: " << xsec << std::endl;
    }

    ///////
    // STEP 1: GENERATE THE JET PROBABILITIES IF NECESSARY
    ///////
    TFile thisFile(path.c_str(),"READ");
    TTree * tree       = (TTree*)(thisFile.Get(treeName.c_str()));    
    TTree * caloHTTree = (TTree*)(thisFile.Get("eventInfo"));    

    // is globalProbabilities were not computed, produce them for each sample
    if(!probProvided) {
      if(debug > -1) std::cout << "Probabilities not provided.." << std::endl;
      if(debug > -1) std::cout << "Constructing localProbabilities..." << std::endl;
      // process the sample by constructing the global probabilities
      globalJetProbabilities * localSampleProb = new globalJetProbabilities(label, stack, isMC, isSig, xsec, tree, jetSel, debug);     

      // parse the json assembled
      Json::Value  sampleJson = localSampleProb->getProbabilitiesJSON();

      // write the json to the stream in the output directory
      std::string jsonOutputName = outputDir+"/"+"json_"+label+".json";
      std::ofstream json_stream;
      json_stream.open(jsonOutputName);      
      // use the style writter
      Json::StyledWriter styledWriter;
      json_stream << styledWriter.write(sampleJson);
      json_stream.close();


      // write the histograms for the file
      if(debug > -1) std::cout << "Writing probability hists..." << std::endl;
      outputFile.cd();      
      localSampleProb->getTaggedHist().Write();
      localSampleProb->getAllHist().Write();
      localSampleProb->getRatioGraph().Write();
      localSampleProb->getCentralEffHist().Write();
      localSampleProb->getCentralEffHistErrUp().Write();
      localSampleProb->getCentralEffHistErrDn().Write();
      
      // set the probabilities to apply to this samples local probabilities
      globalJetProbToApply = localSampleProb; 
    } // end building global probabilities
    else {
      std::cout << "....Loading JSON Probabilities......." << std::endl;
      globalJetProbabilities * loadedProbabilities = new globalJetProbabilities(label, stack, isMC, isSig, xsec, globalProb_root, debug);
      globalJetProbToApply = loadedProbabilities;

      outputFile.cd();
      // write the efficiency histograms
      loadedProbabilities->getCentralEffHist().Write();
      loadedProbabilities->getCentralEffHistErrUp().Write();
      loadedProbabilities->getCentralEffHistErrDn().Write();
      // write the histograms going into the efficiency
      loadedProbabilities->getAllHist().Write();
      loadedProbabilities->getTaggedHist().Write();      

    } // end loading global probaiblities from json

    // check the status of all objects in the global probabilities to by applied
    globalJetProbToApply->printHistStatus();

    ///////////////////
    // STEP 2: APPLY THE GLOBAL / LOCAL PROBABILITIES TO THE SAMPLE 
    // EVENT BY EVENT
    ///////////////////

    if(debug > -1) std::cout << "------ BEGINNING STEP TWO ----- " << std::endl;
    // make one file per sample for the individual tree
    if(debug > -1) std::cout << " Building sample output file  " << std::endl;
    TFile sampleOutputFile((outputDir+"/"+"tree_"+label+".root").c_str(), "RECREATE");    

    // create the tree that will contain the variables from
    // a shallow copy of the original tree
    TTree* jetVariableTree = jetSel.shallowCopyTree(tree);

    // arrays related to the probability tagging    
    int nCat		      = 11;
    int nTagWeight[11]	      = {0,0,0,0,0,0,0,0,0,0};
    double probNTags[11]      = {0,0,0,0,0,0,0,0,0,0};
    double probNTagsErrUp[11] = {0,0,0,0,0,0,0,0,0,0};
    double probNTagsErrDn[11] = {0,0,0,0,0,0,0,0,0,0};
    int   nJetTagArray[11]    = {0,1,2,3,4,5,6,7,8,9,10};
    // jet indexed
    double probJTag[50];
    int   isTagged[50];
    // event indexed / flat number
    int		nTagged = 0;
    long int	evNum	= 0;
    float       caloHT  = 0;
   
    // set the branches for the probabilities
    jetVariableTree->Branch("evNum", &evNum, "evNum/I");
    jetVariableTree->Branch("nCat", &nCat, "nCat/I");
    jetVariableTree->Branch("index", &nJetTagArray, "index[nCat]/I");
    jetVariableTree->Branch("nTagWeight", &nTagWeight, "nTagWeight[nCat]/I");
    jetVariableTree->Branch("probNTags", &probNTags, "probNTags[nCat]/D");
    jetVariableTree->Branch("probNTagsErrUp", &probNTagsErrUp, "probNTagsErrUp[nCat]/D");
    jetVariableTree->Branch("probNTagsErrDn", &probNTagsErrDn, "probNTagsErrDn[nCat]/D");
    // jet indexed
    jetVariableTree->Branch("probJTag", &probJTag, "probJTag[nCaloJets]/D");    
    jetVariableTree->Branch("isTagged", &isTagged, "isTagged[nCaloJets]/I");    
    // event index / flat number
    jetVariableTree->Branch("nTagged", &nTagged, "nTagged/I");    
    jetVariableTree->Branch("caloHT", &caloHT, "caloHT/F");
        
    if(debug > -1) std::cout << "Initializing master probability calculator  " << std::endl;
    // construct the master jet probability calculator to do the n tag calculations
    if(globalJetProbToApply == NULL) std::cout << "global Jet Probabilities not chosen....NULL POINTER" << std::endl;
    jetProbabilityMasterComputer masterJetCalc(globalJetProbToApply, tree, debug);

    // loop event by event
    if(debug > -1) std::cout << " Beginning Event Loop  " << std::endl;
    int nEvents = tree->GetEntries();
    if(maxEvents > 0) nEvents = maxEvents;
    for(long int event = 0; event < nEvents; ++event) {
      if(debug > 2) std::cout << "Checking event selection... "  <<  std::endl;
      bool eventPassSelection = jetSel.doesEventPassSelection(tree, event);
      if(debug > 2) std::cout << "Event passes event selection?? "  << eventPassSelection << std::endl;

      evNum = event;
      if(event % 5000 == 0) std::cout << "Processing Event # --- "  << event << std::endl;
      if(debug > 9) std::cout << "\t\t\tProcessing Event # --- "  << event << std::endl;
      // get the calo HT calculation 
      caloHTTree->GetEntry(event);
      TLeaf *   caloHTLeaf = caloHTTree->GetLeaf("eventCaloHT");
      int       nHT	   = caloHTLeaf->GetNdata();

      if(nHT != 1) {
	std::cout << "caloHT leaf does not contain a single entry" << std::endl;	
	exit(1);
      }
      else {
	caloHT = caloHTLeaf->GetValue(0);
      }

      // get the probability vector for n tagged scenarios
      std::vector<double>    nTagProbVector		    = 
	masterJetCalc.getNJetTaggedVector(event, maxJetTags);
      std::vector<double>    jetProbabilityVector	    = 
	masterJetCalc.getJetProbabilityVector(event);
      std::vector<std::pair<double,double>> nTagProbErrVector = 
	masterJetCalc.getNJetErrorVector(event, maxJetTags);

      if(debug > 2) std::cout << "Getting Vector Size for event: "  << event << std::endl;
      int nJets = jetVariableTree->GetLeaf("nCaloJets")->GetValue(0);//nTagProbVector.size();

      // fill the array with zeros
      for(int jj = 0; jj <= maxJetTags; ++jj) probNTags[jj] = 0;
      
      // look at up to maxJetTags
      for(int jj = 0; (jj <= nJets) && (jj <= maxJetTags); ++jj ) {
	if(debug > 2) std::cout << "Checking for jets jj "  << jj <<  " nJets " << nJets << std::endl;
	double	prob   = nTagProbVector[jj];
	double	probUp = nTagProbErrVector[jj].first;
	double	probDn = nTagProbErrVector[jj].second;

	// fill the arrays given the probabilities make sense
	if(prob == 1 && jj != 0) { 
	  std::cout << "[ERROR] probability is 1?....Exiting....event: " << event << std::endl;
	  std::cout << "[ERROR] nJets: " << nJets << " checking for jet jj = " << jj  << std::endl;
	  std::cout << "[ERROR] probUp: " << probUp << "  probdown: " << probDn  << std::endl;
	  exit(1);
	}

	probNTags[jj]	   = (prob >= 0 && prob <= 1) ? prob : 0;
	probNTagsErrUp[jj] = (probUp >= 0 && probUp <= 1) ? probUp : 0;
	probNTagsErrDn[jj] = (probDn >= 0 && probDn <= 1) ? probDn : 0;

	double weightErrUp = probNTags[jj] + probNTagsErrUp[jj];
	double weightErrDn = probNTags[jj] - probNTagsErrDn[jj];
	// fill the prediciton histogram
	if(debug> 2) {
	  std::cout << "Fill histograms with probabilities p = " << 
	    prob << " pUp =  " << probUp << " pDn" << probDn << 
	    " weight up " << weightErrUp << " weight dn " << weightErrDn << 
	    " jj " << jj << std::endl;
	}

	// fill histograms
	if(eventPassSelection) {
	  nTagHistPred.Fill(jj, probNTags[jj]);
	  nTagHistPredErrUp.Fill(jj, weightErrUp);
	  nTagHistPredErrDn.Fill(jj, weightErrDn);
	}
	
	if(debug > 6) { 
	  std::cout << " histogram status " << std::endl;
	  nTagHistPred.Print();	
	  std::cout << "\tintegral =  " << nTagHistPred.Integral() << std::endl;
	  nTagHistPredErrUp.Print();
	  std::cout << "\tintegral =  " << nTagHistPredErrUp.Integral() << std::endl;
	  nTagHistPredErrDn.Print();
	  std::cout << "\tintegral =  " << nTagHistPredErrDn.Integral() << std::endl;
	}
      }      
      
      // get the per jet tagged vector
      if(debug > 2) std::cout << "Getting the nTagged Jets Vector.." << std::endl;
      std::vector<bool> taggedVector = jetSel.getJetTaggedVector(tree, event);      
      nTagged = 0; // the total number of jets tagged per event      
      if(debug > 2) std::cout << "Filling output Tree Branch.. with tagvector size:" << taggedVector.size() << std::endl;
      if(debug > 5) std::cout << "[jetTagAnalyzer] vector of tags";
      for(int jj = 0; jj < taggedVector.size(); ++jj) {
	if(debug > 5) std::cout << taggedVector[jj];	
	// check for each vector being tagged
	if(taggedVector[jj]) {
	  nTagged++;
	  isTagged[jj] = 1;
	}
	else {
	  isTagged[jj] = 0;
	}
	if(debug > 5) std::cout << std::endl;
	// get the jet probability
	probJTag[jj] = jetProbabilityVector[jj];
      }
      // set the weight vector to the number of true tags 
      nTagWeight[nTagged] = 1;
      // fill the histogram
      if(eventPassSelection) nTagHistTrue.Fill(nTagged,1);
      if(debug > 1) std::cout << " -- Number of jets tagged...." << nTagged << std::endl;

      // fill the Tree
      if(debug > 2) std::cout << "Filling the Tree.." << std::endl;
      jetVariableTree->Fill();

      // reset the weight back to zero
      nTagWeight[nTagged] = 0;
    } // loop over events in a single sample

    sampleOutputFile.cd();

    // build the expectation comparison
    TCanvas canvas("canvas", "prediction vs number of tags per event", 800, 800);

    // lines attributes
    nTagHistPred.SetLineWidth(2);
    nTagHistTrue.SetLineWidth(2);
    nTagHistPredErrUp.SetLineStyle(9);
    nTagHistPredErrDn.SetLineStyle(9);
    nTagHistPred.SetLineColor(kBlue);
    nTagHistTrue.SetLineColor(kBlack);

    // marker  attributes
    nTagHistPred.SetMarkerStyle(25);
    nTagHistTrue.SetMarkerStyle(20);
    nTagHistPred.SetMarkerColor(kBlue);
    nTagHistTrue.SetMarkerColor(kBlack);

    // axes
    nTagHistPredErrUp.GetXaxis()->SetTitle("N Jets Tagged");

    // draw onto the canvas
    nTagHistPredErrUp.Draw();
    nTagHistPredErrDn.Draw("same");
    nTagHistPred.Draw("psame");
    nTagHistTrue.Draw("epsame");

    // write the canvas
    canvas.Write();

    // write the histograms
    nTagHistTrue.Write();
    nTagHistPred.Write();
    nTagHistPredErrUp.Write();
    nTagHistPredErrDn.Write();

    // write out the tree and close 
    jetVariableTree->Write();
    sampleOutputFile.Close();    
  } // loop over samples contained in the the sample json  

  if(debug > -1) std::cout << "Closing outputfile..." << std::endl;
  outputFile.Close();  
} // end of main method



// void processSample(const Json::Value sample) {
  

// }

// void parseSampleJSON() {
// }

// void parseJetSelection() {
// }


// void jetTagAnalyzer::setupOutput() {
// 

// void jetTagAnalyzer::setupProbabilities() {
// }
